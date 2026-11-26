"""
chemwand_mcp.py — MCP server: Claude drives desktop ChemDraw.

Each tool takes SMILES (+ a high-level layout description), builds the figure with
the render-once pipeline, writes a .cdxml, and OPENS it in ChemDraw. Runs on the
Mac; needs rdkit + lxml + numpy + the `mcp` SDK, and the other modules in this
folder (cdxml_core, text_metrics, smiles_to_frag, layout_linear/radial/snaking,
arrows, beautify).

Register it in Claude Desktop, then just ask Claude to draw a molecule / scheme /
radial figure / snaking route.
"""

from __future__ import annotations

import collections
import contextlib
import contextvars
import difflib
import inspect
import json

import descriptors as DESC
import os
import re
import subprocess
import sys
import time
from pathlib import Path
import typing
from typing import Literal

import pydantic

from mcp.server.fastmcp import FastMCP

import beautify as BT
import cdxml_core as C
import decorate as DEC
import graphics as GX
import layout_linear as L
import layout_cycle as CY
import layout_grid as GR
import legend as LG          # aliased: the draw_molecule parameter is also `legend`
import layout_radial as R
# Imported HERE and nowhere else: preflight imports no engine, so the six tools
# validate a spec without a layout module having been touched.
import preflight as PF
import smiles_to_frag as SF
from smiles_to_frag import smiles_to_fragment
import stereo_tools as ST
import table_layout as TL   # only for CFG's key set -- see the schema section

mcp = FastMCP("chemwand")

OUT_DIR = Path(os.environ.get("CHEMWAND_FIG_DIR", str(Path.home() / "chemwand_figures")))
OUT_DIR.mkdir(parents=True, exist_ok=True)
# The bare name "ChemDraw" matches no installed application -- it is installed as
# "ChemDraw 23.1.2.app" -- and the two call sites fail DIFFERENTLY, which is what
# hid it. _open_in_chemdraw catches the failure and falls back to `open <path>`, so
# figures still opened and nothing looked wrong; _clean_up's AppleScript has no
# fallback and is wrapped in `except: pass`, so clean_up=True did nothing at all.
#
# The ".app" is load-bearing and the two callers disagree about it, so it was
# measured rather than guessed:
#   open -a          accepts "ChemDraw 23.1.2.app", REJECTS "ChemDraw 23.1.2"
#   AppleScript      accepts both
# so the suffixed form is the only single value that satisfies both. A site with a
# different build sets CHEMDRAW_APP; include the ".app".
CHEMDRAW_APP = os.environ.get("CHEMDRAW_APP", "ChemDraw 23.1.2.app")

# Structure bottom -> legend top. Chosen off a rendered figure, and it lands
# inside a window the drawing already defines:
#
#   num_gap      10   structure -> its own compound number. The tightest a text
#                     block ever sits under a structure, because that number
#                     BELONGS to it.
#   LEGEND_GAP   16
#   footnote_gap 22   scope table -> its running prose. The loosest, because a
#                     footnote is an aside and reads as a separate block.
#
# A legend is neither: it is a second half, not a caption and not an aside. 22
# was tried first and read as two unrelated objects on the page.
LEGEND_GAP = 16.0

# ---- which build of this code is actually LOADED --------------------------
#
# This server does not hot-reload. Edit a .py and every tool call keeps running the
# code that was imported when Claude Desktop started, silently, with no outward sign
# -- so a fix can easily be tested against a version that was never loaded. So
# the results carry a stamp: the commit they came from, and when the
# module was imported. Cheap insurance against re-litigating a bug that was already
# fixed on disk.
_HERE = Path(__file__).resolve().parent
_IMPORTED_AT = time.time()


def _git(*args: str) -> str:
    """`git ...` in this file's directory, or "" for any reason at all. A version
    stamp must never be the thing that stops a figure from being drawn."""
    # A deployed copy has no .git, and on a Mac without the developer tools
    # /usr/bin/git is a stub that opens an install dialog nobody is there to
    # click. No repository, no call.
    if not (_HERE / ".git").exists():
        return ""
    try:
        r = subprocess.run(("git",) + args, cwd=str(_HERE), capture_output=True,
                           text=True, timeout=5)
        return r.stdout.strip() if r.returncode == 0 else ""
    except Exception:                                            # noqa: BLE001
        return ""


def _commit() -> str:
    rev = _git("rev-parse", "--short", "HEAD")
    if not rev:
        return "no-git"
    # -uno: a tracked file edited but not committed is exactly the state worth
    # flagging; a stray scratch file in the directory is not.
    return rev + "-dirty" if _git("status", "--porcelain", "-uno") else rev


_COMMIT = _commit()          # once, at import: this is the code that is running


def _mmss(seconds: float) -> str:
    m, s = divmod(int(seconds), 60)
    return "%dh%02dm" % divmod(m, 60) if m >= 60 else "%dm%02ds" % (m, s)


def engine_stamp() -> str:
    """The stamp that goes back with a result.

    Resolved at import, EXCEPT the staleness check -- that has to run now, because
    the whole point is to catch a source file that changed after the import it is
    being compared against. Answers the question rather than making the reader
    compare two timestamps by eye.

    EVERY stale module is named, not just the most recently touched one. Reporting
    only the newest file answers "has anything changed?" when the question actually
    being asked is "which of the running modules is not the code just edited?" --
    and edits to the allocator, then the layout, then text_metrics would otherwise
    report text_metrics alone, leaving the other two silently stale with nothing
    to distinguish them from up-to-date.
    """
    loaded = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(_IMPORTED_AT))
    # 2s of grace: a file written during the import itself is not evidence of
    # anything, and clock resolution should not be able to raise a false alarm.
    cutoff = _IMPORTED_AT + 2.0
    stale = []
    for mod in list(sys.modules.values()):
        f = getattr(mod, "__file__", None) or ""
        if not f.endswith(".py"):
            continue
        p = Path(f)
        if p.parent != _HERE:            # only the modules a restart would reload
            continue
        try:
            mt = p.stat().st_mtime
        except OSError:
            continue
        if mt > cutoff:
            stale.append((mt, p.name))
    stamp = "%s loaded %s" % (_COMMIT, loaded)
    if stale:
        stale.sort(key=lambda t: -t[0])          # most recently edited first
        stamp += (" -- STALE: %d loaded module(s) changed after import (%s); quit "
                  "Claude Desktop (Cmd+Q) and reopen it or this is not the code you "
                  "just edited"
                  % (len(stale),
                     ", ".join("%s +%s" % (name, _mmss(mt - _IMPORTED_AT))
                               for mt, name in stale)))
    return stamp


def _open_in_chemdraw(path: Path) -> None:
    if sys.platform == "win32":          # Windows: open with the default .cdxml app
        try:
            os.startfile(str(path))
        except Exception:
            pass
        return
    try:
        subprocess.run(["open", "-a", CHEMDRAW_APP, str(path)], check=True)
    except Exception:
        subprocess.run(["open", str(path)], check=False)   # fall back to default .cdxml handler


def _clean_up() -> None:
    """Optional: Select All + Structure > Clean Up Structure (Cmd+Shift+K).
    Best-effort GUI automation; needs macOS Accessibility permission for the host app."""
    script = (
        f'tell application "{CHEMDRAW_APP}" to activate\n'
        'delay 0.6\n'
        'tell application "System Events"\n'
        '  keystroke "a" using command down\n'
        '  delay 0.2\n'
        '  keystroke "k" using {command down, shift down}\n'
        'end tell\n'
    )
    try:
        subprocess.run(["osascript", "-e", script], check=False, timeout=12)
    except Exception:
        pass


# WHY headless() EXISTS. Every draw_* tool funnels its finished document through
# _write_and_open, so that one call is the whole boundary between "compute the
# CDXML" -- arithmetic, runs anywhere -- and "show it on THIS Mac" -- a subprocess
# `open -a` that needs a GUI and a ChemDraw license. A server has the first half
# and none of the second. The test suite has been cutting in exactly this place
# since the first headless test (eight files monkeypatch _write_and_open or
# _open_in_chemdraw); headless() only makes that cut official and re-entrant.
#
# THE FILE IS STILL WRITTEN, only the opening is skipped. _figure_report re-reads
# the path to run the InChI round trip -- `roundtrip_inchi.check(path, ...)`, and
# note its `not path.startswith("/dev/")` guard, which is what a suppressed write
# would trip. Dropping the write would silently drop the one check that catches
# "the figure came out, but it is not the molecule that was asked for". Servers
# have disks; what they lack is a screen.
#
# ContextVar and not a plain global: FastAPI runs sync endpoints in a thread pool,
# so two users drawing at the same time are two threads inside this module at the
# same time. A module-level flag would leak one request's headlessness into the
# other's, and the symptom -- a figure that silently did not open -- looks like
# nothing at all on a desktop.
_HEADLESS: contextvars.ContextVar[bool] = contextvars.ContextVar(
    "chemdraw_headless", default=False)

# WHY THE DIRECTORY IS ALSO PER-CONTEXT. The filename is Claude's to choose and it
# chooses descriptively, so two people asking for the same reaction get the same
# name -- suzuki_coupling.cdxml, twice, into one directory. The second write wins
# and whoever downloads first gets a figure of someone else's molecule, with no
# error on either side and nothing in the report that looks wrong. One directory
# per caller removes the shared name entirely, which is cheaper than negotiating
# over it and does not depend on Claude naming things uniquely.
_OUT_DIR: contextvars.ContextVar = contextvars.ContextVar(
    "chemdraw_out_dir", default=None)


@contextlib.contextmanager
def headless(out_dir=None):
    """Run draw_* tools without opening ChemDraw (and without clean_up, which is
    GUI automation over the opened window). The .cdxml is still written, and its
    path still comes back in the report's `path` field -- read the CDXML there.

    `out_dir` redirects the write. It is created if missing, so a caller can hand
    over a directory per caller without arranging for it first."""
    tokens = [_HEADLESS.set(True)]
    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        tokens.append(_OUT_DIR.set(out_dir))
    try:
        yield
    finally:
        for t in reversed(tokens):
            t.var.reset(t)


def _footnote_lines(footnotes) -> list:
    """`footnotes` as a list of lines, whatever shape it arrived in.

    THE BUG THIS EXISTS FOR: a finished two-step scheme whose footnote reads,
    verbatim, brackets and quotes and all:

        ["Regioselectivity caveat (step 1): SNAr on 2,4,6-trichloropyrimidine
        is not fully site-selective...", "Step 2 is a routine acid chloride..."]

    The tools take `list | str`, the model writes the list correctly nearly
    always, and occasionally sends the same list SERIALIZED -- a string whose
    content is a JSON array. `[footnotes] if isinstance(footnotes, str)` then
    makes it one footnote, and that one footnote is JSON. Nobody ever meant to
    print brackets under a figure, so a JSON array that arrived as a string is
    read as the array it plainly is.

    `[a] Conditions: ...` IS A REAL FOOTNOTE and must survive untouched. It does,
    and json.loads IS THE GUARD THAT MAKES IT SURVIVE -- not the bracket test.
    Requiring a trailing `]` as well as a leading `[` would change no behaviour,
    because `[a] Conditions: ...` does not parse either way. The bracket test is
    a cheap way to skip json.loads on the ordinary bare-string footnote, and
    nothing more; it is deliberately only the leading `[`, so it claims no credit
    it has not earned.
    Anything reaching json.loads must come back a list of strings or it is left
    exactly as it arrived.

    list() is not called on a bare string anywhere here, which is the hazard of
    the naive per-tool form: list("a note") is ['a', ' ', 'n', ...], one
    footnote per character."""
    if not footnotes:
        return []
    if isinstance(footnotes, str):
        one = footnotes.strip()
        if one.startswith("["):        # cheap skip, not the guard -- see above
            try:
                parsed = json.loads(one)
            except ValueError:
                parsed = None
            if isinstance(parsed, list) and all(isinstance(x, str) for x in parsed):
                return [t for t in (x.strip() for x in parsed) if t]
        return [one] if one else []
    return [str(x) for x in footnotes]


# Characters XML 1.0 does not allow at all. Tab, newline and carriage return are
# allowed and stay. Anything else here makes a file ChemDraw cannot open, and the
# InChI read-back then failed on it and blamed the structures (or, for a figure
# with no structures, said nothing).
_XML_BAD = re.compile("[\x00-\x08\x0b\x0c\x0e-\x1f\ud800-\udfff￾￿]")


def _refuse_bad_chars(doc: str, filename: str) -> None:
    """Raise, naming the text, if the figure holds a character XML forbids.

    Refused rather than replaced: \\x0b is Word's soft line break, so a space
    would quietly change what the caller wrote."""
    m = _XML_BAD.search(filename)
    if m:
        raise ValueError(
            "The filename %r contains an invisible control character (U+%04X). "
            "Nothing was written. Remove it and try again."
            % (filename, ord(m.group())))
    m = _XML_BAD.search(doc)
    if not m:
        return
    # The run of text it sits in: from the tag before it to the tag after --
    # or, inside a tag (an R-group label is an attribute), the quoted value.
    i = m.start()
    if doc.rfind("<", 0, i) > doc.rfind(">", 0, i):
        lo, hi = doc.rfind('"', 0, i) + 1, doc.find('"', m.end())
    else:
        lo, hi = doc.rfind(">", 0, i) + 1, doc.find("<", m.end())
    text = doc[lo:hi if hi >= 0 else len(doc)]
    # Marked before stripping: strip() would eat \x0b / \x0c at either end.
    shown = _XML_BAD.sub(lambda c: "[U+%04X]" % ord(c.group()), text).strip()[:120]
    raise ValueError(
        "A piece of text in this figure contains an invisible control character "
        "(U+%04X): \"%s\". Nothing was written -- ChemDraw cannot open a file "
        "with it. It usually comes from text pasted out of Word or a PDF; remove "
        "it and draw again." % (ord(m.group()), shown))


# A numeric attribute (p, BoundingBox, D, ...) holding nan or inf: a position
# that was never computed. Matched only in values made of numbers, so a word
# such as "ChainAngle" or a caption cannot trip it.
_NONFINITE = re.compile(
    r'\s([A-Za-z0-9]+)="((?:[-+0-9.eE ]|nan|inf)*\b(?:nan|inf)\b(?:[-+0-9.eE ]|nan|inf)*)"')


def _refuse_nonfinite(doc: str) -> None:
    """Raise if any coordinate in the figure is NaN or infinite.

    Such a file used to be written with ok:true, and reading it back for the
    InChI check could run RDKit out of memory."""
    m = _NONFINITE.search(doc)
    if not m:
        return
    tag = doc.rfind("<", 0, m.start())
    raise ValueError(
        "ChemWand could not compute a position in this figure (%s=\"%s\" on "
        "%s). Nothing was written. Please report the input; as a workaround, "
        "write the species that caused it as text instead of as a structure."
        % (m.group(1), m.group(2), doc[tag:tag + 40].split(">")[0] + ">"))


_NAME_MAX = 100     # characters of a file name, before ".cdxml"


def _clean_name(filename: str, tool: str) -> str:
    """A file name that is safe to write and to turn into a download link.

    Only spaces and slashes used to be handled. An empty name wrote a hidden
    ".cdxml"; a newline or tab gave a name the link builder cannot match, so the
    reply carried the server's real path; 300 characters failed with that path
    in the error; "x.cdxml" became "x.cdxml.cdxml"."""
    name = re.sub(r"[\s/\\]+", "_", str(filename or "")).strip("_")
    name = re.sub(r"(?i)\.cdxml$", "", name).lstrip(".")[:_NAME_MAX]
    if name:
        return name
    # Empty after cleaning: the calling tool's own default ("molecule", ...).
    try:
        return inspect.signature(globals()[tool]).parameters["filename"].default
    except Exception:
        return "figure"


def _write_and_open(doc: str, filename: str, clean_up: bool) -> str:
    _refuse_bad_chars(doc, filename)
    _refuse_nonfinite(doc)
    base = _OUT_DIR.get() or OUT_DIR
    stem = _clean_name(filename, sys._getframe(1).f_code.co_name)
    # NEVER OVERWRITE. Every tool has a fixed default name, so two figures drawn
    # without one were both "molecule.cdxml": the second replaced the first, and
    # the link already handed out for the first downloaded the second. A taken
    # name gets _2, _3, ...; "x" mode makes the claim atomic, so two calls in
    # one folder at the same moment cannot both take the same number.
    n = 1
    while True:
        path = base / ("%s.cdxml" % stem if n == 1 else "%s_%d.cdxml" % (stem, n))
        try:
            with open(path, "x", encoding="utf-8") as fh:
                fh.write(doc)
            break
        except FileExistsError:
            n += 1
    if _HEADLESS.get():
        return str(path)
    _open_in_chemdraw(path)
    if clean_up:
        _clean_up()
    return str(path)


# ---- what the SMILES -> structure layer could not do ----------------------
#
# Every one of these has the same shape: the figure came out, but a molecule in it
# is not the molecule that was asked for -- drawn full-length when abbreviation was
# requested, or oriented differently from the rest of the series. None of that is
# visible in a filename, and a table of twenty products is exactly where nobody
# notices. So the tool SAYS it, next to the path, every time.
def _chem_report() -> str:
    """Notes for the molecules this batch could not draw as asked ("" if none).

    Call _new_batch() before converting and this after; the reports live in
    smiles_to_frag and are per-batch, not per-molecule.
    """
    out = []
    deloc = list(DEC.DELOCALIZED)
    if deloc:
        out.append(
            "\nNOTE: %d structure(s) are drawn with a DELOCALIZED charge -- the "
            "double bonds are deliberately absent and a dashed arc plus one "
            "floating symbol stand in for them. That makes the file a DRAWING and "
            "not the molecule, so the structure/inchi-mismatch warning below is "
            "EXPECTED for them and is not suppressed: %s"
            % (len(deloc), ", ".join(str(d) for d in deloc[:5])))
    misses = list(SF.TEMPLATE_MISSES)
    if misses:
        out.append("\nNOTE: %d structure(s) did not match the depiction template and "
                   "were laid out independently, so their scaffold may be oriented "
                   "differently: %s" % (len(misses), ", ".join(misses[:5])))
    bad = list(SF.TEMPLATE_DISTORTED)
    if bad:
        out.append("\nNOTE: %d structure(s) MATCHED the template but came out of the "
                   "constrained layout with uneven bonds, so the template was dropped "
                   "for them and they were laid out independently (longest/shortest "
                   "bond, templated vs on its own): %s"
                   % (len(bad), ", ".join("%s %.2fx vs %.2fx" % (s, t, f)
                                          for s, t, f in bad[:5])))
    drift = list(getattr(SF, "TEMPLATE_DRIFTED", []))
    if drift:
        out.append("\nNOTE: %d structure(s) MATCHED the template and came through "
                   "undeformed, but with the shared scaffold rotated or shifted away "
                   "from where the rest of the series draws it -- each cell looks "
                   "right on its own and crooked next to its neighbours -- so the "
                   "template was dropped for them (drift, in bond lengths): %s"
                   % (len(drift), ", ".join("%s %.2f" % (s, d) for s, d in drift[:5])))
    ezr = list(getattr(SF, "TEMPLATE_EZ_REDRAWN", []))
    if ezr:
        out.append("\nNOTE: %d structure(s) MATCHED the template, but copying its "
                   "layout drew a ring double bond the wrong way round (a trans ring "
                   "alkene squared into a cis one), so the template was dropped for "
                   "them and they were laid out independently: %s"
                   % (len(ezr), ", ".join(ezr[:5])))
    ph = list(getattr(SF, "PLACEHOLDERS", []))
    if ph:
        labels = sorted({l for _s, ls in ph for l in ls})
        out.append("\nNOTE: %d structure(s) carry placeholder labels (%s), so this "
                   "figure is a Markush DRAWING, not a set of molecules -- the "
                   "structures cannot be searched and their formulas are "
                   "meaningless. Give the values in a legend."
                   % (len(ph), ", ".join(labels)))
    miss = list(getattr(SF, "ANCHOR_MISSES", []))
    if miss:
        out.append("\nNOTE: %d highlight anchor(s) matched nothing, so nothing was "
                   "drawn for them -- a figure missing a highlight looks exactly "
                   "like one that was never asked for it: %s"
                   % (len(miss), "; ".join(
                       "%s on %s%s" % (sma, s, " -- " + why if why else "")
                       for s, sma, why in miss[:5])))
    multi = list(getattr(SF, "ANCHOR_MULTI", []))
    if multi:
        out.append("\nNOTE: %d anchor(s) matched several times AND each match names "
                   "a GROUP rather than one atom, so a single spec drew several "
                   "marks per match -- \"the amide\" written as [#6](=[#8])[#7] put "
                   "NINE dots on a tripeptide. Narrow it with a recursive SMARTS "
                   "([$(...)]) or address the atoms with 'maps': %s"
                   % (len(multi), "; ".join(
                       "%s on %s -- %d matches, %d atoms, style %s"
                       % (sma, s, nm, na, st) for s, sma, st, nm, na in multi[:5])))
    nofill = list(getattr(SF, "FILL_NOT_RING", []))
    if nofill:
        out.append("\nNOTE: %d 'fill' anchor(s) matched something with no ring in "
                   "it, and a native fill is an area bounded by the bonds it "
                   "names -- a chain or a lone atom encloses nothing. Those were "
                   "drawn as style 'color' (tinted atoms and bonds) instead, which "
                   "is a different picture from the one asked for: %s"
                   % (len(nofill), "; ".join("%s on %s" % (sma, s)
                                             for s, sma in nofill[:5])))
    sqo = list(getattr(SF, "SQUIGGLE_OVERFLOW", []))
    if sqo:
        out.append("\nNOTE: %d squiggle(s) reach outside the structure's own "
                   "bounding box, by up to %.1f units. That overhang is what "
                   "makes the mark read as a cut, so it is not an error -- but "
                   "the box a packer reserves for a cell is the structure's, so "
                   "the overhang is spent out of the gap to the neighbour. If it "
                   "touches something, widen the gap or lower 'overhang': %s"
                   % (len(sqo), max(d for _s, _sma, d in sqo),
                      "; ".join("%s on %s (%.1f)" % (sma, s, d)
                                for s, sma, d in sqo[:5])))
    flo = list(getattr(SF, "FLOATING", []))
    if flo:
        out.append("\nNOTE: %d structure(s) draw a placeholder on a FLOATING bond "
                   "(%s) -- it crosses a ring edge and stops inside, so the "
                   "picture deliberately does not say which position is "
                   "substituted. The bond in the file no longer runs to the atom "
                   "the SMILES attached it to -- it ends on a bare node inside the "
                   "ring, which ChemDraw counts as a CARBON, so the structure "
                   "carries a methyl that nothing draws. Right on paper, wrong as "
                   "a connection table. Put the real positions in the legend "
                   "values (6-Et, 7-Me, 6,7-OMe)."
                   % (len(flo), ", ".join(sorted({l for _s, ls in flo for l in ls}))))
    noflo = list(getattr(SF, "FLOAT_SKIPPED", []))
    if noflo:
        out.append("\nNOTE: %d placeholder(s) were asked to float and could not, "
                   "so they are drawn on an ordinary bond to the atom named in "
                   "the SMILES: %s"
                   % (len(noflo), ", ".join("%s %s (%s)" % (s, l, why)
                                            for s, l, why in noflo[:5])))
    unused = list(getattr(SF, "RLABEL_UNUSED", []))
    if unused:
        out.append("\nNOTE: %d rlabels key(s) matched no placeholder in the SMILES, "
                   "so the default label was drawn instead -- check for a typo "
                   "(rlabels keys are the map numbers in [*:n]): %s"
                   % (len(unused), ", ".join("%s wanted %r" % (s, k)
                                             for s, k in unused[:5])))
    if getattr(SF, "TEMPLATE_MARKUSH", []) and misses:
        out.append("\nNOTE: the depiction template is itself a Markush structure "
                   "(%s). A placeholder only matches another placeholder, so a "
                   "Markush parent cannot align EXPLICIT products -- that is why "
                   "they missed. Either give every cell the same placeholder, or "
                   "pass template='<an explicit member>'."
                   % ", ".join(SF.TEMPLATE_MARKUSH[:3]))
    if _PALETTE_RESCALED:
        out.append("\nPALETTE: %d colour(s) were given with a component above 1, "
                   "so they were read as 0-255 (bytes) and divided down to the "
                   "0..1 fractions CDXML stores: %s. If you meant fractions, "
                   "those values were already out of range."
                   % (len(_PALETTE_RESCALED),
                      "; ".join("%s -> %s" % (a, b)
                                for a, b in _PALETTE_RESCALED[:5])))
    if _PALETTE_SIZE:
        _added = sum(_PALETTE_SIZE)
        out.append("\nPALETTE: %d custom colour(s) were appended to this "
                   "document's colortable, which now holds %d entries -- the "
                   "eight standard ones at indices 2-9 and yours at %s. That "
                   "index is what a later call has to name to reuse the same "
                   "colour."
                   % (_added, 8 + _added,
                      "index 10" if _added == 1 else "10-%d" % (9 + _added)))
        # NO HARD CAP, because the measurement does not support one. A ladder of
        # test files (32, 64, 128, 256 and 512 appended colours, sixteen
        # hues repeating so a broken row is visible without judging any colour)
        # was opened in ChemDraw 26: every row rendered correctly at every size
        # up to and including 512 appended colours. So <=512 is verified and
        # anything past it is simply unmeasured -- which is a warning, not a
        # refusal. Raise the number here by running one more rung of the ladder,
        # not by reasoning about it.
        if _added > 512:
            out.append("\nWARNING: %d appended colours is more than the 512 that have been "
                       "verified to render correctly. Beyond that it is untested, not known to "
                       "be broken; if colours in this figure come out wrong, this is the first "
                       "thing to suspect."
                       % _added)
    if _GRAPHICS_BOX:
        # The ONE number a caller needs to place anything on the next call, and
        # it cannot be known before the figure is laid out. Reported every time
        # graphics are used, not only on trouble: the workflow is two passes, and
        # the second one is impossible without it.
        out.append("\nGRAPHICS: this figure's content box is %.0f x %.0f pt -- "
                   "what the ENGINE drew, i.e. the rows block. Your graphics "
                   "coordinates are measured from ITS TOP-LEFT corner, so "
                   "[0, 0, %.0f, %.0f] is that block. It does NOT frame everything "
                   "on the page: a slot placed with 'of' / 'side' sits OUTSIDE it, "
                   "and one on the left is at negative x."
                   % (_GRAPHICS_BOX[0], _GRAPHICS_BOX[1],
                      _GRAPHICS_BOX[0], _GRAPHICS_BOX[1]))
    if _SLOTS_MOVED:
        out.append(
            "\n" + "NOTE: %d auto-placed slot(s) needed more room than the usual "
            "gap and were pushed further out: %s. Nothing is wrong with the figure; "
            "if the ligand box now sits oddly far from it, that is because something "
            "else already occupied the near spot."
            % (len(_SLOTS_MOVED),
               "; ".join("%s moved %.0f pt %s (%s)" % (w, d, dirn, how)
                         for w, dirn, d, how in _SLOTS_MOVED[:4])))
    if _PLACED_MOVED:
        out.append(
            "\n" + "NOTE: %d placed graphic(s) did not fit at the usual gap and were pushed "
                   "further out: %s."
            % (len(_PLACED_MOVED),
               "; ".join("%s moved %.0f pt %s (%s)" % (w, d, dirn, how)
                         for w, dirn, d, how in _PLACED_MOVED[:4])))
    if _GRAPHICS_COVERS:
        # Reported, never repositioned: the caller picked that spot on purpose and
        # this layer has no idea what for. Naming WHAT it lands on is the point --
        # "a slot overlaps something" sends you back to counting coordinates, and
        # the owner strings the engines already keep ("cell:3a", "hdr:arrow.above")
        # say it outright.
        out.append(
            "\n" + "NOTE: %d caller-placed graphic(s) land on top of the figure. "
            "The file was written as asked -- nothing was moved -- but a filled or "
            "opaque one HIDES what is under it:\n%s"
            % (len(_GRAPHICS_COVERS),
               "\n".join("      %s covers %s (%s) by %.1f pt"
                          % (w, o, k, d) for w, o, k, d in _GRAPHICS_COVERS[:6])))
    if _GRAPHICS_OUTSIDE:
        box, page = _GRAPHICS_OUTSIDE[0]
        out.append("\nNOTE: a graphic reaches outside the page (drawn "
                   "[%.0f %.0f %.0f %.0f], page [%.0f %.0f %.0f %.0f]). It is in "
                   "the file and ChemDraw will show it, but it will not print. "
                   "Keep inside the content box above, or leave room in `page`."
                   % (tuple(box) + tuple(page)))
    clashes = list(_ID_CLASHES)
    if clashes:
        out.append(
            "\nNOTE: %d cell(s) named an id another cell already held, and a "
            "figure has one structure per id -- so THIS SMILES WAS NOT DRAWN and "
            "the first cell's structure appears in both places: %s. Two cells "
            "SHARING an id is fine when they are the same compound (a struck-out "
            "target and the product both labelled 8a); these differ. Give the "
            "second one its own id, or drop its 'smiles' if it was meant to "
            "repeat the first."
            % (len(clashes),
               "; ".join("id %r kept %s, dropped %s" % (i, k, d)
                         for i, k, d in clashes[:5])))
    none_hit = list(SF.ABBREV_MISSES)
    if none_hit:
        out.append("\nNOTE: abbreviate was on, but %d structure(s) carry no group in "
                   "the abbreviation table, so they are drawn in full: %s"
                   % (len(none_hit), ", ".join(none_hit[:5])))
    expanded = list(getattr(SF, "ABBREV_EXPANDED", []))
    if expanded:
        # Summed over the batch and named BY GROUP. The per-structure report above
        # cannot see this case at all: it fires only on a total strikeout, so a
        # figure whose every cell contracted a Boc and drew a phenyl ring out in
        # full came back clean. Three expanded rings is the difference between a
        # 552 pt column and a 708 pt scheme.
        tally: dict = {}
        for _smi, label, n in expanded:
            tally[label] = tally.get(label, 0) + n
        worst = sorted(tally.items(), key=lambda kv: (-kv[1], kv[0]))
        in_table = {l for l, _s in SF._ABBREV_SMARTS}
        absent = [l for l, _n in worst if l not in in_table]
        out.append("\nNOTE: abbreviate was on, and these groups were drawn in FULL "
                   "anyway: %s. Each expanded group costs the figure width."
                   % ", ".join("%s x%d" % (l, n) for l, n in worst))
        if absent:
            out.append(" Not in the abbreviation table at all, so they cannot be "
                       "contracted yet: %s." % ", ".join(absent))
    dropped = list(getattr(SF, "STEREO_GROUP_DROPPED", []))
    if dropped:
        # THE ONLY REPORT HERE WHOSE FAILURE IS INVISIBLE ON PAPER. A racemic-
        # but-defined product is drawn with ordinary wedges -- the same wedges a
        # single enantiomer gets -- and the AND/OR group is the only thing in the
        # file that says which of the two it is. Lose the group and the drawing
        # does not change by one pixel; what changes is that the structure starts
        # asserting a single enantiomer, which is different chemistry, and the
        # assertion travels wherever the file goes.
        #
        # A REPORT RATHER THAN A FIX, because there is nowhere to put it back:
        # EnhancedStereoType is an attribute on a plain <n>, an atom that
        # contracted into a superatom lives inside that superatom's nested
        # fragment, and an atom drawn as Markush is a different node -- neither
        # can carry one. The wedge on such a bond CAN be carried and now is; this
        # cannot, and the difference is the file format rather than an omission.
        out.append("\nNOTE: %d racemic/relative stereo group(s) (AND / OR) had "
                   "nowhere to land and were NOT written: %s. Those centres now "
                   "read as a SINGLE ENANTIOMER in the file while the drawing is "
                   "unchanged -- the wedges look identical either way, so nothing "
                   "on paper shows this. Keep that centre a plain atom (no "
                   "abbreviation, no placeholder), or put the group on a centre "
                   "that stays one."
                   % (len(dropped),
                      "; ".join("atom %d (%s) %s" % (i, g, why)
                                for _s, i, g, why in dropped[:5])))
        if len(dropped) > 5:
            out.append(" (%d more.)" % (len(dropped) - 5))
    if getattr(SF, "_DL_ERROR", ""):
        out.append("\nNOTE: the hand-drawn depiction library did not load, so every "
                   "molecule was computed -- including any that SMILES cannot express: "
                   "%s" % SF._DL_ERROR)
    DL = getattr(SF, "_DL", None)
    if DL is not None:
        canned = list(DL.CANNED_USED)
        if canned:
            out.append("\nNOTE: %d structure(s) came from the hand-drawn depiction "
                       "library instead of being computed, because SMILES cannot "
                       "express them: %s"
                       % (len(canned), ", ".join("%s -> %s" % (s, n)
                                                 for s, n in canned[:5])))
        # A CANNED CORE IS OFFERED AT LOOKUP TIME AND CAN STILL BE THROWN AWAY.
        # core_lookup records the offer; after that the canned reference faces the
        # same gates as any template (no match, uneven bonds, drift, a ring E/Z
        # drawn the wrong way round), and it does get rejected -- vancomycin's
        # glucosamine matches and is then dropped for dispersion 2.343 against the
        # 1.39 gate, correctly, because keeping it stretches a bond to 24.6 pt. So
        # the offers are netted against those lists here rather than reported as
        # successes.
        _lost = (set(SF.TEMPLATE_MISSES)
                 | {t[0] for t in SF.TEMPLATE_DISTORTED}
                 | {t[0] for t in getattr(SF, "TEMPLATE_DRIFTED", [])}
                 | set(getattr(SF, "TEMPLATE_EZ_REDRAWN", [])))
        core_used = [(s, n) for s, n in DL.CANNED_CORE_USED if s not in _lost]
        core_lost = [(s, n) for s, n in DL.CANNED_CORE_USED if s in _lost]
        if core_used:
            out.append("\nNOTE: %d structure(s) had PART of the drawing taken from "
                       "the hand-drawn depiction library (%s): the molecule is still "
                       "the one asked for, but that piece of it is a chemist's "
                       "picture rather than a computed one -- a sugar arrives as its "
                       "chair, not as a flat ring: %s"
                       % (len(core_used), ", ".join(sorted({n for _, n in core_used})),
                          ", ".join(s for s, _ in core_used[:5])))
        if core_lost:
            out.append("\nNOTE: %d structure(s) CONTAIN a hand-drawn entry (%s) and "
                       "were offered its coordinates, but the layout gates dropped it "
                       "again, so they were computed after all. THE 'DEPICTION "
                       "TEMPLATE' IN THE NOTES ABOVE IS THAT ENTRY -- no template was "
                       "passed for them: %s"
                       % (len(core_lost), ", ".join(sorted({n for _, n in core_lost})),
                          ", ".join(s for s, _ in core_lost[:5])))
        # Minus EVERY molecule the two notes above already covered, kept or lost.
        # Without the "kept" half the same molecule collects two notes that
        # contradict each other, one saying its ring came from the library and one
        # saying it did not -- measured on 6-O-acetyl-D-glucose and maltose, where
        # lookup() sees the whole entry as a substructure AND core_lookup then
        # places it. Without the "lost" half it collects two notes that agree, the
        # second adding nothing: on vancomycin this block repeated, in weaker
        # words, what the gate note directly above had just said with numbers.
        # What is left is the molecule that was never OFFERED a core -- a
        # metallocene, or a sugar whose cell was pinned to a template the caller
        # passed, which outranks the library.
        _got = {s for s, _ in DL.CANNED_CORE_USED}
        core = [(s, n) for s, n in DL.CANNED_CORE_ONLY if s not in _got]
        if core:
            out.append("\nNOTE: %d structure(s) CONTAIN a molecule the depiction "
                       "library has (%s) but carry substituents as well, so they were "
                       "computed by rdkit instead. Whatever the hand drawing had that "
                       "a computed one need not -- a coordination bond for a "
                       "metallocene, a chair for a sugar -- is absent from those "
                       "cells: %s"
                       % (len(core), ", ".join(sorted({n for _, n in core})),
                          ", ".join(s for s, _ in core[:5])))
    return "".join(out)


# Cells whose SMILES was thrown away because another cell already held their id:
# (id, kept smiles, discarded smiles). See _register for why same-id-same-smiles
# is deliberately NOT in here.
_ID_CLASHES: list[tuple[str, str, str]] = []

# How many cells were HANDED to _register this batch, counted before the
# duplicate check. `len(structures)` counts what came out. The two differ
# exactly when a cell was merged into an earlier one, which is the number the
# caller could never see before.
_SUBMITTED: list[int] = [0]


def _layout_report(alloc) -> str:
    """The receipt's layout section, or "" when there is nothing to say.

    TWO REASONS TO SPEAK, NOT ONE. Findings are the obvious one. The other is the
    exemption tally: probe() counts every overlap it decided not to report because
    both sides belong to one object, and allocator.py's own comment at that branch
    says why the number has to be visible -- "0 findings, 0 exemptions" and "0
    findings, 214 exemptions" mean very different things and must not look identical.
    Gating on findings alone cancelled that in all five tools at once: a 20-node
    cycle with 201 measured overlaps returned ok:true and a receipt that said
    nothing whatsoever.

    ONLY THE ENDPOINT KIND IS COUNTED, WHICH IS WHY THERE IS NO THRESHOLD HERE.
    related() exempts on two grounds and they behave nothing alike:

      SAME OBJECT   a thing over its own part -- `cell:H~2~O` over
                    `cell:H~2~O.text` at 8.1 deep, `center.box` over
                    `center.label` at 19.9. Deliberate containment, and it
                    SCALES WITH THE FIGURE: 2 / 4 / 8 text cells in a scheme
                    give 1 / 3 / 7 of them. Any count threshold over this kind
                    fires on a long scheme that is perfectly fine, which is
                    what ruled a threshold out.
      ENDPOINT      an arc over the node it terminates at. Innocent one at a
                    time, a crowding meter in bulk: a plain cycle has 0 up to
                    n=12, then 2 / 13 / 92 / 201 / 363 at n=14 / 16 / 18 / 20 /
                    24 -- and 14 is where render_cycle starts to go (14 arcs
                    just touch, 16 break).

    Splitting on that line separates them with nothing invented. Across a set of
    realistic test figures, written as a caller would write them, every
    zero-finding figure stays silent, including the two that carry a SAME OBJECT
    exemption, and so do schemes up to 8 text cells; the growing cycle speaks from
    n=14 on. The ratio exemptions/reservations was tried first and is worse than
    useless: linear_text_cell_plus sits at 0.111 and a crowded 18-node cycle at
    0.104, so it ranks the healthy figure above the broken one.

    `_exempt` and `_path` are read directly rather than through accessors.
    allocator.py is frozen by a regression test; two tests already read
    `_exempt` from outside; and re-running related()'s first branch here is
    exact rather than approximate, because SAME OBJECT and ENDPOINT are the only
    two grounds related() has, so "not SAME OBJECT" IS "ENDPOINT".

    Returns the BODY only. The four scheme tools splice it after "LAYOUT REPORT --"
    in a prose receipt, the cycle puts it under a `layout_report` key; both shapes
    are unchanged, and whenever findings > 0 the string is byte-identical to what
    each of the five sites built inline before.
    """
    if alloc is None:
        return ""
    if alloc.findings:
        return ("%d finding(s). The file was written anyway; open it and check "
                "these spots:\n%s" % (len(alloc.findings), alloc.report()))
    from allocator import _path as _owner_path
    n = sum(1 for _a, _b, _d in alloc._exempt
            if _owner_path(_a)[0] != _owner_path(_b)[0])
    if not n:
        return ""
    return ("no findings, but %d place(s) where an arrow overlaps the species it "
            "ends at were passed over as self-collisions. One or two of those is "
            "normal; the count climbs steeply once a figure is crowded, and "
            "nothing else in this receipt tracks that. Open it and look:\n%s"
            % (n, alloc.report()))


def _figure_report(path: str, structures: dict, alloc=None,
                   notes: str = "", placement: dict | None = None,
                   gaps: dict | None = None,
                   footnotes: dict | None = None,
                   size: dict | None = None) -> str:
    """The tool's return value: a small machine-readable report, as JSON text.

    WHY THIS EXISTS. Every field here is known at draw time, and without this it
    is thrown away. `draw_linear_scheme` in particular can drop a cell whose id
    collides with an earlier one; a reply that is only a sentence does not
    mention it -- the figure comes back with the right number of positions,
    every one of them filled, and one of them holding the wrong molecule. A
    caller comparing "structures I sent" against "structures I got" cannot,
    unless the second number is in the reply.

    Once `cells.dropped` exists, silence is structurally impossible: the field
    is always present, so an empty list is a positive statement rather than an
    absence of one.

    WHAT IT DELIBERATELY DOES NOT CONTAIN: the CDXML. A caller cannot act on
    50 kB of markup in a reply, and the file path is the actionable form of it.

    `notes` carries the existing prose -- the same wording as before, unchanged,
    because it says what to DO about each condition. It is empty on a clean
    figure, which is what keeps the common case small.

    THE STRONGEST CLAIM THIS SUPPORTS is "the figure faithfully expresses what
    you sent". It says nothing about whether what you sent was what you meant.
    `ok: true` on a scheme built from the wrong SMILES is a true statement about
    a wrong figure.
    """
    dropped = [i for i, _k, _d in _ID_CLASHES]
    findings = list(getattr(alloc, "findings", None) or [])
    rep: dict = {
        "ok": not dropped and not findings,
        "path": path,
        "cells": {"submitted": _SUBMITTED[0],
                  "rendered": len(structures),
                  "dropped": dropped},
        # ALWAYS PRESENT, for the reason this docstring already gives about
        # `cells.dropped`: "the field is always present, so an empty list is a
        # positive statement rather than an absence of one". A dial that appears
        # only when there is bad news cannot be read for good news.
        "warnings": [],
        # RADICAL CENTRES, counted rather than left to be noticed. A radical is
        # rare in the figures this draws and is usually a SMILES mistake: a
        # bracket atom that lost its implicit hydrogens, or a charge put on the
        # wrong atom. Both happen in practice. `[N+]` with three bonds instead of
        # `[P+]` comes out as a nitrogen radical cation and draws the glyph for
        # one, and short of this count only someone counting SymbolTypes by
        # hand would catch it. A number on the receipt is cheaper than that.
        "radicals": 0,
    }
    # HOW BIG THE DRAWING CAME OUT, which every caller has to guess at otherwise.
    # `alloc.bbox()` is the union of everything the engine reserved -- structures,
    # captions, numbers, footnotes -- and it is the same rectangle finalize() is
    # handed as the content box, so it is the INK, not the page and not the budget.
    # (Checked: with `span` 548 on a three-column table the page is 548 wide and
    # this says 451.8, which is also what measuring the fragments gives.)
    #
    # It matters most for `span`, where the two differ by a known law rather than
    # by a rounding: cells fill the span, but the ink stops where the edge cells'
    # own drawings stop, and each cell centres its content. So
    #     ink = span * (cols - 1) / cols + (leftmost + rightmost content) / 2
    # -- measured at 0.500 / 0.6667 / 0.8333 for 2 / 3 / 6 columns. A caller who
    # needs a particular INK width reads this back and solves; two passes are
    # exact, because the law is linear.
    # `size` for the one tool with no ledger to ask: draw_molecule bypasses the
    # registrar on purpose (one structure, no id space to clash in), so it hands
    # the number in directly. A field that five of six tools carry is its own
    # small trap.
    if size is not None:
        rep["size"] = size
    elif alloc is not None:
        try:
            _bb = alloc.bbox()
            rep["size"] = {"w": round(_bb.w, 1), "h": round(_bb.h, 1)}
        except Exception:
            pass
    _rads = _radical_centres()
    rep["radicals"] = sum(r for _sid, _sym, _c, r in _rads)
    if _rads:
        rep["warnings"].append(
            "%d radical centre(s): %s. Rare on purpose -- check the SMILES if "
            "none was meant (a bracket atom carries no implicit hydrogens, so "
            "[C:2] is a radical and [CH:2] is not)."
            % (rep["radicals"],
               ", ".join("cell %s %s%s" % (sid, sym,
                                           "+" if c > 0 else "-" if c < 0 else "")
                         for sid, sym, c, _r in _rads)))
    # THE FALLBACK RATE, as a number rather than a sentence.
    #
    # The three prose notes below already name WHICH molecules lost the depiction
    # template and why. None of them carried the DENOMINATOR, and that is the
    # difference between a footnote and a disaster: "4 structure(s) did not match"
    # is nothing in a series of forty and total failure in a series of four.
    # Without this field, a SMARTS in `template` on a four-cell table comes back
    # `ok: true`, `warnings: []`, and one sentence in `notes`.
    #
    # WHY THIS IS A FIELD AND NOT A WARNING. Whether a template suits a series is
    # only knowable at run time: the schema cannot check it and preflight cannot
    # either, because it depends on the molecules. So the return value is the only
    # place it can be said, and a caller has to be able to READ it rather than
    # grep prose for a count.
    #
    # WIRED TO `warnings` ABOVE HALF, and the threshold is the one thing here that
    # needed choosing. COMPARED TO WHAT: the point of a template is that the
    # series shares one depiction, so the question is whether the figure still has
    # a shared depiction at all. Past half the cells the template is the exception
    # rather than the rule, and the table no longer does the thing a template is
    # for.
    #
    # MEASURED AGAINST REFERENCE FIGURES, from both sides:
    #   seven reference scope tables are 0 fallback out of 4, 4, 6, 8 and 11
    #     (two of them use template="off" and attempt nothing) -- nothing to warn
    #     about, at any threshold. Note in particular that a two-section table is
    #     11/11: sections work, each re-picking its own reference, so the "aliphatic
    #     on an aryl parent" case does NOT fall back in practice.
    #   a four-cell reference table is 3 of 4 -- toluene as the reference for
    #     anisole, chlorobenzene and fluorobenzene, all three MISSES, because the
    #     template's methyl carbon has to map onto a ring-attached O / Cl / F and
    #     does not. The template genuinely helps one cell in four, and that reply
    #     now carries the warning. Its advice was checked on that same figure:
    #     template="c1ccccc1" takes it to 4 of 4.
    #
    # The BLOCK is unconditional; only the warning has a threshold. The rate is a
    # fact and the judgement is separate from it.
    _att = list(getattr(SF, "TEMPLATE_ATTEMPTED", []))
    if _att:
        _pl = len(getattr(SF, "TEMPLATE_PLACED", []))
        _fb = len(_att) - _pl
        rep["template_match"] = {
            "attempted": len(_att),
            "placed": _pl,
            "fallback": _fb,
        }
        if _fb / len(_att) > 0.5:
            rep["warnings"].append(
                "template/fallback: %d of %d structure(s) (%.0f%%) could not be "
                "drawn on the depiction template and were laid out on their own, "
                "so the figure has no shared orientation to speak of -- which is "
                "the one thing a template is for. Either pick a reference the "
                "series actually shares (any member works, and 'auto' takes the "
                "first), or name just the shared core in 'template_patt', which "
                "matches on that part alone -- and write THAT one's ring bonds "
                "as ~, e.g. [#6]1~[#6]~[#6]~[#6]~[#6]~[#6]1, because it is "
                "matched against the kekulized reference and 'c1ccccc1' can "
                "never hit. The notes say which structures missed and why."
                % (_fb, len(_att), 100.0 * _fb / len(_att)))
    # THE PAGE WAS WIDER THAN THE ONE THAT WAS ASKED FOR.
    #
    # finalize() grows a fixed-width page rather than clipping it -- that part is
    # right, a clipped figure is worse -- but it does so in silence, and 'single'
    # is not a preference, it is a column width a journal enforces. Measured: one
    # sotorasib-sized structure at page='single' came back 370 pt wide, 54% over
    # the 240 pt column, `ok: true`, `warnings: []`, nothing in `notes`.
    #
    # HERE rather than in each engine: all six return through this function, and
    # all six reach finalize, so the alternative is the same six lines six times.
    # ONE WARNING PER AXIS, not one per figure: a table can be both too wide and
    # too tall, and the two have OPPOSITE answers -- narrow the content, or widen
    # the page. Collapsing them would send half of that caller to the wrong fix.
    _WIDE_FIX = {
        # wrap breaks BETWEEN cells, so it cannot help a single over-wide one.
        # Say that here rather than let the caller try it and see no change.
        "single": "page='double' (504 pt); abbreviate=True to shorten the "
                  "widest structures; or wrap=True (it breaks rows only "
                  "BETWEEN cells, so it does nothing when one cell alone is "
                  "too wide)",
        "double": "page='letter', abbreviate=True, or split the figure",
        "letter": "split the figure -- it is wider than the paper",
    }
    # A NARROW PAGE IS WHY IT IS TALL, which is why this advice runs the other
    # way. At page='single' a scope table gets ONE column -- two cells do not fit
    # 180 pt of usable width -- so the products come down in a ribbon. Measured
    # on the same eight: 240 x 1114 pt at 'single', 504 x 463 at 'double'.
    _TALL_FIX = {
        "single": "page='double' -- a wider page fits more columns, so the figure comes "
                  "out much shorter",
        "double": "page='letter', or split the figure across two",
        "letter": "split the figure -- it is longer than the paper",
    }
    for _axis, _pm, _ask, _got in sorted(set(C.PAGE_OVERFLOW)):
        rep["warnings"].append(
            # PAGE against PAGE, margins included on both sides. Saying "the
            # figure came out 435 pt" invites measuring the ink instead, which is
            # 370 pt here, and then the two numbers disagree for no reason.
            "page/over-%s: page=%r is %.0f pt (%.2f in) %s and the page came "
            "out %.0f pt (%.2f in), %.0f%% over. It was GROWN to fit rather "
            "than clipped, so nothing is cut off -- but a journal will rescale "
            "it back, shrinking every caption with it. Try: %s."
            % (_axis, _pm, _ask, _ask / 72.0,
               "wide" if _axis == "width" else "tall",
               _got, _got / 72.0, 100.0 * (_got - _ask) / _ask,
               (_WIDE_FIX if _axis == "width" else _TALL_FIX)[_pm]))
    if dropped:
        rep["dropped_reason"] = {
            i: "duplicate id: kept %s, dropped %s" % (k, d)
            for i, k, d in _ID_CLASHES}
    if findings:
        rep["geometry"] = {
            "collisions": len(findings),
            "detail": findings[0].message[:160],
        }
    # Layer 2: read the file back with RDKit's own CDXML parser and check the
    # molecules survived. Best-effort by construction -- a checker that can fail
    # the whole call is worse than no checker.
    try:
        import roundtrip_inchi as _RT
        if _BUILT_FROM and path and not path.startswith("/dev/"):
            # A FLOATING placeholder is excluded rather than compared. Its bond
            # deliberately stops on a bare node inside the ring instead of the
            # atom the SMILES attached it to -- that IS the feature, it is
            # already reported above, and comparing it would fail every such
            # cell for doing exactly what was asked.
            _flo = {_s for _s, _lbls in getattr(SF, "FLOATING", [])}
            _skip = {cid: "floating placeholder -- the connection table is "
                          "deliberately not the SMILES's"
                     for cid, smi in _BUILT_FROM.items() if smi in _flo}
            ic = _RT.check(path, dict(_BUILT_FROM), skip=_skip)
            rep["inchi_match"] = ic
            if ic["failed"]:
                rep["ok"] = False
                # ...AND ONTO THE DIAL PEOPLE ACTUALLY READ.
                #
                # Without this, the shape of this failure is `ok: false` with
                # `warnings: []`: structural integrity reported through
                # `inchi_match`, layout through `warnings`, and the two never
                # speak. A caller reading `warnings` -- which this module's own
                # docstring calls "the one dial to read" -- sees an empty list and
                # believes the figure is clean while a carbene has quietly become
                # a CH2.
                rep["warnings"].append(
                    "structure/inchi-mismatch: %d of %d structure(s) did not "
                    "survive the round trip through the file: %s. The drawing is "
                    "not the molecule that was sent -- read 'inchi_match'."
                    % (len(ic["failed"]), ic["checked"], ", ".join(ic["failed"])))
    except Exception as exc:                      # noqa: BLE001
        rep["inchi_match"] = {"error": "%s: %s" % (type(exc).__name__, exc)}
    # Layer 3: the descriptors the SOURCE printed, against RDKit's CIP labels.
    # Beside inchi_match rather than inside it because it answers the other half
    # of the question: layer 2 asks whether the drawing is the molecule that was
    # sent, this asks whether the molecule that was sent is the one the paper
    # describes. Silent unless a cell carries `cip`.
    try:
        _cm = _cip_report()
        if _cm is not None:
            rep["cip_match"] = _cm
            if _cm["failed"]:
                rep["ok"] = False
                rep["warnings"].append(
                    "structure/cip-mismatch: %d of %d structure(s) do not carry "
                    "the configuration the source prints: %s. The drawing matches "
                    "the SMILES; the SMILES does not match the paper -- read "
                    "'cip_match'."
                    % (len(_cm["failed"]), _cm["checked"],
                       ", ".join("%s (printed %s, drawn %s)"
                                 % (f["cell"],
                                    "".join("%d%s" % (v, k) for k, v in
                                            sorted(f["wanted"].items())) or "nothing",
                                    "".join("%d%s" % (v, k) for k, v in
                                            sorted(f["got"].items())) or "nothing")
                                 for f in _cm["failed"])))
    except Exception as exc:                      # noqa: BLE001
        rep["cip_match"] = {"error": "%s: %s" % (type(exc).__name__, exc)}
    # WHERE THE ENGINE PUT THINGS, in its own LAYOUT frame. Only an engine that
    # fills it reports it, so an unwired tool stays silent rather than reporting
    # an empty structure and being believed. `notes` stays last: it is prose, and
    # prose belongs at the bottom of a machine-readable reply.
    if placement:
        rep["placement"] = placement
    # WHAT THE GAPS CAME OUT AT, and whether the figure or the caller chose them.
    # Same rule as `placement`: only a tool that fills it reports it, so the five
    # engines that do not solve anything stay silent rather than reporting an empty
    # structure and being believed.
    if gaps:
        rep["gaps"] = gaps
    # WHAT THE FOOTNOTE BLOCK DID, on the same terms: only a tool that fills it
    # reports it. A fold is invisible in the file -- the caller sent two strings and
    # cannot tell whether they came out as two lines or six -- and `over` non-empty is
    # the only way a word too wide to break ever reaches them.
    if footnotes:
        rep["footnotes"] = footnotes
    if notes.strip():
        rep["notes"] = [n.strip() for n in notes.split("\n\n") if n.strip()]
    return json.dumps(rep, ensure_ascii=False, indent=1)

# id -> the SMILES it was built from, for this batch. Not part of `structures`
# on purpose: that dict's shape is a contract with five packers.
_BUILT_FROM: dict[str, str] = {}


def _radical_centres():
    """[(cell id, element, charge, n_electrons)] over every SMILES submitted.

    Read off _BUILT_FROM rather than threaded through the build, because the
    build has three entry paths and this needs to see all of them.
    """
    out = []
    try:
        from rdkit import Chem
    except Exception:                                        # pragma: no cover
        return out
    for sid, smi in sorted(_BUILT_FROM.items()):
        try:
            m = Chem.MolFromSmiles(smi, sanitize=True)
        except Exception:
            m = None
        if m is None:
            continue
        for a in m.GetAtoms():
            n = a.GetNumRadicalElectrons()
            if n:
                out.append((sid, a.GetSymbol(), a.GetFormalCharge(), n))
    return out

# id -> the halo / dot / locant objects that belong BESIDE that structure.
#
# A PARALLEL DICT, for the same reason `_BUILT_FROM` is one, and here the reason
# is measurable: `structures[sid]` is the (fragment, bbox) PAIR that five modules
# destructure -- `frag, bb = S[...]` in layout_grid, `(_f, bb)` in
# layout_snaking, and so on. Widening it to a triple breaks all five at once, and
# a sixth reader would be added by whoever writes the next engine.
#
# Being separate also makes it OPTIONAL, which is what lets the engines be wired
# one at a time: an engine that never reads `recipe["decorations"]` behaves
# exactly as it would without it, so the engines do not all have to change at once.
_DECOR: dict[str, str] = {}

# id -> (smiles, index_map, base_id): what an ANCHOR needs to resolve a SMARTS
# against a cell (cross-structure annotations).
#
# A THIRD dict rather than widening _BUILT_FROM. _BUILT_FROM's value is compared
# for equality to decide whether two cells sharing an id are the same compound
# (`kept != smiles`), so making it a tuple would either break that comparison or
# make it read oddly. The other reason -- do not widen the (fragment, bbox) pair
# five packers destructure -- is why this is not in `structures` either.
_ANCHOR_SRC: dict = {}

# id -> the stereo descriptors the SOURCE FIGURE printed, e.g. "(1R,2R,3R,4S,5S)".
#
# A FOURTH parallel dict, for the same reason as the three above -- and this one
# is an ASSERTION rather than data the drawing needs. Nothing reads it to draw; a
# figure with it and a figure without it come out byte-identical.
#
# WHY IT EXISTS. The round trip answers "is the drawing the molecule I sent". It
# cannot answer "is the molecule I sent the right one", and that is where the
# errors actually are: transcribing published figures by hand produces wrong
# stereocentres, and the round trip catches none of them. Most are caught by
# reading a descriptor the figure had printed all along and comparing it by hand.
# This does that comparison in the tool instead. It does NOT see an enantiomer
# whose R and S counts are equal -- see _cip_report.
_CIP_WANTED: dict[str, str] = {}

# An R/S descriptor with its optional locant: "R", "(S)", "1R", "4aS", "aR".
#
# The lookbehind is load-bearing. Without it "trans" ends in an s, "TBS" in an S,
# and both would be read as descriptors -- so the token must not be preceded by a
# letter or digit, and R/S must be UPPERCASE, which is how papers print them and
# is not how "trans" spells its last letter.
_CIP_TOKEN = re.compile(r"(?<![A-Za-z0-9])([0-9]*[a-z]?)([RS])(?![A-Za-z])")


def _cip_parse(printed: str):
    """(Counter of tetrahedral descriptors, [axial ones]) out of what a paper prints.

    Forgiving on purpose: the caller copies the figure's own string -- parentheses,
    "98.7% e.e." and locants included -- because anything else is a transcription
    step, and transcription steps are what this whole feature is about.

    A locant STARTING WITH A LETTER is axial (aR, aS) and is counted separately:
    RDKit perceives no stereogenic axis, so comparing it would fail every
    atropisomer for being one. "4aS" is not axial -- that a is a ring fusion --
    which is why the test is the FIRST character, not "contains an a".
    """
    tet, axial = [], []
    for loc, d in _CIP_TOKEN.findall(printed or ""):
        (axial if loc[:1].isalpha() else tet).append(d)
    from collections import Counter
    return Counter(tet), axial


def _cip_report():
    """Printed descriptors vs RDKit's CIP labels, per cell. None if nobody asked.

    COMPARED AS A MULTISET, not per position. A paper prints IUPAC numbering and
    RDKit answers per ATOM INDEX; no mapping exists between the two without
    parsing the compound's name, so "three R and two S" is what can be asked
    without inventing one.

    Measured before it was written, over every printed descriptor in a set of
    published figures: 11 of 11 correct structures agreed (no false alarm), and
    6 of the 7 transcription errors actually made were caught.

    THE SEVENTH WAS AN ENANTIOMER, and that is the blind spot stated exactly:
    mirror a molecule whose R and S counts are EQUAL and the multiset does not
    move. Getting the relative configuration wrong is a different thing and is
    usually caught -- the trans of that same ring comes back {S:2} against the
    printed {R:1,S:1}. So this does not see "an enantiomer with balanced counts",
    which is narrower than "stereochemistry" and is pinned by a test.
    """
    if not _CIP_WANTED:
        return None
    try:
        from rdkit import Chem
        from rdkit.Chem import rdCIPLabeler
        from collections import Counter
    except Exception:                                        # pragma: no cover
        return None
    checked = matched = 0
    failed, notes = [], []
    for sid in sorted(_CIP_WANTED):
        printed = _CIP_WANTED[sid]
        smi = _BUILT_FROM.get(sid)
        want, axial = _cip_parse(printed)
        mol = Chem.MolFromSmiles(smi) if smi else None
        if mol is None:
            notes.append("%s: could not read the SMILES back, nothing compared" % sid)
            continue
        rdCIPLabeler.AssignCIPLabels(mol)
        got = Counter(a.GetProp("_CIPCode") for a in mol.GetAtoms()
                      if a.HasProp("_CIPCode"))
        checked += 1
        if want == got:
            matched += 1
        else:
            failed.append({"cell": sid, "printed": printed,
                           "wanted": dict(want), "got": dict(got)})
        if axial:
            # Reported, NOT failed: the count legitimately differs by the axis.
            notes.append(
                "%s: %r also names %d axial descriptor(s) (%s). RDKit perceives no "
                "stereogenic axis, so those are not compared -- SMILES cannot hold "
                "one either, and two atropisomers come out as the same drawing."
                % (sid, printed, len(axial), ", ".join("a" + d for d in axial)))
    rep = {"checked": checked, "matched": matched, "failed": failed}
    if notes:
        rep["notes"] = notes
    return rep


def _new_batch() -> None:
    """Reset the per-batch reports (and the template cache) before converting."""
    SF.clear_template_cache()
    _ID_CLASHES.clear()
    _BUILT_FROM.clear()
    _DECOR.clear()
    _ANCHOR_SRC.clear()
    _CIP_WANTED.clear()
    _GRAPHICS_BOX.clear()
    _GRAPHICS_OUTSIDE.clear()
    _GRAPHICS_COVERS.clear()
    _SLOTS_MOVED.clear()
    _PLACED_MOVED.clear()
    _PALETTE_SIZE.clear()
    DEC.DELOCALIZED.clear()
    _PALETTE_RESCALED.clear()
    C.PAGE_OVERFLOW.clear()
    _SUBMITTED[0] = 0


def _register(structures: dict, smiles: str, spec: dict | None, *,
              prefix: str = "s", auto: str | None = None,
              where: str = "", **kw) -> str:
    """Build one structure, file it under its id, and return that id.

    THE ONE PLACE a figure tool turns a cell into an entry in `structures`. Seven
    separate places would be seven chances to forget -- the same count and the
    same shape as the six fragment builders `_build_fragment` replaced.

    Two defects live in per-tool copies, both invisible on paper:

    **`"id": ""` collapses every cell that carries it.** With the id taken from
    `str(spec.get("id", <auto>))`, dict.get only falls back when the key is
    ABSENT. An empty string is a present value, so every item written `"id": ""`
    -- which the docs invite, since id is optional -- lands on the same key and
    the last one wins. Here an id that is empty or blank means "I am not giving
    one", and gets the auto number.

    **A clash replaces rather than drops, and nothing says so.** Both cells
    still draw; they draw the SAME structure. That is worse than losing one,
    because the cell count is right and every position has something in it.

    Reported only when the SMILES DIFFER. Two cells legitimately share an id when
    they are the same compound -- a struck-out failed target and the actual
    product both labelled 8a, which is a common way to draw a scope table -- and
    there the shared structure is the engine working correctly. Warning on every
    duplicate would put a note on that figure every time until somebody turned
    the notes off, which is how a lint with false positives gets disabled whole.

    `structures[sid]` stays the (fragment, bbox) PAIR every packer unpacks --
    layout_grid does `frag, bb = S[...]` and layout_snaking iterates
    `(_f, bb)` -- so what each id was built from is remembered alongside, in
    `_BUILT_FROM`, rather than by widening a tuple five modules destructure.
    """
    # Counted BEFORE the duplicate check, so `submitted` is what the caller
    # handed over and `len(structures)` is what survived. The gap between them
    # is the number that is otherwise unobservable.
    _SUBMITTED[0] += 1
    sid = str((spec or {}).get("id") or "").strip()
    if not sid:
        # `auto` where the caller has a better name than a running count -- the
        # feeder entries encode their position, which is worth more in a report
        # than "s7".
        sid = auto or "%s%d" % (prefix, len(structures) + 1)
    if sid in structures:
        kept = _BUILT_FROM.get(sid)
        if kept != smiles:
            _ID_CLASHES.append((sid, kept or "?", smiles))
        return sid                      # first one wins; it is already built
    # The sink for halo / dot / locants. Passed on every build, not only
    # where the caller asked for one: a sink costs an empty list, and deciding
    # per call site whether to pass it is how tools end up unable to take a
    # highlight in the first place.
    sink: list = []
    # The index_map is what turns a caller's SMARTS into this fragment's node ids,
    # and it is only knowable at build time -- smiles_to_fragment fills it as it
    # converts. Captured for every structure, not only where an anchor was asked
    # for: deciding per call site which structures are anchorable is how tools
    # end up unable to take a highlight.
    imap: dict = {}
    structures[sid] = _build_fragment(smiles, spec, decorations=sink,
                                      index_map=imap,
                                      where=where or ("cell %r" % sid), **kw)
    if imap:
        _ANCHOR_SRC[sid] = (smiles, imap, 1)
    if sink:
        _DECOR[sid] = "\n".join(sink)
    _BUILT_FROM[sid] = smiles
    if (spec or {}).get("cip"):
        _CIP_WANTED[sid] = str(spec["cip"])
    return sid


_GRAPHICS_BOX: list = []          # [w, h] of the last figure's content box

# [appended] -- how many [r, g, b] colours this figure added to the document's
# colortable, on top of the standard eight. A list rather than an int for the
# same reason _GRAPHICS_BOX is one: _new_batch clears it, so "the last figure
# extended nothing" and "no figure has run" are the same empty state, and no
# note is printed for either.
_PALETTE_SIZE: list = []

# Colours that arrived as 0-255 and were divided down. Recorded rather than
# silently converted: rescaling is a GUESS at what the caller meant, and a guess
# that says nothing is the shape this module keeps paying for elsewhere.
_PALETTE_RESCALED: list = []


_SLOT_KEYS = frozenset({"at", "structures", "frame", "fill", "arrange", "cols",
                        "gap", "pad", "label_size", "label_gap", "label_color",
                        # say WHAT it sits beside instead of where it goes.
                        "of", "side", "align", "clear"})
_SLOT_SIDES = ("right", "left", "above", "below")
_SLOT_ALIGNS = ("start", "center", "end")

# ---- a GRAPHIC that says what it sits beside --------------------------------
#
# WHY THIS EXISTS, in one figure. Hand-computed coordinates put a TfOH badge on
# top of the words "Ar-H" and "HFIP", and a TM badge above its reaction row
# instead of on its arrow -- because the caller guessed `w * 0.36, h * 0.17` from
# the figure's reported size. The ledger knows exactly where that arrow is the
# whole time; without this, nothing can ask it.
#
# The four outward sides are slot()'s, arithmetic and all -- `near_rect` plus the
# `free_near` search, the same one, because "find room beside this" has one
# answer and it should not get two implementations. What is new is the two sides
# that do NOT search:
#
#   "center"  -- centred ON the anchor. A catalyst badge sits on its arrow; that
#               overlap is the motif, not a collision, so there is no search to
#               run and nothing is reserved (reserving it would make the NEXT
#               placed graphic step around a badge that is meant to be there).
#   "around"  -- the anchor grown by `pad`. A frame round a block, a tint behind
#               one. Same reasoning: it encloses on purpose.
_PLACE_SIDES = ("above", "below", "left", "right", "center", "around")
# `box` and NOT `size`: a text graphic's `size` is its FONT size, and with
# "size" in here a placed text has its font size stripped as a placement key,
# is measured at the default 10 pt, and lands 0.9 pt off. Two meanings, one
# word, and the figure still draws.
_PLACE_KEYS = frozenset({"of", "side", "align", "gap", "pad", "box", "name"})

# [(which graphic, direction, distance, how)] for placed graphics the search had
# to push further out. Same discipline as _SLOTS_MOVED and CAPTION_MOVED: a
# silent relocation swaps a reported overlap for an unexplained gap.
_PLACED_MOVED: list = []


def _anchor_rect(of, alloc, fig_rect, where):
    """The LAYOUT-frame rectangle a placed graphic hangs off.

    `of` is "figure", an ORDINAL ("arrow#0" -- the first arrow in reading order),
    or an exact ledger owner ("cell:s2").

    THE ORDINAL FORM IS THE USABLE ONE and it is why this is not just
    `alloc.find()`. A linear scheme's arrow is reserved as `arrow:r0.95001` --
    row 0 and an internal id -- which nobody can predict from the recipe they
    wrote. "arrow#0" they can.
    """
    if of == "figure":
        return fig_rect
    res = list(alloc.reservations())
    # Reading order: down the rows first, then left to right inside a row. The
    # y is rounded because two cells on one row differ by fractions of a point
    # (their own heights centre them differently) and an unrounded sort would
    # interleave them with the next row's.
    res.sort(key=lambda r: (round(r.rect.y, 0), r.rect.x))
    if "#" in of:
        role, _, n = of.partition("#")
        same = [r for r in res if r.owner.split(":", 1)[0] == role]
        try:
            return same[int(n)].rect
        except (ValueError, IndexError):
            raise ValueError(
                "%s: 'of' %r -- this figure has %d %r, so the ordinals are "
                "%s. (Or name an owner exactly: %s.)"
                % (where, of, len(same), role,
                   ", ".join("%s#%d" % (role, k) for k in range(len(same))) or "none",
                   ", ".join(r.owner for r in same[:6]) or "none")) from None
    for r in res:
        if r.owner == of:
            return r.rect
    roles = {}
    for r in res:
        roles.setdefault(r.owner.split(":", 1)[0], 0)
        roles[r.owner.split(":", 1)[0]] += 1
    raise ValueError(
        "%s: 'of' %r is not in this figure's ledger. Use \"figure\", an ordinal "
        "(%s), or an exact owner (%s)."
        % (where, of,
           ", ".join("%s#0..%d" % (k, v - 1) if v > 1 else "%s#0" % k
                     for k, v in sorted(roles.items())) or "none",
           ", ".join(r.owner for r in res[:8]) or "none"))


def _place_graphics(specs, alloc, shift, ext):
    """Rewrite every spec that says `of` into one that says `at`.

    Returns a NEW list; specs that already carry `at` come through untouched and
    byte-identical, so a figure that never asks for placement is unaffected.

    THE RESOLUTION HAPPENS HERE AND NOT IN graphics.py, on purpose. graphics.emit
    takes coordinates and knows nothing about ledgers, pages or siblings -- the
    same line slot() is held to by test_slot_cannot_see_outside_itself. Moving it
    would put the spatial ledger inside the object library, which is the refactor
    the design rules out.
    """
    from allocator import Rect as _R
    if not any(isinstance(s, dict) and _PLACE_KEYS & set(s) for s in specs):
        return specs, frozenset()
    if alloc is None or not shift:
        raise ValueError(
            "graphics: 'of' asks the engine to find room, which needs the figure's "
            "space ledger -- and this tool did not hand one over. Use 'at' with "
            "explicit coordinates here.")
    # THE FIGURE AS IT WAS BEFORE ANY GRAPHIC WAS PLACED, captured once. Reading
    # alloc.bbox() per spec would make "right of the figure" mean "right of
    # whatever is currently rightmost", so the second banner would already be
    # clear of the first and the search would never report a thing -- the same
    # trap _auto_slot_at documents for two ligand boxes.
    fig_rect = alloc.bbox()
    ox, oy = shift[0]
    out, placed = [], set()
    for i, spec in enumerate(specs):
        where = "graphics[%d]" % i
        if not isinstance(spec, dict) or not (_PLACE_KEYS & set(spec)):
            out.append(spec)
            continue
        placed.add(i)
        if "at" in spec:
            raise ValueError(
                "%s: give it EITHER 'at' -- you place it -- OR 'of' plus a "
                "'side', and the engine finds room. Got both." % where)
        of = spec.get("of")
        if of is None:
            raise ValueError("%s: 'side' / 'align' / 'gap' only mean something "
                             "with 'of'. Add it, or use 'at'." % where)
        side = spec.get("side", "below")
        if side not in _PLACE_SIDES:
            raise ValueError("%s: 'side' is one of %s, got %r"
                             % (where, ", ".join(_PLACE_SIDES), side))
        align = spec.get("align", "center")
        if align not in _SLOT_ALIGNS:
            raise ValueError("%s: 'align' is one of %s, got %r"
                             % (where, ", ".join(_SLOT_ALIGNS), align))
        anchor = _anchor_rect(of, alloc, fig_rect, where)
        kind = spec.get("kind")
        clean = {k: v for k, v in spec.items() if k not in _PLACE_KEYS}
        pad = float(spec.get("pad", 6.0))

        # -- how big is it -------------------------------------------------
        if kind in ("rect", "line"):
            if side == "around":
                w, h = anchor.w + 2 * pad, anchor.h + 2 * pad
            else:
                box = spec.get("box")
                if not (isinstance(box, (list, tuple)) and len(box) == 2):
                    raise ValueError(
                        "%s: a placed %r needs 'box' [w, h] -- unlike a text or "
                        "a dot it has no content to measure. (side 'around' is "
                        "the exception: it takes the anchor's size.)"
                        % (where, kind))
                # null ON AN AXIS MEANS "as wide / as tall as the anchor".
                # The full-width band -- a title bar across the top of a block --
                # is the single commonest piece of furniture in TOC graphics, and
                # writing its width by hand is the one number that sends the
                # caller back to reading `size` out of the report.
                w = anchor.w if box[0] is None else float(box[0])
                h = anchor.h if box[1] is None else float(box[1])
            clean["at"] = [0.0, 0.0, w, h]
            bx0 = by0 = 0.0
        else:
            # MEASURED BY BUILDING IT ONCE AT THE ORIGIN, the same trick
            # _auto_slot_at uses on a slot. bx0/by0 matter: a dot's `at` is its
            # CENTRE, so its probe box starts at -r and the conversion below has
            # to take that back out or every dot lands one radius off.
            clean["at"] = [0.0, 0.0]
            _probe, pbox = GX.emit([dict(clean)], origin=(0.0, 0.0))
            w, h = pbox[2] - pbox[0], pbox[3] - pbox[1]
            bx0, by0 = pbox[0], pbox[1]

        # -- where does it go ----------------------------------------------
        if side == "around":
            at_layout = (anchor.x - pad, anchor.y - pad)
        elif side == "center":
            at_layout = (anchor.cx - w / 2.0, anchor.cy - h / 2.0)
        else:
            gap = float(spec.get("gap", 8.0))
            # Step by the graphic's OWN extent, not by `gap`: a step smaller than
            # the thing being moved cannot move it clear. _auto_slot_at's note.
            step = (w if side in ("right", "left") else h) + gap

            def cands(_side=side, _w=w, _h=h, _gap=gap, _step=step, _a=anchor):
                for k in range(4):
                    r = alloc.near_rect(_a, _side, _w, _h, _gap + k * _step, align)
                    yield ((r.x, r.y), r,
                           ("further from %s" % of, k * _step, "%d of 3 steps" % k))

            got, _k = alloc.free_near(cands(), "panel", where, _PLACED_MOVED)
            if got is None:
                # Nowhere clear: keep the conventional spot and let _reconcile
                # report the overlap. The convention is the honest answer when
                # there is no better one -- what the caption placers do.
                r = alloc.near_rect(anchor, side, w, h, gap, align)
                got = (r.x, r.y)
            at_layout = got
            # NAMED, so a second graphic can hang off this one: the white
            # title rides on the band it sits in, and "graphic:banner" says that
            # where "graphic:3" would silently follow the list order.
            alloc.reserve(_R(at_layout[0], at_layout[1], w, h),
                          "graphic:%s" % spec.get("name", i), "panel")

        # -- LAYOUT -> PAGE -> CONTENT, the caller's own frame last ---------
        x = at_layout[0] + ox - ext[0] - bx0
        y = at_layout[1] + oy - ext[1] - by0
        clean["at"] = ([x, y, x + w, y + h] if kind in ("rect", "line")
                       else [x, y])
        out.append(clean)
    return out, placed


def _auto_slot_at(sp, si, box, alloc, shift, ext, taken, fig_rect):
    """The corner an auto-placed slot goes at, in the caller's CONTENT frame.

    WITHOUT THIS THE CALLER HAS TO SAY WHERE. Working that out needs to know where
    the figure ended and what room is left beside it, which is the ledger -- and
    the ledger reaches this function. So this is the same slot, drawn by the same
    `graphics.slot()` with the same signature; the only difference is who computes
    the corner.

    THREE FRAMES MEET, as they do in _reconcile: the ledger is LAYOUT, this function
    answers in CONTENT, and PAGE is the bridge.
    `graphics.slot()` is never told about any of them -- it takes a corner, which is
    what keeps slot() unable to see outside itself, as a regression test asserts.
    """
    from allocator import Rect as _R
    of = str(sp.get("of") or "figure")
    side = sp.get("side", "right")
    align = sp.get("align", "start")
    if side not in _SLOT_SIDES:
        raise ValueError("slots[%d]: side is one of %s, got %r"
                         % (si, ", ".join(_SLOT_SIDES), side))
    if align not in _SLOT_ALIGNS:
        raise ValueError("slots[%d]: align is one of %s, got %r"
                         % (si, ", ".join(_SLOT_ALIGNS), align))
    if alloc is None or not shift:
        raise ValueError(
            "slots[%d]: 'of' asks the engine to find room, which needs the figure's "
            "space ledger -- and this tool did not hand one over. Use 'at' with "
            "explicit coordinates here." % si)
    if of == "figure":
        # THE FIGURE AS IT WAS BEFORE ANY SLOT WAS PLACED, captured once by the
        # caller. Using alloc.bbox() here instead would make "right" mean "right of
        # whatever is currently rightmost", so a second slot's first candidate is
        # already clear of the first and the search never reports a thing -- the
        # outward stepping and its note become unreachable code. Two ligand boxes are
        # both "to the right of the scheme"; what separates them is the search.
        anchor = fig_rect
    else:
        try:
            anchor = alloc.find("block:%s" % of).rect
        except Exception:
            blocks = sorted(r.owner[6:] for r in alloc.reservations()
                            if r.owner.startswith("block:"))
            raise ValueError(
                "slots[%d]: 'of' is \"figure\" or one of this figure's blocks (%s), "
                "got %r. Only a scope table has blocks; on a scheme, a cycle or a "
                "radial figure the whole figure is the only anchor there is."
                % (si, ", ".join(blocks) or "none", of))
    w, h = box[2] - box[0], box[3] - box[1]
    gap = float(sp.get("clear", 26.0))          # == table_gap, one block-sized break

    # STEP BY THE SLOT'S OWN EXTENT along the direction it is being pushed, plus the
    # gap -- the same argument give_way_cap makes for its budget. Stepping by `gap`
    # alone (26 pt) cannot clear a 95-pt-wide ligand box in four tries, so the
    # search always fails and falls back to the conventional spot: two boxes exactly
    # on top of each other. A step smaller than the thing being moved cannot move it.
    step = (w if side in ("right", "left") else h) + gap

    def cands():
        for k in range(4):
            r = alloc.near_rect(anchor, side, w, h, gap + k * step, align)
            yield ((r.x, r.y), r,
                   ("further from the figure", k * step, "%d of 3 steps" % k))
    at_layout, _i = alloc.free_near(cands(), "panel", "slots[%d]" % si, taken)
    if at_layout is None:
        # Nowhere clear: keep the conventional spot and let _reconcile report the
        # overlap. Same choice the caption placers make -- the convention is the
        # honest answer when there is no better one.
        r = alloc.near_rect(anchor, side, w, h, gap, align)
        at_layout = (r.x, r.y)
    # OWNER "slot:0", NOT "slots[0]". allocator._path strips a trailing [N] on
    # purpose -- reserve_polyline names an arrow's chunks "owner[0]", "owner[1]",
    # and one chunk is not a different owner from another. So "slots[0]" and
    # "slots[1]" collapse to the same head and related() exempts them from each
    # other: two slots stacked exactly on top of each other would report NOTHING.
    # `[N]` means "part N of one thing"; the Nth of
    # several things takes the role:name form every other owner here uses.
    alloc.reserve(_R(at_layout[0], at_layout[1], w, h), "slot:%d" % si, "panel")
    ox, oy = shift[0]
    # LAYOUT -> PAGE -> CONTENT, and the caller's own frame is the last of those.
    return [at_layout[0] + ox - ext[0] - box[0],
            at_layout[1] + oy - ext[1] - box[1]]

# Every Z already on an object, so the whole body can be re-stamped in one pass.
# Narrow on purpose: only the attribute, never a Z inside a coordinate string.
_Z_STRIP = re.compile(r'\s+Z="\d+"')


def _slot_items(spec, structures, abbreviate, where):
    """Build one slot's structures, through _register like every other cell.

    Not through _build_fragment directly, and the difference is the point: going
    the same way a product cell goes means a slot inherits the duplicate-id check
    and the empty-id fallback that live there. A ligand box repeating a
    catalyst the scheme already drew is then ONE structure under one id, and a
    slot that names an id the figure used for something else says so.
    """
    out = []
    for k, s in enumerate(spec.get("structures") or []):
        if not isinstance(s, dict) or not s.get("smiles"):
            raise ValueError("%s.structures[%d]: needs a 'smiles'" % (where, k))
        sid = _register(structures, s["smiles"], s, prefix="w",
                        abbreviate=abbreviate, where="%s structure" % where)
        frag, bb = structures[sid]
        # The 4th slot is the decoration sink: a slot structure may carry a halo
        # like any other.
        # graphics.slot() stays ignorant of what a SMILES is -- it receives the
        # decoration already built, exactly as it receives the fragment.
        out.append((frag, bb, s.get("label"), _DECOR.get(sid)))
    if not out:
        raise ValueError("%s: a slot needs at least one structure" % where)
    return out


def _reconcile(box, what, alloc, shift):
    """Report what a caller-placed box lands on. Nothing is moved.

    THE HOLE THIS FILLS. Without it `graphics` and `slots` do no space accounting
    at all -- not loose accounting, none: nothing builds or queries an Allocator
    for them, so the only check is against the page edge. A ligand box dropped
    on top of the scheme, or two slots on top of each other, comes out silently
    wrong, and "silently wrong" is what this whole layer exists to remove.

    THREE COORDINATE FRAMES MEET HERE, so they are named:
      LAYOUT   what the engine reserved. Its origin is the engine's own and may be
               negative (layout_cycle puts the ring centre at 0, 0).
      PAGE     what is in the written file. PAGE = LAYOUT + the shift finalize
               applied, which is why that shift is an out-parameter now: it is
               known in exactly one place and re-deriving it here was the
               alternative this function's docstring already refused.
      CONTENT  what the CALLER writes, measured from the figure's content box.
               Converted to PAGE by the caller of this function (`+ ext`), which
               is the two-pass contract the return value reports.

    `box` arrives in PAGE. Reporting only -- never repositioning -- because the
    caller chose that spot on purpose and the engine has no idea what for.
    """
    if alloc is None or not shift:
        return []
    from allocator import Rect as _R
    ox, oy = shift[0]
    r = _R(box[0] - ox, box[1] - oy, box[2] - box[0], box[3] - box[1])
    # tol=0.0: the strict question. The policy tolerance exists to absorb the
    # engine's own prediction drift, and a caller's box is not a prediction -- if
    # it touches a structure at all, the caller wants to know.
    return [(what, res.owner, res.kind, depth)
            for res, depth in alloc.probe(r, tol=0.0)]


def _add_graphics(doc: str, specs, slots=None, structures=None,
                  abbreviate: bool = False, alloc=None, shift=None) -> str:
    """Draw the caller's boxes / rules / blocks onto a FINISHED document.

    Post-processing, deliberately, and the alternative is what makes the case:
    every tool hands its recipe to a render_* in its own layout module and that
    module finalizes. Threading `graphics` through would mean five signatures and
    five insertion points -- writing the same three lines once per engine, which
    is the shape this codebase keeps paying for rather than a fix for it.

    Working on the finished document also makes the coordinate contract true
    rather than approximate: the content's position is MEASURED off the page it
    actually landed on, so the caller's (0, 0) really is the top-left of the box
    the previous call reported, whatever page_mode did with the margins.

    The page is NOT grown. Coordinates inside the reported box are inside the page
    by construction -- the margin lies outside them -- so only a graphic that
    reaches past the margin can escape, and that one is reported rather than
    silently repositioned. Growing the box would mean re-deriving finalize's
    translate here, in a second place, from a worse vantage point.
    """
    if not specs and not slots:
        return doc
    specs = list(specs or [])
    slots = list(slots or [])
    for i, sp in enumerate(slots):
        if not isinstance(sp, dict):
            raise ValueError("slots[%d]: each slot is a dict" % i)
        unknown = set(sp) - _SLOT_KEYS
        if unknown:
            raise ValueError("slots[%d]: %s is not a slot key. Known: %s"
                             % (i, ", ".join(sorted(unknown)),
                                ", ".join(sorted(_SLOT_KEYS))))

    # Colours first, and post-hoc on purpose. A band names its colour by an index
    # into the document's colortable, which is normally written before any body
    # exists -- so without this an [r, g, b] works in draw_linear_scheme and
    # nowhere else (the error message for that names the two tools that CAN
    # extend the palette and then refuses the caller anyway). Extending the
    # finished document's colortable instead makes rgb work in all six,
    # without a palette hookup per tool.
    def _is_rgb(v):
        return isinstance(v, (list, tuple)) and len(v) == 3
    rgb = [g for g in specs if isinstance(g, dict) and _is_rgb(g.get("color"))]
    # (dict, key) pairs, NOT copies. A slot spends its colours on two different
    # keys, and rewriting a copy would leave the spec holding raw rgb while the
    # palette held the index -- colours that silently come out wrong, which is
    # the failure the whole extend-first order of operations exists to avoid.
    wants = [(g, "color") for g in rgb]
    # A text graphic's own background slab takes rgb on the SAME terms its ink
    # does. Leaving it out fails loudly rather than silently only because
    # color_index refuses a list -- were it coerced, every pale slab in the
    # figure would come out in whatever colour index 0.9 rounds to.
    wants += [(g, "fill") for g in specs
              if isinstance(g, dict) and _is_rgb(g.get("fill"))]
    wants += [(sp, k) for sp in slots for k in ("fill", "label_color")
              if _is_rgb(sp.get(k))]
    if wants:
        uniq = list(dict.fromkeys(tuple(d[k]) for d, k in wants))
        doc, idx = BT.extend_colortable(doc, [list(c) for c in uniq])
        at = dict(zip(uniq, idx))
        for d, k in wants:
            d[k] = at[tuple(d[k])]

    i = doc.rindex("</page>")
    body = doc[:i]
    # From AFTER the <page> opening tag. That tag carries a BoundingBox of its
    # own -- the whole sheet -- and measuring it as content would report a letter
    # page's 612 x 792 as the figure's box, which would send every caller's
    # coordinates to the wrong place.
    inner = body[body.index(">", body.index("<page")) + 1:]
    ext = C._content_extent(inner)
    if ext is None:
        raise ValueError("graphics: this figure has no content to place them "
                         "against -- draw something first")
    # `of` -> `at`, before anything measures or emits, so everything downstream
    # -- the union box, the per-spec _reconcile, the page check -- sees one kind
    # of spec. A spec that already carries `at` comes back unchanged.
    specs, _placed = _place_graphics(specs, alloc, shift, ext)
    xml, box = GX.emit(specs, origin=(ext[0], ext[1]))
    _GRAPHICS_BOX[:] = [ext[2] - ext[0], ext[3] - ext[1]]
    # Per graphic, not per list: "something in `graphics` covers a structure" is
    # not actionable when there are nine of them. emit() returns one union box, so
    # each spec is measured on its own here.
    for _gi, _g in enumerate(specs):
        if not isinstance(_g, dict):
            continue
        if _gi in _placed:
            # PLACED BY THE LEDGER, so it is already reserved IN the ledger, and
            # reconciling it would find its own reservation and report the
            # graphic as landing on itself (every title would report "covers
            # graphic:t1"). Exactly the trap the auto-placed slots have; same
            # answer. Its real collisions come out of the ledger instead: one
            # fact, one channel.
            continue
        _one, _b1 = GX.emit([_g], origin=(ext[0], ext[1]))
        if _b1:
            _GRAPHICS_COVERS.extend(
                _reconcile(_b1, "graphics[%d] (%s)" % (_gi, _g.get("kind")),
                           alloc, shift))

    # Slots. The offset happens HERE, not inside slot(): that signature takes its
    # own corner and its own content and nothing else, so the one function that
    # could grow into "just check what is beside me" -- ledger awareness arriving
    # unannounced -- cannot see anything to check.
    _slot_boxes: list = []
    # Captured BEFORE any slot is placed, so every slot is positioned against the
    # same figure rather than against the last slot. See _auto_slot_at.
    _fig_rect = alloc.bbox() if alloc is not None else None
    if slots:
        if structures is None:
            structures = {}
        for si, sp in enumerate(slots):
            at, of = sp.get("at"), sp.get("of")
            if (at is None) == (of is None):
                raise ValueError(
                    "slots[%d]: give it EITHER 'at' [x, y] -- you place it -- OR "
                    "'of' (\"figure\", or a block of a scope table) with an optional "
                    "'side' / 'align' -- the engine finds room beside that. Got "
                    "%s." % (si, "both" if at is not None else "neither"))
            if at is not None and not isinstance(at, (list, tuple)) or (
                    at is not None and len(at) != 2):
                raise ValueError("slots[%d]: 'at' is [x, y], got %r" % (si, at))
            items = _slot_items(sp, structures, abbreviate, "slots[%d]" % si)
            style = {k: sp[k] for k in
                     ("frame", "fill", "arrange", "cols", "gap", "pad",
                      "label_size", "label_gap", "label_color") if k in sp}
            ids = dict(gid=98700 + si * 40, tid=98800 + si * 40,
                       base=970000 + si * 20000)
            if of is not None:
                # MEASURED BY COMPOSING IT ONCE AT THE ORIGIN, then composed again
                # where it goes. slot() measures by building, so asking it twice is
                # cheaper than splitting it -- and splitting it would change the one
                # signature a regression test guards.
                _probe, pbox = GX.slot([0.0, 0.0], items, **style, **ids)
                at = _auto_slot_at(sp, si, pbox, alloc, shift, ext, _SLOTS_MOVED,
                                   _fig_rect)
            sxml, sbox = GX.slot(
                [float(at[0]) + ext[0], float(at[1]) + ext[1]], items,
                **style, **ids)
            xml = (xml + "\n" + sxml) if xml else sxml
            box = GX.union(box, sbox)
            if of is None:
                # ONLY for a slot the CALLER placed. An auto-placed one was put
                # there BY the ledger and is reserved in it, so reconciling it
                # against the ledger finds its own reservation and reports the slot
                # as landing on itself. Its collisions come out of the ledger
                # instead (declared-collision, in the LAYOUT REPORT): one fact, one
                # channel.
                _GRAPHICS_COVERS.extend(
                    _reconcile(sbox, "slots[%d]" % si, alloc, shift))
            # Slot vs slot, which the ledger cannot answer: a slot is never
            # reserved in it (the ledger belongs to the engine, and by here the
            # engine has finished). Two ligand boxes at overlapping "at"s is a
            # plain mistake and worth catching, so it is checked directly.
            for _sj, _prev in enumerate(_slot_boxes):
                if not (sbox[2] <= _prev[0] or sbox[0] >= _prev[2]
                        or sbox[3] <= _prev[1] or sbox[1] >= _prev[3]):
                    _GRAPHICS_COVERS.append(
                        ("slots[%d]" % si, "slots[%d]" % _sj, "slot",
                         min(sbox[2], _prev[2]) - max(sbox[0], _prev[0])))
            _slot_boxes.append(sbox)

    page = re.search(r'<page\b[^>]*BoundingBox="([\d.\- ]+)"', doc)
    if page and box:
        px0, py0, px1, py1 = (float(v) for v in page.group(1).split())
        if box[0] < px0 or box[1] < py0 or box[2] > px1 or box[3] > py1:
            _GRAPHICS_OUTSIDE.append((list(box), [px0, py0, px1, py1]))

    # BEHIND the figure -- and inserting first is NOT what puts it there.
    #
    # "Nothing here carries a Z, so ChemDraw draws in document order" is false.
    # It is false in this repo's own words: stamp_z's docstring records three
    # coloured discs emitted BEFORE a fragment that still covered it, and
    # shaded_band's says in capitals that emitting earlier is not enough.
    # The bug that assumption produces is a grey slot with a phosphine correctly
    # placed inside it and invisible -- data right, picture wrong, and only the
    # ONE slot with a fill affected, which makes it look like a placement fault.
    #
    # So: strip every Z and re-stamp in document order. Stripping rather than
    # topping up, because stamp_z skips anything already stamped -- a table that
    # had stamped its own bands would keep them and the new graphics would be
    # numbered into the same low range, tying with them. Document order is what
    # every emitter here already means by "behind" (fills are emitted first, then
    # frames, then structures), so one pass over the whole body makes that true.
    #
    # Only on figures that have graphics. Every other figure carries no Z and
    # renders correctly, and stamping those would rewrite every golden.
    cut = body.index(">", body.index("<page")) + 1
    merged = body[:cut] + xml + body[cut:]
    merged = _Z_STRIP.sub("", merged)
    merged, _next = C.stamp_z(merged)
    return merged + doc[i:]


_GRAPHICS_OUTSIDE: list = []

# [(what the caller placed, what it lands on, that thing's kind, depth)] for this
# batch. Filled by _reconcile, printed by _chem_report. A list, reset per batch,
# like every other report in this file.
_GRAPHICS_COVERS: list = []

# [(which slot, direction, distance, how)] for auto-placed slots the engine had to
# push further out than the conventional gap. Same discipline as CAPTION_MOVED: a
# silent relocation swaps a reported overlap for an unexplained gap.
_SLOTS_MOVED: list = []


# Keys a STRUCTURE SPEC may carry in any tool. One list, because a cell that
# means the same thing should not be spelled six ways.
_STRUCT_KEYS = ("rlabels", "floating", "highlight")

# ---- what a draw_linear_scheme item may say -------------------------------
#
# THE CHECK THAT CATCHES A SILENTLY DROPPED KEY. A cell read for six keys with
# the rest dropped on the floor -- "highlight" among them -- returns success and
# draws a figure with no highlight in it -- indistinguishable from a figure
# nobody asked to highlight.
#
# Every name below was taken from the code that READS it, not from the docstring,
# because the opposite mistake is just as bad: rejecting a key that works turns a
# working call into an error. Sources:
#   cell     layout_linear pass 1/2 + rowspan.py (_span*/rowspan) + _build_fragment
#   arrow    layout_linear a.get(...) + arrows.py spec.get(...) + this layer's
#            reagent_* rewriting
#   conn     layout_linear con.get(...)
_LINEAR_CELL_KEYS = frozenset({
    "smiles", "id", "label", "sublabels", "bracket", "under", "num", "rowspan",
    "rlabels", "floating", "highlight", "cip",
    # Text cell: a species WRITTEN rather than drawn -- "H-C(sp3)", "+ 4a", or a bare
    # compound number whose structure is defined in a panel further down. It
    # takes the label keys because it stands in a structure's slot; `size`
    # overrides the label size for the odd one that wants to be bigger.
    "text", "size",
})
_LINEAR_ARROW_KEYS = frozenset({
    "above", "below", "yield", "style", "above_color", "below_color",
    "yield_color", "repeat", "reverse", "arcs", "reagentEq",
    "reagent_smiles", "reagent_id", "reagent_highlight", "reagent_rlabels",
})
# One cross-structure annotation. Strict: a typo'd "colour" here would be dropped
# and the arrow would come out black with nothing said.
_ANNOT_KEYS = frozenset({"from", "to", "arrow", "text", "color", "bulge", "size"})
# A connector INSIDE a product grid. No "kind": the two cells it
# links are neighbours in a grid, so the run is a straight line and a routed elbow
# would have nothing to route around.
_GRID_CONN_KEYS = frozenset({"from", "to", "above", "below", "style"})
# One panel. "items" is a header row's vocabulary, because that is literally
# the layout it borrows -- structures, plus signs, one arrow with captions, text
# cells.
_PANEL_KEYS = frozenset({"name", "items", "title", "frame", "of", "side",
                         "align", "gap"})

# The four cells that would otherwise have NO key check at all. Every one of
# them would silently drop a typo -- the same defect, in four more places. Each
# set is read off what the code actually consumes: the tool's own keys plus the
# three _build_fragment spends (rlabels / floating / highlight), never a guess at
# what "ought" to be allowed.
#
# DOCUMENTED BEFORE ENFORCED, deliberately. A validator that lands before the
# docs turns a working key into an "unknown key" error, which is a false
# positive.
_FRAG_KEYS = frozenset({"rlabels", "floating", "highlight"})
_PRODUCT_KEYS = frozenset({"smiles", "id", "yield", "sublabels", "label", "foot", "cip",
                           "legend"}) | _FRAG_KEYS
_HDR_ITEM_KEYS = frozenset({"smiles", "id", "templated", "plus", "arrow",
                            "text", "size"}) | _FRAG_KEYS
_RADIAL_PRODUCT_KEYS = frozenset({"smiles", "id", "yield", "foot",
                                  # annotate the arrow that POINTS AT this
                                  # product; each is a list of lines
                                  "above", "below",
                                  # "in" | "out"; default "out". Which way THIS
                                  # spoke points -- per slot, because a real figure
                                  # is one arrow in and several out.
                                  "direction"}) | _FRAG_KEYS
_SNAKING_CELL_KEYS = frozenset({"smiles", "id", "align", "p"}) | _FRAG_KEYS
_LINEAR_CONN_KEYS = frozenset({
    "from", "to", "kind", "via", "style", "above", "below",
})

# Near-misses worth a sentence instead of a bare "unknown key". check_keys takes
# these as its `gone` argument -- the mechanism is already there for renamed cycle
# keys, and this is the same problem: a name that reads as obviously right.
#
# "yield" is the one with evidence behind it. Every a.get("yield") in
# layout_linear belongs to the ARROW dict; no cell-level yield is read anywhere,
# so one written on a structure has always been dropped in silence.
_LINEAR_CELL_HINTS = {
    "yield": ("a yield printed under a STRUCTURE goes in 'sublabels' (['82%']); "
              "'yield' belongs to an arrow, where it sits under the below-caption."),
    "color": "a structure is coloured through 'highlight', e.g. "
             "[{'on': <SMARTS>, 'style': 'color', 'color': 'red'}]",
    "yields": "did you mean 'sublabels'?",
}


def _check_linear_item(it: dict, where: str) -> None:
    """Reject an unknown key on one cell / arrow / connector, naming the path.

    Uses layout_cycle's checker so there is ONE message shape in the codebase --
    that function is public for exactly this reason (see its docstring).
    """
    if not isinstance(it, dict):
        raise ValueError(
            "%s: a cell is a dict like {\"smiles\": \"...\"}, {\"arrow\": {...}} "
            "or {\"plus\": true}, got %s" % (where, type(it).__name__))
    if "arrow" in it:
        CY.check_keys(it, frozenset({"arrow"}), where)
        CY.check_keys(it["arrow"] or {}, _LINEAR_ARROW_KEYS, where + ".arrow")
    elif "plus" in it:
        CY.check_keys(it, frozenset({"plus"}), where)
    else:
        CY.check_keys(it, _LINEAR_CELL_KEYS, where, _LINEAR_CELL_HINTS)
        if "smiles" not in it and "text" not in it:
            raise ValueError(
                "%s: needs 'smiles' (a structure), 'text' (a species written out, "
                "e.g. {\"text\": \"H-C(sp3)\"}), 'arrow' (a reaction arrow) or "
                "'plus'. Got keys: %s" % (where, ", ".join(sorted(it)) or "none"))
        if "smiles" in it and "text" in it:
            raise ValueError(
                "%s: 'smiles' and 'text' are mutually exclusive -- a cell is "
                "either a drawn structure or a piece of type, not both. Same rule "
                "a catalytic-cycle node follows." % where)


def _rgb01(triple):
    """[r, g, b] in whichever convention the caller used -> the 0..1 one CDXML wants.

    CDXML stores a component as a fraction of full scale. Every colour in the
    hand-drawn reference files obeys that, ChemDraw's own export included, all
    inside 0.000..1.000. But 0-255 is what a colour picker, a stylesheet or a
    paper's caption hands you, and without this that spelling is formatted
    straight through: `<color r="217" .../>`, out of range, ok:true, no warning,
    and no measurement of what ChemDraw then displays.

    THE TEST IS THE RANGE, NOT THE SPELLING. "the components look like integers"
    is the rule that suggests itself and it is wrong: a pale yellow halo written
    as [1, 1, 0.75] is two thirds integers, which that rule would
    read as (0.004, 0.004, 0.003) -- black. Only a component ABOVE 1 is
    impossible in the 0..1 convention, so only that is evidence.

    The ambiguous region is therefore the cube from (0,0,0) to (1,1,1) read as
    0-255: eight colours, every one of them indistinguishable from black. Nobody
    asks for those, so nothing is lost by reading them as 0..1.

    Out of BOTH conventions -- negative, or past 255 -- is not a convention
    question and is refused, naming the value. Guessing needs two candidates.
    """
    try:
        vals = [float(v) for v in triple]
    except (TypeError, ValueError):
        raise ValueError(
            "highlight colour %r: an [r, g, b] must be three numbers. Give them "
            "as fractions (0..1, what CDXML stores) or as bytes (0..255, what a "
            "colour picker shows) -- either is read." % (list(triple),))
    bad = [v for v in vals if v < 0 or v > 255]
    if bad:
        raise ValueError(
            "highlight colour %r: %s is outside both conventions -- 0..1 (a "
            "fraction of full scale, what CDXML stores) and 0..255 (bytes, what "
            "a colour picker shows). There is no reading of it to guess at."
            % (list(triple), ", ".join("%g" % v for v in bad)))
    if max(vals) <= 1.0:
        return vals
    out = [round(v / 255.0, 4) for v in vals]
    _PALETTE_RESCALED.append((list(triple), out))
    return out


# {#RRGGBB} -- a custom colour inside a text markup.
#
# The inline colour markup `{name}text{}` already exists (cdxml_core.styled_text),
# but it only takes the palette's eight names or a **colortable index** -- and an
# index is document-level, so the caller cannot compute it. Yields, compound
# numbers, section titles and footnotes could then only use those eight saturated
# colours, while the accent colours in published figures are all custom (a blue
# e.r., coloured labels).
#
# The hex values are scanned out here and appended to the colortable **once**,
# together with the spec's [r,g,b], and the markup is rewritten in place to the
# index. The scan walks arbitrarily nested dicts/lists -- because these strings
# are scattered over rows / header.items / arrow captions / legend rows /
# sublabels / footnotes, seven or eight places; a hand-written list of them
# would miss one, and the one missed raises "unknown colour" quietly.
_HEX_MARK = re.compile(r"\{#([0-9a-fA-F]{6})\}")


def _hex_in(objs):
    """Every six-digit hex of a {#RRGGBB} in objs, in order of first appearance."""
    seen, out = set(), []

    def walk(o):
        if isinstance(o, str):
            for m in _HEX_MARK.finditer(o):
                k = m.group(1).lower()
                if k not in seen:
                    seen.add(k)
                    out.append(k)
        elif isinstance(o, dict):
            for v in o.values():
                walk(v)
        elif isinstance(o, (list, tuple)):
            for v in o:
                walk(v)
    walk(objs)
    return out


def _rewrite_hex(o, idx):
    """Rewrite {#RRGGBB} to {index} in place; a string (immutable) comes back new."""
    if isinstance(o, str):
        return _HEX_MARK.sub(lambda m: "{%d}" % idx[m.group(1).lower()], o)
    if isinstance(o, dict):
        for k in list(o):
            o[k] = _rewrite_hex(o[k], idx)
        return o
    if isinstance(o, list):
        for i in range(len(o)):
            o[i] = _rewrite_hex(o[i], idx)
        return o
    return o


class _PaletteIndex(int):
    """A colortable index that REMEMBERS the [r, g, b] it was resolved from.

    _extend_palette rewrites a spec's `color` in place, from a triple to an index.
    Its docstring justifies that with "the specs come from one JSON tool call, so
    there is nothing else holding a reference to them" -- true of one tool call,
    FALSE of a composed figure, where the script is Python and one module-level
    spec dict is naturally shared by several parts:

        GREEN_RING = {"style": "color", "on": ARENE, "color": [0.06, 0.40, 0.29]}
        c.table("a", rows=[{..., "highlight": [GREEN_RING]}])      # -> index 10
        c.linear("c", items=[{..., "highlight": [GREEN_RING, PINK]}])

    By the second call the triple is gone and `color` is the int 10 -- an index
    into part a's colortable. Part c's table is a different document whose index
    10 is a DIFFERENT colour, so the whole molecule comes out pink. Nothing
    complains: the round trip checks connectivity, collisions() checks boxes, and
    an index is a perfectly valid thing for a spec to carry (a caller may name one
    on purpose to reuse a colour inside ONE document -- see the PALETTE note).

    Downstream this IS an int -- color_index, the f-strings, json.dumps all see an
    int and the bytes are unchanged. The only reader that looks further is
    _extend_palette itself, which re-resolves .rgb against the new document.
    """

    # No __slots__: an int subtype cannot have a nonempty one (TypeError at class
    # creation). The attribute lives in a per-instance __dict__, which is fine --
    # there are a handful of these per figure, not thousands.

    def __new__(cls, index, rgb):
        self = int.__new__(cls, index)
        self.rgb = tuple(rgb)
        return self


def _extend_palette(header: str, spec_lists, keys=("color",), texts=None):
    """Append every [r, g, b] a highlight spec asked for to the document palette.

    A named colour is a SATURATED one -- the ACS table holds primaries and nothing
    else -- and a saturated fill reads as a highlighter pen rather than as a tint.
    So a spec may give [r, g, b] instead, which is appended to this document's own
    colortable and referred to by index. Same move a hand-drawn ChemDraw file
    makes: it can carry an orange at index 10 that the ACS eight do not contain.

    Takes ALL the figure's spec lists at once because the colortable belongs to the
    document, not to a cell -- sixteen scope cells asking for the same pale yellow
    must end up at the same index, and each list is rewritten in place of its
    colours with those indices.

    Returns the extended header and REWRITES h["color"] in place in the spec dicts
    it was given. In place because the alternative -- handing back copies -- means
    every call site has to remember to use them, and the failure mode when one
    forgets is an index written against a palette that never got the colour: a
    silently wrong colour, or an out-of-range index. The specs come from one JSON
    tool call, so there is nothing else holding a reference to them.

    A colour asked for twice is appended ONCE and both specs point at that index.

    A 0-255 TRIPLE IS ACCEPTED AND DIVIDED DOWN. See _rgb01.
    """
    # MORE THAN ONE KEY CAN HOLD A COLOUR. A highlight has `color`; a section's
    # pill badge has `color` for the fill AND `text_color` for the type on it,
    # and a black badge without white type is unreadable. The colortable does not
    # care which key a triple came from -- it is one document-level list either
    # way -- so the caller names the keys and the rest is unchanged.
    # isinstance(h, dict) FIRST. This now runs over `boxes` too, and a box spec
    # is the caller's -- `boxes: ["grid"]` (a bare string, which _emit_boxes
    # refuses with a sentence about what to write instead) would otherwise
    # reach here and die on `.get` with an AttributeError. A malformed spec
    # must still get the
    # engine's message, not this function's stack trace.
    # DEDUPED BY OBJECT IDENTITY, and it is not fussiness. The rewrite below is
    # IN PLACE -- h[k] goes from a triple to an int index -- so visiting the same
    # dict twice does tuple(int) and raises TypeError. One dict lands here twice
    # the moment a caller shares it, which is the natural way to write a figure:
    #
    #     NEW = [{"on": "[CX3:1]=[CX3:2]", "style": "color", "color": [200,40,40]}]
    #     def cell(...): return {..., "highlight": NEW + (extra or [])}
    #
    # `NEW + extra` is a new LIST holding the SAME dict, so every cell's highlight
    # list points at one object. A single-key version has the same line and the
    # same hole; it just takes a shared spec to fall in.
    _seen, flat = set(), []
    for lst in spec_lists:
        for h in (lst or []):
            if not isinstance(h, dict):
                continue
            for k in keys:
                _v = h.get(k)
                # A _PaletteIndex is collected TOO: it is a colour this module
                # already resolved once, against a colortable that is not this
                # document's. Re-resolving is the whole fix; skipping it here is
                # what silently carried a stale index into the next part.
                if not (isinstance(_v, _PaletteIndex)
                        or (isinstance(_v, (list, tuple)) and len(_v) == 3)):
                    continue
                if (id(h), k) in _seen:
                    continue
                _seen.add((id(h), k))
                flat.append((h, k))
    for h, k in flat:
        h[k] = (list(h[k].rgb) if isinstance(h[k], _PaletteIndex)
                else _rgb01(h[k]))
    hexes = _hex_in(texts) if texts is not None else []
    if not flat and not hexes:
        return header
    uniq, order = {}, []

    def want(rgb):
        key = tuple(rgb)
        if key not in uniq:
            uniq[key] = None
            order.append(list(key))
        return key

    for h, k in flat:
        want(h[k])
    hex_key = {}
    for hx in hexes:
        hex_key[hx] = want([int(hx[i:i + 2], 16) / 255.0 for i in (0, 2, 4)])
    header, idx = BT.extend_colortable(header, order)
    for key, i in zip(list(uniq), idx):
        uniq[key] = i
    for h, k in flat:
        h[k] = _PaletteIndex(uniq[tuple(h[k])], h[k])
    if hexes:
        _rewrite_hex(texts, {hx: uniq[hex_key[hx]] for hx in hexes})
    # Recorded, not merely returned: the colortable is a DOCUMENT-level resource
    # that a caller building a figure in two passes has no other way to see, and
    # the ceiling below is the only thing anyone has measured about it.
    _PALETTE_SIZE.append(len(order))
    return header


# ----------------------------------------------------------------------
# Appended to every tool that takes a SMILES. ONE copy, six tools.
#
# It is here rather than written out per tool because the feature is INVISIBLE
# from the signature: `smiles` accepts CXSMILES, no new parameter announces it,
# and a caller who does not know reaches for `@`/`@@`, gets a figure that looks
# right, and ships a file whose embedded data says "single enantiomer". The whole
# point of this note is that this is the one defect the picture cannot show you.
# ----------------------------------------------------------------------
_CXSMILES_DOC = """
    RACEMIC / RELATIVE STEREOCHEMISTRY (dr products). Any `smiles` field here
    accepts CXSMILES, so this needs no extra parameter -- append an extension
    block and the enhanced-stereo group is written into the file:

        "C[C@H](O)[C@@H](C)Cl |&1:1,3|"   AND group 1: both centres are as drawn
                                          OR both inverted -- i.e. racemic with
                                          the RELATIVE configuration defined.
                                          This is what a dr product needs.
        "...                  |o1:1,3|"   OR group 1 (either, not a mixture)
        "...                  |a:1|"      Absolute (the default when unstated)

      The numbers after the colon are ATOM INDICES into that SMILES, 0-based.
      Plain `@`/`@@` alone asserts a SINGLE ENANTIOMER, which is what a figure of
      dr products should NOT say; the wedges look identical either way, so this is
      invisible on paper and wrong only in the data.

      NOTHING IS INFERRED. A structure labelled "dr 3:1" in the caption is not
      guessed into a group -- the engine cannot tell dr from er by looking, and
      guessing would state the wrong chemistry while looking right. You have to
      write the group.

      `&` and `o` REQUIRE a group number (`|&1:...|`, not `|&:...|`). `a` takes
      none.
"""


_BOND_STYLE_DOC = """
    BOND APPEARANCE. `highlight` is not only about colour -- five of its styles
    change what a BOND looks like, and like "squiggle" their anchor names a bond,
    so map the two atoms with :1 and :2:

        {"on": "[C:1]-[C:2]",   "style": "bold"}     heavier: a bond just formed
        {"on": "[C:1]-[C:2]",   "style": "hash"}     evenly hashed
        {"on": "[C:1]-[C:2]",   "style": "dash"}     dashed
        {"on": "[C:1]-[C:2]",   "style": "wavy"}     the bond drawn AS a wave --
                                                     a scope table's attachment
                                                     point, where the shared
                                                     fragment was cut away
        {"on": "[#7:1]-[#5:2]", "style": "dative"}   coordination arrow, N->B.
                                                     The anchor's atom ORDER is
                                                     the arrow's direction.

      "wavy" is NOT "squiggle": that one draws a separate wave ACROSS an intact
      bond (the retrosynthetic cut mark). Same English word, different object.

      Restyling a bond that already said something -- a wedge, or a double bond --
      replaces it and loses that chemistry.

      Full detail, including hash-vs-dash (the engine takes no view) and the
      colour options, is in draw_molecule's `highlight`.
"""


_TEXT_CELL_DOC = """
    A SPECIES YOU WRITE RATHER THAN DRAW. Some things in a scheme are chemistry
    without being a structure, and a cell may carry "text" instead of "smiles":

        {"text": "H–C(sp~3~)"}                  the reactant of a C-H
                                                functionalisation -- a bond type,
                                                which has no structure to draw
        {"text": "+", "id": "4a"}               a co-product the figure will not
                                                redraw: the number goes exactly
                                                where every other product's does
        {"text": "2a"}                          a compound defined in a panel
                                                further down, cited by number
        {"text": "≡"}                           a connective between structures

      It stands in a STRUCTURE's slot, so it takes the same furniture: "id" (the
      bold compound number underneath), "label", "sublabels" (dr / er stacked
      below the number). "size" overrides the type size for the odd one that
      wants to be bigger; the default is a structure's weight, not the small
      size a "+" rides at.

      TWO SPELLINGS THAT LOOK ALIKE AND ARE NOT. Copy the one you want:

        a compound cited BY NUMBER, nothing above it
            {"text": "#2a#"}                one bold line, no number row.
                                            #...# is the bold markup; the
                                            number IS the cell.
        a CONNECTIVE with a number under it
            {"text": "+", "id": "4a"}       "+" in the species position, "4a"
                                            bold underneath where every other
                                            compound number sits.

      Both put "2a"/"4a" on the page, at different sizes and different heights,
      and the picture does not say which you asked for -- so pick from the two
      lines above rather than reasoning it out. In particular {"text": "",
      "id": "2a"} is NOT the first one: it draws an empty species with a number
      hanging below the row.

      This works in a scope table's `header` items as well as in a scheme's
      cells, which is where most of them turn up: the reactant of a C-H
      functionalisation, and a substrate the table cites by number.
"""


_SLOT_DOC = """
    SIDE SLOTS -- structures parked BESIDE the reaction (`slots`). Say EITHER where
    it goes or WHAT it goes beside:
      "of": "figure"        the engine finds room beside the whole figure. On a
                            scope table "of" may also name a block -- "scheme",
                            "grid", "table" -- and a panel is placed the same way.
      "side": "right"       right | left | above | below   (default right)
      "align": "start"      start | center | end
      "clear": 26           the gap; a second slot on the same side is pushed out by
                            its own width and the return value says so.
      "at": [x, y]          you place it, in the content-box frame the return value
                            reports.
    LIGANDS. A framed box beside the arrow holding a ligand or catalyst (L8,
    Xantphos, 4CzIPN, a DES panel) is a slot.

        slots=[{"at": [260, 0],                 # top-left corner, content-box frame
                "frame": "RoundEdge",           # omit for NO frame
                "fill": [0.93, 0.93, 0.93],
                "structures": [{"smiles": "...", "id": "L8", "label": "L8"}]}]

      Each structure's name goes UNDERNEATH it, where a product grid puts a
      compound number -- one rule, not two. "label_color" for the green names
      some papers use.

      "arrange": "row" (default) | "column" | "grid" with "cols". So the four
      shapes in the literature are:
        one structure, one frame     one slot, one structure
        two boxes side by side       TWO slots, each with its own "at"
        one frame round three,
          each with its name         one slot, three structures, arrange "row"
        2x2 staggered                one slot, arrange "grid", cols 2

      "frame" omitted draws no frame at all -- which is a real layout (a panel set off by a vertical dashed rule is a frameless slot plus a
      `graphics` line), not a degenerate one.

      PLACING WITH "at": draw once, read the content box out of the GRAPHICS line in
      the return value, then place against it -- the same two-pass workflow the
      `graphics` boxes use, and the same coordinate frame. With "of" the engine
      finds the room itself.
"""


_ANNOTATE_DOC = """
    ANNOTATIONS INSIDE A STRUCTURE (`highlight` style "annotate"). One primitive
    covering what a mechanism or a scope figure draws ON a molecule -- a curve
    with an optional arrowhead and an optional caption. Which of the five you get
    is which keys you write:

        {"on": "[OX1]", "style": "annotate", "text": "α"}
            a Greek letter beside the atom. No offset -> no leader line.
        {"on": "[CH3]", "style": "annotate", "text": "95%D",
         "offset": [26, -18], "arrow": "full"}
            a caption out at the offset with a leader and an arrowhead.
        {"on": "[CX3:1](=O)[CX4:2]", "style": "annotate",
         "text": "Mukaiyama aldol", "offset": [34, 22]}
            a leader off a BOND to a name. Map the two atoms :1 and :2.
        {"on": "[OX1:1]=[CX3:2]", "style": "annotate", "arrow": "full"}
            an electron-pushing arrow between two atoms. They do NOT have to be
            bonded -- a lone pair to an empty orbital is the usual case.
        {"on": "[OX1:1]=[CX3:2]", "style": "annotate", "arrow": "half"}
            the same with ONE barb: the fishhook a single-electron (homolytic)
            arrow draws.
        {"on": "[cH:1]1[cH][cH:2]ccc1", "style": "annotate"}
            a bare arc, no head and no caption -- what a figure draws to mark
            where a minor regioisomer connects.

      "offset" DECIDES THE LEADER BY BEING THERE, not by being big enough: give
      it and a line is drawn out to the caption, leave it off and the caption
      hugs the atom. "bulge" (0.35) sets how far the curve swells sideways,
      "arrow" is "full" | "half" | "none", "color" and "size" as elsewhere.

      Every match is annotated, not just the first -- the same SMARTS re-run over
      a scope table lands on each molecule's own copy.

      BOTH ENDS MUST BE IN ONE STRUCTURE. An arrow from one cell's atom to
      another cell's is not expressible through `highlight`: the anchors resolve
      inside a single molecule. draw_linear_scheme's own `annotations` does cross
      two cells, and still draws that one as an ARC.

    ADDRESS BY THE SUBJECT'S OWN :n -- "maps", and for a mechanism it is the
    better address. Write the map numbers in the cell's SMILES and name them:

        {"smiles": "[Br-:1].[CH:2]1(CCCCC1)[O:3][P+:4]C",
         "highlight": [{"style": "annotate", "arrow": "full",
                        "tail": [1], "head": [2], "head_dir": 3},
                       {"style": "color", "color": "red", "maps": [2, 3]}]}

      * an annotate arrow with "tail"/"head" needs no "on" at all.
      * "maps" is the same address for every style that takes an atom LIST --
        color, fill, halo, dot, locants, lonepair.
      * WHY BOTHER: a SMARTS cannot say WHICH copy -- all four bromines of CBr4
        come back as ONE match, so `nth` cannot reach the near one. A map number
        is unique.
      * A BRACKET ATOM HAS NO IMPLICIT HYDROGENS and writing :n forces brackets,
        so `[C:2]` in a ring is a RADICAL. Write `[CH:2]`.

      "tail_dir" / "head_dir": one map number each, meaning "this end's tangent
      points at that atom". This is the Sn2 case -- `head_dir` naming the leaving
      group makes the arrow arrive anti to it, i.e. on the backside. Without
      them, each end's direction is chosen automatically.

    LONE PAIRS -- style "lonepair", with "n" for how many (default 1):

        {"on": "[OX2]", "style": "lonepair", "n": 2}
        {"style": "lonepair", "maps": [3], "n": 1}

      Placed automatically in the widest gaps between the atom's bonds.

      DRAW ONLY THE PAIR THAT REACTS: a neutral nucleophile attacks
      through its lone pair, so draw it; an anion attacks as a charge and the
      circled minus already says so.

      RADICAL DOTS ARE AUTOMATIC -- every atom the SMILES gives radical
      electrons to gets one, and there is no style for it. An atom carrying BOTH a charge and
      one unpaired electron gets a single RadicalCation / RadicalAnion glyph
      instead of a dot and a charge separately.
"""


_INLINE_COLOUR_DOC = """
    COLOUR PART OF A LINE. Every piece of type this engine writes -- a section
    title, a table cell, an arrow caption, a sublabel, a footnote -- takes the
    same inline markers, and one of them sets colour for a SPAN:

        "{red}94{}:6 er"        the figure red, the rest of the line as it was
        "see the {blue}SI{}"    one word in a footnote
        "{10}grey{} rest"       a raw palette index, for a colour you added

      **{name} OPENS AND A BARE {} CLOSES.** It is the only asymmetric marker
      here -- the other five need no argument -- so `{red}...{red}` is wrong and
      will render literally.

      A whole line still takes a colour the usual way (above_color,
      below_color, ...); the marker is for spans inside it.

      MARKERS DO NOT NEST, colour included. A green line with bold words in it is
      a block colour plus #bold#, which works; "{green}#94#{}" gives you a green
      run reading "#94#".

      Write {{ and }} for a literal brace.

      The others, unchanged: ~x~ subscript, ^x^ superscript, *x* italic,
      #x# bold, _x_ underline.
"""


_CELL_LEGEND_DOC = """
    A VALUE TABLE UNDER A CELL (`legend` on a product cell). A Markush structure
    is half a figure; the other half says what R stood for and what each gave.

        {"smiles": "[*:1]c1ccccc1", "rlabels": {"1": "R"},
         "legend": {"rows": [["3a", "R = H",   "72%"],
                             ["3b", "R = OMe", "81%"],
                             ["3c", "R = Cl",  "65%"]]}}

      A row is an ORDERED LIST OF CELLS and the order is yours; columns line up
      across rows. Same object draw_molecule's `legend` takes.

      OMIT `id` WHEN THE NUMBERS ARE IN THE TABLE. That is the usual shape -- one
      drawn structure, its variants listed underneath by number -- and there is no
      switch for it: leave `id` out and the cell has no number of its own, so the
      table is the only numbering.

      A LEGEND CELL MAY BE A STRING OR A DICT:

          "R = OMe"                                a plain cell
          {"swatch": "red", "text": "this work"}   a filled dot, then the words

      The swatch is the colour key a figure puts beside structures tinted to
      match. **Its dot is always LEFT-ALIGNED, whatever the column's `align`
      says**; the words beside it follow the dot.

      IT WIDENS THE WHOLE GRID. Every column is as wide as the widest thing in
      any cell: the neighbouring cells stretch with the wide one.
"""


def _cx_doc(fn):
    """Append the shared CXSMILES, bond-style, text-cell, slot, annotate,
    inline-colour and cell-legend notes.

    Before FastMCP reads the docstring.

    Order matters: this must sit BELOW @mcp.tool() in the source so it runs
    FIRST. FastMCP snapshots __doc__ into the tool description at decoration
    time, so appending afterwards would change nothing anyone can see -- which is
    the same class of silent no-op this file guards against elsewhere.
    A regression test asserts the text reaches list_tools().
    """
    fn.__doc__ = ((fn.__doc__ or "") + _CXSMILES_DOC + _BOND_STYLE_DOC
                  + _TEXT_CELL_DOC + _SLOT_DOC + _ANNOTATE_DOC
                  + _INLINE_COLOUR_DOC + _CELL_LEGEND_DOC)
    return fn


def _build_fragment(smiles: str, spec: dict | None = None, *,
                    abbreviate: bool = False, template=None,
                    template_patt: str | None = None,
                    where: str | None = None, decorations: list | None = None,
                    index_map: dict | None = None):
    """(fragment, bbox) for ONE structure, however the tool spells its cell.

    Every figure tool calling smiles_to_fragment itself, with a slightly
    different argument list each time, is exactly how one tool ends up with
    highlight support and another drops the key on the floor. Six call sites,
    six chances to forget one. This is the one place
    that decides what a structure spec may say, so the tools cannot drift apart
    again: adding a key here adds it everywhere at once.

    `spec` is the caller's cell dict; the decoration/markush keys are read off it
    by name. Pass None for a structure that has no dict of its own (the radial
    centre, a bare reagent).
    """
    spec = spec or {}
    return smiles_to_fragment(
        smiles, abbreviate=abbreviate, template=template,
        template_patt=template_patt,
        rlabels=spec.get("rlabels"), floating=spec.get("floating"),
        highlight=spec.get("highlight"),
        where=where, decorations=decorations, index_map=index_map)


@mcp.tool()
@_cx_doc
def draw_molecule(smiles: str, filename: str = "molecule", clean_up: bool = False,
                  page: str = "letter", abbreviate: bool = False,
                  rlabels: dict = None, legend: dict = None,
                  floating: list | bool = None, highlight: list = None) -> str:
    """Draw a single chemical structure from a SMILES string as a ChemDraw (CDXML) file.
    Use UPPERCASE Kekule SMILES for aromatics (e.g. C1=CC=CC=C1).
    page: 'letter' (default) | 'single' | 'double' | 'fit' | 'auto'.

    MARKUSH (generic) STRUCTURES. Write the placeholder as a dummy atom in the
    SMILES and it is drawn as a labelled node: a bare * is "R", and [*:1] [*:2] are
    "R1" "R2" with the digit superscripted. So the usual one-R figure needs nothing
    but the SMILES. Placeholders work anywhere an atom does -- hanging off a ring,
    in the middle of a chain, or AS a ring vertex -- and ChemDraw trims every bond
    that meets one, however many there are.
      rlabels: {"1": "Ar", "2": "X"} renames them, keyed by the map number in
      [*:n]. A value in brackets whose symbols are all real elements ("[N]",
      "[Pd]", "[N,O]") becomes a ChemDraw element-list query instead of a generic
      label. Use ^ ^ for a superscript ("R^1^") and ~ ~ for a subscript; an
      unmarked label goes through ChemDraw's auto-formula mode, which subscripts
      trailing digits (so "CF3" gives CF-sub-3, and "R1" would give R-sub-1).
    The return value says when a figure contains placeholders: it is then a
    drawing, not a molecule, and nothing about it can be searched.

    floating: [1, 2] -- map numbers to draw on a FLOATING bond, the one that
      crosses a ring edge and stops inside because the paper is not saying which
      position is substituted. Use it when the legend values carry their own
      positions ("6-Et", "7-Me", "6,7-OMe") or the numbering does (3b-o / 3b-m /
      3b-p); a bare list of groups ("Et", "OMe") means the drawn position is the
      real one, so leave it off. Per placeholder, not per structure -- one
      molecule routinely mixes a floating R with a fixed one. The placeholder
      still attaches to a ring atom in the SMILES; that attachment stops meaning
      "here" and starts meaning "in from this side". Pass floating=true for all
      of them.

    highlight: say something about a named part of the structure. TWO KINDS of
      thing live here, and the name only suggests the first:
        DECORATIVE  "color" "fill" "halo" "dot" "locants" -- laid over or behind
                    the structure. Remove them and the molecule is unchanged.
        STRUCTURAL  "squiggle" and the five BOND STYLES below (bold / hash / dash
                    / wavy / dative). These change what the drawing SAYS: a wavy
                    bond marks an attachment point, a dative bond is a
                    coordination arrow, a squiggle is where the bond breaks.
      So reach for `highlight` when you want a bond drawn differently, not only
      when you want something coloured.

        [{"on": "c1ccc2c(c1)ccc(=O)o2", "style": "halo",  "color": "green"},
         {"on": "[#8]=[#6]",            "style": "dot",   "color": "red"},
         {"on": "c1ccccc1", "style": "locants", "labels": ["5","6","7","8"]},
         {"on": "[#7]",                 "style": "color", "color": "blue"}]
      "on" is a SMARTS, not a list of atom numbers and not coordinates -- so the
      same spec re-runs over a whole series and lands on each molecule's own copy
      of the group. Every match is decorated, not just the first.

      WHEN YOU MEAN ONE OF SEVERAL IDENTICAL GROUPS, "on" cannot say it and no
      SMARTS can. The two CO2Et on a malonate carbon are ENANTIOTOPIC -- replace
      one and you get R, the other and you get S -- so no constitutional pattern
      separates them, which is precisely why a lipase desymmetrising them is
      worth a figure. Colour them with "on" and BOTH go green.
      Address one with "maps": write the number into the cell's own SMILES.
          "smiles": "CC(C)CC(C#N)C(C(=O)OCC)[C:1](=[O:2])[O:3][CH2:4][CH3:5]"
          "highlight": [{"maps": [1,2,3,4,5], "style": "color", "color": "green"}]
      A map number is unique by construction. "maps" is accepted by every style
      that takes an atom list -- color, fill, halo, dot, locants, lonepair -- and
      the mechanism notes below document it in full.
      A MAPPED ATOM IS A BRACKET ATOM, so it carries no implicit hydrogens: write
      [CH2:4], not [C:4]. A bare [C:4] in a chain is a carbon with three unpaired
      electrons and is refused by name.
      style: "color" (default -- tints the atoms and the bonds between them)
             A WHOLE SKELETON IN ONE COLOUR is `"on": "*"`, and it is the base
             style of most figures that use colour at all. Put the reactive site
             in a SECOND spec after it: the later spec wins on the overlap. All
             of a pattern's matches count as ONE region, so a single-atom pattern
             tints the bonds between them too -- `"[#6]"` gives the carbon
             skeleton and no label, a carbon having none. A PLACEHOLDER's label
             takes colour like any other, by `"[#0]"` or by its map.
             "fill"    ChemDraw's OWN highlight: the ring is filled with colour,
                       the way Color > Highlight does it by hand. Native, so it
                       survives editing, moves with the structure and stays
                       selectable -- prefer it over "halo" for a ring.
                       NEEDS A RING: a chain or a single atom falls back to
                       "color" and the return value says which ones did.
             "halo"    a soft disc behind the group
             "dot"     a filled disc on each matched atom
             "locants" position numbers outside a ring; "labels" go on the
                       matched atoms IN SMARTS ORDER, so the query fixes the
                       numbering. Give null or "" to skip a position.
             "squiggle" the bond-cleavage mark: a wavy line drawn ACROSS one
                       bond, for a retrosynthesis or mechanism figure. This is
                       the ONE style whose anchor names a BOND, so say which
                       two atoms with the map numbers :1 and :2 --
                       {"on": "[CX3:1](=O)[CX3:2]", "style": "squiggle"}.
                       Without them the first two atoms IN SMARTS ORDER are
                       taken, and for "[CX3](=O)[CX3]" those are the C and the
                       O, so the mark lands on the C=O and not on the C-C you
                       meant. Optional "amp" (how far it swings, 0.25),
                       "wavelength" (0.20) and "overhang" (how far past each
                       end of the bond, 0.30) are FRACTIONS OF A BOND LENGTH.
                       Unlike halo / dot / locants it is drawn inside the
                       structure, so it works in every figure tool and moves
                       with the cell it belongs to.
                       NOT "wavy" -- see below.

      BOND STYLES. `highlight` is not only about colour: these five change
      what a BOND LOOKS LIKE. They take an anchor and no colour, and like
      "squiggle" the anchor names a BOND, so map the two atoms with :1 and :2.
             "bold"    a heavier line -- the new bond formed in this step
             "hash"    evenly hashed
             "dash"    dashed. Hash and dash are BOTH flat, non-wedge bonds and
                       the engine assigns no meaning to either. (The tapered stereo hash is a
                       different thing again and is drawn automatically from the
                       SMILES -- you never ask for it.)
             "wavy"    THE BOND ITSELF is drawn as a wave: the attachment-point
                       mark a scope table puts where the shared fragment was cut
                       away. **This is not "squiggle"**, which draws a separate
                       wave ACROSS an intact bond as a cleavage mark. Same word
                       in English, different object and different meaning.
             "dative"  a coordination arrow, N->B, L->M. DIRECTION IS THE ANCHOR'S
                       ORDER: {"on": "[#7:1]-[#5:2]", "style": "dative"} draws
                       N->B and swapping the maps draws B->N.
      A bond can carry only ONE Display and ONE Order, so restyling a bond that
      already said something -- a wedge, or a double bond -- REPLACES it, and the
      return value says so per bond. That is a real loss of chemistry (a wedge is
      a configuration; Order 2 is a double bond) and it is reported rather than
      refused, because cutting a shared fragment away at a stereocentre is a
      legitimate reason to do it.

      Colours: black red yellow green cyan blue magenta, a colortable index, or
      [r, g, b] with each 0..1, which is appended to the document's palette.
      USE [r, g, b] FOR A HALO. The named colours are the ACS primaries and a
      saturated glow reads as a highlighter pen over the structure rather than
      as light behind it -- [1, 1, 0.75] is a pale yellow, [0.85, 0.95, 1] a
      pale blue.
      WATCH THE AROMATICITY. The engine perceives more rings as aromatic than a pen
      does -- coumarin comes back with its lactone ring aromatic, carbonyl and
      all -- so "O=C-O" matches nothing there while "[#8]=[#6][#8]" matches. When
      an anchor hits nothing the return value says so, and says this if it looks
      like the reason.

    legend: the table of values under the structure -- the other half of a Markush
      figure. Each row is an ORDERED LIST OF CELLS and the order is yours.
      Columns line up across rows.
        {"rows": [["#23#,", "R^1^ = OMe,", "72%,", ">19:1 *dr*"],
                  ["#24#,", "R^1^ = Me,",  "68%,", ">19:1 *dr*"],
                  "from Benzocaine"],           # a bare string is a note line
         "header"?: [null, "R^1^", null, null], # a heading over one column only
         "align"?:  ["left", "left", "right", "left"],
         "cols"?: 2}                            # two blocks side by side
      Same text markers as everywhere else: #bold# for the compound number,
      ^sup^ ~sub~ *italic*.
    abbreviate=True collapses recognised pendant groups (Boc, Ts, OTBS, NO2, ...) to
    contracted ChemDraw superatom labels (expand/contract in the app). If nothing in
    the table matches, the structure is drawn in full and the return value says so.

    Molecules SMILES cannot express (ferrocene: the metal is not BONDED to the ring
    carbons) are drawn from a prepared depiction instead -- write the SMILES any
    reasonable way and the return value reports the swap."""
    # REFUSE A BAD SPEC BEFORE ANYTHING IS BUILT. First statement in the
    # tool, ahead of _new_batch, so a rejected call leaves no per-batch state
    # behind either. See preflight.py for what it checks and why the engines'
    # own key whitelists are not enough.
    PF.preflight_molecule({"smiles": smiles, "page": page, "rlabels": rlabels,
                           "highlight": highlight, "legend": legend})
    _new_batch()
    # A named colour is a SATURATED one -- the ACS table holds primaries and
    # nothing else -- and a saturated glow is a highlighter pen, not a glow. So a
    # spec may give [r, g, b] instead, which is appended to the document's own
    # colortable and used by index. Same move a hand-drawn ChemDraw file makes:
    # it can carry an orange at index 10 that the ACS eight do not contain.
    header = _extend_palette(C.acs_header_prefix(), [highlight])
    # draw_molecule is the ONE tool that never repositions its structure, so it is
    # the one that can take halo / dot / locants -- it hands in a list to receive
    # them. Everywhere else that list is absent and those styles raise; see
    # smiles_to_fragment.
    deco_parts: list = []
    frag, fb = _build_fragment(smiles, {"rlabels": rlabels, "floating": floating,
                                        "highlight": highlight},
                               abbreviate=abbreviate, decorations=deco_parts)
    deco = "\n".join(deco_parts)
    parts = [deco, frag]
    if legend and legend.get("rows"):
        lx, lb = LG.compose_legend(legend, gid=1)
        fw, fh = fb[2] - fb[0], fb[3] - fb[1]
        w = max(fw, lb[2])
        dx = (w - fw) / 2.0
        # Both centred on the wider one. base 100/200 keeps the two id ranges
        # apart the same way the packers do it; the decorations already live at
        # 7000+ so they need no base of their own.
        #
        # The decorations move by the SAME dx as the structure, not by place():
        # place() shifts a piece by ITS OWN bbox, which would slide a halo off
        # the group it belongs to. The fragment is normalised to (0,0), so dx is
        # the whole of the shift.
        parts = [C.translate_objects(deco, dx, 0.0) if deco else "",
                 LG.place_legend(lx, 200, (w - lb[2]) / 2.0, fh + LEGEND_GAP),
                 C.place(frag, 100, dx, 0.0)]
    body = "\n".join(p for p in parts if p)
    if deco:
        # THE STRUCTURE IS THE FRONT LAYER. Written last and stamped last, so it
        # takes the highest Z -- and Z, not file order, is what ChemDraw draws
        # by. Only stamped when there ARE decorations: every other figure here
        # carries no Z and renders correctly, and adding one would rewrite every
        # golden to say something it never needed to say.
        body, _next = C.stamp_z(body)
    doc = C.document(header, body, page_mode=page)
    # draw_molecule is the one tool that legitimately bypasses _register (one
    # structure, no id space to clash in), so nothing recorded what it was
    # built from and the round-trip had nothing to compare against. Record it
    # here rather than routing the tool through the registrar, which would be a
    # behaviour change for a reporting fix.
    _BUILT_FROM.setdefault("1", smiles)
    # The union of everything drawn -- a halo reaches further than the atoms do,
    # so the fragment's own box is not the answer.
    # NORMALISED, because not every BoundingBox is written min-then-max: a disc
    # writes centre-then-edge, so a halo's box has x0 > x2 and a plain min()/max()
    # over the raw numbers drops it silently -- and the halo'd molecule then
    # reports the same width as the bare one.
    _bxs = [[float(v) for v in m.split()]
            for m in re.findall(r'BoundingBox="([\d.\- ]+)"', body)]
    _bxs = [[min(b[0], b[2]), min(b[1], b[3]), max(b[0], b[2]), max(b[1], b[3])]
            for b in _bxs if len(b) == 4]
    return _figure_report(_write_and_open(doc, filename, clean_up),
                          {"1": None}, notes=_chem_report(),
                          size={"w": round(max(b[2] for b in _bxs)
                                           - min(b[0] for b in _bxs), 1),
                                "h": round(max(b[3] for b in _bxs)
                                           - min(b[1] for b in _bxs), 1)}
                          if _bxs else None)


@mcp.tool()
@_cx_doc
def draw_linear_scheme(items: list, filename: str = "scheme", clean_up: bool = False,
                       page: str = "letter", connectors: list = None,
                       abbreviate: bool = False, align: str = "left",
                       wrap: bool = False, wrap_style: str = "straight",
                       max_arrow: float = None,
                       graphics: list = None, slots: list = None,
                       annotations: list = None,
                       # ONE NOTE MAY BE WRITTEN BARE. The body takes a string
                       # -- `[footnotes] if isinstance(footnotes, str)` -- and
                       # the tests exercise it, but a plain `list` annotation
                       # would refuse it at the pydantic boundary, so the shape
                       # would work in every direct call and in no real one.
                       # draw_radial_figure has the same branch and
                       # the same widening; draw_scope_table does NOT, and must
                       # not get it -- see the note at its own `footnotes`.
                       footnotes: list | str | None = None,
                       reserve: list = None, span: list = None,
                       template: str = "off", template_patt: str = None,
                       footnote_width: float | None = None) -> str:
    """Draw a LINEAR reaction scheme (one row OR several rows) as a ChemDraw (CDXML) file.

    items: a single row (list of cells) OR a list of rows (each a list of cells) for a
    MULTI-ROW figure in one call. A cell is:
      {"smiles": "...", "id": "4", "label"?: "2HCl", "sublabels"?: ["dr 3:1"], "bracket"?: true,
   "cip"?: "(1R,2R,3R,4S,5S)",   # THE DESCRIPTOR THE SOURCE FIGURE PRINTS, copied
                          # verbatim -- parentheses, "98.7% e.e." and locants and
                          # all. Checked against the CIP labels of your SMILES and
                          # reported as `cip_match`; a mismatch sets ok=false,
                          # because the drawing then matches your SMILES while
                          # your SMILES does not match the paper. THIS IS THE ONLY
                          # CHECK THAT LOOKS AT WHETHER THE INPUT WAS RIGHT.
                          # It does not see an enantiomer whose R and S counts are
                          # equal; it does see a wrong relative configuration,
                          # which moves the counts. Axial descriptors (aR, aS) are
                          # reported and not compared -- SMILES cannot hold a
                          # stereogenic axis.
       "highlight"?: [...],   # THE WHOLE spec draw_molecule takes, all of it: "fill"
                              # (native ChemDraw highlight), "color", "squiggle" (the
                              # bond-cleavage mark), the five BOND STYLES (bold / hash
                              # / dash / wavy / dative), AND "halo" / "dot" /
                              # "locants". The last three are
                              # separate objects beside the structure rather than
                              # inside it, and travel with the cell.
                              # "wavy" makes the BOND a wave (a scope table's attachment
                              # point) -- that is not "squiggle", which draws a wave
                              # ACROSS an intact bond.
       "rlabels"?: {...}, "floating"?: [...],
       "under"?: "30",   # x-align this structure under one in an EARLIER row, which is
                         # what makes a connectors arrow between them come out VERTICAL
       "rowspan"?: 3}    # ONE STRUCTURE FACING SEVERAL ROWS -- a fork (or a merge).
                         # Put it at slot 0 of its row and the rows BELOW it become
                         # its branches: each carries its own arrow and product, and
                         # a spine plus a stub is drawn joining this structure to all
                         # of their arrows. Put it anywhere else and the same key
                         # means a MERGE: the covered rows END at it, and their last
                         # arrows are stretched so they all arrive together.
                         #   fork    [{"smiles": ..., "id": "1", "rowspan": 3},
                         #             {"arrow": {...}}, {"smiles": ...}],   <- branch 1
                         #            [{"arrow": {...}}, {"smiles": ...}],   <- branch 2
                         #            [{"arrow": {...}}, {"smiles": ...}]    <- branch 3
                         # The value counts the OWNING row too, so 3 = 3 branches.
                         # NOT equalized: each branch arrow is still sized by its own
                         # caption, so unequal captions give a ragged product column.
      {"arrow": {"above"?:[...], "below"?:[...], "yield"?:"82%", "style"?:{"kind":"equilibrium"},
                 "above_color"?:"green", "below_color"?:["green",null], "yield_color"?:"blue",
                 "repeat"?:{"n":2,"dx":10,"dy":10}, "reverse"?:true,
                 "arcs"?:[{"at":0.6,"off":[-50,-48]}],
                 "reagent_smiles"?:"...", "reagent_id"?:"29",
                 "reagent_highlight"?:[...], "reagent_rlabels"?:{...}}}
      {"plus": true}

    MULTI-ROW / spanning a step across rows — use ordinary STRAIGHT arrows, never a routed
    turn. Split one step at the row edge in one of two ways:
      • trailing: end a row with the step's {"arrow": ...} cell; its product structure then
        starts the next row (arrow points off toward the margin).
      • leading:  end a row on a structure; start the next row with the step's {"arrow": ...}
        cell, which leads from the left margin into the first structure of that row.
    Both just place the arrow as the last / first cell of a row — no special syntax.

    align: 'left' (default — rows share a left margin and read top-to-bottom like text) |
      'center' (rows centred on a common axis).
    wrap: if True AND page is 'letter'/'single'/'double', break every row that is wider
      than the column into rows that fit, BETWEEN cells; a "+" never starts a row. With
      several rows in `items`, an arrow that ENDS a row (its product starts the next
      one) is never left on a row by itself. The pieces take the row's place, in order;
      a row that fits is left exactly as written, so your own rows survive wherever
      they fit. A rowspan limits it: a row a span COVERS, and a row that owns a FORK, is
      never broken; a row that owns a MERGE breaks only before the merge cell, which
      stays in the last piece -- never at its slot 0 -- with the covered rows right
      after it. The reply says which rows were broken and names every row still wider
      than the column, and why. No effect on 'fit'/'auto' pages.
    wrap_style: 'straight' (default — every wrapped row reads left-to-right, joined at the
      margin) | 'snake' (rows alternate direction and a 90 degree turn carries the route
      down, so the last product of one row sits directly above the first substrate of the
      next). The snake breaks a single flat row only: with several rows in `items`
      nothing is broken. Only meaningful with wrap=True and only
      once the row actually breaks. Every row but the last is JUSTIFIED to the widest one
      so the turns drop plumb, which lengthens arrows on the shorter rows; max_arrow caps
      that, and a row that cannot reach keeps its turn as an elbow.
    max_arrow: optional cap (pts) on arrow length; an over-wide reagent-above-arrow is
      shrunk to fit instead of stretching the arrow. Omit to size every arrow to its caption.
    page: 'letter' (full page, default) | 'double' | 'single' (ACS column widths) | 'fit' | 'auto'.
    Arrow style.kind: normal|retro|equilibrium|resonance|harpoon|fishhook|noGoCross|noGoHash;
    style.line: dashed|bold. ('wavy' is accepted but ChemDraw ignores it on arrows;
      use lineWidth for a thicker shaft than bold.) Text markers: ~sub~ ^sup^ *italic* #bold# _underline_.
    arrow.repeat: MULTI-STEP twin arrow — n copies of the same arrow, each offset
      (dx, dy) down-right from the last, the group centered on the row baseline.
      {"n":2,"dx":10,"dy":10} is the conventional "several steps" double arrow. The
      layout reserves the overhang, so captions and neighbours stay clear.
    arrow.reverse: draw the arrow RIGHT-TO-LEFT (head on the left end). The packer
      still lays the row out left-to-right, so use this when the product sits to the
      LEFT of its starting material - e.g. a branch that folds back along a row.
    COLOUR. Names: white black red yellow green cyan blue magenta (or a raw colortable
      index; custom [r, g, b] colours given elsewhere in the figure are appended from 10).
      * an ARROW is coloured through its style: {"style": {"color": "red"}}
      * a CAPTION through above_color / below_color / yield_color. One value colours the
        whole block, or pass a list parallel to the lines - a caption routinely mixes them
        by line, e.g. below: ["B~2~(pin)~2~", "*via C-C coupling*"] with
        below_color: ["green", null].
    Other verified arrow style keys: dipole (true), hashSpacing (dash pitch), lineWidth
      (thicker shaft AND head; the layout reserves the extra ink), boldWidth,
      equilibriumRatio (asymmetric equilibrium - shows which way it lies; 50 is the
      known working value).
    arrow.arcs: curved arrows attached to this arrow's shaft, for a reagent joining the
      flow or a by-product peeling off it. Each entry: "at" 0..1 along the shaft (fractional,
      so it survives re-lengthening), "off" [dx,dy] to the free end, "dir" "in" (arrives at
      the shaft, default) or "out" (departs from it). By default the arc meets the shaft
      TANGENTIALLY, which is what reads as merging into it; give "bulge" instead for a plain
      bow that crosses at an angle. "head"/"tail" set the arrowheads.
      The layout reserves the arc's footprint, but captions still clear only the shaft.
    arrow.yield: a yield line placed under the below-caption. To put the yield next
      to the PRODUCT instead (the more common convention), use that structure's
      "sublabels" and leave arrow.yield off.
    connectors: from/to arrows between PLACED structures — the only way to express a
      BRANCH POINT (one compound with several outgoing arrows), since rows are plain
      left-to-right sequences with no from/to.
      • kind "straight" is the RIGHT tool for a branch: an ordinary straight arrow at
        any angle, including vertical. On a mostly-vertical arrow the captions rotate
        with it — "above" goes to its right, "below" to its left (conditions-right,
        yield-left, as drawn in the literature). Use this to hang a second or third
        branch off a compound instead of restating that compound on a new row.
      • kind "elbow_vh" / "elbow_hv" / "via" are the routed/turning arrows and stay
        DISCOURAGED — for spanning rows prefer trailing/leading straight arrows.
      • kind "turn" picks between them from where the two structures actually landed:
        a straight drop when they stack, an elbow when the target is off to one side.
        This is what wrap_style='snake' generates for itself; name it by hand when you
        are joining rows whose widths you have not worked out.
      NOTE: between ADJACENT rows the row gap grows to fit a connector's captions.
      A connector spanning more than one row runs through the rows in between, and
      the side captions of a vertical connector get no extra horizontal room --
      leave those rows far enough apart, or check the result.
      [{"from": id, "to": id, "kind": ..., "via"?: [[x,y],...],
      "style"?: {...}, "above"?: [...], "below"?: [...]}]; from/to match the ids above.
    annotations: an electron-pushing arrow, a leader or a bare arc whose two ends are
      in DIFFERENT cells -- the piperidine N attacking the pyridinium ring. Each end
      names a cell and a SMARTS, never atom numbers (abbreviate=True changes the atom
      count, and you cannot know the engine's numbering anyway):
        [{"from": {"cell": "1", "on": "[#7;H1]"},
          "to":   {"cell": "2", "on": "[#6]1[#6][#6][#7+][#6][#6]1"},
          "arrow"?: "full" | "half" | "none",   # half = the single-electron fishhook
          "text"?: "*addition*", "color"?: "red",
          "bulge"?: 0.35,        # sideways swell as a fraction of the chord
          "size"?: 9}]
      A whole-ring SMARTS arrives at the ring centre (what "attacks the ring" means)
      and a single-atom one lands on that atom.
      WITHIN one structure, use the cell's own `highlight` with style "annotate"
      instead: that goes inside the fragment and so moves if the cell moves. This one
      belongs to neither cell, so it is drawn on the page and will NOT follow either
      structure if you drag one in ChemDraw afterwards. It does reserve space, so
      nothing else lands on it.

    reserve: leave room BELOW a row for something a CALLER will draw there -- the
      catalytic cycle that belongs under one step's arrow, drawn by its own engine
      call and parked in the gap.
        [{"after": "8-2", "h": 208}]
      "after" names a cell; the room goes under THAT CELL'S ROW. Height only.
    template: draw every row structure on ONE reference depiction, so a scaffold
      that runs through the route keeps its orientation instead of being laid out
      afresh at each step. "off" (the default) | "auto" (the figure's own first
      structure) | a SMILES of your own.
      OPT-IN HERE, unlike draw_scope_table where it is on. Reach for it when the
      same skeleton runs the length of the scheme and comes out flipped halfway
      along.
      Row cells only, not a reagent above an arrow. The return says how many
      structures matched and how many fell back.
    template_patt: a SMARTS naming WHICH part of `template` to match on, when the
      whole reference is more than the shared core.
    span: fix HOW FAR APART two cells sit, in points -- `reserve` answers "leave
      this much room below", this answers "put these two exactly this far apart".
        [{"from": "III", "to": "II", "w": 308.7}]
      CENTRE TO CENTRE; "from"/"to" name a pair, not a direction. The
      distance is paid for by lengthening the arrow(s) between the two cells, and
      it is split evenly when there is more than one.
      WHEN YOU NEED IT: this scheme has to match a distance fixed by another figure
      (e.g. two intermediates of a separately drawn catalytic cycle). Anchoring one
      end is not enough -- the distance inside this figure has to be set.
      An arrow is never shortened -- its caption sets its length.
    footnotes: running prose under the whole figure, wrapped to the figure's width,
      left-aligned, in 8 pt. One string per note; each is broken across lines as
      needed, and ~sub~ ^sup^ *it* #bold# {color} work as they do everywhere else.
      Sits under EVERYTHING, including a connector routed below the last row.
    footnote_width: WHICH FIGURE the footnotes run under, in points. Omitted, they
      fold to THIS figure's own width -- what "running prose under this figure"
      means, and right for every ordinary call. Pass a number only when this
      drawing is one PART of a bigger one (set beside another tool's output),
      because then the figure the prose runs under is the composed one."""
    # REFUSE A BAD SPEC BEFORE ANYTHING IS BUILT. First statement in the
    # tool, ahead of _new_batch, so a rejected call leaves no per-batch state
    # behind either. See preflight.py for what it checks and why the engines'
    # own key whitelists are not enough.
    PF.preflight_linear({"items": items, "connectors": connectors, "page": page,
                         "align": align, "reserve": reserve, "span": span,
                         "template": template,
                         "slots": slots})
    # A MISSPELLED STYLE MUST NOT BE A QUIET STRAIGHT WRAP. draw_snaking_scheme
    # takes Lx/Rx and never reads them -- four wildly different pairs, one
    # identical file -- and nothing says so. A new mode word that silently means
    # "no" would be the same trap with a new name.
    if wrap_style not in ("straight", "snake"):
        raise ValueError("wrap_style is 'straight' or 'snake', got %r" % (wrap_style,))
    _new_batch()
    rows_in = items if (items and isinstance(items[0], list)) else [items]
    # BEFORE anything is built. This loop rebuilds each item into the engine's own
    # vocabulary ("s" where the caller writes "smiles"), so an unknown key would be
    # dropped here and layout_linear would never see it to complain -- which is
    # exactly how a cell can carry "highlight" and draw nothing.
    _multirow = bool(items) and isinstance(items[0], list)
    for _ri, _row in enumerate(rows_in):
        for _ci, _it in enumerate(_row):
            _check_linear_item(_it, "items[%d][%d]" % (_ri, _ci) if _multirow
                               else "items[%d]" % _ci)
    for _ci, _con in enumerate(connectors or []):
        CY.check_keys(_con, _LINEAR_CONN_KEYS, "connectors[%d]" % _ci)
    # The colortable belongs to the DOCUMENT, so every cell's and every reagent's
    # custom colours are appended once, up front, and the specs rewritten to the
    # indices they landed on.
    _specs = [it.get("highlight") for r in rows_in for it in r if it.get("highlight")]
    _specs += [it["arrow"].get("reagent_highlight") for r in rows_in for it in r
               if "arrow" in it and (it["arrow"] or {}).get("reagent_highlight")]
    _hdr = _extend_palette(C.acs_header_prefix(), _specs)
    # ONE REFERENCE FOR THE WHOLE ROUTE. Same two sentinels the scope table reads
    # ("auto" -> the figure's own first structure, "off" -> none), resolved here
    # for the same reason it resolves them before building anything: a template is
    # a set of COORDINATES to donate, so it has to be a real drawable structure by
    # the time the first cell is built.
    #
    # DEFAULT "off", NOT "auto" LIKE THE TABLE. A scope table is n drawings of one
    # scaffold and matching them is the whole point; a scheme is a ROUTE, where
    # consecutive compounds differ by a bond and the caller may well want each one
    # drawn the way its own chemistry reads. So this is opt-in, and every figure
    # drawn before it existed is byte-identical.
    #
    # ROW CELLS ONLY. A reagent above an arrow is shrunk to fit the arrow
    # (see DECOR_UNSCALED), so donating coordinates to it buys a shape that is
    # then scaled anyway.
    _tmpl = None
    if str(template).lower() == "auto":
        _tmpl = next((it["smiles"] for r in rows_in for it in r if it.get("smiles")),
                     None)
    elif str(template).lower() != "off":
        _tmpl = template

    structures: dict = {}
    spec_rows: list = []
    for row_items in rows_in:
        row: list = []
        for it in row_items:
            if "smiles" in it:
                sid = _register(structures, it["smiles"], it,
                                abbreviate=abbreviate, template=_tmpl,
                                template_patt=template_patt)
                cell = {"s": sid}
                if not str(it.get("id") or "").strip():
                    cell["num"] = False        # same rule as the grid's
                if it.get("num") is False:
                    # SAME HAZARD AS `rowspan` BELOW: whitelisted in
                    # _LINEAR_CELL_KEYS, consumed by layout_linear
                    # (`num=cell.get("num", True)`), and so it must be copied
                    # across in between -- otherwise `num: False` is accepted and
                    # does nothing, and the number prints anyway.
                    #
                    # What needs it: a resonance pair. A published mechanism
                    # figure draws [A <-> B] under ONE compound number, and the
                    # two forms cannot share an
                    # id (preflight refuses it, correctly -- `structures` is keyed
                    # by id, so the second form would be silently replaced by the
                    # first). Giving B its own id and suppressing its number is
                    # how you write that.
                    #
                    # Only False is honoured. `num: True` is already the default,
                    # and on a cell with no id there is nothing to print -- so
                    # there is no case where True means anything the default does
                    # not already do.
                    cell["num"] = False
                if it.get("label"):
                    cell["label"] = it["label"]
                if it.get("sublabels"):
                    cell["sublabels"] = it["sublabels"]
                if it.get("bracket"):
                    cell["bracket"] = True
                if it.get("under"):
                    cell["under"] = str(it["under"])
                if it.get("rowspan"):
                    # THE KEY IS WHITELISTED, SO IT MUST BE COPIED HERE.
                    # _LINEAR_CELL_KEYS is built by reading the code that
                    # CONSUMES each key -- rowspan.py does consume this one -- so
                    # the validator waves it through, and a loop that did not
                    # copy it would give the caller ok:true and a figure with no
                    # fan: three unrelated reactions where a divergence was
                    # asked for.
                    # Validated in rowspan.plan_spans (>= 2, must not run past the
                    # last row, root must be a structure), not here: one rule, one
                    # home, and that one already says it with row/slot indices.
                    cell["rowspan"] = it["rowspan"]
                row.append(cell)
            elif "text" in it:
                # A species written rather than drawn. Same keys a structure cell
                # takes, because on paper it occupies the same slot: "+ 4a" is a
                # product with a number under it and no structure above.
                cell = {"text": str(it["text"])}
                for k in ("id", "label", "sublabels", "size"):
                    if it.get(k) is not None:
                        cell[k] = it[k]
                row.append(cell)
            elif "arrow" in it:
                arrow = dict(it["arrow"])
                if arrow.get("reagent_smiles"):
                    # a reagent above an arrow is a structure like any other, so it
                    # takes the same keys -- under "reagent_" names, because it
                    # shares its dict with the arrow's own.
                    # It never carries a caller id, so `prefix` is what names it.
                    _rsmi = arrow.pop("reagent_smiles")
                    rid = _register(
                        structures, _rsmi,
                        {"highlight": arrow.pop("reagent_highlight", None),
                         "rlabels": arrow.pop("reagent_rlabels", None)},
                        prefix="_r", abbreviate=abbreviate,
                        where="reagent above arrow (%s)" % _rsmi)
                    arrow["reagent"] = rid
                    if arrow.get("reagent_id"):
                        arrow["reagent_label"] = str(arrow.pop("reagent_id"))
                row.append({"arrow": arrow})
            elif it.get("plus"):
                row.append({"plus": True})
        spec_rows.append(row)
    recipe = {"structures": structures, "rows": spec_rows}
    if footnotes:
        recipe["footnotes"] = _footnote_lines(footnotes)
    if _DECOR:
        recipe["decorations"] = dict(_DECOR)
    if connectors:
        recipe["connectors"] = connectors   # straight = branch points; elbow/via discouraged
    if reserve:
        recipe["reserve"] = reserve
    if span:
        recipe["span"] = span
    cfg = {"row_align": align}
    if footnote_width:
        cfg["footnote_width"] = float(footnote_width)
    if max_arrow:
        cfg["max_arrow"] = max_arrow
    wrap_width = C.usable_width(page, 30.0) if wrap else None
    # HAND THE ENGINE AN ALLOCATOR. Without one, layout_linear takes the
    # `alloc is None` path through every put_box and through the whole connector
    # block, so not one reservation is made and the collision report is empty --
    # not "clean", EMPTY. The figure comes out exactly the same either way, which
    # is what makes the difference invisible: the overlaps are still there, only
    # nobody is measuring them. Same pattern as draw_scope_table.
    if annotations:
        for _ai, _an in enumerate(annotations):
            CY.check_keys(_an, _ANNOT_KEYS, "annotations[%d]" % _ai)
            for _end in ("from", "to"):
                _e = _an.get(_end)
                if not isinstance(_e, dict) or not _e.get("cell") or not _e.get("on"):
                    raise ValueError(
                        "annotations[%d].%s: needs {\"cell\": <an id in this "
                        "figure>, \"on\": <a SMARTS>}, got %r"
                        % (_ai, _end, _e))
                if str(_e["cell"]) not in structures:
                    raise ValueError(
                        "annotations[%d].%s: cell %r is not in this figure. The "
                        "ids are %s." % (_ai, _end, _e["cell"],
                                          ", ".join(map(repr, sorted(structures)))))
        recipe["annotations"] = [
            {**a, "from": {**a["from"], "cell": str(a["from"]["cell"])},
             "to": {**a["to"], "cell": str(a["to"]["cell"])}}
            for a in annotations]
        recipe["anchor_src"] = dict(_ANCHOR_SRC)
    from allocator import Allocator as _Alloc
    _alloc = _Alloc()
    L.DECOR_UNSCALED.clear()
    L.ACCOUNT_DIFF.clear()
    L.CAPTION_MOVED.clear()
    L.CONN_DETOURED.clear()
    L.SPAN_SHORT.clear()
    L.WRAP_ROWS.clear()
    _shift: list = []
    _placement: dict = {}
    import anchors as _ANC
    _frame = _ANC.Anchors()
    doc = L.render_linear(recipe, _hdr, cfg=cfg, page_mode=page,
                          wrap_width=wrap_width,
                          serpentine=(wrap_style == "snake"),
                          alloc=_alloc, shift=_shift,
                          frame=_frame, placement=_placement)
    doc = _add_graphics(doc, graphics, slots, structures, abbreviate,
                        alloc=_alloc, shift=_shift)
    _path = _write_and_open(doc, filename, clean_up)
    out = _wrap_note(wrap_width)
    out += _span_short_note()
    if L.CONN_DETOURED:
        # SAME TERMS AS CAPTION_MOVED. The arrow now bends, and a reader who is not
        # told why sees an arrow that swings wide for no reason. Naming what was in
        # the way is also the actionable part: moving THAT cell is what gets the
        # straight line back.
        out += ("\n\nNOTE: %d connector(s) could not run straight without crossing "
                "something and were routed around it: %s. The figure reads "
                "correctly and nothing is drawn over, but a bent connector is not "
                "the straight arrow this engine's conventions ask for -- move what "
                "is named out of the way, or put the two cells in adjacent rows, "
                "to get it back."
                % (len(L.CONN_DETOURED),
                   "; ".join("%s bent aside to x=%.0f, past %s" % (o, x, ", ".join(w))
                             for o, x, w in L.CONN_DETOURED[:4])))
    if L.CAPTION_MOVED:
        out += ("\n\nNOTE: %d connector caption(s) could not sit where the "
                "convention puts them (a gap from the shaft) and were moved clear "
                "of what was in the way: %s. The figure reads correctly, but the "
                "caption now sits further from its arrow than it should -- usually "
                "because the arrow itself runs through a structure, which the "
                "LAYOUT REPORT above will say if so. Widening the gap between the "
                "cells this connector bridges puts it back."
                % (len(L.CAPTION_MOVED),
                   "; ".join("%s moved %.0f pt %s (%s)" % (o, d, dirn, how)
                             for o, dirn, d, how in L.CAPTION_MOVED[:4])))
    if L.DECOR_UNSCALED:
        out += ("\n\nNOTE: %d reagent structure(s) above an arrow were shrunk to fit "
                "the arrow, and a halo / dot / locant cannot be shrunk with them "
                "(the disc keeps its size in absolute coordinates), so THOSE MARKS "
                "WERE NOT DRAWN: %s. The structure itself is there and correct. To "
                "get the marks, give that step a longer arrow -- a longer caption, "
                "or raise max_arrow -- so the reagent fits at full size."
                % (len(L.DECOR_UNSCALED),
                   ", ".join("%s (shrunk to %.0f%%)" % (i, f * 100)
                             for i, f in L.DECOR_UNSCALED[:4])))
    _lrep = _layout_report(_alloc)
    if _lrep:
        out += "\n\nLAYOUT REPORT -- " + _lrep
    return _figure_report(_path, structures, alloc=_alloc,
                          notes=out + _chem_report(), placement=_placement)


# A section marker in `rows`. Strict, because a typo here would otherwise turn a
# section heading into a product cell with no SMILES and vanish.
_SECTION_KEYS = frozenset({"section", "template", "color", "style", "letter",
                           "text_color"})


def _scope_cell(c, where: str) -> dict:
    """A product cell, or a ValueError saying which one and why it is not.

    THESE USED TO BE DROPPED, NOT REFUSED. Only cells with a 'smiles' were kept,
    and the key check runs later, on what was kept -- so {"smile": ...},
    {"SMILES": ...} or a cell with an id and a yield but no structure vanished
    from the table while the reply said dropped: []. A misspelt 'yield' was
    caught; a misspelt 'smiles' never was."""
    if isinstance(c, dict) and isinstance(c.get("smiles"), str) and c["smiles"].strip():
        return c
    if not isinstance(c, dict):
        raise ValueError("%s: a product cell is an object like {\"smiles\": ..., "
                         "\"id\": ...}, got %r" % (where, c))
    name = " (id %r)" % c["id"] if c.get("id") not in (None, "") else ""
    if "smiles" in c:
        raise ValueError("%s%s: 'smiles' is empty, so there is nothing to draw. "
                         "Give the structure, or leave the cell out." % (where, name))
    near = difflib.get_close_matches("smiles", [str(k).lower() for k in c], n=1, cutoff=0.6)
    hint = (" Did you mean 'smiles'? (got %r)" % next(k for k in c if str(k).lower() == near[0])
            if near else "")
    raise ValueError("%s%s has no 'smiles', so it would be left out of the table. "
                     "Its keys are: %s.%s" % (where, name, ", ".join(map(str, c)), hint))


def _split_sections(rows: list):
    """Cut `rows` into sections at every {"section": ...} marker.

    Returns (groups, colours). A group is
        {"title": str|None, "template": str|None, "rows": [[cell, ...]],
         "explicit": bool}
    where explicit means the caller wrote the rows out literally. The first group
    has title None, so a call with no marker in it comes back as exactly one
    untitled group -- which is what keeps the no-section path byte-identical.

    colours is [(title, colour)] for markers that asked for one, which is not
    wired yet; the tool says so rather than dropping it.
    """
    seq = list(rows or [])
    groups = [{"title": None, "template": None, "rows": [], "explicit": False}]
    colours: list = []
    loose: list = []

    def flush():
        if loose:
            groups[-1]["rows"].append(list(loose))
            loose.clear()

    def is_marker(x):
        return isinstance(x, dict) and x.get("section") is not None

    for i, it in enumerate(seq):
        if is_marker(it):
            bad = set(it) - _SECTION_KEYS
            if bad:
                raise ValueError(
                    "rows[%d]: a section marker takes %s, got unknown key(s) %s"
                    % (i, sorted(_SECTION_KEYS), sorted(bad)))
            flush()
            if it.get("color"):
                colours.append((str(it["section"]), str(it["color"])))
            groups.append({"title": str(it["section"]),
                           "template": it.get("template"),
                           "style": it.get("style"),
                           "letter": it.get("letter"),
                           "color": it.get("color"),
                           "text_color": it.get("text_color"),
                           "rows": [], "explicit": False})
        elif isinstance(it, list):
            if any(is_marker(x) for x in it):
                raise ValueError(
                    "rows[%d]: a section marker goes BETWEEN rows, not inside one "
                    "-- move it out of the list" % i)
            cells = [_scope_cell(c, "rows[%d][%d]" % (i, j)) for j, c in enumerate(it)]
            if cells:
                flush()
                groups[-1]["rows"].append(cells)
                groups[-1]["explicit"] = True
        else:
            loose.append(_scope_cell(it, "rows[%d]" % i))
    flush()
    # A section whose products arrived as one loose run gets PACKED (row shape
    # derived from the page); one that wrote its rows out keeps them literally,
    # and any loose run in the same section becomes one more literal row.
    groups = [g for g in groups if g["rows"]]
    return (groups or [{"title": None, "template": None, "rows": [],
                        "explicit": False}]), colours


def _section_template(group: dict, table_tmpl):
    """This section's depiction reference.

    Omitted -> inherit whatever the table set, which is what makes one shared
    scaffold across sections work. "auto" -> this section's OWN first product,
    which is the point of the feature: an aliphatic section templated on an aryl
    parent matches nothing and every cell falls back.
    """
    t = group.get("template")
    if t is None:
        return table_tmpl
    s = str(t)
    if s.lower() == "off":
        return None
    if s.lower() == "auto":
        first = next((c for r in group["rows"] for c in r), None)
        return first["smiles"] if first else None
    return t


def _span_short_note() -> str:
    """The `span` shortfall sentence, for both tools that take one.

    SHARED, unlike the four notes beside it, because THIS one is the request whose
    result is invisible. Every other note is about something a reader can see went
    wrong -- a bent connector, a displaced caption. A span that came up short draws
    a perfectly ordinary figure, and its caller is a PARENT that has already
    anchored the other end of this chain to another engine's drawing, so what they
    will actually see is the OTHER figure's structure apparently misplaced. Two
    tools saying that in two slightly different ways is two chances to say it less
    well.
    """
    if not L.SPAN_SHORT:
        return ""
    # TWO KINDS OF SHORTFALL AND ONLY ONE REMEDY EACH, so the explanation is
    # attached to the kind that is actually present rather than to both.
    why = ""
    if any(g is not None for *_, g in L.SPAN_SHORT):
        why = (" An arrow is only ever LENGTHENED for a span -- its own caption sets "
               "its natural length, and a shorter shaft prints the conditions past the "
               "head -- so a distance under the one these cells already sit at is "
               "reported rather than granted. To close it, narrow the cells themselves, "
               "or move the other part of the composite figure instead.")
    if any(g is None for *_, g in L.SPAN_SHORT):
        why += (" A span whose cells were auto-wrapped into different rows is not "
                "a horizontal distance any more; switch `wrap` off, or break the "
                "rows yourself so the pair stays together.")
    return ("\n\nNOTE: %d span(s) were not reached: %s.%s"
            % (len(L.SPAN_SHORT),
               "; ".join(("%s..%s asked %.0f pt and the two were auto-wrapped into "
                          "different rows" % (a, b, w)) if g is None else
                         ("%s..%s asked %.0f pt, came to %.0f" % (a, b, w, g))
                         for a, b, w, g in L.SPAN_SHORT[:4]), why))


def _wrap_note(wrap_width) -> str:
    """What `wrap` did to the rows, read off layout_linear.WRAP_ROWS, in the
    caller's row numbers.

    REPLACES "Auto-wrapped to page width.", which was printed whenever wrap was on
    and the page had a width -- for several rows the wrap had not touched, and for
    one row that already fitted. Silent when nothing was broken and nothing is
    left over the column: a note saying that nothing happened is the blanket line
    again.

    A row whose ONLY piece is still too wide was not broken at all, so it is
    listed with the rows left whole, under its own reason, rather than as "part 1
    after breaking". The sentence about spans is added only when a span is one of
    the reasons.
    """
    if not wrap_width or not L.WRAP_ROWS:
        return ""

    def name(r):
        return "the row" if r is None else "items[%d]" % r

    broken = {e[1] for e in L.WRAP_ROWS if e[0] == "broke"}
    snake, broke, left, after, by_span = [], [], [], [], False
    for e in L.WRAP_ROWS:
        if e[0] == "snake":
            snake.append("WRAP: wrap_style='snake' breaks a single flat row only, so "
                         "the %d rows of items were left as they are; "
                         "wrap_style='straight' breaks the over-wide ones." % e[1])
        elif e[0] == "broke":
            broke.append("%s was broken into %d rows" % (name(e[1]), e[2]))
        elif e[0] == "whole":
            left.append("%s, %.0f pt (%s)" % (name(e[1]), e[2], e[3]))
            by_span = True
        else:                                                    # "over"
            _k, r, j, w, m = e
            why = ("everything from the merge into %s on stays together" % m if m
                   else "no break is allowed inside it -- wrap breaks only "
                        "between cells")
            if r in broken:
                after.append("%s, part %d, %.0f pt (%s)" % (name(r), j + 1, w, why))
            else:
                left.append("%s, %.0f pt (%s)" % (name(r), w, why))
    out = list(snake)
    if broke or left or after:
        s = "WRAP (%.0f pt column):" % wrap_width
        if broke:
            s += " %s." % "; ".join(broke)
        if left:
            s += " Left whole although wider than the column: %s." % "; ".join(left)
            if by_span:
                s += (" A row a rowspan covers, or one that owns a fork, is never "
                      "broken -- the span counts rows -- and a covered row is "
                      "measured together with the span cell it faces.")
        if after:
            s += " Still wider than the column after breaking: %s." % "; ".join(after)
        out.append(s)
    return "\n\n".join(out)


@mcp.tool()
@_cx_doc
def draw_scope_table(rows: list, filename: str = "scope", clean_up: bool = False,
                     page: str = "letter", header: dict = None, title: str = None,
                     footnotes: list = None, cols: int = None,
                     label_align: str = "center", struct_valign: str = "bottom",
                     template: str = "auto", template_patt: str = None,
                     abbreviate: bool = False,
                     last_row: str = "auto", fill_width: bool = True,
                     span: float | None = None,
                     table: dict = None, table_notes: list = None,
                     table_style: dict = None, order: list = None,
                     panel_pos: str = "above", frame: str = "none",
                     graphics: list = None, slots: list = None,
                     boxes: list = None,
                     connectors: list = None, panels: list = None,
                     footnote_width: float | None = None) -> str:
    """Draw a SUBSTRATE-SCOPE TABLE (the n x m product grid) as a ChemDraw (CDXML) file.

    The other common methodology-paper figure besides a scheme: one generic
    reaction on top, then a grid of products each with its number, yield and any
    er / dr / footnote marker.

    rows: a flat list of product cells, OR a list of rows for a RAGGED grid.
      Give explicit rows whenever the shape is uneven -- real tables are often
      3,3,3,4,4 because the first rows draw the aryl group in full and the last
      ones abbreviate it. No single column count can express that.
      A product cell is:
        {"smiles": "...", "id": "3aa", "yield"?: "92% yield,",
         "sublabels"?: ["96:4 er & 7:1 dr^[b]^"], "label"?: "2HCl", "foot"?: "b",
         "rlabels"?: {"1": "Ar"}, "floating"?: [1], "legend"?: {...},
         "highlight"?: [...]}
      "highlight" is the WHOLE spec draw_molecule takes, all of it -- "color",
      "fill", "squiggle", the five bond styles, and "halo" / "dot" / "locants".
      "id" is the bold compound number under the structure. Put er/dr and any
      footnote marker in sublabels -- one string per line, in order.
      A cell may be a MARKUSH structure: write the placeholder as * or [*:n] in the
      SMILES (see draw_molecule for the full rules) and name it with "rlabels".

      SECTIONS. Anywhere in the rows sequence you may put a marker instead of a
      cell, which ends the current row and starts a new block under a heading:
        {"section": "(het)aryl alkynes", "template"?: "<SMILES>"}
      ("color" fills the title badge of style "pill", and "text_color"
      colours that title.) Each section sizes its OWN cells, so two sections may have
      different column counts and widths, and short rows are set flush left.
      "template" per section is the real point: omitted it inherits the table's,
      "auto" re-picks from that section's own first product. An aliphatic section
      templated on an aryl parent matches nothing and every cell falls back --
      give such a section its own reference.

    template: THE thing that makes a scope table readable. Fixes the shared
      scaffold to one depiction so every cell draws it identically and the eye
      lands on the substituent instead of on a rotated core.
        "auto" (default) -- use the FIRST product as the reference
        "<SMILES>"       -- an explicit reference
        "off"            -- independent layout per molecule
      SMILES, NOT SMARTS: the template has to be a structure that can be drawn.
      Any series member works as the reference: as a substructure query a phenyl
      matches a tolyl and a methyl matches the first carbon of an ethyl. Products
      that fail to match, or that the template did not actually help, fall back
      to independent layout and are named in the return value.
      HOW MANY fell back is `template_match` in the return value:
      {"attempted", "placed", "fallback"}. The prose notes name which molecules and
      why.
    template_patt: a SMARTS narrowing WHICH PART of `template` has to match, and
      the answer whenever the shared core is a pattern rather than a molecule.
      Two arguments, two jobs, and they are not interchangeable:
        template      the drawable structure whose coordinates get copied  (SMILES)
        template_patt the part of it that has to match each member         (SMARTS)
      So a sulfoximine series aligns on a real member with the reactive core named
      as the pattern -- template="CS(=NC=O)C",
      template_patt="[#6]S(=N[#6]=O)[#6]" -- and members the plain template misses
      (an aryl one, say) come through placed. Omitted, the whole template must
      match.
      Applies to whatever reference is in force, including a section that re-picks
      its own with "auto" -- if the pattern is not inside that section's reference
      those cells fall back, and `template_match` is where that shows.

    header: the generic reaction above the table. Same cell shapes as
      draw_linear_scheme, in one row:
        {"items": [{"smiles": "...", "id": "1a"}, {"plus": true},
                   {"smiles": "...", "id": "2"},
                   {"arrow": {"above"?: [...], "below"?: [...]}},
                   {"smiles": "...", "id": "3", "templated"?: true}]}
      A header structure takes "rlabels", "floating" and "highlight" as well, on the
      same terms as a product cell. {"text": ...} is a written species (see A SPECIES YOU WRITE RATHER THAN DRAW).
      A rule is drawn under it. A header item may be a generic structure: write the
      placeholder as * or [*:n] and add "rlabels" to name it, same as a product
      cell. Set "templated": true on the product to depict it on the series
      template -- but not on a generic one, which will not match explicit cells.
      THIS WORKS WITH NO PRODUCT GRID, which is the optimisation figure: one
      reaction over a table of conditions. Give `template` a SMILES -- normally
      the substrate -- and the templated species are drawn on its scaffold, so
      the product does not come out as the substrate's mirror image. With
      `template: "auto"` and no products the reference is the FIRST header
      structure. When the two ends differ in aromaticity (an arene oxidised to a
      quinone) no SMILES matches both: name the shared skeleton in
      `template_patt` with `~` bonds, e.g.
      "[#6]1~[#6]~[#6]~[#7]~[#6]2~[#6]~1~[#6]~[#6]~[#6]~[#6]~2".

    title: the figure's heading, set above everything. Bold the label inline:
      "#Table 2:# Substrate scope of alkenes".
    footnotes: running prose under the table, wrapped to the table width. One
      string per note; each is broken across lines as needed.
    footnote_width: WHICH FIGURE the footnotes run under, in points. Omitted, they
      fold to this table's own width -- right for every ordinary call. Pass a
      number only when this table is one PART of a bigger figure, because then the
      thing the prose runs under is the composed figure, not the table.

    table: a DATA TABLE (the optimisation table: entry / catalyst / base /
      solvent / yield), drawn under -- or over -- the structure panel. Text
      only: a
      cell refers to a structure by its short code (`L1`, `Co-5`) and the panel
      holds the structure itself.
        {"cols": [{"header": "yield (%)", "align"?: "center|left|right",
                   "note"?: "b", "id"?: "y_3a"}, ...],
         "rows": [{"cells": ["1", "Cs~2~CO~3~", "59"],
                   "note"?: "c",              # mark on the entry number
                   "mark"?: "best" | "group:<name>"}, ...],
         "groups"?: [{"header": "yield (%)", "note": "b",
                      "cols": ["y_3a", "y_4a"]}],   # spans ADJACENT columns
         "split"?: 2,                # side by side, header repeated in full
         "style"?: "default" | "acs",
         "width_mode"?: "natural" | "match"}
      A cell is a string, or {"v": "80", "sec": "77", "note": "c"} for a value
      with a bracketed secondary reading -> `80 (77)^c`. "v" may be None: a
      bracket with no main value is legal.
      "mark" keeps two things apart on purpose: "best" says this row is the
      answer, "group:x" says these adjacent rows belong together (they share one
      band, and different groups alternate colours).
      Cell text is normalised -- `--`, `–`, blank all become an em dash, `nd` /
      `N.D.` become `n.d.`, `Trace` becomes `trace`, `> 95:5` closes up.
      "width_mode" is how wide the table comes out:
        "natural" (default) -- every column is as wide as its own widest string
          and the table is centred under whatever else is in the figure.
        "match" -- the columns absorb the slack until the table spans the WHOLE
          figure width, so its bottom rule lines up with the scheme above it.
      "match" is CAPPED at 1.30x the natural width (table_style["stretch_max"]),
      so a narrow table under a wide scheme may still come out short of the
      figure. Nothing ever compresses: a target narrower than the
      natural width is ignored.
    table_notes: the footnote block under the data table, as
      [{"mark": "a", "text": "..."}, {"mark": None, "text": "n.d., not detected."}].
      mark=None is the unmarked abbreviation line.
    table_style: overrides for the table's own layout, e.g. {"zebra": [0.95, 0.95,
      0.95]} or {"vrule": True}. Keys: size, header_size, col_pad, row_pad, rule_top, rule_header,
      rule_bottom, rule_group, rule_vertical, vrule, gutter, notes_size,
      notes_gap, notes_lead, header_bold, header_fill, best_fill, group_fills,
      zebra, stretch_max, width_mode. Two of these belong to
      width_mode "match": {"stretch_max": 1.6} raises the 1.30 cap on how far the
      table may be pulled past its natural width (>= 1.0; 1.0 means "never
      stretch"), and "width_mode" itself is accepted here too.
    order: the VERTICAL ORDER of the figure's blocks, as an explicit sequence of
      "scheme" (the header reaction), "grid" (the product/ligand panel) and
      "table" (the data table):
        ["scheme", "grid", "table"]   default -- reaction, products, data table
        ["scheme", "table", "grid"]   reaction, data table, ligand panel below it
        ["scheme", "table"]           the optimisation figure: no product grid
        ["table", "scheme"]           table first
      A name for a block you did not pass is ignored. Omit `order` and it comes
      from panel_pos.
    panel_pos: DEPRECATED, kept working -- use `order`. "above" (structures above
      the data table, the usual layout) = ["scheme","grid","table"];
      "below" = ["scheme","table","grid"]. Ignored when `order` is given.
    frame: "none" | "enclose" -- raise the table's top rule ABOVE the reaction
      scheme so the scheme sits inside the table's frame.
    panels: a lettered, titled, framed box holding a small reaction, placed
      relative to a block or to an earlier panel -- e.g. panels A/B/C in a row, one
      panel tucked into another's bottom-right corner, a ligand panel beside a data
      table.
        [{"name": "A",                      # the letter drawn as "(A)", and how a
                                            #   later panel refers to this one
          "title": "Aroyl chlorides",
          "items": [...],                   # SAME cells `header` takes: structures,
                                            #   {"plus": true}, one {"arrow": {...}}
                                            #   with captions, {"text": ...}
          "frame": "Plain",                 # omit for no frame
          "of": "grid",                     # "scheme" | "grid" | "table", or the
                                            #   name of a panel placed EARLIER
          "side": "below",                  # right | left | above | below
          "align": "end",                   # start | center | end -- "end" is what
                                            #   lines up two panels' bottom edges
          "gap": 26}]
      Panels are placed IN ORDER, so "of" may only name a panel before it.
      A panel holds one reaction row; it cannot contain a product grid or a data
      table.
    connectors: an arrow BETWEEN two product cells, for a small transformation drawn
      inside the grid (e.g. 20 -[O]-> 21, where each side keeps its own number
      and yield). Straight only -- two neighbours in a grid have nothing to route
      around, so there is no "kind":
        [{"from": "20", "to": "21", "above"?: ["[O]"], "below"?: ["88%"],
          "style"?: {...}}]
      Captions ROTATE with the run, as they do on a scheme's connectors: above/below
      for a left-to-right arrow, right/left for a vertical one (conditions-right,
      yield-left). A vertical arrow leaves BELOW the cell's compound number and
      yield, not below its structure.
      THE COLUMNS ARE NOT WIDENED for it. Two structures sit 40 pt apart, and a long
      caption in that gap prints over the cells. Give the caption fewer words, or
      fewer columns.

    boxes: a frame round something the figure already drew --
        [{"around": "grid", "style"?: "Dashed", "pad"?: 8}]
      `around`: "figure" | a block ("scheme" | "grid" | "table") | "section:0" |
      ONE header species ("hdr:L4") | ONE product cell ("cell:3a"). The last two
      use the id you gave that structure, and they frame it WITH its compound
      number and yield. A box round SEVERAL cells is not built -- put those cells
      in their own `section` and box that.
      `style` is an OUTLINE ("Plain" default, "Dashed", "RoundEdge", "RoundEdge
      Dashed") or a FILL ("Filled", "Shaded", "RoundEdge Filled", "RoundEdge
      Shaded"). A fill takes `color` and is drawn UNDER everything, so it is the
      tinted panel behind a block or a section:
        [{"around": "section:1", "style": "RoundEdge Shaded", "color": "yellow"}]
      `pad` defaults to 6 pt.
      THE LIGAND BESIDE THE ARROW is `hdr:<id>`, and put the ligand in
      `header.items` rather than in `panels`: as a header item it shares the
      reaction row's width and alignment, where a panel costs the figure a column
      of its own.
    section style: a section marker also takes "style" ("rule" default: lead
      line, title, rule out | "plain": title alone | "centered": line, title
      centred, line) and "letter", the bold a / b / c at the left margin --
        {"section": "aliphatic amines", "style": "centered", "letter": "b"}
    cols: force a column count. Omitted, it is derived from the page width, then
      the products are spread over the fewest rows that hold them.
    span: how wide the figure MAY be, in points, instead of asking the page.
      It answers three questions: how many columns the
      products fall into (when `cols` is not given), how far `fill_width` spreads
      the cells, and when the figure is reported as too wide. Pass it when this
      table is one PART of a larger figure -- a part has no page of its own, and
      the width it has to live within is its sibling's. Omitted, the page decides.
      A BUDGET, NOT A TARGET: cells are widened to fill the span, but the ink
      stops where the edge cells' own drawings stop, so it comes out narrower:
          ink = span * (cols - 1) / cols + (leftmost + rightmost content) / 2
      The return's `size` reports the ink, so a caller who needs a particular INK
      width reads it back and solves -- two passes are exact.
      Outside the middle regime the law does not apply: a span under the natural
      width moves nothing, and one past max_spread (2x natural) is capped. A cell
      is only ever WIDENED, never narrowed.
    last_row: "auto" (centre a short final row) | "center" | "left".
    label_align: "center" (numbers centred under each structure, the usual house
      style) | "left" (structure and label share a left edge).
    struct_valign: "bottom" (default -- every structure in a row sits on one
      baseline) | "center".
    fill_width: widen the cells so the grid spans the column instead of huddling
      to the left. Capped at 2x the natural cell width.
    page: 'letter' (default) | 'double' | 'single' | 'fit' | 'auto'.
    Text markers: ~sub~ ^sup^ *italic* #bold# _underline_.
    """
    # REFUSE A BAD SPEC BEFORE ANYTHING IS BUILT. First statement in the
    # tool, ahead of _new_batch, so a rejected call leaves no per-batch state
    # behind either. See preflight.py for what it checks and why the engines'
    # own key whitelists are not enough.
    PF.preflight_scope_table({"rows": rows, "header": header, "page": page,
                              "span": span,
                              "struct_valign": struct_valign,
                              "connectors": connectors, "panels": panels,
                              "slots": slots,
                              "template": template})
    _new_batch()

    groups, colours = _split_sections(rows)
    _title_box = [title]          # see the texts= note at _extend_palette below
    flat = [c for g in groups for r in g["rows"] for c in r]
    # BEFORE the first _register. _register builds the fragment immediately, so a
    # highlight's colour is spent at that moment -- an [r, g, b] rewritten after
    # it would name an index the drawing never saw. THREE sources here, and
    # missing one is silent: the products, the header block's own structures, and
    # the panels'. (Same reason the tool validates all three: a header structure
    # takes the whole highlight spec, and a panel is a reaction drawn inside the
    # figure, not decoration.)
    _hdr = _extend_palette(C.acs_header_prefix(), (
        [c.get("highlight") for c in flat if c.get("highlight")]
        + [it.get("highlight") for it in ((header or {}).get("items") or [])
           if isinstance(it, dict) and it.get("highlight")]
        + [it.get("highlight") for p in (panels or [])
           for it in ((p or {}).get("items") or [])
           if isinstance(it, dict) and it.get("highlight")]
        # A BOX'S FILL IS A COLOUR LIKE ANY OTHER. The eight names in the ACS
        # table are SATURATED, and a saturated panel reads as a highlighter pen;
        # tinted panels in published figures are pastels. _extend_palette already
        # appends an [r, g, b] to this document's own colortable and rewrites the
        # spec in place -- boxes have to be handed to it too, or a fill could
        # only ever be one of the eight.
        + ([boxes] if boxes else [])
        # GROUPS, NOT THE RAW MARKERS. _split_sections has already run (nine lines
        # up) and COPIED each marker's colour into its group dict, so rewriting
        # the marker here would rewrite a value nobody reads again. The groups are
        # what becomes `sections`, so they are the copy that reaches the packer.
        + [groups]),
        keys=("color", "text_color"),
        # {#RRGGBB} inside text. What is passed is the **container**: the rewrite
        # is in place, and title is a bare string, so it only comes back if it is
        # wrapped in a one-element list. The whole rows/header/panels/table tree
        # is handed over rather than listing fields one by one -- whichever field
        # is missed raises "unknown colour" quietly.
        texts=[rows, header, panels, table, table_notes, footnotes, boxes,
               _title_box,
               # GROUPS TOO, and for the same reason they are in spec_lists
               # above: _split_sections has already run, and each section title
               # is **copied** into its group, so rewriting the marker in rows
               # cannot reach that copy. The section `color` has the same trap.
               groups])
    # A SCHEME is a block too. A guard demanding products or a table would refuse
    # `rows: []` plus a header -- one generic reaction, drawn alone -- although
    # the engine can draw it. Raise only when NOTHING resolves.
    hdr_items = [it for it in ((header or {}).get("items") or [])
                 if it.get("smiles") or it.get("plus") or it.get("arrow")]
    if not flat and not table and not hdr_items:
        return ("nothing to draw: no product cells with a 'smiles' key, no "
                "'table', and no 'header' items")

    # THE REFERENCE every "templated" structure is drawn on. Gating this on
    # `flat` would make the mechanism unreachable on the figure that needs it
    # most: an OPTIMISATION table is one reaction over a table of conditions and
    # has NO product grid, so `flat` is empty, `tmpl` stays None, and a header
    # item's "templated": true is accepted and does nothing at all.
    #
    # WHAT THAT COSTS, measured on a quinoline -> quinoline-5,8-dione pair:
    # untemplated, RDKit puts the ring N at dx +12.5 on the substrate and -14.5
    # on the product, so the product reads as the substrate's MIRROR IMAGE and
    # the eye has to re-map the scaffold across the arrow. Templated on the
    # substrate it flips back to +14.5 and the two oxygens stay symmetric.
    #
    # "auto" still needs a series to read the reference off. With no products it
    # takes the FIRST HEADER STRUCTURE, which on this figure shape is the
    # substrate -- the thing every other species in the row should be drawn like.
    hdr_smiles = [it["smiles"] for it in hdr_items if it.get("smiles")]
    # ASKED FOR A TEMPLATE AND THERE IS NONE. Reported rather than raised, and
    # reported rather than passed over: "templated": true with no reference to
    # draw on is the same silence this block was written to remove, one level up.
    # It is reachable on purpose -- `template: "off"` plus a templated item is a
    # contradiction the caller wrote, not a shape the tool should refuse.
    _no_ref: list = []
    tmpl = None
    if template and str(template).lower() != "off":
        if template != "auto":
            tmpl = template
        elif flat:
            tmpl = flat[0]["smiles"]
        elif hdr_smiles:
            tmpl = hdr_smiles[0]

    # WHICH PART of that reference has to match. `auto` takes the first product,
    # and a product carries its own substituent -- so a member whose substituent
    # is not the same kind of atom cannot match the whole of it. Measured on a
    # 17-member para-substituted benzamide series: 7 placed, and the ten misses
    # were exactly the members with a non-carbon para substituent (F, Cl, Br,
    # OMe, NO2). The series' largest common core is the part they DO share:
    # 16 of 17.
    #
    # ONLY FOR `auto`, and only when the caller named no pattern of their own. An
    # explicit template is a decision; this is filling in the one the tool made.
    def _patt_for(members, ref):
        if template_patt or template != "auto" or ref is None:
            return template_patt
        # No products means the series IS the header row, which is where `ref`
        # came from in that case; handing auto_core_patt an empty member list
        # would ask it for the common core of nothing.
        smis = [c["smiles"] for c in members] or hdr_smiles
        return SF.auto_core_patt(smis, ref)

    structures = {}

    def _convert(row, sec_tmpl, sec_patt=None):
        out = []
        for it in row:
            CY.check_keys(it, _PRODUCT_KEYS, "a product cell (id %r)"
                          % (it.get("id") or "?"), _LINEAR_CELL_HINTS)
            sid = _register(structures, it["smiles"], it,
                            abbreviate=abbreviate, template=sec_tmpl,
                            template_patt=sec_patt)
            cell = {"s": sid}
            # No id from the caller -> no compound number. The id that reaches
            # the packer is the structures-dict KEY, invented when none was
            # given, and printing that would put "s1" on the page as a compound
            # number. Leaving `id` out IS the switch -- the shape a Markush cell
            # takes when its numbers live in the value table below it.
            if not str(it.get("id") or "").strip():
                cell["num"] = False
            if it.get("legend"):
                cell["legend"] = it["legend"]
            for k in ("yield", "sublabels", "label", "foot"):
                if it.get(k):
                    cell[k] = it[k]
            out.append(cell)
        return out

    sectioned = any(g["title"] is not None for g in groups)
    recipe = {"structures": structures}
    if sectioned:
        secs = []
        for g in groups:
            sec_tmpl = _section_template(g, tmpl)
            cells = [c for r in g["rows"] for c in r]
            # PER SECTION. Each section sizes and templates its own cells -- an
            # aliphatic section under an aryl reference is the case the feature
            # exists for -- so the core it shares is its own, not the table's.
            conv = _convert(cells, sec_tmpl, _patt_for(cells, sec_tmpl))
            sec = {"title": g["title"], "products": conv}
            for _k in ("style", "letter", "color", "text_color"):
                if g.get(_k) is not None:
                    sec[_k] = g[_k]
            if g["explicit"]:
                # explicit rows are taken literally, same contract as recipe["rows"]
                sec["row_sizes"] = [len(r) for r in g["rows"]]
            secs.append(sec)
        recipe["sections"] = secs
    else:
        _flat_patt = _patt_for(flat, tmpl)
        spec_rows = [_convert(r, tmpl, _flat_patt) for r in groups[0]["rows"]]
        spec_rows = [r for r in spec_rows if r]
        if not spec_rows:
            pass                      # table-only figure: no panel at all
        elif len(spec_rows) > 1:
            recipe["rows"] = spec_rows
        else:
            recipe["products"] = spec_rows[0]
    title = _title_box[0]
    if title and str(title).strip():
        recipe["title"] = str(title)
    if footnotes:
        # NO HAZARD IN WIDENING THE ANNOTATION. Widening it to `list | str` with a
        # naive list() below would shred a bare string into one footnote per
        # CHARACTER, because list("a note") is ['a', ' ', 'n', ...].
        # _footnote_lines never calls list() on a string, so widening the
        # annotation is a decision about the schema rather than a trap. It is
        # still declared `list` here; that is not changed, only safe to change.
        recipe["footnotes"] = _footnote_lines(footnotes)
    if cols:
        recipe["cols"] = int(cols)

    if header and header.get("items"):
        hitems = []
        for _hi, it in enumerate(header["items"]):
            CY.check_keys(it, _HDR_ITEM_KEYS, "header.items[%d]" % _hi,
                          _LINEAR_CELL_HINTS)
            if it.get("text") is not None:
                # A text cell in the header, which is where such cases usually
                # live: the reactant of a C-H functionalisation, a compound
                # cited by number, and the "≡" between two species.
                hitems.append({k: it[k] for k in ("text", "id", "size")
                               if it.get(k) is not None})
            elif it.get("plus"):
                hitems.append({"plus": True})
            elif "arrow" in it:
                # Membership, not truthiness: {"arrow": {}} is a bare arrow with
                # no captions, which a truthiness test would drop.
                hitems.append({"arrow": dict(it["arrow"] or {})})
            elif it.get("smiles"):
                if it.get("templated") and tmpl is None:
                    _no_ref.append(str(it.get("id") or it["smiles"])[:40])
                # THE case for the empty-id rule: two header items both written
                # "id": "" would collapse onto one entry and draw the second
                # structure twice. An empty id means "not given".
                hid = _register(
                    structures, it["smiles"], it, prefix="h",
                    abbreviate=abbreviate,
                    template=tmpl if it.get("templated") else None,
                    # The SAME pattern the products got: "templated": true means
                    # "draw this on the series' depiction", and a header product
                    # drawn on a narrower match than its own table is a header
                    # that no longer matches the table it heads.
                    template_patt=(_patt_for(flat, tmpl)
                                   if it.get("templated") else None),
                    where="header structure")
                hitems.append({"s": hid})
        if hitems:
            recipe["header"] = {"items": hitems}

    if boxes:
        recipe["boxes"] = list(boxes)
    if table:
        recipe["table"] = table
        if table_notes:
            recipe["table_notes"] = list(table_notes)
        if table_style:
            recipe["table_cfg"] = dict(table_style)

    cfg = {"footnote_width": float(footnote_width)} if footnote_width else {}
    if span:
        cfg["span"] = float(span)
    cfg |= {"label_align": label_align, "struct_valign": struct_valign,
           "last_row": last_row, "fill_width": bool(fill_width),
           "panel_pos": panel_pos, "frame": frame}
    stale = None
    if order:
        cfg["order"] = [str(b) for b in order]
        if panel_pos != "above":
            # Both given and they disagree. `order` wins because it is the more
            # specific statement, and the loser is named rather than dropped.
            stale = ("panel_pos=%r was ignored: `order` says %s."
                     % (panel_pos, cfg["order"]))

    # The Allocator is load-bearing in this engine, so READ IT. Silently throwing
    # the report away is how an overlap ships: the instrument fires and nobody
    # looks. Findings are reported alongside the file, not raised -- a figure with
    # a flagged overlap is still worth opening and looking at.
    if connectors:
        # Ids are checked HERE, where the cell ids are known. The engine asks the
        # ledger for "cell:3a" and a typo would come back as an AllocatorError
        # naming every reservation in the figure -- true, but not a sentence a
        # chemist should have to read.
        for _ci, _con in enumerate(connectors):
            CY.check_keys(_con, _GRID_CONN_KEYS, "connectors[%d]" % _ci)
            for _end in ("from", "to"):
                if str(_con.get(_end) or "") not in structures:
                    raise ValueError(
                        "connectors[%d].%s: %r is not a cell in this figure. The "
                        "ids are %s." % (_ci, _end, _con.get(_end),
                                          ", ".join(map(repr, sorted(structures)))))
        recipe["connectors"] = [{**c, "from": str(c["from"]), "to": str(c["to"])}
                                for c in connectors]
    if panels:
        _seen: set = set()
        for _pi, _p in enumerate(panels):
            CY.check_keys(_p, _PANEL_KEYS, "panels[%d]" % _pi)
            if not str(_p.get("name") or "").strip():
                raise ValueError(
                    "panels[%d]: needs a 'name' -- it is the letter drawn on the "
                    "panel AND how a later panel refers to it in 'of'" % _pi)
            if _p["name"] in _seen:
                raise ValueError(
                    "panels[%d]: two panels are both named %r, so 'of' could not "
                    "say which one it meant" % (_pi, _p["name"]))
            _seen.add(_p["name"])
            if not _p.get("items"):
                raise ValueError("panels[%d]: needs 'items' to draw" % _pi)
            for _ii, _it in enumerate(_p["items"]):
                _check_linear_item(_it, "panels[%d].items[%d]" % (_pi, _ii))
        _pspec = []
        for _p in panels:
            _items = []
            for _it in _p["items"]:
                if _it.get("text") is not None:
                    _items.append({k: _it[k] for k in ("text", "id", "size")
                                   if _it.get(k) is not None})
                elif _it.get("plus"):
                    _items.append({"plus": True})
                elif "arrow" in _it:
                    _items.append({"arrow": dict(_it["arrow"] or {})})
                elif _it.get("smiles"):
                    _items.append({"s": _register(
                        structures, _it["smiles"], _it, prefix="p",
                        abbreviate=abbreviate, where="panel structure")})
            _pspec.append({**{k: v for k, v in _p.items() if k != "items"},
                           "items": _items})
        recipe["panels"] = _pspec
    # LAST, after every _convert and after the header loop: those are what call
    # _register, so this is the first moment _DECOR is complete. Same trap as the
    # cycle's -- see the note there.
    if _DECOR:
        recipe["decorations"] = dict(_DECOR)

    from allocator import Allocator as _Alloc
    _alloc = _Alloc()
    _notes: list = []
    if _no_ref:
        _notes.append(
            "header item(s) %s asked for \"templated\": true and there is no "
            "reference to draw them on, so each was laid out on its own. "
            "`template` is %r: give it a SMILES, or add a product cell for "
            "\"auto\" to read the series off."
            % (", ".join(_no_ref), template))
    _shift: list = []
    doc = GR.render_grid(recipe, _hdr, cfg=cfg, page_mode=page,
                         alloc=_alloc, notes=_notes, shift=_shift)
    doc = _add_graphics(doc, graphics, slots, structures, abbreviate,
                        alloc=_alloc, shift=_shift)
    _path = _write_and_open(doc, filename, clean_up)
    out = ""

    if stale:
        _notes.insert(0, stale)
    if _notes:
        out += "\n\nNOTE: " + "\n      ".join(_notes)

    _lrep = _layout_report(_alloc)
    if _lrep:
        out += "\n\nLAYOUT REPORT -- " + _lrep

    # `color` ON A SECTION IS THE PILL'S FILL, so no "colour is not wired yet,
    # drawn in black" note is needed here. Accepting and reporting a key that is
    # not yet served is the right way to hold a request you cannot serve -- it
    # is discoverable, and the day it lands there is one place to delete.
    return _figure_report(_path, structures, alloc=_alloc,
                          notes=out + _chem_report())


# What a cycle node takes AT THIS LAYER -- "smiles" here, which becomes the engine's
# "s" (a key into the fragments this layer builds) a few lines into the tool.
_CYCLE_NODE_KEYS = frozenset({"smiles", "text", "label", "id", "as",
                              "highlight", "rlabels", "floating",
                              # node kind and placement
                              "kind", "at", "pairs",
                              # `content` -- a slot holding a SEQUENCE (a metal
                              # centre plus its drawn chain, a pair of products)
                              "content",
                              # `frame` -- the square-bracket form only
                              "frame",
                              # schema only for now: used -> raises
                              "graph"})
# Edge keys.
_CYCLE_EDGE_KEYS = frozenset({"as", "from", "to", "via", "ring", "name",
                              "feeders", "style", "reverse",
                              # accepted so a recipe reads whole; set -> raises
                              "crossing", "muted",
                              "bulge"})
# Ring keys.
_CYCLE_RING_KEYS = frozenset({"as", "nodes", "centre", "radius", "start_angle",
                              "direction", "label",
                              # Mirrors layout_cycle's own
                              # RING_KEYS -- FOUR whitelists, not three: this set,
                              # the engine's, this set again as the published JSON
                              # schema below, and `_cycle_ring` passing the dict
                              # through. The engine checks what these two mean.
                              "tangent_to", "at"})
_CYCLE_FEEDERS_KEYS = frozenset({"in", "out", "style", "at", "aside"})
# Path keys, the fourth table.
_CYCLE_PATH_KEYS = frozenset({"as", "edges", "name", "muted"})
# Keys that exist in the schema and are NOT implemented yet. Refused when
# actually used rather than ignored: never draw a figure that is quietly
# missing a piece the caller asked for.
_PHASE1_UNBUILT = {
    "graph": "a nested subgraph; only one level of nodes can be drawn",
}
# ONE ITEM of a node's `content`. The species vocabulary and nothing else: `at` /
# `kind` / `pairs` / `frame` / `as` describe the SLOT, and an item lives inside the
# slot rather than being one. `abbreviate` is absent for the same reason it is
# absent from a node -- it is a whole-figure argument here, unlike a feeder entry,
# which is one molecule and may differ.
_CYCLE_CONTENT_KEYS = frozenset({"smiles", "text", "label", "id",
                                 "highlight", "rlabels", "floating"})
# ...and a feeder entry, which is a species in the same sense: the node vocabulary
# plus the two keys only a feeder has, plus its own "abbreviate" (nodes take that as
# a tool argument for the whole figure; a feeder is one molecule and can differ).
_CYCLE_FEEDER_KEYS = frozenset({"smiles", "text", "label", "id", "at", "angle",
                                "abbreviate", "as",
                                "highlight", "rlabels", "floating"})


def _cycle_content_item(it, where: str, structures: dict, abbreviate: bool):
    """ONE item of a node's `content` slot, in the engine's words.

    Same job as `_cycle_feeder_entry` one level down: this layer turns a caller's
    `smiles` into a key in `structures`, and the engine sees only `s`. The
    connective passes straight through as the bare string it was -- it names no
    species, so there is nothing to build -- and the engine is what decides which
    connectives it can draw.

    Halos, R-labels and pushing arrows are per ITEM, not per slot: the drawn chain
    beside a written metal centre is the thing that gets coloured, and `_register`
    already takes the whole spec, so passing the item is the whole of it.
    """
    if isinstance(it, str):
        return it
    if not isinstance(it, dict):
        raise ValueError(
            "%s: an item of 'content' is a dict -- {\"text\": ...} or "
            "{\"smiles\": ...} -- or the bare string \"+\". Got %s."
            % (where, type(it).__name__))
    CY.check_keys(it, _CYCLE_CONTENT_KEYS, where)
    has_t, has_s = bool(it.get("text")), bool(it.get("smiles"))
    if has_t and has_s:
        raise ValueError(
            "%s: 'text' and 'smiles' are mutually exclusive -- an item is either a "
            "drawn structure or a piece of type, not both. Two things in one slot "
            "is what the sequence is for: write them as two items." % where)
    if not (has_t or has_s):
        raise ValueError(
            "%s: needs either 'smiles' (draw it) or 'text' (write it out)" % where)
    out = {}
    if it.get("label"):
        out["label"] = it["label"]
    if has_t:
        out["text"] = it["text"]
        return out
    out["s"] = _register(structures, it["smiles"], it,
                         auto=None, abbreviate=abbreviate, where="content item")
    return out


def _cycle_feeder_entry(it, where: str, structures: dict, abbreviate: bool) -> dict:
    """One `in` / `out` species, turned into something the engine can place.

    Lifted out of _cycle_step so the old `steps` form and the new `edges` form
    build a feeder entry through ONE function -- two copies would be two chances
    for the forms to disagree about what a feeder is, and that disagreement would
    show up as a figure rather than as an error.
    """
    if not isinstance(it, dict):
        raise ValueError(
            "%s: a feeder entry is a dict like {\"text\": \"H~2~O\"} or "
            "{\"smiles\": \"C#Cc1ccccc1\", \"id\": \"2a\"}, got %s"
            % (where, type(it).__name__))
    CY.check_keys(it, _CYCLE_FEEDER_KEYS, where)
    if not it.get("smiles"):
        return dict(it)
    if it.get("text"):
        raise ValueError(
            "%s: 'text' and 'smiles' are mutually exclusive -- a feeder entry is "
            "either a drawn structure or a piece of type, not both" % where)
    e = {k: v for k, v in it.items()
         if k not in ("smiles", "id", "abbreviate", "highlight", "rlabels",
                      "floating")}
    e["s"] = _register(structures, it["smiles"], it,
                       auto=where.replace("[", "_").replace("]", "").replace(
                           ".", "_").replace(" ", ""),
                       abbreviate=bool(it.get("abbreviate", abbreviate)),
                       where="%s entry" % where)
    return e


def _cycle_edge(ed, i: int, structures: dict, abbreviate: bool) -> dict:
    """One entry of the NEW `edges` table.

    Passes the edge's own fields through and turns any feeder SMILES into a
    fragment, exactly as _cycle_step does for the old form. `crossing` and `bulge`
    are built; `muted` belongs on the `paths` row and is refused here, so a recipe
    that writes it on an edge is told where it went rather than quietly drawing a
    black arrow.
    """
    where = "edges[%d]" % i
    if not isinstance(ed, dict):
        raise ValueError("%s: an edge is a dict with 'from' and 'to', got %s"
                         % (where, type(ed).__name__))
    CY.check_keys(ed, _CYCLE_EDGE_KEYS, where)
    for end in ("from", "to"):
        if not ed.get(end):
            raise ValueError(
                "%s: an edge needs '%s' -- the address of the object it starts or ends "
                "at." % (where, end))
    if "muted" in ed:
        raise ValueError(
            "%s: 'muted' is not an edge key -- it moved to the `paths` table, "
            "because greying is a property of a SET: the edge, its label and its "
            "decorations go grey together. To grey this one edge, give it an 'as' "
            "and write a one-edge path: paths=[{\"edges\": [\"<that as>\"], "
            "\"muted\": true}]." % where)
    out = {k: ed[k] for k in ("as", "from", "to", "via", "ring", "name",
                              "style", "reverse", "crossing", "bulge")
           if k in ed}
    feed = ed.get("feeders")
    if feed is not None:
        if not isinstance(feed, dict):
            raise ValueError("%s.feeders: a dict with 'in' / 'out' / 'style', got "
                             "%s" % (where, type(feed).__name__))
        CY.check_keys(feed, _CYCLE_FEEDERS_KEYS, where + ".feeders")
        nf = {k: feed[k] for k in ("style", "at", "aside") if k in feed}
        for tag in ("in", "out"):
            if feed.get(tag) is None:
                continue
            nf[tag] = [_cycle_feeder_entry(it, "%s.feeders.%s[%d]" % (where, tag, k),
                                           structures, abbreviate)
                       for k, it in enumerate(feed[tag])]
        out["feeders"] = nf
    return out


def _cycle_ring(rg, i: int) -> dict:
    """One entry of the NEW `rings` table.

    ANY number of rings, not just two -- triple-ring cycles occur in published
    figures. Each ring carries its
    own `centre`, `radius`, `start_angle`, `direction` and `label`; the engine
    validates the contents, this layer only checks the shape and the keys.
    """
    where = "rings[%d]" % i
    if not isinstance(rg, dict):
        raise ValueError("%s: a ring is a dict with 'nodes', got %s"
                         % (where, type(rg).__name__))
    CY.check_keys(rg, _CYCLE_RING_KEYS, where)
    if not rg.get("nodes"):
        raise ValueError(
            "%s: a ring needs 'nodes' -- the addresses of its members in cyclic "
            "order." % where)
    return dict(rg)


# The grey a muted path is drawn in, as [r, g, b] 0..1. Measured on published
# hand-drawn figures, not chosen: of the five greys found, 0.60 is the most common
# (13 uses) and the only one that appears on a LINE rather than as a pale table
# fill. Overridable per figure through cfg["muted_rgb"].
MUTED_RGB = (0.6, 0.6, 0.6)


def _cycle_path(pa, i: int) -> dict:
    """One entry of the `paths` table.

    A path is a NAMED SET OF EDGES -- the host of a group-level name (a TM, an
    ROH, path A / path B) and the unit of muting.
    Shape and keys here; the engine checks that every named edge exists.
    """
    where = "paths[%d]" % i
    if not isinstance(pa, dict):
        raise ValueError("%s: a path is a dict with 'edges', got %s"
                         % (where, type(pa).__name__))
    CY.check_keys(pa, _CYCLE_PATH_KEYS, where)
    if not pa.get("edges"):
        raise ValueError(
            "%s: a path needs 'edges' -- the 'as' names of the edges it groups, "
            "e.g. {\"edges\": [\"e3\", \"e4\"], \"name\": {\"above\": [\"path A\"]}}."
            % where)
    if not isinstance(pa["edges"], (list, tuple)):
        raise ValueError("%s.edges: a list of edge names, got %s"
                         % (where, type(pa["edges"]).__name__))
    return dict(pa)


def _cycle_step(st, i: int, structures: dict, abbreviate: bool) -> dict:
    """One step, with any feeder SMILES turned into a fragment the engine can place.

    The same translation the node loop above does, and it has to happen here for the
    same reason: this layer speaks "smiles" and the engine speaks "s" (a key into the
    fragments it was handed), so an entry rebuilt anywhere else would have its unknown
    keys dropped on the floor before layout_cycle could complain about them. Steps
    are not passed through untouched, so they get checked.

    Everything else about the step -- captions, style, "at", feeder_style -- is copied
    across unread and validated by the engine, which owns that vocabulary.
    """
    st = dict(st or {})
    for tag in ("in", "out"):
        if not st.get(tag):
            continue
        out = []
        for k, it in enumerate(st[tag]):
            where = "steps[%d].%s[%d]" % (i, tag, k)
            if not isinstance(it, dict):
                raise ValueError(
                    "%s: a feeder entry is a dict like {\"text\": \"H~2~O\"} or "
                    "{\"smiles\": \"C#Cc1ccccc1\", \"id\": \"2a\"}, got %s"
                    % (where, type(it).__name__))
            CY.check_keys(it, _CYCLE_FEEDER_KEYS, where)
            if not it.get("smiles"):
                out.append(dict(it))
                continue
            if it.get("text"):
                raise ValueError(
                    "%s: 'text' and 'smiles' are mutually exclusive -- a feeder "
                    "entry is either a drawn structure or a piece of type, not both"
                    % where)
            # The keys THIS LAYER consumes, dropped before the engine sees them.
            # highlight / rlabels / floating are spent by _register just below,
            # exactly as "smiles" is -- and leaving them in is a real bug, not a
            # tidiness point: _CYCLE_FEEDER_KEYS accepts "highlight" here and the
            # engine's FEEDER_KEYS does not, so a feeder that asked for one would
            # raise `steps[0].in[0]: unknown key 'highlight'` from layout_cycle.
            # Both checks are right; the seam between them is what this closes.
            e = {k2: v for k2, v in it.items()
                 if k2 not in ("smiles", "id", "abbreviate",
                               "highlight", "rlabels", "floating")}
            # The structures dict is shared with the nodes, so a reagent drawn as a
            # feeder here and as a species elsewhere is built once under one id --
            # which is exactly the legitimate reuse _register stays quiet about
            # when the SMILES agree, and reports when they do not.
            e["s"] = _register(
                structures, it["smiles"], it,
                auto="%s%d_%d" % (tag, i, k),
                abbreviate=bool(it.get("abbreviate", abbreviate)),
                where="%s entry" % tag)
            out.append(e)
        st[tag] = out
    return st


@mcp.tool()
@_cx_doc
def draw_catalytic_cycle(nodes: list, steps: list = None, filename: str = "cycle",
                         clean_up: bool = False, page: str = "fit",
                         # A ONE-LINE CENTRE NOTE MAY BE WRITTEN BARE, same shape
                         # as draw_radial_figure's `center_label: str | list`. The
                         # body takes either (`list(center) if
                         # isinstance(center, (list, tuple)) else [center]`), and
                         # the string form would otherwise only work for a caller
                         # that calls the function directly and never meets
                         # pydantic.
                         center: str | list | None = None, radius: float = None,
                         direction: str = "cw", start_angle: float = 270.0,
                         label_pos: str = "auto", arc_style: dict = None,
                         abbreviate: bool = False, graphics: list = None,
                         slots: list = None,
                         edges: list = None, rings: list = None,
                         paths: list = None) -> dict:
    """Draw a CATALYTIC CYCLE / mechanism ring as a CDXML file.

    Species sit on a circle and the reaction arrows ARE the ring.

    nodes: the species IN RING ORDER, starting at the top and running clockwise
      by default. Each is EITHER a structure:
        {"smiles": "...", "id": "A"?, "label"?: "L~n~Pd(0)",
         "rlabels"?: {..}, "floating"?: [..], "highlight"?: [...]}
      "highlight" is the whole spec draw_molecule takes -- "color", "fill",
      "squiggle", the five bond styles, "halo" / "dot" / "locants". A feeder entry
      (see below) takes it too.
      "label" is the text drawn with the species -- normally the formula, e.g.
      "L~n~Pd(0)" or "Ar-Pd(II)-X". It falls back to "id" if omitted.
      OR a species with no structure to draw, written as type:
        {"text": "Cu(I)", "id": "Cu"?}
      Use "text" for a species you would write rather than draw -- "Cu(I)",
      "L~n~Pd(0)", "Ar-Pd(II)-X". That text IS the node; no structure is drawn and
      no label is added unless you pass "label" explicitly, so the name appears
      exactly once. Do NOT fake this with {"smiles": "[Cu]", "label": "Cu(I)"} --
      that draws a stray lone Cu atom next to the words.
      "id" is a key only, never drawn.

    steps: one entry per arrow, in the same order, len(nodes) to close the ring
      (the usual case) or len(nodes)-1 to leave it open. Each is:
        {"above"?: ["oxidative addition", "+ ArX"],
         "inside"?: true, "reverse"?: true, "style"?: {...}}
      Species entering or leaving a step can simply be caption lines -- "+ ArX",
      "- HBr". Captions sit OUTSIDE the ring; "inside": true puts one in the
      middle instead.

    FEEDER ARROWS, per step -- the literature's other way of showing what enters
    and leaves: a small arc TANGENT to the ring arc, with the reagent/by-product
    written at its free end.
        {"in":    [{"text": "PhC≡CH (#2a#)"}],
         "out":   [{"text": "H~2~O"}],
         "aside": ["*t*-BuOK"],
         "at": 0.5,
         "feeder_style": "split"}
      feeder_style: how the pair meets the ring.
        "split" (default) is two arcs sharing one tangency point, so they read as
                one continuous S through the ring: the reagent arrives from
                outside headless, the by-product peels away into the middle of the
                cycle carrying the only arrowhead. Any number of in/out entries;
                successive pairs fan along the arc.
        "single" is one arc outside the ring, touching it at the middle of its
                sweep, tail = reactant and head = product, both labels outside.
      in / out: lists. Each entry is a species, WRITTEN or DRAWN:
        {"text": "PhC≡CH (#2a#)"}                       -- type at the arc's free end
        {"smiles": "C#Cc1ccccc1", "id": "2a",
         "label"?: "#2a#", "abbreviate"?: true}         -- the molecule itself
      A drawn feeder gets its "id" under it as a label unless you give a "label" of
      your own; the structures are shared with the nodes, so the same "id" means
      the same molecule built once.
      Plus two optional placement keys:
        "at":   where on the ring arc this one touches, strictly between 0 and 1,
                same scale as the step's "at" below. Overrides it -- use it to pull
                an in/out pair apart so they no longer share a point.
        "angle": degrees to swing this entry away from the arc's own tangent, for a
                nudge you can see and the geometry cannot. It moves WHERE the species
                sits; a drawn structure is never rotated by it, or by anything else.
      Which end carries the arrowhead is NOT a per-entry choice: it follows from
      feeder_style.
      aside: the base or additive, e.g. "*t*-BuOK", set beside the tangency point:
        outside the ring under "split",
        inside it under "single". Type only -- it is a reagent note, not a species.
      at: where along the step's arrow the feeder meets it, strictly between 0 and 1.
        DEFAULT 0.5 = THE GEOMETRIC MIDPOINT OF THE ARROW: a feeder joins its step
        halfway along, and moving it off centre is a special case, not a knob to
        reach for. Reach for it when something about that particular step asks for
        it -- pulling an in/out pair onto two separate meetings with the ring, or
        steering the motif away from a crowded quarter of the figure -- and leave it
        alone otherwise.
        Measured along the arrow AS DRAWN: 0 and 1 are the two ends of the visible
        arc, where it starts and stops clear of the two species, NOT the species
        centres, so 0.5 is the midpoint of the ink. It is the middle of THIS step's
        own arrow, not of a 360/n slice.
      Omit all of these and the step has no feeder.

    ----------------------------------------------------------------------------
    THE EDGE TABLE -- `edges` (+ optional `rings`), an ALTERNATIVE to `steps`.

    `steps` says "one arrow per gap, in node order", so it cannot express an arrow
    that skips a node, or one that leaves the circle. `edges` states each arrow's
    two ends by name, which is the whole difference:

        nodes=[{"text": "L~n~Pd(0)", "as": "Pd0"},
               {"text": "Ar-Pd^II^-X", "as": "ArPdX"},
               {"smiles": "...", "as": "prod", "at": [320, 40], "label": "3a"}],
        edges=[{"from": "Pd0",   "to": "ArPdX", "name": {"above": ["ArX"]}},
               {"from": "ArPdX", "to": "Pd0"},
               {"from": "ArPdX", "to": "prod", "name": {"below": ["- HX"]}}],
        rings=[{"nodes": ["Pd0", "ArPdX"]}]

    `from` / `to` name nodes by their `as`. Give one, or the engine names them
    `node#0`, `node#1`, ... An edge also takes `via`, `name` (`above` / `below`),
    `feeders` (`in` / `out` / `style` / `at` / `aside`), `style`, `reverse`, `as`.

    `rings` is optional: leave it out and the ring is every node WITHOUT an `at`,
    in the order written. A ring takes `nodes` (cyclic order), `as`, `centre`,
    `radius`, `start_angle`, `direction`, `label`.

    `at` ON A NODE puts it at those coordinates instead of on the circle -- the
    species drawn inside the ring and the one hanging off it are the same thing
    here, at different coordinates. A pinned node spends no circumference, so it
    does not push the radius out.

    THE SHAPE OF AN EDGE FOLLOWS FROM THE DATA, it is not declared: two adjacent
    ring members is the ring's ARC, two non-adjacent members is a CHORD, and an
    edge with an end off the ring is an OFF-RING edge. Arcs are drawn as arcs;
    chords and off-ring edges are drawn as straight arrows. Say `via: "line"` to
    draw a ring step straight as well.

    `steps` and `edges` are mutually exclusive.

    RETURNS a dict:
        {"path": "<download link to cycle.cdxml>",
         "warnings": [...],          # EVERYTHING wrong with this figure: see below
         "feeder_extents": {...},    # degrees each motif subtends at the centre,
                                     #   or null -- see feeder_extent_undefined
         "feeder_extent_undefined": [...],   # motifs with no angular width to quote
         "ring": {...},              # how the circle was sized: see below
         "layout_report": "..."}     # the full layout findings, when there are any
      "warnings" is the one dial to read. It carries every problem the layout found
      -- a caption that could not clear a feeder, any known overlap -- so it is
      never [] on a figure with a known overlap in it. "layout_report" is the same
      findings in full, grouped, for when a warning is not enough to find the spot.
      "feeder_extents" gives each motif's angular width at the ring centre, and null
      for a motif holding a box that COVERS that centre; those motifs are named in
      "feeder_extent_undefined".
      "ring" is why the circle came out the size it did: "radius", "gaps" (the
      degrees each step was allotted) and "need" (of 360).
      A feeder's step is given room so the feeder does not run into the species
      beside it or into the centre note. What this cannot fix is a
      motif whose OWN two labels are too long for one arc to carry (a 30-character
      reagent and a 30-character by-product on the same step), so still read
      "warnings" and "layout_report" rather than assuming it came out clean.
      LONG LINES FOLD. Any line of type wider than 10 bond lengths (144 units,
      roughly 32 characters at the caption size) is broken at a word boundary into
      as many lines as it needs. Nothing in the recipe
      turns this on or off; give a feeder an "at" of its own, or shorten the name, if
      a folded pair still reports an overlap.

    ARROW STYLES, per step in "style" (or for every step at once via arc_style):
        {"head": "Full", "tail": "Full"}   heads on BOTH ends: a reversible step
        {"head": null}                     plain headless arc
        {"head": "HalfLeft"} / "HalfRight" single-barb head, for one-electron steps
        {"line": "dashed"}                 dashed shaft, e.g. a proposed step
        {"head_type": "Solid"|"Hollow", "head_size": n, "head_width": n}
      Plus "reverse": true on the step to run that arrow AGAINST the ring, which
      swaps the heads and leaves the geometry alone.

    center: text for the middle of the ring, one string per line -- usually the
      net transformation. The ring radius grows to keep it clear of the species,
      so a long net reaction pushes the circle out rather than crowding it.

    label_pos: "auto" (default) puts each species label BELOW its structure, as
      the other engines do, and falls back to radially outward where below is
      occupied. "below" or "outward" force one or the other.
    direction: "cw" (default) | "ccw".
    start_angle: 270 puts the first species at the TOP. 0 = right, 90 = bottom,
      180 = left (CDXML y grows down).
    radius: force the ring radius; omitted, it is derived from the species sizes,
      the neighbour spacing and the centre note.
    page: 'fit' (default -- a ring wants a tight page) | 'letter' | 'double' |
      'single' | 'auto'.
    Text markers: ~sub~ ^sup^ *italic* #bold# _underline_.
    """
    if not nodes:
        # raises rather than returning a message: this tool now returns a dict,
        # and a string that looks like a result is worse than a clear failure
        raise ValueError("nodes is empty -- a cycle needs at least one species")

    # Validate rather than skip. Dropping any node without a "smiles" key would
    # lose the species AND silently misassign every caption after it:
    # `steps` is passed through untouched and the engine pairs steps[i] with node
    # i, so removing one node shifts the whole ring's captions by one.
    # REFUSE A BAD SPEC BEFORE ANYTHING IS BUILT. First statement in the
    # tool, ahead of _new_batch, so a rejected call leaves no per-batch state
    # behind either. See preflight.py for what it checks and why the engines'
    # own key whitelists are not enough.
    PF.preflight_cycle({"nodes": nodes, "steps": steps, "edges": edges,
                        "rings": rings, "arc_style": arc_style, "page": page,
                        "direction": direction, "label_pos": label_pos})
    _new_batch()
    # BEFORE the first _register (see the note in draw_scope_table). TWO sources:
    # the ring's own nodes, and every feeder entry hanging off a step -- those go
    # through _cycle_step, which registers them further down and so
    # is also downstream of this.
    # THE MUTED GREY, and it rides the same mechanism. The ACS eight hold no grey
    # at all (white black red yellow green cyan blue magenta), so a greyed path
    # needs a colour appended to this document's colortable -- which is exactly
    # what _extend_palette does for a highlight's [r, g, b].
    #
    # 0.60 is NOT a taste call: it is the grey published figures already use. A
    # sweep of hand-drawn ChemDraw files finds five distinct greys, and 0.60 is both
    # the most common (13 uses) and the only one used on LINES rather than as a
    # table-row tint (it appears on a <graphic GraphicType="Line">);
    # the others sit at 0.85-0.91, which is fill territory.
    #
    # CONDITIONAL. Appended only when some path is actually muted, so a figure with
    # no muting keeps the eight-entry colortable it has always had -- which is what
    # keeps the byte-identity regression gate green.
    _mute_spec = [{"color": list(MUTED_RGB)}] if any(
        isinstance(x, dict) and x.get("muted") for x in (paths or [])) else []
    _hdr = _extend_palette(C.acs_header_prefix(), (
        [nd.get("highlight") for nd in (nodes or [])
         if isinstance(nd, dict) and nd.get("highlight")]
        # ...and a highlight on an ITEM of a `content` slot. A halo lives on
        # the species, and inside a slot the species are the items -- the chain
        # beside a written metal centre is the thing that gets coloured. Left out,
        # its colour would simply not be in the colortable and the halo would come
        # out as whatever index 0-7 happened to sit there.
        + [it.get("highlight") for nd in (nodes or [])
           if isinstance(nd, dict) and isinstance(nd.get("content"), list)
           for it in nd["content"]
           if isinstance(it, dict) and it.get("highlight")]
        + [e.get("highlight") for st in (steps or []) if isinstance(st, dict)
           for side in ("in", "out") for e in (st.get(side) or [])
           if isinstance(e, dict) and e.get("highlight")]
        + [_mute_spec]))
    structures, spec_nodes = {}, []
    for i, nd in enumerate(nodes):
        if not isinstance(nd, dict):
            raise ValueError(
                "nodes[%d]: a node is a dict like {\"smiles\": \"...\"} or "
                "{\"text\": \"Cu(I)\"}, got %s" % (i, type(nd).__name__))
        # Strict, and it has to happen HERE: this loop rebuilds each node into the
        # engine's own vocabulary, so an unknown key would be dropped on the floor
        # before layout_cycle could complain about it. `steps` is passed through
        # untouched, so the engine checks that level itself.
        CY.check_keys(nd, _CYCLE_NODE_KEYS, "nodes[%d]" % i)
        for _k, _what in _PHASE1_UNBUILT.items():
            if nd.get(_k) is not None:
                raise ValueError(
                    "nodes[%d]: '%s' is not supported -- it is %s. Refused rather than "
                    "ignored, so the figure you get is never quietly missing what you asked "
                    "for."
                    % (i, _k, _what))
        # `kind` may be stated. Inferred from smiles / text when
        # it is not, so the same `nodes` list serves both input forms -- see the
        # decision list. A stated kind must AGREE, because a disagreement means one
        # of the two is a typo and guessing which would be guessing.
        _kind = nd.get("kind")
        if _kind is not None and _kind not in ("structure", "text", "marker"):
            raise ValueError(
                "nodes[%d]: kind is 'structure', 'text' or 'marker' (a nested 'subgraph' "
                "is not supported). Got %r."
                % (i, _kind))
        has_t, has_s = bool(nd.get("text")), bool(nd.get("smiles"))
        _content = nd.get("content")
        if _content is not None:
            if has_t or has_s:
                raise ValueError(
                    "nodes[%d]: 'content' IS the slot's contents, so it replaces "
                    "'smiles' / 'text' rather than sitting beside them. Put the %s "
                    "into the list: {\"content\": [{\"%s\": ...}, ...]}."
                    % (i, "structure" if has_s else "type",
                       "smiles" if has_s else "text"))
            # STATED KIND REFUSED, and this is the one place in the node loop where
            # that is the right answer. The schema splits keys into "falls out of
            # the data" and "must be declared": for a slot holding a sequence, the
            # structure/text question falls out (is anything in it drawn?), so a
            # stated word could only agree or be a typo. And `marker` -- the word
            # that genuinely cannot be inferred -- has no meaning here: the pairing
            # rule turns on one mark with edges continuing THROUGH it, and a
            # sequence is not one mark.
            if _kind is not None:
                raise ValueError(
                    "nodes[%d]: a 'content' slot does not take a 'kind'. It is a structure "
                    "slot if anything in the sequence is drawn and a text slot otherwise, "
                    "which the sequence already says. ('marker' is not available either: a "
                    "marker is ONE mark that edges pair through, not a row of things.)" % i)
            # A non-list, or an empty one, goes THROUGH untouched so that the
            # engine's `_content_layout` is the single place that complains about
            # the shape of the sequence -- two spellings of that message is two
            # chances for them to drift.
            if isinstance(_content, list) and _content:
                _content = [
                    _cycle_content_item(it, "nodes[%d].content[%d]" % (i, j),
                                        structures, abbreviate)
                    for j, it in enumerate(_content)]
        if _kind == "structure" and not has_s:
            raise ValueError("nodes[%d]: kind 'structure' needs a 'smiles'" % i)
        if _kind in ("text", "marker") and not has_t:
            raise ValueError(
                "nodes[%d]: kind %r is written out, so it needs a 'text' -- the SET box, "
                "the e- circle, the kcal/mol cloud."
                % (i, _kind))
        if nd.get("pairs") is not None and _kind != "marker":
            raise ValueError(
                "nodes[%d]: 'pairs' declares which edges continue into which "
                "THROUGH A MARKER, so it belongs on a 'marker' node; this one is "
                "%r." % (i, _kind or ("structure" if has_s else "text")))
        if has_t and has_s:
            raise ValueError(
                "nodes[%d]: 'text' and 'smiles' are mutually exclusive -- a node "
                "is either a drawn structure or a piece of type, not both" % i)
        if not (has_t or has_s or _content is not None):
            raise ValueError(
                "nodes[%d]: needs either 'smiles' (draw a structure) or 'text' "
                "(write the species out, e.g. {\"text\": \"Cu(I)\"}), or "
                "'content' for a slot holding several things" % i)
        auto = "n" + str(i + 1)
        if _content is not None:
            # A content slot's own id still has to exist -- it is the allocator owner
            # and the report's name for this node -- and there is no single
            # structures key to take it from, so it works exactly as a text node's
            # does: what was given, or the running number.
            cell = {"content": _content,
                    "id": str(nd.get("id") or "").strip() or auto}
        elif has_t:
            # A text node holds no structure, so it is not registered -- but an
            # empty id has to mean the same thing on both branches.
            cell = {"text": nd["text"],
                    "id": str(nd.get("id") or "").strip() or auto}
        else:
            cell = {"s": _register(structures, nd["smiles"], nd, auto=auto,
                                   abbreviate=abbreviate, where="node")}
            # The key is still `auto` ("n1") when no id was given -- the engine
            # needs a name -- but it is ours, not a compound number, so say
            # outright that nothing is printed under the structure.
            if not str(nd.get("id") or "").strip():
                cell["label"] = ""
        if nd.get("label"):
            cell["label"] = nd["label"]
        if nd.get("as"):
            # The ADDRESS name, passed through untouched. Deliberately NOT folded
            # into `id`: "id" is spent here on the structures dict (chemical
            # identity, shared on purpose), and an address is per instance. Two
            # namespaces, neither giving anything up to the other -- see
            # layout_cycle._as_id.
            cell["as"] = nd["as"]
        if _kind:
            # Passed through so the engine sees the caller's own word. Without it a
            # marker would be indistinguishable from a text node, and the pairing
            # rule keys off exactly that distinction.
            cell["kind"] = _kind
        if nd.get("pairs") is not None:
            cell["pairs"] = nd["pairs"]
        if nd.get("at") is not None:
            # Explicit coordinates. "Inside the ring / outside the
            # ring" is not a property of a node -- an `at` node is placed where the
            # caller says and a node without one is placed by its ring, and the
            # centre-of-the-circle species and the off-ring chain of a figure are the
            # SAME thing in the model at different coordinates.
            cell["at"] = nd["at"]
        if nd.get("frame") is not None:
            # `frame`. Passed through; the engine validates the
            # shape, because it is the engine that knows which forms it can draw.
            cell["frame"] = nd["frame"]
        spec_nodes.append(cell)

    recipe = {"structures": structures, "nodes": spec_nodes}
    # `edges` (+ optional `rings`) is the NEW input form. `steps` is the old one.
    # Exactly one of them, and old-form recipes take the identical path they
    # always did -- which is what keeps the byte-identity regression gate green
    # (every recipe written in the old form must still produce identical bytes).
    if edges is not None and steps:
        raise ValueError(
            "'steps' and 'edges' are two alternative input forms; give one, not "
            "both. `edges` states each arrow's two endpoints, which `steps` cannot.")
    if edges is not None:
        recipe["edges"] = [_cycle_edge(e, i, structures, abbreviate)
                           for i, e in enumerate(edges)]
        if rings is not None:
            recipe["rings"] = [_cycle_ring(r, i) for i, r in enumerate(rings)]
        if paths is not None:
            recipe["paths"] = [_cycle_path(x, i) for i, x in enumerate(paths)]
    elif rings is not None:
        raise ValueError(
            "'rings' describes the circles the `edges` table's arcs belong to, so "
            "it needs `edges`. The old `steps` form has one implicit ring and "
            "takes its geometry from `radius` / `start_angle` / `direction`.")
    elif paths is not None:
        raise ValueError(
            "'paths' groups edges by name, so it needs `edges`. The old `steps` form "
            "has no way to name a subset of its arrows -- that is what `paths` is "
            "for.")
    if steps:
        recipe["steps"] = [_cycle_step(s, i, structures, abbreviate)
                           for i, s in enumerate(steps)]
    # AFTER the steps, not before: _cycle_step registers each feeder's structure,
    # so a feeder's halo does not exist in _DECOR until that loop has run. Copying
    # it early would silently drop exactly the decorations that FEEDER_KEYS goes out
    # of its way to accept.
    if _DECOR:
        recipe["decorations"] = dict(_DECOR)
    # AFTER the steps for the same reason `decorations` is: _cycle_step registers
    # each feeder's structure, so a feeder has no _ANCHOR_SRC entry until that loop
    # has run, and a feeder's atoms would be unaddressable.
    #
    # UNCONDITIONAL, unlike layout_linear, which sets this key only when the caller
    # passed `annotations`. Two reasons it can be: it is a shallow copy of
    # references, so it costs nothing; and the expensive half of an atom lookup
    # (node_points + RDKit) is deferred inside address.py rather than here, which
    # is where the cost actually is. Making the CHEAP half conditional is how
    # linear ends up unable to address a reagent at all.
    recipe["anchor_src"] = dict(_ANCHOR_SRC)
    if center:
        recipe["center"] = list(center) if isinstance(center, (list, tuple)) else [center]
    if radius:
        recipe["radius"] = float(radius)

    cfg = {"direction": direction, "start_angle": float(start_angle),
           "label_pos": label_pos}
    if arc_style:
        cfg["arc_style"] = dict(arc_style)
    # The colortable INDEX the grey landed on, not the rgb: _extend_palette
    # rewrote the spec in place when it appended it. Set only when something is
    # muted, so `cfg` is byte-for-byte what it always was otherwise.
    if _mute_spec:
        cfg["muted_color"] = _mute_spec[0]["color"]

    from allocator import Allocator as _Alloc
    _alloc = _Alloc()
    _rep: dict = {}
    _shift: list = []
    doc = CY.render_cycle(recipe, _hdr, cfg=cfg,
                          page_mode=page, alloc=_alloc, report=_rep,
                          shift=_shift)
    doc = _add_graphics(doc, graphics, slots, structures, abbreviate,
                        alloc=_alloc, shift=_shift)

    _path = _write_and_open(doc, filename, clean_up)
    # ONCE. _figure_report re-reads the written file and runs the InChI round
    # trip, so calling it twice would do that work twice.
    _base = json.loads(_figure_report(_path, structures, alloc=_alloc,
                                      notes=_chem_report()))
    out = {**_base,
           "engine": engine_stamp(),
           # CONCATENATED, not replaced: _figure_report puts a structural failure
           # on this same key, and assigning here would drop it on the floor.
           "warnings": (_base.get("warnings") or []) + (_rep.get("warnings") or []),
           "feeder_extents": _rep.get("feeder_extents") or {},
           "feeder_extent_undefined": _rep.get("feeder_extent_undefined") or [],
           "ring": _rep.get("ring") or {},
           # one entry per path.
           "paths": _rep.get("paths") or {},
           # one entry per ring. `ring` stays ring 0 so no
           # existing caller moves.
           "rings": _rep.get("rings") or {}}
    # Read the instrument. Throwing the report away is how an overlap ships.
    _lrep = _layout_report(_alloc)
    if _lrep:
        out["layout_report"] = _lrep
    chem = _chem_report()
    if chem:
        out["structure_report"] = chem.strip()
    return out


def _gap_value(name: str, v, *, keywords: frozenset = frozenset()):
    """One gap argument, coerced to a positive float -- or a keyword left as itself.

    THE BUG THIS CLOSES. `row_gap: float | str` is a pydantic SMART union: a JSON
    string matches the `str` member EXACTLY, so no float coercion is ever attempted,
    and row_gap=200 sent over the wire as "200" arrives here as the string '200' and
    would be refused as "not a positive number". Widening the annotation so "auto"
    can reach the engine admits every other string with it, and nothing narrows it
    back. A plain `float` annotation has the OTHER half of the same problem: "auto"
    is rejected by pydantic with "unable to parse string as a number" BEFORE this
    file is entered, where no message from this file can say what the forms are.

    So both annotations are deliberately loose (`float | str`) and THIS is the
    gate. A string that spells a number is a number -- there is no reading of "200"
    that means anything else -- and anything that is neither that nor a listed keyword
    is refused here, in a message that names both accepted forms.

    `keywords` is per-parameter rather than a shared constant because the parameters
    do not accept the same ones: what "auto" means is a property of the quantity being
    solved, not of gaps in general.

    A KEYWORD IS SPELLED EXACTLY. A `.strip()` first would accept `" auto "` here
    and nowhere else -- a direct call takes it and a validating client refuses it,
    which is the worst possible split: the form that works in a test is the one
    that fails in the field. Not stripping makes this gate agree with
    what the schema publishes. It does NOT make `" 200 "` an error -- `float()` skips
    surrounding whitespace on its own, and that has never been this function's doing.
    """
    if v is None:
        return None
    forms = "a positive number" + "".join(
        ' or "%s"' % k for k in sorted(keywords))
    if isinstance(v, str):
        if v in keywords:
            return v
        try:
            v = float(v)
        except ValueError:
            raise ValueError("%s must be %s (got %r)" % (name, forms, v)) from None
    # bool BEFORE the number test: True is an int in Python and would otherwise sail
    # through as a 1.0-unit gap.
    if isinstance(v, bool) or not isinstance(v, (int, float)) or v <= 0:
        raise ValueError("%s must be %s (got %r)" % (name, forms, v))
    return float(v)


@mcp.tool()
@_cx_doc
def draw_radial_figure(center_smiles: str, products: dict = None, box: bool = False,
                       center_label: str | list = "", filename: str = "radial",
                       clean_up: bool = False, page: str = "letter",
                       abbreviate: bool = False, highlight: list = None,
                       graphics: list = None, slots: list = None,
                       # The same branch draw_linear_scheme has, so the same
                       # widening -- these two share FOOTNOTE_DEFAULTS and are
                       # meant to be one object (a regression test holds them
                       # together).
                       footnotes: list | str | None = None, trunk: bool = False,
                       # BOTH SAY EXACTLY THE TWO FORMS, and _gap_value is still
                       # the gate. The two obvious annotations are each wrong in
                       # their own direction: plain `float` makes pydantic refuse
                       # "auto" before this function is entered, with a message
                       # that cannot say what the accepted forms are, and
                       # `float | str` is a SMART union that matches the `str`
                       # member EXACTLY, so row_gap sent over the wire as "200"
                       # arrives as the STRING '200' and every other string arrives
                       # with it. Naming the keyword as a Literal fixes both at
                       # once: "200" no longer matches the string member, so it is
                       # coerced to 200.0; "auto" matches; "wide" matches neither
                       # and pydantic's own error names both forms.
                       #
                       # _gap_value is NOT redundant after this. It is what refuses
                       # a negative, a zero and a bool -- pydantic passes -5 and
                       # True straight through as floats -- and what says so in one
                       # message that names the keyword too.
                       col_gap: float | Literal["auto"] | None = None,
                       row_gap: float | Literal["auto"] | None = None,
                       rotate_captions: bool = False,
                       spokes: dict = None,
                       footnote_width: float | None = None) -> str:
    """Draw a RADIAL (one center -> many products) figure as a CDXML file.

    products: {slot: {"smiles": .., "id": .., "yield"?: .., "foot"?: ..,
                      "above"?: [..], "below"?: [..],
                      "rlabels"?: {..}, "floating"?: [..], "highlight"?: [...]}}
      where slot is one of TL, T, TR, L, R, BL, B, BR (the 3x3 grid around the
      center). "highlight" is the same spec draw_molecule takes, all of it -- the
      figure-level `highlight` argument decorates the CENTRE, and a product carries
      its own.

    above / below: reaction conditions on the arrow that POINTS AT this product,
      one list entry per line ("above": ["Pd(OAc)~2~ (5 mol%)", "K~2~CO~3~"]).
      They are the arrow's caption, not the product's.
      On a diagonal spoke "above" means the UPPER side of the shaft (the two
      vertical slots T and B have no upper side, so there it means the right).
      Text is not rotated unless rotate_captions=True (below).

    footnotes: lines set under the whole figure, left-aligned, in 8 pt. This is
      the other half of the lettered-caption idiom -- write "(a)" on each arrow
      and spell the conditions out down here.

      FOLDED TO THE FIGURE'S OWN WIDTH, so a long note comes back as several lines
      instead of running off the side. ~sub~ ^sup^ *it* #bold# spans are kept whole
      across the break. See footnote_width for the one case where a caller knows
      better: this drawing is PART of a bigger figure.
    footnote_width: WHICH FIGURE the footnotes run under, in points. Omitted, they
      fold to THIS figure's own width -- what "running prose under this figure"
      means, and right for every ordinary call. Pass a number only when this
      drawing is one PART of a bigger one (set beside another tool's output),
      because then the figure the prose runs under is the composed one.

      THE REPLY SAYS WHAT HAPPENED, in a "footnotes" field beside "gaps": `width`
      the cap that was used, `lines` how many lines each note folded into, and
      `over` the lines still too wide AFTER folding. `over` is normally empty and is
      not empty when a single word is wider than the whole figure -- such a word
      gets its own line rather than being hyphenated, so the figure is
      honestly too wide and the reply says which note did it.

    col_gap / row_gap: distance from the centre to a product column / row, in page
      units. BOTH ARE SOLVED FROM THE FIGURE BY DEFAULT -- pass nothing and each
      gap is derived from the structures you actually sent. Pass a number to pin
      one; pass "auto" to ask for the solved value explicitly. The reply's "gaps"
      field says what was chosen and how ("solved" / "pinned" / "default").

    rotate_captions=True: set the conditions ALONG a diagonal spoke instead of
      horizontally beside it, which is what published radial figures do. OFF by
      default.

      DIAGONALS ONLY. The four corner slots turn; T, B, L and R keep their captions
      horizontal, and so does every slot on a trunk (a trunk turns its diagonals into
      horizontal arrows).

      IT CHANGES THE ROW PITCH: the solved row_gap is re-derived for it and the
      reply's `gaps` field reports the value actually used. Turning this on and
      pinning the old row_gap is how you get a collision.

      Type always reads left-to-right, so on a left-hand spoke it runs back toward
      the hub. That is the convention, not a bug.

    trunk=True: the slots down one side (TL/L/BL, and TR/R/BR) stop radiating
      separately -- one stub leaves the hub, meets a vertical line, and each
      product hangs off that line on a horizontal arrow. This is the common
      literature drawing for a diversification figure. It applies per side and
      only where that side holds TWO OR MORE products; a lone slot, and T and B
      (which sit on no side), keep their own spoke. Products do not move: the
      3x3 grid is untouched and only the arrows change shape. Default False.

    box=true draws a box around the center, enclosing the structure and its label.

    center_label: the centre compound's own label, set UNDER the structure and
      centred -- the same place a product's "id, yield" goes. One string, or a list
      for a two- or three-line label:
      ["#2a#", "70%", "from 2.0 mmol of #1a#"]. Markup is the usual ~sub~ / ^sup^ /
      #bold#, so a counterion still reads "BF~4~^-^".
    abbreviate=True collapses recognised pendant groups to contracted superatom labels.

    spokes: an ALIAS for `products`, and the better name once a figure has mixed
      directions. Identical in every other way; pass one or the other, never both.

    direction (per slot, inside products/spokes): "out" (the default) points the
      arrow FROM the hub AT that slot; "in" points it from that slot AT the hub --
      the reaction's entry.

      PER SLOT, NOT PER FIGURE: the figures this exists for are "one arrow in,
      several out".

      An inbound slot is never folded into a trunk.

      IT MOVES ONLY THAT SPOKE: its arrow and that arrow's caption shift 2.0 units
      toward the hub; every other spoke, structure, label and the page stay put."""
    # REFUSE A BAD SPEC BEFORE ANYTHING IS BUILT. First statement in the
    # tool, ahead of _new_batch, so a rejected call leaves no per-batch state
    # behind either. See preflight.py for what it checks and why the engines'
    # own key whitelists are not enough.
    PF.preflight_radial({"products": products, "center_smiles": center_smiles,
                         "page": page, "highlight": highlight})
    _new_batch()
    # ONE OF THE TWO, and saying which. Silently preferring one when both are given
    # would draw a figure the caller did not ask for and there is nothing in the file
    # to notice it by.
    if products is not None and spokes is not None:
        raise ValueError(
            "pass products= or spokes=, not both -- they are the same argument "
            "under two names (spokes reads better when one of the slots is the "
            "reaction's entry rather than a product).")
    products = products if products is not None else spokes
    if products is None:
        raise ValueError("draw_radial_figure needs products= (or spokes=): a dict "
                         "keyed by slot, one of TL T TR L R BL B BR.")
    # THROUGH THE REGISTRAR, like every other tool. Bypassing it would cost more
    # than tidiness: no duplicate-id check, no _BUILT_FROM, no anchor source -- on
    # the one tool whose products are keyed by slot, so two products both labelled
    # "2" would draw the same structure twice in silence, which is precisely the
    # defect the duplicate-id check exists to report.
    #
    # THE ENGINE IS NOT TOUCHED. layout_radial still receives {"frag": (frag, bbox)}
    # per slot; `structures` is built alongside and the entry is handed over from it.
    # Unifying that shape is a change to the engine plus four of its tests, it buys
    # only consistency, and it is deliberately NOT done here.
    structures: dict = {}
    # BEFORE the first _register (see the note in draw_scope_table). TWO sources:
    # each product's own spec, and the figure-level `highlight`, which decorates
    # the CENTRE. The centre's spec dict is built fresh below as
    # {"highlight": highlight}, but `highlight` is the caller's own list and
    # _extend_palette rewrites the dicts INSIDE it, so the rewrite still lands.
    _hdr = _extend_palette(C.acs_header_prefix(), (
        [p.get("highlight") for p in (products or {}).values()
         if isinstance(p, dict) and p.get("highlight")]
        + ([highlight] if highlight else [])))
    # slot -> decoration, which is how layout_radial keys them (a radial species has
    # no name but its slot). _DECOR is keyed by structure id, so this maps across.
    _rdeco: dict = {}

    def _rd(slot_key, smiles, spec, where, auto):
        sid = _register(structures, smiles, spec, auto=auto,
                        abbreviate=abbreviate, where=where)
        if sid in _DECOR:
            _rdeco[slot_key] = _DECOR[sid]
        return structures[sid]

    center = {"frag": _rd("center", center_smiles, {"highlight": highlight},
                          "centre", "center"),
              "box": box}
    if center_label:
        center["label"] = center_label
    prod: dict = {}
    for slot, d in products.items():
        CY.check_keys(d, _RADIAL_PRODUCT_KEYS, "products[%r]" % slot)
        # No id, no label. The slot ("TL") is a position, not a compound number,
        # and it used to be printed in bold under the product as if it were one.
        e = {"frag": _rd(slot, d["smiles"], d, "product %r" % slot, slot),
             "id": d.get("id", "")}
        if d.get("yield"):
            e["yield"] = d["yield"]
        if d.get("foot"):
            e["foot"] = d["foot"]
        # PASSED THROUGH UNVALIDATED, on purpose: layout_radial._inbound is the one
        # place the string is read and the one place it is refused, so a bad value
        # gets the engine's message rather than a second, differently worded one
        # here. Absent stays absent, so a call that does not mention direction hands
        # the engine an entry it cannot tell from the ones it got before this.
        if d.get("direction"):
            e["direction"] = d["direction"]
        # A bare string is taken as one line -- block_box normalises it anyway, so
        # rejecting it would be a rule with no reason behind it.
        for _k in ("above", "below"):
            if d.get(_k):
                e[_k] = [d[_k]] if isinstance(d[_k], str) else list(d[_k])
        prod[slot] = e
    # Same gap as draw_linear_scheme without its Allocator, one step less severe.
    # layout_radial builds its own Allocator when handed none -- there the map is
    # load-bearing page-extent bookkeeping, not a shadow -- so the measuring DOES
    # happen; the findings just die with the local variable. Detection alive,
    # reporting dead. Pass one in and the report has somewhere to go.
    from allocator import Allocator as _Alloc
    _alloc = _Alloc()
    _recipe = {"center": center, "products": prod}
    if _rdeco:
        _recipe["decorations"] = _rdeco
    if footnotes:
        _recipe["footnotes"] = _footnote_lines(footnotes)
    if trunk:
        _recipe["trunk"] = True
    _shift: list = []
    # Built only from what the caller actually passed, so a call that names
    # neither hands render_radial cfg=None and cannot take a different path.
    # ---- the two gaps, solved unless the caller pinned them --------------------
    # None and "auto" are THE SAME REQUEST. There is no third behaviour hiding
    # between "I did not say" and "work it out": both mean the figure decides, and
    # keeping them distinct would only be useful for a default that is not solved.
    # They part company in exactly one place, below, where "auto" is an explicit
    # request for something unsound and silence is not.
    _AUTO = frozenset({"auto"})
    _col_in = _gap_value("col_gap", col_gap, keywords=_AUTO)
    _row_in = _gap_value("row_gap", row_gap, keywords=_AUTO)
    _cfg: dict = {"footnote_width": float(footnote_width)} if footnote_width else {}
    _how: dict = {}
    # SET BEFORE THE SOLVERS RUN, because the row solver reads it: a rotated caption
    # is displaced by cap_gap + h/2 instead of cap_gap + _support(w, h, ..), and
    # solving the pitch for the wrong displacement is solving nothing. Threaded
    # through cfg to both solvers below for the same reason.
    if rotate_captions:
        _cfg["cap_rotate"] = True
    # THE COLUMN PITCH is solvable for any figure -- closed form, no search and no
    # allocator probing. See layout_radial.auto_col_gap for what sets it.
    if isinstance(_col_in, float):
        _cfg["col_gap"], _how["col"] = _col_in, "pinned"
    else:
        _cfg["col_gap"], _how["col"] = R.auto_col_gap(center, prod, _cfg, trunk=trunk), "solved"
    # THE ROW PITCH is solved AFTER the column, and it needs it: a diagonal spoke's
    # direction is set by both pitches at once, so where its caption lands depends on
    # the column that was just chosen. Hence col_gap is passed down rather than the
    # row term guessing at it.
    #
    # THE GATE IS NOT "ANY CAPTIONED DIAGONAL". Refusing whenever a DIAGONAL spoke
    # carries a caption would rest on the fact that the floor there goes to 149.4
    # where a formula without the caption term says 115.5. That evidence is right
    # and the conclusion is not -- such a formula is INCOMPLETE. The missing term is
    # the caption block itself, displaced along the shaft normal, whose vertical
    # component carries it into the very gap being solved; with it the floor is
    # predicted to 0.00000 on all 44 profile x caption-line-count points where that
    # family binds. Since typical radial figures carry multi-line conditions on the
    # diagonals, this is the difference between a solved pitch of about 50-170 and
    # the flat 285.
    #
    # WHAT STILL DECLINES is a caption on a VERTICAL spoke with a diagonal beside it:
    # that caption sits in the column corridor and the diagonal's shaft sweeps through
    # it, colliding continuously from row 112 to 496, so no pitch clears it. Those
    # figures keep the historical 285 and SAY SO: reporting a declined solve as
    # "pinned" would hide which of the two happened.
    if isinstance(_row_in, float):
        _cfg["row_gap"], _how["row"] = _row_in, "pinned"
    elif prod and R.row_gap_is_sound(prod, trunk):
        _cfg["row_gap"], _how["row"] = R.solve_row_gap(
            center, prod, _cfg, trunk=trunk, col_gap=_cfg["col_gap"]), "solved"
    elif _row_in == "auto":
        raise ValueError(
            "row_gap=\"auto\" cannot be solved for this figure: a VERTICAL spoke "
            "(the T or B slot) carries a caption while a DIAGONAL spoke stands "
            "beside it. That caption sits in the column corridor, and the diagonal "
            "spoke's SHAFT sweeps through the corridor as the row pitch changes, so "
            "no row spacing clears it and there is nothing here to solve for. A "
            "caption on a DIAGONAL spoke is solved, and so is a captioned T and B "
            "with no diagonal present. Pass a number instead, set trunk=True so the "
            "diagonals become horizontal, or move the caption off the vertical slot.")
    else:
        _cfg["row_gap"], _how["row"] = R.CFG["row_gap"], "default"
    _gaps = {"col": round(_cfg["col_gap"], 2), "row": round(_cfg["row_gap"], 2),
             "solver": _how}
    _fold: dict = {}
    doc = R.render_radial(_recipe, _hdr, cfg=_cfg, page_mode=page,
                          alloc=_alloc, shift=_shift, footnote_report=_fold)
    # `structures` rather than None now that radial has one: a slot repeating a
    # ligand the figure already drew is then ONE structure under one id, and a
    # slot that names an id the figure used for something else says so.
    doc = _add_graphics(doc, graphics, slots, structures, abbreviate,
                        alloc=_alloc, shift=_shift)
    _path = _write_and_open(doc, filename, clean_up)
    out = ""
    _lrep = _layout_report(_alloc)
    if _lrep:
        out += "\n\nLAYOUT REPORT -- " + _lrep
    return _figure_report(_path, structures, alloc=_alloc,
                          notes=out + _chem_report(), gaps=_gaps,
                          footnotes=_fold)


@mcp.tool()
@_cx_doc
def draw_snaking_scheme(rows: list, filename: str = "snaking", Lx: float = 40.0,
                        Rx: float = 600.0, clean_up: bool = False, page: str = "letter",
                        abbreviate: bool = False, graphics: list = None,
                        slots: list = None, conditions: dict = None,
                        reserve: list = None, span: list = None) -> str:
    """Draw an EQUAL-WIDTH SNAKING multi-row scheme as a CDXML file.

    rows: list of rows (top to bottom); each row = list of cells
      {"smiles": .., "id": .., "align"?: "L"|"R",
       "rlabels"?: {..}, "floating"?: [..], "highlight"?: [...]}
      "highlight" is the same spec draw_molecule takes, all of it -- "color",
      "fill", "squiggle", the five bond styles, "halo" / "dot" / "locants".
    align says which end a row starts at: a row whose FIRST cell is 'R' (or whose
    LAST cell is 'L') reads right-to-left. For a snake, give row 2's first cell 'R'
    and its last 'L', and alternate. With no backwards row there is no snake: rows
    wrap at the margin and the step between rows is an arrow at the start of the
    next row. Every row but the last is justified to the widest one. Lx / Rx and
    p are not read; the row edges come from `page`.
    abbreviate=True collapses recognised pendant groups to contracted superatom labels.

    conditions: a caption over the arrow BETWEEN two cells, keyed "from->to" by the
      ids above -- {"1->2": ["NaBH~4~", "MeOH, 0 °C"]}. One string per line.
      Ids are the ones you gave (or the running
      numbers the engine invents when you give none), and a key naming a pair with
      no arrow between them draws nothing -- the arrows follow cell order within a
      row, plus one turn from each row's last cell to the next row's first.
    Text markers: ~sub~ ^sup^ *italic* #bold# _underline_.

    reserve: leave room BELOW a row for something a CALLER will draw there -- the
      catalytic cycle that belongs under one step's arrow, drawn by its own engine
      call and parked in the gap.
        [{"after": "8-2", "h": 208}]
      "after" names a cell; the room goes under THAT CELL'S ROW. Height only.
    span: fix HOW FAR APART two cells sit, in points -- `reserve` answers "leave
      this much room below", this answers "put these two exactly this far apart".
        [{"from": "III", "to": "II", "w": 308.7}]
      CENTRE TO CENTRE; "from"/"to" name a pair, not a direction. The
      distance is paid for by lengthening the arrow(s) between the two cells, and
      it is split evenly when there is more than one.
      WHEN YOU NEED IT: this scheme has to match a distance fixed by another figure
      (e.g. two intermediates of a separately drawn catalytic cycle). Anchoring one
      end is not enough -- the distance inside this figure has to be set.
      An arrow is never shortened -- its caption sets its length.
      ON THIS ENGINE, one extra thing to know: the row carrying a span is no
      longer flush with its neighbours, and its turn becomes a diagonal elbow.
    """
    # REFUSE A BAD SPEC BEFORE ANYTHING IS BUILT. First statement in the
    # tool, ahead of _new_batch, so a rejected call leaves no per-batch state
    # behind either. See preflight.py for what it checks and why the engines'
    # own key whitelists are not enough.
    PF.preflight_snaking({"rows": rows, "page": page, "reserve": reserve,
                          "span": span})
    _new_batch()
    # BEFORE the first _register (see the note in draw_scope_table). ONE source:
    # a snaking figure is cells and nothing else.
    _hdr = _extend_palette(C.acs_header_prefix(),
                           [c.get("highlight") for r in (rows or []) for c in r
                            if isinstance(c, dict) and c.get("highlight")])
    structures: dict = {}
    spec_rows: list = []
    for _ri, r in enumerate(rows):
        cells = []
        for _ci, c in enumerate(r):
            CY.check_keys(c, _SNAKING_CELL_KEYS, "rows[%d][%d]" % (_ri, _ci))
            sid = _register(structures, c["smiles"], c,
                            auto=str(len(structures) + 1),
                            abbreviate=abbreviate)
            cell = {"id": sid}
            if c.get("align"):
                cell["align"] = c["align"]
            elif "p" in c:
                cell["p"] = c["p"]
            cells.append(cell)
        spec_rows.append(cells)
    _recipe = {"structures": structures, "rows": spec_rows}
    if reserve:
        _recipe["reserve"] = reserve
    if span:
        _recipe["span"] = span
    if _DECOR:
        _recipe["decorations"] = dict(_DECOR)
    if conditions:
        # Keys are checked here rather than in the engine, which reads the dict with
        # .get() and so cannot tell a typo'd pair from a pair with no arrow. An id
        # that is not in the figure at all is a typo every time.
        for _k in conditions:
            _bad = [p for p in str(_k).split("->") if p not in structures]
            if len(str(_k).split("->")) != 2 or _bad:
                raise ValueError(
                    "conditions key %r: write it as \"from->to\" using the cell ids "
                    "in this figure (%s)%s"
                    % (_k, ", ".join(sorted(structures)),
                       "; %s is not one of them" % ", ".join(map(repr, _bad))
                       if _bad else ""))
        _recipe["conditions"] = {k: (v if isinstance(v, list) else [v])
                                 for k, v in conditions.items()}

    # ---- ONE PIPELINE. Translate the snake into the linear engine's own rows ----
    # layout_snaking placed by two lines of arithmetic -- x from align/p, y from a
    # row gap fixed at 130 pt -- and measured nothing: no fold_box, no row_gap_y,
    # no struct_low, no ink_extra, not one collision check. A structure taller than
    # 130 pt ran into the row below, and the report for it was not clean, it was
    # absent.
    #
    # The snaking convention is that a row is written in REACTION order and its
    # cells are pinned by align/p; the linear convention is that layout order IS
    # reading order. So a row whose reaction starts at the right is handed over
    # reversed, with its arrows drawn head-first the other way -- the same
    # transform wrap_style='snake' performs on a row it broke itself.
    _cond = _recipe.get("conditions", {})
    _lin_rows, _turns, _back_rows = [], [], set()
    for _ri, _cells in enumerate(spec_rows):
        _ids = [c["id"] for c in _cells]
        # A row reads backwards when its FIRST cell is pinned to the right edge or
        # its LAST to the left -- which is how the docstring tells callers to write
        # one ("put the row's first/last structures at align 'L'/'R'").
        _back = (_cells[0].get("align") == "R") or (_cells[-1].get("align") == "L")
        _seq = list(reversed(_ids)) if _back else _ids
        if _back:
            _back_rows.add(_ri)
        _row = []
        for _k, _sid in enumerate(_seq):
            _row.append({"s": _sid})
            if _k == len(_seq) - 1:
                continue
            # in REACTION terms the step still runs the way the caller wrote it
            _a, _b = (_seq[_k + 1], _seq[_k]) if _back else (_seq[_k], _seq[_k + 1])
            _arw: dict = {}
            if _cond.get("%s->%s" % (_a, _b)):
                _arw["above"] = _cond["%s->%s" % (_a, _b)]
            if _back:
                _arw["reverse"] = True
            _row.append({"arrow": _arw})
        _lin_rows.append(_row)
        if _ri + 1 < len(spec_rows):
            _a, _b = _ids[-1], spec_rows[_ri + 1][0]["id"]
            _t = {"from": _a, "to": _b, "kind": "turn"}
            if _cond.get("%s->%s" % (_a, _b)):
                _t["above"] = _cond["%s->%s" % (_a, _b)]
            _turns.append(_t)
    # NO ROW READS BACKWARDS, SO THERE IS NO SNAKE TO DRAW -- and the honest answer
    # is the one the straight wrap already gives: break the row, let the reader jump
    # the margin, and leave the step at the break as an ARROW AT THE START OF THE
    # NEXT ROW. `align` is the only way a caller says which end a row begins at, so
    # a figure with none of it has not asked for a boustrophedon.
    #
    # THE OTHER TWO ANSWERS ARE BOTH WORSE, measured on one six-compound route:
    #
    #   * Keep the turn and let auto_under place the next row under the compound it
    #     points at. Each row then starts where the last one
    #     ENDED: a descending staircase, 1056 pt wide, 72% over letter.
    #   * Keep the turn and switch auto_under off. The rows line up at the margin
    #     and the turn becomes a 390 pt carriage return sweeping back along row 2's
    #     own baseline, through a label in it. That is the defect the old
    #     layout_snaking had, and the reason it is not worth having back.
    #
    # Inventing the alternation instead -- reversing every other row so it IS a
    # snake -- is the one thing not on the table: a route reads left to right and a
    # new row starts again at the left unless the caller says otherwise.
    #
    # The "NOT A SNAKE" note below still fires, and names the one-word change.
    if not _back_rows:
        for _ri, _t in enumerate(_turns):
            _arw = {}
            if _t.get("above"):
                _arw["above"] = _t["above"]
            _lin_rows[_ri + 1].insert(0, {"arrow": _arw})
        _turns = []
    _recipe["rows"] = _lin_rows
    _recipe.pop("conditions", None)
    if _turns:
        _recipe["connectors"] = _turns

    from allocator import Allocator as _Alloc
    _alloc = _Alloc()
    _shift: list = []
    L.DECOR_UNSCALED.clear()
    L.ACCOUNT_DIFF.clear()
    L.CAPTION_MOVED.clear()
    L.CONN_DETOURED.clear()
    L.SPAN_SHORT.clear()
    doc = L.render_linear(_recipe, _hdr, page_mode=page, snake_rows=_back_rows,
                          alloc=_alloc, shift=_shift)
    doc = _add_graphics(doc, graphics, slots, structures, abbreviate,
                        alloc=_alloc, shift=_shift)
    _path = _write_and_open(doc, filename, clean_up)
    out = ""
    # NOT A SNAKE, SAID IN WORDS. Rows that all read left-to-right make a
    # DESCENDING STAIRCASE: each one starts under where the last one ended. That is
    # a perfectly good figure and it is drawn cleanly -- but it is not what this
    # tool's name promises, and the caller asked for the name. The old engine said
    # so only by failing: it swept the turn back across the whole of row 2, through
    # every structure in it, and left the overlaps in the report for someone to
    # read backwards into a cause.
    # THREE PARAMETERS THAT NO LONGER DECIDE ANYTHING, said out loud. Lx/Rx never
    # did -- the tool wrote them into the recipe and layout_snaking read its own
    # cfg, so four wildly different pairs produced one identical file and nothing
    # mentioned it. `p` did, and now does not: it existed because the old engine had
    # no layout at all and had to be told every x. The linear pipeline measures the
    # row and justifies it, and honouring an arbitrary proportion on top of that can
    # only put structures back on top of one another.
    _pinned = [c for r in (rows or []) for c in r
               if isinstance(c, dict) and "p" in c]
    if _pinned:
        out += ("\n\n%d cell(s) carry `p`, and it is no longer read. Positions come from "
                "the measured width of each structure and the arrows between them, with "
                "every row but the last justified to the widest one. Use align 'L'/'R' "
                "to say which end a row starts at; to line a structure up under one in "
                "an earlier row, use draw_linear_scheme, whose cells take `under`." % len(_pinned))
    if Lx != 40.0 or Rx != 600.0:
        out += ("\n\nLx/Rx are ignored: the row edges come from `page` ('letter' | "
                "'double' | 'single').")
    if len(spec_rows) > 1 and not _back_rows:
        out += ("\n\nNOT A SNAKE: every row reads left-to-right, so the rows wrap at the "
                "margin and each new row starts with an arrow at its left end. For a "
                "snake, pin the NEXT row's first structure to the edge the one above "
                "ended at: align 'R' on row 2's first cell and 'L' on its last.")
    out += _span_short_note()
    _lrep = _layout_report(_alloc)
    if _lrep:
        out += "\n\nLAYOUT REPORT -- " + _lrep
    return _figure_report(_path, structures, alloc=_alloc,
                          notes=out + _chem_report())


# ============================ stereochemistry tools ============================
# Read + audit the stereochemistry you drew (never invent it). RDKit reads SDF/MOL,
# so ingest takes an SDF exported from ChemDraw (File > Save As > SDF) OR isomeric
# SMILES (Copy-as-SMILES). See stereo_tools.py.

@mcp.tool()
def merge_figures(files: list, filename: str = "merged", arrange: str = "column",
                  gap: float = 36.0) -> str:
    """Put figures you have ALREADY DRAWN into ONE ChemDraw file, in the order
    given -- e.g. the schemes of a paper, drawn one call at a time. Nothing is
    redrawn: each figure keeps its layout exactly and is only moved.

    files: the figures, by the `filename` you drew them with ("scheme1" or
      "scheme1.cdxml"); a download link or a local path to one works too. Only
      figures drawn under your own name are found, and on the shared server a
      figure is kept 8 days.
    arrange: "column" (top to bottom) or "row" (left to right), centred on
      each other.
    gap: space between neighbouring figures, in points (36 = half an inch);
      0 or more. The merged page may be at most 40000 pt on either side."""
    from urllib.parse import unquote
    if arrange not in ("column", "row"):
        raise ValueError('arrange is "column" or "row", not %r' % (arrange,))
    if not isinstance(files, list) or not 2 <= len(files) <= 40:
        raise ValueError("files: give a list of 2 to 40 figure names")
    # A negative gap silently reversed the figures' order (-500 put Scheme 2
    # above Scheme 1); an infinite one has no page at all.
    gap = 36.0 if gap is None else float(gap)
    if not 0 <= gap < float("inf"):
        raise ValueError("gap must be 0 or more points, not %r" % (gap,))
    base = _OUT_DIR.get() or OUT_DIR
    # The same name mangling _write_and_open applied when each figure was written,
    # so the name the caller drew with finds the file; a link or a path is cut to
    # its last component first, which also leaves no way out of `base`.
    def stem(f):
        f = unquote(str(f).split("?")[0].rstrip("/").rsplit("/", 1)[-1])
        f = f[:-6] if f.lower().endswith(".cdxml") else f
        return f.replace(" ", "_").replace("/", "_")
    stems = [stem(f) for f in files]
    out = filename.replace(" ", "_").replace("/", "_")
    if out in stems:
        raise ValueError("filename %r is one of the figures being merged and would "
                         "overwrite it -- pick another name" % filename)
    missing = [f for f, s in zip(files, stems) if not (base / (s + ".cdxml")).is_file()]
    if missing:
        have = sorted((p for p in base.glob("*.cdxml")),
                      key=lambda p: p.stat().st_mtime, reverse=True)[:30]
        raise ValueError("not found among your figures: %s. Most recent ones: %s"
                         % (", ".join(map(str, missing)),
                            ", ".join(p.stem for p in have) or "none"))

    figdir = str(Path(__file__).resolve().parent / "figures")
    if figdir not in sys.path:
        sys.path.insert(0, figdir)
    import cdxdoc as CDOC
    import compose as CP
    comp = CP.Composer(margin=48.0, workdir=str(base))
    for i, s in enumerate(stems):
        name = "f%d" % i
        comp.parts.append(CP.Part(name, CDOC.Doc((base / (s + ".cdxml")).read_text(), name), {}))
        if i:
            if arrange == "column":
                comp.below(name, "f%d" % (i - 1), gap=float(gap))
            else:
                comp.right_of(name, "f%d" % (i - 1), gap=float(gap))
    doc, _shift, size = comp.build()
    # 40000 pt (14 m, about 930 benzene rings in a row) is the largest side
    # measured to open in ChemDraw, in under half a second; no real figure is near it.
    if max(size) > 40000:
        raise ValueError("the merged page would be %d x %d pt; the limit is 40000 pt "
                         "on either side -- merge fewer figures, or use a smaller gap"
                         % (round(size[0]), round(size[1])))
    # Same exit check as Composer.write: each figure gets its own id range, so a
    # repeat here means one input already carried it.
    body = doc[doc.index("<page"):]
    ids = re.findall(r'\bid="(\d+)"', body[body.index(">") + 1:])
    # One pass: ids.count() per id was quadratic -- 16.8 s on a 3.8 MB merge.
    dup = sorted((i for i, n in collections.Counter(ids).items() if n > 1), key=int)
    path = _write_and_open(doc, filename, False)
    rep = {"ok": not dup, "path": path, "merged": stems,
           "page": [round(size[0]), round(size[1])]}
    if dup:
        rep["warnings"] = ["repeated object ids %s inside one of the inputs; "
                           "ChemDraw may misplace those objects" % dup[:10]]
    if any(n.startswith("字体") for n in comp.notes):
        rep.setdefault("warnings", []).append(
            "two figures define the same font id differently; the first one's was kept")
    return json.dumps(rep, indent=1, ensure_ascii=False)


@mcp.tool()
def stereo_ingest(sdf_path: str = "", smiles: list = None) -> str:
    """Read the stereochemistry ChemDraw perceived from your wedges. Pass a list of
    isomeric SMILES via smiles (ChemDraw: Edit > Copy As > SMILES). sdf_path is not
    available on this server. Returns canonical isomeric SMILES
    plus a stereocenter census (defined R/S + any undefined centers) per structure.
    Never invents stereo — only reads what you drew."""
    recs = ST.ingest(sdf_path=(sdf_path or None), smiles=(smiles or None))
    return ST.format_ingest(recs)


@mcp.tool()
def stereo_lint(entries: list) -> str:
    """Audit a scheme's stereochemistry. entries = [[id, smiles], ...]. Flags parse
    failures and UNDEFINED stereocenters (centers still needing a source). Reliable
    per-structure report; to carry stereo from one structure to another, use
    stereo_propagate."""
    pairs = [(str(e[0]), e[1]) for e in entries]
    report, n = ST.lint(pairs)
    return report + f"\n\n{n} structure(s) flagged (parse error or undefined center)."


@mcp.tool()
@_cx_doc                 # draws through smiles_to_fragment like the six above
def stereo_gallery(entries: list, filename: str = "stereo_gallery",
                   clean_up: bool = False) -> str:
    """Render every structure WITH WEDGES into one ChemDraw file for a visual diff
    against your drawn scheme. entries = [[id, smiles], ...]."""
    # _new_batch FIRST. The reports in smiles_to_frag are
    # per-batch; reading them without resetting would have this tool describe
    # whatever the PREVIOUS call left behind -- the same defect the depiction
    # library had, and the MCP server never restarts to wash it out.
    _new_batch()
    pairs = [(str(e[0]), e[1]) for e in entries]
    doc, skipped = ST.gallery_cdxml(pairs, C.acs_header_prefix())
    # THE COUNT IS WHAT WAS DRAWN, NOT len(pairs), WHICH IS ONLY WHAT WAS ASKED FOR.
    # Every other loss in this pipeline is at least visible on the page; a
    # structure missing from a proofreading gallery looks exactly like one that
    # passed proofreading.
    out = ("Drew stereo gallery (%d of %d structures): %s."
           % (len(pairs) - len(skipped), len(pairs),
              _write_and_open(doc, filename, clean_up)))
    if skipped:
        out += ("\n\nNOTE: %d entr(y/ies) could not be built and are NOT in the "
                "file: %s. Nothing on the page says so -- check the SMILES."
                % (len(skipped),
                   "; ".join("%s (%s)" % (cid, why)
                             for cid, _s, why in skipped[:5])))
        if len(skipped) > 5:
            out += " (%d more.)" % (len(skipped) - 5)
    return out + _chem_report()


@mcp.tool()
def stereo_propagate(anchors: list, targets: list) -> str:
    """Draw the anchors once -> fill @/@@ across the whole scheme. anchors = list of
    ISOMERIC SMILES (defined-stereo reference structures, e.g. from stereo_ingest).
    targets = [[id, flat_smiles], ...]. Transfers perceived chirality onto each
    target's shared scaffold (MCS match), leaves NEW or substituent-changed centers
    UNDEFINED (never guesses), and AUTO-RUNS stereo_lint on the result. Centers still
    flagged undefined need their own source anchor."""
    tpairs = [(str(t[0]), t[1]) for t in targets]
    results = ST.propagate(list(anchors), tpairs)
    return ST.format_propagate(results)


@mcp.tool()
def molecular_properties(smiles: list) -> str:
    """The numbers a paper reports about a molecule: formula, average MW, exact
    (monoisotopic) mass, tPSA, Crippen logP, HBD, HBA, rotatable bonds, heavy
    atoms, rings. Takes a LIST of SMILES; a bad one comes back as an error entry
    beside the others rather than losing them.

    WHEN QUOTING THESE BACK, carry the method with the number:

    - `logp_crippen` is RDKit's Wildman-Crippen ESTIMATE. It is not ChemDraw's
      cLogP and not an experimental logD -- those are different fitted models and
      routinely disagree by 0.3-0.5. Never call it "cLogP".
    - `tpsa` is RDKit's. It can differ from PubChem's value: caffeine
      is 61.82 here and 58.44 there, entirely from how the two perceive the amide
      nitrogens of the ring. Say it is RDKit's if the user is comparing.
    - `mw` goes on a bottle and into a stoichiometry calculation; `exact_mass` is
      what "HRMS calcd for" is compared against. They are not interchangeable.

    Nothing is drawn and no file is written."""
    import descriptors as DESC
    return json.dumps({"molecules": DESC.describe_all(list(smiles or []))},
                      ensure_ascii=False, indent=1)


@mcp.tool()
def lookup_structure(names: list) -> str:
    """The structure of a NAMED compound, from PubChem -- instead of writing its
    SMILES from memory.

    USE THIS BEFORE DRAWING ANYTHING THAT HAS A NAME: a reagent, a catalyst, a
    drug, a natural product. Then draw the `smiles` it returns, copied exactly;
    re-typing a SMILES is how chirality flips. A remembered SMILES reads exactly
    like a correct one and nothing downstream can tell: the InChI check confirms
    the drawing matches the SMILES you gave, not that it is the compound.

    names: a LIST of names or CAS numbers, e.g. ["vancomycin", "DIAD",
    "2172794-56-6"]. One that fails comes back as an error entry beside the rest.

    ONLY PUBLISHED NAMES. Each name is sent to PubChem (NCBI) and nothing else
    is. Never send the user's compound numbers ("3a"), internal codes or names
    of unpublished compounds: they leave the machine and are not found anyway.
    A generic structure with R groups cannot be looked up.

    READ `title` AND `warnings` ON EVERY RECORD, and pass the warnings on:
    - "PubChem took X to mean Y": say it, and let the user confirm. TMS came
      back as a tetramethoxystilbene and lutidine as 2,6-lutidine; for an
      abbreviation that expanded correctly (DMAP) it is a one-line confirmation.
    - undefined stereocenters: PubChem does not know them either. Look up a name
      that states the isomer ("(-)-menthol"). NEVER add @/@@ yourself.
    - a name asking for a stereoisomer whose record has no stereo: (S)-BINAP and
      (R)-BINOL come back racemic, because axial chirality has no SMILES
      notation. Say the drawing cannot show it.
    - "Contains <metal>": PubChem often breaks metal-ligand bonds (n-BuLi comes
      back as a Li atom beside a butyl radical). In a scheme, write such a
      reagent as text on the arrow rather than drawing it.
    - several records under one name (DIAD): do not pick one silently.

    NOT FOUND IS NOT A REASON TO GUESS. Many abbreviations are not in PubChem
    (DMP, TEA, PCC): try the full name or the CAS number, then ask the user. If
    `unreachable`, say PubChem could not be reached -- do not fall back to a
    SMILES from memory.

    Nothing is drawn and no file is written."""
    import pubchem as PC
    return json.dumps({"source": "PubChem PUG-REST",
                       "results": PC.lookup_all(list(names or []))},
                      ensure_ascii=False, indent=1)


@mcp.tool()
def reaction_amounts(scale: dict, reagents: list) -> str:
    """How much of each reagent to weigh out, for one step at a given scale.

    DO NOT DO THIS ARITHMETIC YOURSELF. Call this tool and quote the numbers it
    returns. Your multiplication errors read exactly like your correct answers,
    and "weigh 1.38 g" in someone's supporting information is not a place to be
    approximately right.

    scale: {"mmol": 5.0} OR {"mass_mg": 1000, "smiles": "Cc1ccc(Br)cc1"} -- give
        the limiting reagent's SMILES and the molecular weight is looked up.
    reagents: one object per reagent.
        {"name": "K2CO3", "equiv": 2.0, "form": "solid", "smiles"|"mw": ...}
        {"name": "Et3N",  "equiv": 3.0, "form": "neat",   "mw": ..., "density": 0.726}
        {"name": "n-BuLi","equiv": 1.2, "form": "solution", "molarity": 2.5}
        {"name": "THF",   "form": "solvent", "concentration": 0.2}
        optional "purity": 0.60 for a dispersion or an aqueous reagent.
        `smiles` may replace `mw` anywhere; it is resolved here.

    IF `ask` COMES BACK NON-EMPTY, ASK THE USER THOSE QUESTIONS. Do not supply the
    missing density, molarity or concentration from memory, and do not pick a
    "typical" value. The only correct source for what is in the user's bottle is
    the user: a table saying n-BuLi is 2.5 M, beside a bottle that says 2.2 M,
    produces a wrong number with no sign that anything went wrong. Measured, the
    same n-BuLi is 2.40 mL as a 2.5 M solution and 0.57 mL treated as a neat
    liquid -- a factor of 4.2, and the reaction is dead.

    Ask for the scale too if it was not given, rather than assuming one.

    WRITING IT UP: ACS style, as in a J. Org. Chem. / Org. Lett. SI. Use the
    `display` string for each reagent -- it is already in the units and precision
    a procedure is written in. Molecular weights come back as `scale_mw` and each
    reagent's `mw`; quote those rather than recalling one, and call
    molecular_properties for anything else you need a number for. A remembered
    molecular weight is the one number in a procedure that nobody checked. Addition order, temperature and time are NOT
    computed here; if the user has not said, ask or state plainly that you are
    proposing them.

    Nothing is drawn and no file is written."""
    import descriptors as DESC
    import stoichiometry as STOICH

    scale = dict(scale or {})
    rs = [dict(r) for r in (reagents or [])]
    bad = [e for e in [DESC.fill_mw(d) for d in [scale] + rs] if e]
    out = STOICH.compute(scale, rs)

    if bad:
        out["ok"] = False
        out["ask"] = ["SMILES could not be parsed: " + b for b in bad] + list(out.get("ask", []))
    return json.dumps(out, ensure_ascii=False, indent=1)


@mcp.tool()
def check_disconnection(target: str, precursors: list) -> str:
    """Check ONE retrosynthetic disconnection: does every heavy atom of the target
    come from one of the proposed precursors?

    USE THIS ON YOUR OWN PROPOSALS, EVERY TIME, BEFORE SHOWING THEM. A
    retrosynthesis is a guess, and a wrong one reads exactly like a right one.
    This settles the single part of the claim that is mechanically falsifiable:
    atoms cannot appear out of nothing. If it comes back with atoms unaccounted
    for, the disconnection is wrong regardless of how good the reaction names
    sound, and you should say so and try again rather than present it.

    ONE DISCONNECTION AT A TIME. For a multi-step route the user has asked for,
    check each step separately -- the target of step n is the precursor of step
    n+1.

    READ `shape` BEFORE YOU READ `ok`. With two or more precursors ("coupling"),
    an unaccounted atom is a defect: say the disconnection is wrong. With ONE
    ("fgi"), it may equally be the reagent's -- Vilsmeier formylation is short a
    carbon and an oxygen because DMF brings them, and calling that impossible
    would be wrong. Code cannot separate a reagent-supplied atom from a forgotten
    precursor, so on a single-precursor step report BOTH readings and let the
    conditions decide, rather than announcing a verdict this tool did not reach.

    WHAT A PASS DOES NOT MEAN, and you must pass this on rather than imply
    otherwise. `not_checked` comes back on every reply, including failures, and it
    is not decoration. The check is blind to connectivity: 4-methylbiphenyl from
    2-BROMOtoluene plus phenylboronic acid accounts for all thirteen atoms and
    passes, and it is wrong -- the methyl ends up ortho. It is equally blind to
    stereochemistry, chemoselectivity, protecting groups, whether the reaction has
    any literature precedent, and whether the conditions exist.

    So the reply to the user is "the atoms account for, and here is what that does
    and does not establish" -- never "this route works". When you name the
    disconnection, name the reaction type and how common the precedent is for a
    substrate like this one, and say which step you would worry about. That is
    what lets a chemist reject it in thirty seconds, which is the useful thing.

    Nothing is drawn and no file is written."""
    import retro as RETRO
    return json.dumps(RETRO.check_disconnection(str(target or ""),
                                                list(precursors or [])),
                      ensure_ascii=False, indent=1)


@mcp.tool()
def write_procedure(scale: dict, reagents: list, limiting: str,
                    names: dict = None, order: list = None,
                    catalysts: list = None, temperature: str = None,
                    time: str = None, addition_temperature: str = None,
                    atmosphere: str = None, workup: str = None,
                    purification: str = None, product: dict = None,
                    actions: list = None) -> str:
    """Write the experimental paragraph for one step, ACS supporting-information style.

    YOU SUPPLY THE WORDS. THIS SUPPLIES THE NUMBERS. Do not put a mass, a volume
    or a millimole figure anywhere in what you pass in, and do not work one out to
    check: every quantity in the returned paragraph is computed from `scale` and
    `reagents` by the same calculator reaction_amounts uses. Your arithmetic is
    wrong in the same voice as your arithmetic is right, and "weigh 1.38 g" in
    someone's supporting information is not a place to be approximately right.

    `scale` and `reagents` are exactly reaction_amounts' arguments. If anything is
    missing -- a density, a molarity -- this refuses to write and hands back the
    same questions. ASK THEM. A procedure with a hole where a mass belongs is
    worse than no procedure, because it reads complete.

    `limiting` IS REQUIRED AND NEVER GUESSED. Usually the substrate, but an
    expensive coupling partner or a scarce ligand can be the limiting one. If the
    user has not said, ask.

    WHAT IS YOURS TO DECIDE, and it is the useful half: `order` (addition order),
    `temperature`, `addition_temperature` (when the charge goes in cold and is
    heated after), `time`, `atmosphere`, and then EITHER `actions` (preferred) OR
    `workup` + `purification` as free prose.

    PREFER `actions`. It is a list of steps drawn from a controlled vocabulary,
    and because it is structured it can be CHECKED. Free prose cannot: a
    misspelling, an invented field or two reactions run together all read fine.
    The 25 types are:

      Add(material, dropwise?, temperature?, atmosphere?, duration?)
      MakeSolution(materials[])   SetTemperature(temperature)
      Stir(duration?, temperature?, atmosphere?)
      Reflux(duration?, dean_stark?, atmosphere?)
      Wait(duration, temperature?)   Microwave(...)   Sonicate(...)
      Irradiate(wavelength?, light_source?, duration?, temperature?, atmosphere?)
      Degas(gas?, duration?)         Quench(material, dropwise?, temperature?)
      PH(material?, ph?, dropwise?, temperature?)
      Extract(solvent, repetitions?) Wash(material, repetitions?)
      Partition(material_1, material_2)   PhaseSeparation()
      CollectLayer(layer)           Filter(phase_to_keep?)
      DrySolution(material?)        DrySolid(duration?, temperature?, atmosphere?)
      Concentrate()   Purify()   Recrystallize(solvent)   Triturate(solvent)
      Yield(material)

    Each is {"action": "Add", "material": {"name": "K2CO3"}, ...}. A CHEMICAL IS A
    NAME AND NOTHING ELSE -- do not give it a quantity. Amounts are substituted
    from the calculator, and a name matching a reagent row picks up its computed
    amount automatically.

    EACH REAGENT IS WRITTEN ONCE. The opening "To a solution of X in Y was added
    Z" is generated from the amounts, and any row your sequence puts in -- by Add,
    MakeSolution, Quench or PH -- is taken out of it and written where your action
    adds it, with its computed amount: a catalyst after degassing, a base dropwise
    at 0 °C. Leave what simply goes in at the start to the opening sentence. Name
    a chemical by its key or its written name; the prose prints the written name.

    A workup step (Extract, Wash, Partition, Recrystallize, Triturate) never
    carries an amount, even when its solvent shares a name with the reaction
    solvent: the computed volume is the reaction's, not what you extracted with.

    IF THE SEQUENCE IS REJECTED, `problems` says why and nothing is written. Fix it
    rather than falling back to prose.

    Two things the vocabulary cannot express, reported when they come up:
    `Purify` has no fields at all -- it can say "purified" and not how -- and a
    chemical you name that is not in the reagent table gets no amount. The second
    is normal for a workup solvent and a missing table row for anything weighed.

    Free prose still works, and is right when the USER dictated the workup. The
    volumes in it are then theirs, not calculated.

    `names` maps each reagent key to its written name ({"K2CO3": "potassium
    carbonate"}); `catalysts` lists the ones to report in mol% instead of equiv,
    as names or {"name": ...}. One that is not a row (no molecular weight, say)
    comes back in `catalysts_not_in_the_table` -- ask for its structure or MW.
    `product` is {"name", "mass_mg"?, "yield_pct"?} for the closing clause -- those
    are the user's measured numbers, not computed ones.

    TELL THE USER WHICH HALF IS WHICH. `note` and `model_wrote` come back on every
    reply for that purpose. A paragraph that reads uniformly authoritative hides
    that its addition order and workup volumes were invented by you.

    Nothing is drawn and no file is written; call draw_linear_scheme for the
    scheme."""
    import procedure as PROC
    import stoichiometry as STOICH

    scale = dict(scale or {})
    rs = [dict(r) for r in (reagents or [])]
    bad = [e for e in [DESC.fill_mw(d) for d in [scale] + rs] if e]
    amounts = STOICH.compute(scale, rs)
    if bad:
        amounts["ok"] = False
        amounts["ask"] = ["SMILES could not be parsed: " + b
                          for b in bad] + list(amounts.get("ask", []))

    out = PROC.render(amounts, limiting=str(limiting or ""), order=order,
                      names=names, catalysts=catalysts,
                      temperature=temperature, time=time,
                      addition_temperature=addition_temperature,
                      atmosphere=atmosphere, workup=workup,
                      purification=purification, product=product,
                      actions=actions)
    out["amounts"] = amounts
    return json.dumps(out, ensure_ascii=False, indent=1)


# ======================== the published input schema ==========================
#
# WHAT FastMCP PUBLISHES WITHOUT THIS. Every tool's inputSchema is generated from
# the signature, and the signature says `items: list`, `header: dict`. That comes
# out as {"type": "array", "items": {}} and {"type": "object"} -- true, and empty.
# It names not one nested key, and JSON Schema's default for `additionalProperties`
# is TRUE, so the published contract for all six actively stated that an unknown
# key is fine. It is not fine anywhere in this file: `check_keys` refuses unknown
# keys at the top of each tool, and without this the schema says the opposite of
# the code.
#
# THE PROPERTY NAMES ARE THE FROZENSETS THEMSELVES, never a second hand-typed list.
# `_obj` takes the whitelist a tool already enforces at runtime and REFUSES to
# narrow a property that whitelist does not have -- at import, so the server does
# not start. A key renamed in one place cannot leave a stale entry in the other.
#
# WHAT THIS LAYER DOES NOT DO, said plainly, because the gap is the whole point of
# the next one. FastMCP validates an incoming call against a pydantic model built
# from the SIGNATURE (`Tool.run` -> `fn_metadata.call_fn_with_arg_validation`), not
# against `Tool.parameters`, and pydantic's default for an extra field is to IGNORE
# it. So `additionalProperties: false` here is a statement to the caller and a
# fixture for the CI invariant; a client that does not validate can still send an
# unknown top-level kwarg and have it dropped in silence. Closing that is the
# preflight layer's job, not this one's.
#
# CLOSED RECORD vs CLOSED MAP. `additionalProperties: false` is only expressible
# for an object with a FIXED key set. A map -- `rlabels` keyed by placeholder
# number, `products` keyed by slot, `conditions` keyed by "from->to" -- has no such
# set, so it is closed the only way a map can be: `propertyNames` on the KEY, where
# the key set is known, and one schema for every value.

def _obj(keys, narrow: dict | None = None) -> dict:
    """A CLOSED object: exactly the keys `keys` names, and nothing else.

    `narrow` types individual properties; everything else is left as `{}` (any),
    because a property this file cannot state a type for honestly is one it should
    not be stating a type for. A name in `narrow` that `keys` does not have raises
    AT IMPORT -- that is the guard: the property list comes from the frozenset the
    tool enforces, so a key renamed there cannot quietly leave a stale narrowing
    behind here.
    """
    props: dict = {k: {} for k in sorted(keys)}
    for k, v in (narrow or {}).items():
        if k not in props:
            raise ValueError(
                "schema: narrowed %r, which is not one of the keys this level "
                "accepts (%s)" % (k, ", ".join(sorted(keys))))
        props[k] = v
    return {"type": "object", "properties": props, "additionalProperties": False}


def _arr(item: dict, **kw) -> dict:
    return dict({"type": "array", "items": item}, **kw)


def _map(value: dict, key_enum=None) -> dict:
    """A map: the caller's own keys, one schema for every value."""
    out: dict = {"type": "object", "additionalProperties": value}
    if key_enum is not None:
        out["propertyNames"] = {"enum": sorted(key_enum)}
    return out


# ---- colour, in the three positions this codebase actually has ---------------
#
# The eight names are cdxml_core.COLORS, i.e. the ACS header's own colortable,
# READ from it rather than retyped. There is no fourth spelling: nothing in this
# repository parses "#rrggbb" -- an RGB colour is an [r, g, b] triple of 0..1
# floats, which is what beautify.extend_colortable appends.
_S_COLOR_NAME = {"type": "string", "enum": sorted(C.COLORS)}
_S_COLOR_INDEX = {"type": "integer", "minimum": 0}
_S_RGB = {"type": "array", "minItems": 3, "maxItems": 3,
          "items": {"type": "number", "minimum": 0.0, "maximum": 1.0}}
# 1. A position `cdxml_core.color_index` reads: a NAME from the document palette,
#    or a raw colortable index. It refuses an [r, g, b] there BY NAME, because
#    only the tool layer holds the header to append one to.
_S_COLOR = {"oneOf": [_S_COLOR_NAME, _S_COLOR_INDEX]}
# 2. A position the tool layer palette-extends FIRST (_extend_palette for
#    highlights, _add_graphics for graphics and slots), where an [r, g, b] is
#    accepted and rewritten to the index it landed on.
_S_COLOR_RGB = {"oneOf": [_S_COLOR_NAME, _S_COLOR_INDEX, _S_RGB]}
# 3. A CAPTION colour, where a LIST already means one colour per line and null
#    means "that line keeps the block colour". An [r, g, b] cannot be spelled here
#    at all -- cdxml_core.color_index rejects a float in this position by name.
_S_CAPTION_COLOR = {"oneOf": [
    _S_COLOR_NAME, _S_COLOR_INDEX,
    _arr({"oneOf": [_S_COLOR_NAME, _S_COLOR_INDEX, {"type": "null"}]})]}


# ---- key sets this file did not already have --------------------------------
#
# Every name below is read off the code that CONSUMES it, the standard
# _LINEAR_CELL_KEYS was built to. These are the hand-derived ones, listed together
# so it is obvious which they are: each names a dict its own module reads with
# .get(), so there was no frozenset to point `_obj` at until one was written here.
#
# FIVE CLOSED KEY SETS LIVE ELSEWHERE and are aliases here.
#
# NOT A CLAIM ABOUT THE REST OF THIS SECTION. The sets still defined below are
# not audited, and two spot checks say at least some of them have the same hole:
# a typo inside `slots[].structures[]` and an unknown key on `header` itself are
# BOTH accepted in silence. They are left as they are rather than fixed here,
# because "which sets are enforced" needs the same check-every-cell pass the
# four above got, not a guess.
#
# ALL FIVE LIVE IN preflight.py: `_HIGHLIGHT_KEYS` and the other four, each the
# same way: the definition goes where it is ENFORCED at runtime, and this file
# points at it so `_obj` below and the check the tools run cannot come apart.
# One list, not two.
#
# WHAT THAT MEASURED: twelve places publish one of these four objects, and in
# all twelve an unenforced typo'd key produced a file BYTE-IDENTICAL to the
# correct spec. Four of the twelve were worse than that -- the scope table's
# header and panel arrows read neither `style` nor `arcs` at all, so a
# correctly spelled key did nothing either. Those two are refused whole rather
# than spell-checked (`preflight.refuse_unread_arrow`).
_HIGHLIGHT_KEYS = PF.HIGHLIGHT_KEYS
_ARROW_STYLE_KEYS = PF.ARROW_STYLE_KEYS
_ARC_ENTRY_KEYS = PF.ARC_ENTRY_KEYS
_LEGEND_KEYS = PF.LEGEND_KEYS
_LEGEND_CELL_KEYS = PF.LEGEND_CELL_KEYS
# One structure inside a `slots` box: _slot_items reads "smiles" and "label",
# _register reads "id", _build_fragment reads the three fragment keys.
_SLOT_STRUCT_KEYS = frozenset({"smiles", "id", "label"}) | _FRAG_KEYS
# An annotation endpoint, checked by hand in draw_linear_scheme.
_ANNOT_END_KEYS = frozenset({"cell", "on"})
# A scope table's `header` is read for ONE key and nothing else.
_SCOPE_HEADER_KEYS = frozenset({"items"})
# table_layout._reject_unknown's five call sites, which pass inline literals.
_TABLE_KEYS = frozenset({"cols", "rows", "groups", "split", "style", "width_mode"})
_TABLE_COL_KEYS = frozenset({"header", "align", "note", "id"})
_TABLE_ROW_KEYS = frozenset({"cells", "note", "mark"})
_TABLE_CELL_KEYS = frozenset({"v", "sec", "note"})
_TABLE_GROUP_KEYS = frozenset({"header", "note", "cols"})
_TABLE_NOTE_KEYS = frozenset({"mark", "text"})


# ---- the pieces every engine shares -----------------------------------------
#
# `rlabels` stays an open map on purpose: its keys are the caller's placeholder
# numbers, and a key that matches no placeholder is ALREADY REPORTED in the return
# value (smiles_to_frag.RLABEL_UNUSED). There is nothing silent there to close.
_S_RLABELS = {"type": "object"}
_S_FLOATING = {"oneOf": [{"type": "boolean"},
                         {"type": "array", "items": {"type": "integer"}}]}
_S_HIGHLIGHT = _arr(_obj(_HIGHLIGHT_KEYS, {"color": _S_COLOR_RGB}))
# The three keys _build_fragment spends, wherever a structure spec is taken.
_S_FRAG = {"rlabels": _S_RLABELS, "floating": _S_FLOATING,
           "highlight": _S_HIGHLIGHT}
# The shaft styles an arrow can actually have, published as an enum
# rather than left open. THREE SURFACES, ONE OBJECT: `_S_ARROW_STYLE` is what
# `items[].arrow.style` and both tools' `connectors[].style` point at, so
# narrowing it once closes all three; `arcs[].line` is the fourth and gets the
# same enum below. The set is preflight's, not a second hand-typed list -- the
# tools enforce it at runtime and this publishes the same names.
_S_ARROW_LINE = {"type": "string", "enum": sorted(PF.ARROW_LINE_TYPES)}
_S_ARROW_STYLE = _obj(_ARROW_STYLE_KEYS, {"color": _S_COLOR,
                                          "line": _S_ARROW_LINE})
_S_LEGEND = _obj(_LEGEND_KEYS, {
    "rows": _arr({"oneOf": [
        {"type": "string"},                            # a bare NOTE row
        _arr({"oneOf": [{"type": "string"},
                        _obj(_LEGEND_CELL_KEYS, {"swatch": _S_COLOR})]})]})})
# `graphics` and `slots` are the same two arguments in all six tools, so they get
# one schema each. The graphics keys are graphics._KEYS itself, per kind.
#
# `kind` IS PINNED PER BRANCH, and it has to be: the three kinds share "kind",
# "at" and "color", so without the const a plain {"kind": "line", "at": [...]}
# satisfies the rect branch as well and `oneOf` -- which means EXACTLY one --
# refuses a graphic that is perfectly legal.
_S_GRAPHICS = _arr({"oneOf": [
    _obj(GX._KEYS[_k], {"kind": {"const": _k}, "color": _S_COLOR_RGB})
    for _k in GX.KINDS]})
_S_SLOTS = _arr(_obj(_SLOT_KEYS, {
    "fill": _S_COLOR_RGB, "label_color": _S_COLOR_RGB,
    "side": {"type": "string", "enum": sorted(_SLOT_SIDES)},
    "align": {"type": "string", "enum": sorted(_SLOT_ALIGNS)},
    "structures": _arr(_obj(_SLOT_STRUCT_KEYS, _S_FRAG))}))


# ---- 1. draw_molecule --------------------------------------------------------
_S_MOLECULE = dict(_S_FRAG, legend=_S_LEGEND)

# ---- 2. draw_linear_scheme ---------------------------------------------------
_S_LINEAR_ARROW = _obj(_LINEAR_ARROW_KEYS, {
    "style": _S_ARROW_STYLE,
    "above_color": _S_CAPTION_COLOR,
    "below_color": _S_CAPTION_COLOR,
    "yield_color": _S_CAPTION_COLOR,
    "arcs": _arr(_obj(_ARC_ENTRY_KEYS, {"line": _S_ARROW_LINE})),
    "reagent_rlabels": _S_RLABELS,
    "reagent_highlight": _S_HIGHLIGHT,
})
# One item of a row: a structure/text cell, a bare plus, or the arrow. Three
# shapes, and _check_linear_item picks between them the same way -- by which of
# "arrow" / "plus" is present.
_S_LINEAR_ITEM = {"oneOf": [
    _obj(_LINEAR_CELL_KEYS, _S_FRAG),
    _obj({"plus"}),
    _obj({"arrow"}, {"arrow": _S_LINEAR_ARROW}),
]}
_S_LINEAR = {
    # ONE row (a list of cells) or SEVERAL (a list of rows). Both, because
    # draw_linear_scheme's first line accepts both.
    "items": _arr({"oneOf": [_S_LINEAR_ITEM, _arr(_S_LINEAR_ITEM)]}),
    "connectors": _arr(_obj(_LINEAR_CONN_KEYS, {"style": _S_ARROW_STYLE})),
    "annotations": _arr(_obj(_ANNOT_KEYS, {
        "color": _S_COLOR,
        "from": _obj(_ANNOT_END_KEYS),
        "to": _obj(_ANNOT_END_KEYS)})),
    "graphics": _S_GRAPHICS,
    "slots": _S_SLOTS,
}

# ---- 3. draw_scope_table -----------------------------------------------------
_S_PRODUCT = _obj(_PRODUCT_KEYS, dict(_S_FRAG, legend=_S_LEGEND))
# A section marker's colour is the one place the scope table names a colour by
# itself. It is NOT wired -- the tool says so in the return value instead of
# dropping it -- so the honest schema is the eight names the palette has and no
# [r, g, b], which nothing downstream would know what to do with.
_S_SECTION = _obj(_SECTION_KEYS, {"color": _S_COLOR_NAME})
_S_SCOPE_ROW_ITEM = {"oneOf": [_S_PRODUCT, _S_SECTION]}
_S_SCOPE = {
    "rows": _arr({"oneOf": [_S_SCOPE_ROW_ITEM, _arr(_S_PRODUCT)]}),
    # ONE branch, not the linear item's three: _HDR_ITEM_KEYS already holds
    # "plus" and "arrow" alongside the structure keys, and header.items is checked
    # against that one whitelist. Splitting it into three would make {"plus": true}
    # match two of them, which `oneOf` refuses.
    "header": _obj(_SCOPE_HEADER_KEYS, {
        "items": _arr(_obj(_HDR_ITEM_KEYS,
                           dict(_S_FRAG, arrow=_S_LINEAR_ARROW)))}),
    "table": _obj(_TABLE_KEYS, {
        "cols": _arr(_obj(_TABLE_COL_KEYS)),
        "rows": _arr(_obj(_TABLE_ROW_KEYS, {
            "cells": _arr({"oneOf": [{"type": "string"}, {"type": "null"},
                                     _obj(_TABLE_CELL_KEYS)]})})),
        "groups": _arr(_obj(_TABLE_GROUP_KEYS))}),
    "table_notes": _arr(_obj(_TABLE_NOTE_KEYS)),
    # Every CFG key is honoured from here (table_layout._cfg does `out.update`),
    # so the schema is CFG's own keys. The four that carry a colour carry an
    # [r, g, b]: table_layout writes them into a fill, not through color_index.
    "table_style": _obj(TL.CFG, {
        "zebra": _S_RGB, "header_fill": _S_RGB, "best_fill": _S_RGB,
        "group_fills": _arr(_S_RGB)}),
    "connectors": _arr(_obj(_GRID_CONN_KEYS, {"style": _S_ARROW_STYLE})),
    "panels": _arr(_obj(_PANEL_KEYS, {
        "items": _arr({"oneOf": [
            _obj(_LINEAR_CELL_KEYS, _S_FRAG),
            _obj({"plus"}),
            _obj({"arrow"}, {"arrow": _S_LINEAR_ARROW})]})})),
    "graphics": _S_GRAPHICS,
    "slots": _S_SLOTS,
}

# ---- 4. draw_catalytic_cycle -------------------------------------------------
# minItems 2, not 3, and the number is the ENGINE's: layout_cycle says "a cycle
# needs at least two species on the ring" in three places, and a checked-in
# reference figure is a two-node ring.
_S_CYCLE_MIN_NODES = 2
# `in` / `out` are LISTS of species (a step can consume two reagents), `aside` is
# type and only type -- a reagent note on the step, one string per line, never a
# structure. Three different shapes under three sibling keys, and the code says so
# at three different places: the `for k, it in enumerate(st[tag])` loop for the
# first two, and layout_cycle's `[aside] if isinstance(aside, str) else list(aside)`
# for the third.
_S_CYCLE_FEEDER = _arr(_obj(_CYCLE_FEEDER_KEYS, _S_FRAG))
_S_CYCLE_ASIDE = {"oneOf": [{"type": "string"},
                            _arr({"type": "string"})]}
_S_CYCLE = {
    "nodes": _arr(_obj(_CYCLE_NODE_KEYS, dict(
        _S_FRAG,
        frame=_obj(CY.FRAME_KEYS),
        content=_arr({"oneOf": [{"type": "string"},
                                _obj(_CYCLE_CONTENT_KEYS, _S_FRAG)]}))),
                  minItems=_S_CYCLE_MIN_NODES),
    "steps": _arr({"oneOf": [
        {"type": "null"},
        _obj(CY.STEP_KEYS, {"style": _obj(CY._ARC_KEYS),
                            "in": _S_CYCLE_FEEDER, "out": _S_CYCLE_FEEDER,
                            "aside": _S_CYCLE_ASIDE})]}),
    "edges": _arr(_obj(_CYCLE_EDGE_KEYS, {
        "style": _obj(CY._ARC_KEYS),
        "feeders": _obj(_CYCLE_FEEDERS_KEYS, {
            "in": _S_CYCLE_FEEDER, "out": _S_CYCLE_FEEDER,
            "aside": _S_CYCLE_ASIDE})})),
    "rings": _arr(_obj(_CYCLE_RING_KEYS, {
        "nodes": _arr({}, minItems=_S_CYCLE_MIN_NODES)})),
    "paths": _arr(_obj(_CYCLE_PATH_KEYS)),
    "arc_style": _obj(CY._ARC_KEYS),
    "graphics": _S_GRAPHICS,
    "slots": _S_SLOTS,
}

# ---- 5. draw_radial_figure ---------------------------------------------------
# The slot names are layout_radial.POS, read from it. An unknown slot already
# raises -- as a bare KeyError, which says nothing -- so naming them here is what
# lets a validating caller be told which eight there are.
_S_RADIAL_PRODUCTS = _map(_obj(_RADIAL_PRODUCT_KEYS, _S_FRAG), key_enum=R.POS)
_S_RADIAL = {
    "products": _S_RADIAL_PRODUCTS,
    "spokes": _S_RADIAL_PRODUCTS,
    "highlight": _S_HIGHLIGHT,
    "graphics": _S_GRAPHICS,
    "slots": _S_SLOTS,
}

# ---- 6. draw_snaking_scheme --------------------------------------------------
_S_SNAKING = {
    "rows": _arr(_arr(_obj(_SNAKING_CELL_KEYS, _S_FRAG))),
    "graphics": _S_GRAPHICS,
    "slots": _S_SLOTS,
}


_TOOL_SCHEMAS = {
    "draw_molecule": _S_MOLECULE,
    "draw_linear_scheme": _S_LINEAR,
    "draw_scope_table": _S_SCOPE,
    "draw_catalytic_cycle": _S_CYCLE,
    "draw_radial_figure": _S_RADIAL,
    "draw_snaking_scheme": _S_SNAKING,
}


def _install_schema(name: str, narrow: dict) -> None:
    """Close one tool's published schema and narrow the parameters named.

    The GENERATED schema is the base, not a thing to replace: it already carries
    every top-level parameter's name, type, default and the `required` list, all
    read off the signature, and re-typing those by hand is the drift this whole
    section is arranged to avoid. So `narrow` only replaces the nested shape of
    the parameters it names -- keeping the generated `title` and `default` -- and
    the one thing added to the top level is the closure.

    A name in `narrow` the signature does not have raises at import, same guard as
    `_obj`'s: a renamed parameter stops the server rather than losing its schema.
    """
    tool = mcp._tool_manager.get_tool(name)
    if tool is None:
        raise ValueError("schema: no tool named %r is registered" % name)
    schema = dict(tool.parameters)
    props = dict(schema.get("properties") or {})
    for k, v in narrow.items():
        if k not in props:
            raise ValueError(
                "schema: %s has no parameter %r to narrow (it takes %s)"
                % (name, k, ", ".join(sorted(props))))
        keep = {j: props[k][j] for j in ("title", "default") if j in props[k]}
        props[k] = dict(v, **keep)
    schema["properties"] = props
    schema["additionalProperties"] = False
    tool.parameters = schema


def tool_schema(name: str) -> dict:
    """The published inputSchema for one tool, for a caller that wants to
    validate against it -- which is what the test suite does."""
    tool = mcp._tool_manager.get_tool(name)
    if tool is None:
        raise ValueError("no tool named %r is registered" % name)
    return tool.parameters


for _name, _narrow in _TOOL_SCHEMAS.items():
    _install_schema(_name, _narrow)


def _accept_explicit_null(tool) -> None:
    """Let an explicit null mean "not given" for every parameter whose default
    is None.

    The published schema says `"default": null` for these, and a client may
    send exactly that -- on 2026-10-29 a draw_molecule call through Claude Code
    arrived with `floating: null` and was refused, because the validator is
    built from the signature and `floating: list = None` admits a list only.
    58 parameters across the tools had the same mismatch. Only the validator
    changes; the published schema is left as it is."""
    am = tool.fn_metadata.arg_model
    loose = {k: (f.annotation | None, f) for k, f in am.model_fields.items()
             if f.default is None and type(None) not in typing.get_args(f.annotation)}
    if loose:
        tool.fn_metadata.arg_model = pydantic.create_model(
            am.__name__, __base__=am, **loose)


for _tool in mcp._tool_manager._tools.values():
    _accept_explicit_null(_tool)


if __name__ == "__main__":
    mcp.run()
