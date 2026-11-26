"""text_metrics.py — ChemDraw (Arial) text-box model, measured not fitted.

Width and height are computed from the font itself rather than from a fitted
per-character table. ChemDraw uses the font's own hmtx advance widths, with no
rounding (FractionalWidths="yes") and exactly linear scaling with point size.
So text size is arithmetic, not estimation.

    width      = sum(advance), NO kerning        <- the default
    ink width  = the same, minus lsb(first glyph) and rsb(last glyph)
    ink height = glyph bounds; ChemDraw adds no vertical padding at all

KERNING IS OFF, and that is not an oversight: ChemDraw uses the font's own
hmtx advance widths but NOT its kern pairs. Over the 5498 comparable text runs
in ChemDraw's own 652
sample files, 91 contain a non-zero Arial kern pair worth -0.0737 em on
average, yet their mean residual differs from the un-kerned runs by +0.00018
em. Least squares puts the fraction of Arial kerning ChemDraw actually applies
at -0.045, i.e. zero. Turning kerning ON moves the bias from +0.0012 to
+0.0749 em and the RMS from 0.019 to 0.086.

This is the one place where matching HarfBuzz and matching ChemDraw are not
the same thing: HarfBuzz applies kern/GPOS by default, ChemDraw does not.

ANCHOR CONVENTION. The <t> p= attribute is the FIRST LINE'S BASELINE, not the
box bottom. Measured against ChemDraw's own saved boxes over 69 test strings:

    bbox.left   - p.x   median  +0.54   (the left side bearing)
    bbox.bottom - p.y   median  +0.12   (descender depth, 0 when none)
    bbox.top    - p.y   median  -7.16   (ascender height)

Baseline and box bottom coincide only for strings with no descender, which is
why reading p as the box bottom-left looks right until a descender appears.
cdxml_core.styled_text writes `baseline = top + asc`, i.e. this convention.

FACES. Bold and italic have their own advance tables; measuring them with the
regular table costs up to 6.7 pt on a bold string and 1.5 pt on an italic one.
Pass bold=/italic= and the matching font file is used.

MEASURED CONSTANTS
    SUBSCRIPT_W   0.7489  subscript/superscript size factor (a ChemDraw setting).
                          Measured against ChemDraw's own files rather than
                          against round-trip test files. Least squares on the
                          3377 face=96 runs gives k = 0.74888 (RMS 0.00243 em);
                          the explicitly marked runs agree -- face=32 subscripts
                          give s_sub = 0.7525, face=64 superscripts s_sup =
                          0.7524. All three within 0.4% of each other, so one
                          constant covers auto-formula and manual markup alike.
                          A value of 0.7729, from round-trip test files, is
                          high by 0.024, which costs 0.013 em per shrunk digit.
    SUB_DROP_EM   0.2120  subscripts sit one descender below the baseline;
                          this equals the font's hhea descent, so it should
                          follow the font rather than being Arial-specific
    PITCH_EM      1.0604  see the caveat below

KNOWN OPEN, both currently erring on the safe side:

  * SUPERSCRIPT. A lone superscript rises 0.5207 em (close to the font's
    x-height). But a block holding BOTH a subscript and a superscript measured
    12.86 where those constants predict about 9.7, so the combination is not
    yet understood and no rise is applied. There is also no API for it: a
    superscript can only be passed in subscript_chars, which would push it
    DOWN. Callers should leave superscripts out of that bag until this is
    measured properly.

  * LINE PITCH. CaptionLineHeight is literally the string "variable", and it
    is: pitch measured 0.63 em between two x-height-only lines and 1.03 em
    between a descender line and an ascender line. Same-content lines follow
    roughly (that line's ink height + 0.12 em). PITCH_EM below is a single
    constant at the TOP of the observed range, so multi-line blocks reserve
    slightly more room than they need rather than less.

Rollback / A-B testing:
    TEXT_METRICS_FITTED=1   force the old fitted path (also restores the old
                            PITCH_EM, so a whole run reverts together)
    CHEMDRAW_FONT=/path     override regular-face autodetection
"""

import os
import sys
from functools import lru_cache

# ----------------------------------------------------------------------
# Original fitted constants. Kept verbatim: other modules may import them,
# and the fitted path below is still the fallback when the font is absent.
# ----------------------------------------------------------------------
INTERCEPT_EM = -0.10861576113418797
H0, HP, HD, HC = 0.5364094463242136, 0.19933708630241334, 0.04586395882463665, 0.18856138572745962
SUBSCRIPT_W = 0.7489

# Characters ChemDraw's auto-formula mode (face="96") renders at SUBSCRIPT_W.
# Digits were expected; the charge signs were forced by the data. 'N+' was
# mispredicted by a constant -0.1489 em, and 0.25 * advance('+') = 0.1460 --
# i.e. the '+' was being measured full size when ChemDraw had shrunk it.
# Both the ASCII hyphen and the Unicode minus occur in real labels, and the
# superscript forms are included because a caller may paste them in.
FORMULA_SHRINK_CHARS = frozenset("0123456789+-−⁺⁻")


def should_shrink(char, prev_char):
    """Does `char` shrink, given the character before it, inside a face=96 run?

    Positional, which is the whole reason a character SET cannot express it:
    the '2' in 'NH2' shrinks and the '2' in 'HBr (2 eq)' does not. Requiring a
    preceding letter, digit or ')' took the residual RMS on face=96 runs from
    0.0176 em to 0.0032 em; adding the charge signs to the set above then took
    the count of residuals over 0.03 em from 109 to zero.

    prev_char is None (or "") for the first character of a run. A leading
    character has no predecessor, so it never shrinks -- which is exactly the
    'HBr (2 eq)' case once the run is split on the space.
    """
    if char not in FORMULA_SHRINK_CHARS:
        return False
    if not prev_char:
        return False
    return prev_char.isalpha() or prev_char.isdigit() or prev_char == ")"

# Measured from ChemDraw, not from a font metric: the caption line pitch
# matches neither the hhea line box (1.1499) nor the OS/2 typo box. The hhea
# line box would over-reserve 7.9% per line.
PITCH_EM_MEASURED = 1.060413
PITCH_EM_FITTED = 1.1516958346149935
PITCH_EM = (PITCH_EM_FITTED if os.environ.get("TEXT_METRICS_FITTED")
            else PITCH_EM_MEASURED)

# Subscripts sit exactly one descender below the baseline: the fitted 0.2120
# matches Arial's hhea descent 0.2119.
SUB_DROP_EM = 0.212

ADVANCE_EM = {" ": 0.26231485831965984, "%": 0.8902490125473141, "'": 0.20861576113418656, "(": 0.3562490125473119, ")": 0.35624901254731167, "+": 0.5806157611341888, ",": 0.21461576113419034, "-": 0.37861576113419027, ".": 0.19466154790715479, "/": 0.3866157611341856, "0": 0.5628263620918337, "1": 0.44562636209183437, "2": 0.561388962875861, "3": 0.5636815915435246, "4": 0.6036157611341868, "5": 0.5826157611341886, "6": 0.5806157611341887, "7": 0.5726157611341889, "8": 0.5796157611341914, "9": 0.5786157611341887, "=": 0.5806157611341873, "A": 0.7704993416982067, "B": 0.6496157611341907, "C": 0.7343431761028966, "D": 0.9391411516388746, "E": 0.6286815915435248, "H": 0.68224901254731, "I": 0.20361576113418706, "M": 0.8192052537113961, "N": 0.6586815915435296, "O": 0.7832172759241244, "P": 0.6556157611341882, "R": 0.8059341695906598, "S": 0.671343176102893, "T": 0.6264395479729813, "W": 0.9722385870447295, "X": 0.7656157611341889, "a": 0.5726615479071561, "c": 0.5582947993202788, "d": 0.34717345432097413, "e": 0.5679911819357603, "f": 0.35195893723708055, "i": 0.2134615107247564, "l": 0.21337353293748323, "m": 0.8232490125473123, "n": 0.5462628642992792, "o": 0.5969367229480971, "p": 0.5592756249759648, "q": 0.6265820504207562, "r": 0.371866090285084, "t": 0.3337273783164941, "u": 0.7056271752530026, "°": 0.37961576113419004, "·": 0.20861576113418867}
DEFAULT_ADVANCE = 0.5

# ----------------------------------------------------------------------
# Font faces
# ----------------------------------------------------------------------

# (bold, italic) -> candidate paths, best first
FONT_CANDIDATES = {
    (False, False): ["/System/Library/Fonts/Supplemental/Arial.ttf",
                     "/Library/Fonts/Arial.ttf",
                     os.path.expanduser("~/Library/Fonts/Arial.ttf"),
                     r"C:\Windows\Fonts\arial.ttf"],
    (True, False):  ["/System/Library/Fonts/Supplemental/Arial Bold.ttf",
                     "/Library/Fonts/Arial Bold.ttf",
                     os.path.expanduser("~/Library/Fonts/Arial Bold.ttf"),
                     r"C:\Windows\Fonts\arialbd.ttf"],
    (False, True):  ["/System/Library/Fonts/Supplemental/Arial Italic.ttf",
                     "/Library/Fonts/Arial Italic.ttf",
                     os.path.expanduser("~/Library/Fonts/Arial Italic.ttf"),
                     r"C:\Windows\Fonts\ariali.ttf"],
    (True, True):   ["/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf",
                     "/Library/Fonts/Arial Bold Italic.ttf",
                     os.path.expanduser("~/Library/Fonts/Arial Bold Italic.ttf"),
                     r"C:\Windows\Fonts\arialbi.ttf"],
}

_FACES = {}          # (bold, italic) -> (hb, font, upem, path) or False
_WARNED = set()


def _warn(msg):
    if msg not in _WARNED:
        print("text_metrics: " + msg, file=sys.stderr)
        _WARNED.add(msg)


def _face(bold=False, italic=False):
    """Load and cache one face. Falls back to regular, then to the fitted model."""
    key = (bool(bold), bool(italic))
    if key in _FACES:
        return _FACES[key] or None

    if os.environ.get("TEXT_METRICS_FITTED"):
        _FACES[key] = False
        return None
    try:
        import uharfbuzz as hb
    except ImportError:
        _warn("uharfbuzz missing, using fitted widths (pip install uharfbuzz)")
        _FACES[key] = False
        return None

    path = None
    if key == (False, False):
        env = os.environ.get("CHEMDRAW_FONT")
        if env and os.path.exists(env):
            path = env
    if path is None:
        path = next((p for p in FONT_CANDIDATES[key] if os.path.exists(p)), None)

    if path is None:
        if key == (False, False):
            _warn("Arial not found, using fitted widths (set CHEMDRAW_FONT)")
            _FACES[key] = False
            return None
        # A styled face is missing but the regular one may not be. Measuring a
        # bold string with the regular table costs up to 6.7 pt, so say so
        # once rather than drifting silently.
        _warn("no font file for bold=%s italic=%s; measuring it with the "
              "regular face (up to ~6.7 pt too narrow)" % key)
        got = _face(False, False)
        _FACES[key] = got or False
        return got

    face = hb.Face(open(path, "rb").read())
    got = (hb, hb.Font(face), face.upem, path)
    _FACES[key] = got
    return got


# ----------------------------------------------------------------------
# Exact path
# ----------------------------------------------------------------------

def _shrunk_flags(text, subscript_chars, formula):
    """Per character: does it render at SUBSCRIPT_W of the nominal size?

    Two different questions, deliberately kept apart:

      subscript_chars   the CALLER asserts these characters are subscript.
                        Set membership, position-independent. This is right for
                        text the caller has already marked up (~x~ in a caption):
                        the marker is a statement of intent, not something to
                        re-derive.

      formula=True      ChemDraw's auto-formula mode, i.e. what it does to any
                        run written with face="96" -- which is every atom and
                        superatom label this engine emits. The rule is
                        positional, so a set cannot express it.

    Passing both applies both; a character shrinks if either says so.
    """
    prev = None
    out = []
    for ch in text:
        flag = ch in subscript_chars
        if formula and not flag:
            flag = should_shrink(ch, prev)
        out.append(flag)
        prev = ch
    return out


def _runs(text, subscript_chars, formula=False):
    """Split into consecutive runs of same subscript-ness.

    Kerning applies only within a run, which is also what a renderer does
    across a style change. (Kerning is off by default -- see the module
    docstring -- but the run split still matters for the scale factor.)
    """
    if not text:
        return []
    flags = _shrunk_flags(text, subscript_chars, formula)
    out, cur, flag = [], text[0], flags[0]
    for c, f in zip(text[1:], flags[1:]):
        if f == flag:
            cur += c
        else:
            out.append((cur, flag))
            cur, flag = c, f
    out.append((cur, flag))
    return out


@lru_cache(maxsize=8192)
def _exact_width(text, size, subscript_chars, bold, italic, trim=True,
                 formula=False):
    """Width in points, or None if the exact path is unavailable.

    trim=False -> the ADVANCE box. THE DEFAULT, and what ChemDraw records in
                  the overwhelming majority of its own files. Also the right
                  choice when summing several runs into one line: ink boxes
                  lose a pair of bearings at every style join, and they lose
                  it in the narrow direction.
    trim=True  -> the INK box, outer side bearings removed.

    WHICH BOX CHEMDRAW ACTUALLY WRITES (measured over the 652 CDXML
    files shipped in ChemDraw 23.1.2 -- 6495 <t> objects, all with a
    BoundingBox, zero parse failures):

        advance box   435 of the 443 files that contain text   (98%)
        ink box         8 files                                 ( 2%)

    The ink box is the 2% case, and defaulting to it costs 60% of this
    module's total width error -- more than the subscript factor and the
    kerning setting combined.

    The choice is a DOCUMENT-level property, not a per-<t> one. Within a file
    every label follows the same convention, and no <t> attribute predicts it:
    LabelAlignment, LabelJustification, Justification, LabelSize,
    InterpretChemically, Z and the node's bond count all score at the
    all-one-class baseline (purity 0.871, Cramer's V <= 0.08), while `file`
    scores 0.999. It correlates with the writing build but is NOT monotone in
    version -- 23.1 writes advance and 23.5 writes ink; 25.0.2 writes advance
    and 25.0.0 writes ink -- so a version test is not a safe discriminator
    either. A reader that needs to know must decide per file, e.g. by measuring
    that file's own single-character atom labels against both references.
    """
    got = _face(bold, italic)
    if not got:
        return None
    hb, font, upem, _ = got

    total = 0.0
    first = last = None          # (glyph_id, scale) / (glyph_id, scale, adv)
    for chunk, is_sub in _runs(text, subscript_chars, formula):
        scale = SUBSCRIPT_W if is_sub else 1.0
        buf = hb.Buffer()
        buf.add_str(chunk)
        buf.guess_segment_properties()
        hb.shape(font, buf, {"kern": False})
        infos, poss = buf.glyph_infos, buf.glyph_positions
        if not infos:
            continue
        total += sum(p.x_advance for p in poss) / upem * size * scale
        if first is None:
            first = (infos[0].codepoint, scale)
        last = (infos[-1].codepoint, scale, poss[-1].x_advance)

    if first is None:
        return 0.0

    if not trim:
        return total

    # trim to the ink box: drop the outer side bearings
    try:
        g, scale = first
        e = font.get_glyph_extents(g)
        if e is not None:
            total -= e.x_bearing / upem * size * scale
        g, scale, adv = last
        e = font.get_glyph_extents(g)
        if e is not None:
            total -= (adv - (e.x_bearing + e.width)) / upem * size * scale
    except Exception:
        # extents unavailable: return the advance box instead. Always >= ink,
        # so layout over-reserves slightly rather than colliding.
        pass

    return total


@lru_cache(maxsize=4096)
def _vbounds(ch, bold, italic):
    """(yMin, yMax) of one glyph's ink, in em. None if the glyph has no ink."""
    got = _face(bold, italic)
    if not got:
        return None
    hb, font, upem, _ = got
    buf = hb.Buffer()
    buf.add_str(ch)
    buf.guess_segment_properties()
    hb.shape(font, buf, {"kern": False})
    if not buf.glyph_infos:
        return None
    try:
        e = font.get_glyph_extents(buf.glyph_infos[0].codepoint)
    except Exception:
        return None
    if e is None or not e.height:
        return None
    a = e.y_bearing / upem
    b = (e.y_bearing + e.height) / upem
    return (min(a, b), max(a, b))


@lru_cache(maxsize=8192)
def _exact_extremes(text, size, subscript_chars, bold, italic, formula=False):
    """(yMin, yMax) of one line's ink, in points relative to its baseline.

    Kept separate from the height because a multi-line block's box runs from
    the FIRST line's top to the LAST line's bottom -- taking max(height) over
    lines is a different quantity, and differs whenever the lines differ in
    ascenders and descenders.
    """
    if not _face(bold, italic):
        return None
    hi, lo = [], []
    flags = _shrunk_flags(text, subscript_chars, formula)
    for ch, shrunk in zip(text, flags):
        b = _vbounds(ch, bold, italic)
        if b is None:            # spaces and other blank glyphs
            continue
        if not shrunk:
            hi.append(b[1])
            lo.append(b[0])
        elif formula and ch in "+-\u2212\u207a\u207b" and ch not in subscript_chars:
            # A charge sign in formula mode is a SUPERscript: it shrinks but it
            # rises, it does not drop. The rise is still unmeasured (see the
            # module docstring), so apply the shrink and leave the baseline
            # alone -- over-reserving slightly rather than placing it too low.
            hi.append(b[1] * SUBSCRIPT_W)
            lo.append(b[0] * SUBSCRIPT_W)
        else:
            hi.append(b[1] * SUBSCRIPT_W - SUB_DROP_EM)
            lo.append(b[0] * SUBSCRIPT_W - SUB_DROP_EM)
    if not hi:
        return None
    return (min(lo) * size, max(hi) * size)


@lru_cache(maxsize=8192)
def _exact_height1(text, size, subscript_chars, bold, italic, formula=False):
    e = _exact_extremes(text, size, subscript_chars, bold, italic, formula)
    return None if e is None else (e[1] - e[0])


def _fitted_width(text, size, subscript_chars="", formula=False):
    """The original per-character fitted model. Retained as fallback."""
    total = INTERCEPT_EM
    for c, shrunk in zip(text, _shrunk_flags(text, subscript_chars, formula)):
        adv = ADVANCE_EM.get(c, DEFAULT_ADVANCE)
        total += adv * (SUBSCRIPT_W if shrunk else 1.0)
    return size * total


def _feats(text):
    return (1.0,
            1.0 if any(c in "()[]{}" for c in text) else 0.0,
            1.0 if any(c in "gjpqy,;()" for c in text) else 0.0,
            1.0 if any(c.isupper() or c.isdigit() for c in text) else 0.0)


def _fitted_height1(text, size):
    f = _feats(text)
    return size * (H0 * f[0] + HP * f[1] + HD * f[2] + HC * f[3])


# ----------------------------------------------------------------------
# Public API — original signatures, plus optional bold=/italic=
# ----------------------------------------------------------------------

def width(text, size, subscript_chars="", bold=False, italic=False,
          trim=False, formula=False):
    """Rendered ADVANCE width in points.

    subscript_chars: the characters in `text` rendered as subscript.
    bold / italic: optional, default off, so two-argument callers are
    unaffected. They matter: a bold string measured with the regular table
    comes out up to 6.7 pt narrow, an italic one up to 1.5 pt.

    trim=True switches to the ink box. It is not the default: of the 443
    ChemDraw sample files that contain text, 435 (98%) record the advance box.
    See _exact_width for which files are the exception and how to tell.
    """
    w = _exact_width(text, size, subscript_chars, bool(bold), bool(italic),
                     bool(trim), bool(formula))
    return w if w is not None else _fitted_width(text, size, subscript_chars,
                                                 bool(formula))


def height1(text, size, subscript_chars="", bold=False, italic=False,
            formula=False):
    """Rendered ink height of a single line, in points."""
    h = _exact_height1(text, size, subscript_chars, bool(bold), bool(italic),
                       bool(formula))
    return h if h is not None else _fitted_height1(text, size)


def height(lines, size, subscript_chars="", bold=False, italic=False,
           formula=False):
    """lines: str or list of strings. Total rendered height incl. line pitch.

    Each line sits one pitch below the one before, so line i's ink occupies
    [yMin_i - i*pitch, yMax_i - i*pitch]. The block spans the highest top and
    the lowest bottom across the lines that actually carry ink.

    Two things this gets right that a max(height1) formula cannot:

      * a descender in a MIDDLE line falls between baselines and does not
        deepen the block, while the same descender in the last line does;
      * a blank line carries no ink but still advances the baseline, so
        "a\\n\\nb" is three lines tall, and a trailing newline adds nothing.
    """
    if isinstance(lines, str):
        lines = lines.split("\n")
    lines = list(lines) or [""]
    pitch = size * PITCH_EM
    b, i_ = bool(bold), bool(italic)

    tops, bots = [], []
    for i, line in enumerate(lines):
        e = _exact_extremes(line, size, subscript_chars, b, i_, bool(formula))
        if e is None:
            continue                     # blank line: advances, adds no ink
        tops.append(e[1] - i * pitch)
        bots.append(e[0] - i * pitch)
    if tops:
        return max(tops) - min(bots)

    return (max(height1(l, size, subscript_chars, b, i_, formula) for l in lines)
            + (len(lines) - 1) * pitch)


def box(lines, size, subscript_chars="", bold=False, italic=False,
        formula=False):
    """(w, h) rendered box for a (multi-line) text block.

    A single string containing newlines is split, so passing either a list of
    lines or one '\\n'-joined string gives the same answer.
    """
    if isinstance(lines, str):
        lines = lines.split("\n")
    lines = list(lines) or [""]
    w = max(width(l, size, subscript_chars, bold, italic, formula=formula)
            for l in lines)
    return w, height(lines, size, subscript_chars, bold, italic, formula)


# ----------------------------------------------------------------------
# Diagnostics and guards
# ----------------------------------------------------------------------

def backend():
    """('exact', path) for the regular face, or ('fitted', reason)."""
    got = _face(False, False)
    if got:
        return ("exact", got[3])
    if os.environ.get("TEXT_METRICS_FITTED"):
        return ("fitted", "TEXT_METRICS_FITTED is set")
    return ("fitted", "uharfbuzz or Arial unavailable")


def faces():
    """{(bold, italic): path or None} -- which faces are actually measurable."""
    out = {}
    for key in FONT_CANDIDATES:
        got = _face(*key)
        out[key] = got[3] if got else None
    return out


def assert_font(name):
    """Fail loudly if the document asks for a font that is not being measured.

    The worst failure mode for this module is not being wrong -- it is being
    wrong SILENTLY and self-consistently, which is exactly what a fitted
    table does. If a document's fonttable names something other than the face
    loaded here, every width and height is wrong while every assertion still
    passes. Call this once per document with the fonttable's font name.
    """
    kind, detail = backend()
    if kind != "exact":
        raise RuntimeError(
            "text_metrics is on the fitted fallback (%s); measurements are "
            "estimates, not the font's own values" % detail)
    got = os.path.basename(detail)
    norm = lambda s: "".join(s.lower().split()).replace("-", "")
    if norm(name) not in norm(got):
        raise RuntimeError(
            "document uses font %r but text_metrics loaded %r -- every width "
            "and height would be silently wrong" % (name, got))


def missing_chars(text):
    """Characters with no glyph in the regular face.

    These are where any width model silently lies. Arial has no glyph for the
    equilibrium arrows U+21CC and U+21C4, for instance, so a caption using one
    both prints a placeholder and measures wrong.
    """
    got = _face(False, False)
    if not got:
        return [c for c in sorted(set(text)) if c not in ADVANCE_EM]
    hb, font, upem, _ = got
    out = []
    for c in sorted(set(text)):
        if not c.strip():
            continue
        buf = hb.Buffer()
        buf.add_str(c)
        buf.guess_segment_properties()
        hb.shape(font, buf, {"kern": False})
        if buf.glyph_infos and buf.glyph_infos[0].codepoint == 0:
            out.append(c)
    return out


if __name__ == "__main__":
    kind, detail = backend()
    print(f"backend: {kind}  ({detail})")
    print(f"PITCH_EM = {PITCH_EM:.6f}   SUB_DROP_EM = {SUB_DROP_EM}")
    print("\nfaces:")
    for (b, i), p in sorted(faces().items()):
        tag = ("bold " if b else "") + ("italic" if i else "") or "regular"
        print(f"  {tag:<12} {p or 'NOT FOUND'}")

    print(f"\n{'':<22}{'width':>22}{'height':>10}")
    print(f"{'':<22}{'regular':>8}{'bold':>7}{'ital':>7}{'exact':>10}")
    for s, sub in [("Pd(dppf)Cl2", "2"), ("K2CO3", "23"), ("35b", ""),
                   ("n-BuLi", ""), ("then", ""), ("gy", "")]:
        print(f"  {s!r:<20}{width(s, 10, sub):>8.2f}"
              f"{width(s, 10, sub, bold=True):>7.2f}"
              f"{width(s, 10, sub, italic=True):>7.2f}"
              f"{height1(s, 10, sub):>10.2f}")

    bad = missing_chars("\u03b1\u03bc\u00b0\u00b7\u2212\u2192\u21cc\u21c4\u2264")
    print(f"\nno glyph in Arial: {bad if bad else 'none of the probed set'}")
