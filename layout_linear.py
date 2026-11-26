"""
layout_linear.py — analytical linear reaction-scheme layout (render once).

Follows the earlier JavaScript implementation for LINEAR routes, but every position
is COMPUTED from text_metrics instead of discovered by rendering + reading back.
Arrows are sized to their widest caption; captions, yields, compound numbers, salt
labels and optional reagent structures are placed with exact clearance; multi-row
schemes are centered to a common axis. Output is a finalized CDXML string (positive
quadrant, tight page) ready to open once in desktop ChemDraw.

Recipe:
  {
    "structures": { id: (fragment_str, bbox), ... },   # bbox = [x0,y0,x1,y1]
    "rows": [
      [ {"s": id, "label": "2HCl"?, "under": id?},   # "under" x-aligns this structure
                                                     # under one in an EARLIER row, so a
                                                     # connector between them is VERTICAL
        {"arrow": {"above":[lines]?, "below":[lines]?, "yield":str?,
                   "reagent": id?, "reagentEq": str?}},
        {"s": id}, ... ],
      ...
    ],
    "footnotes": ["^a^ 0.2 mmol scale.", ...]?    # under the whole figure, wrapped
                                                  # to its width, left-aligned
  }
Caption line strings may use ~x~ for subscripts (e.g. "CuSO~4~\u00b75H~2~O").
"""

from __future__ import annotations

import re

import cdxml_core as C
import text_metrics as M
import arrows as A
import nested as _N
import rowspan as _RS

# Arial cap height / ascent as a fraction of point size. The <t> anchor is the
# FIRST LINE's baseline; the box top sits ASCENT_EM*size above it. Measured
# constant from the calibration readback (top_off/size = 0.716 across all blocks).
ASCENT_EM = 0.716

CFG = dict(
    row_gap_y=130.0,     # MINIMUM baseline-to-baseline row spacing (short rows keep this)
    row_clear=28.0,      # clearance between adjacent rows' ACTUAL content extents
    gap=30.0,            # horizontal gap between a structure and an adjacent arrow
    arrow_pad=12.0,      # arrow extends this far beyond the widest caption on each side
    min_arrow=60.0,
    max_arrow=None,      # optional cap; when set, an over-wide reagent is shrunk instead
    txt_gap_above=6.0,   # caption BOTTOM sits this far above the arrow line
    txt_gap_below=8.0,   # caption TOP sits this far below the arrow line
    txt_gap_side=4.0,    # VERTICAL arrow: caption edge sits this far to either SIDE of
                         # the shaft. One symmetric value - txt_gap_above/below are an
                         # asymmetric pair tuned for text sitting over/under a shaft,
                         # which is the wrong axis once the arrow turns vertical.
    yield_gap=5.0,       # gap between the below-caption and a yield line
    num_gap=10.0,        # structure bottom -> compound-number TOP
    num_size=11,
    # Chords per attached arc when the ledger reserves one. NOT a new number:
    # layout_cycle's cfg["arc_samples"] is 48 for the same job (ring arcs and
    # feeders through arc_points -> reserve_polyline). Read off the existing answer
    # rather than picked -- and it only affects what is RESERVED, never what is drawn.
    arc_samples=48,
    cap_size=9,
    label_size=10,
    label_dx=12.0,       # salt label offset to the right of a structure
    group_salt=False,    # True -> emit compound+salt as one ChemDraw <group> (Ctrl+G).
                         # Confirm the <group> tag against a manual Ctrl+G export first.
    reagent_gap=14.0,    # arrow line -> reagent structure BOTTOM
    reagent_num_gap=4.0,
    sublabel_gap=3.0,    # gap between stacked labels (number -> dr -> ee ...)
    sublabel_size=9,
    plus_size=12,        # size of the "+" between products
    plus_gap=14.0,       # space on each side of a "+"
    # A species written rather than drawn. NOT plus_size -- a "+" is a
    # connective and rides small, while "H-C(sp3)" is a chemical species standing
    # in a structure's slot and reads at the weight of the things beside it.
    text_cell_size=11,
    text_cell_pad=6.0,   # breathing room each side, so type does not touch a bond
    arrow_margin=16.0,   # gap between a routed connector arrow and the boxes it links
    row_align="left",    # multi-row alignment: "left" reads like wrapped text (reading
                         # order); "center" = common-axis centering
    margin=30.0,
)

# ---- footnotes: NOT in CFG, on purpose --------------------------------------
# Merged into cfg at the one call site that draws them, so a figure with no
# footnotes carries no footnote keys and CFG stays the set of things every linear
# figure has. The four values are copied VERBATIM from layout_grid.GRID_DEFAULTS
# (footnote_size / _gap / _lead / _width) so a footnote under a scheme and a
# footnote under a scope table are the same object at the same size in the same
# place. If one of them is ever re-measured, both move -- change them in both
# dicts or say in the commit why they diverged.
FOOTNOTE_DEFAULTS = {
    "footnote_size": 8,
    "footnote_gap": 22.0,
    "footnote_lead": 3.0,
    "footnote_width": None,   # None -> wrap to the figure's own width
}

# ---- ~x~-aware text measuring (uses the fitted advance table) ----
# THE SAME OBJECT the writer splits on, not a copy of it. A copy -- the identical
# regex written out twice, here and in cdxml_core -- goes stale when a marker is
# added, and that costs: {red}94{} measured 53 pt against the 28 pt it draws,
# because the measuring side counted the markup as letters. A slot
# reserved at nearly twice the width of its content, with a gap on the page and
# nothing to say why.
#
# One grammar, one definition. A marker added to the writer is understood here
# the moment it exists; the only thing left to decide below is whether it changes
# the ADVANCE (bold and italic do, colour and underline do not).
_MARK = C._STYLE_MARK


# --- shadow allocator ------------------------------------------------------
# Enable with CDXML_SHADOW=1; report is written to /tmp/cdxml_shadow.txt at
# exit. Disabled by default: nothing below runs unless the env var is set.
import os as _os
import atexit as _atexit
try:
    from allocator import (Allocator as _Allocator, Rect as _Rect,
                           polyline_hull as _polyline_hull,
                           run_shapes as _run_shapes,
                           ARROW_INK_WIDTH as _INK)
    _SHADOW = bool(_os.environ.get("CDXML_SHADOW"))
except Exception:
    _Allocator = _Rect = _polyline_hull = _run_shapes = None
    _INK = 4.0
    _SHADOW = False

_SHADOW_LOG = []

# [(reagent id, scale factor)] for reagents whose decoration had to be dropped
# because the reagent itself was shrunk to fit its arrow. Read and cleared by the
# tool layer, the same way layout_grid's _WIDEST_LEGEND is -- a list rather than a
# raise, because the figure is still worth opening with one mark missing.
DECOR_UNSCALED: list = []

# [(owner, direction, distance, how)] for every connector caption the placer had
# to move off the spot the drawing convention asks for.
#
# WHY THIS EXISTS. A silent search costs information: without the search, a
# caption sitting on a structure is REPORTED as an overlap, so the figure's real
# defect is visible. With the search, the caption dodges and
# nothing says why it is now three line-heights from its arrow -- so the reader
# sees an oddly floating caption and no explanation. Avoiding a collision is an
# improvement; doing it quietly moves the figure away from the convention this
# engine documents, and that is a thing the caller is entitled to know.
CAPTION_MOVED: list = []

# [(owner, corridor_x, [what was in the way])] for every connector whose SHAFT had
# to leave the straight line to miss something.
#
# WHY THE SHAFT NEEDS THIS AT ALL. The captions have free_near (see the note
# above); without a detour, a connector spanning non-adjacent rows is drawn
# straight through whatever sits between them -- measured on a convergent scheme,
# where the A2 -> C arrow ran down the middle of the azide in row 2, and the engine
# reported the collision and drew it anyway. The connector-caption note in
# render_linear says why row spacing cannot fix that one ("a connector from row 0
# to row 3 runs vertically THROUGH rows 1 and 2"): the answer has to be the route.
#
# REPORTED for exactly the reason CAPTION_MOVED is. A detour is an improvement, but
# it is also a departure from the straight line this engine's docs promise, and the
# reader is owed the sentence "it bends because B2 was in the way" rather than an
# arrow that inexplicably swings wide.
CONN_DETOURED: list = []

# [(from, to, wanted, got)] for every `span` the stretch below could not reach.
#
# WHY IT EXISTS. A span is the one request in this engine whose result is
# invisible in the output: the figure is drawn either way, and the only difference
# between "308.7 as asked" and "181.2, because the arrows had nothing to give" is a
# number nobody can read off the page. The caller is a PARENT that has already
# anchored one end of this chain somewhere else, so a span that quietly came up
# short does not look like a short span -- it looks like the other engine put its
# structure in the wrong place.
#
# Two ways to come up short, both here: the distance asked for is under the natural
# one (an arrow is never shrunk -- its length is set by its own caption, and a
# shorter arrow prints the conditions past the head), or `max_arrow` capped the
# growth. Neither is refusable in preflight, which has no widths.
SPAN_SHORT: list = []

# What `wrap_width` did to the rows, for the reply -- so it does not say a blanket
# "Auto-wrapped to page width." whenever wrap is on, including for rows it has not
# touched. Filled by render_linear's wrap block, read and cleared by the
# tool layer like SPAN_SHORT. A row is numbered as the CALLER wrote it (None for
# one flat row), because those are the only numbers a caller can act on.
#   ("broke", row, n)                row became n rows
#   ("whole", row, width, why)       wider than the column and NOT broken; `why`
#                                    names the span that forbids it
#   ("over", row, piece, width, m)   piece `piece` of a broken row is still wider
#                                    than the column; m is the id of the merge
#                                    cell it holds, or None
#   ("snake", n_rows)                wrap_style='snake' on several rows: nothing
#                                    was broken
WRAP_ROWS: list = []

# --- the two accounts, and the check that holds them together ----------------
#
# This engine keeps TWO records of where things are. `boxes` (a plain list, whose
# union is what finalize is handed as the page extent) and the Allocator (whose
# bbox() is what grid, radial, snaking and cycle all hand to finalize instead).
# Every put_box writes to both -- except the connector block, which writes to
# `boxes` with reserve=False and then reserves a DIFFERENT rect: the polyline's
# chunks for a straight connector, the routed hull for an elbow.
#
# That is the same object counted twice, so the two need a consistency check.
#
# THIS IS A MEASUREMENT, NOT A FIX. The page extent comes from `boxes`, so the
# comparison never changes the output. Whether a disagreement means "boxes is
# reserving air" or "alloc is under-reserving, add to it" is decided from the
# list below, not guessed. The tempting inference is "alloc.bbox() is smaller,
# so alloc lost something", and the inverse is just as available: `boxes` may be
# reserving air.
#
# [(descriptor, boxes_bbox, alloc_bbox, per-side deltas)] for this batch. Recorded
# ALWAYS -- one tuple per figure, and the tool layer clears it the way it clears
# DECOR_UNSCALED -- so the agreement is a standing invariant a test can assert
# rather than something only an env var can see. CDXML_ACCOUNTS=1 additionally
# writes the table at exit.
ACCOUNT_DIFF: list = []
_ACCOUNTS = bool(_os.environ.get("CDXML_ACCOUNTS"))


def _account_diff(recipe, boxes, content_bbox, alloc) -> None:
    """Compare the two accounts. Records only; never changes what is drawn."""
    if alloc is None:
        return
    ab = alloc.bbox()
    if ab is None:
        return
    al = (ab.x, ab.y, ab.right, ab.bottom)
    # POSITIVE = alloc reaches further out on that side. Signs are per-side rather
    # than absolute so "alloc is smaller everywhere" and "alloc is bigger on one
    # side" cannot be read as the same finding.
    d = (content_bbox[0] - al[0], content_bbox[1] - al[1],
         al[2] - content_bbox[2], al[3] - content_bbox[3])
    rows = recipe.get("rows") or []
    desc = "rows=%d cells=%d arrows=%d conn=%d%s" % (
        len(rows),
        sum(1 for r in rows for c in r if "s" in c or "text" in c),
        sum(1 for r in rows for c in r if "arrow" in c),
        len(recipe.get("connectors") or []),
        "".join(" " + k for k in ("decorations",) if recipe.get(k)))
    ACCOUNT_DIFF.append((desc, tuple(round(v, 2) for v in content_bbox),
                         tuple(round(v, 2) for v in al),
                         tuple(round(v, 2) for v in d)))


def _accounts_dump():
    if not (_ACCOUNTS and ACCOUNT_DIFF):
        return
    path = _os.environ.get("CDXML_ACCOUNTS_OUT", "/tmp/cdxml_accounts.txt")
    agree = sum(1 for _d, _b, _a, dd in ACCOUNT_DIFF if not any(abs(v) > 0.005 for v in dd))
    try:
        with open(path, "w") as fh:
            fh.write("layout_linear: boxes-union vs alloc.bbox(), %d figures, "
                     "%d agree exactly\n" % (len(ACCOUNT_DIFF), agree))
            fh.write("delta signs: POSITIVE = alloc reaches FURTHER OUT on that side"
                     " (left, top, right, bottom)\n\n")
            for desc, b, a, d in ACCOUNT_DIFF:
                fh.write("%-46s boxes=%-34s alloc=%-34s d=%s\n"
                         % (desc, b, a, d))
    except Exception:
        pass


_atexit.register(_accounts_dump)


def _shadow_dump():
    if not _SHADOW_LOG:
        return
    path = _os.environ.get("CDXML_SHADOW_OUT", "/tmp/cdxml_shadow.txt")
    try:
        with open(path, "w") as fh:
            for i, al in enumerate(_SHADOW_LOG, 1):
                fh.write("===== figure %d =====\n%s\n\n" % (i, al.report()))
    except Exception:
        pass


_atexit.register(_shadow_dump)
# --- end shadow allocator --------------------------------------------------


def _runs(line: str, marks=frozenset()):
    """yield (text, marks); marks a subset of {'sub','sup','bold','italic'}.

    A SET rather than one name because the markers nest: `#{red}3a{}#` is bold
    and red, `*n*Bu~4~N` puts a subscript inside nothing but `#H~2~O#` puts one
    inside a bold span. A single name would measure the nested case with the
    INNER markers still in the string -- `#{blue}Table 2{}#`
    came out 70.6 pt against the 34.5 pt ChemDraw draws, because "{blue}" and the
    two "#" were counted as seven glyphs. Colour and underline add no mark at
    all: they change no advance (measured identical to 0.00 pt), and what matters
    is only that their MARKUP is not measured.

    Bold and italic are NOT folded into 'plain', though they look to have
    about the same advances. Measured against ChemDraw they do not: a
    bold string comes out up to 6.7 pt wider than the regular table predicts,
    an italic one up to 1.5 pt. Underline does not change the box at all
    (measured: identical to within 0.00 pt), so it stays folded in.

    Anything that is not a string is a caller mistake, and it gets named here
    because this is the one funnel every text field on every tool flows
    through -- one guard covers all six engines and 218 call sites. Raw, the
    failure is `_MARK.split` saying "expected string or bytes-like object, got
    'dict'", which names neither the field that was filled in nor the shape it
    wanted. {"text": ...} is the commonest wrong guess: that IS the shape a
    cycle node takes, so reusing it for `center=` reads as the obvious thing to
    write. Instrumented over 119 stress figures the guard is a no-op for valid
    input -- 3925 of 3944 calls arrive as `str`, and every one of the 19 that
    did not came from a deliberately malformed recipe that raises.
    """
    if not isinstance(line, str):
        raise ValueError(
            "a text line has to be a string, got %s (%.120r). Every text "
            "field on every tool takes either ONE string or a LIST of "
            "strings, one per line -- e.g. center=\"Pd catalysis\" or "
            "above=[\"NaBH~4~\", \"MeOH, 0 \u00b0C\"]. A dict like "
            "{\"text\": ...} is the shape a NODE takes, not a line of text."
            % (type(line).__name__, line))
    for p in _MARK.split(line):
        if p == "":
            continue
        if p.startswith("{") and p.endswith("{}"):
            yield from _runs(p[1:].split("}", 1)[1][:-2], marks)
        elif p.startswith("~") and p.endswith("~"):
            yield from _runs(p[1:-1], marks | {"sub"})
        elif p.startswith("^") and p.endswith("^"):
            yield from _runs(p[1:-1], marks | {"sup"})
        elif len(p) >= 2 and p[0] == "*" and p[-1] == "*":
            yield from _runs(p[1:-1], marks | {"italic"})
        elif len(p) >= 2 and p[0] == "#" and p[-1] == "#":
            yield from _runs(p[1:-1], marks | {"bold"})
        elif len(p) >= 2 and p[0] == "_" and p[-1] == "_":
            yield from _runs(p[1:-1], marks)   # underline: no effect on the box
        else:
            # Un-double a brace exactly as the writer does, or an escaped one
            # measures two characters wide and draws one. Every recursion ends
            # here: a span's inside is two characters shorter than the span.
            yield (p.replace("{{", "{").replace("}}", "}"), marks)


def _line_width_legacy(line: str, size: float) -> float:
    """The inlined fitted table. Kept for LAYOUT_LEGACY_METRICS=1 comparisons."""
    total = M.INTERCEPT_EM
    for txt, kind in _runs(line):
        scale = (M.SUBSCRIPT_W if "sub" in kind
                 else (size - 2) / size if "sup" in kind else 1.0)
        for c in txt:
            total += M.ADVANCE_EM.get(c, M.DEFAULT_ADVANCE) * scale
    return size * total


def _line_width(line: str, size: float, bold: bool = False) -> float:
    """Rendered width of one ~x~-marked line, in points.

    Asks text_metrics rather than re-implementing it. Reading M.ADVANCE_EM
    directly would stop every improvement to the model at the module boundary,
    where no figure would ever see it.

    trim=False on each run: run widths are SUMMED, and ink boxes drop a pair
    of side bearings at every join. Advance boxes add up correctly and are
    never short, which is the direction a reservation should err in.

    `bold` is the BLOCK's weight, for a caller that sets a whole line bold
    rather than marking a run with #...#: a table's header band is bold by
    style, and its text may already carry markup of its own, so wrapping the
    string in #...# to get the weight would swallow the markup inside it.
    Default False: a caller that does not ask measures with the regular table.
    """
    if _os.environ.get("LAYOUT_LEGACY_METRICS"):
        return _line_width_legacy(line, size)
    total = 0.0
    for txt, kind in _runs(line):
        if "sub" in kind or "sup" in kind:
            total += M.width(txt, size * M.SUBSCRIPT_W,
                             bold=(bold or "bold" in kind), trim=False)
        else:
            total += M.width(txt, size, bold=(bold or "bold" in kind),
                             italic=("italic" in kind), trim=False)
    return total


def _plain(line: str) -> str:
    return "".join(t for t, _ in _runs(line))


def block_box(lines, size: float, bold: bool = False):
    """(w, h) rendered box for a (multi-line) ~x~-marked text block.

    `bold` sets the whole block bold; see _line_width. It matters rather than
    being cosmetic -- a bold string measured with the regular table comes out
    up to 6.7 pt narrow, and a table's column width is often set by its header.
    """
    if isinstance(lines, str):
        lines = [lines]
    w = max(_line_width(l, size, bold) for l in lines)
    if _os.environ.get("LAYOUT_LEGACY_METRICS"):
        h = (max(M.height1(_plain(l), size) for l in lines)
             + (len(lines) - 1) * size * M.PITCH_EM)
    else:
        # THE ONE HEIGHT ROUTINE, for the same reason fold_lines is the one
        # folding routine. A formula copied here would take max(height1) over
        # the lines, which is a different quantity from the block's ink: a
        # descender in a middle line sits between baselines and does not
        # deepen the block, while one in the last line does.
        h = M.height([_plain(l) for l in lines], size, bold=bold)
    return w, h


def _text_cell_box(cell, cfg):
    """(w, h) of a written-out species cell, padding included.

    ONE measurer, called from three places -- the wrap check, the item builder
    and, through the item's stored w/h, the renderer. The `plus` it sits next to
    measures itself twice in two of those and gets away with it because "+" is a
    constant; a caller's string is not, and two spellings of the same measurement
    is how a cell ends up drawn wider than the slot reserved for it.

    The padding is on the WIDTH only. Vertically a species is centred on the row
    like a structure, and the row's height already comes from the tallest thing
    in it.
    """
    size = cell.get("size", cfg["text_cell_size"])
    w, h = block_box([str(cell["text"])], size)
    return w + 2 * cfg["text_cell_pad"], h


def _words(line: str) -> list[str]:
    """`line` cut at the spaces BETWEEN words, with ~sub~ ^sup^ *it* #bold# spans
    kept whole.

    Atomic spans because a fold inside `~2~` would strand a tilde at the end of one
    line and open an unterminated subscript on the next -- ChemDraw would print both
    marks. `*t*-BuOK` and `Cs~2~CO~3~,` therefore come back as ONE word each: the
    marker pieces glue onto whatever plain text they touch, and only a real space
    starts a new word.
    """
    words = [""]
    for p in _MARK.split(line):
        if p == "":
            continue
        if _MARK.fullmatch(p):
            words[-1] += p
            continue
        bits = p.split(" ")
        words[-1] += bits[0]
        words.extend(bits[1:])
    return [w for w in words if w]


def fold_lines(lines, size: float, cap: float) -> list[str]:
    """`lines` re-broken at word boundaries so none renders wider than `cap` points.

    THE ONE FOLDING ROUTINE for this pipeline. Whoever reserves space for a text
    block and whoever writes that block MUST fold it the same way, or the reserved
    rectangle is not the rectangle that lands on the page -- so both go through
    here, and in practice through fold_box, which returns the lines and their box
    together so that a caller cannot take one without the other.

    A line that already fits is passed through UNTOUCHED, before any tokenising
    happens. That is what keeps a figure whose captions all fit byte-identical to
    one drawn without folding: no round trip through _words can normalise a
    double space or an odd marker.

    A single word wider than `cap` gets a line of its own rather than being
    hyphenated. The box is then measured from that line, so an unfoldable label
    makes the figure honestly too wide instead of quietly disagreeing with its own
    reservation -- which is the failure this routine exists to remove.
    """
    if isinstance(lines, str):
        lines = [lines]
    out: list[str] = []
    for line in lines:
        if _line_width(line, size) <= cap:
            out.append(line)
            continue
        cur = ""
        for w in _words(line):
            trial = "%s %s" % (cur, w) if cur else w
            if cur and _line_width(trial, size) > cap:
                out.append(cur)
                cur = w
            else:
                cur = trial
        if cur:
            out.append(cur)
    return out


def block_leading(lines, size: float) -> float:
    """The blank space one line of this block carries: pitch less line height.

    This is the gap ChemDraw itself leaves between two lines of the SAME text block,
    which makes it the natural floor for how close two SEPARATE blocks may sit -- any
    closer and a reader takes them for one paragraph. 3.43 units at size 9, and it
    needs no constant of its own because the type size already implies it.

    Independent of the line count: n lines of pitch less a box that is one line plus
    (n-1) pitches leaves exactly one line's leading, whatever n is.
    """
    if isinstance(lines, str):
        lines = [lines]
    _w, h = block_box(lines, size)
    return size * M.PITCH_EM * len(lines) - h


def fold_box(lines, size: float, cap: float):
    """(folded lines, (w, h)) -- ask for both at once and they cannot disagree.

    The height is the FOLDED height: n lines of pitch, not one. Callers that turn a
    box into a radius (a cycle's `reach`, a feeder's disk) get the taller box, which
    is the whole point -- folding a label trades width for height, and a packer told
    only about the width would pack the figure on the unfolded number.
    """
    folded = fold_lines(lines, size, cap)
    return folded, block_box(folded, size)


def render_linear(recipe: dict, header_prefix: str, cfg: dict | None = None,
                  page_mode: str = "fit", wrap_width: float | None = None,
                  serpentine: bool = False, snake_rows: set | None = None,
                  alloc=None, shift: list | None = None, frame=None,
                  placement: dict | None = None) -> str:
    cfg = {**CFG, **(cfg or {})}
    # ARE THESE ROWS A SNAKE? Two ways to be one, and only these two: the wrap
    # below BUILT one (it knows, because it chose the breaks and did the
    # reversing), or a caller handed rows that already alternate and said so.
    #
    # NOT "the caller asked for wrap_style='snake'". A caller who splits their own
    # rows has said where each one starts, and their rows read forwards -- hanging
    # every other one at the far edge, on the strength of an alternation that is
    # not there, moved a left-to-right row 123 pt to the right.
    #
    # AND IT CARRIES WHICH ROWS READ BACKWARDS, as a set of row indices, rather
    # than "the odd ones". A wrap builds strict alternation, but a caller placing
    # rows by hand may write two that both read forwards -- a staircase -- and
    # "every other row hangs at the far edge" is false for that figure.
    _snaked = snake_rows is not None
    _back_rows = set(snake_rows or ())
    # Compose any {"stack": [...]} cell into a synthetic S entry before
    # the packer sees it. IDENTITY (same object back) when no cell uses one, so
    # a recipe without stacks cannot take a different path -- see nested.py.
    recipe = _N.expand_stacks(recipe, cfg)
    if alloc is None and _SHADOW:
        alloc = _Allocator()
    # AND _SHADOW, or this leaks. The list is read by exactly one thing -- the
    # atexit dump below, which is off unless CDXML_SHADOW is set -- and gated only
    # on "did anyone pass an allocator", the append would always run, because the
    # tool layer always passes one. A long-running server would then keep one
    # Allocator, with every rect of every figure in it, for every figure it had
    # ever drawn. Measured: 5 draws, 5 kept.
    if alloc is not None and _SHADOW:
        _SHADOW_LOG.append(alloc)
    S = recipe["structures"]

    def cap_w(lines, size):
        return block_box(lines, size)[0] if lines else 0.0

    # ---------- pass 1: size arrows, lay out each row's x, find centers ----------
    def measure_arrow(a):
        widths = []
        for k in ("above", "below"):
            if a.get(k):
                widths.append(cap_w(a[k], cfg["cap_size"]))
        if a.get("yield"):
            widths.append(cap_w([a["yield"]], cfg["cap_size"]))
        if a.get("reagentEq") and not a.get("reagent"):
            widths.append(cap_w([a["reagentEq"]], cfg["cap_size"]))
        if a.get("reagent"):
            rb = S[a["reagent"]][1]
            _rwid = rb[2] - rb[0]
            if a.get("reagentEq"):
                # the equivalents caption sits BESIDE the reagent, so the arrow
                # must span both plus the gap -- not the wider of the two, which
                # would let the caption run past the arrow into the next cell
                _rwid += (cfg.get("reagent_eq_gap", 8.0)
                          + cap_w([a["reagentEq"]], cfg["cap_size"]))
            widths.append(_rwid)
        widest = max(widths) if widths else 0.0
        length = max(widest + 2 * cfg["arrow_pad"], cfg["min_arrow"])
        if cfg.get("max_arrow"):
            length = min(length, cfg["max_arrow"])
        return length

    # ---------- optional: auto-wrap rows to a target usable width ----------
    # Straight arrows only, reading order. The overflowing cell picks the break style:
    # a STRUCTURE that won't fit means the arrow before it stayed on this row (trailing)
    # and its product starts the next row; an ARROW that won't fit LEADS the next row.
    # A "+" never starts a row. (No effect when wrap_width is None.)
    #
    # SEVERAL ROWS: every row wider than the column is broken by that same rule and
    # its pieces take its place, in order; a row that fits is left exactly as the
    # caller wrote it. Rows split upstream by COUNTING cells need this, because only
    # this layer knows how wide a cell comes out: such rows measured 656-867 pt on a
    # 552 pt column. Three things differ from one flat row:
    #
    #   * A ROWSPAN IS READ OFF ITS SLOT -- by plan_spans, which runs AFTER this:
    #     slot 0 is a fork, anything else a merge -- and it counts the rows it
    #     covers. So a row a span covers, and a row that owns a fork, is not broken
    #     at all; a row that owns a merge breaks only BEFORE the merge cell, which
    #     then lands in the last piece with the covered rows straight after it. A
    #     break that would put the merge cell at slot 0 moves back one cell, so the
    #     arrow into it leads the piece instead: at slot 0 the same cell would draw
    #     as a FORK, and nothing downstream would object.
    #   * A MERGE CELL COSTS ITS CONVERGENCE ARROW -- see rowspan.merge_lead.
    #   * AN ARROW THAT ENDS A ROW NEVER SITS ALONE. A caller's row often ends on the
    #     step whose product starts the NEXT row, and when that arrow is the cell
    #     that overflows, "it leads the next row" gives it a row to itself with its
    #     product one row further down. That was 3 of the 4 over-wide rows without a
    #     span in the route figures measured. The break moves back one cell, so the
    #     arrow trails its substrate. One flat row keeps the single-row loop
    #     untouched.
    #
    # The snake stays single-row: on several rows nothing is broken. WRAP_ROWS
    # records what happened here, for the reply.
    rows_src = recipe["rows"]
    if wrap_width and serpentine and len(rows_src) != 1:
        WRAP_ROWS.append(("snake", len(rows_src)))
    if wrap_width and (len(rows_src) == 1 or not serpentine):
        def _cell_adv(cell, length=None):
            if "s" in cell:
                bb = S[cell["s"]][1]
                w = bb[2] - bb[0]
                if cell.get("label"):
                    w += cfg["label_dx"] + block_box([cell["label"]], cfg["label_size"])[0]
                return w
            if "plus" in cell:
                return block_box(["+"], cfg["plus_size"])[0] + 2 * cfg["plus_gap"]
            if "text" in cell:
                return _text_cell_box(cell, cfg)[0]
            # `length`: a fork branch's arrow, which fork_arrow_length forces to the
            # longest branch's -- the length it is actually drawn at
            aw = measure_arrow(cell["arrow"]) if length is None else length
            ae = A.ink_extra(cell["arrow"], aw)
            return aw + ae["box_left"] + ae["box_right"] + 2 * cfg["gap"]

        _merge_extra = _RS.merge_lead(cfg) - cfg["gap"]

        def _break(cells, merges=(), hold_last=False):
            """[(piece, width)]: `cells` broken at wrap_width by the rule above.

            `merges` are the slots of this row's merge cells, ascending: no break
            after the first of them, and a break falling ON it moves back one
            cell. `hold_last` makes the same move for an arrow that ends the row.
            Where moving back would leave nothing behind, or start a row with a
            "+", there is no break and the piece stays wide -- its width says so.
            With neither argument this is the single-row loop, cell for cell.
            """
            keep = merges[0] if merges else None
            out, cur, wx = [], [], 0.0
            for k, cell in enumerate(cells):
                adv = _cell_adv(cell) + (_merge_extra if k in merges else 0.0)
                if (cur and (wx + adv) > wrap_width and "plus" not in cell
                        and (keep is None or k <= keep)):
                    if k == keep or (hold_last and k == len(cells) - 1
                                     and "arrow" in cell):
                        if len(cur) > 1 and "plus" not in cur[-1]:
                            prev = cur.pop()
                            pw = _cell_adv(prev)
                            out.append((cur, wx - pw)); cur, wx = [prev], pw
                    else:
                        out.append((cur, wx)); cur, wx = [], 0.0
                cur.append(cell); wx += adv
            if cur:
                out.append((cur, wx))
            return out

        if len(rows_src) == 1:
            _pieces = _break(rows_src[0])
            rows_src = [_pc for _pc, _pw in _pieces]
            if len(_pieces) > 1:
                WRAP_ROWS.append(("broke", None, len(_pieces)))
            WRAP_ROWS.extend(("over", None, _j, _pw, None)
                             for _j, (_pc, _pw) in enumerate(_pieces)
                             if _pw > wrap_width)
        else:
            # plan_spans ON THE ROWS AS WRITTEN: a refusal then names the caller's
            # row numbers, and these spans are the ones the rules above are about.
            # It runs again below, on the broken rows; that pass is the one that
            # lays out.
            _grows, _sp0 = _RS.plan_spans(rows_src)
            _flen = _RS.fork_arrow_length(_grows, _sp0, measure_arrow) if _sp0 else {}
            _whole, _merges, _faces = {}, {}, {}
            for _sp in _sp0:
                if _sp["fork"]:
                    _whole.setdefault(_sp["owner"], "owns the fork from %s" % _sp["id"])
                else:
                    _merges.setdefault(_sp["owner"], []).append(_sp["slot"])
                for _r in range(_sp["owner"] + 1, _sp["last"] + 1):
                    _whole.setdefault(_r, "covered by the %s %s" % (
                        "fork from" if _sp["fork"] else "merge into", _sp["id"]))
                    _faces.setdefault(_r, []).append(_sp)
            _out = []
            for _ri, _row in enumerate(rows_src):
                if _ri in _whole:
                    # MEASURED WITH THE SPAN IT FACES. A covered row is drawn with the
                    # span cell in its slot -- a merge's convergence arrow and product
                    # after it, a fork's root in front -- and that is the width it
                    # takes on the page. Measured on a published route: a branch row
                    # of 484 pt alone fits; with the merge it runs into it is 722,
                    # and that is the whole of why that figure is over-wide.
                    _w = sum(_cell_adv(_c, _flen.get((_ri, _ci)))
                             for _ci, _c in enumerate(_grows[_ri]) if "ghost" not in _c)
                    _w += _merge_extra * len(_merges.get(_ri, ()))
                    for _sp in _faces.get(_ri, ()):
                        _root = rows_src[_sp["owner"]][_sp["slot"]]
                        if _sp["fork"]:
                            # align_x starts the branches at x0 + the root's WIDTH,
                            # not its advance: a label on the root is not counted
                            _bb = S[_root["s"]][1]
                            _w += _bb[2] - _bb[0]
                        else:
                            _w += _merge_extra + _cell_adv(_root)
                    if _w > wrap_width:
                        WRAP_ROWS.append(("whole", _ri, _w, _whole[_ri]))
                    _out.append(_row)
                    continue
                _mg = sorted(_merges.get(_ri, ()))
                _pieces = _break(_row, _mg, hold_last=True)
                if len(_pieces) > 1:
                    WRAP_ROWS.append(("broke", _ri, len(_pieces)))
                    _out.extend(_pc for _pc, _pw in _pieces)
                else:
                    _out.append(_row)
                for _j, (_pc, _pw) in enumerate(_pieces):
                    if _pw > wrap_width:
                        _held = bool(_mg) and any(_c is _row[_mg[0]] for _c in _pc)
                        WRAP_ROWS.append(("over", _ri, _j, _pw,
                                          str(_row[_mg[0]]["s"]) if _held else None))
            rows_src = _out

        # ---- serpentine: the same break, read boustrophedon ----
        # A SNAKE IS A READING ORDER, NOT A LAYOUT. Every row is still broken,
        # measured, spaced and reserved by the pipeline below; all that changes is
        # that odd rows are handed over right-to-left, with their arrows drawn
        # head-first the other way and one elbow carrying the route down.
        #
        # The turns are computed BEFORE the reversal, in REACTION order: the turn
        # leaves a row's last product and arrives at the next row's first
        # substrate, and after the reversal that substrate sits at the same edge
        # the turn left from. Reading the reversed list would aim it across the
        # whole page at the row's far end.
        if serpentine and len(rows_src) > 1:
            # THE ARROW AT THE BREAK *BECOMES* THE TURN. A straight wrap can leave
            # the step arrow lying at the end of one row (trailing) or the start of
            # the next (leading) and let the reader jump the margin; a snake has a
            # real arrow going down, so that step has to BE the turn. Adding a turn
            # ALONGSIDE it would draw the step twice -- "AcCl" on a horizontal arrow
            # and a bare elbow into the same compound. Its caption and style come
            # along: the conditions belong to the step, not to the row it ended up on.
            _turns = []
            for _ri in range(len(rows_src) - 1):
                _a = None
                if rows_src[_ri] and "arrow" in rows_src[_ri][-1]:
                    _a = rows_src[_ri].pop()["arrow"] or {}
                elif rows_src[_ri + 1] and "arrow" in rows_src[_ri + 1][0]:
                    _a = rows_src[_ri + 1].pop(0)["arrow"] or {}
                _f = next((c["s"] for c in reversed(rows_src[_ri]) if "s" in c), None)
                _t = next((c["s"] for c in rows_src[_ri + 1] if "s" in c), None)
                if _f is None or _t is None:
                    continue
                _turn = {"from": _f, "to": _t, "kind": "turn"}
                # `reverse` is deliberately NOT carried: a turn's direction is the
                # from/to pair, and a flag meaning "head on the left end" has no
                # reading on a route that goes down.
                for _k in ("above", "below", "style"):
                    if (_a or {}).get(_k):
                        _turn[_k] = _a[_k]
                _turns.append(_turn)
            for _ri in range(1, len(rows_src), 2):
                # Cells are the caller's own dicts; a reversed arrow is a NEW one.
                # `reverse` is toggled rather than set, so a retro step the caller
                # had already flipped still points the way they drew it.
                rows_src[_ri] = [
                    dict(c, arrow=dict(c["arrow"],
                                       reverse=not (c["arrow"] or {}).get("reverse")))
                    if "arrow" in c else c
                    for c in reversed(rows_src[_ri])]
            _snaked = True
            _back_rows = set(range(1, len(rows_src), 2))
            if _turns:
                recipe = dict(recipe,
                              connectors=list(recipe.get("connectors") or []) + _turns)

    # Find rowspan cells and insert ghosts so each covered row's x
    # advances past the span's slot. Same object back when nothing spans.
    rows_src, _spans = _RS.plan_spans(rows_src)

    # One length for every branch of a fork, decided BEFORE any x is
    # accumulated below -- see rowspan.fork_arrow_length for why this cannot be
    # align_x's stretch. Empty unless the figure has a fork, and then nothing
    # here changes.
    _fork_len = _RS.fork_arrow_length(rows_src, _spans, measure_arrow) if _spans else {}

    plans = []
    for _ri, row in enumerate(rows_src):
        items, x = [], 0.0
        first_left = last_right = None
        for _ci, cell in enumerate(row):
            if "s" in cell:
                bb = S[cell["s"]][1]
                w, h = bb[2] - bb[0], bb[3] - bb[1]
                label = cell.get("label")
                # compound + salt label = ONE layout unit (à la Ctrl+G): reserve the
                # label's width in the cell's advance so the next cell / "+" / arrow
                # can't overlap it. (w stays the MOLECULE width — pass 2 still uses it
                # to place both the fragment and the label; only the advance grows.)
                adv = w
                if label:
                    lw, _lh = block_box([label], cfg["label_size"])
                    adv = w + cfg["label_dx"] + lw
                items.append(dict(type="s", id=cell["s"], x=x, w=w, h=h,
                                  label=label, sublabels=cell.get("sublabels"),
                                  bracket=cell.get("bracket"), under=cell.get("under"),
                                  num=cell.get("num", True)))
                first_left = x if first_left is None else first_left
                last_right = x + adv
                x += adv
            elif "ghost" in cell:
                # a slot occupied by a rowspan cell owned by an earlier row.
                # Width is patched in by align_x once the span's x is known;
                # here it only has to exist so `post` cells land after it.
                items.append(dict(type="ghost", x=x, w=0.0, h=0.0,
                                  _span=cell["_span"]))
                last_right = x
            elif "plus" in cell:                         # "+" between products
                x += cfg["plus_gap"]
                pw = block_box(["+"], cfg["plus_size"])[0]
                items.append(dict(type="plus", x=x, w=pw))
                first_left = x if first_left is None else first_left
                last_right = x + pw
                x += pw + cfg["plus_gap"]
            elif "text" in cell:                         # a species written out
                # Typed like a structure, not like a "+": it stands in a
                # structure's slot and takes the same label furniture, so its
                # advance is measured the same way and it keeps `type="s"`'s
                # neighbours honest about where the row got to.
                tw, th = _text_cell_box(cell, cfg)
                items.append(dict(type="text", x=x, w=tw, h=th,
                                  text=str(cell["text"]),
                                  size=cell.get("size", cfg["text_cell_size"]),
                                  id=cell.get("id"), label=cell.get("label"),
                                  sublabels=cell.get("sublabels")))
                first_left = x if first_left is None else first_left
                last_right = x + tw
                x += tw
            elif "arrow" in cell:
                a = cell["arrow"]
                # A fork branch arrow is handed the longest branch's natural
                # length; every other arrow measures its own. `in` rather than
                # .get(... ) or a truthiness test: a forced length is a real
                # number that must be used even if it equalled the natural one.
                length = (_fork_len[(_ri, _ci)] if (_ri, _ci) in _fork_len
                          else measure_arrow(a))
                # A composite arrow (multi-step twins, no-go cross) puts ink
                # outside the shaft span. Reserve it here or the overhang lands on
                # the neighbouring structure. Cached so pass 2 and the row-spacing
                # pre-pass use the SAME numbers.
                ext = A.ink_extra(a, length)
                cx = x + cfg["gap"] + ext["box_left"] + length / 2
                items.append(dict(type="arrow", a=a, cx=cx, length=length, ext=ext))
                x = cx + length / 2 + ext["box_right"] + cfg["gap"]
                last_right = x
        center = ((first_left or 0.0) + (last_right if last_right is not None else x)) / 2
        plans.append(dict(items=items, right=x, center=center))

    # NOTHING TO PLACE, SAID IN WORDS. `min(lefts)` on no rows raises "min()
    # iterable argument is empty", and that sentence would go straight out to whoever
    # called the tool -- a chemist on the web app reading Python's internals about
    # a figure they asked for. Every other refusal in this pipeline names the
    # offending key and shows a working shape; this one has to as well.
    if not plans:
        raise ValueError(
            "nothing to draw: this figure has no rows. `items` (or a snaking "
            "figure's `rows`) needs at least one row holding at least one cell, "
            "e.g. [[{\"smiles\": \"CCO\", \"id\": \"1\"}]].")

    # A row's TRUE left edge, arrow ink overhang included. Hoisted out of the
    # row_align branch below because `placement` reports the same quantity once
    # the shifts are settled, and a second copy of it would drift.
    def _row_left(p):
        xs = [(it["cx"] - it["length"] / 2 - it["ext"]["box_left"])
              if it["type"] == "arrow" else it["x"]
              for it in p["items"]]
        return min(xs) if xs else 0.0

    def _row_right(p):
        """The mirror of _row_left, and NOT p["right"]. That one is the packing
        cursor: it carries a trailing `gap` after an arrow and none after a
        structure, so comparing two rows through it calls them different widths for
        ending in different kinds of cell."""
        xs = [(it["cx"] + it["length"] / 2 + it["ext"]["box_right"])
              if it["type"] == "arrow" else it["x"] + it["w"]
              for it in p["items"]]
        return max(xs) if xs else 0.0

    # ---------- serpentine: every row must FINISH where the next one STARTS -----
    # A turn leaves a row's last product and arrives at the next row's first
    # substrate, and a snake puts those two at the same edge by construction -- so
    # the turn is vertical exactly when the two rows share that edge. Rows are
    # packed independently and then left-aligned, so a short row's right edge can
    # come up 123 pt shy of the row above and the turn then comes down as a diagonal
    # elbow. Sharing the LEFT edge is free (they are left-aligned); sharing the
    # RIGHT edge means sharing a width.
    #
    # JUSTIFIED TO THE WIDEST ROW, NOT TO THE COLUMN. The widest row is the one the
    # breaker already judged full; the column is 552 pt at letter and justifying to
    # it would have put 232 pt onto the single arrow of a row whose content never
    # asked to reach that far. The stretch itself is align_x's, the one rowspan
    # already uses to make a covered row reach its span.
    #
    # THE LAST ROW IS NOT STRETCHED, for the reason the last line of a paragraph is
    # not: nothing follows it, so it has no edge to meet. It only has to START
    # where the turn above comes down -- already true when it reads forwards, and
    # a shift when it reads backwards.
    if _snaked and len(plans) > 1:
        _widths = [_row_right(p) - _row_left(p) for p in plans]
        _target = max(_widths)
        for _ri in range(len(plans) - 1):
            _short = _target - _widths[_ri]
            _arrows = [it for it in plans[_ri]["items"] if it["type"] == "arrow"]
            if _short <= 1e-9 or not _arrows:
                continue          # a row with nothing to stretch stays short, and
                                  # its turn degrades to the elbow, which is drawn
            _each = _short / len(_arrows)
            _cap = cfg.get("max_arrow")
            _off = 0.0
            for _it in plans[_ri]["items"]:
                if _it["type"] != "arrow":
                    _it["x"] += _off
                    continue
                _grow = _each if not _cap else max(0.0, min(_each,
                                                            _cap - _it["length"]))
                _it["cx"] += _off + _grow / 2.0
                _it["length"] += _grow
                _off += _grow
            plans[_ri]["right"] += _off
            plans[_ri]["center"] += _off / 2.0

    # ---- the parent fixes how far apart two cells sit --------------------------
    # THE SAME OPERATION AS THE BLOCK ABOVE, ON A SHORTER RULER. That one is handed
    # a row and a target WIDTH and pays for it out of every arrow in the row; this
    # one is handed a PAIR OF CELLS and a target DISTANCE and pays for it out of the
    # arrows between them. The body below is that body with the range narrowed.
    #
    # NOT FOLDED INTO A SHARED HELPER, and this is a decision rather than an
    # oversight: the two differ in what they measure (row edges vs cell centres), in
    # who they skip (the last row vs nothing), and in what they do when the money
    # runs out (leave the row short and let the turn degrade vs report). What is
    # left over is eight lines of arithmetic, and a helper taking four flags to
    # cover both would be longer than the copy and read worse.
    #
    # AFTER the serpentine justifier, not before, because a span is a REQUEST and
    # that is an automatic tidy-up: running first would let the tidy-up grow the row
    # further and silently break the number the caller gave. It does mean a snaked
    # figure with a span does not have flush rows -- said out loud in the docstring,
    # since the alternative is honouring neither.
    #
    # WHY A PARENT NEEDS IT, measured on a published figure: a catalytic cycle puts
    # II and III 308.7 pt apart, and the chain a SECOND engine call draws between
    # them comes to 181.2, so anchoring one end leaves the other 127.5 pt out.
    # Without a span there is no way to say so: FORTY SPACES padded into the arrow's
    # caption buy 96 of the 127 and begin printing over the ring at sixty.
    #
    # CENTRE TO CENTRE, because that is what the composer measures when it matches
    # two drawings of the same molecule (`compose.join_shape` reads the fragment
    # BoundingBox's midpoint). Left edge to left edge would differ from it by half
    # the width of two different structures, and could not be checked against
    # anything the parent can see.
    for _sv in (recipe.get("span") or []):
        _ends = [(_ri, _ci, _it) for _ri, _p in enumerate(plans)
                 for _ci, _it in enumerate(_p["items"])
                 if _it["type"] in ("s", "text")
                 and _it.get("id") is not None
                 and str(_it["id"]) in (str(_sv.get("from")), str(_sv.get("to")))]
        _want = float(_sv.get("w") or 0.0)
        if len(_ends) == 2 and _ends[0][0] != _ends[1][0]:
            # THE ONE PREFLIGHT CANNOT SEE. It reads the rows the CALLER wrote, and
            # `wrap_width` breaks one of them here, long after that -- so two cells
            # that were in one row when the spec was checked can still arrive in
            # two. Everything else in this guard is a belt on a refusal that has
            # already happened; this branch is the only reader of it.
            SPAN_SHORT.append((str(_sv.get("from")), str(_sv.get("to")),
                               _want, None))
            continue
        if len(_ends) != 2 or _want <= 0:
            continue                     # preflight refuses both; this is the belt
        (_ri, _ai, _a), (_, _bi, _b) = sorted(_ends, key=lambda t: t[2]["x"])
        _have = (_b["x"] + _b["w"] / 2.0) - (_a["x"] + _a["w"] / 2.0)
        _arrows = [_it for _it in plans[_ri]["items"][_ai + 1:_bi]
                   if _it["type"] == "arrow"]
        _short = _want - _have
        if _short <= 1e-9 or not _arrows:
            SPAN_SHORT.append((str(_sv.get("from")), str(_sv.get("to")),
                               _want, _have))
            continue
        _each = _short / len(_arrows)
        _cap = cfg.get("max_arrow")
        _off = 0.0
        for _ci, _it in enumerate(plans[_ri]["items"]):
            if _ci <= _ai:
                continue                          # in front of the span: untouched
            if _ci < _bi and _it["type"] == "arrow":
                _grow = _each if not _cap else max(0.0, min(_each,
                                                            _cap - _it["length"]))
                _it["cx"] += _off + _grow / 2.0
                _it["length"] += _grow
                _off += _grow
            elif _it["type"] == "arrow":
                _it["cx"] += _off
            else:
                _it["x"] += _off
        plans[_ri]["right"] += _off
        plans[_ri]["center"] += _off / 2.0
        if _off < _short - 0.01:                  # max_arrow took the rest
            SPAN_SHORT.append((str(_sv.get("from")), str(_sv.get("to")),
                               _want, _have + _off))

    if cfg.get("row_align", "center") == "left":
        # left-align: line every row's TRUE left edge (leading arrows included) up to a
        # common margin, so wrapped rows read top-to-bottom like text and a short final
        # row doesn't float to center.
        lefts = [_row_left(p) for p in plans]
        base = min(lefts)
        shifts = [base - l for l in lefts]
        # A BACKWARDS LAST ROW HANGS AT THE FAR EDGE. Its reaction starts at its
        # RIGHT end -- that is what reading it backwards means -- and that is the
        # end the turn above comes down onto. Left-aligning it leaves the turn
        # reaching across the whole row to find it.
        _li = len(plans) - 1
        if _snaked and _li >= 1 and _li in _back_rows:
            shifts[_li] = ((_row_right(plans[_li - 1]) + shifts[_li - 1])
                           - _row_right(plans[_li]))
    else:
        common = max(p["center"] for p in plans)
        shifts = [common - p["center"] for p in plans]

    # Per-row override — a cell may carry "under": <id of a structure in an EARLIER
    # row>, and that row is shifted so the two structures share a centre x. Without this
    # a branch arrow drawn with connectors comes out as a long diagonal, because rows are
    # packed independently and nothing lines a branch product up under its precursor.
    # Rows are visited in order, so the target row's shift is already final.
    # Auto-"under". A connector already names the structure a later row hangs
    # from, so derive the alignment instead of making every caller hand-write it.
    # Two guards, both needed:
    #   * one row is pulled at most once  -- a converge has several sources and can
    #     only sit under one of them
    #   * one source pulls at most one row -- a FORK to two rows must not align both
    #     under the same structure, or the second connector is driven straight down
    #     through the first branch
    #   * only a row's FIRST structure may anchor it -- aligning on a later structure
    #     drags the whole row sideways by everything in front of it
    # An explicit "under" on any cell in a row always beats the derived one.
    if recipe.get("connectors") and cfg.get("auto_under", True):
        _row_of, _cell_of = {}, {}
        for _ri, _p in enumerate(plans):
            for _it in _p["items"]:
                if _it["type"] == "s":
                    _row_of.setdefault(_it["id"], _ri)
                    _cell_of.setdefault(_it["id"], _it)
        _first = {}          # row index -> id of that row's FIRST structure
        for _ri, _p in enumerate(plans):
            for _it in _p["items"]:
                if _it["type"] == "s":
                    _first[_ri] = _it["id"]
                    break
        _taken, _used_src = set(), set()
        for _ri, _p in enumerate(plans):
            for _it in _p["items"]:
                if _it["type"] == "s" and _it.get("under"):
                    _taken.add(_ri)
                    _used_src.add(_it["under"])
        for _con in recipe["connectors"]:
            _f, _t = _con.get("from"), _con.get("to")
            if _f not in _row_of or _t not in _row_of:
                continue
            _tr = _row_of[_t]
            if _row_of[_f] >= _tr or _tr in _taken or _f in _used_src:
                continue
            if _first.get(_tr) != _t:
                continue          # a row hangs from its FIRST structure only
            _cell_of[_t]["under"] = _f
            _taken.add(_tr)
            _used_src.add(_f)

    if _spans:
        # a covered row's x is fixed by its span, not by row_align -- see rowspan.py
        _RS.align_x(plans, shifts, _spans, cfg, S)

    placed_cx = {}                       # structure id -> final centre x
    for ri, plan in enumerate(plans):
        for it in plan["items"]:
            if it["type"] != "s":
                continue
            if it.get("under") and it["under"] in placed_cx:
                mine = it["x"] + it["w"] / 2 + shifts[ri]
                shifts[ri] += placed_cx[it["under"]] - mine
        for it in plan["items"]:         # record only AFTER this row's shift is settled
            if it["type"] == "s":
                placed_cx.setdefault(it["id"], it["x"] + it["w"] / 2 + shifts[ri])

    # ---------- vertical extents per row -> adaptive baselines ----------
    # Each row reaches UP and DOWN from its arrow baseline by different amounts. A
    # reagent drawn ABOVE an arrow (plus the caption stacks) genuinely takes vertical
    # space, which a fixed row_gap_y ignores -> a tall reagent could collide with the
    # row above. Space rows by their real extents (never tighter than row_gap_y, so
    # single-row and short multi-row schemes keep exactly row_gap_y).
    def _cap_h(lines):
        return block_box(lines, cfg["cap_size"])[1]

    def _struct_up(it):
        u = it["h"] / 2
        if it.get("bracket"):
            u = max(u, it["h"] / 2 + 9.0)
        if it.get("label"):
            u = max(u, block_box([it["label"]], cfg["label_size"])[1] / 2)
        return u

    def _struct_down(it):
        d = it["h"] / 2 + cfg["num_gap"] + cfg["num_size"]
        for sub in (it.get("sublabels") or []):
            d += cfg["sublabel_gap"] + block_box([sub], cfg["sublabel_size"])[1]
        if it.get("bracket"):
            d = max(d, it["h"] / 2 + 9.0)
        if it.get("label"):
            d = max(d, block_box([it["label"]], cfg["label_size"])[1] / 2)
        return d

    def _arrow_up(a, length):
        # matches pass-2 stacking: above caption, then the reagent (scaled to the arrow
        # span, as pass 2 shrinks it) with its optional id label, sitting on top of
        # that caption
        above = (cfg["txt_gap_above"] + _cap_h(a["above"])) if a.get("above") else 0.0
        if a.get("reagent"):
            rb = S[a["reagent"]][1]
            rw, rh = rb[2] - rb[0], rb[3] - rb[1]
            if rw > length:
                rh *= length / rw
            lbl_h = (cfg["num_size"] + cfg["reagent_num_gap"]) if a.get("reagent_label") else 0.0
            if a.get("reagentEq"):
                rh = max(rh, _cap_h([a["reagentEq"]]))
            return above + cfg["reagent_gap"] + lbl_h + rh
        if a.get("reagentEq"):
            # equivalents with no reagent structure still occupy the reagent slot
            return above + cfg["reagent_gap"] + _cap_h([a["reagentEq"]])
        return above

    def _arrow_down(a):
        if a.get("below"):
            d = cfg["txt_gap_below"] + _cap_h(a["below"])
            if a.get("yield"):
                d += cfg["yield_gap"] + _cap_h([a["yield"]])
        elif a.get("yield"):
            d = cfg["txt_gap_below"] + _cap_h([a["yield"]])
        else:
            d = 0.0
        return d

    ups, downs = [], []
    for plan in plans:
        up = down = 0.0
        for it in plan["items"]:
            if it.get("_spanning"):
                continue          # belongs to no single row; align_y handles it
            if it["type"] == "s":
                up = max(up, _struct_up(it)); down = max(down, _struct_down(it))
            elif it["type"] == "arrow":
                e = it["ext"]                 # twins / no-go / attached arcs
                up = max(up, e["box_up"], e["up"] + _arrow_up(it["a"], it["length"]))
                down = max(down, e["box_down"], e["down"] + _arrow_down(it["a"]))
        ups.append(up); downs.append(down)
    # ---- a connector's caption is content of the gap between two rows -----------
    # Rows spaced without knowing connectors exist leave a branch arrow's
    # conditions nowhere to go, and pass 2's placer shoves them sideways (visible
    # in CAPTION_MOVED). Asking for
    # the room HERE puts the caption back where the convention wants it: beside the
    # shaft, centred on it.
    #
    # NO NEW CONSTANT. The caption is simply another thing between the two rows, so
    # it joins the chain with the SAME row_clear on each side that already separates
    # them: down + clear + caption + clear + up. Inventing a threshold here would
    # need something to measure it against, and there is nothing -- whereas
    # row_clear is already the answer to "how far apart do two pieces of row
    # content sit".
    #
    # ADJACENT ROWS ONLY. A connector from row 0 to row 3 runs vertically THROUGH
    # rows 1 and 2, and no amount of row spacing fixes that -- it is a real
    # collision and the ledger reports it. Widening a gap it does not sit in would
    # move the figure and fix nothing.
    _row_of = {it["id"]: ri for ri, plan in enumerate(plans)
               for it in plan["items"] if it["type"] == "s" and it.get("id")}
    _conn_gap = [0.0] * len(plans)
    for _con in (recipe.get("connectors") or []):
        _ri, _rj = _row_of.get(_con.get("from")), _row_of.get(_con.get("to"))
        if _ri is None or _rj is None:
            continue
        _lo, _hi = min(_ri, _rj), max(_ri, _rj)
        if _hi != _lo + 1:
            continue
        _h = max([_cap_h(_con[k]) for k in ("above", "below") if _con.get(k)]
                 or [0.0])
        if _h:
            _conn_gap[_hi] = max(_conn_gap[_hi], _h + cfg["row_clear"])

    # ---- a sibling figure parked in the gap is content of the gap too ----------
    # Same shape as the connector-caption block above, and for the same reason. A
    # composer wants to drop a catalytic cycle under one row's arrow -- published
    # figures do exactly that -- and without `reserve` the packer leaves 84 pt for
    # a cycle that is 208 tall, with no way to say so: on `draw_snaking_scheme` not
    # even a blank `sublabels` to fake it with, because its cells take no such key.
    #
    # NO NEW CONSTANT, for the block above's reason and not by analogy to it: the
    # reserved block IS one more thing between the two rows, so it joins the same
    # chain with row_clear on each side -- down + clear + reserved + clear + up.
    #
    # ONE MECHANISM, TWO ENGINES. draw_snaking_scheme renders through this
    # function (it hands over `snake_rows`), which is why both would leave the same
    # 84 and why both are handled here.
    #
    # `_conn_gap` holds two kinds of thing. The name comes from the connector-caption
    # block and stays: it is read in one place, four lines down, and renaming it
    # would touch five sites for no reader's benefit.
    for _rv in (recipe.get("reserve") or []):
        _ri = _row_of.get(_rv.get("after"))
        _h = float(_rv.get("h") or 0.0)
        if _ri is not None and _h > 0 and _ri + 1 < len(plans):
            _conn_gap[_ri + 1] = max(_conn_gap[_ri + 1], _h + cfg["row_clear"])

    row_base = [0.0] * len(plans)
    for ri in range(1, len(plans)):
        needed = downs[ri - 1] + cfg["row_clear"] + _conn_gap[ri] + ups[ri]
        row_base[ri] = row_base[ri - 1] + max(cfg["row_gap_y"], needed)

    if _spans:
        _RS.align_y(plans, _spans, ups, downs, row_base, cfg)

    # ---------- placement: the geometry a caller has to be able to read back ----
    # An OUT-PARAMETER, filled in place, exactly like `shift` -- and for the same
    # reason: the return value is a document, and these numbers are not in it.
    #
    # FRAME IS "layout": the engine's own, the one every Allocator reservation is
    # in, the one `shift` converts to PAGE. Not the caller's CONTENT frame, which
    # does not exist until finalize has run.
    #
    # FILLED HERE because every number is final by this line and no earlier one:
    # align_x has stretched whatever it stretches and align_y has moved whatever
    # it moves, while pass 2 below only consumes these. Reading them out of pass 2
    # instead would mean reading them once per emitted element.
    #
    # `arrows` IS A LIST. A linear row is A -> B -> C as often as it is one step,
    # so reporting "the row's arrow" would be silently wrong for every multi-step
    # row -- and silently is the failure mode this whole layer exists to remove.
    if placement is not None:
        _rows_out = []
        for _ri, _plan in enumerate(plans):
            _sh, _arrows = shifts[_ri], []
            for _it in _plan["items"]:
                if _it["type"] != "arrow":
                    continue
                _x0 = _it["cx"] + _sh - _it["length"] / 2
                _x1 = _it["cx"] + _sh + _it["length"] / 2
                _a = _it["a"] or {}
                # A reversed arrow is LAID OUT left-to-right and DRAWN with its
                # head on the left end (emit_arrow turns `reverse` into
                # heads="start"). Naming the ends by geometry would put "head"
                # where there is no head, so the two swap.
                if (_a.get("reverse")
                        or (_a.get("style") or {}).get("heads") == "start"):
                    _x0, _x1 = _x1, _x0
                _arrows.append({"tail": [_x0, row_base[_ri]],
                                "head": [_x1, row_base[_ri]]})
            _rows_out.append({"arrows": _arrows})
        placement.clear()
        placement.update({
            "frame": "layout",
            "rows": _rows_out,
            # The leftmost INK of the rows block. Not necessarily a tail x: an
            # arrow whose composite ink overhangs its shaft reaches further left
            # than its own tail does, and a row led by a structure has no tail at
            # all. A caller that wants a tail should read one out of `rows`.
            "left_margin": min((_row_left(_p) + shifts[_i]
                                for _i, _p in enumerate(plans)), default=0.0),
        })

    # ---------- pass 2: emit + track true bboxes (incl. text) for the finalize ----------
    parts, boxes = [], []
    # Halo / dot / locants for the cells that asked for them, keyed by structure
    # id. `behind` collects them as they are placed and is joined BEFORE `parts`,
    # because a shaded orbital is OPAQUE -- measured on ChemDraw output: a green
    # halo drawn over a benzene hid the ring outright. Empty
    # for every figure that asks for no decoration, and then nothing below changes.
    D = recipe.get("decorations") or {}
    behind: list = []
    # Cross-structure: {id: (smiles, index_map, base_id)}, what an anchor needs
    # to resolve a SMARTS against a cell. Empty unless the caller asked for one.
    SRC = recipe.get("anchor_src") or {}
    tid, aid, base, gid = 90000, 95000, 1000, 97000
    struct_pos = {}          # id -> (cx, cy, hw, hh) of each placed structure (for connectors)
    # id -> y of the BOTTOM of everything hanging under that structure (compound
    # number, salt label, sublabels). struct_pos carries a symmetric half-height,
    # which cannot express a stack that exists only below -- and a connector
    # leaving downward trimmed to the structure box would drive straight through
    # the cell's own compound number. Measured on a branch figure: 3.9 pt into
    # "caption:2", every time.
    struct_low = {}

    def put_box(x0, y0, x1, y1, owner=None, kind="structure", reserve=True):
        boxes.append((x0, y0, x1, y1))
        if alloc is not None and reserve:
            alloc.reserve(_Rect(x0, y0, x1 - x0, y1 - y0),
                          owner or "%s:%d" % (kind, len(boxes)), kind)

    def make_text(cx, base_y_target, lines, size, *, anchor="top", face=None,
                  color=None, owner=None):
        """Build a centered text-block <t> string and record its box (for the page
        extent), but DON'T append it. Returns (xml, top, bottom). Lets a caller group
        the text with a fragment before emitting (see the salt label in pass 2).

        `owner` NAMES the box this reserves, and defaults to None -- which is what
        every caller but one passes, and which leaves put_box to generate the usual
        "caption:N". It exists for the one caller that draws text
        INSIDE a box it also reserves itself (the written-out species cell), where
        two anonymous strangers on one rectangle is a collision report against
        nobody. See that branch in pass 2.

        THE NUMBERING DOES NOT MOVE when a caller names its box: put_box appends to
        `boxes` before it looks at the owner, so len(boxes) -- and therefore every
        OTHER caption's auto-name -- runs the same sequence either way."""
        nonlocal tid
        w, h = block_box(lines, size)
        asc = ASCENT_EM * size
        if anchor == "top":
            top = base_y_target
        else:  # bottom
            top = base_y_target - h
        baseline = top + asc                 # <t> p.y is the first-line baseline
        # Center justification anchors every line about cx, so multi-line captions
        # read centered (not ragged-left). The block bbox is the same [cx-w/2, cx+w/2],
        # so nothing shifts — only the internal per-line alignment changes.
        xml = C.styled_text(tid, cx, baseline, lines, size=size,
                            bold=(face == 1), align="center", color=color)
        tid += 1
        put_box(cx - w / 2, top, cx + w / 2, top + h, owner=owner, kind="caption")
        return xml, top, top + h

    def put_text(cx, base_y_target, lines, size, *, anchor="top", face=None,
                 color=None, owner=None):
        """Place a centered text block. anchor='top' -> box top at base_y_target;
        'bottom' -> box bottom at base_y_target. Returns (top, bottom)."""
        xml, top, bottom = make_text(cx, base_y_target, lines, size, anchor=anchor,
                                     face=face, color=color, owner=owner)
        parts.append(xml)
        return top, bottom

    for ri, (plan, sh) in enumerate(zip(plans, shifts)):
        base_y = row_base[ri]
        for it in plan["items"]:
            if it["type"] == "ghost":
                continue                      # a slot held by a span from above
            if it["type"] == "s":
                x = it["x"] + sh
                # a spanning cell is centred on the rows it covers, not on the
                # baseline of the row that happens to declare it
                y_top = (it["_span_cy"] if it.get("_span_cy") is not None
                         else base_y) - it["h"] / 2
                _cy = it["_span_cy"] if it.get("_span_cy") is not None else base_y
                struct_pos[it["id"]] = (x + it["w"] / 2, _cy, it["w"] / 2, it["h"] / 2)
                frag = S[it["id"]][0]
                frag_xml, _dxml = C.place_with(frag, D.get(it["id"]), base, x, y_top)
                if frame is not None and it["id"] in SRC:
                    # After place(), with the placed string: the coordinates in it
                    # are the ones actually emitted, so an anchor cannot drift from
                    # the drawing. See anchors.Anchors.record.
                    _sm, _im, _bi = SRC[it["id"]]
                    frame.record(it["id"], frag_xml, base,
                                 smiles=_sm, index_map=_im, base_id=_bi)
                base += 1000
                if _dxml:
                    behind.append(_dxml)
                put_box(x, y_top, x + it["w"], y_top + it["h"],
                        owner="cell:%s" % it["id"], kind="structure")
                # salt / annotation label to the right, vertically centered on the row.
                # compound + salt are ONE layout unit; the packer already reserved the
                # label width in pass 1 so nothing overlaps it. With cfg["group_salt"],
                # also emit them as a single ChemDraw <group> (the Ctrl+G equivalent)
                # so they select / move / align together in the app.
                if it["label"]:
                    lw, lh = block_box([it["label"]], cfg["label_size"])
                    label_cx = x + it["w"] + cfg["label_dx"] + lw / 2
                    if cfg.get("group_salt"):
                        label_xml, lt, lb = make_text(label_cx, base_y - lh / 2,
                                                      [it["label"]], cfg["label_size"],
                                                      anchor="top")
                        gbb = (x, min(y_top, lt),
                               x + it["w"] + cfg["label_dx"] + lw,
                               max(y_top + it["h"], lb))
                        parts.append(C.group(gid, [frag_xml, label_xml], bbox=gbb)); gid += 1
                    else:
                        parts.append(frag_xml)
                        put_text(label_cx, base_y - lh / 2, [it["label"]],
                                 cfg["label_size"], anchor="top")
                else:
                    parts.append(frag_xml)
                # compound number, top just under the structure.
                # "num": False suppresses it -- a COMPOSED cell (nested.py: a
                # vertical stack) numbers its members internally and has only a
                # synthetic id to print. Default True, so every existing recipe
                # takes the identical path and the goldens cannot move.
                if it.get("num", True):
                    _, nbot = put_text(x + it["w"] / 2, y_top + it["h"] + cfg["num_gap"],
                                       [str(it["id"])], cfg["num_size"], anchor="top", face=1)
                else:
                    nbot = y_top + it["h"] + cfg["num_gap"]
                # Stacked sub-labels (dr, ee, ...) beneath the number
                yy = nbot
                for sub in (it.get("sublabels") or []):
                    _, yy = put_text(x + it["w"] / 2, yy + cfg["sublabel_gap"],
                                     [sub], cfg["sublabel_size"], anchor="top")
                struct_low[it["id"]] = yy
                if it.get("bracket"):                     # native ChemDraw square brackets
                    pad = 9.0
                    parts.append(C.native_bracket(gid, x - pad, y_top - pad,
                                                  x + it["w"] + pad, y_top + it["h"] + pad))
                    gid += 2
                    put_box(x - pad - 3, y_top - pad, x + it["w"] + pad + 3,
                            y_top + it["h"] + pad,
                            owner="cell:%s.bracket" % it["id"], kind="decoration")
            elif it["type"] == "plus":                   # "+" between products
                pw, ph = block_box(["+"], cfg["plus_size"])
                put_text(it["x"] + sh + pw / 2, base_y - ph / 2, ["+"], cfg["plus_size"], anchor="top")
            elif it["type"] == "text":                   # a species written out
                # Centred on the row like a structure, THEN the same number and
                # sublabel stack underneath -- "+ 4a" is a product whose number
                # sits where every other product's number sits, and a reader
                # should not be able to tell that this one has no picture.
                x = it["x"] + sh
                tw, th = it["w"], it["h"]
                y_top = base_y - th / 2
                # ONE OBJECT, TWO PARTS, AND THEY HAVE TO SAY SO. This branch is
                # the only place in the engine that draws text inside a rectangle it
                # ALSO reserves: the glyph's own box and the cell's padded slot are
                # concentric, so they overlap by the whole glyph, every time. Left
                # anonymous the glyph comes back as "caption:N", allocator.related()
                # compares "cell:X" against "caption:N", finds nothing in common and
                # reports the cell as sitting on top of a stranger -- 7.3 units deep
                # on a bare [{"text": "X"}], with no neighbour in the figure at all.
                #
                # NOT A PLACEMENT BUG AND NOT FIXED LIKE ONE. The slot contains the
                # glyph exactly as intended (text_cell_pad each side), so nothing
                # moves: the two rectangles, the `boxes` ledger and the page extent
                # are what they were, and only the NAME changes. Same device the
                # radial engine uses for a caption on its own shaft -- extend the
                # owner's path so related() sees one object.
                _own = "cell:%s" % (it.get("id") or it["text"])
                put_text(x + tw / 2, y_top, [it["text"]], it["size"], anchor="top",
                         owner="%s.text" % _own)
                put_box(x, y_top, x + tw, y_top + th, owner=_own, kind="structure")
                nbot = y_top + th + cfg["num_gap"]
                if it.get("id"):
                    _, nbot = put_text(x + tw / 2, nbot, [str(it["id"])],
                                       cfg["num_size"], anchor="top", face=1)
                yy = nbot
                for sub in (it.get("sublabels") or []):
                    _, yy = put_text(x + tw / 2, yy + cfg["sublabel_gap"],
                                     [sub], cfg["sublabel_size"], anchor="top")
            else:
                a = it["a"]; cx = it["cx"] + sh; length = it["length"]
                ext = it["ext"]
                # One logical arrow -> one or more <arrow> elements
                els, aid = A.emit_arrow(aid, cx, base_y, length, a)
                parts.extend(els)
                _aown = "arrow:r%s.%s" % (ri, aid)
                # The same condition ink_extra uses to decide whether box_* differs
                # from *, so the two cannot disagree about whether this arrow has
                # arcs at all. No arcs -> the single put_box below, into both
                # accounts, and nothing under it runs.
                _no_arcs = not a.get("arcs")
                # ONE entry in `boxes`, the whole footprint -- so the
                # page extent does not move and neither does every later caption's
                # auto-name, which put_box derives from len(boxes).
                #
                # THE LEDGER GETS SOMETHING ELSE. `box_*` includes attached arcs, an
                # arc reaches 46 units off the shaft, and declaring that as ONE
                # rectangle tells every caption beside the arrow it is sitting on
                # the arrow. Measured on a two-arc figure: two ERRORs, 7.8
                # and 6.6 units "deep", while the captions were 4 pt clear of the shaft
                # and 0 of 181 sampled points on either curve came within 2 pt of them.
                # The connector block avoids the same mistake -- the axis-aligned
                # box of a curve is mostly the part with no ink in it.
                #
                # ink_extra hands back both pairs and says which is which in its own
                # docstring ("up / down / left / right -- the SHAFT's own ink. Captions
                # clear this, and only this"). This site reads the shaft pair.
                put_box(cx - length / 2 - ext["box_left"], base_y - 2 - ext["box_up"],
                        cx + length / 2 + ext["box_right"], base_y + 2 + ext["box_down"],
                        owner=_aown, kind="arrow", reserve=_no_arcs)
                if alloc is not None and not _no_arcs:
                    # the SHAFT's own ink...
                    _shaft = alloc.reserve(
                        _Rect(cx - length / 2 - ext["left"], base_y - 2 - ext["up"],
                              length + ext["left"] + ext["right"],
                              4 + ext["up"] + ext["down"]), _aown, "arrow")
                    # ...and each arc as the CURVE it is, sampled into chords the way
                    # layout_cycle reserves its ring arcs and feeders.
                    #
                    # THE TWO ACCOUNTS DIVERGE HERE, by up to 2 pt, in one place: the
                    # element box that went into `boxes` reaches the arrowhead's TIP and
                    # the sampled curve stops at its base (measured 1.73 and 1.99 pt on
                    # that figure's two arcs). The ledger is the short one; the page
                    # extent, which is what gets drawn, is still the head-inclusive box.
                    # That is the trade for not over-reserving by 46 -- see
                    # test_two_accounts, which pins both the size and the direction.
                    for _ai, _el in enumerate(els):
                        _apts = C.arc_points(_el, n=cfg["arc_samples"])
                        if not _apts:                # not an arc: no AngularSize
                            continue
                        # An arc MEETS the shaft -- tangentially, for an "in" arc, which
                        # is the whole point of the motif rather than a collision. Same
                        # `allow` the cycle feeders use for the ring arc they feed.
                        alloc.reserve_polyline(_apts, "%s#arc%d" % (_aown, _ai),
                                               "arrow", allow=[_shaft.rid])
                # caption gaps are measured from the OUTERMOST shaft axis (for a
                # plain arrow that is base_y, so nothing shifts), which keeps a twin
                # or no-go arrow from pushing its own text through the shaft
                top_axis = base_y - ext["up"]
                bot_axis = base_y + ext["down"]
                above_top = top_axis                      # running top of the stack above the arrow
                if a.get("above"):
                    t, _ = put_text(cx, top_axis - cfg["txt_gap_above"], a["above"],
                                    cfg["cap_size"], anchor="bottom",
                                    color=a.get("above_color"))
                    above_top = t
                below_bottom = None
                if a.get("below"):
                    _, below_bottom = put_text(cx, bot_axis + cfg["txt_gap_below"], a["below"],
                                               cfg["cap_size"], anchor="top",
                                               color=a.get("below_color"))
                if a.get("yield"):
                    ytop = (below_bottom + cfg["yield_gap"]) if below_bottom is not None \
                        else bot_axis + cfg["txt_gap_below"]
                    put_text(cx, ytop, [a["yield"]], cfg["cap_size"], anchor="top",
                             color=a.get("yield_color"))
                if a.get("reagent"):
                    rf, rb = S[a["reagent"]]
                    rw, rh = rb[2] - rb[0], rb[3] - rb[1]
                    _rdeco = D.get(a["reagent"])
                    if rw > length:                      # shrink to fit the arrow span
                        f = length / rw
                        rf = C.scale_fragment(rf, f)
                        rw, rh = rw * f, rh * f
                        if _rdeco:
                            # C.scale_fragment covers p= and BoundingBox and NOT the
                            # absolute 3D endpoints a disc carries, so a halo cannot
                            # be shrunk with its reagent. Dropped rather than drawn
                            # at full size: the halo is OPAQUE (measured on
                            # ChemDraw output), so an unshrunk glow over
                            # a shrunk reagent would hide the reagent completely.
                            # That is the one outcome worse than losing the mark.
                            DECOR_UNSCALED.append((a["reagent"], f))
                            _rdeco = None
                    # stack the reagent ABOVE the conditions text (no overlap);
                    # only draw a number if the caller gave a real reagent id (never the internal _r id)
                    lbl = a.get("reagent_label")
                    lbl_h = (cfg["num_size"] + cfg["reagent_num_gap"]) if lbl else 0.0
                    # reagent + equivalents caption = ONE group, centred on the
                    # arrow. Centring the reagent and hanging the caption off its
                    # right edge would push the caption past the arrow end.
                    _eq = a.get("reagentEq")
                    _eqw, _eqh = (block_box([str(_eq)], cfg["cap_size"])
                                  if _eq else (0.0, 0.0))
                    _gapw = cfg.get("reagent_eq_gap", 8.0) if (_eq and rw) else 0.0
                    _bandw = rw + _gapw + _eqw
                    _bandh = max(rh, _eqh)
                    ry_top = above_top - cfg["reagent_gap"] - lbl_h - _bandh
                    _bx = cx - _bandw / 2
                    _sy = ry_top + (_bandh - rh) / 2
                    _rxml, _rd = C.place_with(rf, _rdeco, base, _bx, _sy)
                    parts.append(_rxml); base += 1000
                    if _rd:
                        behind.append(_rd)
                    put_box(_bx, _sy, _bx + rw, _sy + rh,
                            owner="arrow:r%s.reagent" % ri, kind="structure")
                    if _eq:
                        _ex = _bx + rw + _gapw
                        put_text(_ex + _eqw / 2,
                                 ry_top + (_bandh - _eqh) / 2,
                                 [str(_eq)], cfg["cap_size"], anchor="top")
                    if lbl:
                        put_text(cx, ry_top + _bandh + cfg["reagent_num_gap"],
                                 [str(lbl)], cfg["num_size"], anchor="top", face=1)
                elif a.get("reagentEq"):
                    # equivalents given without a reagent structure: draw them in
                    # the reagent slot rather than dropping them on the floor
                    _eq = str(a["reagentEq"])
                    _eqh = block_box([_eq], cfg["cap_size"])[1]
                    put_text(cx, above_top - cfg["reagent_gap"] - _eqh,
                             [_eq], cfg["cap_size"], anchor="top")

    if _spans:
        # join each spanning cell to its arrows -- without this the arrows share
        # a start x but touch nothing, and the fan reads as unrelated reactions
        _bus, _bus_boxes, gid = _RS.emit_bus(_spans, plans, shifts, row_base, cfg, gid)
        parts.extend(_bus)
        for _b in _bus_boxes:
            put_box(*_b, owner="rowspan.bus", kind="decoration")

    # ---------- connectors: turning / routed arrows between placed structures ----------
    # The "layout representation" for arrows.py's elbow_hv / elbow_vh / via routes:
    # name a source + target structure id and a route kind, and the arrow is drawn
    # between their placed boxes (like the snaking engine's inter-row turns) once
    # everything is positioned. Lets a linear figure carry recycle / branch / down /
    # up-and-over arrows that a single horizontal segment can't express.
    #   recipe["connectors"] = [{"from": id, "to": id,
    #       "kind": "elbow_vh"|"elbow_hv"|"via"|"straight", "via": [[x,y],...]?,
    #       "style": {...}?, "above": [lines]?, "below": [lines]?}, ...]
    # Connector captions land after every row is placed, so nothing has moved out
    # of their way. Test the spot the convention asks for; if it is taken, try
    # the mirrored side, then further out. The FIRST candidate is the conventional
    # spot, so a caption that fits keeps its exact position.
    def _cap_free(rect):
        return alloc is None or not alloc.probe(rect, kind="caption")

    # Both searches go through alloc.free_near, which owns the probe (tol=0.0) and
    # the "it had to move" record. The CANDIDATE SEQUENCES stay here, per
    # free_near's own docstring: mirror-then-outward is not the same idea as whole
    # line heights, and merging them would be a layout change dressed as a cleanup.
    def _free_cap_x(lines, mx, top, sg, owner="?"):
        """x for a SIDE caption on a vertical-ish connector."""
        w, h = block_box(lines, cfg["cap_size"])
        base = cfg["txt_gap_side"] + w / 2
        if alloc is None:
            return mx + sg * base

        def cands():
            for trial in (sg, -sg):
                for extra in (0.0, w * 0.75, w * 1.5):
                    cx = mx + trial * (base + extra)
                    yield (cx, _Rect(cx - w / 2, top, w, h),
                           ("sideways", abs(cx - (mx + sg * base)),
                            "mirrored" if trial != sg else "outward"))
        cx, _i = alloc.free_near(cands(), "caption", owner, CAPTION_MOVED)
        # nowhere clear: keep the convention and let the overlap report itself
        return mx + sg * base if cx is None else cx

    def _free_cap_y(lines, mx, top0, sg, owner="?"):
        """top y for an above/below caption on a horizontal-ish connector."""
        w, h = block_box(lines, cfg["cap_size"])
        if alloc is None:
            return top0

        def cands():
            for extra in (0.0, h, h * 2.0, h * 3.0):
                top = top0 + sg * extra
                yield (top, _Rect(mx - w / 2, top, w, h),
                       ("away from the shaft", extra,
                        "%d of 3 steps" % round(extra / h)))
        top, _i = alloc.free_near(cands(), "caption", owner, CAPTION_MOVED)
        return top0 if top is None else top

    # ---- the SHAFT gets what the captions already had ------------------------
    # Test the route the caller asked for; if it is clear, that is the route, so
    # every figure that is already clear is byte-identical. If it is not, step
    # sideways into a clear lane and say so. Same three moves free_near makes, and the
    # same ending: nothing clear means keep the convention and let the overlap report.
    #
    # NOT free_near itself. That one probes ONE rect per candidate and hands back a
    # key; a route is a polyline whose every segment has to miss everything, and the
    # candidate is a whole re-route rather than a coordinate. Its docstring already
    # draws this line ("WHAT IS SHARED IS THE PROBE AND THE REPORT, NOT THE
    # CANDIDATES") -- what is shared here is alloc.probe and CONN_DETOURED.
    def _touches(pts):
        """[(owner, depth)] the ledger would REPORT this run against, worst first.

        Default tolerance, not the strict tol=0.0 a free-spot search normally
        wants, and that is the point: the trigger has to be exactly "the ledger
        would have complained". A stricter question would bend arrows that draw
        cleanly today and that nothing was ever going to report.
        """
        if alloc is None or _run_shapes is None:
            return []
        worst: dict = {}
        for _rect, _shape in _run_shapes(pts):
            for _res, _d in alloc.probe(_rect, kind="arrow", shape=_shape):
                worst[_res.owner] = max(worst.get(_res.owner, 0.0), _d)
        return sorted(worst.items(), key=lambda kv: -kv[1])

    def _lanes(y0, y1, direct_x):
        """Clear vertical lanes crossing the band y0..y1, nearest the straight line
        first.

        The lanes are the GAPS BETWEEN what is already reserved in that band, so
        the candidate set is read off the page rather than invented: there is no
        step size to choose and no search to bound. A gap has to be 2 * _INK wide
        to count, which leaves the shaft half its own width of air on each side --
        the one number here, and it comes from the ink, not from taste.
        """
        spans = []
        for _res in alloc.reservations():
            q = _res.rect
            if q.bottom > y0 and q.y < y1:
                spans.append([q.x, q.right])
        spans.sort()
        merged: list = []
        for a, b in spans:
            if merged and a <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        out = []
        for i in range(len(merged) - 1):
            lo, hi = merged[i][1], merged[i + 1][0]
            if hi - lo >= 2 * _INK:
                out.append((0, (lo + hi) / 2.0))
        if merged:                       # and the open page on either side
            out.append((1, merged[0][0] - _INK))
            out.append((1, merged[-1][1] + _INK))
        # INSIDE FIRST, then by how far the arrow has to lean. A lane between two
        # things that are already there costs nothing; one past the outermost of
        # them widens the figure, and width is a promise here -- page='single' is
        # a column an editor enforces. So an outside lane is a last resort rather
        # than simply a farther one, which a single distance key would make it
        # whenever the figure is lopsided.
        out.sort(key=lambda t: (t[0], abs(t[1] - direct_x)))
        return [x for _side, x in out]

    def _detour(spec, r, aid, owner):
        """`r`, or a route through a clear lane when `r` would be reported."""
        if alloc is None or _run_shapes is None or spec.get("via"):
            return r                      # the caller named this route; it is theirs
        hit = _touches(r["pts"])
        if not hit:
            return r
        (fx, fy), (_fhw, fhh) = spec["f"], spec["fh"]
        (tx, ty), (_thw, thh) = spec["t"], spec["th"]
        if abs(ty - fy) <= abs(tx - fx):
            # MOSTLY SIDEWAYS, and a lane is a COLUMN. Stepping a horizontal run
            # into a vertical corridor is not a sidestep, it is a different arrow
            # -- and the thing it would have to clear is a row, not a column. Not
            # handled: the ledger reports it as an ordinary overlap.
            return r
        m = cfg["arrow_margin"]
        # Leave through the edge the run already leaves through, so the detour is a
        # sidestep rather than a different arrow.
        sg = 1.0 if ty >= fy else -1.0
        y_out, y_in = fy + sg * (fhh + m), ty - sg * (thh + m)
        if abs(y_in - y_out) < 2 * _INK:
            return r                      # the two rows touch: no band to route in
        for cx in _lanes(min(y_out, y_in), max(y_out, y_in), (fx + tx) / 2.0):
            cand = dict(spec, kind="via", via=[[cx, y_out], [cx, y_in]])
            rc = A.route(cand, aid=aid, margin=m)
            if not _touches(rc["pts"]):
                CONN_DETOURED.append((owner, cx, [o for o, _d in hit[:3]]))
                return rc
        return r                          # nowhere clear: let the overlap report

    for con in recipe.get("connectors", []):
        fpos, tpos = struct_pos.get(con.get("from")), struct_pos.get(con.get("to"))
        if fpos is None or tpos is None:
            continue
        # LEAVE BELOW THE LABEL STACK, not below the structure box.
        #
        # `fh` is a symmetric half-height, so passing the structure's own is the same
        # as claiming nothing hangs underneath -- and a compound number always does.
        # The route then trims to the structure edge and drives through it: measured
        # 3.9 pt into "caption:2" on branch figures.
        #
        # Only when the connector actually GOES DOWN. A sideways one leaves through
        # the left or right edge, where the label stack is not, and inflating the
        # half-height there would push the arrow away from a structure for no reason.
        _fh = (fpos[2], fpos[3])
        _low = struct_low.get(con.get("from"))
        if _low is not None and tpos[1] > fpos[1]:
            _fh = (fpos[2], max(fpos[3], _low - fpos[1]))
        _th = (tpos[2], tpos[3])
        _tlow = struct_low.get(con.get("to"))
        if _tlow is not None and fpos[1] > tpos[1]:
            # ...and symmetrically: arriving from BELOW means arriving through the
            # target's own label stack.
            _th = (tpos[2], max(tpos[3], _tlow - tpos[1]))
        # "turn": PICK THE ROUTE FROM WHERE THE TWO ACTUALLY LANDED, which is not
        # known when a snake's turns are generated -- that happens at wrap time,
        # before any x has been accumulated. elbow_vh (down, then sideways into the
        # target's flank) is right when the target is off to one side; when the two
        # structures stack, its sideways leg has nowhere to go and it drove 17 pt
        # down INTO the lower compound and 16 pt out through its left edge.
        _kind = con.get("kind", "elbow_vh")
        _f_at, _t_at = (fpos[0], fpos[1]), (tpos[0], tpos[1])
        if _kind == "turn":
            _stacked = abs(fpos[0] - tpos[0]) < fpos[2] + tpos[2]
            _kind = "straight" if _stacked else "elbow_vh"
            if _stacked:
                # AND IT DROPS PLUMB. "straight" runs centre to centre, and two
                # stacked compounds' centres differ by however much the two
                # molecules differ in width -- 7 pt for the Boc carbamate over the
                # amine -- so the turn leant 7 pt across a 58 pt drop. A scheme's
                # turn between rows is drawn vertical; split the difference and go
                # straight down.
                _mx = (fpos[0] + tpos[0]) / 2.0
                _f_at, _t_at = (_mx, fpos[1]), (_mx, tpos[1])
        spec = {"kind": _kind,
                "f": _f_at, "fh": _fh,
                "t": _t_at, "th": _th,
                "style": con.get("style")}
        if con.get("via"):
            spec["via"] = con["via"]
        # ONE copy of the owner name, for put_box, the reservation, the caption
        # placer and the detour report below.
        _cown = "conn:%s->%s" % (con.get("from"), con.get("to"))
        r = A.route(spec, aid=aid, margin=cfg["arrow_margin"])
        r = _detour(spec, r, aid, _cown)
        parts.extend(r["segs"])
        aid += len(r["segs"]) + 1
        hx0, hy0, hx1, hy1 = r["head"]
        # THE SAME MODEL AS THE ALLOCATOR'S, asked for rather than restated. This
        # box is the `boxes` ledger's account of the connector -- page extent only,
        # with reserving turned off (spelled out rather than quoted, because
        # test_only_one_call_site_can_make_them_diverge COUNTS that keyword in this
        # function's source and a comment mentioning it reads as a new call site).
        # test_two_accounts holds it against the allocator's account to
        # 0.005. Padding a point hull by a hand-written 2 (half of an isotropic
        # width) is wrong at every corner of an oriented strip: measured 0.97 units
        # too far right and 0.18 too high on the via fixture. polyline_hull IS the
        # allocator's arithmetic, so there is no second copy to keep in step.
        #
        # The annotation block below asks the same way. Any further restatement
        # wants this function, not another constant.
        #
        # UNIONED with the unpadded anchor hull because the ledger counts the
        # cell anchors and the head segment, and the routed run does not
        # necessarily reach all of them.
        pts = [(fpos[0], fpos[1]), (tpos[0], tpos[1]), (hx0, hy0), (hx1, hy1)]
        pts += [tuple(p) for p in con.get("via", [])]
        _bx0, _by0 = min(p[0] for p in pts), min(p[1] for p in pts)
        _bx1, _by1 = max(p[0] for p in pts), max(p[1] for p in pts)
        _chull = _polyline_hull(r["pts"]) if _polyline_hull else None
        if _chull is not None:
            _bx0, _by0 = min(_bx0, _chull.x), min(_by0, _chull.y)
            _bx1, _by1 = max(_bx1, _chull.right), max(_by1, _chull.bottom)
        put_box(_bx0, _by0, _bx1, _by1, owner=_cown,
                kind="arrow", reserve=False)
        if alloc is not None:
            # THE RUN, not the box drawn around it. A routed connector is an L
            # (or a via chain), and the axis-aligned box of an L is mostly the
            # empty corner: reserving that box makes the connector collide with
            # whatever sits in the corner, without any ink ever going near it --
            # measured: five ERRORs on one figure, every one of them untrue.
            #
            # `pts` comes back from route() rather than being rebuilt here, for
            # the reason the two-account note at the top of this file gives: the
            # corner points are route()'s arithmetic, and a second copy of them
            # would be a second account of the same geometry.
            #
            # "straight" is unaffected by this -- its pts ARE its two
            # endpoints.
            alloc.reserve_polyline(r["pts"], _cown, "arrow")
        # Caption placement ROTATES with the arrow. Rotating a left-to-right arrow
        # 90 deg clockwise sends it downward and maps up->right, down->left, so a
        # mostly-VERTICAL arrow puts "above" on its right and "below" on its left --
        # the conditions-right / yield-left convention of a branch-point arrow.
        # Captions always above/below the shaft would, on a vertical arrow, land on
        # top of the source and target structures.
        mx, my = (hx0 + hx1) / 2, (hy0 + hy1) / 2
        if abs(hy1 - hy0) > abs(hx1 - hx0):                     # vertical-ish
            # Both side captions share ONE top (hence one baseline, since the ascent
            # depends only on the point size) so they read as horizontally aligned.
            # Centring each block on its OWN height instead puts the baselines ~1 unit
            # apart, because text_metrics' height model depends on which glyphs are
            # present ("79%" has no brackets or descenders; "n) BF3.OEt2, Me2S" has
            # both). The pair is centred on the shaft midpoint by the taller block.
            sides = [(k, sg) for k, sg in (("above", +1), ("below", -1)) if con.get(k)]
            if sides:
                boxes_ = {k: block_box(con[k], cfg["cap_size"]) for k, _ in sides}
                top = my - max(h for _, h in boxes_.values()) / 2
                for key, sg in sides:
                    put_text(_free_cap_x(con[key], mx, top, sg, _cown), top,
                             con[key], cfg["cap_size"], anchor="top")
        else:                                                   # horizontal-ish
            if con.get("above"):
                _ch = block_box(con["above"], cfg["cap_size"])[1]
                _t0 = min(hy0, hy1) - cfg["txt_gap_above"] - _ch
                put_text(mx, _free_cap_y(con["above"], mx, _t0, -1, _cown),
                         con["above"], cfg["cap_size"], anchor="top")
            if con.get("below"):
                _t0 = max(hy0, hy1) + cfg["txt_gap_below"]
                put_text(mx, _free_cap_y(con["below"], mx, _t0, +1, _cown),
                         con["below"], cfg["cap_size"], anchor="top")

    # ---------- annotations whose two ends live in two different cells ------------
    # AFTER the connectors, and after every cell has been placed, because an anchor
    # cannot be resolved until the atom it names is somewhere. Page-level objects:
    # a curve from the piperidine N to the pyridinium ring belongs to neither
    # fragment, so it cannot ride inside one (see anchors.between, which says what
    # that costs).
    if recipe.get("annotations"):
        # Deferred: anchors -> decorate -> smiles_to_frag, and smiles_to_frag is
        # what builds the fragments this module is handed. A top-level import
        # would close that loop.
        import anchors as _AN
    for _ai, _an in enumerate(recipe.get("annotations") or []):
        _axml, _aink, _aboxes = _AN.between(
            frame, _an, gid, tid, where="annotations[%d]" % _ai,
            # anchors speaks plain (x, y, w, h) tuples, as it does for the boxes it
            # hands back -- keeping allocator types out of a module that has no
            # allocator. The adapter lives here, where the Rect does.
            free=(lambda r: _cap_free(_Rect(*r))) if _Rect else None,
            moved=CAPTION_MOVED)
        if not _axml:
            continue               # unresolved; already in ANCHOR_MISSES
        parts.append(_axml)
        gid += 1
        tid += len(_aboxes)         # one text at most, but keep the ranges apart
        _aown = "annot:%s->%s" % (_an["from"]["cell"], _an["to"]["cell"])
        if _aink:
            # TWO ACCOUNTS, on purpose, exactly as the connector block does it:
            #   alloc  a chain of short boxes. One hull would be 85 x 26 pt of
            #          mostly air for a hairline and would swallow whatever sits
            #          under the bow -- measured: three compound numbers on one
            #          figure.
            #   boxes  the hull, for the PAGE EXTENT, which is all `boxes` is for.
            #          Without it the bow reaches past the page and gets clipped.
            # Padded the way the allocator pads it, by ASKING, for the reason the
            # connector-hull note above gives: a hand-written `_hw = 2.0` (half of
            # an isotropic width) is the wrong number for a strip.
            # test_two_accounts holds the two accounts to the bit.
            _bh = _polyline_hull(_aink) if _polyline_hull else None
            if _bh is not None:
                put_box(_bh.x, _bh.y, _bh.right, _bh.bottom,
                        owner=_aown, kind="decoration", reserve=False)
            if alloc is not None:
                alloc.reserve_polyline(_aink, _aown, "decoration")
        for _rect, _own, _kind in _aboxes:
            put_box(_rect[0], _rect[1], _rect[0] + _rect[2], _rect[1] + _rect[3],
                    owner=_own, kind=_kind)

    # ---------- footnotes: running prose under the whole figure -------------------
    # LAST, and anchored to the union of `boxes` rather than to the bottom row --
    # because "under everything" has to include the things that are not rows. By
    # here `boxes` already holds the connector hulls (the put_box above that skips
    # the ledger) and the annotation hulls, so a connector routed below the last row
    # pushes the footnotes down with it.
    #
    # The token is spelled out in prose nowhere here on purpose:
    # test_two_accounts.test_only_one_call_site_can_make_them_diverge counts the
    # bypass keyword in this function's source, and a comment mentioning it reads
    # to that test exactly like a new call site. layout_radial places its notes last for the same
    # reason, and says so.
    #
    # The anchor is taken ONCE, before the first line is emitted: put_box appends
    # to `boxes`, so reading the union inside the loop would let each footnote
    # push the next one down.
    _notes = recipe.get("footnotes") or []
    if isinstance(_notes, str):
        _notes = [_notes]
    if _notes:
        _fcfg = {**FOOTNOTE_DEFAULTS, **cfg}
        _fsz = _fcfg["footnote_size"]
        _fx = min(b[0] for b in boxes)
        _fcap = _fcfg["footnote_width"] or (max(b[2] for b in boxes) - _fx)
        _fy = max(b[3] for b in boxes) + _fcfg["footnote_gap"]
        # fold_lines, not a second wrapper: it is THE folding routine for this
        # pipeline, it keeps ~sub~ ^sup^ *it* #bold# spans whole, and it measures
        # with the same _line_width that block_box reserves with -- so the folded
        # line and its rectangle cannot disagree.
        for _fi, _fl in enumerate(fold_lines([str(n) for n in _notes], _fsz, _fcap)):
            _fw, _fh = block_box([_fl], _fsz)
            # LEFT, flush on the figure's own left edge, which is why this does not
            # go through make_text: that one is centred by construction. Centring a
            # footnote and back-solving for its left edge leaves a ragged margin,
            # because the engine's width and ChemDraw's differ slightly and every
            # line differs differently -- layout_grid, layout_radial and table_layout
            # each carry the same note.
            parts.append(C.styled_text(tid, _fx, _fy + ASCENT_EM * _fsz, [_fl],
                                       size=_fsz, align="left"))
            tid += 1
            # reserve=True (the default): the block is real ink, so a caption or a
            # connector landing in it must surface as a collision rather than as a
            # figure nobody measured. put_box already skips the Allocator when
            # there is none, so this is the same conditional the rest of the file
            # uses -- and the tool layer always passes one.
            put_box(_fx, _fy, _fx + _fw, _fy + _fh,
                    owner="footnote:%d" % _fi, kind="caption")
            _fy += _fh + _fcfg["footnote_lead"]

    # ---------- finalize: positive-quadrant shift + page (fit / letter / …) ----------
    # Delegate to cdxml_core.finalize so page_mode is actually honored. It never
    # clips (any page dimension smaller than the content expands to content+margin),
    # and page_mode="fit" gives the tight page.
    # Decorations FIRST, then everything else, then one stamp_z over the lot --
    # the same order and the same reason as draw_molecule's. Document order is not
    # drawing order in ChemDraw (stamp_z's docstring, measured twice), so emitting
    # early is necessary and not sufficient; the Z is what decides.
    #
    # stamp_z ONLY when there is something to layer. A figure with no decoration
    # carries no Z at all and renders correctly, and stamping those would rewrite
    # every golden to say something it never needed to say. So a figure with no
    # decoration takes a byte-identical path -- structurally, rather than by a
    # check.
    body = "\n".join(behind + parts) if behind else "\n".join(parts)
    if behind:
        body, _z = C.stamp_z(body)
    # The decoration is deliberately NOT in `boxes`: a halo overhangs its atoms by
    # 0.8 bond on every side and ChemDraw recomputes a fragment's BoundingBox from
    # the atoms anyway (measured), so there is no honest way to make the packer
    # leave room via the bbox. The overhang is reported instead -- same choice
    # squiggle made, see SQUIGGLE_OVERFLOW.
    # ...AND THE SAME AGAIN ONE ROW DOWN. A row can exist and still place nothing
    # -- `items=[]` becomes one EMPTY row, which gets past the guard above and
    # arrives here with no boxes at all. Same leak, same sentence out of Python,
    # so the same answer.
    if not boxes:
        raise ValueError(
            "nothing to draw: the rows are there but every one of them is empty. "
            "A row is a list of cells, e.g. [{\"smiles\": \"CCO\", \"id\": \"1\"}, "
            "{\"arrow\": {\"above\": [\"NaBH~4~\"]}}, {\"smiles\": \"CCO\", "
            "\"id\": \"2\"}].")
    xs = [b[0] for b in boxes] + [b[2] for b in boxes]
    ys = [b[1] for b in boxes] + [b[3] for b in boxes]
    content_bbox = (min(xs), min(ys), max(xs), max(ys))
    _account_diff(recipe, boxes, content_bbox, alloc)
    return C.finalize(header_prefix, body, content_bbox,
                      page_mode=page_mode, margin=cfg["margin"], shift=shift)
