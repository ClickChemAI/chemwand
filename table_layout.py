"""table_layout.py — the data table: a grid of strings, and nothing else.

The data-table engine, minus the features deliberately deferred:
column widths and row heights from text_metrics, three horizontal rules,
per-column alignment, markup passed through untouched, folded headers that
bottom-align, footnote marks, spanning groups with a short rule covering only
the columns they span, the {"v","sec","note"} value object, row marks, the
notes block, one-spelling-per-thing, side-by-side panels, shaded bands,
and the two style presets.

NOT here, and each waiting on something named: automatic folding to width (there
is no width budget in the default "natural" mode -- it wants "match"), a
serif face (text_metrics measures Arial and its three constants were measured on
Arial), coloured rules (no coloured Line exists in any ChemDraw file to copy),
and row bold. Composition with a reaction scheme and a structure panel lives in
layout_grid.

WHY THIS IS NOT layout_grid. A scope table sizes ONE cell for the whole figure,
because a scope table only reads as a table if the compound numbers land on the
same line across a row (layout_grid.measure). A data table is the opposite: every
column is sized by its own widest string, and the columns are deliberately
unequal. The two share their measuring and their footnotes, not their packing.

NO STRUCTURES. None of the surveyed published tables draws a structure inside a
data cell, so a cell is a string, row height comes from
text_metrics alone, and no fragment bbox is involved anywhere below. That is what
makes this file short.

TWO HEADER BANDS, NOT A TREE. A group is one level deep and that is all the
surveyed published tables have: the few that span columns use a single group over
2 or 3 adjacent columns, and none has a third level. So the header is two
bands -- group titles above, column headers below -- and every column header
lives in the LOWER band whether it is grouped or not, which is where published
tables put "entry" relative to grouped headers such as "3a"/"4a".
"""

from __future__ import annotations

import re

import beautify as BT
import cdxml_core as C
import layout_linear as LL
import text_metrics as M

# THE measuring routine -- not text_metrics directly. block_box is the one that
# understands ~sub~ ^sup^ *italic* #bold#; handing a marked-up string to the raw
# metrics counts the tildes as characters and over-measures every formula.
block_box = LL.block_box

# The <t> anchor is the first line's BASELINE, and the box top sits this far
# above it. Same constant the rest of the pipeline places text with, which is
# also what makes every cell in a row share one baseline: the offset does not
# depend on the string, so a column of digits and a column with descenders sit
# on the same line rather than on their own ink tops.
ASCENT_EM = LL.ASCENT_EM


CFG = {
    # Type size. 9 is what layout_grid already sets a yield at (yield_size), so
    # a data table and a scope table under the same header set numbers alike.
    "size": 9,
    "header_size": 9,

    # PADDING IS SET BY EYE, against the type. Three drafts were read on paper --
    # 4.5/3.15 too tight, 9.0/4.5 still crowded relative to the font -- and this
    # sits at the top of the acceptable window (col 6..10, row 4..5.75):
    #
    #     col_pad 10.0  -> 20 pt (2.2 em) between the ink of adjacent columns
    #     row_pad  5.0  -> 10 pt (1.1 em) between the ink of adjacent rows
    #
    # col_pad is exactly 2 x row_pad, which is the only relation between them
    # that survived all three drafts: a table needs more air across than down,
    # because down already has the baseline grid holding it together.
    #
    # THESE ARE POINTS, NOT RATIOS. 10.0 is 1.11 x the 9 pt type, not a whole
    # multiple of it, so if "size" ever moves these two do not follow and want
    # looking at again on paper.
    #
    # NOT anchored to BOND_LENGTH. 14.4 is the structure engine's unit of
    # "short"; a bond length either side would put 28.8 pt between two columns
    # of digits, which reads as two tables standing side by side.
    "col_pad": 10.0,
    "row_pad": 5.0,

    # The three rules. WIDTHS ONLY: <graphic GraphicType="Line"> in this
    # pipeline takes no colour attribute, and there is no coloured Line in any
    # ChemDraw file to copy a spelling from, so "rule colour" is not a knob this
    # engine can honestly offer. A few published tables want one; that is a
    # separate question, and a hand-drawn ChemDraw file answers it a different
    # way (a coloured BOND, not a graphic).
    "rule_top": 1.2,            # None -> not drawn (frame:"enclose")
    "rule_header": 0.6,
    "rule_bottom": 1.2,
    # The short rule under a group title. Same weight as the rule under the
    # header, because it is the same kind of object: a rule separating a header
    # from what it heads. Published tables do not distinguish them.
    "rule_group": 0.6,

    # --- the notes block ----------------------------------------------------
    # 8 pt against the body's 9, same step down layout_grid takes from its
    # yield_size to its footnote_size.
    "notes_size": 8,
    # ONE ROW GAP below the bottom rule. layout_grid puts 22 between a scope
    # table and its prose and calls it the loosest gap in the figure, because a
    # scope table ends in a ragged row of structures and needs the white to read
    # as a break. A data table ends in a RULE, and the rule is already the
    # break; the gap only has to keep the type off it.
    "notes_gap": 10.0,          # 2 x row_pad
    "notes_lead": 3.0,          # between wrapped lines, as layout_grid sets it

    # --- width --------------------------------------------------------------
    # "natural" -- the columns are as wide as their own widest string and the
    # table is centred in whatever the figure width turns out to be. "match" --
    # the table is stretched to `target_w`, which the FIGURE hands down (see
    # layout_grid: the parent owns the total width, a block never decides it).
    #
    # Published tables do both and the choice belongs to the caller, so this
    # is a knob and not a rule. Nothing here ever COMPRESSES: there is no
    # per-column folding yet (see _LATER["wrap"]), so a column cannot give width
    # back, and a target narrower than the natural width is ignored.
    #
    # THE CALLER SETS THIS ON THE TABLE OBJECT (`table["width_mode"]`, resolved
    # by _cfg_for), which is where the data model puts it: it is a statement
    # about this table, like `split`. Reaching it through table_style also
    # works, because every key here does.
    "width_mode": "natural",    # "natural" | "match"
    "target_w": None,           # set by the parent packer, not by the caller
    # How far "match" may pull the table past its natural width. AESTHETIC, not
    # measured. It exists so that a narrow
    # three-column table under a wide scheme does not end up as three columns of
    # digits with an inch of air between them; at 1.3 the columns stay recognisably
    # sized by their content. Prefer the result that looks right over the one that
    # hits the number, and report the miss.
    #
    # REPORTING THE MISS: layout_grid._TableBlock.fit appends a note
    # when this cap leaves the table short of the figure's W, because a caller
    # who asked for "match" and got a narrower table would otherwise have to
    # measure the output to find out. Raise it through table_style (any value
    # >= 1.0) when the figure needs the full span more than it needs the
    # content-sized columns.
    "stretch_max": 1.30,

    # --- side-by-side panels ------------------------------------------------
    # TWICE the gap between two columns of the same panel. The break between
    # panels has to be one level louder than the break between columns, or the
    # eye reads two 7-column panels as one 14-column table.
    "gutter": 40.0,             # 4 x col_pad
    # The rule down the gutter, rare in published tables. Same weight as the
    # header rule.
    "rule_vertical": 0.6,
    "vrule": False,

    # --- style --------------------------------------------------------------
    # Header weight. Surveyed published tables set the header roman somewhat
    # more often than bold, and it is one of the six parameters that vary
    # INDEPENDENTLY of the others -- which is why it is a knob rather than
    # something bundled into "the ACS look".
    "header_bold": True,

    # --- fills (drawn by cdxml_core.shaded_band) ----------------------------
    # Every one of these is [r, g, b] in 0..1, appended to the document's own
    # colortable, because the ACS eight are all saturated and a saturated band
    # is a highlighter pen rather than a tint.
    #
    # TWO OF THESE FOUR ARE MEASURED, and they are the two chosen by eye in
    # ChemDraw and read off a hand-drawn file. The other two are derived from
    # that blue rather than invented: converted to HLS it is
    # h=213°, l=0.823, s=0.995, so the group colours below hold that lightness
    # and saturation and only turn the hue. They are the numbers to argue with
    # if a band reads too strong on paper.
    "header_fill": None,                          # acs sets the measured grey
    "best_fill": [0.6463, 0.8026, 0.9992],        # MEASURED: the pale blue
    "group_fills": [[0.9992, 0.8345, 0.6463],     # derived, hue 32°
                    [0.6463, 0.9992, 0.7051]],    # derived, hue 130°
    "zebra": None,                                # or a very light grey
}


# ---------------------------------------------------------------------------
# Presets -- a named snapshot of parameter values, not a bundle
# ---------------------------------------------------------------------------
# A survey of published tables is the reason these are thin. Published tables
# do NOT fall into two style clusters; they show at least six parameters that
# take their values independently, and two apparent "clusters" in a small sample
# turn out to be two publishers each being internally consistent. So a preset is
# a convenient starting point, never a thing you are meant to stay inside.
#
# WHAT IS NOT IN HERE, and why -- every one of these is a style key that has no
# working implementation behind it yet, and putting it in a preset would mean
# the preset silently does not do what it says:
#
#   font (Arial | Times)   Blocked. text_metrics measures Arial and only Arial
#       (FONT_CANDIDATES is four Arial paths), and its SUBSCRIPT_W, PITCH_EM and
#       SUB_DROP_EM were each measured against ChemDraw round trips of Arial. A
#       Times preset that kept those numbers would lay out every column with the
#       wrong widths while every assertion still passed -- which is the exact
#       failure text_metrics.assert_font exists to prevent.
#   rule colour   Blocked. <graphic GraphicType="Line"> carries no colour in
#       this pipeline and no ChemDraw file has a coloured one to copy; a
#       hand-drawn file draws its coloured lines as BONDS instead.
#
# So today the two presets differ in two parameters, header_bold and
# header_fill (the row fills are style keys a caller sets). That is honest
# rather than impressive: the rest are named above with what they are waiting on.
STYLES = {
    # unshaded three-rule table + bold header (one publisher's look)
    "default": {"header_bold": True, "header_fill": None},
    # grey-shaded header + roman type (the other publisher's look). The grey
    # was picked by eye in ChemDraw: 0.8780 on all three channels, a 12% grey,
    # which lands in the 5-15% window the style calls for.
    "acs": {"header_bold": False, "header_fill": [0.8780, 0.8780, 0.8780]},
}

# Rule widths are deliberately NOT split between the presets. Both use
# 1.2 / 0.6 / 1.2, and line widths are one of the things published images are
# too coarse to measure (the top/bottom rules can look thicker than the rule
# under the header, but the difference would need measuring). Inventing a
# difference to make the two presets look more different would be a mistake.

_ALIGN = ("center", "left", "right")
_WIDTH_MODES = ("natural", "match")

# Everything the data model has that this engine does not honour. Named, so the
# error can say WHICH feature it is rather than "unknown key".
_LATER = {
    "wrap":   "automatic column folding",
    "mark":   "a shading option",
    "bold":   "an emphasis option",
    "split":  "a side-by-side column option",
    "style":  "a style preset",
    "width":  "a fixed column width",
}


# ---------------------------------------------------------------------------
# One spelling for each thing, chosen once
# ---------------------------------------------------------------------------
# Published tables write "no result" five ways and "not detected" three, and the
# differences carry no information. So the engine picks one of each.
#
# WHOLE CELL ONLY. Every rewrite below tests the cell's entire contents, never a
# substring, because the same characters mean something else inside a value:
# "73-78" is a range and "45/35" is a pair, and a substring rule would eat the
# dash out of the first one.
_NONE_FORMS = {"", "-", "--", "---", "–", "—", "−"}
_ONE_SPELLING = {
    "nd": "n.d.", "n.d": "n.d.", "n.d.": "n.d.",
    "ND": "n.d.", "N.D": "n.d.", "N.D.": "n.d.",
    "trace": "trace", "Trace": "trace", "TRACE": "trace",
    "traces": "trace", "Traces": "trace",
    "nr": "n.r.", "n.r": "n.r.", "n.r.": "n.r.",
    "NR": "n.r.", "N.R": "n.r.", "N.R.": "n.r.",
}
# "> 95:5" and ">99" are the same statement (published tables use both). Close
# it up.
_COMPARATOR = re.compile(r"^([<>≤≥])\s+")


def normalize(s: str) -> str:
    """The one spelling of `s`. Anything unlisted passes through."""
    t = s.strip()
    if t in _NONE_FORMS:
        return "—"
    if t in _ONE_SPELLING:
        return _ONE_SPELLING[t]
    return _COMPARATOR.sub(r"\1", t)


def _cfg(user: dict | None, style=None) -> dict:
    """CFG, then the named preset, then the caller's own overrides.

    That order is the point of a preset: it is a starting point you can argue
    with one key at a time, not a mode you are locked into.
    """
    out = dict(CFG)
    if style is not None:
        if style not in STYLES:
            raise ValueError("table['style'] = %r; the presets are %s."
                             % (style, ", ".join(sorted(STYLES))))
        out.update(STYLES[style])
    out.update(user or {})

    _check_width(out)
    return out


def _check_width(out: dict) -> None:
    """The two width knobs, checked where every path resolves through.

    Both are CALLER-REACHABLE -- `width_mode` from the table object and
    either of them from table_style -- so they are checked rather than trusted.
    """
    if out["width_mode"] not in _WIDTH_MODES:
        raise ValueError("width_mode = %r; expected one of %s: the columns are as wide as their "
                         "content, or they stretch to the figure's width."
                         % (out["width_mode"], ", ".join(_WIDTH_MODES)))
    if float(out["stretch_max"]) < 1.0:
        raise ValueError(
            "stretch_max = %r, which is below 1.0 -- that would ask the table to "
            "come out NARROWER than its own content. Nothing here can compress "
            "(there is no column folding), so 1.0 is the floor and it means 'natural "
            "width, never stretched'."
            % (out["stretch_max"],))


def _cfg_for(table: dict, user: dict | None) -> dict:
    """`_cfg`, plus the one layout key that lives on the TABLE OBJECT.

    `width_mode` is the caller's switch, so the data model puts it on the table
    itself and not in the style overrides -- it is a statement about this table,
    like `split`, not a tweak to the house style. It is still honoured from
    table_style, because every CFG key is; when both are given and they DISAGREE
    this raises instead of picking one, since either answer would be a figure the
    caller did not ask for and both spellings are equally deliberate.
    """
    out = _cfg(user, table.get("style"))
    wm = table.get("width_mode")
    if wm is not None:
        other = (user or {}).get("width_mode")
        if other is not None and other != wm:
            raise ValueError(
                "width_mode is %r on the table and %r in table_style. Those are "
                "two different figures; drop one." % (wm, other))
        out["width_mode"] = wm
        _check_width(out)       # the table's spelling gets the same check
    return out


def _reject_unknown(d: dict, allowed: set, what: str) -> None:
    """Refuse a key this engine will not honour, instead of dropping it.

    The lesson is draw_linear_scheme's silent discard of duplicate ids: a key
    that is quietly ignored produces a figure that is wrong in a way neither a
    test nor a glance at the output can catch, because the thing that is wrong
    is something that ISN'T THERE. A table that refuses to draw is cheap; a
    table that draws without the shading you asked for is expensive.
    """
    for k in d:
        if k in allowed:
            continue
        raise ValueError(
            "%s: %r is %s. It is not implemented -- and honouring it silently is not "
            "on offer, because you would get a figure that simply does not have it." % (what, k, _LATER.get(k, "not part of the data model")))


def _header_lines(text: str, note=None) -> list[str]:
    r"""A header's lines, with its footnote mark on the LAST one.

    \n is an explicit fold and the only fold there is here: an automatic
    one measures against a "table width budget", and in the default
    natural width mode there is no budget to measure against -- the width is
    whatever the columns come to.
    """
    lines = str(text).split("\n")
    if note:
        lines[-1] = "%s^%s^" % (lines[-1], note)
    return lines


def cell_text(cell, where: str = "cell") -> str:
    """One cell, as the string that gets drawn (one spelling per thing).

        "59"                                  -> 59
        {"v": "80", "sec": "77", "note": "c"} -> 80 (77)^c
        {"v": None, "sec": "91", "note": "d"} -> (91)^d
        {"v": "45/35"}                        -> 45/35

    THE MARK GOES OUTSIDE THE BRACKET, which is the majority spelling in
    published tables: `0 (63)^c` and `(91)^d` are common, `80 (77^c)` is rare.
    Either is correct -- the difference carries no information -- so this is a
    house style, not a fact about the data.
    """
    if isinstance(cell, str):
        return normalize(cell)
    if not isinstance(cell, dict):
        raise ValueError("%s is %s; a cell is a string or a {'v','sec','note'} "
                         "value object" % (where, type(cell).__name__))

    _reject_unknown(cell, {"v", "sec", "note"}, where)
    v, sec, note = cell.get("v"), cell.get("sec"), cell.get("note")
    if v is None and sec is None:
        raise ValueError("%s: a value object needs 'v' or 'sec' (or both); "
                         "with neither there is nothing to draw" % where)

    mark = ("^%s^" % note) if note else ""
    if sec is None:
        return normalize(str(v)) + mark
    bracket = "(%s)" % normalize(str(sec))
    if v is None:
        # A bracket with no main value is legal input, not a mistake: published
        # tables contain exactly this.
        return bracket + mark
    return "%s %s%s" % (normalize(str(v)), bracket, mark)


def row_cells(row: dict, where: str = "row") -> list[list[str]]:
    """A row's cells, each already a list of lines.

    The row's own note lands on the FIRST cell, because that is the entry
    number -- a row note is placed on the entry number, and most surveyed
    published tables lead with an entry column. A table without one puts the mark on whatever its
    first column is, which is the only reading of "first" available here.
    """
    out = [cell_text(c, "%s.cells[%d]" % (where, j))
           for j, c in enumerate(row["cells"])]
    if row.get("note"):
        out[0] = "%s^%s^" % (out[0], row["note"])
    # \n in a data cell is an explicit fold, same as in a header. Published
    # tables fold a catalyst cell in two, and a row is sized from the height of
    # each cell after folding, so a folded cell has to be able to reach the row
    # height.
    return [t.split("\n") for t in out]


def normalizations(table: dict) -> list:
    """[(row, col, before, after)] for every cell the one-spelling rule rewrote.

    Rewriting the author's own strings is a silent transformation, which is the
    shape of bug a pipeline like this one is prone to. So it is countable: a
    caller can print this, and "0 rewrites" and "31 rewrites" are different
    figures.
    """
    out = []
    for i, row in enumerate(table["rows"]):
        for j, c in enumerate(row["cells"]):
            if not isinstance(c, str):
                continue
            after = normalize(c)
            if after != c:
                out.append((i, j, c, after))
    return out


def group_order(table: dict) -> list:
    """Group names in order of first appearance, over the WHOLE table.

    Computed before the table is cut into panels, so a group living entirely in
    the right-hand panel still gets the colour its position in the table earns
    it. Working panel by panel would restart the alternation at each gutter.
    """
    seen = []
    for r in table["rows"]:
        mk = r.get("mark")
        if isinstance(mk, str) and mk.startswith("group:") and mk not in seen:
            seen.append(mk)
    return seen


def bands(rows: list, cfg: dict, order: list | None = None) -> list:
    """[(first_row, last_row, rgb)] -- every shaded band over the data rows.

    THREE SEPARATE THINGS, deliberately not merged (published tables confuse
    them and this engine should not):

      mark "best"        one row, the emphasis colour, e.g. the optimal entry.
      mark "group:<n>"   a RUN of adjacent rows sharing a name, one band across
                         all of them, and adjacent groups in different colours,
                         e.g. two adjacent control entries.
      zebra              every other unmarked row, a very light tint. This one
                         is decoration and carries no meaning at all, which is
                         why a marked row keeps its own colour instead.

    A group's colour follows its NAME, in order of first appearance, so the same
    group coming back later in the table comes back the same colour.
    """
    out, seen, i, n = [], list(order or []), 0, len(rows)
    while i < n:
        mk = rows[i].get("mark")
        if mk == "best":
            out.append((i, i, cfg["best_fill"]))
            i += 1
        elif isinstance(mk, str) and mk.startswith("group:"):
            j = i
            while j + 1 < n and rows[j + 1].get("mark") == mk:
                j += 1
            if mk not in seen:
                seen.append(mk)
            fills = cfg["group_fills"]
            out.append((i, j, fills[seen.index(mk) % len(fills)]))
            i = j + 1
        else:
            # Odd rows, so the first data row is clean under a shaded header.
            if cfg["zebra"] and i % 2 == 1:
                out.append((i, i, cfg["zebra"]))
            i += 1
    return [b for b in out if b[2]]


def fills_used(table: dict, cfg: dict | None = None) -> list:
    """Every distinct [r,g,b] the figure needs, in a stable order.

    Asked BEFORE anything is drawn, because a band names its colour by index
    into the document's colortable and the colortable lives in the header --
    which is written before the body. Same order of operations draw_molecule
    uses for its highlight colours.
    """
    cfg = _cfg_for(table, cfg)
    out = []
    if cfg["header_fill"]:
        out.append(tuple(cfg["header_fill"]))
    order = group_order(table)
    for p in panels(table):
        for _a, _b, rgb in bands(p["rows"], cfg, order):
            if tuple(rgb) not in out:
                out.append(tuple(rgb))
    return list(dict.fromkeys(out))


def panels(table: dict) -> list[dict]:
    """The table cut into `split` side-by-side panels.

    NOT page-break continuation. The published answer to a table too long for the
    page is to set it twice side by side with the header repeated in full -- none
    of the surveyed tables continues onto a second block -- so each panel below is
    a COMPLETE table, header and all three rules included, and each measures its
    own column widths. Published two-panel tables genuinely disagree about how
    wide a shared column is, and that is the published layout.

    ceil(n/2) so the LEFT panel takes the extra row, as in a published 9|8.
    """
    n = int(table.get("split", 1) or 1)
    if n < 1:
        raise ValueError("table['split'] = %r; expected 1 or more" % n)
    rows = table["rows"]
    if n == 1:
        return [table]
    if n > len(rows):
        raise ValueError("table['split'] = %d but there are only %d rows; one "
                         "panel would come out empty" % (n, len(rows)))
    per = -(-len(rows) // n)
    base = {k: v for k, v in table.items() if k not in ("rows", "split")}
    return [dict(base, rows=rows[i * per:(i + 1) * per]) for i in range(n)]


def _validate(table: dict) -> None:
    _reject_unknown(table, {"cols", "rows", "groups", "split", "style",
                            "width_mode"}, "table")
    cols, rows = table.get("cols"), table.get("rows")
    if not cols:
        raise ValueError("table needs 'cols'")
    if not rows:
        raise ValueError("table needs 'rows' (a table with no data rows has no "
                         "row height to measure)")

    for j, col in enumerate(cols):
        _reject_unknown(col, {"header", "align", "note", "id"}, "cols[%d]" % j)
        if not isinstance(col.get("header"), str):
            raise ValueError("cols[%d]['header'] must be a string" % j)
        al = col.get("align", "center")
        # styled_text falls back to Center on an unknown alignment, so an
        # unchecked typo here is another silent wrong answer.
        if al not in _ALIGN:
            raise ValueError("cols[%d]['align'] = %r; expected one of %s"
                             % (j, al, ", ".join(_ALIGN)))

    for i, row in enumerate(rows):
        _reject_unknown(row, {"cells", "note", "mark"}, "rows[%d]" % i)
        mk = row.get("mark")
        if mk is not None and mk != "best" and not (
                isinstance(mk, str) and mk.startswith("group:") and len(mk) > 6):
            raise ValueError(
                "rows[%d]['mark'] = %r; expected None, 'best', or 'group:<name>'. These "
                "are two different statements and are kept apart on purpose: 'best' says "
                "this row is the answer, 'group:x' says these rows belong together." % (i, mk))
        cells = row.get("cells")
        if not isinstance(cells, (list, tuple)):
            raise ValueError("rows[%d]['cells'] must be a list" % i)
        if len(cells) != len(cols):
            raise ValueError("rows[%d] has %d cells, but there are %d columns"
                             % (i, len(cells), len(cols)))
        row_cells(row, "rows[%d]" % i)          # raises on a malformed cell


def _groups(table: dict) -> list:
    """[(group, first_col, last_col)] -- validated and in column order.

    A group is addressed by the ids of the columns it spans, which means the
    columns it spans need ids; the column spec has no other key for it, so 'id'
    is the one key added to the model here rather than inventing a positional
    form the data model does not use.
    """
    cols = table["cols"]
    ids: dict = {}
    for j, c in enumerate(cols):
        cid = c.get("id")
        if cid is None:
            continue
        if cid in ids:
            raise ValueError("cols[%d] and cols[%d] share id %r"
                             % (ids[cid], j, cid))
        ids[cid] = j

    out, claimed = [], {}
    for gi, g in enumerate(table.get("groups") or []):
        _reject_unknown(g, {"header", "note", "cols"}, "groups[%d]" % gi)
        names = list(g.get("cols") or [])
        if len(names) < 2:
            raise ValueError("groups[%d] spans %d column(s); a group that spans "
                             "one column is just that column's header"
                             % (gi, len(names)))
        missing = [n for n in names if n not in ids]
        if missing:
            raise ValueError("groups[%d]['cols'] names %s, which is not any "
                             "column's 'id' (known ids: %s)"
                             % (gi, missing, sorted(ids) or "none"))
        idx = [ids[n] for n in names]
        if idx != list(range(idx[0], idx[0] + len(idx))):
            raise ValueError(
                "groups[%d] spans %s, which are columns %s -- a group has to "
                "cover adjacent columns, in order, or the short rule under it "
                "would not sit over the columns it names" % (gi, names, idx))
        for j in idx:
            if j in claimed:
                raise ValueError("column %r is in groups[%d] and groups[%d]"
                                 % (cols[j].get("id"), claimed[j], gi))
            claimed[j] = gi
        out.append((g, idx[0], idx[-1]))
    return sorted(out, key=lambda t: t[1])


def _notes_lines(notes, cfg: dict, width: float) -> list[str]:
    """The notes block, already folded to the table's width.

    One note per line before folding -- surveyed published tables set them that
    way somewhat more often than running them together as prose, and per-note is
    the more general of the two: joining them is a style decision for a preset.
    """
    out = []
    for i, n in enumerate(notes or []):
        _reject_unknown(n, {"mark", "text"}, "notes[%d]" % i)
        # mark=None is the unmarked abbreviation line: "n.d.,
        # not detected." belongs under the table without claiming a letter.
        head = ("^%s^ " % n["mark"]) if n.get("mark") else ""
        out.extend(LL.fold_lines([head + str(n.get("text", ""))],
                                 cfg["notes_size"], width))
    return out


def _measure_one(table: dict, cfg: dict) -> dict:
    """One panel's plan: column widths, both header bands, row heights, size.

    Column width is the widest of the header and EVERY cell -- not the widest of
    some typical sample. A published table can have one ligand name that sets its
    whole column, or a full sentence in a numeric column, so there is no typical
    content to measure instead.
    """
    _validate(table)
    groups = _groups(table)
    cols, rows = table["cols"], table["rows"]
    hs, bs, rp = cfg["header_size"], cfg["size"], cfg["row_pad"]

    hb = cfg["header_bold"]
    hdr = [_header_lines(c["header"], c.get("note")) for c in cols]
    grp = [_header_lines(g["header"], g.get("note")) for g, _a, _b in groups]
    body = [row_cells(r, "rows[%d]" % i) for i, r in enumerate(rows)]

    col_w = []
    for j in range(len(cols)):
        w = block_box(hdr[j], hs, hb)[0]
        for cells in body:
            w = max(w, block_box(cells[j], bs)[0])
        col_w.append(w + 2 * cfg["col_pad"])

    # A group title can be wider than the columns under it -- "conversion (%)"
    # over two three-character columns -- and a title hanging out past its own
    # short rule is simply a broken drawing. Widen the spanned columns until it
    # fits, evenly, so the sub-headers stay centred under their share of it.
    for (lines, (_g, a, b)) in zip(grp, groups):
        need = block_box(lines, hs, hb)[0] + 2 * cfg["col_pad"]
        have = sum(col_w[a:b + 1])
        if need > have:
            add = (need - have) / (b - a + 1)
            for j in range(a, b + 1):
                col_w[j] += add

    row_h = [max(block_box(c, bs)[1] for c in cells) + 2 * rp for cells in body]

    col_x, x = [], 0.0
    for w in col_w:
        col_x.append(x)
        x += w

    # Band 2 holds EVERY column's header, grouped or not. Band 1 holds the group
    # titles and exists only when there are groups.
    h2 = max(block_box(h, hs, hb)[1] for h in hdr)
    k2 = max(len(h) for h in hdr)
    h1 = max((block_box(l, hs, hb)[1] for l in grp), default=0.0)
    k1 = max((len(l) for l in grp), default=0)

    if groups:
        # The short rule sits in the MIDDLE of an ordinary row gap: row_pad of
        # air above it and row_pad below, so the distance from the group title's
        # line to the sub-header's line is 2*row_pad, exactly what separates two
        # body rows. One spacing rule, used twice.
        band1_y = rp
        rule_y = rp + h1 + rp
        band2_y = rule_y + rp
    else:
        band1_y = rule_y = None
        band2_y = rp

    head_h = band2_y + h2 + rp
    return {"col_w": col_w, "col_x": col_x, "row_h": row_h,
            "head_h": head_h, "w": sum(col_w),
            "table_h": head_h + sum(row_h),
            "hdr": hdr, "grp": grp, "groups": groups, "body": body,
            "h1": h1, "k1": k1, "h2": h2, "k2": k2,
            "band1_y": band1_y, "band2_y": band2_y, "rule_y": rule_y}


def measure(table: dict, notes=None, cfg: dict | None = None) -> dict:
    """The whole figure: every panel, the gutters between them, and the notes."""
    return _plan(table, notes, _cfg_for(table, cfg))


def _stretch(plans: list, slack: float) -> None:
    """Hand `slack` to the columns IN PROPORTION to their natural widths.

    The rule for width_mode "match", and not equal shares: a column that
    holds `4,4'-dibromo-2,2'-bipyridine` and one that holds `59` do not want the
    same number of points added to them. Proportional keeps the table looking
    like its own content after the stretch; equal shares make every column drift
    towards the same width, which is the scope-table look (one cell size for the
    whole figure) rather than the data-table one.

    In place, across ALL panels: with split:2 the two panels share one figure
    width, so the slack is one pool and not one pool each.
    """
    nat = sum(sum(p["col_w"]) for p in plans)
    if nat <= 0 or slack <= 0:
        return
    for p in plans:
        p["col_w"] = [w + slack * (w / nat) for w in p["col_w"]]
        x, col_x = 0.0, []
        for w in p["col_w"]:
            col_x.append(x)
            x += w
        p["col_x"], p["w"] = col_x, x


def _plan(table: dict, notes, cfg: dict) -> dict:
    ps = panels(table)
    plans = [_measure_one(p, cfg) for p in ps]
    order = group_order(table)
    for p, plan in zip(ps, plans):
        plan["bands"] = bands(p["rows"], cfg, order)

    gutter = cfg["gutter"] if len(plans) > 1 else 0.0

    # width_mode "match": grow into the width the parent handed down. Measured
    # BEFORE panel_x and the notes fold, so both see the stretched width -- the
    # notes belong to the whole figure and fold to its width.
    target = cfg.get("target_w")
    if cfg.get("width_mode") == "match" and target:
        nat_total = sum(p["w"] for p in plans) + gutter * (len(plans) - 1)
        # The cap, and the reason there is no floor: a target BELOW natural is
        # ignored rather than honoured, because no column can fold to give width
        # back (_LATER["wrap"]). So table.min == table.natural.
        reach = min(float(target), nat_total * cfg["stretch_max"])
        _stretch(plans, reach - nat_total)

    panel_x, x = [], 0.0
    for p in plans:
        panel_x.append(x)
        x += p["w"] + gutter
    w = x - gutter
    table_h = max(p["table_h"] for p in plans)

    # The notes fold to the WHOLE figure's width, not one panel's: they sit
    # under both panels and belong to both.
    nl = _notes_lines(notes, cfg, w)
    ns, lead = cfg["notes_size"], cfg["notes_lead"]
    notes_h = (sum(block_box([l], ns)[1] for l in nl) + lead * (len(nl) - 1)
               if nl else 0.0)

    return {"plans": plans, "panel_x": panel_x, "gutter": gutter,
            "w": w, "table_h": table_h,
            "h": table_h + (cfg["notes_gap"] + notes_h if nl else 0.0),
            "notes": nl, "notes_h": notes_h,
            # One panel is the common case, so its numbers stay reachable
            # without indexing into plans.
            **({k: plans[0][k] for k in
                ("col_w", "col_x", "row_h", "head_h", "h1", "h2")}
               if len(plans) == 1 else {})}


def _anchor(x: float, w: float, align: str, pad: float) -> float:
    """The x that <t> is anchored at, which is the justification's own edge."""
    if align == "left":
        return x + pad
    if align == "right":
        return x + w - pad
    return x + w / 2.0


def _render_one(table: dict, m: dict, cfg: dict, x0: float, y0: float,
                rule, put, band) -> float:
    """Emit one complete panel. Returns the y of its bottom rule."""
    cols = table["cols"]
    hs, bs = cfg["header_size"], cfg["size"]
    pad, rp = cfg["col_pad"], cfg["row_pad"]
    hb = cfg["header_bold"]
    pitch = hs * M.PITCH_EM
    x1 = x0 + m["w"]

    def bottom_baseline(band_y, kmax):
        """Where the LAST line of every cell in this band sits.

        Bottom alignment, not top: a published table that folds a header to
        three lines sets it flush with the bottom of the band, which is what keeps a one-line
        neighbour on the same line as the fold's last line instead of floating
        above it.
        """
        return y0 + band_y + ASCENT_EM * hs + (kmax - 1) * pitch

    # FILLS FIRST. Not because document order draws them first -- it does not,
    # ChemDraw sorts by Z -- but because stamp_z hands out Z in document order,
    # so being first in the string is what puts a band behind its own text.
    if cfg["header_fill"]:
        band(x0, y0, x1, y0 + m["head_h"], cfg["header_fill"])
    ry = y0 + m["head_h"]
    for a, b, rgb in m["bands"]:
        top = ry + sum(m["row_h"][:a])
        band(x0, top, x1, top + sum(m["row_h"][a:b + 1]), rgb)

    rule(y0, x0, x1, cfg["rule_top"])

    if m["groups"]:
        b1 = bottom_baseline(m["band1_y"], m["k1"])
        for lines, (_g, a, b) in zip(m["grp"], m["groups"]):
            gx0 = x0 + m["col_x"][a]
            gx1 = x0 + m["col_x"][b] + m["col_w"][b]
            put(lines, (gx0 + gx1) / 2.0, b1 - (len(lines) - 1) * pitch,
                hs, "center", hb)
            # ONLY the spanned columns. This is the whole point of a group rule:
            # it says which columns the title above it belongs to.
            rule(y0 + m["rule_y"], gx0, gx1, cfg["rule_group"])

    b2 = bottom_baseline(m["band2_y"], m["k2"])
    for j, col in enumerate(cols):
        al = col.get("align", "center")
        lines = m["hdr"][j]
        put(lines, _anchor(x0 + m["col_x"][j], m["col_w"][j], al, pad),
            b2 - (len(lines) - 1) * pitch, hs, al, hb)

    y = y0 + m["head_h"]
    rule(y, x0, x1, cfg["rule_header"])
    for i, cells in enumerate(m["body"]):
        # Data cells align on their FIRST line, the opposite of a header band.
        # A folded catalyst name has to start level with its own entry number;
        # a folded header has to END level with its one-line neighbours.
        base = y + rp + ASCENT_EM * bs
        for j, (col, lines) in enumerate(zip(cols, cells)):
            al = col.get("align", "center")
            put(lines, _anchor(x0 + m["col_x"][j], m["col_w"][j], al, pad),
                base, bs, al)
        y += m["row_h"][i]
    rule(y, x0, x1, cfg["rule_bottom"])
    return y


def render(table: dict, notes=None, cfg: dict | None = None, colors=None,
           x0: float = 0.0, y0: float = 0.0, tid: int = 70000,
           gid: int = 76000):
    """(body, bbox) -- the figure's objects and the box they occupy.

    Returned rather than finalized so the caller can put the table under a
    reaction scheme later without this function knowing about pages.
    The bbox is also the ONE rectangle an Allocator needs for a table: the body
    is all text of known size, so there is nothing inside it for a collision
    check to discover.

    `notes` is a sibling of the table in the object model, not part of it,
    which is why it is a separate argument rather than a key in `table`.
    """
    cfg = _cfg_for(table, cfg)
    m = _plan(table, notes, cfg)
    ps = panels(table)
    parts, fills = [], []

    def band(xa, ya, xb, yb, rgb):
        if colors is None or tuple(rgb) not in colors:
            raise ValueError(
                "no colortable index for the fill %s. A band names its colour "
                "by index, and the colortable is written into the header before "
                "the body exists -- so call document(), which asks fills_used() "
                "first and extends the header, or pass `colors` yourself."
                % (list(rgb),))
        fills.append(C.shaded_band(band.gid, xa, ya, xb, yb, colors[tuple(rgb)]))
        band.gid += 10          # four atoms, four bonds, the fragment, the area
    band.gid = gid + 1000

    def rule(y, xa, xb, w):
        # None means "somebody else is drawing this one". frame:"enclose" moves
        # the top rule up above the reaction scheme, where it encloses the
        # scheme and the table together (some published tables draw it that
        # way).
        if not w:
            return
        parts.append(C._line(rule.gid, xa, y, xb, y, w))
        rule.gid += 1
    rule.gid = gid

    def put(lines, ax, baseline, size, align, bold=False):
        parts.append(C.styled_text(put.tid, ax, baseline, lines,
                                   size=size, align=align, bold=bold))
        put.tid += 1
    put.tid = tid

    bottoms = [_render_one(p, plan, cfg, x0 + dx, y0, rule, put, band)
               for p, plan, dx in zip(ps, m["plans"], m["panel_x"])]
    y = max(bottoms)

    if len(ps) > 1 and cfg["vrule"]:
        # Down the middle of each gutter, top rule to the LOWEST bottom rule.
        # With 9|8 rows the two panels do not end together, and a rule stopping
        # at the shorter one would look like it had been cut off.
        for dx, plan in zip(m["panel_x"][1:], m["plans"][1:]):
            gx = x0 + dx - m["gutter"] / 2.0
            parts.append(C._line(rule.gid, gx, y0, gx, y, cfg["rule_vertical"]))
            rule.gid += 1

    if m["notes"]:
        ns = cfg["notes_size"]
        ny = y + cfg["notes_gap"]
        for line in m["notes"]:
            # Left, always. layout_grid found that centring footnotes and
            # back-solving for the left edge leaves a ragged margin, because
            # every note is a different width and the width carries some error.
            put([line], x0, ny + ASCENT_EM * ns, ns, "left")
            ny += block_box([line], ns)[1] + cfg["notes_lead"]
        y = ny - cfg["notes_lead"]

    # Bands ahead of everything, so stamp_z gives them the lowest Z in the
    # figure and they end up behind the rules and the type.
    return "\n".join(fills + parts), [x0, y0, x0 + m["w"], y]


def document(table: dict, notes=None, cfg: dict | None = None,
             header_prefix: str | None = None, page_mode: str = "fit") -> str:
    """A standalone CDXML document holding the table, its bands and its notes.

    Colour first, then body, then Z. The order is forced: a band names its
    colour by an index into the colortable, the colortable is part of the
    header, and the header is written before the body.
    """
    cfg = _cfg_for(table, cfg)
    header = header_prefix or C.acs_header_prefix()

    rgbs = fills_used(table, cfg)
    colors = None
    if rgbs:
        header, idx = BT.extend_colortable(header, [list(c) for c in rgbs])
        colors = dict(zip(rgbs, idx))

    body, bbox = render(table, notes, cfg, colors)
    if colors:
        # ONLY when there are bands. stamp_z's own docstring is clear that
        # every other figure in this repo carries no Z at all and renders
        # correctly, and stamping those would rewrite goldens to say something
        # they never needed to say.
        body, _ = C.stamp_z(body)
    return C.finalize(header, body, bbox, page_mode=page_mode)
