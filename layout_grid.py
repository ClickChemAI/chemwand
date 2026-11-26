"""layout_grid.py — substrate-scope grid (the n x m product table).

The second most common figure in a methodology paper and the one the desktop
pipeline did not have. Topologically it is the easy case: the products have no
edges between them, so there is no arrow routing at all, only packing.

Difference from layout_linear: the Allocator is LOAD-BEARING here, not a shadow.
There is no second `boxes` list. Every primitive reserves space first, carries
its reservation id into emit, and alloc.bbox() is what finalize() is handed. If
something is drawn without being reserved, the page extent is wrong and you find
out immediately instead of three bug reports later.

Recipe (fragments already built, same convention as layout_linear):

    {
      "structures": { id: (fragment_str, bbox), ... },     # bbox = [x0,y0,x1,y1]
      "products": [
        {"s": id, "yield": "82%"?, "sublabels": ["94% *ee*"]?,
         "foot": "b"?, "label": "2HCl"?},
        ...
      ],
      "cols": 4?,                      # omitted -> derived from the page width
      "footnotes": ["*a* Conditions: ..."]?,
    }

Text may use the usual markers: ~sub~ ^sup^ *italic* #bold#.
"""

from __future__ import annotations

import math
import re

import arrows as A
import beautify as BT
import cdxml_core as C
import layout_linear as LL
import legend as _LG
import table_layout as TL
from allocator import Allocator, Policy, Rect

block_box = LL.block_box
ASCENT_EM = getattr(LL, "ASCENT_EM", 0.72)

# Base id for a cell's composed legend. Its own high range because compose_legend
# hands back a <group> whose children number upward from it, and C.place adds a
# per-cell base on top -- overlapping the grid's own ids would renumber a legend
# onto a structure.
_LEGEND_GID = 60000

# Base id for a cell's Ctrl+G <group> (one per cell, counting up). ITS OWN RANGE,
# above every other pool this engine hands out: legend groups start at 60000 and
# number upward, the grid's text ids at 90000, the section bands at 96100, panels
# at 97500. A group id colliding with a child's id is the one collision that does
# not show up as a wrong drawing -- ChemDraw resolves it silently and the group
# ends up holding the wrong object.
_CELL_GID = 98000

# A section letter's own text id. ITS OWN POOL because the band already spends
# `tid` on the title and its caller advances tid by exactly one -- a second text
# object taken from the same counter would land on the next section's title.
# Indexed by section, so it cannot collide with itself either.
_SECTION_LETTER_TID = 95000

# The pill badge's graphic id. Its own pool too, for the same reason: not to
# shift anyone else's numbers. The section bands' gids start at 96100 and take
# two per section (lead line + right line), so squeezing the badge into that
# sequence would push every later section's lines up by one -- the drawing is
# identical, but every golden file changes. 96900 is eight hundred ids above
# 96100, and no figure has eight hundred sections.
_SECTION_PILL_GID = 96900

# [(cell id, width)] when a cell's value table is what set the grid's column
# width. Filled by measure(), read by the page-overflow note. A list because
# the reports elsewhere in this pipeline are lists and reset the same way.
_WIDEST_LEGEND: list = []

# [(cell id, table width, structure width, cause)] for every value table drawn
# more than legend_fit_hi times its own structure's width. Filled by measure(),
# read by _finish UNCONDITIONALLY -- unlike _WIDEST_LEGEND above, which only
# speaks when the figure also overflowed the page WIDTH. That is the wrong
# trigger for this: a wide value table makes the figure narrow and TALL, because
# it inflates the cell size until the column count collapses. Measured -- six
# products plus one Markush cell came out 322 x 1189 pt, 469 pt past the page,
# and not one line named the table that did it.
_LEGEND_WIDE: list = []

# Only the keys the grid adds. Everything else (num_size, sublabel_gap, ...) is
# inherited from layout_linear.CFG so the two engines set type identically.
GRID_DEFAULTS = {
    "cell_pad_x": 20.0,      # breathing room either side of the widest thing in a cell
    "cell_pad_y": 14.0,
    "yield_size": 9,
    "footnote_size": 8,
    "footnote_gap": 22.0,
    "footnote_lead": 3.0,
    "fill_width": True,
    "max_spread": 2.0,
    "last_row": "auto",      # "auto" | "center" | "left"
    "label_align": "center", # "center" | "left"
    "struct_valign": "center",   # "center" | "bottom"
    "footnote_width": None,  # None -> wrap to the grid width

    # THE FIGURE'S OWN TITLE -- "Table 2: Substrate scope of alkenes." Every
    # surveyed figure has one and this engine had nowhere to put it.
    #
    # SAME PROSE BAND AS THE FOOTNOTES, TURNED UPSIDE DOWN, and that is where both
    # numbers come from rather than from taste: one block of running text folded to
    # the figure's own width, set flush left, one gap away from the drawing. So the
    # gap above equals the gap below (footnote_gap) and the line lead is the same.
    # Only the SIZE differs, and it differs for a reason -- a heading is read before
    # the figure and a footnote after it, so the heading is set at body size while
    # the footnote stays smaller.
    #
    # BOLD IS INLINE, NOT A FLAG. Published figures set a bold label and regular
    # prose -- "**Table 2:** Substrate scope of alkenes" -- so the caller writes
    # `#Table 2:# Substrate scope of alkenes` and the existing marker does it.
    # Measured: block_box reads that inline marker to within 0.01 pt of real Arial
    # (163.97 vs 163.96), so fold_lines breaks it in the right place. A whole-line
    # bold flag would NOT be safe here: fold_lines takes no `bold`, so a bold title
    # folded with the regular table would break too late -- the same defect the
    # section band guards against (see _emit_section_band).
    "title_size": 10,
    "title_gap": None,        # None -> footnote_gap, i.e. symmetric with the footnotes
    "title_lead": None,       # None -> footnote_lead
    "header_gap": 26.0,
    "cols": None,
    "span": None,            # None -> the page decides how wide the grid may be

    # A header species written rather than drawn. Its own two keys because
    # the grid's type scale is its own -- layout_linear's text_cell_size is set
    # against a scheme's structures, and a header sits above a product grid.
    "hdr_text_size": 11,
    "hdr_text_pad": 6.0,

    # The break between a cell's label stack and the value table under it.
    "legend_gap": 6.0,
    # HOW MUCH WIDER THAN ITS OWN STRUCTURE a value table may be drawn before the
    # report says so. The number is read off published figures: the table and
    # the Markush structure are drawn about the same width so the pair reads as
    # one object -- and in those figures the generic structure is drawn BIGGER
    # than an ordinary product for exactly this reason.
    #
    # A SOFT TARGET: going past it is acceptable and the figure is still right.
    # Nothing enforces it; it only decides when a note is worth writing.
    #
    # ONLY AN UPPER BOUND. A table NARROWER than its structure changes nothing
    # else in the figure, so a note about it would be one nobody acts on.
    "legend_fit_hi": 1.2,

    # CTRL+G: emit each cell's structure, compound number, yield and value table
    # as ONE ChemDraw <group>, so they select, move and align together in the app
    # instead of being loose objects that have to be re-selected every time.
    #
    # ON by default, which rests on two checks:
    #
    #   THE SPELLING. See cdxml_core.group()'s docstring for both directions; the
    #   short version is that this repo's own <group> survives a ChemDraw round
    #   trip AS A GROUP with its children still nested.
    #
    #   THE READERS. Wrapping a cell puts its objects one level deeper, so a
    #   reader that asks for them by direct-child path sees almost nothing: the
    #   same two-cell figure rendered to 45 drawing elements ungrouped and 4
    #   grouped -- a near-blank preview, silently. A reader has to descend through
    #   groups (or use .iter(), which is recursive); then the two renders come out
    #   byte-identical.
    #
    # WHAT IT COSTS: one <group> wrapper per cell and nothing else. Measured on a
    # four-cell table -- 156 objects before, 156 after, every tag, every attribute
    # and every string identical. RDKit and the InChI gate read it the same either
    # way (13/14 atoms back, checked 2 / matched 2).
    "group_cells": True,

    # --- section bands -------------------------------------------------------
    # PROVISIONAL. The shape is taken from the source figure -- a short lead line,
    # the title, then a line out to the right edge, NOT one rule struck through
    # the text -- and the numbers are to be tuned against a real one on paper. So
    # these are not pinned to measured windows the way fold_bond or the template
    # guards are; each one is only anchored to the nearest thing already in the
    # figure, named here so the next reader knows what to argue with:
    "section_size": None,    # None -> cfg["num_size"], i.e. compound-number size
    "section_lead": 14.0,    # one bond length: the drawing's own unit of "short"
    "section_text_gap": 5.0, # the break either side of the title
    "section_lw": 0.8,       # same weight as the rule under a header (see below)
    "section_gap_above": 20.0,   # < header_gap 26: a section is a smaller break
    "section_gap_below": 10.0,   # == num_gap, so title-to-row reads like num-to-row
    # HOW THE BAND IS DRAWN. Published figures use at least six shapes, and with
    # only one of them every sectioned figure would come out in the same dress
    # whatever the source did.
    #   "rule"      lead line, title, line to the right edge   (this engine's own)
    #   "plain"     the title and nothing else                 ("Condition A" / "B")
    #   "centered"  line, title centred, line                  (the most common
    #               shape in published figures)
    # A section may override it; this is the figure-wide default. NOT DONE HERE:
    # the filled pill badge and the tinted band both need a fill behind text,
    # which is a different primitive -- see the tool's docstring.
    "section_style": "rule",
    # Pill badge: the padding round the title and the corner radius. The radius
    # is in C.rect's units (1/100 pt); 600 is the value ChemDraw itself writes
    # in a file it saved.
    "pill_pad_x": 7.0,
    "pill_pad_y": 2.5,
    "pill_radius": 600,

    # A BOX DRAWN AROUND SOMETHING THE FIGURE ALREADY CONTAINS. Published figures
    # use a frame as their universal "this is not a product cell" mark -- the
    # catalyst or ligand beside the header, a whole section, the unsuccessful
    # substrates, the whole figure -- and 30 of 102 surveyed figures carry one.
    # `graphics` cannot draw it: it takes absolute coordinates, and nobody knows
    # where the grid landed until after it is packed.
    "box_pad": 6.0,           # breathing room between the box and what it encloses
    "box_gid": 97000,         # its own id pool; panels take 97500, cells 98000

    # --- the data table (combination layer) -----------------------------------
    # "above" puts the STRUCTURE PANEL above the table, which is 16/18 of the
    # surveyed figures; "below" is the other two. The name refers to the panel,
    # not to the table.
    "panel_pos": "above",
    # "none" | "enclose". enclose raises the table's top rule above the reaction
    # scheme so the scheme sits inside the table's frame -- "top rule drawn above
    # the reaction scheme", 5/40 surveyed figures; a full box round both is a
    # different thing and not this.
    "frame": "none",
    # Between the panel and the table. Same break as header-to-grid, because it
    # is the same kind of break: one whole block ending and another starting.
    "table_gap": 26.0,           # == header_gap

    # Clearance between a within-grid connector and the two cells it links.
    # 8.0, not layout_linear's arrow_margin of 16.0, and the number is chosen
    # against the space that exists rather than by taste: the gutter between two
    # structures is 2 * cell_pad_x = 40 pt, and 16 either side would leave 8 pt of
    # shaft. Same value snaking uses for the same reason (its rows are close too).
    "conn_margin": 8.0,
}


def _wrap(paragraph: str, max_w: float, size: float):
    """Greedy word wrap. Footnotes under a scope table are running prose spanning
    the full column, not one line per note, so they have to be broken to width.

    THE ONE FOLDING ROUTINE. This delegates to `fold_lines` -- which calls itself
    "THE ONE FOLDING ROUTINE for this pipeline" -- rather than running a second
    greedy loop that cuts on a bare `paragraph.split()`. The two agree on the
    break width to the last
    decimal (`block_box(s, n)[0]` IS `_line_width(s, n)`), so plain prose folds
    identically; they part company on ONE thing, and it is the thing fold_lines'
    tokeniser exists for.

    A MARKER SPAN WITH A SPACE IN IT GETS CUT IN HALF by `.split()`, which sees
    "*combined mass of both diastereomers*" as five words and will happily end a
    line after "*combined", which strands one asterisk on each line: ChemDraw then
    prints both marks and the italic is LOST -- measured, face="2" became no face
    at all. Spans without spaces (`~2~`, `^1^`, `*i*-PrOH`) tokenise the same
    either way, which is why it is easy to miss: it only bites when the fold
    happens to land inside a multi-word span, so it comes and goes with the
    column width.

    Italic is not decoration here. *trans*, *cis*, *i*-Pr, *N*-, *o*-/*m*-/*p*-,
    *R*/*S* all mean something, and a footnote that prints them upright between two
    asterisks says something different from what the author wrote.
    """
    return LL.fold_lines([paragraph], size, max_w)


def _cfg(user: dict | None) -> dict:
    out = dict(LL.CFG)
    out.update(GRID_DEFAULTS)
    out.update(user or {})
    return out


def _stack(p: dict, cfg: dict):
    """What sits under a structure, top to bottom.

    Each entry: (kind, lines, size, gap_before, bold).
    """
    out = []
    ident = str(p.get("s", ""))
    if p.get("foot"):
        ident = "%s^%s^" % (ident, p["foot"])
    # "num": False -> the cell has no compound number of its own. `s` is the
    # STRUCTURES-DICT KEY, which is a caller's id when they gave one and an
    # invented "s1"/"s2" when they did not -- so printing it unconditionally puts
    # the engine's internal bookkeeping on the page as a compound number: a
    # Markush cell whose numbers live in its value table would come out labelled
    # "s1". layout_linear reads the same key.
    if p.get("num", True):
        out.append(("num", [ident], cfg["num_size"], cfg["num_gap"], True))
    if p.get("label"):
        out.append(("label", [str(p["label"])], cfg["label_size"],
                    cfg["sublabel_gap"], False))
    if p.get("yield"):
        out.append(("yield", [str(p["yield"])], cfg["yield_size"],
                    cfg["yield_gap"], False))
    for sl in (p.get("sublabels") or []):
        out.append(("sub", [sl], cfg["sublabel_size"], cfg["sublabel_gap"], False))
    return out


def _legend_widest_column(spec: dict):
    """(column number, its width, its longest entry) for the column that set this
    value table's width -- or None when the table has no columns.

    NAMES THE CAUSE, not the symptom. "40% wider than its structure" is a fact
    the caller cannot act on; "column 2 is 92 of the 146 pt, and its longest
    entry is 'R = 4-morpholinoC6H4,'" is the same fact with the edit in it. The
    page-overflow note one function down follows the same rule.

    Called ONLY when the note is about to fire. _LG.measure re-reads the recipe,
    and a second measurement taken for layout is the trap measure() warns about
    three lines up -- this one is read for prose and nothing is sized from it.
    """
    m = _LG.measure(spec)
    colw = m.get("colw") or []
    if not colw:
        return None
    j = max(range(len(colw)), key=lambda k: colw[k])
    lcfg = m.get("cfg") or {}
    size = lcfg.get("size", 9)
    text, wide = "", 0.0
    for row in (spec.get("rows") or []):
        if isinstance(row, str) or j >= len(row):
            continue                      # a bare string is a note row: no columns
        w = _LG._cell_width(row[j], size, lcfg)
        if w > wide:
            text, wide = _LG._cell_text(row[j]), w
    return j + 1, colw[j], text


def measure(recipe: dict, cfg: dict):
    """Pass 1: size every product, then fix ONE cell size for the whole grid.

    A scope table only reads as a table if the cells are congruent -- the
    compound numbers must land on the same line across a row. So the structure
    band and the label band are each sized to their own maximum, rather than
    each cell shrink-wrapping itself.
    """
    S = recipe["structures"]
    meas = []
    for p in recipe["products"]:
        frag, bb = S[p["s"]]
        stack = _stack(p, cfg)
        widths = [block_box(l, sz)[0] for _k, l, sz, _g, _b in stack]
        # The value table under a Markush cell. Composed HERE, once, and
        # carried in the record -- not rebuilt in _emit_cell, for the reason the
        # header's text box is not: composing twice is two chances for the width
        # reserved and the width drawn to disagree.
        #
        # It is NOT another _stack entry. A stack entry is one line of type with
        # a size; a legend is a block with its own column alignment and its own
        # box, and it is routinely wider than the structure above it.
        leg_xml = leg_box = None
        if p.get("legend"):
            leg_xml, leg_box = _LG.compose_legend(p["legend"], gid=_LEGEND_GID)
            # ...AND HOW IT CAME OUT NEXT TO ITS OWN STRUCTURE. Report only; the
            # bytes are unchanged. Folding cannot enforce the target and was
            # measured before it was written off: 1 column is the narrowest a
            # table gets, 2 columns is 2.39x wider (two blocks plus the block
            # gap), so a table already past the target has no fold that helps.
            # Short entries land at 0.83-0.91 of the structure on their own,
            # which is inside the target with nothing done.
            _sw = bb[2] - bb[0]
            _lw = leg_box[2] - leg_box[0]
            if _sw > 0 and _lw > _sw * cfg["legend_fit_hi"]:
                _LEGEND_WIDE.append((str(p.get("s", "?")), _lw, _sw,
                                     _legend_widest_column(p["legend"])))
        meas.append({
            "p": p,
            "frag": frag,
            # The halo / dot / locants for this cell, looked up HERE for the
            # same reason the legend is composed here -- the record is the one
            # place a cell's pieces are collected, so _emit_cell cannot get a
            # different answer than the measurement did.
            "deco": (recipe.get("decorations") or {}).get(p["s"]),
            "sw": bb[2] - bb[0],
            "sh": bb[3] - bb[1],
            "stack": stack,
            "stack_w": max(widths) if widths else 0.0,
            "stack_h": sum(g + block_box(l, sz)[1]
                           for _k, l, sz, g, _b in stack),
            "legend": leg_xml,
            "legend_w": (leg_box[2] - leg_box[0]) if leg_box else 0.0,
            "legend_h": ((leg_box[3] - leg_box[1]) + cfg["legend_gap"]
                         if leg_box else 0.0),
        })

    cell_w = (max(max(m["sw"], m["stack_w"], m["legend_w"]) for m in meas)
              + 2 * cfg["cell_pad_x"])

    # WHICH CELL MADE THE TABLE THIS WIDE, when the answer is a legend.
    # One cell size for the whole grid means a wide value table widens EVERY
    # column -- which is what the published figures do, so it is right -- but it
    # makes "the table suddenly got wider" an effect with a cause the caller
    # cannot see. Told only "17 pt too wide", they reach for the page size; told
    # which legend did it, they can split it. A report has to point at the
    # cause, not just the symptom.
    _WIDEST_LEGEND.clear()
    widest = max(meas, key=lambda m: m["legend_w"])
    if widest["legend_w"] and widest["legend_w"] >= max(
            max(m["sw"], m["stack_w"]) for m in meas):
        _WIDEST_LEGEND.append((str(widest["p"].get("s", "?")),
                               widest["legend_w"]))
    return meas, cell_w


def _flatten(recipe: dict):
    """Accept either a flat product list or explicit rows.

    Real scope tables are often ragged -- 3,3,3,4,4 when the first rows draw the
    aryl group in full and the last rows abbreviate it to Ph. No single column
    count expresses that, so recipe["rows"] is taken literally and every
    heuristic below is skipped.
    """
    if recipe.get("rows"):
        rows = [list(r) for r in recipe["rows"]]
        return [p for r in rows for p in r], [len(r) for r in rows]
    return list(recipe["products"]), None



def measure_header(recipe: dict, cfg: dict):
    """Pass 1 for the SCHEME block: what each item is, how wide, and how far the
    row reaches above and below its axis. Emits nothing.

    Lifted out of _emit_header without changing an expression, because the
    figure's total width cannot be resolved until every block's width is known
    and the scheme was the only block whose width was unknowable until it had
    already been drawn. The grid has measure() and the table has
    TL.measure(); this is the third.

    Returns None when there is no header -- including a header whose "items" is
    empty, which is not a block. (Otherwise the max() below would see an empty
    sequence and raise; the tool layer only ever sets recipe["header"] when at
    least one item survived, so an ordinary tool call never takes that path.)
    """
    h = recipe.get("header")
    if not h or not h.get("items"):
        return None
    S = recipe["structures"]

    # width of every item, and how far the tallest thing reaches up/down.
    # `boxes` carries the measurement of a text item forward to _emit_header, by
    # position in `sized`. The other kinds re-measure themselves down there and
    # get away with it -- "+" is a constant and a structure's box comes from the
    # fragment -- but a caller's string does not, and two spellings of one
    # measurement is a cell drawn wider than the slot reserved for it, with
    # nothing in the picture to say so.
    sized, boxes = [], {}
    for it in h["items"]:
        if it.get("text") is not None:
            size = it.get("size", cfg["hdr_text_size"])
            tw, th = block_box([str(it["text"])], size)
            tw += 2 * cfg["hdr_text_pad"]
            # A number under it only when one was asked for. "#2a#" is a bold
            # compound number standing alone -- it IS the cell -- while {"text":
            # "+", "id": "4a"} is a connective with a number beneath, and the two
            # must not render alike or nobody can tell which they wrote.
            lbl = (cfg["num_gap"] + block_box([str(it["id"])], cfg["num_size"])[1]
                   if it.get("id") else 0.0)
            boxes[len(sized)] = (tw, th, size)
            sized.append(("text", it, tw, th / 2.0, th / 2.0 + lbl))
        elif it.get("plus"):
            w, hh = block_box(["+"], cfg["plus_size"])
            sized.append(("plus", it, w + 2 * cfg["plus_gap"], hh, 0.0))
        # MEMBERSHIP, not truthiness -- the same rule nested.py follows: an arrow
        # with no captions is still an arrow. `it.get("arrow")` on {"arrow": {}}
        # is falsy, so a truthiness test would send a bare arrow to the structure
        # branch to die on KeyError 's'.
        elif "arrow" in it:
            a = it["arrow"] or {}
            caps = [a[k] for k in ("above", "below") if a.get(k)]
            cw = max([block_box(c, cfg["cap_size"])[0] for c in caps] or [0.0])
            w = max(cfg["min_arrow"], cw + cfg["arrow_pad"] * 2)
            up = (cfg["txt_gap_above"] + block_box(a["above"], cfg["cap_size"])[1]
                  if a.get("above") else 0.0)
            dn = (cfg["txt_gap_below"] + block_box(a["below"], cfg["cap_size"])[1]
                  if a.get("below") else 0.0)
            sized.append(("arrow", it, w + 2 * cfg["gap"], up, dn))
        else:
            bb = S[it["s"]][1]
            sw, sh = bb[2] - bb[0], bb[3] - bb[1]
            lbl = cfg["num_gap"] + block_box([str(it["s"])], cfg["num_size"])[1]
            sized.append(("s", it, sw + cfg["gap"], sh / 2.0, sh / 2.0 + lbl))

    return {"sized": sized, "boxes": boxes,
            "natural": sum(w for _k, _i, w, _u, _d in sized),
            "up": max(u for _k, _i, _w, u, _d in sized),
            "down": max(d for _k, _i, _w, _u, d in sized)}


def _emit_header(recipe, cfg, alloc, put, grid_w, tid0, base0, m=None, y0=0.0,
                 behind=None, x0=0.0):
    """The generic reaction that sits above a scope table.

    A deliberately small linear row -- structures separated by "+" or by one
    arrow with captions -- rather than a call into layout_linear, which returns a
    finished document and cannot be embedded as a body fragment.

    `m` is measure_header()'s record. Passed in by a caller that has already
    measured (which is every caller once the figure width is resolved before
    emit); measured here when omitted, so the two existing call sites stay as
    they are.

    Returns (height, next_tid, next_base). `y0` is the top of the block: 0.0 for
    a scheme at the top of the figure, which is every arrangement but the ones
    where `order` puts something above it.
    """
    m = measure_header(recipe, cfg) if m is None else m
    if m is None:
        return 0.0, tid0, base0
    S = recipe["structures"]
    tid, base = tid0, base0
    sized, up, down = m["sized"], m["up"], m["down"]
    axis = y0 + up
    # `x0` is where this block's BAND starts, and `grid_w` how wide the band is.
    # Both are 0.0 and the figure width on every figure whose scheme spans the
    # whole thing, which is all of them but one: a panel on the scheme's left
    # takes the left-hand strip, and the scheme then centres in what is LEFT
    # rather than in the page. Centring in the page is what put the reaction
    # under the panel instead of beside it.
    x = x0 + max(0.0, (grid_w - m["natural"]) / 2.0)

    for _i, (kind, it, w, _u, _d) in enumerate(sized):
        cx = x + w / 2.0
        if kind == "text":
            # Read, not re-measured: measure_header put the box in m["boxes"].
            tw, th, size = m["boxes"][_i]
            ty = axis - th / 2.0
            put(C.styled_text(tid, cx, ty + ASCENT_EM * size,
                              [str(it["text"])], size=size, align="center"),
                Rect(cx - tw / 2, ty, tw, th), "hdr:text", "caption")
            tid += 1
            if it.get("id"):
                nw, nh = block_box([str(it["id"])], cfg["num_size"])
                ny = ty + th + cfg["num_gap"]
                put(C.styled_text(tid, cx, ny + ASCENT_EM * cfg["num_size"],
                                  [str(it["id"])], size=cfg["num_size"],
                                  bold=True, align="center"),
                    Rect(cx - nw / 2, ny, nw, nh), "hdr:%s.num" % it["id"], "label")
                tid += 1
        elif kind == "plus":
            pw, ph = block_box(["+"], cfg["plus_size"])
            put(C.styled_text(tid, cx, axis - ph / 2.0 + ASCENT_EM * cfg["plus_size"],
                              ["+"], size=cfg["plus_size"], align="center"),
                Rect(cx - pw / 2, axis - ph / 2, pw, ph), "hdr:plus", "caption")
            tid += 1
        elif kind == "arrow":
            a = it["arrow"]
            ln = w - 2 * cfg["gap"]
            put(C.arrow(90500 + tid, cx, axis, ln),
                Rect(cx - ln / 2, axis - 2, ln, 4), "hdr:arrow", "arrow")
            for key, sign in (("above", -1), ("below", 1)):
                if not a.get(key):
                    continue
                cw, ch = block_box(a[key], cfg["cap_size"])
                gap = cfg["txt_gap_above"] if sign < 0 else cfg["txt_gap_below"]
                ty = axis - gap - ch if sign < 0 else axis + gap
                put(C.styled_text(tid, cx, ty + ASCENT_EM * cfg["cap_size"], a[key],
                                  size=cfg["cap_size"], align="center"),
                    Rect(cx - cw / 2, ty, cw, ch), "hdr:arrow." + key, "caption")
                tid += 1
        else:
            frag, bb = S[it["s"]]
            sw, sh = bb[2] - bb[0], bb[3] - bb[1]
            sx, sy = cx - sw / 2.0, axis - sh / 2.0
            _fx, _dx = C.place_with(
                frag, (recipe.get("decorations") or {}).get(it["s"]),
                base, sx, sy)
            put(_fx, Rect(sx, sy, sw, sh), "hdr:%s" % it["s"], "structure")
            if _dx and behind is not None:
                behind.append(_dx)
            base += 1000
            nw, nh = block_box([str(it["s"])], cfg["num_size"])
            ny = sy + sh + cfg["num_gap"]
            put(C.styled_text(tid, cx, ny + ASCENT_EM * cfg["num_size"],
                              [str(it["s"])], size=cfg["num_size"], bold=True,
                              align="center"),
                Rect(cx - nw / 2, ny, nw, nh), "hdr:%s.num" % it["s"], "label")
            tid += 1
        x += w
    return up + down, tid, base


def _group_cell(parts, start, box, cfg, gid):
    """Ctrl+G: fold everything the cell just emitted into one <group>.

    `parts[start:]` is exactly what this cell wrote, because `put` appends in
    emit order and nothing else runs between the two call sites. Taking the slice
    rather than threading a sink through _emit_cell keeps that function byte-for-
    byte the one the scope goldens were cut from.

    RESERVATIONS ARE UNTOUCHED. Each piece is still reserved on its own, so the
    Allocator keeps finding overlaps at the granularity it always did -- a group
    is a fact about the FILE, not about the packing. Nothing is reserved for the
    group itself either: its box is the union of boxes already reserved, and
    reserving it again would report every cell as overlapping itself.

    DECORATIONS STAY OUT, and this is a real limit rather than an oversight. A
    halo/dot/locant goes to `behind`, which _finish writes BEFORE parts so stamp_z
    gives every glow a lower Z than every structure. Pulling one into the cell's
    group would move it to the group's place in document order and raise its Z
    over the neighbouring cell it overhangs -- trading a correct print for a
    convenient drag. So a decorated cell drags without its glow.

    A cell that emitted ONE object is left alone: Ctrl+G on a single object does
    nothing in ChemDraw either.
    """
    if not cfg.get("group_cells") or box is None or len(parts) - start < 2:
        return gid
    chunk = parts[start:]
    del parts[start:]
    parts.append(C.group(gid, chunk, bbox=[round(v, 2) for v in box]))
    return gid + 1


def _emit_cell(m, cell_x, cell_w, top, band_h, cfg, put, tid, base, lalign,
               behind=None):
    """One product: its structure, then the label stack under it.

    Lifted out of render_grid verbatim so the sectioned packer can reuse it. The
    scope goldens are what prove the lift changed nothing -- the order of the two
    put() calls and of the id handouts is the output, byte for byte.

    Returns (tid, base, box), where box is the union of every rect this cell
    reserved -- the one thing _group_cell cannot work out for itself, and read
    off the same numbers that were reserved rather than measured a second time.
    """
    _bx: list = []

    def _seen(x, y, w, h):
        _bx.append((x, y, x + w, y + h))
    cx = cell_x + cell_w / 2.0
    sid = m["p"]["s"]

    # When the labels are left-anchored the structure is too, so the two
    # share one left edge instead of drifting apart.
    sx = (cell_x + cfg["cell_pad_x"] if lalign == "left"
          else cx - m["sw"] / 2.0)
    # "bottom" sits every structure in a row on one baseline, which is the
    # usual scope-table look and kills the ragged gap a short structure gets
    # when it is centred inside a band sized by a tall neighbour.
    sy = (top + band_h - m["sh"] if cfg.get("struct_valign") == "bottom"
          else top + (band_h - m["sh"]) / 2.0)
    _fx, _dx = C.place_with(m["frag"], m.get("deco"), base, sx, sy)
    put(_fx, Rect(sx, sy, m["sw"], m["sh"]), "cell:%s" % sid, "structure")
    _seen(sx, sy, m["sw"], m["sh"])
    if _dx and behind is not None:
        behind.append(_dx)
    base += 1000

    # label stack: starts at the SAME y in every cell, so the numbers line up
    y = top + band_h
    for kind, lines, size, gap, bold in m["stack"]:
        w, h = block_box(lines, size)
        y += gap
        if lalign == "left":
            lx = cell_x + cfg["cell_pad_x"]
            anchor, ax = "left", lx
        else:
            anchor, ax, lx = "center", cx, cx - w / 2.0
        put(C.styled_text(tid, ax, y + ASCENT_EM * size, lines,
                          size=size, bold=bold, align=anchor),
            Rect(lx, y, w, h),
            "cell:%s.%s" % (sid, kind),
            "label" if kind == "num" else "caption")
        _seen(lx, y, w, h)
        tid += 1
        y += h
    if m.get("legend"):
        # BELOW the stack, as published figures do: a Markush cell that carries a
        # value table usually has no number of its own -- the numbers are IN the table
        # -- so a caller omits `id` and the stack above is empty. No switch for
        # that; leaving the key out is the switch.
        lw = m["legend_w"]
        lx = (cell_x + cfg["cell_pad_x"] if lalign == "left" else cx - lw / 2.0)
        ly = y + cfg["legend_gap"]
        lh = m["legend_h"] - cfg["legend_gap"]
        put(_LG.place_legend(m["legend"], base, lx, ly), Rect(lx, ly, lw, lh),
            "cell:%s.legend" % sid, "caption")
        _seen(lx, ly, lw, lh)
        base += 1000
    box = ((min(b[0] for b in _bx), min(b[1] for b in _bx),
            max(b[2] for b in _bx), max(b[3] for b in _bx)) if _bx else None)
    return tid, base, box


def _emit_footnotes(recipe, cfg, alloc, put, tid):
    """Running prose under the table, wrapped to the grid width, one left edge.

    Centring them and back-solving for the left edge makes the margin ragged,
    because block_box's width estimate carries a little error and every footnote
    has a different width -- so anchor left.
    """
    notes = recipe.get("footnotes") or []
    if not notes:
        return tid
    grid = alloc.bbox()
    left = grid.x
    y = grid.bottom + cfg["footnote_gap"]
    wrapped = []
    for fn in notes:
        wrapped.extend(_wrap(fn, cfg["footnote_width"] or grid.w,
                             cfg["footnote_size"]))
    for j, fn in enumerate(wrapped):
        w, h = block_box([fn], cfg["footnote_size"])
        put(C.styled_text(tid, left, y + ASCENT_EM * cfg["footnote_size"],
                          [fn], size=cfg["footnote_size"], align="left"),
            Rect(left, y, w, h), "footnote:%d" % j, "caption")
        tid += 1
        y += h + cfg["footnote_lead"]
    return tid


# The two families of box, kept apart because they are emitted to TWO DIFFERENT
# PLACES: outlines are drawn last (on top of the content, which is harmless --
# it is only a ring of line), fills are drawn first (under the content,
# otherwise they are a block of colour covering the structure). A file saved by
# ChemDraw itself contains both families, spelled RectangleType="Filled" +
# color=<index>; there is no FillType attribute.
_BOX_OUTLINE = ("Plain", "Dashed", "RoundEdge", "RoundEdge Dashed")
_BOX_FILL = ("Filled", "Shaded", "RoundEdge Filled", "RoundEdge Shaded")


def _box_targets(alloc, cap: int = 8) -> str:
    """What this figure will accept as a box's `around`, read off its own ledger.

    A FIXED VOCABULARY WOULD BE A LIE HERE. Which names exist depends on what the
    figure turned out to be: a table-only figure has no 'grid', an unsectioned one
    has no 'section:1', and the cells are named after the caller's own compound
    numbers. So this reads the reservations rather than reciting the documentation,
    and a caller who mistypes '3a' as '3A' is told the ids that are actually there.

    CELLS ARE CAPPED, because a 40-product scope would otherwise bury the two
    block names under forty compound numbers. Header species are not: there are
    never many, and the ligand beside a scheme is exactly the one people look for.
    """
    # kind == "structure", NOT just the name prefix. The header also reserves
    # `hdr:arrow`, `hdr:plus` and `hdr:rule` -- its furniture -- and offering those
    # as "header species" would send a caller looking for a ligand to the plus
    # sign. (They still RESOLVE if someone asks; a box round the arrow is a
    # coherent thing to want. They are just not what this list is advertising.)
    res = list(alloc.reservations())
    regions = sorted({r.owner for r in res
                      if r.kind == "region" and "." not in r.owner})
    hdrs = sorted({r.owner for r in res
                   if r.owner.startswith("hdr:") and "." not in r.owner
                   and r.kind == "structure"})
    cells = sorted({r.owner for r in res
                    if r.owner.startswith("cell:") and "." not in r.owner
                    and r.kind == "structure"})
    bits = []
    if regions:
        bits.append("blocks and sections: " + ", ".join(map(repr, regions)))
    if hdrs:
        bits.append("header species: " + ", ".join(map(repr, hdrs)))
    if cells:
        shown = ", ".join(map(repr, cells[:cap]))
        bits.append("cells: %s%s" % (shown, " and %d more" % (len(cells) - cap)
                                     if len(cells) > cap else ""))
    if not bits:
        return "This figure has nothing on the ledger yet."
    return ("This figure has %s -- plus 'figure' for the lot. A cell or a header "
            "species is framed WITH its compound number and yield."
            % "; ".join(bits))


def _emit_boxes(recipe, cfg, alloc, put, behind=None):
    """Frames drawn AROUND things the figure has already placed.

    `around` names a region the packer registered rather than a rectangle the
    caller measured, and that is the whole point: the catalyst box in a paper sits
    around the header whatever the header turned out to be, and a caller who had
    to supply coordinates would be measuring the engine's output by hand. This is
    the same lookup `panels` already does to place itself next to "block:grid".

    Names, all of them owner paths the packer registers:
        figure                 everything on the ledger  (alloc.bbox)
        scheme | grid | table  one block
        section:N              the Nth section, counting from 0
        hdr:<key>              ONE species in the generic reaction
        cell:<id>              ONE product cell

    THE LAST TWO ARE THE EASY ONES TO MISS. The ledger is keyed by owner path and
    `find` takes any of them, so `hdr:L4` resolves like any block name -- and they
    cover the one thing published figures ask for most after a block (the ligand
    beside the arrow; one cell ringed as special).

    <key> IS THE CALLER'S OWN ID when the header item gave one -- `{"smiles": ...,
    "id": "L4"}` is `hdr:L4`. With no id the engine falls back to its own structure
    key and the figure prints THAT as the compound number, so a caller who wants to
    name a species has already had to name it for the figure to read right.

    A CELL OR A SPECIES IS FRAMED WITH ITS LABELS, which is why the resolution
    unions `find_all` instead of taking `find`: the structure is `cell:3a` and the
    compound number and yield hang off `cell:3a.num` / `cell:3a.yield`. Blocks and
    sections are unaffected -- a region covers its own contents by construction, so
    the union of a region and its children is the region.

    STILL NOT BUILT: a box round SEVERAL cells (the "unsuccessful examples"
    sub-box of published figures). That is a range, not a name, and the answer
    for now is to put those cells in their own `section` and box that.

    DRAWN AFTER, SO IT SITS ON TOP, and that is safe for the kinds offered here --
    they are outlines. It is NOT safe for "Filled" or "Shaded", which is why those
    are refused: a filled box emitted last would cover the thing it was drawn
    around, and putting it behind instead needs the decorations list the cells use.
    """
    boxes = recipe.get("boxes") or []
    if not boxes:
        return
    gid = cfg["box_gid"]
    for i, b in enumerate(boxes):
        if not isinstance(b, dict) or not b.get("around"):
            raise ValueError("boxes[%d]: needs 'around' -- what the box goes "
                             "round. Try 'figure', 'grid', 'scheme', 'table' or "
                             "'section:0'." % i)
        what = str(b["around"])
        kind = str(b.get("style") or "Plain")
        if kind not in _BOX_OUTLINE and kind not in _BOX_FILL:
            raise ValueError(
                "boxes[%d].style %r: an outline (%s) or a fill (%s)."
                % (i, kind, ", ".join(map(repr, _BOX_OUTLINE)),
                   ", ".join(map(repr, _BOX_FILL))))
        fill = kind in _BOX_FILL
        if fill and behind is None:
            raise ValueError(
                "boxes[%d].style %r: this packer cannot draw a fill -- it has no "
                "`behind` list, so the panel would land ON TOP of what it is "
                "meant to sit under. Use an outline here." % (i, kind))
        pad = float(b.get("pad", cfg["box_pad"]))
        if what == "figure":
            r = alloc.bbox()
        else:
            owner = what if ":" in what else "block:" + what
            # find_all, NOT find, and the difference is the whole point for a
            # single cell. The ledger's owner grammar is dotted: a product cell
            # reserves its structure as `cell:3a` and hangs its compound number and
            # yield off `cell:3a.num` / `cell:3a.legend`. find() returns the
            # structure alone, so a box asked for one cell would frame the drawing
            # and leave the number and the yield hanging outside it -- which is not
            # what anyone means by "box that cell".
            #
            # For a block or a section the two agree: a region covers its own
            # contents by construction (Allocator.region), so the union of the
            # region and its children IS the region, and a box round a block comes
            # out the same either way. Checked by a regression test, not assumed.
            parts = alloc.find_all(owner)
            if not parts:
                raise ValueError(
                    "boxes[%d].around %r: nothing on the ledger is called %r. %s"
                    % (i, what, owner, _box_targets(alloc)))
            r = parts[0].rect
            for p in parts[1:]:
                r = r.union(p.rect)
        x0, y0 = r.x - pad, r.y - pad
        x1, y1 = r.right + pad, r.bottom + pad
        xml = C.rect(gid + i, x0, y0, x1, y1, kind=kind, color=b.get("color"))
        if fill:
            # A FILL GOES IN `behind`, NOT THROUGH put(). Everything put() emits
            # lands after the cells, and _finish stamps Z in emission order -- so
            # a filled rectangle pushed through put() would be a solid colour laid
            # OVER the structures it was meant to sit under. `behind` is the list
            # the halos already use for exactly this reason, and it is the same
            # place because it is the same requirement: under the ink, not over.
            #
            # Reserved directly rather than through put(), because put() also
            # appends to `parts`. The reservation still has to happen or the page
            # will not grow to hold a panel that sticks out past the cells.
            behind.append(xml)
            alloc.reserve(Rect(x0, y0, x1 - x0, y1 - y0), "box:%d" % i, "decoration")
        else:
            put(xml, Rect(x0, y0, x1 - x0, y1 - y0), "box:%d" % i, "decoration")


def _emit_title(recipe, cfg, alloc, put, tid):
    """The figure's own title, ABOVE everything it has drawn so far.

    WRITTEN LAST AND PLACED FIRST, which is what makes it simple: the title folds
    to the FIGURE's width, and the figure's width is not known until the figure is
    drawn. Emitting it at the end and walking UP from alloc.bbox().y sidesteps the
    circularity entirely -- negative y is fine, finalize() translates the whole
    document into the positive quadrant afterwards. It is exactly the footnote
    emitter with the sign flipped, and it is placed next to it for that reason.

    LEFT, like the footnotes, and for the footnotes' stated reason: block_box's
    width estimate carries a little error, so centring and back-solving for the
    left edge makes the margin ragged when lines differ in length.
    """
    title = recipe.get("title")
    if not (title and str(title).strip()):
        return tid
    size = cfg["title_size"]
    lead = cfg["title_lead"] if cfg["title_lead"] is not None else cfg["footnote_lead"]
    gap = cfg["title_gap"] if cfg["title_gap"] is not None else cfg["footnote_gap"]
    box = alloc.bbox()
    lines = _wrap(str(title), box.w, size)
    heights = [block_box([l], size)[1] for l in lines]
    total = sum(heights) + (len(lines) - 1) * lead
    y = box.y - gap - total
    for j, ln in enumerate(lines):
        w, h = block_box([ln], size)
        put(C.styled_text(tid, box.x, y + ASCENT_EM * size, [ln],
                          size=size, align="left"),
            Rect(box.x, y, w, h), "title:%d" % j, "caption")
        tid += 1
        y += h + lead
    return tid


def _emit_section_band(title, x0, x1, y, cfg, put, tid, gid, idx,
                       style=None, letter=None, fill=None, text_color=None,
                       behind=None):
    """The section band. Returns its height; the title's vertical centre is where
    any lines sit.

    THREE SHAPES, out of at least six that published figures use:

      rule      a short lead line, the title, then a line to the right edge.
                Deliberately NOT one rule struck through the text -- three objects
                with a break either side. This engine's original and the default,
                so every existing figure is byte-identical.
      plain     the title alone, e.g. two blocks headed "Condition A" /
                "Condition B" with no rules at all.
      centered  a line, the title centred between x0 and x1, a line. The most
                common shape in published figures.

    `letter` is the bold a / b / c published figures put at the left margin,
    drawn at x0 ahead of whatever the style does. It is a SEPARATE argument rather than part of
    the title because it is anchored to the margin while the title may be centred --
    gluing them would drag the letter into the middle.
    """
    size = cfg["section_size"] or cfg["num_size"]
    # BOLD ON BOTH SIDES OF THE MEASUREMENT. The title is DRAWN bold four lines
    # down, so it has to be MEASURED bold too -- the right-hand rule starts at
    # tx + w. Measured regular on "with other (hetero)aryl groups" at 11 pt:
    # 143.68 pt reserved against 156.48 pt actually drawn, so the rule would
    # begin 12.8 pt early and strike through the last two characters. Two
    # independent readings of that 156: this pipeline's own
    # text_metrics.width(..., bold=True) says 156.45, and PIL against the system
    # Arial Bold says 156.48 -- so the metrics model is right as long as the call
    # asks for bold. It scales with the string, ~8.9%, so a longer title loses more.
    #
    # The docstring above says this band is deliberately NOT a rule struck through
    # the text. Without this keyword it would be one, on every sectioned figure --
    # and block_box's own docstring warns that a bold string measured with the
    # regular table comes out up to 6.7 pt narrow.
    #
    # THE ALLOCATOR CANNOT SEE IT: the title would reserve the 143.68 box and the
    # rule the space after it, so the two never overlap on the ledger. The ink
    # overlaps because the reservation is too small.
    w, h = block_box([title], size, bold=True)
    mid = y + h / 2.0
    lw = cfg["section_lw"]
    gap = cfg["section_text_gap"]
    style = style or cfg.get("section_style") or "rule"
    if style not in ("rule", "plain", "centered", "pill"):
        raise ValueError(
            "section style %r: known styles are 'rule' (lead line + title + rule), "
            "'plain' (title only), 'centered' (line, centred title, line) and "
            "'pill' (the title on a filled rounded badge, lines either side). "
            "A pill takes its colour from the marker's `color`." % (style,))
    if fill is not None and style != "pill":
        raise ValueError(
            "section %r asked for a `color` but style is %r. A colour on a section "
            "is the PILL's fill -- there is nowhere else on a rule or a plain "
            "heading to put one. Set style 'pill', or drop the colour."
            % (title, style))
    if style == "pill" and behind is None:
        raise ValueError(
            "section style 'pill': this packer has no `behind` list, so the badge "
            "would be painted OVER its own title. Use 'centered'.")

    # THE LETTER FIRST, at the margin, whatever the style does with the title.
    if letter:
        lw_, lh_ = block_box([str(letter)], size, bold=True)
        put(C.styled_text(_SECTION_LETTER_TID + idx, x0, y + ASCENT_EM * size,
                          [str(letter)], size=size, bold=True, align="left"),
            Rect(x0, y, lw_, lh_), "section:%d.letter" % idx, "label")
        x0 = x0 + lw_ + gap
        h = max(h, lh_)
        mid = y + h / 2.0

    if style == "plain":
        put(C.styled_text(tid, x0, y + ASCENT_EM * size, [title],
                          size=size, bold=True, align="left"),
            Rect(x0, y, w, h), "section:%d.title" % idx, "label")
        return h

    if style == "pill":
        # Same placement as centered; the only difference is a filled rounded
        # rectangle under the title. The rectangle goes to `behind`, for the same
        # reason as the boxes' fills and the halos: drawn over the ink it is a
        # block of colour covering the text. The text itself still goes through
        # put(), so it lands after the badge.
        px, py = cfg["pill_pad_x"], cfg["pill_pad_y"]
        tx = x0 + (x1 - x0 - w) / 2.0
        bx0, by0, bx1, by1 = tx - px, y - py, tx + w + px, y + h + py
        behind.append(C.rect(_SECTION_PILL_GID + idx, bx0, by0, bx1, by1,
                             kind="RoundEdge Filled", color=fill,
                             radius=cfg["pill_radius"]))
        put(C.styled_text(tid, tx, y + ASCENT_EM * size, [title],
                          size=size, bold=True, align="left", color=text_color),
            Rect(bx0, by0, bx1 - bx0, by1 - by0), "section:%d.title" % idx, "label")
        # The lines stop at the badge's edge, not the title's -- the badge is
        # wider than the text, so a line drawn up to the text would run into it.
        if bx0 - x0 > lw + gap:
            put(C._line(gid, x0, mid, bx0 - gap, mid, lw),
                Rect(x0, mid - lw / 2.0, bx0 - gap - x0, lw),
                "section:%d.lead" % idx, "decoration")
        if x1 - (bx1 + gap) > lw:
            put(C._line(gid + 1, bx1 + gap, mid, x1, mid, lw),
                Rect(bx1 + gap, mid - lw / 2.0, x1 - bx1 - gap, lw),
                "section:%d.rule" % idx, "decoration")
        return by1 - by0

    if style == "centered":
        tx = x0 + (x1 - x0 - w) / 2.0
        if tx - x0 > lw + gap:
            put(C._line(gid, x0, mid, tx - gap, mid, lw),
                Rect(x0, mid - lw / 2.0, tx - gap - x0, lw),
                "section:%d.lead" % idx, "decoration")
        put(C.styled_text(tid, tx, y + ASCENT_EM * size, [title],
                          size=size, bold=True, align="left"),
            Rect(tx, y, w, h), "section:%d.title" % idx, "label")
        rx = tx + w + gap
        if x1 - rx > lw:
            put(C._line(gid + 1, rx, mid, x1, mid, lw),
                Rect(rx, mid - lw / 2.0, x1 - rx, lw),
                "section:%d.rule" % idx, "decoration")
        return h

    lead = min(cfg["section_lead"], max(0.0, x1 - x0))
    put(C._line(gid, x0, mid, x0 + lead, mid, lw),
        Rect(x0, mid - lw / 2.0, lead, lw), "section:%d.lead" % idx, "decoration")
    tx = x0 + lead + gap
    put(C.styled_text(tid, tx, y + ASCENT_EM * size, [title],
                      size=size, bold=True, align="left"),
        Rect(tx, y, w, h), "section:%d.title" % idx, "label")
    rx = tx + w + gap
    if x1 - rx > lw:                      # a table narrower than its own title
        put(C._line(gid + 1, rx, mid, x1, mid, lw),
            Rect(rx, mid - lw / 2.0, x1 - rx, lw),
            "section:%d.rule" % idx, "decoration")
    return h


def _pack_section(sec, recipe, cfg, usable):
    """One section's own measurements, column count and rows.

    Each section sizes its own cells -- two sections may end up with different
    column counts and different widths, which is deliberate: the eye compares
    within a section, and forcing one width on the whole table makes an aliphatic
    section of small chains carry the gutters of an aryl section.

    Rows are LEFT-ALIGNED, never centred: a short row at the end of a section (and
    at the end of the table) is set flush left in print.
    """
    # DECORATIONS COME ALONG. measure() looks a cell's halo/dot up in
    # `recipe["decorations"]` (see the "deco" key it builds), so this sub-recipe
    # has to carry it as well as structures + products -- without it, on the
    # sectioned path every cell's deco comes back None and every glow and dot
    # vanishes, silently, with no warning anywhere. Measured: two haloed cells
    # emit 2 <graphic GraphicType="Orbital"> with no section marker and, without
    # the key, 0 with one.
    sub = {"structures": recipe["structures"], "products": sec["products"],
           "decorations": recipe.get("decorations")}
    meas, nat_w = measure(sub, cfg)
    if sec.get("row_sizes"):
        sizes = list(sec["row_sizes"])
        cols = max(sizes)
    else:
        fits = max(1, int((usable or 540.0) // nat_w))
        n_rows = max(1, -(-len(meas) // fits))
        cols = max(1, min(-(-len(meas) // n_rows), len(meas)))
        sizes = []
        left = len(meas)
        while left > 0:
            sizes.append(min(cols, left))
            left -= cols
    cell_w = nat_w
    if cfg.get("fill_width", True) and cols > 1 and usable:
        cell_w = min(max(nat_w, usable / cols), nat_w * cfg["max_spread"])
    rows, i = [], 0
    for n in sizes:
        rows.append(list(range(i, i + n)))
        i += n
    return {"title": sec.get("title"), "meas": meas, "rows": rows,
            "cell_w": cell_w, "cols": cols,
            "style": sec.get("style"), "letter": sec.get("letter"),
            "color": sec.get("color"), "text_color": sec.get("text_color")}


def _render_sectioned(recipe, header_prefix, cfg, page_mode, alloc, usable, hm,
                      notes=None, shift=None):
    """The grid with section bands in it.

    A separate packer from the plain path on purpose. Sections change the packing
    model -- per-section cell widths, left-aligned short rows -- and the plain path
    has six byte-exact goldens on it, so the two do not share a code path where a
    change to one could move the other.

    THE TABLE IS BUILT HERE TOO. Without a _TableBlock on this path, `sections` +
    `table` would return a figure with no table in it and say nothing. Every rule
    in the figure spans W = the wider of the section
    grid and the table; the cells stay flush left, which is what _pack_section
    already decided about short rows.
    """
    plans = [_pack_section(s, recipe, cfg, usable) for s in recipe["sections"]]
    grid_w = max(p["cell_w"] * p["cols"] for p in plans)

    parts, behind = [], []   # decorations first, see _finish
    tid, base, gid = 90000, 1000, 96100
    cgid = _CELL_GID

    def put(xml, rect, owner, kind):
        res = alloc.reserve(rect, owner, kind)
        alloc.check_emit(res.rid, owner)
        parts.append(xml)
        return res

    order = _resolve_order(cfg, {"scheme": hm is not None, "grid": True,
                                 "table": bool(recipe.get("table"))})
    enclose = (cfg.get("frame") == "enclose" and "scheme" in order
               and "table" in order
               and order.index("scheme") < order.index("table"))
    blk = (_TableBlock(recipe, cfg, header_prefix, enclose=enclose)
           if recipe.get("table") else None)
    if blk is not None:
        header_prefix = blk.header_prefix
    W = max(grid_w, hm["natural"] if hm is not None else 0.0,
            blk.w if blk is not None else 0.0)
    if blk is not None:
        blk.fit(W, notes)
    below = blk is not None and order.index("table") < order.index("grid")

    y = 0.0
    if hm is not None:
        hh, tid, base = _emit_header(recipe, cfg, alloc, put, W, tid, base, hm,
                                     behind=behind)
        ry = hh + cfg["header_gap"] * 0.5
        put(C._line(96000, 0.0, ry, W, ry, 0.8),
            Rect(0.0, ry - 0.4, W, 0.8), "hdr:rule", "decoration")
        y = hh + cfg["header_gap"]
        if enclose:
            fw = TL._cfg(blk.cfg, recipe["table"].get("style"))["rule_top"] or 1.2
            fy = -cfg["header_gap"] * 0.5
            put(C._line(96001, 0.0, fy, W, fy, fw),
                Rect(0.0, fy - fw / 2, W, fw), "frame:top", "decoration")

    # "below" means the PANEL goes under the table, with the scheme staying on
    # top -- the same reading the plain path gives it.
    table_y = y
    if below:
        y += blk.h + cfg["table_gap"]

    lalign = cfg.get("label_align", "center")
    # THE BLOCK REGIONS. The main packer publishes block:scheme / block:grid /
    # block:table so panels have somewhere to hang off, and the sectioned packer
    # -- a separate function on purpose, see its docstring -- needs the same three
    # lines. Without them, anything that asks the ledger where a block is fails:
    # `boxes` around "grid" would work on an unsectioned figure and raise on a
    # sectioned one.
    if hm is not None:
        alloc.region("block:scheme", Rect(0.0, 0.0, W, hm["up"] + hm["down"]))
    _grid_y0 = y

    for si, plan in enumerate(plans):
        _sec_y0 = y
        if plan["title"]:
            if si or recipe.get("header"):
                y += cfg["section_gap_above"]
            y += _emit_section_band(plan["title"], 0.0, W, y, cfg,
                                    put, tid, gid, si,
                                    style=plan.get("style"),
                                    letter=plan.get("letter"),
                                    fill=plan.get("color"),
                                    text_color=plan.get("text_color"),
                                    behind=behind)
            tid += 1
            gid += 2
            y += cfg["section_gap_below"]
        meas, cw = plan["meas"], plan["cell_w"]
        for row in plan["rows"]:
            band_h = max(meas[i]["sh"] for i in row)
            # legend_h included: the value table is part of what hangs below
            # the structure, and a band sized without it overlaps the next row.
            stack_h = max(meas[i]["stack_h"] + meas[i]["legend_h"] for i in row)
            top = y + cfg["cell_pad_y"]
            for col, i in enumerate(row):
                _n0 = len(parts)
                tid, base, _cbox = _emit_cell(meas[i], col * cw, cw, top,
                                              band_h, cfg, put, tid, base,
                                              lalign, behind=behind)
                cgid = _group_cell(parts, _n0, _cbox, cfg, cgid)
            y += band_h + stack_h + 2 * cfg["cell_pad_y"]
        # ONE REGION PER SECTION, so a box can be asked for "section:0" without the
        # caller measuring anything. The main path already does this for its three
        # blocks; a section is the same kind of thing one level down. kind="region"
        # is tolerant against everything (see Allocator.region), so registering it
        # cannot turn a clean figure into a reported overlap, and a region emits no
        # ink -- every existing figure comes out byte-identical.
        alloc.region("section:%d" % si, Rect(0.0, _sec_y0, W, max(0.0, y - _sec_y0)))
    alloc.region("block:grid", Rect(0.0, _grid_y0, W, max(0.0, y - _grid_y0)))

    # AFTER every cell: a connector asks the ledger where they landed.
    _aid, tid = _emit_connectors(recipe, cfg, alloc, put, 90600, tid)

    # Panels, placed relative to a block region or to an earlier panel.
    if recipe.get("panels"):
        import panels as _PN                 # deferred: _PN imports this module
        _PN.place(recipe["panels"], recipe, cfg, alloc, parts.append,
                  notes if notes is not None else [],
                  gid=97500, tid=tid + 400, base=base + 40000)

    if blk is not None:
        blk.emit((W - blk.w) / 2.0, table_y if below else y + cfg["table_gap"], put)

    tid = _emit_footnotes(recipe, cfg, alloc, put, tid)
    # AFTER the footnotes, and the order is load-bearing rather than tidy: the
    # title folds to the FIGURE's width, so everything that can set that width
    # has to be on the ledger first. Footnotes cannot widen it (they fold to the
    # grid), but reading the box after them is what makes that a fact rather
    # than an assumption.
    # BEFORE the title and AFTER everything else: a box asked for "figure" should
    # enclose the drawing, and the heading belongs outside it -- which is where
    # published figures put it.
    _emit_boxes(recipe, cfg, alloc, put, behind)
    tid = _emit_title(recipe, cfg, alloc, put, tid)
    # _finish, not a bare finalize: a table with a shaded band needs stamp_z, and
    # this path can carry a table.
    return _finish(header_prefix, parts, alloc, page_mode, W, cfg["margin"],
                   notes, behind=behind, shift=shift, span=cfg.get("span"))


def _emit_connectors(recipe, cfg, alloc, put, aid, tid):
    """An arrow BETWEEN two product cells.

    A product grid has no connector concept of its own -- the topology is "no edges
    between the cells", which is true of a scope table and false of a figure where,
    say, 20 -[O]-> 21 is a small transformation drawn inside the grid, each side
    carrying its own number and yield.

    WHERE THE CELLS ARE is asked of the ledger rather than tracked in a second
    variable: alloc.find("cell:3a") is exactly the query the ledger exists to
    answer. layout_linear keeps a parallel `struct_pos`
    for the same purpose, which is the duplication this avoids.

    THE GUTTER IS NOT WIDENED. cell_pad_x is 20 either side, so two structures sit
    40 pt apart and an arrow with a caption may not fit. Making the column solver
    account for it is a change to how every scope table is sized, for a rare
    figure -- so this draws in the room that exists and RESERVES it, and the
    allocator says so when it does not fit. Report, not silence, and not a rewrite
    of the packer either.
    """
    def _half(res, going_down):
        """(hw, hh) for one end, taking the LABEL STACK into account when the run
        leaves downward through it.

        Same correction layout_linear makes (without it a connector drives
        through the compound number of the cell it leaves), and here it costs nothing
        to get right: every label is already reserved as `cell:3a.num` /
        `.yield` / `.sub`, so find_all on the cell's path IS the cell's full
        content box. linear keeps a parallel `struct_low` dict for the same answer.

        Only downward, for the same reason as there: a sideways run leaves through
        the left or right edge, where the stack is not.
        """
        hw, hh = res.rect.w / 2.0, res.rect.h / 2.0
        if going_down:
            low = max(r.rect.bottom for r in alloc.find_all(res.owner))
            hh = max(hh, low - res.rect.cy)
        return hw, hh

    for i, con in enumerate(recipe.get("connectors") or []):
        a, b = alloc.find("cell:%s" % con["from"]), alloc.find("cell:%s" % con["to"])
        own = "conn:%s->%s" % (con["from"], con["to"])
        r = A.route({"kind": "straight",
                     "f": (a.rect.cx, a.rect.cy),
                     "fh": _half(a, b.rect.cy > a.rect.cy),
                     "t": (b.rect.cx, b.rect.cy),
                     "th": _half(b, a.rect.cy > b.rect.cy),
                     "style": con.get("style") or {}},
                    aid=aid, margin=cfg["conn_margin"])
        aid += len(r["segs"]) + 1
        hx0, hy0, hx1, hy1 = r["head"]
        # ONE segment: kind "straight" returns exactly one, and there is no routed
        # kind here (see _GRID_CONN_KEYS -- two neighbours in a grid have nothing to
        # route around). Indexing says that; iterating with a break only hid it.
        put(r["segs"][0], Rect(min(hx0, hx1), min(hy0, hy1) - 2.0,
                               abs(hx1 - hx0) or 4.0, abs(hy1 - hy0) + 4.0),
            own, "arrow")
        mx, my = (hx0 + hx1) / 2.0, (hy0 + hy1) / 2.0
        # THE CAPTION ROTATES WITH THE ARROW, exactly as layout_linear's caption
        # placement note says: turning a left-to-right arrow 90 degrees maps
        # up->right and down->left, so a mostly-VERTICAL run puts "above" on its right and "below"
        # on its left -- the conditions-right / yield-left convention the literature
        # draws. Without this, "above" on a vertical connector lands in the label
        # stack of the cell it just left, which the ledger duly reported.
        _vertical = abs(hy1 - hy0) > abs(hx1 - hx0)
        _sides = [(k, s) for k, s in (("above", +1), ("below", -1)) if con.get(k)]
        if _vertical and _sides:
            # One shared top, so the two read as horizontally aligned -- centring
            # each on its own height puts the baselines a unit apart, because the
            # height model depends on which glyphs are present.
            _boxes = {k: block_box(con[k], cfg["cap_size"]) for k, _s in _sides}
            _top = my - max(h for _w, h in _boxes.values()) / 2.0
            for key, sgn in _sides:
                w, h = _boxes[key]
                cx = mx + sgn * (cfg["txt_gap_side"] + w / 2.0)
                put(C.styled_text(tid, cx, _top + ASCENT_EM * cfg["cap_size"],
                                  con[key], size=cfg["cap_size"], align="center"),
                    Rect(cx - w / 2, _top, w, h), own + "." + key, "caption")
                tid += 1
        else:
            for key, sgn in (("above", -1), ("below", +1)):
                if not con.get(key):
                    continue
                lines = con[key]
                w, h = block_box(lines, cfg["cap_size"])
                gap = cfg["txt_gap_above"] if sgn < 0 else cfg["txt_gap_below"]
                top = (min(hy0, hy1) - gap - h) if sgn < 0 else (max(hy0, hy1) + gap)
                put(C.styled_text(tid, mx, top + ASCENT_EM * cfg["cap_size"], lines,
                                  size=cfg["cap_size"], align="center"),
                    Rect(mx - w / 2, top, w, h), own + "." + key, "caption")
                tid += 1
    return aid, tid


# The three blocks a combined figure can hold, in the order they default to.
_BLOCKS = ("scheme", "grid", "table")


def _resolve_order(cfg, present):
    """The vertical block order, as a list of the names actually present.

    `order` replaces `panel_pos`, which could only say whether the
    structure panel came before or after the table, and so could not express a
    scheme above the table with a ligand grid below it.

    A name that is not in the figure is dropped: asking for a grid you did not
    pass is not an error. A block that IS in the figure and is NOT named is an
    error -- dropping a block silently is exactly the failure `order` exists to
    prevent, and an `order` that quietly loses a table would reintroduce it at the
    API level.
    """
    want = cfg.get("order")
    if want:
        seq = [str(b) for b in want]
        bad = [b for b in seq if b not in _BLOCKS]
        if bad:
            raise ValueError("order: unknown block(s) %s; the blocks are %s"
                             % (sorted(set(bad)), ", ".join(_BLOCKS)))
        if len(set(seq)) != len(seq):
            raise ValueError("order: %s names a block twice" % (seq,))
        missing = [b for b in _BLOCKS if present.get(b) and b not in seq]
        if missing:
            raise ValueError(
                "order: %s leaves out %s, which this figure has. A block left "
                "out of `order` would be dropped without a word; name it, or "
                "take it out of the recipe." % (seq, ", ".join(missing)))
        return [b for b in seq if present.get(b)]
    # panel_pos, the deprecated alias: it said where the PANEL went relative to
    # the table, with the scheme always on top.
    seq = (["scheme", "table", "grid"] if cfg.get("panel_pos") == "below"
           else ["scheme", "grid", "table"])
    return [b for b in seq if present.get(b)]


def _block_gap(a, b, cfg):
    """The break between two stacked blocks.

    header_gap between a scheme and the grid under it, table_gap anywhere else.
    The two are both 26.0 and GRID_DEFAULTS says so on purpose ("the
    same break as header-to-grid, because it is the same kind of break"), so this
    preserves every existing figure's arithmetic exactly while keeping the two
    names attached to the two ideas.
    """
    return cfg["header_gap"] if (a, b) == ("scheme", "grid") else cfg["table_gap"]


class _TableBlock:
    """recipe['table'], measured, with its colours already in the colortable.

    Built before anything is emitted, because a shaded band names its colour by
    an index into the colortable and the colortable is part of the header --
    which is written before the body exists. Same order draw_molecule uses for
    its highlight colours.
    """

    def __init__(self, recipe, cfg, header_prefix, enclose=None):
        self.table = recipe["table"]
        self.notes = recipe.get("table_notes")
        self.cfg = dict(recipe.get("table_cfg") or {})
        # `enclose` None means "read it off cfg". It is passed explicitly by a
        # caller that knows whether there is actually a SCHEME to enclose:
        # frame:"enclose" with no scheme would drop rule_top and leave nobody to
        # redraw it, so the table would come out with no top rule at all -- a
        # silent drop, one level down.
        if enclose is None:
            enclose = cfg.get("frame") == "enclose"
        if enclose:
            # The scheme and the table share ONE top rule, drawn above the
            # scheme so the scheme sits inside the table's frame.
            self.cfg["rule_top"] = None
        self.colors = None
        rgbs = TL.fills_used(self.table, self.cfg)
        if rgbs:
            header_prefix, idx = BT.extend_colortable(
                header_prefix, [list(c) for c in rgbs])
            self.colors = dict(zip(rgbs, idx))
        self.header_prefix = header_prefix
        m = TL.measure(self.table, self.notes, self.cfg)
        self.w, self.h = m["w"], m["h"]
        self.natural_w = m["w"]

    def fit(self, W, notes=None):
        """Pass 2 for the table: lay it out against the width the parent resolved.

        A no-op under width_mode "natural" (the default):
        the table keeps its own width and the parent centres it. Under "match" the
        columns absorb the slack, so both w and h can move -- h because the notes
        fold to the figure width and a wider table may need fewer lines of them.

        AND IT REPORTS A SHORT REACH. `stretch_max` caps the stretch at 1.30 x
        natural, so a narrow table under a wide scheme can ask for "match" and
        still come out narrower than the figure -- the cap wins and the table
        stays centred in the leftover white. That is the right LOOK (see
        table_layout.CFG), but it is not what the caller asked for, and the only
        way to notice it in the output is to measure the file. So it goes in
        `notes`, the same channel the page-overflow note uses.
        """
        cfg = TL._cfg_for(self.table, self.cfg)
        if cfg["width_mode"] != "match":
            return
        self.cfg["target_w"] = W
        m = TL.measure(self.table, self.notes, self.cfg)
        self.w, self.h = m["w"], m["h"]
        if notes is not None and self.natural_w > 0 and self.w < W - 0.01:
            notes.append(
                "width_mode 'match' asked the table to span the figure's %.0f pt, "
                "but stretch_max %.2f caps it at %.2f x its natural %.0f pt = "
                "%.0f pt, so it comes out %.0f pt short and stays centred. To "
                "reach %.0f pt the cap has to be at least %.2f: pass "
                "table_style={'stretch_max': %.2f}."
                % (W, cfg["stretch_max"], cfg["stretch_max"], self.natural_w,
                   self.w, W - self.w, W, W / self.natural_w,
                   math.ceil(W / self.natural_w * 100.0) / 100.0))

    def emit(self, x, y, put):
        """The whole table as ONE reservation.

        The body is text of known size with no structure in it, so
        there is nothing inside for a collision check to find. Registering the
        block keeps the footnotes and the page extent honest; registering every
        cell would add three hundred boxes that cannot collide with anything.
        """
        body, bb = TL.render(self.table, self.notes, self.cfg, self.colors,
                             x0=x, y0=y)
        put(body, Rect(bb[0], bb[1], bb[2] - bb[0], bb[3] - bb[1]),
            "table", "caption")


def _by_area(xmls):
    """`behind` sorted biggest-first, so a big panel cannot bury a small badge."""
    def area(x):
        m = re.search(r'BoundingBox="([-\d.eE ]+)"', x)
        if not m:
            return 0.0
        v = [float(t) for t in m.group(1).split()]
        return abs(v[2] - v[0]) * abs(v[3] - v[1]) if len(v) >= 4 else 0.0
    return sorted(xmls, key=area, reverse=True)


def _finish(header_prefix, parts, alloc, page_mode, W=None, margin=None,
            notes=None, behind=None, shift=None, span=None):
    """Join, stamp Z if there are bands, and finalize.

    stamp_z ONLY when a band exists: its own docstring says every other figure
    here carries no Z at all and renders correctly, and stamping those would
    rewrite twelve goldens to say something they never needed to say.

    Also the one place the page-overflow line is written. REPORT ONLY: every
    fixed page mode already grows rather than clips (cdxml_core.finalize, "the
    page NEVER clips"), and that invariant is shared by all six engines, so this
    says so instead of changing it.
    """
    if notes is not None and W and margin is not None:
        # WHOEVER SET THE LIMIT DECIDES WHAT THE SENTENCE SAYS. Under a page mode
        # the limit is a column and the page is worth naming; under `span` the
        # limit is a number the caller chose, and naming the page there sends them
        # to fix something that is not the constraint. Same overflow, same
        # channel, two different pieces of news.
        u = float(span) if span else C.usable_width(page_mode, margin)
        if u and W > u + 0.01:
            notes.append(
                ("the figure comes to %.0f pt wide, %.0f pt past the span of "
                 "%.0f pt it was given. Nothing is clipped, but it will overhang "
                 "whatever it was meant to sit inside. A cell is only ever "
                 "WIDENED to fill a span, never narrowed, so the grid is this "
                 "wide because its widest product is: the remedy is fewer "
                 "columns, or a smaller drawing in that cell."
                 % (W, W - u, u))
                if span else
                ("the figure comes to %.0f pt wide, %.0f pt past the usable width "
                 "of page %r (%.0f pt). The page grows to fit rather than clip, so "
                 "nothing is lost -- but it will not sit in that column. The "
                 "table's own remedy is split: 2; the panel's is fewer columns."
                 % (W, W - u, page_mode, u)))
            if _WIDEST_LEGEND:
                sid, lw = _WIDEST_LEGEND[0]
                notes.append(
                    "the widest thing in the grid is the value table under cell "
                    "%r (%.0f pt) -- every column is as wide as it, because a "
                    "scope table only reads as one when the cells are congruent. "
                    "Splitting that legend into more rows narrows the whole grid."
                    % (sid, lw))
    # SAID WHETHER OR NOT THE PAGE OVERFLOWED, which is the whole difference from
    # _WIDEST_LEGEND above -- see the note on _LEGEND_WIDE for the figure that
    # went 469 pt past the page without a word about the table that did it.
    if notes is not None:
        for _sid, _lw, _sw, _cause in _LEGEND_WIDE:
            _where = ""
            if _cause:
                _col, _cw, _txt = _cause
                _where = (" Column %d is %.0f of those points%s."
                          % (_col, _cw,
                             (", its longest entry %r" % _txt) if _txt else ""))
            notes.append(
                "the value table under cell %r came out %.0f pt wide against a "
                "%.0f pt structure (%.2fx). Drawn about as wide as its own "
                "structure the two read as one object, and EVERY column in the "
                "grid is at least as wide as this cell -- so a wide table widens "
                "the whole figure and can cost it a column.%s Nothing was "
                "changed and the drawing is right either way: the published fix "
                "is to draw the generic structure larger than an ordinary "
                "product, or to shorten that entry."
                % (_sid, _lw, _sw, _lw / _sw, _where))
    # Decorations first, so stamp_z gives every glow a lower Z than every
    # structure -- not just than its own. A halo overhangs its atoms by 0.8 bond
    # and cell_pad_x is 20, so it can reach into the NEXT cell, and "just before
    # its own structure" would then leave it on top of the neighbour.
    # BIGGEST FIRST, INSIDE `behind`. Everything here is under the ink, but they
    # are also under EACH OTHER in emission order -- and the order things get
    # appended in is the packer's convenience, not a stacking intent. Measured:
    # a section's pill badge is appended during the section loop and a `boxes`
    # panel round that same section is appended afterwards, so the panel (982 x
    # 414) was painted over the badge (139 x 15) and the badge vanished.
    #
    # Area, descending, stable. A halo's BoundingBox is degenerate by design (see
    # cdxml_core.disc) so it measures ~0 and lands last -- on top of the panels,
    # under the structures, which is exactly where a glow belongs.
    body = (("\n".join(_by_area(behind)) + "\n") if behind else "") + "".join(parts)
    if behind or "<ColoredMolecularArea" in body:
        body, _ = C.stamp_z(body)
    bb = alloc.bbox()
    return C.finalize(header_prefix, body, [bb.x, bb.y, bb.right, bb.bottom],
                      page_mode=page_mode, shift=shift)


def _render_no_grid(recipe, header_prefix, cfg, page_mode, alloc, hm, notes=None,
                    shift=None):
    """A figure with NO product grid: a scheme, a data table, or both.

    This is the optimisation figure -- a generic reaction over a table of
    conditions, with no products drawn between them -- and it is also a table
    with no scheme at all. Both are legal input; `rows` may be empty.

    EVERY BLOCK IS READ HERE. Emitting the table and returning would hand a caller
    who passed a header AND a table a document with the scheme silently missing,
    and `footnotes` would go the same way. Nothing would fail and nothing would be
    reported.

    NO SCHEME RULE IS DRAWN HERE. In the plain path that rule separates the
    scheme from the product grid under it. With no grid, the next thing below the
    scheme is the table's own top rule, which is already that separator; drawing
    both would put two rules 26 pt apart with nothing between them.
    """
    parts, behind = [], []   # decorations first, see _finish

    def put(xml, rect, owner, kind):
        res = alloc.reserve(rect, owner, kind)
        alloc.check_emit(res.rid, owner)
        parts.append(xml)
        return res

    order = _resolve_order(cfg, {"scheme": hm is not None, "grid": False,
                                 "table": bool(recipe.get("table"))})
    enclose = (cfg.get("frame") == "enclose" and "scheme" in order
               and "table" in order
               and order.index("scheme") < order.index("table"))
    blk = (_TableBlock(recipe, cfg, header_prefix, enclose=enclose)
           if recipe.get("table") else None)
    # A SIDE PANEL JOINS ITS BLOCK'S BAND, and is measured here because a panel
    # measured after the layout can only be hung off the outside. Measured on a
    # published figure hung that way: both blocks ended at 545 pt and the ligand
    # box ran 557 -> 673, so the figure was 673 pt wide with everything in the left
    # 545 and no two right edges alike; it also ran 72 pt below the scheme,
    # alongside the table's first rows. A panel is part of the figure, not a
    # column standing next to it.
    #
    # `pads[name]` is [left, right, tallest panel] in points, and it is what the
    # block GIVES UP rather than what the figure gains: the band is W minus the
    # pads, the panel sits in the strip, and the two share the figure's edge.
    #
    # LEFT/RIGHT ONTO A BLOCK, and nothing else. A panel anchored to another panel
    # stacks outside the band by construction, and above/below panels do not take a
    # block's width away -- both still overhang, and both are their own piece of
    # work. Two panels on the SAME side take the max, not the sum: they would land
    # on each other, and that is a collision for the ledger to report, not a
    # spacing this loop should quietly make room for.
    pads = {"scheme": [0.0, 0.0, 0.0], "grid": [0.0, 0.0, 0.0],
            "table": [0.0, 0.0, 0.0]}
    if recipe.get("panels"):
        import panels as _PN                 # deferred: _PN imports this module
        for _p in recipe["panels"]:
            _of, _side = _p.get("of", "grid"), _p.get("side", "right")
            if _of not in order or _side not in ("left", "right"):
                continue
            _pw, _ph = _PN.measure(_p, recipe, cfg)
            _i = 0 if _side == "left" else 1
            pads[_of][_i] = max(pads[_of][_i],
                                float(_p.get("gap", cfg["table_gap"])) + _pw)
            pads[_of][2] = max(pads[_of][2], _ph)

    # W, for the two blocks that can be here. Each is centred in its own BAND, so
    # the widest band-plus-pads sets the figure width and the narrower one sits in
    # the middle of what is left over.
    natural = {"scheme": hm["natural"] if hm else 0.0,
               "table": blk.w if blk is not None else 0.0}
    W = max([natural[nm] + pads[nm][0] + pads[nm][1] for nm in order] or [0.0])
    if blk is not None:
        blk.fit(W - pads["table"][0] - pads["table"][1], notes)

    heights = {"scheme": (hm["up"] + hm["down"]) if hm is not None else 0.0,
               "table": blk.h if blk is not None else 0.0}
    # A band is as tall as the taller of the block and its side panel, so the NEXT
    # block starts below both. Without this the panel is level with the scheme and
    # the table starts under the SCHEME, which on that same figure put the table's
    # first three rows alongside the ligand box.
    for name in order:
        heights[name] = max(heights[name], pads[name][2])
    top, y = {}, 0.0
    for i, name in enumerate(order):
        if i:
            y += _block_gap(order[i - 1], name, cfg)
        top[name] = y
        y += heights[name]
    # Each block becomes a REGION, so a panel can be placed relative to one. The
    # extents are the ones this loop just computed -- not measured again, which
    # would be two accounts of one number.
    for name in order:
        if heights[name] > 0.0:
            alloc.region("block:%s" % name,
                         Rect(pads[name][0], top[name],
                              W - pads[name][0] - pads[name][1], heights[name]))

    tid, base = 90000, 1000
    if hm is not None:
        hh, tid, base = _emit_header(
            recipe, cfg, alloc, put,
            W - pads["scheme"][0] - pads["scheme"][1], tid, base, hm,
            y0=top["scheme"], behind=behind, x0=pads["scheme"][0])
        if enclose:
            # The table's top rule, moved up so the scheme sits inside it. Same
            # idiom as the plain path, so the two agree on weight and offset.
            # Across the TABLE's band, because it is the table's rule: a panel
            # beside the scheme must not lengthen it.
            fw = TL._cfg(blk.cfg, recipe["table"].get("style"))["rule_top"] or 1.2
            ry = top["scheme"] - cfg["table_gap"] * 0.5
            rx0, rx1 = pads["table"][0], W - pads["table"][1]
            put(C._line(96001, rx0, ry, rx1, ry, fw),
                Rect(rx0, ry - fw / 2, rx1 - rx0, fw), "frame:top", "decoration")

    # Panels, placed relative to a block region or to an earlier panel. The SAME
    # three lines the other two packers end with: the block regions are published
    # forty lines up so a panel can be placed relative to one, and `boxes`, the
    # footnotes and the title are wired the same way.
    #
    # WHY IT MATTERS HERE. An optimisation table IS this path: one generic
    # reaction over a table of conditions, no product grid -- and 17 of 32
    # surveyed optimisation figures carry a ligand or catalyst panel. Without this
    # call `panels` would be accepted, drawn nowhere and reported nowhere, and the
    # catalysts would have to be smuggled through `rows` instead, which then needs
    # `template: "off"` on every one to stop the packer templating a set of
    # unrelated catalysts on each other.
    #
    # `of` STILL DEFAULTS TO "grid", which this path does not have, so a panel
    # written without one raises rather than landing somewhere arbitrary. That is
    # the right answer here: on a figure with no product grid there is no default
    # worth guessing, and the message already names the blocks this figure does
    # have.
    if recipe.get("panels"):
        _PN.place(recipe["panels"], recipe, cfg, alloc, parts.append,
                  notes if notes is not None else [],
                  gid=97500, tid=tid + 400, base=base + 40000)

    if blk is not None:
        blk.emit(pads["table"][0]
                 + (W - pads["table"][0] - pads["table"][1] - blk.w) / 2.0,
                 top["table"], put)
    _emit_footnotes(recipe, cfg, alloc, put, tid)
    # AFTER the footnotes, and the order is load-bearing rather than tidy: the
    # title folds to the FIGURE's width, so everything that can set that width
    # has to be on the ledger first. Footnotes cannot widen it (they fold to the
    # grid), but reading the box after them is what makes that a fact rather
    # than an assumption.
    # BEFORE the title and AFTER everything else: a box asked for "figure" should
    # enclose the drawing, and the heading belongs outside it -- which is where
    # published figures put it.
    _emit_boxes(recipe, cfg, alloc, put, behind)
    tid = _emit_title(recipe, cfg, alloc, put, tid)
    return _finish(blk.header_prefix if blk is not None else header_prefix,
                   parts, alloc, page_mode, W, cfg["margin"], notes,
                   behind=behind, shift=shift, span=cfg.get("span"))


def render_grid(recipe: dict, header_prefix: str, cfg: dict | None = None,
                page_mode: str = "letter", alloc: Allocator | None = None,
                policy: Policy | None = None,
                notes: list | None = None,
                shift: list | None = None) -> str:
    """`notes` is an out-parameter: pass a list and the packer appends anything the
    caller has to be told (currently the page-overflow note). Same shape as `alloc`,
    which the tool layer already reads back after the call -- and the same reason,
    that a figure worth opening should not be turned into an exception."""
    cfg = _cfg(cfg)
    # PER FIGURE, not per measure(): the sectioned packer calls measure() once for
    # every section, so clearing it in there would keep only the last section's
    # findings. (_WIDEST_LEGEND does clear per call and so loses the earlier
    # sections; it is a different note with a different trigger.)
    _LEGEND_WIDE.clear()
    # Measured once, here, and handed to whoever emits it. The scheme is the only
    # block whose width is otherwise known only after it is drawn.
    hm = measure_header(recipe, cfg)
    has_grid = bool(recipe.get("products") or recipe.get("rows")
                    or recipe.get("sections"))
    # A SCHEME counts as a block. Testing the grid and the table only would refuse
    # a figure that is one reaction and nothing else while allowing a figure that
    # is one table.
    if not has_grid and hm is None and not recipe.get("table"):
        raise ValueError("recipe needs 'products', 'rows', 'sections', 'header' "
                         "or 'table'")
    alloc = alloc if alloc is not None else Allocator(policy)

    if not has_grid:
        return _render_no_grid(recipe, header_prefix, cfg, page_mode, alloc,
                               hm, notes, shift=shift)

    def _usable():
        """Usable page width, or None. 'fit' / 'auto' size the page to the
        content, so they legitimately have no fixed width -- usable_width
        RETURNS None rather than raising, which is easy to walk into.

        `span` ANSWERS THE SAME QUESTION FROM A DIFFERENT AUTHORITY. Every number
        below is asking one thing -- how wide may this be? -- and without `span`
        the only voice that can answer is the page. A grid that is one PART of a
        larger figure has no page: it sits next to a sibling 548 pt wide, and that is the
        width it has to live within. So a caller who knows it says so and the page
        is not consulted. This is footnote_width's move one level up: not a second
        definition of the width, the same definition pointed at the right figure.

        Deliberately NOT clamped to the page. A composed part is routinely wider
        than any single column -- the parent scales or splits the whole figure --
        and quietly shrinking it to letter's 552 would make the caller's number a
        lie in exactly the case it was passed for."""
        if cfg.get("span"):
            return float(cfg["span"])
        try:
            u = C.usable_width(page_mode, cfg["margin"])
        except Exception:
            return None
        return u if u else None

    # Sections are a different packing model, so they get their own packer and
    # leave this one alone -- six byte-exact goldens ride on the path below.
    if recipe.get("sections"):
        return _render_sectioned(recipe, header_prefix, cfg, page_mode, alloc,
                                 _usable(), hm, notes, shift=shift)

    flat, row_sizes = _flatten(recipe)
    recipe = dict(recipe, products=flat)
    meas, cell_w = measure(recipe, cfg)

    if row_sizes:
        # Explicit rows are taken literally; none of the width arithmetic below
        # applies, and running it on a page mode with no fixed width crashes.
        cols = max(row_sizes)
    else:
        cols = recipe.get("cols") or cfg.get("cols")
        if not cols:
            usable = _usable() or 540.0
            fits = max(1, int(usable // cell_w))
            # Packing as many columns as fit leaves ugly remainders -- 8 products
            # in a 7-wide page becomes 7 + 1. Take the fewest ROWS that hold
            # them, then spread the products over those rows: 8 -> 4 x 2.
            rows = max(1, -(-len(meas) // fits))
            cols = -(-len(meas) // rows)
        cols = max(1, min(int(cols), len(meas)))

    # A scope table spans the column. Once the column count is fixed, widen the
    # cells to fill the usable width so the products spread out instead of
    # huddling in the top-left. Capped so a 2-column grid of small rings does not
    # end up with absurd gutters.
    if cfg.get("fill_width", True) and cols > 1:
        usable = _usable()
        if usable:
            _nat, _want = cell_w, max(cell_w, usable / cols)
            cell_w = min(_want, _nat * cfg["max_spread"])
            # A SHORT REACH IS REPORTED, NOT SWALLOWED -- the rule the value
            # table's width_mode "match" already follows 200 lines up, and for
            # its reason: the cap is the right LOOK, but when it wins the figure
            # comes out narrower than was asked for and the only way to notice is
            # to measure the file.
            #
            # ONLY UNDER `span`. Under a page mode the cap winning means "the page
            # is wider than this grid wants to be", which is nobody's request. A
            # span is a request, so missing it is news.
            #
            # max_spread is not a tool argument, so the remedy named is the one
            # the caller actually has: the widest this grid can reach is
            # cols x natural x max_spread, and only `cols` is theirs to move.
            if notes is not None and cfg.get("span") and cell_w < _want - 0.01:
                notes.append(
                    "span %.0f pt was not reached: max_spread caps a cell at "
                    "%.2f x its natural %.0f pt, so %d columns stop at %.0f pt "
                    "and the grid sits to the left of the rest. More columns "
                    "reach further -- a cell's own width is set by the widest "
                    "product in it, not by the span."
                    % (usable, cfg["max_spread"], _nat, cols,
                       _nat * cfg["max_spread"] * cols))

    parts, behind = [], []   # decorations first, see _finish
    tid, base = 90000, 1000
    order = _resolve_order(cfg, {"scheme": hm is not None, "grid": True,
                                 "table": bool(recipe.get("table"))})
    # enclose is only meaningful when the scheme is ABOVE the table: it means "the
    # table's top rule is drawn above the scheme instead". With the table first
    # there is nothing to raise it over, and suppressing the table's own rule
    # would leave the figure with no top edge at all.
    enclose = (cfg.get("frame") == "enclose" and "scheme" in order
               and "table" in order
               and order.index("scheme") < order.index("table"))
    blk = (_TableBlock(recipe, cfg, header_prefix, enclose=enclose)
           if recipe.get("table") else None)
    if blk is not None:
        header_prefix = blk.header_prefix

    # A short final row reads as an accident when it is left-hanging under a full
    # grid. Journals centre it, which makes any remainder look deliberate. Only
    # worth doing when the row is clearly short -- nudging a 3-of-4 row by half a
    # cell just breaks the column alignment for no gain.
    if row_sizes:
        n_rows = len(row_sizes)
        n_last = row_sizes[-1]
    else:
        n_rows = -(-len(meas) // cols)
        n_last = len(meas) - (n_rows - 1) * cols
    mode = cfg.get("last_row", "auto")
    centre_short = (mode == "center" or
                    (mode == "auto" and n_last * 2 <= cols))
    last_dx = (cols - n_last) * cell_w / 2.0 if centre_short else 0.0

    if row_sizes:
        slots = []
        for r, k in enumerate(row_sizes):
            pad = (cols - k) * cell_w / 2.0 if centre_short else 0.0
            slots += [(r, c, pad) for c in range(k)]
    else:
        slots = [(i // cols, i % cols,
                  last_dx if i // cols == n_rows - 1 else 0.0)
                 for i in range(len(meas))]

    grid_w = cols * cell_w

    # Row heights are computed PER ROW. Sizing the structure band from the
    # tallest structure in the WHOLE table pushes every compound number down to
    # clear it, so a short structure ends up with its label stranded far below.
    members = [[] for _ in range(n_rows)]
    for i, (r, _c, _p) in enumerate(slots):
        members[r].append(i)
    band = [max(meas[i]["sh"] for i in mem) for mem in members]
    stack = [max(meas[i]["stack_h"] + meas[i]["legend_h"] for i in mem)
             for mem in members]
    grid_h = sum(band[r] + stack[r] + 2 * cfg["cell_pad_y"]
                 for r in range(n_rows))

    # ---- W, the figure width: resolved HERE, then handed down --------------
    # One rule, no branches: the widest natural wins and everything else centres
    # in it. grid_w is already the grid's page-aware fill target, so a
    # figure with neither a scheme nor a table resolves W to exactly the width it
    # drew before -- which is what keeps the panel-only goldens byte-identical.
    W = max(hm["natural"] if hm is not None else 0.0, grid_w,
            blk.w if blk is not None else 0.0)
    if blk is not None:
        blk.fit(W, notes)       # width_mode "match" only; "natural" is a no-op

    # ---- y: one accumulator over `order`, no per-arrangement special cases ---
    heights = {"scheme": (hm["up"] + hm["down"]) if hm is not None else 0.0,
               "grid": grid_h, "table": blk.h if blk is not None else 0.0}
    top, y = {}, 0.0
    for i, name in enumerate(order):
        if i:
            y += _block_gap(order[i - 1], name, cfg)
        top[name] = y
        y += heights[name]
    # Each block becomes a REGION, so a panel can be placed relative to one. The
    # extents are the ones this loop just computed -- not measured again, which
    # would be two accounts of one number.
    for name in order:
        if heights[name] > 0.0:
            alloc.region("block:%s" % name, Rect(0.0, top[name], W, heights[name]))

    def put(xml, rect, owner, kind):
        res = alloc.reserve(rect, owner, kind)
        alloc.check_emit(res.rid, owner)
        parts.append(xml)
        return res

    # Emitted scheme -> cells -> table -> footnotes regardless of the VERTICAL
    # order, which only decides the y values above. Two reasons to keep the
    # sequence fixed: the default order then produces byte-identical output, and
    # stamp_z hands out Z in document order, so a table's shaded bands keep
    # landing behind the type instead of depending on where the table sits.
    if hm is not None:
        hh, tid, base = _emit_header(recipe, cfg, alloc, put, W, tid, base, hm,
                                     y0=top["scheme"], behind=behind)
        # The rule under the scheme separates it from the GRID. With the table
        # next instead, the table's own top rule is already that separator and a
        # second rule 26 pt above it would be a line with nothing on either side.
        nxt = order[order.index("scheme") + 1] if order[-1] != "scheme" else None
        if nxt == "grid":
            ry = top["scheme"] + hh + cfg["header_gap"] * 0.5
            try:
                put(C._line(96000, 0.0, ry, W, ry, 0.8),
                    Rect(0.0, ry - 0.4, W, 0.8), "hdr:rule", "decoration")
            except Exception:
                pass
        if enclose:
            # ABOVE the scheme, and as heavy as the table's own top rule: this
            # IS the table's top rule, moved up so the scheme sits inside it.
            fw = TL._cfg(blk.cfg, recipe["table"].get("style"))["rule_top"] or 1.2
            fy = top["scheme"] - cfg["header_gap"] * 0.5
            put(C._line(96001, 0.0, fy, W, fy, fw),
                Rect(0.0, fy - fw / 2, W, fw), "frame:top", "decoration")

    row_top, y_acc = [], top["grid"]
    for r in range(n_rows):
        row_top.append(y_acc + cfg["cell_pad_y"])
        y_acc += band[r] + stack[r] + 2 * cfg["cell_pad_y"]

    # The grid is at grid_w, which may be short of W when the scheme or the table
    # is the widest block, or when max_spread capped the fill. The residue is
    # whitespace either side of it, never wider cells.
    gdx = (W - grid_w) / 2.0
    lalign = cfg.get("label_align", "center")
    cgid = _CELL_GID
    for i, m in enumerate(meas):
        row, col, pad = slots[i]
        _n0 = len(parts)
        tid, base, _cbox = _emit_cell(m, gdx + col * cell_w + pad, cell_w,
                                      row_top[row], band[row], cfg, put, tid,
                                      base, lalign, behind=behind)
        cgid = _group_cell(parts, _n0, _cbox, cfg, cgid)

    # AFTER every cell: a connector asks the ledger where they landed.
    _aid, tid = _emit_connectors(recipe, cfg, alloc, put, 90600, tid)

    # Panels, placed relative to a block region or to an earlier panel.
    if recipe.get("panels"):
        import panels as _PN                 # deferred: _PN imports this module
        _PN.place(recipe["panels"], recipe, cfg, alloc, parts.append,
                  notes if notes is not None else [],
                  gid=97500, tid=tid + 400, base=base + 40000)

    if blk is not None:
        blk.emit((W - blk.w) / 2.0, top["table"], put)

    tid = _emit_footnotes(recipe, cfg, alloc, put, tid)
    # AFTER the footnotes, and the order is load-bearing rather than tidy: the
    # title folds to the FIGURE's width, so everything that can set that width
    # has to be on the ledger first. Footnotes cannot widen it (they fold to the
    # grid), but reading the box after them is what makes that a fact rather
    # than an assumption.
    # BEFORE the title and AFTER everything else: a box asked for "figure" should
    # enclose the drawing, and the heading belongs outside it -- which is where
    # published figures put it.
    _emit_boxes(recipe, cfg, alloc, put, behind)
    tid = _emit_title(recipe, cfg, alloc, put, tid)
    return _finish(header_prefix, parts, alloc, page_mode, W, cfg["margin"],
                   notes, behind=behind, shift=shift, span=cfg.get("span"))
