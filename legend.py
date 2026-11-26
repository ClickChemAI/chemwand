"""legend.py — the table of values that sits under a Markush structure.

A generic structure is only half a figure. The other half is the list that says
what the placeholders stand for, and in published figures that list is written
five different ways:

    23, R = OMe, 72%, >19:1 dr        number first, variable repeated per line
    R = H, 2a, 73% yield              variable first, number in the middle
    3aa (X = Me, 94%, E/Z = 3.6/1)    number first, the rest inside brackets
    R2 = Me, 3aa: 74%                 variable first, number last, colon
            R1                        variable name lifted into a heading
    4a =    H, 97%, 3 h               so the line itself is only values

WHY A LIST OF CELLS AND NOT A RECORD

The obvious model is {id, value, yield, dr}. It breaks on the fourth line above.
Compare it with the fifth: both use "=", and in one the "=" joins a VARIABLE to
its value while in the other it joins a COMPOUND NUMBER to its value. One symbol,
two unrelated jobs. A typed record has to grow a field or a flag for every new
arrangement, and the arrangements belong to each figure's author, not here.

So a row is an ORDERED LIST OF CELLS and the order is the caller's. This module
never asks what a cell means. It asks how wide it is.

WHAT IT ACTUALLY ADDS

Column alignment across rows, which is the one thing the pipeline did not have.
layout_grid stacks a number, a yield and any sublabels under each structure and
centres each stack on its own cell -- nothing anywhere lines a column up across
rows. Everything else a legend needs was already here: ~sub~ ^sup^ *italic*
#bold# all work, and block_box already measures them.

Deliberately NOT merged into layout_grid. A scope-table cell holds a molecule and
a legend cell holds a word; the only thing they share is "a column is as wide as
its widest member", and that is one line, not an engine. If the scope table ever
wants aligned columns it can call this.

RECIPE

    {
      "rows": [
        ["#23#,", "R = OMe,", "72%,", ">19:1 *dr*"],
        ["#24#,", "R = Me,",  "68%,", ">19:1 *dr*"],
        "from Benzocaine",              # a bare string is a NOTE row: no columns,
                                        # flush left, for the asides that sit
                                        # between entries in a real legend
      ],
      "header": [None, "R^1^", None, None],   # optional, per column; None = none
      "align":  ["left", "left", "right", "left"],   # optional, per column
      "cols": 2,                        # optional: two blocks side by side
    }

Returns (object_string, bbox) -- the same pair `S` holds and `nested.compose_stack`
returns, so a legend can be registered as a structure and placed by any packer, or
dropped straight under a molecule with C.place().

It composes at the origin and does NOT reserve anything: a composite that is about
to be translated cannot reserve its final rect, and the caller who places it knows
where it lands. Call measure() to get the box, reserve that, then place.
"""

from __future__ import annotations

import math

import cdxml_core as C
# A leaf module: nothing in the layout engines imports this one, so the plain
# import is safe. block_box is THE marker-aware measurer for this pipeline
# (~sub~ ^sup^ *it* #bold#) and a legend must measure with the same ruler the
# captions beside it use, or the two sizes of "9pt" diverge.
from layout_linear import ASCENT_EM, block_box

LEGEND_DEFAULTS = {
    "size": 9,             # == layout_linear's cap_size: a legend IS a caption
    "header_size": None,   # None -> size
    "note_size": None,     # None -> size
    "header_bold": False,  # published headings are plain; the NUMBERS are bold,
                           # and those are written #23# in the cell itself
    "row_gap": 3.0,        # == sublabel_gap: line-to-line inside one label stack
    "col_gap": 10.0,       # between columns of one block
    "block_gap": 30.0,     # between the blocks of a two-column legend
    "header_gap": 4.0,     # extra space under the heading row
    "note_indent": 0.0,
    "align": "left",       # default per-column alignment
    "fill": "column",      # "column": down the first block then the second
                           # "row":    alternate, so reading order goes across
}


def _cfg(user: dict | None) -> dict:
    out = dict(LEGEND_DEFAULTS)
    out.update(user or {})
    return out


def _normalize(recipe: dict) -> list:
    """Rows as ("cells", [str, ...]) or ("note", str)."""
    out = []
    for r in recipe.get("rows") or []:
        if isinstance(r, str):
            out.append(("note", r))
        else:
            # A dict cell survives as a dict. str() on it would give the
            # repr -- "{'swatch': 'red', ...}" -- which measures 174 pt wide and
            # draws itself out on the page.
            out.append(("cells", ["" if c is None
                                  else (c if isinstance(c, dict) else str(c))
                                  for c in r]))
    return out


def _line_height(texts, size: float) -> float:
    """The step every line of this legend takes: the TALLEST line in it.

    One height for all of them, for the reason layout_grid sizes all its cells
    together -- a list only reads as a list if the lines step evenly, and a row
    measured on its own content steps short whenever it happens to have no
    descender. "24," is taller than "24" by the depth of the comma, and half a
    legend's rows end in one.

    Taking the max rather than a nominal line height is not cosmetic: the box
    returned to the caller is what the Allocator reserves and nobody re-measures
    the contents, so a band shorter than its own ink is an overflow that never
    gets reported. With "X" alone setting this, the commas hang out the bottom.
    """
    return max([block_box(["X"], size)[1]]
               + [block_box([_cell_text(t)], size)[1] for t in texts
                  if _cell_text(t)])


# --- a colour-swatch cell -----------------------------------------------------
# {"swatch": "red", "text": "previous work"} -- a filled dot, then the words.
# The colour key legends carry down one side of a late-stage-functionalisation
# figure, pairing each tint on the structures with what it means.
#
# A DICT RATHER THAN A SEVENTH INLINE MARKER, and the reason is a kind
# difference rather than a count. ~x~ ^x^ *x* #x# _x_ {red} all answer "what does
# this text look like"; a swatch answers "put an object here". Folding it into
# that grammar would invite `@red@` in a sublabel or an arrow caption, where a
# dot means nothing and would be drawn anyway. A dict is visibly a legend cell's
# own shape and nobody carries it somewhere it does not belong.
SWATCH_R = 3.0          # dot radius
SWATCH_GAP = 5.0        # dot to text


def _cell_text(cell) -> str:
    """The words of a legend cell, whichever shape it was written in."""
    return str(cell.get("text", "")) if isinstance(cell, dict) else str(cell or "")


def _cell_width(cell, size, cfg) -> float:
    w = block_box([_cell_text(cell)], size)[0] if _cell_text(cell) else 0.0
    if isinstance(cell, dict) and cell.get("swatch") is not None:
        w += 2 * SWATCH_R + SWATCH_GAP
    return w


def _widths(rows: list, header, cfg: dict) -> list:
    """Column widths, over EVERY row in every block.

    Shared across blocks on purpose. Two blocks of a two-column legend are two
    halves of one list, and the eye reads down both; sizing them apart would put
    the same column at two different widths on one line.
    """
    n = max([len(c) for k, c in rows if k == "cells"]
            + [len(header or [])] + [0])
    w = [0.0] * n
    for kind, cells in rows:
        if kind != "cells":
            continue
        for j, cell in enumerate(cells):
            if cell:
                w[j] = max(w[j], _cell_width(cell, cfg["size"], cfg))
    for j, head in enumerate(header or []):
        if head:
            w[j] = max(w[j], block_box([head], cfg["header_size"] or cfg["size"])[0])
    return w


def _blocks(rows: list, cols: int, fill: str) -> list:
    if cols <= 1 or not rows:
        return [rows]
    if fill == "row":
        return [rows[k::cols] for k in range(cols)]
    per = math.ceil(len(rows) / cols)
    return [rows[k * per:(k + 1) * per] for k in range(cols)]


def _anchor(x0: float, w: float, align: str):
    """(anchor x, justification) for a cell in a column of width w."""
    if align == "right":
        return x0 + w, "right"
    if align == "center":
        return x0 + w / 2.0, "center"
    return x0, "left"


def _aligns(recipe: dict, ncols: int, cfg: dict) -> list:
    spec = recipe.get("align", cfg["align"])
    if isinstance(spec, str):
        return [spec] * ncols
    out = list(spec) + [cfg["align"]] * ncols
    return out[:ncols]


def measure(recipe: dict, cfg: dict | None = None) -> dict:
    """Everything the layout needs, without emitting anything.

    Returns {"w", "h", "colw", "colx", "block_w", "bands", "blocks", "head_h"}.
    Separate from compose_legend() so a caller can reserve the right rect before
    the legend exists -- the Allocator wants the box first, and a composite that
    is about to be translated cannot reserve its own.
    """
    cfg = _cfg(cfg)
    rows = _normalize(recipe)
    header = recipe.get("header")
    colw = _widths(rows, header, cfg)

    colx, acc = [], 0.0
    for w in colw:
        colx.append(acc)
        acc += w + cfg["col_gap"]
    block_w = max(0.0, acc - cfg["col_gap"])

    # Read off the RECIPE, with cfg as the fallback, the same way "align" is:
    # how many blocks and which way they fill are properties of this legend's
    # shape, not house style. (cfg["fill"] alone would silently ignore the recipe.)
    blocks = _blocks(rows, int(recipe.get("cols") or 1),
                     recipe.get("fill", cfg["fill"]))

    # One band per row INDEX, tall enough for that row in every block, so the two
    # halves of a two-column legend step down together instead of drifting apart.
    depth = max((len(b) for b in blocks), default=0)
    body_h = _line_height([c for k, cells in rows if k == "cells" for c in cells],
                          cfg["size"])
    note_h = _line_height([p for k, p in rows if k == "note"],
                          cfg["note_size"] or cfg["size"])
    bands = []
    for i in range(depth):
        kinds = {b[i][0] for b in blocks if i < len(b)}
        bands.append(note_h if kinds == {"note"} else body_h)

    head_h = 0.0
    if header and any(h for h in header):
        head_h = (_line_height(header, cfg["header_size"] or cfg["size"])
                  + cfg["header_gap"])

    ncols = len(blocks)
    total_w = ncols * block_w + (ncols - 1) * cfg["block_gap"] if block_w else 0.0
    total_h = head_h + sum(bands) + max(0, len(bands) - 1) * cfg["row_gap"]
    return {"w": total_w, "h": total_h, "colw": colw, "colx": colx,
            "block_w": block_w, "bands": bands, "blocks": blocks,
            "head_h": head_h, "cfg": cfg, "header": header,
            "aligns": _aligns(recipe, len(colw), cfg)}


def place_legend(xml: str, base: int, x: float, y: float) -> str:
    """Move a composed legend to (x, y). USE THIS, NOT C.place().

    compose_legend builds at the origin, so the shift is exactly (x, y) -- but
    C.place() translates with C.translate(), which moves `p=` and BoundingBox and
    NOT the Center3D / MajorAxisEnd3D a disc carries. Placed that way, a
    swatch legend comes out with its words in the cell and its coloured dots
    stranded at the top-left corner of the page.

    The general rule -- "a decoration reaches every engine for free if and only
    if it is expressible inside the fragment", see also translate_objects' own
    docstring -- is easy to miss here, because a legend LOOKS like text and
    text moves fine. It only bites when a legend has a swatch in it.

    One function, used by every caller, so the next one cannot reach for place().
    """
    return C.translate_objects(C.renumber(xml, base), x, y)


def compose_legend(recipe: dict, cfg: dict | None = None, *, gid: int = 1):
    """Build the legend as one placeable, groupable object at the origin."""
    m = measure(recipe, cfg)
    cfg = m["cfg"]
    header, colx, colw, aligns = m["header"], m["colx"], m["colw"], m["aligns"]
    hsize = cfg["header_size"] or cfg["size"]
    nsize = cfg["note_size"] or cfg["size"]

    children, tid = [], gid + 1

    def put(text, x0, width, align, y, size, bold=False):
        nonlocal tid
        ax, just = _anchor(x0, width, align)
        children.append(C.styled_text(tid, ax, y + ASCENT_EM * size, [text],
                                      size=size, bold=bold, align=just))
        tid += 1

    for k, block in enumerate(m["blocks"]):
        bx = k * (m["block_w"] + cfg["block_gap"])
        y = 0.0
        if m["head_h"]:
            for j, head in enumerate(header or []):
                if head and j < len(colw):
                    put(head, bx + colx[j], colw[j], aligns[j], y, hsize,
                        cfg["header_bold"])
            y += m["head_h"]
        for i, band in enumerate(m["bands"]):
            if i < len(block):
                kind, payload = block[i]
                if kind == "note":
                    if payload:
                        put(payload, bx + cfg["note_indent"], m["block_w"],
                            "left", y, nsize)
                else:
                    for j, cell in enumerate(payload):
                        if not cell or j >= len(colw):
                            continue
                        if isinstance(cell, dict) and cell.get("swatch") is not None:
                            # dot first, then the words beside it -- the pairing
                            # a colour key reads left to right. LEFT-aligned
                            # whatever the column says: a dot centred in a wide
                            # column stops lining up with the dots above it,
                            # which is the one thing a key has to get right.
                            x = bx + colx[j]
                            children.append(C.disc(
                                tid, x + SWATCH_R, y + SWATCH_R * 2.0,
                                SWATCH_R, cell["swatch"]))
                            tid += 1
                            txt = _cell_text(cell)
                            if txt:
                                put(txt, x + 2 * SWATCH_R + SWATCH_GAP,
                                    colw[j] - 2 * SWATCH_R - SWATCH_GAP,
                                    "left", y, cfg["size"])
                        else:
                            put(_cell_text(cell), bx + colx[j], colw[j],
                                aligns[j], y, cfg["size"])
            y += band + (cfg["row_gap"] if i < len(m["bands"]) - 1 else 0.0)

    # Ceil to 2dp for the same reason nested.compose_stack does: C.group writes
    # the bbox with %.2f, so a raw float would let the number the caller reserves
    # and the number ChemDraw reads differ by up to half a hundredth -- half the
    # time in the direction that under-reserves. Round UP, always.
    w = math.ceil(m["w"] * 100.0) / 100.0
    h = math.ceil(m["h"] * 100.0) / 100.0
    bbox = [0.0, 0.0, w, h]
    if not children:
        return C.group(gid, [], bbox=bbox), bbox
    return C.group(gid, children, bbox=bbox), bbox
