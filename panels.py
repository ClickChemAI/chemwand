"""panels.py — a lettered, titled, framed box holding a small reaction.

Published figures put one or more of these beside or under the main content:
lettered A/B/C sub-schemes, framed bands, a panel C tucked into the bottom-right
of panel B. What they hold is a transformation -- nearly always a reaction in at
least one sub-region -- plus a letter and a title.

WHY THIS IS NOT A NEW ENGINE, which the layout design rules out. Two decisions,
and both come from the evidence rather than from taste:

  A PANEL IS COMPOSED, NOT LAID OUT. `C.group()` wraps already-built objects and
  gives the result its own BoundingBox, and C.place() then treats it as one
  ordinary cell -- nested.py does the same for a vertical stack. Arrows inside
  such a group are supported: place() moves the absolute 3D endpoints an <arrow>
  carries, and ChemDraw (measured) keeps such a group, selects it as one object
  and moves the arrow with it.

  A PANEL IS PLACED RELATIVE TO A BLOCK, not on a canvas of its own. Every
  two-dimensional need seen in published figures is of the form "panel C at B's
  bottom-right" or "the ligand panel to the RIGHT of the data table" -- none
  needs absolute coordinates. So a panel names what it sits beside, and
  Allocator.reserve_near does the arithmetic it was written for.

THE ROW INSIDE IS layout_grid's, borrowed rather than rewritten. A panel's content
is a small reaction row -- structures, "+", one arrow with captions, text cells
-- which is exactly what `measure_header` / `_emit_header` already build for the
scheme above a scope table. That pair is called here with a THROWAWAY allocator and
a collecting `put`, so the row arithmetic has one home. Writing a third row layout
(after layout_linear and _emit_header) is the shape this codebase keeps paying for.

WHAT IS DELIBERATELY NOT HERE (waiting for a second use): a panel holding a whole
product grid or data table. That needs render_grid to be embeddable as a body
fragment, which is the "shared canvas" refactor the layout design excludes.
"""

from __future__ import annotations

import cdxml_core as C
from allocator import Allocator, Rect

# Panel furniture. Anchored to numbers already in the figure rather than chosen:
PANEL_PAD = 12.0          # frame -> content. == legend_gap 6 doubled; see below.
PANEL_TITLE_GAP = 6.0     # title baseline block -> content
PANEL_LETTER_SIZE = 12    # one step above num_size 11: the letter labels a BLOCK
PANEL_TITLE_SIZE = 10     # == label_size, which is what a titled thing uses


def compose(panel: dict, recipe: dict, cfg: dict, *, gid: int, tid: int,
            base: int):
    """(xml, bbox) for one panel, built at the origin.

    `panel` is {"name", "items", "title"?, "frame"?}; `recipe` is the figure's own,
    because the row builder reads `structures` out of it.

    Built at the origin and placed later, which is what makes it a cell: the panel
    cannot know where it will go, and nothing in here asks.
    """
    import layout_grid as GR              # deferred: GR imports nothing from here

    # The row, through the grid's own header builder. A throwaway allocator: this
    # measures and emits into panel-local coordinates, and the FIGURE's ledger
    # learns about the panel as one region plus its frame (see place()). Reserving
    # the row's insides twice -- once here, once shifted -- would be two accounts of
    # one row.
    local = Allocator()
    parts: list = []

    def put(xml, rect, owner, kind):
        parts.append(xml)
        return local.reserve(rect, "panel:%s.%s" % (panel["name"], owner), kind)

    sub = {"structures": recipe["structures"], "header": {"items": panel["items"]}}
    m = GR.measure_header(sub, cfg)
    if m is None:
        raise ValueError("panel %r has no items to draw" % panel.get("name"))
    _h, tid, base = GR._emit_header(sub, cfg, local, put, m["natural"], tid, base,
                                    m, y0=0.0)
    inner = local.bbox()

    # letter and title ABOVE the row, in the reading order a paper uses: "(A)" then
    # the title then the chemistry.
    head: list = []
    y = inner.y
    if panel.get("title") or panel.get("name"):
        lines = []
        if panel.get("name"):
            lines.append("#(%s)#" % panel["name"])
        if panel.get("title"):
            lines.append(str(panel["title"]))
        from layout_linear import block_box, ASCENT_EM
        w, h = block_box(lines, PANEL_TITLE_SIZE)
        y = inner.y - PANEL_TITLE_GAP - h
        head.append(C.styled_text(tid, inner.x, y + ASCENT_EM * PANEL_TITLE_SIZE,
                                  lines, size=PANEL_TITLE_SIZE, align="left"))
        tid += 1
        local.reserve(Rect(inner.x, y, w, h),
                      "panel:%s.title" % panel["name"], "caption")
        inner = local.bbox()

    box = [inner.x - PANEL_PAD, inner.y - PANEL_PAD,
           inner.right + PANEL_PAD, inner.bottom + PANEL_PAD]
    # THE COMPOSITION CLAIM, checked where the data is. Everything above was
    # reserved in `local` under "panel:<name>.<something>", so declaring the frame as
    # the region and asking check_region is a real question with a real answer: is
    # anything composed here sticking out of the box about to be drawn round it? The
    # figure's own ledger cannot answer it -- the insides never go in there (see
    # place: the panel arrives as ONE reservation), which is the point of composing.
    local.region("panel:%s" % panel["name"],
                 Rect(box[0], box[1], box[2] - box[0], box[3] - box[1]))
    local.check_region("panel:%s" % panel["name"])
    escapes = [f for f in local.findings if f.code == "region-escape"]
    frame = panel.get("frame")
    if frame is not None:
        if frame not in C.RECTANGLE_TYPES:
            raise ValueError(
                "panel %r: frame %r has no precedent in any ChemDraw file. Known: %s"
                % (panel.get("name"), frame, ", ".join(C.RECTANGLE_TYPES)))
        head.insert(0, C.rect(gid, *box, kind=frame))

    xml = C.group(gid + 1, head + parts,
                  bbox=[round(v, 2) for v in box])
    return xml, [round(v, 2) for v in box], tid, base, escapes


def measure(panel: dict, recipe: dict, cfg: dict) -> tuple:
    """(w, h) of one panel, composed and thrown away.

    A packer has to know how wide a SIDE panel is before it resolves the figure
    width -- a panel measured after the fact can only be hung off the outside,
    which is what happens: on one optimisation figure the scheme and the
    table both ended at 545 pt and the ligand box ran from 557 to 673, so the
    figure was 673 pt wide with its contents in the left 545 and nothing sharing
    a right-hand edge.

    Nothing predicts that width from the spec: the row goes through the grid's own
    header builder, so the only thing that knows is compose(). Composing twice is
    the price. The measuring pass emits into compose's own throwaway allocator and
    its xml, its escapes and its id counters are all discarded here, so it cannot
    reach the document -- which is also why it can reuse the same ids.
    """
    _xml, box, _tid, _base, _escapes = compose(panel, recipe, cfg,
                                               gid=0, tid=0, base=0)
    return box[2] - box[0], box[3] - box[1]


# --- placement ---------------------------------------------------------------
_SIDES = ("right", "left", "above", "below")
_ALIGNS = ("start", "center", "end")


def place(panels, recipe, cfg, alloc: Allocator, emit, notes: list, *,
          gid: int, tid: int, base: int):
    """Compose and place every panel. Returns (gid, tid, base).

    Placement is `reserve_near` against the region the panel names -- a block
    ("grid", "table", "scheme") or an earlier panel. That function was written
    for exactly this.

    ORDER MATTERS AND IS THE CALLER'S: a panel may sit beside a panel, so they are
    placed in the order given and a forward reference is an error rather than a
    silent 0,0.
    """
    for i, p in enumerate(panels):
        of = p.get("of", "grid")
        anchor_name = of if of.startswith("panel:") else "block:%s" % of
        if of in {q.get("name") for q in panels[:i]}:
            anchor_name = "panel:%s" % of
        # BEFORE THE LOOKUP, because a failed lookup is not the evidence: the
        # sectioned packer publishes block:grid and block:scheme, because `boxes`
        # needs somewhere to hang a frame on a sectioned figure. The anchor
        # therefore RESOLVES, and without this check a panel on a sectioned table
        # would place itself without a word -- which is the thing this refusal
        # exists to prevent.
        #
        # THE REASON is not "there is nothing to sit beside"; it is that this
        # pair has never been laid out or checked. One render says it comes out
        # plausibly -- panel to the right of the whole sectioned grid, level with
        # the first row -- and one render is not a specification: which section a
        # panel should anchor to, and what happens when there are three of them,
        # are open questions. Turning the combination on is its own piece of work
        # with its own tests.
        if of == "grid" and (recipe or {}).get("sections"):
            raise ValueError(
                "panel %r: this figure is built from `sections`, and a panel "
                "beside a SECTIONED table is not built -- the two features have "
                "never been laid out together, so which section the panel lines "
                "up with is undefined. Either drop `sections` (a flat `rows` list "
                "takes panels), or place the panel with an explicit 'at'."
                % (p.get("name"),))
        try:
            anchor = alloc.find(anchor_name)
        except Exception:
            raise ValueError(
                "panel %r: 'of' names %r, which is not a block or an earlier panel "
                "in this figure. Blocks are scheme / grid / table (whichever this "
                "figure has); a panel may also name one placed BEFORE it."
                % (p.get("name"), of))
        side = p.get("side", "right")
        if side not in _SIDES:
            raise ValueError("panel %r: side is one of %s, got %r"
                             % (p.get("name"), ", ".join(_SIDES), side))
        align = p.get("align", "start")
        if align not in _ALIGNS:
            raise ValueError("panel %r: align is one of %s, got %r"
                             % (p.get("name"), ", ".join(_ALIGNS), align))

        xml, box, tid, base, escapes = compose(p, recipe, cfg, gid=gid, tid=tid,
                                              base=base)
        for f in escapes:
            notes.append("panel %r: %s" % (p["name"], f.message))
        w, h = box[2] - box[0], box[3] - box[1]
        gap = float(p.get("gap", cfg["table_gap"]))
        # ONE reservation for the whole panel, kind "panel" -- NOT "region".
        # A region is tolerant of everything, because it covers content named
        # elsewhere; a panel is INK, and a panel sitting on the grid or on another
        # panel is the one thing this layer exists to prevent. ink_tol["panel"] is
        # 0.0 for the same reason a caption's is.
        res = alloc.reserve_near(anchor, side, w, h,
                                 "panel:%s" % p["name"], kind="panel",
                                 gap=gap, align=align)
        # `emit` only APPENDS: reserve_near already made the one reservation,
        # and reserving the same rect again would be two accounts of one panel.
        emit(C.place(xml, base, res.rect.x, res.rect.y))
        alloc.check_emit(res.rid, "panel:%s" % p["name"])
        base += 1000
        gid += 2
    return gid, tid, base
