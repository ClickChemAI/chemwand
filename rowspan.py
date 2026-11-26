"""rowspan.py -- one cell occupying ONE SLOT across SEVERAL ROWS.

Two shapes a scheme has to say directly:

    one thing facing several      23 fans out to 25/26/27/28   (fork)
    several things acting as one  A and B both feed C          (merge, cross-row)

`nested.py` says the in-row version of the second shape. This says the cross-row
version of both.

THE COMMITTED DESIGN: THE FAN IS HORIZONTAL
-------------------------------------------
Arrows leave a spanning cell's right edge at each covered row's own baseline.
Four rows, four ordinary horizontal arrows -- NOT four lines radiating from one
point. The reason that settles it is not aesthetic:

    a horizontal fan arrow IS an ordinary arrow cell

so `above` / `below` / `yield` / `reagent` / `reagentEq`, composite arrows and
`ink_extra` overhang reservation all work on it the day it is written. A
radiating arrow is a new object type and every one of those would have to be
re-derived along a diagonal. That is the difference between reusing the arrow
subsystem and rebuilding it.

It also makes diagonal caption placement and caption rotation MOOT rather
than unsolved, and removes connector-crossing-connector collisions by
construction.

THE MODEL
---------
A spanning cell occupies slot k in rows i..j. Every covered row splits into

    pre   cells at slots < k        post  cells at slots > k

    fork    pre is empty            (23 at slot 0, products to its right)
    merge   post is empty           (A and B to the left, C at the last slot)

and the two cases are then the SAME computation, which is why this is one
mechanism and not two:

    1. lay out each covered row's `pre`
    2. the span's left edge is the RIGHTMOST of those (plus the gap), so no
       row's pre is ever cut short
    3. every other row's shortfall is absorbed by lengthening the last arrow
       in its pre -- a row with an empty pre has nothing to lengthen and needs
       nothing, which is the fork falling out of the merge
    4. every covered row's `post` starts at the span's right edge plus the gap

WHY THE ARROW HAS TO STRETCH
----------------------------
`measure_arrow` derives a length from ONE thing -- the widest caption:

    length = max(widest + 2*arrow_pad, min_arrow)

In a merge, two rows' arrows must both reach the same spanning cell while the
cells in front of them have different widths. So the length is now the larger
of two relations, not one. This is the "compared to what" pattern again: a
quantity that was determined by one thing is now determined by two.

WHY row_align STOPS APPLYING TO A COVERED ROW
---------------------------------------------
`row_align: left` puts every row's true left edge on a common margin. A covered
row in a fork begins with an arrow that must start at the SPAN's right edge, not
at the margin. So a covered row's x is fixed by the span, and the alignment rule
only governs rows the span does not touch.

SCOPE OF THIS VERSION
---------------------
The span must sit at the same slot index in every row it covers, and a covered
row must be either all-pre or all-post (i.e. fork or merge, not a span with
content on both sides in a covered row). The owning row may have both. That
covers the fork and merge schemes examined; a covered row straddling the span
has not been seen in published schemes and is rejected loudly rather than
mis-drawn.
"""

from __future__ import annotations

import cdxml_core as C


def _cell_kind(cell) -> str:
    """What a cell IS, named the way the packer names it (layout_linear pass 1),
    so an error message and the code use one vocabulary."""
    if "s" in cell:
        # Named here too, because the merge-direction message below reports the
        # cell that SHOULD have been an arrow and "cell with keys ['s']" is a
        # worse way to say "structure". The raise above can never see this branch
        # -- it only fires when "s" is absent.
        return "structure cell"
    for k in ("text", "plus", "arrow", "stack", "ghost"):
        if k in cell:
            return "%s cell" % k
    return "cell with keys %s" % sorted(cell)


def plan_spans(rows_src):
    """Find `rowspan` cells and insert ghost cells so pass 1's x advances.

    Returns (rows, spans). `rows` is a NEW list when anything spans, and the
    SAME object when nothing does -- so a recipe without rowspan cannot take a
    different path, the same opt-in guarantee nested.py gives.
    """
    if not any(c.get("rowspan") for row in rows_src for c in row):
        return rows_src, []

    rows = [list(r) for r in rows_src]
    spans = []
    for ri, row in enumerate(rows):
        for k, cell in enumerate(row):
            n = cell.get("rowspan")
            if not n:
                continue
            if n < 2:
                raise ValueError("rowspan must be >= 2 (got %r)" % n)
            last = ri + n - 1
            if last >= len(rows):
                raise ValueError("rowspan %d at row %d runs past the last row"
                                 % (n, ri))
            # The root is fetched back later as a laid-out STRUCTURE -- align_x
            # and align_y both look for `it["type"] == "s"` with this id. A cell
            # that has no "s" never matches, and the failure arrives as a bare
            # StopIteration out of a generator: no message, no row number, and
            # nothing pointing at the cell that caused it. Say it here, where the
            # indices are still in hand. "stack" is not an exception to the rule:
            # nested.expand_stacks has already rewritten it into an "s" cell by
            # the time this runs.
            if "s" not in cell:
                raise ValueError(
                    'rows[%d][%d]: the root of a fork/merge must be a STRUCTURE '
                    'cell ("s", or a "stack", which becomes one) -- got a %s. '
                    'The spanning cell is what every covered row\'s arrow leaves '
                    'from or points at, so there has to be a drawn structure '
                    'there.' % (ri, k, _cell_kind(cell)))
            spans.append(dict(owner=ri, slot=k, rows=n, last=last,
                              id=cell.get("s")))
    # Where does the span attach in the rows it COVERS? Not at the same slot
    # index -- that model fails the moment the covered row has a different
    # number of cells than the owning one: the span can be the owner row's
    # cell 2 of 3 while the covered row's cell 2 of 5 is the middle of a
    # route.
    #
    # The side is what is fixed, not the index:
    #
    #     slot 0 in the owning row  -> a FORK  -> the covered rows START at it
    #     anywhere else             -> a MERGE -> the covered rows END at it
    #
    # A span at slot 0 cannot be a merge: there is no room to the left of it for
    # a covered row's cells to be.
    for si, sp in enumerate(spans):
        sp["fork"] = (sp["slot"] == 0)
        ghost = {"ghost": sp["slot"], "_span": si}
        for r in range(sp["owner"] + 1, sp["last"] + 1):
            # A MERGE'S COVERED ROW HAS TO POINT AT THE SPAN. Its cells are all
            # `pre`, so the LAST of them is what reaches the spanning cell -- and
            # if that is a structure rather than an arrow, nothing does.
            #
            # The row still draws, which is the whole problem. It comes out as
            # that structure FEEDING the span; ok stays true, the InChI round trip
            # passes, the geometry report finds nothing, and the only thing wrong
            # is that every arrow means the opposite of what was intended. Seen
            # once on a real endgame, where the span was a natural product and the
            # covered row its phosphate prodrug: the drawing said the prodrug made
            # the natural product.
            #
            # Rejected rather than guessed at, and the guess is named instead:
            # "ends with a structure" is exactly what a FORK's branch looks like.
            if not sp["fork"] and rows[r] and "arrow" not in rows[r][-1]:
                raise ValueError(
                    "rows[%d] is covered by the MERGE rooted at rows[%d][%d]: its "
                    "cells all sit to the LEFT of the spanning cell, so the last "
                    "one is what points at it -- and that is a %s, not an arrow. "
                    "Nothing reaches the span. It would still DRAW, reading as "
                    "\"this row feeds the spanning cell\", and every check would "
                    "pass. If you meant the spanning cell to branch OUT into these "
                    "rows, that is a FORK, and a fork roots at SLOT 0: move the "
                    "rowspan cell to the front of its row and give each covered "
                    "row its own arrow followed by its product."
                    % (r, sp["owner"], sp["slot"], _cell_kind(rows[r][-1])))
            rows[r] = ([ghost] + rows[r]) if sp["fork"] else (rows[r] + [ghost])

    # TWO GHOSTS IN A ROW IS NORMAL; TWO **FORK** GHOSTS IS NOT.
    # A row can have a fork's ghost at its front and a merge's at its end --
    # the row leaves one spanning cell and arrives at another, and both ghosts
    # are where they belong. Two forks is the shape that cannot work: both are
    # PREPENDED, so the second one pushes the first off index 0, and index 0 is
    # the whole of how align_x recognises "this row has no `pre`". The row is
    # then treated as a merge, align_x goes looking for an arrow in front of the
    # spanning cell to lengthen, finds a ghost, and raises -- naming a missing
    # arrow, which is the symptom two steps downstream of the cause.
    #
    # And two forks can only share a covered row when one fork's OWNER row lies
    # inside the other's span: that is precisely a branch whose content is
    # another fork. So this test names the unsupported shape exactly.
    for r, row in enumerate(rows):
        owners = [spans[c["_span"]]["owner"] for c in row
                  if "ghost" in c and spans[c["_span"]]["fork"]]
        if len(owners) > 1:
            raise ValueError(
                "rows[%d] is covered by %d forks at once (rooted in rows %s): a "
                "fork whose branch content is itself a fork is not supported. "
                "Both spanning cells would have to sit at the start of this row, "
                "and only one can. Draw the second divergence as its own figure, "
                "or flatten it into sibling branches of the first."
                % (r, len(owners), sorted(owners)))
    return rows, spans


def _row_extent(plan, upto):
    """Right edge of a row's items before index `upto`, in row-local x."""
    right = 0.0
    for it in plan["items"][:upto]:
        if it["type"] == "arrow":
            right = max(right, it["cx"] + it["length"] / 2 + it["ext"]["box_right"])
        else:
            right = max(right, it["x"] + it["w"])
    return right


def merge_lead(cfg) -> float:
    """How far a MERGE's spanning cell sits past the spine: the ONE arrow a
    convergence ends in (see align_x).

    Its own function because layout_linear's wrap needs it too, before any x
    exists: a merge cell costs this arrow on top of its own width. Measured as a
    plain cell, one route's [arrow, 4] merge came out 30 pt wider than its cells
    add up to. Two spellings of one length is how the width the wrap reserves
    stops being the width align_x spends.
    """
    return cfg.get("merge_arrow", cfg["min_arrow"])


def align_x(plans, shifts, spans, cfg, S):
    """Make every covered row agree with its span's x, in FINAL coordinates.

    Mutates `shifts` and, where a row's pre must stretch, the arrow item that
    absorbs it. Returns the spans annotated with x0/x1 for pass 2.
    """
    gap = cfg["gap"]
    for si, sp in enumerate(spans):
        owner = plans[sp["owner"]]
        span_it = next(it for it in owner["items"]
                       if it["type"] == "s" and it["id"] == sp["id"])
        covered = list(range(sp["owner"], sp["last"] + 1))

        span_it["_spanning"] = True

        # --- where each covered row's `pre` ends, in FINAL coords -------------
        pre_ends = {}
        for r in covered:
            plan = plans[r]
            idx = next((n for n, it in enumerate(plan["items"])
                        if it.get("_span") == si or
                        (r == sp["owner"] and it is span_it)), None)
            if idx is None:
                continue
            pre_ends[r] = (_row_extent(plan, idx) + shifts[r], idx)

        # --- the span's left edge: the RIGHTMOST pre, so none is cut short ----
        # A row with an EMPTY pre (the fork case) says nothing about where the
        # span goes -- it has no content to the left of it. Only the owning row
        # and any row that actually reaches the span constrain x0. Letting an
        # empty pre contribute its bare shift would pin x0 to whatever
        # row_align had already decided.
        # A fork's branches start one `gap` beyond the cell. A merge is not the
        # mirror image: convergence ends in ONE arrow pointing at the product,
        # so the spine has to sit a whole arrow-length back from it, not a gap.
        lead = gap if sp.get("fork", True) else merge_lead(cfg)
        sp["lead"] = lead
        cands = [e + lead for r, (e, idx) in pre_ends.items() if idx]
        own_e, own_idx = pre_ends[sp["owner"]]
        cands.append(own_e + (lead if own_idx else 0.0))
        x0 = max(cands)

        # --- absorb every other row's shortfall -------------------------------
        for r, (end, idx) in pre_ends.items():
            if idx == 0:                       # empty pre (fork): nothing to stretch
                continue
            arrows = [it for it in plans[r]["items"][:idx] if it["type"] == "arrow"]
            if not arrows:
                raise ValueError(
                    "row %d must reach the spanning cell but has no arrow in "
                    "front of it to lengthen" % r)
            a = arrows[-1]                     # the one nearest the span
            deficit = x0 - lead - end
            if deficit > 1e-9:
                a["length"] += deficit
                a["cx"] += deficit / 2.0       # left edge stays put
            # A converging route carries no arrowhead of its own: the head
            # belongs to the single arrow at the point of convergence, pointing
            # at the product. Two headed arrows aiming at a bare line is not how
            # a convergent scheme is drawn. Copied, not mutated in place, so the
            # caller's recipe dict is untouched.
            a["a"] = {**a["a"],
                      "style": {**(a["a"].get("style") or {}), "heads": "none"}}

        # --- the span itself, and every covered row's `post` ------------------
        x1 = x0 + span_it["w"]
        sp["x0"], sp["x1"] = x0, x1
        span_it["x"] = x0 - shifts[sp["owner"]]        # pass 2 re-adds the shift
        for r in covered:
            if r not in pre_ends:
                continue
            end, idx = pre_ends[r]
            if idx != 0:
                continue                       # merge: nothing after the span
            first = plans[r]["items"][idx + 1:]
            if not first:
                continue
            lead = first[0]
            cur = ((lead["cx"] - lead["length"] / 2 - lead["ext"]["box_left"])
                   if lead["type"] == "arrow" else lead["x"]) + shifts[r]
            shifts[r] += (x1 + gap) - cur
    return spans


def fork_arrow_length(rows, spans, measure_arrow):
    """{(row, cell): length} forcing every branch of one fork to one arrow length.

    Arrow length tracks the widest caption line, so unequal branch captions give
    a ragged product column: measured at 127.64 / 246.21 / 127.64 for three
    branches whose captions differ, a spread of 118.57 pt. The tails are already
    common -- they all leave the spine -- so it is the HEADS that scatter, and a
    product sits at its arrow's head.

    WHY THIS IS A PRE-PASS AND NOT align_x's STRETCH. That stretch is two lines
    (`length += deficit`, `cx += deficit/2`) and the arithmetic is right for a
    fork too, but it cannot be reached by deleting the `idx == 0` skip:

      * It runs AFTER pass 1, and it only gets away with that because the arrow
        it lengthens is the LAST item of a merge's `pre` -- nothing follows it in
        the row, so nothing downstream is left stale. A fork's branch arrow is
        followed by its product, whose x pass 1 has already accumulated from the
        old length. Reusing the stretch there would need a third operation that
        pushes the rest of the row along, and `post` cannot supply it: it shifts
        the whole row, which moves arrow and product together and cannot change
        the distance between them.
      * That path also forces `heads: "none"`, because a converging route hands
        its arrowhead to the single arrow at the point of convergence. A fork's
        branches each keep their own head.

    Decided here instead, before pass 1 accumulates any x, which is where it
    has to happen: `measure_arrow` is the same function, the length
    is simply the max of what it returns, and no retrofit is needed at all.

    THE FIRST arrow of each covered row, not every arrow. A branch may be a chain
    of several steps, and its later steps are steps within the
    branch, sized by their own captions like any other row. Only the arrow leaving
    the spine is a branch arrow.

    One fork at a time: two unrelated forks in one figure each get their own max,
    because equalizing across them would tie two divergences that share nothing.
    A row with no arrow at all contributes nothing and receives nothing -- it has
    no connection to the spine either way, which is a separate defect if it
    happens. Empty dict when the figure has no fork, so pass 1 is untouched.
    """
    # A FORK WHOSE ROWS ALSO END IN A MERGE IS LEFT ALONE. When one row is
    # [fork root, arrow, merge root], that single arrow is both the
    # fork's branch arrow and the merge's `pre` arrow -- and align_x will then
    # stretch it to reach the merge root. Measured on such a figure:
    # equalizing raised its natural length 23.12 pt, the stretch took it to
    # 634.07 while the other branch stayed at 98.12, and the figure grew 23 pt
    # wider. The equality does not survive, so it is not bought; it is only paid
    # for. Skipping keeps that figure byte-identical and keeps this function from
    # promising something the next pass revokes.
    merge_rows: set = set()
    for other in spans:
        if not other.get("fork", True):
            merge_rows.update(range(other["owner"], other["last"] + 1))

    out: dict = {}
    for sp in spans:
        if not sp.get("fork", True):
            continue
        if any(r in merge_rows for r in range(sp["owner"], sp["last"] + 1)):
            continue
        picks = []
        for r in range(sp["owner"], sp["last"] + 1):
            ci = next((k for k, c in enumerate(rows[r]) if "arrow" in c), None)
            if ci is None:
                continue
            picks.append((r, ci, measure_arrow(rows[r][ci]["arrow"])))
        if len(picks) < 2:
            continue                    # nothing to agree with
        longest = max(p[2] for p in picks)
        for r, ci, _nat in picks:
            out[(r, ci)] = longest
    return out


def fork_stub_y(shaft_ys):
    """Y of the stub -- the run from the spanning cell out to the spine.

    ODD N SNAPS TO THE MIDDLE SHAFT. The stub must be COLINEAR with that arrow,
    because that is what it reads as: that arrow, extended back to the cell. A
    few points off does not read as "slightly off", it reads as a mistake --
    measured at 5.95 pt on one fork, which is plainly visible.

    DO NOT SIMPLIFY THIS TO AN AVERAGE: the first/last midpoint is exactly
    that mistake. Those two agree only when the rows are evenly
    spaced, and covered rows are deliberately NOT evenly spaced: align_y below
    gives each covered gap `downs + row_clear + ups` from its own two rows, so
    one branch carrying a taller product or an extra sublabel moves the midpoint
    and leaves the middle arrow behind. On that fork the gaps are 121.1 and
    133.0 pt, and the midpoint lands 5.95 pt below the middle shaft.

    EVEN N has no middle arrow to be colinear with, so there is nothing to snap
    to and the midpoint of the middle PAIR is the honest answer -- it sits
    symmetrically between the two arrows nearest the centre. Note this is the
    middle pair, not the outer pair: with four unevenly spaced branches those
    differ too.
    """
    n = len(shaft_ys)
    mid = n // 2
    if n % 2:
        return shaft_ys[mid]
    return 0.5 * (shaft_ys[mid - 1] + shaft_ys[mid])


def align_y(plans, spans, ups, downs, row_base, cfg):
    """Two-pass row heights.

    Pass 1 (already done by the caller) heights every row from its own cells.
    A spanning cell was NOT counted there -- it belongs to no single row -- so
    it may be taller than the rows it covers. Distribute that excess over the
    covered gaps, then centre the cell on the rows it now spans.

    "Centre the cell over the rows it spans" is circular until the excess has
    been distributed, which is why this is two passes and not one.
    """
    for sp in spans:
        owner = plans[sp["owner"]]
        span_it = next(it for it in owner["items"]
                       if it["type"] == "s" and it["id"] == sp["id"])
        i, j = sp["owner"], sp["last"]

        # --- the rows a span covers are BRANCHES, not independent rows -------
        # `row_gap_y` (130) is a floor with no relation to content; `row_clear`
        # (28) is the real relation, applied to the rows' actual extents. The
        # floor is kept only so existing figures do not move. A fan's covered
        # rows carry one arrow and one
        # product, so they NEED about 77 and the floor charges 130 -- a third of
        # the figure comes out as dead space, and branches spread that far apart
        # read as four separate reactions instead of one divergence, which is
        # the exact complaint the fan was drawn to fix.
        #
        # So covered rows are spaced by what they need. Rows the span does not
        # touch keep the floor, so nothing without a rowspan can move.
        inc = [0.0] + [row_base[r] - row_base[r - 1]
                       for r in range(1, len(row_base))]
        for r in range(i + 1, j + 1):
            inc[r] = downs[r - 1] + cfg["row_clear"] + ups[r]
        for r in range(1, len(row_base)):
            row_base[r] = row_base[r - 1] + inc[r]

        # A fan reads as a fan when the spanning cell lines up with the ARROWS,
        # so it is centred on the baseline midpoint -- NOT on the content extent.
        # Those two differ whenever the first row's captions above and the last
        # row's sublabels below are not symmetric, which is the normal case:
        # centred on the extent, the cell sags 11.9 units towards whichever row
        # carries more text.
        #
        # Centred there, the cell reaches min(ups[i], downs[j]) past the outer
        # baselines before it runs out of room, so that -- doubled -- is what
        # the covered rows actually offer.
        have = (row_base[j] - row_base[i]) + 2 * min(ups[i], downs[j])
        need = span_it["h"] + 2 * cfg.get("rowspan_pad", 0.0)
        if need > have:
            extra = (need - have) / (j - i)          # per covered gap
            for r in range(i + 1, len(row_base)):
                row_base[r] += extra * min(r - i, j - i)
        # NOT (row_base[i] + row_base[j]) / 2 -- see fork_stub_y. Read off the
        # DISTRIBUTED row_base, so the snap follows the rows to wherever the
        # excess above put them; that is the second half of why this is two
        # passes rather than one.
        sp["cy"] = fork_stub_y([row_base[r] for r in range(i, j + 1)])
        span_it["_span_cy"] = sp["cy"]
    return row_base


def emit_bus(spans, plans, shifts, row_base, cfg, gid, ink=False):
    """The line joining a spanning cell to the tails (or heads) of its arrows.

    WHY THIS EXISTS -- found by looking at the render, not by reasoning
    ------------------------------------------------------------------
    Horizontal arrows leave a spanning cell's right edge at each covered row's
    OWN baseline. The cell is centred on those baselines. So in a four-row fan
    the outer arrows start 115 units above and below the cell and NOTHING
    touches it: the figure comes out as a molecule beside four unrelated
    reactions, which is the exact reading the fan was drawn to replace.

    A test that checks the four arrows share a start x misses it -- that is
    true, and irrelevant. Sharing an x is not being connected.

    The fix keeps every arrow an ordinary horizontal arrow (so `above`,
    `below`, `yield`, `reagent`, composite arrows and ink_extra all still work
    on them) and adds the connection as two plain line segments:

        fork              merge
        [cell]--+          +--[cell]
                |          |
                +-->       <--+
                |          |
                +-->       <--+

    A line is not an arrow, so this needs no arrow machinery at all -- which is
    why the horizontal decision survives the defect it turned up.

    WHAT IT HANDS BACK TO THE ALLOCATOR (`ink`)
    -------------------------------------------
    The shape is an L (or a T): a tall thin spine, plus a short thin stub joining
    it to the cell. Its axis-aligned HULL is therefore mostly the empty corner
    between the two, and reserving the hull declares that whole corner occupied.

    With ink=True the two segments are reserved as themselves instead. That
    matters wherever something else is legitimately placed in the corner -- in a
    radial figure the corner IS the corridor between the hub and the trunk, only
    `lead` (24) units wide, and a hull reservation there reports anything sitting
    in it as a collision with a line that is 24 units away. Measured on one
    radial figure: the hull came back 26.0 x 320.7 for a drawing whose ink
    is two 2-unit-wide segments, and the centre label was reported 9.6 deep into
    it -- 9.6 being the LABEL's height, not any distance to the trunk, which is
    why cutting the label's length by three quarters moved the reading by 0.0.

    Default False, and the linear engine keeps it: there the bus sits in the
    gutter between two columns with nothing else in it, so the hull costs nothing
    and switching it is a change to that engine's fixtures for no gain here.

    The pad is +-1 either way -- the same figure the hull already used, so
    nothing new is measured.
    """
    out, boxes = [], []
    gap = cfg["gap"]
    for sp in spans:
        i, j = sp["owner"], sp["last"]
        y0, y1 = row_base[i], row_base[j]
        cy = sp["cy"]
        fork = sp.get("fork", True)
        # the arrows tails sit one gap outside the cell on the side they leave
        lead = sp.get("lead", gap)
        bx = (sp["x1"] + lead) if fork else (sp["x0"] - lead)
        ex = sp["x1"] if fork else sp["x0"]
        lw = cfg.get("bus_width", 1.0)
        out.append(C._line(gid, bx, y0, bx, y1, lw)); gid += 1          # the spine
        if abs(bx - ex) > 1e-9:
            if fork:
                # the branches carry the heads; this only joins the cell to them
                out.append(C._line(gid, ex, cy, bx, cy, lw)); gid += 1
            else:
                # the ONE arrow of a convergence, head on the product
                out.append(C.arrow(gid, (bx + ex) / 2.0, cy, ex - bx)); gid += 1
        if ink:
            boxes.append((bx - 1, y0 - 1, bx + 1, y1 + 1))            # the spine
            if abs(bx - ex) > 1e-9:                    # ...only if one was drawn
                boxes.append((min(bx, ex) - 1, cy - 1, max(bx, ex) + 1, cy + 1))
        else:
            boxes.append((min(bx, ex) - 1, y0 - 1, max(bx, ex) + 1, y1 + 1))
    return out, boxes, gid
