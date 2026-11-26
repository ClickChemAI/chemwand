"""nested.py — a cell that holds a small layout INSIDE it, stacked VERTICALLY.

Lets a row position hold a vertical stack and, through it, draw two inputs with
ONE arrow (rather than across two rows) and a bracket around a whole mixture
rather than around one structure.

WHY IT TOUCHES THE PACKER ALMOST NOWHERE
----------------------------------------
`layout_linear`'s pass 1 never asks a CELL how big it is.  It asks `S`:

    bb = S[cell["s"]][1]                       # layout_linear L355
    w, h = bb[2]-bb[0], bb[3]-bb[1]

`S` is just `id -> (object_string, bbox)`.  Nothing requires that object string
to be a single <fragment>.  `C.group()` wraps several already-built objects and
gives the result its own BoundingBox, and `C.place()` finds that BoundingBox
first (it precedes the children in the string), renumbers, and translates.

So a stack is composed BEFORE render_linear runs and registered in `S` under a
synthetic id.  The packer then treats it as an ordinary structure, which means
`bracket`, `label`, `sublabels` and `under` all keep working with no new code:
the bracket in particular is drawn around `it["h"]` and does not care whether
that height is one molecule or three.

TWO CONSTRAINTS, both discovered from the source and both load-bearing
---------------------------------------------------------------------
1. ARROWS INSIDE A STACK ARE ALLOWED ONLY BECAUSE `place()` MOVES THEIR 3D
   ENDPOINTS.  `C.translate()` moves p=, BoundingBox= and CurvePoints= and NOT
   the 3D endpoints an <arrow> carries, so a group placed with it alone would
   come out with its structures moved and its arrow stranded where it was
   composed.

   `place()` moves the 3D endpoints too, and that is byte-neutral for everything
   without an arrow -- a <fragment> has no such attributes, so for a fragment
   the two are the same function.

   Two things were measured in ChemDraw:
   a <group> holding an arrow IS a real group there -- it selects and drags as one,
   arrow included -- and it SURVIVES a re-save with the arrow still nested, which is
   NOT what happens to a <graphic> inside a <fragment> (that one gets flattened out).
   So "ChemDraw flattens nesting" is not a general rule; it depends on the container
   and the contained.

   The tear was also worse than it looked: ChemDraw reads a stranded arrow as a
   reaction step whose product is the whole group, so the embedded reaction data
   lied along with the picture.

   This is what opens horizontal nesting, and with it panels.py: a panel is
   structures + arrows + "+" + type, composed and placed as one cell.  The
   constraints below still hold and are what a panel has to obey.

2. THE COMPOSITE BBOX IS THE UNION OF THE MEMBERS' OWN bboxes, never re-derived
   from coordinates.  `fragment()` writes the corrected extent (atom-label
   overhang included) into the BoundingBox attribute; scanning p="" instead
   would silently drop that overhang again.

ID BUDGET
---------
pass 2 places structures with `base += 1000`, so everything inside one stack
must fit in fewer than 1000 ids or the next cell collides with it.  Members are
packed tightly and the total is asserted.

DEPTH
-----
A member is a structure or a "+".  A stack inside a stack is NOT accepted --
nesting is capped at 2 and this is level 1; the cap is enforced here rather
than discovered later as a pathological figure.
"""

from __future__ import annotations

import math
import re

import cdxml_core as C

_ID_RE = re.compile(r'\bid="(\d+)"')

# The <t> anchor is the FIRST LINE's baseline; the box top sits ASCENT_EM*size
# above it.  Same measured constant layout_linear uses -- duplicated rather than
# imported to keep this module free of a circular dependency on the engine that
# will call it.
ASCENT_EM = 0.716


def _max_id(s: str) -> int:
    ids = [int(m) for m in _ID_RE.findall(s)]
    return max(ids) if ids else 0


def _default_measure(lines, size):
    from layout_linear import block_box          # deferred: avoids the cycle
    return block_box(lines, size)


def compose_stack(members, cfg, *, measure=None, gid: int = 1):
    """Build one vertical stack as a single placeable object.

    members: list, in top-to-bottom order, of either
        {"frag": <object string>, "bbox": [x0,y0,x1,y1], "num": <label or None>}
        {"plus": True}

    Returns (object_string, bbox) -- exactly the shape `S` holds, so the caller
    can drop it straight into recipe["structures"].
    """
    measure = measure or _default_measure
    plus_gap = cfg["plus_gap"]
    plus_size = cfg["plus_size"]
    num_gap = cfg["num_gap"]
    num_size = cfg["num_size"]

    # ---- measure every member; the stack's width is the widest of them ----
    meas = []
    for m in members:
        # membership, NOT truthiness: an EMPTY stack is still a stack, and
        # `m.get("stack")` would wave it through to a KeyError on "bbox" later
        if "stack" in m:
            raise ValueError("a stack inside a stack is not supported (nesting is capped at depth 2)")
        if m.get("plus"):
            w, h = measure(["+"], plus_size)
            meas.append(dict(kind="plus", w=w, h=h))
            continue
        bb = m["bbox"]
        mw, mh = bb[2] - bb[0], bb[3] - bb[1]
        num = m.get("num")
        if num not in (None, False):
            nw, nh = measure([str(num)], num_size)
        else:
            nw = nh = 0.0
        meas.append(dict(kind="s", frag=m["frag"], num=num,
                         mw=mw, mh=mh, nw=nw, nh=nh,
                         w=max(mw, nw),                       # a long number can
                         h=mh + ((num_gap + nh) if nh else 0.0)))
    if not meas:
        raise ValueError("nested.compose_stack: empty stack")

    W = max(x["w"] for x in meas)

    # ---- place, top to bottom, every member centred on the stack axis ----
    # The gap between two members is plus_gap on each side of a "+", which is
    # the same constant the HORIZONTAL "+" already uses (layout_linear L373/378).
    # Borrowed from the case it mirrors rather than invented -- see the note at
    # the bottom of this file about whether a plus-free stack wants its own.
    children = []
    nid = 1
    y = 0.0
    for i, mm in enumerate(meas):
        if i:
            y += plus_gap
        if mm["kind"] == "plus":
            children.append(C.styled_text(nid, W / 2.0, y + ASCENT_EM * plus_size,
                                          ["+"], size=plus_size, align="center"))
            nid += 1
            y += mm["h"]
            continue
        # molecule
        placed = C.place(mm["frag"], nid, (W - mm["mw"]) / 2.0, y)
        nid = _max_id(placed) + 1
        children.append(placed)
        y += mm["mh"]
        # its own compound number, INSIDE the stack -- the source numbers each
        # member of a mixture individually, and the outer cell has only a
        # synthetic id to print (see expand_stacks, which suppresses it)
        if mm["nh"]:
            children.append(C.styled_text(
                nid, W / 2.0, y + num_gap + ASCENT_EM * num_size,
                [str(mm["num"])], size=num_size, bold=True, align="center"))
            nid += 1
            y += num_gap + mm["nh"]

    if nid >= 1000:
        raise ValueError("this stack holds too many objects (%d ids); split it into smaller stacks" % nid)

    # `C.group` writes the bbox with %.2f, so the number the PACKER reserves and
    # the number ChemDraw reads would otherwise differ by up to 0.005 -- and
    # %.2f rounds to nearest, i.e. half the time in the UNSAFE direction.
    # Ceil both to 2dp and hand the same value back, so the two agree exactly
    # and the error, where there is one, is always over-reservation.
    # (The asymmetric acceptance rule for a bbox: at or ABOVE, never
    # "as close as possible".)
    W = math.ceil(W * 100.0) / 100.0
    y = math.ceil(y * 100.0) / 100.0
    bbox = [0.0, 0.0, W, y]
    return C.group(gid, children, bbox=bbox), bbox


def expand_stacks(recipe: dict, cfg: dict, *, measure=None, prefix: str = "__stack") -> dict:
    """Rewrite `{"stack": [...]}` cells into synthetic `{"s": id}` cells.

    IDENTITY when no cell carries "stack" -- the same object comes back, so a
    recipe that does not use the feature cannot take a different code path and
    existing output cannot drift.  That acceptance test is made structural
    instead of checked by eye.

    Recipe syntax:

        {"stack": [{"s": 12}, {"plus": true}, {"s": 13}], "bracket": true}

    A member's compound number defaults to its id (the packer's own convention
    that a structure's id IS its number); "num": false suppresses it.
    """
    rows = recipe.get("rows") or []
    if not any("stack" in cell for row in rows for cell in row):
        return recipe

    S = dict(recipe["structures"])
    out_rows = []
    n = 0
    for row in rows:
        out_row = []
        for cell in row:
            if "stack" not in cell:
                out_row.append(cell)
                continue
            members = []
            for m in cell["stack"]:
                if m.get("plus"):
                    members.append({"plus": True})
                    continue
                frag, bb = S[m["s"]]
                members.append(dict(frag=frag, bbox=bb,
                                    num=(m["s"] if m.get("num", True) is True
                                         else m.get("num"))))
            obj, bb = compose_stack(members, cfg, measure=measure)
            sid = "%s%d" % (prefix, n)
            n += 1
            S[sid] = (obj, bb)
            # everything a cell may carry rides along unchanged; only the
            # compound number is turned off, because the members carry theirs
            new = {k: v for k, v in cell.items() if k != "stack"}
            new["s"] = sid
            new["num"] = False
            out_row.append(new)
        out_rows.append(out_row)
    return {**recipe, "structures": S, "rows": out_rows}


# ---------------------------------------------------------------------------
# OPEN:  the gap between two members with NO "+" between them currently
# uses plus_gap (14) as well.  There is no measured constant for that case -- it
# does not occur horizontally -- so it is borrowed, not derived.  If a plus-free
# stack ever appears, ask "compared to what" before fixing a number to it.
# ---------------------------------------------------------------------------
