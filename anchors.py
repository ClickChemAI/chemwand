"""anchors.py — "the pyridinium ring of cell 2" -> a point on the page.

Everything decorative in this pipeline already names
atoms by SMARTS rather than by index or coordinate (see match_node_ids: an index
dies when the SMILES is rewritten, a coordinate dies when a substituent changes),
but that naming only ever reached INSIDE one fragment. `decorate.apply` resolves
against one structure's index_map and knows nothing about where a packer put it,
so an electron-pushing arrow whose two ends live in two different cells could not
be drawn without it -- e.g. a piperidine N in one cell attacking a pyridinium ring
in another.

WHAT IS NEEDED is not coordinates. `C.node_points()` reads a fragment's own atom
positions, and every engine knows the (dx, dy) it handed to C.place(). This
module writes the two down together. An engine's own position record (such as
layout_linear's `struct_pos`) is per-CELL and local, not per-ATOM, and not
reachable from outside.

    record()  called once per cell, right after C.place(), with the PLACED string
    points()  a SMARTS -> page points, resolved through the same match_node_ids
              that every highlight goes through

SMARTS, NOT INDICES, and this is not a preference. `abbreviate=True` collapses a
pendant group into one superatom node, so the atom COUNT changes; and a caller
cannot know the engine's numbering anyway. The interface has to be the one thing
that survives both, which is the same argument match_node_ids' docstring makes at
length.

COORDINATE FRAME: LAYOUT, the engine's own, the one every Allocator reservation
is in -- NOT the page. It is the frame the caller of record() was working in, and
converting here would need finalize's shift, which does not exist yet when a cell
is placed.
"""

from __future__ import annotations

import math

import cdxml_core as C
import decorate as _DEC
import smiles_to_frag as S


class Anchors:
    """Per-figure: which atoms of which cell ended up where.

    One instance per render, handed in by the tool layer the way `alloc` and
    `shift` are. Empty unless somebody records into it, so an engine that has not
    been wired costs nothing and answers nothing -- rather than answering wrongly.
    """

    def __init__(self):
        self._pts: dict = {}      # key -> {placed node id: (x, y)}
        self._src: dict = {}      # key -> (smiles, index_map, base_id, base)
        self._els: dict = {}      # key -> {placed node id: {attribute: value}}

    # ---- filling it in --------------------------------------------------

    def record(self, key, placed_xml: str, base: int, *,
               smiles: str, index_map: dict, base_id: int = 1) -> None:
        """File one placed cell.

        `placed_xml` is what C.place() returned, so the coordinates in it are
        already the engine's, and `base` is the renumbering offset place() applied
        -- a matched id from the pre-place numbering becomes `id + base` here.

        Reading the PLACED string rather than translating the pre-place points is
        deliberate: it cannot drift from what was emitted, because it IS what was
        emitted. The same reason node_points exists at all instead of recomputing a
        layout.
        """
        if not smiles or index_map is None:
            return
        self._pts[key] = C.node_points(placed_xml)
        # AND the attributes, because an arrow END needs to know what it landed
        # on, not only where: an O or an H takes a wider setback than a bare
        # carbon vertex. Read off the SAME placed string as the points, so the
        # two cannot disagree about which node is which.
        self._els[key] = C.node_props(placed_xml)
        self._src[key] = (smiles, index_map, base_id, base)

    def __contains__(self, key) -> bool:
        return key in self._pts

    def keys(self):
        return tuple(self._pts)

    # ---- asking it questions --------------------------------------------

    def points(self, key, smarts: str, where: str | None = None) -> list:
        """[[(x, y), ...], ...] -- one list of page points per SMARTS match.

        A miss goes through S.ANCHOR_MISSES, the SAME channel a highlight's miss
        goes through, so "your query matched nothing" reads identically wherever
        it came from. Silence here would be the worst case: an arrow that is simply
        absent looks exactly like one nobody asked for.
        """
        src = self._src.get(key)
        if src is None:
            S.ANCHOR_MISSES.append(
                (where or ("cell %r" % (key,)), smarts,
                 "no structure was recorded under id %r in this figure; the ids "
                 "are %s" % (key, ", ".join(map(repr, sorted(map(str, self._pts))))
                             or "(none)")))
            return []
        pts = self._pts.get(key, {})
        return [[pts[i] for i in g] for g in self.ids_of(key, smarts, where)]

    def ids_of(self, key, smarts: str, where: str | None = None) -> list:
        """[[placed node id, ...], ...] -- one list per SMARTS match.

        points() alone throws the ids away. An arrow end needs them: how far
        short of its target it stops is decided by the ELEMENT it landed on,
        and only the id reaches that. Same resolver, same miss channel, so a
        miss reads the same from either.
        """
        src = self._src.get(key)
        if src is None:
            S.ANCHOR_MISSES.append(
                (where or ("cell %r" % (key,)), smarts,
                 "no structure was recorded under id %r in this figure; the ids "
                 "are %s" % (key, ", ".join(map(repr, sorted(map(str, self._pts))))
                             or "(none)")))
            return []
        smiles, imap, base_id, base = src
        pts = self._pts[key]
        out = []
        for ids in S.match_node_ids(smiles, smarts, imap, base_id=base_id,
                                    where=where or ("cell %r" % (key,))):
            here = [i + base for i in ids if (i + base) in pts]
            if here:
                out.append(here)
        return out

    def centre_kind(self, key, smarts: str, where: str | None = None):
        """(point, ANNOTATE_INSET bucket) for the FIRST match, else (None, None).

        centre() plus the one extra fact an arrow end needs: what it has to
        clear. One atom -> its element decides; two -> the bond between them;
        three or more -> a ring centre, where the head belongs ON the centroid
        and nothing is pulled back.
        """
        gs = self.ids_of(key, smarts, where=where)
        if not gs:
            return None, None
        g = gs[0]
        pts = self._pts[key]
        pt = (sum(pts[i][0] for i in g) / len(g),
              sum(pts[i][1] for i in g) / len(g))
        if len(g) != 1:
            return pt, _DEC.endpoint_kind(natoms=len(g))
        pr = self._els.get(key, {}).get(g[0], {})
        return pt, _DEC.endpoint_kind(int(pr.get("Element", 6)),
                                      int(pr.get("Charge") or 0),
                                      node_type=pr.get("NodeType"))

    def centre(self, key, smarts: str, where: str | None = None):
        """The mean of the FIRST match, or None.

        First rather than all: two ends of one arrow are two points, and a query
        that lands twice in one cell is a query that has not said which one it
        means. The miss report above names it either way.
        """
        got = self.points(key, smarts, where=where)
        if not got:
            return None
        g = got[0]
        return (sum(p[0] for p in g) / len(g), sum(p[1] for p in g) / len(g))


# --- the cross-structure half ------------------------------------------------
# The within-structure variants (a pointing label, a bond leader, a pushing arrow,
# a regiochemistry arc) all live in decorate.apply and go INSIDE the fragment, so
# C.place() carries them. This one cannot: it belongs to neither fragment.
#
# CONSEQUENCE, STATED ONCE: a page-level curve does NOT move if somebody later
# drags either structure in ChemDraw. That is acceptable for a figure written in
# one pass and it is the only place in this pipeline where a mark does not travel
# with its owner -- written down so it is not rediscovered as a bug. It IS
# reserved, so nothing else lands on it.

ARROWS = _DEC.ANNOTATE_ARROWS          # "full" | "half" | "none", one definition
DEFAULT_BULGE = _DEC.ANNOTATE_BULGE
DEFAULT_SIZE = _DEC.ANNOTATE_TEXT_SIZE


def between(frame: "Anchors", spec: dict, gid: int, tid: int, where: str = "?",
            free=None, moved: list | None = None):
    """(xml, ink, [(rect, owner, kind), ...]) for one cross-structure annotation.

    `ink` is the curve flattened to a polyline, for the engine to hand to
    reserve_polyline. NOT one bounding box: the axis-aligned hull of a bow is
    mostly air -- measured, a hull reservation came to 85 x 26 pt for a
    hairline and swallowed three compound numbers underneath it. Same lesson
    reserve_polyline was written for and radial's spokes cite.

    Returns boxes to reserve rather than reserving them: this module has no
    allocator and should not grow one. The engine owns its ledger; this is a guest.

    `free(rect) -> bool` lets the caption step clear of what is already placed.
    Needed because THE FAR END IS OFTEN INSIDE SOMETHING: a pushing arrow arrives
    at a ring, so its target point is the ring's centroid, and "just past the
    point, along the chord" is still in the middle of the ring. Without this the
    label prints over the structure it is naming. `moved` receives a record when
    the label had to leave its conventional spot, for the same reason
    layout_linear.CAPTION_MOVED exists -- a silent relocation swaps a reported
    overlap for an unexplained gap.
    """
    a, kind_a = frame.centre_kind(spec["from"]["cell"], spec["from"]["on"],
                                  where="%s.from" % where)
    b, kind_b = frame.centre_kind(spec["to"]["cell"], spec["to"]["on"],
                                  where="%s.to" % where)
    if a is None or b is None:
        return "", [], []      # already reported through ANCHOR_MISSES
    head = ARROWS.get(str(spec.get("arrow", "full")).lower(), "MISSING")
    if head == "MISSING":
        raise ValueError(
            "%s: arrow=%r; one of %s. 'half' is the fishhook a single-electron "
            "arrow draws." % (where, spec.get("arrow"), ", ".join(sorted(ARROWS))))
    colour = spec.get("color", "red")
    bulge = float(spec.get("bulge", DEFAULT_BULGE))
    pts = _DEC.bow(a, b, bulge)
    # NATIVE <arrow>, not <curve>, the same as the in-structure annotate: a
    # <curve> carrying an arrowhead draws the HEAD ALONE with no line -- measured
    # three ways, see the list in decorate._push_arrow -- so a CROSS-CELL
    # electron-pushing arrow would come out as a lone triangle or as nothing at
    # all, while the tool reported ok with no warnings. `_push_arrow` negates the
    # bulge on its way into arc_chord, which is exactly what keeps this bowing to
    # the same side as bow().
    # SAME SETBACKS AS AN IN-STRUCTURE ARROW. One table of clearances
    # (_DEC.ANNOTATE_INSET) for both paths -- two sets of aesthetics for the same
    # object would be the worse bug.
    out = [_DEC._push_arrow(gid, a, b, bulge, head, colour,
                            back=_DEC.endpoint_backs(kind_a, kind_b))]
    own = "annot:%s->%s" % (spec["from"]["cell"], spec["to"]["cell"])
    # The BOW flattened, not the chord: the bulge is the half a reader notices, and
    # reserving the chord would leave it unaccounted.
    ink = C.curve_polyline(pts)
    boxes = []
    if spec.get("text"):
        import layout_linear as _LL       # deferred: text metrics only
        size = int(spec.get("size", DEFAULT_SIZE))
        w, h = _LL.block_box([str(spec["text"])], size)
        # Outward along the chord from the far end, stepping until clear. The label
        # of a pushing arrow belongs where the arrow ARRIVES -- that is what it
        # names -- but the arrival point is frequently a ring centroid, i.e. inside
        # the structure, so the first candidate is routinely occupied.
        ux, uy = b[0] - a[0], b[1] - a[1]
        d = math.hypot(ux, uy) or 1.0
        ux, uy = ux / d, uy / d
        base = max(w, h) / 2.0 + _DEC.ANNOTATE_TEXT_PAD
        step = max(w, h) * 0.6
        tx = ty = None
        for k in range(6):
            cx = b[0] + ux * (base + k * step)
            cy = b[1] + uy * (base + k * step)
            r = (cx - w / 2.0, cy - h / 2.0, w, h)
            if free is None or free(r):
                tx, ty = cx, cy
                if k and moved is not None:
                    moved.append((own + ".text", "along the arrow",
                                  k * step, "%d of 5 steps" % k))
                break
        if tx is None:                      # nowhere clear: keep the convention
            tx = b[0] + ux * base
            ty = b[1] + uy * base           # and let the allocator report it
        out.append(C.styled_text(tid, tx, ty - h / 2.0 + 0.716 * size,
                                 [str(spec["text"])], size=size, align="center",
                                 color=C.color_index(colour)))
        boxes.append(((tx - w / 2.0, ty - h / 2.0, w, h), own + ".text", "caption"))
    return "\n".join(out), ink, boxes
