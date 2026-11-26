"""Page-space allocator for the ChemWand CDXML pipeline.

Design stance: tolerant by default.

Two kinds of failure are deliberately kept apart:

  structural   someone emitted a primitive without reserving space for it,
               or two reservations were declared on top of each other.
               These are caller bugs and have nothing to do with geometric
               precision, so they are always detectable in pure Python.

  drift        the reserved rect turned out to be a bit smaller than what
               beautify / text_metrics / ChemDraw actually produced.
               This is prediction error, not a bug. It is absorbed silently
               unless it grows large enough to actually collide with a
               neighbour.

Nothing raises unless you opt in (Policy.strict / Policy.enforce_rid).
Run in reporting mode first, read Allocator.report(), then tighten.

TWO SHAPES, and the second one is narrow on purpose. `Rect` is axis-aligned and is
what almost everything is -- structures, labels, captions, hulls, page extents.
`OBox` is an oriented strip and exists for one thing: a SLANTED run, whose bounding
box is mostly empty and whose chunked approximation over-claims 6.4x
perpendicular to itself. See the note above OBox for the measurement. Horizontal and
vertical runs stay Rects, because there a strip and its bounding box are the same
rectangle.

Units are CDXML page units; BOND is the pipeline's target bond length, and
ARROW_INK_WIDTH is how wide a reserved run is (measured -- see the constant).
"""

from __future__ import annotations

import math
import os
import re
from dataclasses import dataclass, field
from enum import IntEnum
from typing import Iterable, Sequence

BOND = 30.0

# HOW THICK A RESERVED RUN IS, and it is read off a measurement rather than picked.
# ChemDraw's OWN recomputed BoundingBox, measured on its output for every arrow
# style at HeadSize 1000:
#
#   headless (the bare line)  reach 0.30   = half of the document's LineWidth 0.60
#   plain arrow               reach 2.50   = max single-side travel of the ink
#   plain arrow               half-height 2.00, which is what ChemDraw reports
#   plain arrow               end overhang 0.00 -- the head is drawn entirely
#                                                 BEHIND its Head3D point
#
# So 4.00 is 2 x the measured plain-arrow HALF-HEIGHT. That is its provenance, and
# it is NOT "a line width" -- LineWidth is 0.60, 6.7x smaller.
#
# WHY NOT THE 0.60 INK. A shaft driven clean through a structure scores the strip's
# own thickness as its depth, and ink_tol holds an arrow at 2.00. A 0.60-wide strip
# would score 0.60, land under the tolerance, and a total crossing would go silent
# -- the exact regression test_hairline_through_a_structure_is_not_invisible exists
# to prevent. Measured: deflating this to 2.0 fails 11 tests, two of them
# detection-alive tests.
#
# WHY NOT THE 2.50 REACH, which is the more conservative and arguably truer number:
# the same measurement puts the ink about 2.55 one side of the axis and 1.50 the other,
# so a symmetric 2.00 under-reserves roughly half a unit on the fat side. Raising
# this to 5.00 was measured and it MOVES SHIPPED FIGURES -- three further tests
# fail, one of them a cycle CDXML golden, because layout_cycle places by probing and
# a thicker strip pushes a caption. That is a separate decision with its own
# evidence.
#
# The end overhang being 0.00 is why a strip is NOT extended along its own axis.
# Isotropic inflation would add width/2 past each endpoint, which no measured
# arrow style does.
ARROW_INK_WIDTH = 4.00

Side = str      # "left" | "right" | "above" | "below"
Align = str     # "start" | "center" | "end"


class Severity(IntEnum):
    OK = 0
    NOTE = 1
    WARN = 2
    ERROR = 3


@dataclass(frozen=True)
class Rect:
    x: float
    y: float
    w: float
    h: float

    @property
    def right(self) -> float:
        return self.x + self.w

    @property
    def bottom(self) -> float:
        return self.y + self.h

    @property
    def cx(self) -> float:
        return self.x + self.w / 2.0

    @property
    def cy(self) -> float:
        return self.y + self.h / 2.0

    def union(self, other: "Rect") -> "Rect":
        x = min(self.x, other.x)
        y = min(self.y, other.y)
        return Rect(x, y,
                    max(self.right, other.right) - x,
                    max(self.bottom, other.bottom) - y)

    def overlap(self, other: "Rect") -> tuple[float, float]:
        """Interpenetration along each axis. Non-positive means no overlap."""
        dx = min(self.right, other.right) - max(self.x, other.x)
        dy = min(self.bottom, other.bottom) - max(self.y, other.y)
        return dx, dy

    def depth(self, other: "Rect") -> float:
        """How far you would have to push to separate. 0.0 if disjoint."""
        dx, dy = self.overlap(other)
        if dx <= 0.0 or dy <= 0.0:
            return 0.0
        return min(dx, dy)

    def overlap_phrase(self, other: "Rect") -> str:
        """How deep, AND the shape of it: "1.5 units deep (35 wide x 1.5 tall)".

        The depth alone misleads, which is the whole reason this exists. It is a
        minimum translation distance -- the shorter way out -- so a pair of captions
        printed through each other across 35 units of shared baseline reports "1.5
        units deep" and reads as a rounding error. Anyone triaging needs the span
        too, and the span is free to compute.
        """
        dx, dy = self.overlap(other)
        return "%.1f units deep (%.0f wide x %.0f tall)" % (
            self.depth(other), max(0.0, dx), max(0.0, dy))

    def excess_over(self, other: "Rect") -> float:
        """How far self pokes out beyond other on the worst side."""
        return max(0.0,
                   other.x - self.x,
                   other.y - self.y,
                   self.right - other.right,
                   self.bottom - other.bottom)


# ---------------------------------------------------------------------------
# Oriented strips
# ---------------------------------------------------------------------------
# A DIAGONAL RUN IS NOT ITS BOUNDING BOX. Rect is axis-aligned, so the only way it
# could describe a slanted shaft was to chop the run into short pieces and let each
# piece's own box be nearly tight. That works, and it costs two things:
#
#   * the pieces are inflated ISOTROPICALLY, width/2 on all four sides, so the
#     clearance a run really claims PERPENDICULAR to itself is
#         L*|sin|*|cos| + (width/2)*(|sin| + |cos|)
#     which is 12.83 units at 45 degrees for L=20, width=4 -- against 2.00 for the
#     same ink laid flat. Measured, and it is what manufactures the phantom 175-180
#     error band that test_radial_auto_gaps pins as an artefact rather than ink.
#   * one contact is reported once per piece, so two near-collinear connectors
#     (0.35 degrees apart, measured) come back as fifty findings of one crossing.
#
# An OBox is the strip itself: a centre, a unit direction, a half-length along it
# and a half-thickness across it, with collisions decided by an exact
# separating-axis test. The perpendicular clearance is then `ht` at EVERY angle and
# nothing is chunked.
#
# WHAT THIS IS NOT. It does not replace Rect and Rect is not touched. Structures,
# labels, captions and every hull stay axis-aligned boxes, because that is what
# they are. An axis-aligned RUN also stays a Rect: for ux or uy on an axis the
# strip and its own bounding box are the same rectangle, so an OBox there would add
# a shape that answers every question identically. Both properties are pinned, in
# test_allocator, rather than asserted here.
@dataclass(frozen=True)
class OBox:
    """A rectangle rotated to lie along (ux, uy). Ink, not a hull.

    (cx, cy) is the centre, (ux, uy) a UNIT direction, `hl` the half-extent along
    it, `ht` the half-extent across it.

    Degenerate strips are deliberately not representable: a run with no length has
    no direction, so there is no axis to orient and reserve_polyline hands that
    case to Rect instead of inventing one.
    """

    cx: float
    cy: float
    ux: float
    uy: float
    hl: float
    ht: float

    def corners(self) -> list:
        px, py = -self.uy, self.ux
        return [(self.cx + self.ux * a + px * b, self.cy + self.uy * a + py * b)
                for a, b in ((-self.hl, -self.ht), (self.hl, -self.ht),
                             (self.hl, self.ht), (-self.hl, self.ht))]

    def aabb(self) -> Rect:
        """The axis-aligned hull -- what gets RESERVED as the reservation's `rect`.

        The reservation keeps a Rect so bbox() and the page extent are untouched by
        this whole change: "how far does the ink reach in x and y" is a question
        about the hull, and only the COLLISION question consults the strip.
        """
        cs = self.corners()
        xs = [p[0] for p in cs]
        ys = [p[1] for p in cs]
        return Rect(min(xs), min(ys), max(xs) - min(xs), max(ys) - min(ys))


def _shape_axes(shape):
    """(corners, separating axes) for a Rect or an OBox. Axes are already unit."""
    if isinstance(shape, OBox):
        return shape.corners(), ((shape.ux, shape.uy), (-shape.uy, shape.ux))
    return ([(shape.x, shape.y), (shape.right, shape.y),
             (shape.right, shape.bottom), (shape.x, shape.bottom)],
            ((1.0, 0.0), (0.0, 1.0)))


def sat_depth(a, b) -> float:
    """Minimum translation distance separating two convex quads. 0.0 if disjoint.

    THE SAME ANSWER AS Rect.depth FOR TWO RECTS, and the reason is worth stating
    because the whole change rests on it: both rectangles contribute the axes (1,0)
    and (0,1), so the two candidate overlaps ARE dx and dy, and the smaller of them
    is min(dx, dy) -- which is Rect.depth's own expression. Argued here, pinned
    against 2000 random pairs in test_allocator, because "provably identical" is
    worth exactly as much as the test that says so.
    """
    ca, aa = _shape_axes(a)
    cb, ab = _shape_axes(b)
    best = math.inf
    for ax, ay in aa + ab:
        pa = [p[0] * ax + p[1] * ay for p in ca]
        pb = [p[0] * ax + p[1] * ay for p in cb]
        ov = min(max(pa), max(pb)) - max(min(pa), min(pb))
        if ov <= 0.0:
            return 0.0
        best = min(best, ov)
    return 0.0 if best is math.inf else best


def run_shapes(pts, width: float = ARROW_INK_WIDTH):
    """[(rect, shape), ...] for the segments of a polyline. `shape` may be None.

    THE ONE PLACE that decides what a run occupies. reserve_polyline reserves these,
    polyline_hull unions them, and layout_radial's row solver tests captions against
    them -- so the allocator's account of a run and every other account of the same
    run cannot drift. That drift is not hypothetical: restating the padding by hand
    goes wrong as soon as the shape changes, and restating the shaft as a bare strip
    under-answers by 1.95 units because the real thing reaches half a width past
    each endpoint.

    PUBLIC for that third caller: a solver that has to measure the very same run
    needs the shapes themselves, and the alternative is a fourth restatement.
    """
    half = width / 2.0
    pts = [(float(a), float(b)) for a, b in pts]
    out = []
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        span = math.hypot(x1 - x0, y1 - y0)
        if span == 0.0:
            # A ZERO-LENGTH RUN. Nothing currently shipping reaches this branch
            # -- draw_snaking_scheme goes through layout_linear, whose turns are
            # routed and trimmed -- but it stays DECIDED anyway:
            # arrows.route can hand back a degenerate segment for any two points
            # that trim to the same place, and defaulting is how that becomes a
            # silent zero-area reservation.
            #
            # DECIDED, NOT DEFAULTED, and the decision is Rect. A run with no length
            # has no direction, so there is no axis to orient a strip along -- and
            # Rect is precisely the orientation-free shape, so this is not a
            # fallback so much as the right type. `width x width` about the point is
            # also the established reservation for this case, to the unit, so the
            # choice moves no figure.
            #
            # The two alternatives were weighed and rejected. A DEGENERATE OBox
            # (hl=0, some invented axis) reserves a zero-area sliver and would stop
            # reporting a contact that is reported today -- silence bought with an
            # invented direction. REFUSING would raise on a drawing that renders
            # today, turning a legal degenerate figure into an exception on no
            # evidence. Pinned by test_a_zero_length_run_stays_a_plain_rect.
            out.append((Rect(x0 - half, y0 - half, width, width), None))
            continue
        ux, uy = (x1 - x0) / span, (y1 - y0) / span
        if ux == 0.0 or uy == 0.0:
            # AXIS-ALIGNED, so the strip and its own bounding box are the SAME
            # rectangle. An OBox here would carry an orientation that changes no
            # answer, so a horizontal or vertical run stays a plain Rect and the
            # only reservations carrying oriented ink are the ones that need it.
            # Pinned by test_an_axis_aligned_run_stays_a_plain_rect.
            out.append((Rect(min(x0, x1) - half, min(y0, y1) - half,
                             abs(x1 - x0) + width, abs(y1 - y0) + width), None))
            continue
        shape = OBox((x0 + x1) / 2.0, (y0 + y1) / 2.0, ux, uy, span / 2.0, half)
        out.append((shape.aabb(), shape))
    return out


def polyline_hull(pts, width: float = ARROW_INK_WIDTH) -> "Rect | None":
    """The axis-aligned hull of everything reserve_polyline would reserve for `pts`.

    For callers keeping a SECOND, purely geometric account of the page extent beside
    the allocator's. There is one such caller (layout_linear's `boxes` ledger) and
    test_two_accounts holds the two together to 0.005, so it has to pad a run the
    way the allocator pads it. Restating the number instead goes wrong the moment
    the shape changes.

    None for a polyline with no segments, which is the same nothing
    reserve_polyline reserves for it.
    """
    boxes = [r for r, _s in run_shapes(pts, width)]
    if not boxes:
        return None
    hull = boxes[0]
    for b in boxes[1:]:
        hull = hull.union(b)
    return hull


@dataclass(frozen=True)
class Reservation:
    rid: int
    rect: Rect
    owner: str
    kind: str
    slack: float
    # The ORIENTED ink, when this reservation is a slanted run. None means the rect
    # IS the shape -- every structure, label, caption and axis-aligned run -- so
    # nothing that was a Rect has quietly become something else.
    shape: "OBox | None" = None


@dataclass(frozen=True)
class Finding:
    severity: Severity
    code: str
    message: str
    rids: tuple[int, ...] = ()
    amount: float = 0.0


@dataclass
class Policy:
    """Tolerances, in page units. Defaults are deliberately loose."""

    slack: dict[str, float] = field(default_factory=lambda: {
        "structure":  0.25 * BOND,   # 7.5
        "caption":    0.30 * BOND,   # 9.0
        "arrow":      0.20 * BOND,   # 6.0
        "label":      0.30 * BOND,   # 9.0
        "decoration": 0.40 * BOND,   # 12.0
    })
    default_slack: float = 0.25 * BOND

    # Two rects must interpenetrate by MORE than this before anyone complains.
    # A few units of clipping is invisible in ChemDraw; this is not a bug.
    overlap_tol: float = 0.25 * BOND

    # ...but a thin object can never interpenetrate by more than its own width,
    # so a hairline driven straight through a structure would score barely any
    # depth at all and slip under the tolerance. Cap the tolerance at this
    # fraction of the thinner object's short side.
    #
    # WRITTEN FOR ARROWS, and it must not be read as applying to text. A reserved
    # run is ARROW_INK_WIDTH across -- 4 units, and see that constant for where the
    # 4 comes from: it is twice ChemDraw's measured plain-arrow half-height, NOT "a
    # line width" (LineWidth is 0.60, 6.7x smaller).
    # A caption's short side, by contrast, is a LINE HEIGHT, and "half the thinner
    # side" then says half a line of overlap is acceptable. It is not: see
    # text_kinds below.
    thin_frac: float = 0.5

    # Two TEXT boxes tolerate NO interpenetration at all.
    #
    # Not strictness for its own sake -- a text box is an under-estimate of its own
    # ink. text_metrics fits a box 6.94 units tall for a 9-point line, where the same
    # line of Arial actually spans about 8.3 from ascender to descender. So two boxes
    # that merely touch are already printing through each other, and by the time the
    # boxes overlap by the 3.47 units `thin_frac` would allow, the glyphs have
    # collided for a third of their height.
    #
    # Under thin_frac a 27-character split feeder passes in silence: two labels
    # overlapping 35 units along their shared baseline and 1.5 units across it score
    # min(dx, dy) = 1.5 against a 3.47 tolerance. min() is the right question for two
    # chunky boxes clipping at a corner -- 4 units between two 60-unit cells really is
    # invisible, and that case is pinned by a test -- so the metric stays and the
    # TOLERANCE is what is wrong, for one kind of pair.
    text_kinds: set = field(default_factory=lambda: {"caption", "label"})
    text_overlap_tol: float = 0.0

    # How deep this kind may be interpenetrated before it shows. This is a
    # property of what the box CONTAINS, not of how big the box is:
    #
    #   caption/label  solid ink across the whole box. Any overlap prints
    #                  through, so zero.
    #   arrow          a 4-unit shaft inside a wide box; the box is mostly
    #                  air, and it can never interpenetrate deeper than its
    #                  own width anyway.
    #   decoration     same, a hairline.
    #   structure      bond lines are sparse; a few units clipped at a corner
    #                  between two 60-unit cells really is invisible, and a
    #                  test pins that.
    #
    # A pair takes min(), so adding a kind means answering ONE question --
    # how dense is its ink -- not one question per existing kind.
    ink_tol: dict = field(default_factory=lambda: {
        "caption": 0.0,
        "label": 0.0,
        "arrow": 2.0,
        "decoration": 2.0,
        "structure": 7.5,
        # A region is a SUMMARY rect covering content named elsewhere, so it
        # overlaps that content by construction -- see tolerates(). What this number
        # governs is region-vs-region.
        "region": 0.0,
        # A PANEL is different from a region and the two must not be confused: it is
        # real ink (a frame with a reaction inside), so it is strict against
        # everything. A panel dropped on the grid, or on another panel, is the one
        # failure the panel layer exists to prevent.
        "panel": 0.0,
    })

    # Kind pairs that are allowed to sit on top of each other by design.
    #
    # {arrow, caption} and {arrow, label} are NOT here. Both would stand in for
    # ownership -- "the arrow it annotates" -- which a kind pair cannot express,
    # so they would license a caption to sit on ANY arrow at any depth. Layout
    # already places an arrow's caption a gap above the shaft, and ownership is
    # threaded through the reservations instead (see related() below).
    #
    # The decoration pairs stay: a bracket really is drawn over the thing it
    # brackets, and that IS a property of the kind, not of a relation.
    tolerant_pairs: set = field(default_factory=lambda: {
        frozenset({"decoration", "structure"}),
        frozenset({"decoration", "decoration"}),
    })

    strict: bool = False        # raise on ERROR findings
    enforce_rid: bool = False   # raise when a primitive emits without a rid

    def slack_for(self, kind: str) -> float:
        return self.slack.get(kind, self.default_slack)

    def tolerates(self, a: str, b: str) -> bool:
        # A REGION IS DRAWN AROUND ITS OWN CONTENTS, so it overlaps every one of them
        # by construction -- enumerating pairs would mean listing every kind that
        # exists and re-listing it whenever a kind is added. Region-vs-REGION stays
        # strict: two panels on top of each other is the one thing this whole layer
        # is for. Written as a rule rather than a table because it is a property of
        # what a region IS, the same argument the decoration pairs make.
        if "region" in (a, b):
            return not (a == b == "region")
        return frozenset({a, b}) in self.tolerant_pairs

    def tol_for(self, a: str | None, b: str, ra: Rect, rb: Rect) -> float:
        """How deep two boxes of these kinds may interpenetrate before anyone minds.

        One place, asked by probe() and by nothing else, so "what counts as an
        overlap" has a single answer that can be read in one sitting.
        """
        if os.environ.get("ALLOCATOR_LEGACY_TOL"):
            if a is not None and a in self.text_kinds and b in self.text_kinds:
                return self.text_overlap_tol
            thin = min(ra.w, ra.h, rb.w, rb.h)
            return min(self.overlap_tol, thin * self.thin_frac)

        # ra and rb are deliberately unread: tolerance is a property of the
        # kinds, not of the sizes. Reading the boxes is what made the same
        # pair get different verdicts in different figures.
        return min(self.ink_tol.get(a, self.overlap_tol),
                   self.ink_tol.get(b, self.overlap_tol))


# ---------------------------------------------------------------------------
# Ownership
# ---------------------------------------------------------------------------
# Without the {arrow, caption} and {arrow, label} kind exemptions, a plain
# overlap test reports a flood of "caption over its OWN arrow" -- 8 of them
# across the cycle tests -- so ownership has to be threaded through the
# reservations.
#
# It does not have to be threaded very far. The owner strings ALREADY encode
# the relation:
#
#     step:n1->n2[0]                 names n1 and n2
#     node:n1.label                  belongs to n1
#     step:n0->n1.in[0].curve[0]     and .in[0].label share a path prefix
#
# so two reservations can be asked whether they are parts of one thing without
# any caller passing anything new. Two relations, and only two:
#
#   PREFIX   one owner's dotted path is a prefix of the other's. A caption and
#            the arrow it annotates; a bracket and the cell it brackets.
#
#   ENDPOINT a `step:` naming a species, versus that species' own `node:`.
#            An arrow legitimately reaches into the label of the node it starts
#            from or ends at.
#
# Step-vs-step is deliberately NOT related even when the two steps share a
# node. Two different arrows meeting at one species is exactly the crowding the
# packer exists to resolve, and exempting it would hide the thing worth seeing.
_SEG_IX = re.compile(r"\[\d+\]")


def _path(owner: str) -> list:
    """Dotted path with segment indices dropped: reserve_polyline appends [n] per
    SEGMENT, and one segment of an arrow is not a different owner from another.

    Unaffected by oriented strips, and it is worth saying why, because the
    indices are scarce: a straight run is one reservation rather than ten chunks,
    so most owners arrive as plain `[0]`. Stripping an index that is there and
    stripping one that is not are the same no-op, so `arrow:a->b` and
    `arrow:a->b[0]` reduce to the same component 0 either way -- which is what keeps
    a caption related to its own shaft. Checked by running the suite with the index
    dropped entirely: no new findings anywhere."""
    return _SEG_IX.sub("", owner).split(".")


def _species(owner: str) -> tuple:
    """(role, {names}) for the head of an owner path.

    'step:n1->n2' -> ('step', {'n1','n2'});  'node:n1' -> ('node', {'n1'})
    Anything else -> (head, set()), which relates to nothing.
    """
    head = _path(owner)[0]
    role, _, rest = head.partition(":")
    if role == "step":
        return "step", set(rest.split("->"))
    if role == "node":
        return "node", {rest}
    return head, set()


def related(a: str, b: str) -> bool:
    """Are these two reservations parts of the same object?

    SAME HEAD is the whole of the first relation. `a.b.label` and `a.b.curve`
    are SIBLINGS, not prefixes of one another -- comparing prefixes gets the
    feeder case wrong -- and every owner already carries its top-level object in
    component 0, so sharing that component is exactly "parts of one thing".
    """
    if _path(a)[0] == _path(b)[0]:             # SAME OBJECT
        return True
    ra, sa = _species(a)
    rb, sb = _species(b)
    if {ra, rb} == {"step", "node"} and sa & sb:   # ENDPOINT
        return True
    return False


class AllocatorError(RuntimeError):
    pass


def _require_owner(owner: str, kind: str, what: str) -> None:
    """Refuse an unnamed reservation, loudly.

    Deliberately an exception and not a Finding. A Finding would go into the same
    report the missing name has just made unreadable, and reserve() is called from
    layout code that is not looking at the report while it runs. The caller is
    wrong, it is wrong every time, and it is wrong at the moment of the call.

    It matters more because of ownership: related() reads the owner STRING to
    decide what not to report, so an empty owner is not merely unhelpful in a
    message -- it is an input to the exemption logic.
    """
    if owner is None or not str(owner).strip():
        raise AllocatorError(
            "%s() called with owner=%r (kind=%r). Every reservation must name "
            "what it belongs to -- an anonymous one occupies the map, absorbs "
            "collisions and reports them against nobody." % (what, owner, kind))


class Allocator:
    """Owns the occupancy map for one page region.

    It does NOT decide row width, pagination or gutters. It only guarantees
    that within the region it was handed, nothing meaningfully collides, and
    that bbox() reports the honest extent of everything that was reserved.
    """

    def __init__(self, policy: Policy | None = None):
        self.policy = policy or Policy()
        self._res: dict[int, Reservation] = {}
        self._next_rid = 1
        self.findings: list[Finding] = []
        self._drift: list[tuple[str, float]] = []   # (kind, excess) calibration log
        # (owner, other_owner, depth) for every overlap ownership decided not to
        # report. Surfaced by report(); see the branch in probe() for why.
        self._exempt: list[tuple[str, str, float]] = []

    # ---- allocation -----------------------------------------------------

    def reserve(self, rect: Rect, owner: str, kind: str = "structure",
                slack: float | None = None,
                allow: Iterable[int] = (),
                shape: "OBox | None" = None) -> Reservation:
        """`shape` is the ORIENTED ink, for a slanted run; `rect` stays its hull.

        Callers other than reserve_polyline pass no shape: with shape=None the rect
        is the shape.

        The MESSAGE keeps measuring the hulls on purpose. depth comes from the
        strip, which is the honest number, but overlap_phrase's "(4 wide x 35
        tall)" is what tells a reader the overlap runs for 35 units rather than
        nicking a corner -- and that span is a property of the boxes. It also
        carries the information a chunked model would convey by reporting the same
        contact fifty times, once per chunk.
        """
        _require_owner(owner, kind, "reserve")
        s = self.policy.slack_for(kind) if slack is None else slack
        for other, depth in self.probe(rect, kind=kind, ignore=allow,
                                       owner=owner, shape=shape):
            self._record(Severity.ERROR, "declared-collision",
                         f"{owner!r} ({kind}) declared on top of {other.owner!r} "
                         f"({other.kind}), {rect.overlap_phrase(other.rect)}",
                         (other.rid,), depth)
        res = Reservation(self._next_rid, rect, owner, kind, s, shape)
        self._res[res.rid] = res
        self._next_rid += 1
        return res

    @staticmethod
    def near_rect(anchor, side: Side, w: float, h: float, gap: float = 0.0,
                  align: Align = "center") -> Rect:
        """Where a w x h box goes beside `anchor`. Reserves NOTHING.

        Split out of reserve_near because a placer has to ASK before it commits:
        calling reserve_near per candidate and then deleting the reservation it did
        not want is a hack that reaches into `_res` from outside. Pure arithmetic
        has no such problem.

        `anchor` may be a Reservation or a bare Rect. The whole figure's extent is a
        Rect with no owner -- there is no reservation for "everything" -- and
        requiring one was incidental rather than meaningful.
        """
        a = anchor.rect if isinstance(anchor, Reservation) else anchor
        if side == "right":
            x, y = a.right + gap, Allocator._align_cross(a.y, a.h, h, align)
        elif side == "left":
            x, y = a.x - gap - w, Allocator._align_cross(a.y, a.h, h, align)
        elif side == "above":
            x, y = Allocator._align_cross(a.x, a.w, w, align), a.y - gap - h
        elif side == "below":
            x, y = Allocator._align_cross(a.x, a.w, w, align), a.bottom + gap
        else:
            raise ValueError(f"unknown side {side!r}")
        return Rect(x, y, w, h)

    def reserve_near(self, anchor, side: Side, w: float, h: float,
                     owner: str, kind: str = "caption", gap: float = 0.0,
                     align: Align = "center",
                     slack: float | None = None) -> Reservation:
        return self.reserve(self.near_rect(anchor, side, w, h, gap, align),
                            owner, kind, slack)


    def reserve_polyline(self, pts, owner: str, kind: str = "arrow",
                         width: float = ARROW_INK_WIDTH,
                         allow: Iterable[int] = ()) -> list[Reservation]:
        """Reserve a line or polyline as its own ink: ONE strip per segment.

        The axis-aligned bounding box of a DIAGONAL line is mostly empty space --
        that is why one cross-row connector could look like it collides with half
        the page. Chopping the run into short pieces so each piece's own box is
        nearly tight only approximates the answer; here it is answered exactly, by
        reserving the oriented strip and testing it with sat_depth. See the OBox
        note above for what the chopping costs.

        SEGMENTS, NOT RUNS, so an arc still arrives as many reservations: cycle
        feeders and ring arcs come through here as a chain of sampled chords, and
        each chord gets its own strip. A straight shaft is one segment and so one
        reservation.

        Pieces never flag each other: each is reserved with the previous ones on its
        allow-list.

        `allow` names reservations this run is MEANT to cross, and is added to that
        same list. A catalytic-cycle feeder arrow crossing the ring arc it feeds is
        the whole point of the motif, not a collision -- without this the report
        fills with dozens of by-design crossings and buries the real ones. Default
        empty.
        """
        # Checked HERE and not left to reserve(): each piece is reserved as
        # f"{owner}[i]", so a blank owner would arrive as the perfectly truthy
        # "[0]" and sail straight past the guard.
        _require_owner(owner, kind, "reserve_polyline")
        out: list[Reservation] = []
        base_allow = list(allow)
        for rect, shape in run_shapes(pts, width):
            out.append(self.reserve(
                rect, f"{owner}[{len(out)}]", kind,
                allow=base_allow + [r.rid for r in out], shape=shape))
        return out

    @staticmethod
    def _align_cross(origin: float, span: float, size: float, align: Align) -> float:
        if align == "start":
            return origin
        if align == "end":
            return origin + span - size
        return origin + (span - size) / 2.0

    # ---- queries --------------------------------------------------------

    def probe(self, rect: Rect, kind: str | None = None,
              ignore: Iterable[int] = (), tol: float | None = None,
              pad: float = 0.0, owner: str | None = None,
              shape: "OBox | None" = None
              ) -> list[tuple[Reservation, float]]:
        """Reservations this rect meaningfully collides with, worst first.

        `tol` overrides the policy tolerance for this question only. Pass 0.0 to ask
        the strict question -- "does this touch anything at all" -- which is what a
        caller SEARCHING for a free spot wants: the policy tolerance exists to absorb
        prediction drift, and a placer that aims at it has spent the entire error
        budget before any drift has happened. Reporting still uses the policy value,
        so this changes what a search accepts, not what counts as a finding.

        `pad` grows `rect` on every side before testing, which asks for CLEARANCE
        rather than mere separation: two boxes that miss each other by a hair are
        disjoint, so no measure of interpenetration can see them, and two captions a
        hair apart are one paragraph as far as a reader is concerned. The padding
        affects the question only -- what gets reserved is the caller's real box.
        """
        skip = set(ignore)
        probe_rect = (rect if not pad else
                      Rect(rect.x - pad, rect.y - pad,
                           rect.w + 2.0 * pad, rect.h + 2.0 * pad))
        # WHICH SHAPE THIS QUESTION IS ABOUT. `rect` is still the argument, because
        # every caller but reserve_polyline asks about a box and the padding above
        # is theirs. When a shape is given it is the thing being tested, and `pad`
        # grows it the same way -- on both of its own axes, not on the page's.
        mine = probe_rect if shape is None else (
            shape if not pad else
            OBox(shape.cx, shape.cy, shape.ux, shape.uy,
                 shape.hl + pad, shape.ht + pad))
        out = []
        for r in self._res.values():
            if r.rid in skip:
                continue
            if kind is not None and self.policy.tolerates(kind, r.kind):
                continue
            # parts of one object are not colliding with each other. Only the
            # REPORTING path passes an owner; a placer searching for a free spot
            # calls probe() without one and still sees everything, because a
            # caption should be steered off its own arrow even though sitting on
            # it is not a finding.
            if (owner is not None and related(owner, r.owner)
                    and not (kind in self.policy.text_kinds
                             and r.kind in self.policy.text_kinds)):
                # ...but TWO PIECES OF TEXT belonging to one object still must
                # not overlap: a caption above an arrow and the yield below it
                # are siblings, and them landing on each other is a real defect,
                # not a self-collision. Ownership excuses text over ink, not
                # text over text.
                #
                # COUNTED, NOT SWALLOWED. This is the branch that decides a real
                # geometric overlap is not worth reporting, which is exactly the
                # shape of thing that hides a broken engine: the report goes
                # quiet and the quiet reads as clean. So every exemption is
                # tallied and report() prints the tally. It is not a Finding --
                # nothing here is wrong -- but the number has to be visible,
                # because "0 findings, 0 exemptions" and "0 findings, 214
                # exemptions" mean very different things and would otherwise
                # look identical.
                d = sat_depth(mine, r.shape or r.rect)
                if d > 0.0:
                    self._exempt.append((owner, r.owner, d))
                continue
            d = sat_depth(mine, r.shape or r.rect)
            if d <= 0.0:
                continue
            t = (self.policy.tol_for(kind, r.kind, rect, r.rect)
                 if tol is None else tol)
            if d > t:
                out.append((r, d))
        out.sort(key=lambda t: -t[1])
        return out

    def bbox(self) -> Rect | None:
        """Honest extent of everything DRAWN. This is what the packer reads.

        REGIONS ARE EXCLUDED, because a region is bookkeeping and not ink. It covers
        its own contents by construction, so leaving it out loses nothing -- and
        leaving it IN breaks a golden: a block region spans the full figure width, so
        registering one would widen the page of every scope table whose content
        happened to be narrower than its widest block. The page must be the extent of what is
        on it, not of what was declared about it.

        A panel placed outside the existing content still grows the page, because its
        structures and its frame are each reserved in their own right; only the
        summary rectangle is silent.
        """
        rects = [r.rect for r in self._res.values() if r.kind != "region"]
        if not rects:
            return None
        box = rects[0]
        for r in rects[1:]:
            box = box.union(r)
        return box

    def reservations(self) -> Sequence[Reservation]:
        return tuple(self._res.values())

    def find(self, owner: str) -> Reservation:
        """THE reservation named `owner`. Raises on none and on more than one.

        Not `-> Reservation | None`, and the reason is what the callers are for: a
        grid asking "where did cell 3a end up" so it can draw an arrow to cell 3b
        cannot continue without an answer, and "no such cell" and "two cells with
        that name" are both states it must not paper over. Returning None only moves
        the failure to the next line, as an AttributeError with less to say.

        `find_all` is the question that may legitimately come back empty.
        """
        hits = [r for r in self._res.values() if r.owner == owner]
        if len(hits) == 1:
            return hits[0]
        near = sorted({r.owner for r in self._res.values()})
        raise AllocatorError(
            "find(%r) matched %d reservations; it must match exactly one. "
            "Reserved owners are: %s"
            % (owner, len(hits), ", ".join(near) if near else "(nothing yet)"))

    def region(self, name: str, rect: Rect) -> Reservation:
        """Declare that `name` occupies `rect`, and that its parts live under it.

        A REGION IS NOT A NEW OBJECT, it is an owner-path prefix with a rectangle.
        `related()` already reads dotted paths to decide what belongs to what, so a
        parallel parent/child structure would be a second source of truth for the
        same question. Everything whose owner starts with `name + "."` is inside.

        kind="region" is tolerant against everything by design: a region is drawn
        over -- indeed around -- its own contents, and against a NEIGHBOUR the thing
        that must not overlap is the neighbour's region, which is the one pair left
        strict. Same argument the decoration/structure exemption makes: a property of
        the kinds, not of a relation.
        """
        return self.reserve(rect, name, "region")

    def check_region(self, name: str) -> Severity:
        """Report any part of `name` that pokes outside it. Nothing is moved.

        The point of a region is that its contents are IN it, and the way that stops
        being true is a panel composed at one size and filled with something wider --
        which is invisible on the page until it collides with the next panel. Uses
        Rect.excess_over, which already answers "how far does this stick out".

        NO THRESHOLD. The rule is to ask "compared to what" before fixing a number,
        and there is nothing here to compare to -- so this reports the measurement
        and lets the reader decide. A tolerance can be added when there is a figure
        to calibrate it against.
        """
        outer = self.find(name)
        worst = Severity.OK
        for r in self.find_all(name):
            if r.rid == outer.rid or r.kind == "region":
                continue
            out = r.rect.excess_over(outer.rect)
            if out > 0.0:
                self._record(Severity.WARN, "region-escape",
                             "%r reaches %.1f units outside region %r"
                             % (r.owner, out, name), (r.rid,), out)
                worst = max(worst, Severity.WARN)
        return worst

    def free_near(self, candidates, kind: str, owner: str, moved: list,
                  ignore: Iterable[int] = (), pad: float = 0.0):
        """First candidate whose rect touches nothing. Returns (key, index).

        `candidates` is an iterable of (key, rect, label). The KEY is whatever the
        caller wants back -- an x, a top, an offset -- because the four searches that
        use it each hand a different thing to their own arithmetic, and
        layout_cycle's give_way argues at length for returning an offset rather than
        a rectangle (one piece of arithmetic, evaluated once). Returns (None, None)
        when nothing was free; the caller decides whether that means "keep the
        convention and let the overlap report" or something else.

        WHAT IS SHARED IS THE PROBE AND THE REPORT, NOT THE CANDIDATES. The four
        sequences differ -- mirror-then-outward for a side caption, whole line
        heights for an above/below one, bond-length steps for a feeder, along-the-
        chord for an annotation label -- and each was arrived at separately. Merging
        them would be a layout change dressed as a cleanup; the rule is to measure
        the sequences before deciding they are one.

        tol=0.0, ALWAYS. The policy tolerance exists to absorb prediction error after
        the fact, and a placer aiming at it has spent the whole budget before any
        drift has happened. Left at the policy value, give_way stops one step early
        and settles for a 4.5-unit overlap that the allocator then declines to
        report.

        `moved` IS REQUIRED, and not for symmetry. A search that relocates in silence
        trades a reported overlap for an unexplained gap -- the caption ends up three
        line heights from its arrow and nothing says why, and the findings a silent
        placer stops producing read as an improvement in a count. Making the sink a
        positional argument means "search without a way to report it" cannot be
        written.
        """
        for i, (key, rect, label) in enumerate(candidates):
            if self.probe(rect, kind=kind, ignore=ignore, tol=0.0, pad=pad):
                continue
            if i:
                moved.append((owner,) + (tuple(label) if isinstance(label, tuple)
                                         else (label,)))
            return key, i
        return None, None

    def find_all(self, prefix: str) -> Sequence[Reservation]:
        """Every reservation whose owner path starts with `prefix`.

        Path-aware, not string-aware: "cell:3" does NOT match "cell:3a", because
        the owner grammar is dotted (`cell:3a.bracket`) and a prefix that cuts a
        component in half is a typo rather than a query. `related()` reads the same
        paths for the same reason.
        """
        out = []
        for r in self._res.values():
            if r.owner == prefix:
                out.append(r)
            elif r.owner.startswith(prefix) and r.owner[len(prefix)] in ".[":
                out.append(r)
        return tuple(out)

    # ---- drift ----------------------------------------------------------

    def reconcile(self, rid: int, measured: Rect) -> Severity:
        """Feed back what a primitive actually measured after rendering.

        Within slack: absorbed silently, the stored rect grows so bbox stays
        honest. Beyond slack but hitting nobody: logged, nothing escalates.
        Beyond slack and colliding: that is the one case that produces a bad
        figure, so it is an ERROR.
        """
        res = self._res.get(rid)
        if res is None:
            self._record(Severity.WARN, "unknown-rid",
                         f"reconcile() for rid {rid} which was never reserved")
            return Severity.WARN

        excess = measured.excess_over(res.rect)
        self._drift.append((res.kind, excess))
        grown = res.rect.union(measured)

        if excess <= res.slack:
            self._res[rid] = Reservation(rid, grown, res.owner, res.kind, res.slack)
            return Severity.OK

        hits = self.probe(measured, kind=res.kind, ignore=(rid,))
        self._res[rid] = Reservation(rid, grown, res.owner, res.kind, res.slack)

        if not hits:
            self._record(Severity.WARN, "drift-clear",
                         f"{res.owner!r} ({res.kind}) grew {excess:.1f} units past "
                         f"its {res.slack:.1f} slack but hits nothing",
                         (rid,), excess)
            return Severity.WARN

        other, depth = hits[0]
        self._record(Severity.ERROR, "drift-collision",
                     f"{res.owner!r} ({res.kind}) grew {excess:.1f} units and now "
                     f"overlaps {other.owner!r} by {depth:.1f}",
                     (rid, other.rid), depth)
        return Severity.ERROR

    # ---- emit guard -----------------------------------------------------

    def check_emit(self, rid: int | None, owner: str = "?") -> None:
        """Call from the CDXML writer for every primitive."""
        if rid is None or rid not in self._res:
            self._record(Severity.ERROR, "unreserved-emit",
                         f"{owner!r} emitted with no live reservation "
                         f"(rid={rid!r})")
            if self.policy.enforce_rid:
                raise AllocatorError(f"{owner!r} emitted without a reservation")

    # ---- reporting ------------------------------------------------------

    def note(self, severity: Severity, code: str, message: str,
             rids: tuple[int, ...] = (), amount: float = 0.0) -> None:
        """Report something this layer could not have measured for itself.

        THE ENGINE SOMETIMES KNOWS SOMETHING THE GEOMETRY DOES NOT SAY. A radial
        spoke whose two inward trims have crossed is drawn pointing back at its own
        source, and both endpoints are legal points, the polyline pieces are legal
        boxes in legal places, and a reversed arrow overlaps nothing -- so nothing
        in here can see it. The layout module can, at the one point where the
        intended source->target sense and the drawn sense are both in hand, and
        this is how it says so.

        DELEGATES TO _record RATHER THAN APPENDING TO .findings, which is the whole
        reason it exists. An outside caller appending its own Finding would skip
        the policy.strict raise and make a second way for one to come into being --
        two channels to keep in step by hand, differing exactly when it matters.
        One channel; this is its front door.
        """
        self._record(severity, code, message, rids, amount)

    def _record(self, severity: Severity, code: str, message: str,
                rids: tuple[int, ...] = (), amount: float = 0.0) -> None:
        self.findings.append(Finding(severity, code, message, rids, amount))
        if severity >= Severity.ERROR and self.policy.strict:
            raise AllocatorError(message)

    def worst(self) -> Severity:
        return max((f.severity for f in self.findings), default=Severity.OK)

    def drift_summary(self) -> dict[str, dict[str, float]]:
        """Per-kind observed drift. Use this to tune slack from real figures
        instead of guessing, and as a calibration set for text_metrics."""
        out: dict[str, dict[str, float]] = {}
        for kind, excess in self._drift:
            b = out.setdefault(kind, {"n": 0.0, "max": 0.0, "sum": 0.0})
            b["n"] += 1
            b["sum"] += excess
            b["max"] = max(b["max"], excess)
        for b in out.values():
            b["mean"] = b["sum"] / b["n"] if b["n"] else 0.0
        return out

    def report(self, max_lines: int = 20) -> str:
        """Grouped, capped and worst-first.

        A polyline reservation splits one arc or connector into dozens of pieces,
        so a single bad figure can raise hundreds of near-identical findings --
        778 of them on a four-species cycle drawn with a negative gap. An
        un-collapsed dump is unreadable, and an unreadable report is one nobody
        looks at, which defeats the whole point of having the instrument. Collapse
        by (severity, code, message with piece indices and depths normalised) and
        print the worst groups first.
        """
        lines = [f"reservations: {len(self._res)}   worst: {self.worst().name}"]
        groups: dict[tuple, dict] = {}
        for f in self.findings:
            msg = re.sub(r"\[\d+\]", "[*]", f.message)
            # depth AND its span, so a run of polyline pieces still collapses to one
            # group -- the numbers differ piece by piece, the complaint does not
            msg = re.sub(r"-?[\d.]+ units deep"
                         r"(?: \(-?[\d.]+ wide x -?[\d.]+ tall\))?",
                         "units deep", msg)
            g = groups.setdefault((f.severity, f.code, msg),
                                  {"n": 0, "amount": 0.0})
            g["n"] += 1
            g["amount"] = max(g["amount"], f.amount)
        ranked = sorted(groups.items(),
                        key=lambda kv: (-kv[0][0], -kv[1]["amount"], -kv[1]["n"]))
        for (sev, code, msg), g in ranked[:max_lines]:
            tail = ""
            if g["amount"]:
                tail += " (max %.1f" % g["amount"]
                tail += ", x%d)" % g["n"] if g["n"] > 1 else ")"
            elif g["n"] > 1:
                tail += " (x%d)" % g["n"]
            lines.append(f"  [{sev.name:5}] {code}: {msg}{tail}")
        if len(ranked) > max_lines:
            lines.append("  ... and %d more distinct finding(s); %d in total"
                         % (len(ranked) - max_lines, len(self.findings)))
        if os.environ.get("CDXML_SHADOW_VERBOSE"):
            lines.append("  reservations (x, y, w, h):")
            for r in sorted(self._res.values(), key=lambda r: (r.rect.y, r.rect.x)):
                lines.append("    %-26s %-11s %7.1f %7.1f %7.1f %7.1f"
                             % (r.owner, r.kind, r.rect.x, r.rect.y,
                                r.rect.w, r.rect.h))
        if self._exempt:
            deepest = max(self._exempt, key=lambda t: t[2])
            lines.append("  ownership exemptions: %d overlap(s) not reported "
                         "because both sides belong to one object "
                         "(deepest %.1f, %s over %s)"
                         % (len(self._exempt), deepest[2],
                            deepest[0], deepest[1]))
        ds = self.drift_summary()
        if ds:
            lines.append("  drift by kind (n / mean / max, page units):")
            for kind, b in sorted(ds.items()):
                lines.append(f"    {kind:11} {int(b['n']):4d}  "
                             f"{b['mean']:6.2f}  {b['max']:6.2f}  "
                             f"(slack {self.policy.slack_for(kind):.1f})")
        return "\n".join(lines)
