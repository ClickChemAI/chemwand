"""layout_cycle.py — catalytic cycle / mechanism ring.

The sixth topology, and the one the desktop pipeline had no engine for. Unlike a
scheme it is polar: species sit on a circle and the reaction arrows ARE the ring.

Allocator is load-bearing, as in layout_grid: no parallel `boxes` list, every
primitive reserves, alloc.bbox() is what finalize() gets. Arc arrows reserve the
ink itself via cdxml_core.arc_points -> Allocator.reserve_polyline, the same
primitive connectors and radial spokes use -- the axis-aligned hull of an arc is
almost entirely empty space.

Recipe (fragments already built, same convention as the other engines):

    {
      "structures": { id: (fragment_str, bbox), ... },
      "nodes": [ {"s": id, "label"?: "L~n~Pd(0)"}, ... ],   # in ring order
                 # or, for a species with no structure to draw:
                 #   {"text": "Cu(I)", "id"?: key}
                 # "s" and "text" are mutually exclusive. A text node renders as
                 # that type and nothing else -- no label unless one is given
                 # explicitly, or the name would print twice.
      "steps": [ {"above"?: ["oxidative addition", "+ ArX"]}, ... ],
      "center": ["net:", "ArX + HNR~2~ -> ArNR~2~"]?,
    }

len(steps) may be len(nodes) (closes the ring) or len(nodes) - 1 (open arc).
Species entering or leaving a step are text lines in that step's caption --
"+ ArX", "- HBr" -- which is how most papers draw them. Or they get a FEEDER, the
other literature convention, which a step asks for with

    {"in": [{"text": "PhC≡CH (#2a#)"}], "out": [{"text": "H~2~O"}],
     "aside": ["*t*-BuOK"], "at": 0.5, "feeder_style": "split"}

A feeder is a CIRCLE ARC TANGENT TO THE RING -- see _feeder_circle. An `in` or `out`
entry may DRAW A STRUCTURE instead of writing type, by naming a fragment exactly as a
node does: {"s": "2a", "label"?: "#2a#", "angle"?: -20}. It is drawn upright at the
free end of its arc and never rotated; `angle` turns the direction the box is set out
along, which is the manual nudge and the only one. A structure always sits OUTSIDE the
ring, so a structural by-product is refused on the "split" style, whose by-product arc
is tangent on the inside -- refused rather than quietly restyled, because silently
turning two arcs into one is a different drawing (_check_steps).

THE RING IS PACKED, not divided into n equal slices. 360 degrees is a budget:
every species spends some of it on its own box, every step spends some on its arc
and on whatever feeder hangs off it, and the radius is the smallest one at which
the bill comes to 360 or less (solve_radius, by bisection). Whatever is left over
is shared out evenly, so a step with a feeder gets the room it asked for and the
others keep equal shares of the rest -- the polar form of layout_linear's adaptive
gutter. A ring with no feeder on it spends exactly what the pre-packer formula
gave it and comes out byte for byte identical; that equality is asserted, not
hoped for (legacy_radius).

Text may use the usual markers: ~sub~ ^sup^ *italic* #bold#.

EVERY piece of type here -- species name, label, step caption, feeder label, aside,
centre note -- goes through _fold, which returns the lines to write and their box in
one call. Measuring a block one way and drawing it another lets two 40-character
feeder labels print through each other with a clean layout report, so
there is exactly one text-box path and no way to take half of it.
"""

from __future__ import annotations

import math
from typing import NamedTuple

import address as ADDR
import cdxml_core as C
import graphics as _GX          # straight_head only: two points -> a collinear cubic
import layout_linear as LL
from allocator import Allocator, Policy, Rect, Severity


# --- THE ONE-OPERATION INVARIANT: placement and registration are ONE operation -
#
# Every wrapper below that appends XML takes a REQUIRED `addr_as`, so a new
# placement cannot be written without deciding whether it is addressable. Python
# raises TypeError on a forgotten argument, which is a stronger guard than a test:
# two `place_with` sites with only one of them recording is possible precisely
# when registration is a separate optional line that one path can simply not
# have.
#
# NOT_ADDRESSED is not a way out, it is the other half of the same declaration:
# it carries a reason, and a test enumerates every use.
class Addressed(NamedTuple):
    """This placement IS an object with a house number."""
    aid: str
    kind: str
    handle: object = None
    ink: tuple = ()


class NotAddressed(NamedTuple):
    """This placement is NOT an addressable object, and here is why.

    Captions, labels and the arcs inside a feeder motif are marks ON an object
    rather than objects -- they are reserved in the ledger (so nothing lands on
    them) and have no id, because nothing can point at them.
    """
    why: str


def _ink_rect(pts) -> Rect:
    """The axis-aligned box round a polyline, for an edge's `rect`.

    An edge is a line, so its rect is derived rather than given -- and it is the
    HULL of the ink, which for an arc is mostly air. That is fine here and would
    not be in the ledger: this rect is the record's `rect` field, used for
    host_box (where a caption hangs) and never for collision, which the allocator
    already accounts chunk by chunk through reserve_polyline. Same two-accounts
    split layout_linear's annotation block documents.
    """
    xs = [float(p[0]) for p in pts]
    ys = [float(p[1]) for p in pts]
    return Rect(min(xs), min(ys), max(xs) - min(xs), max(ys) - min(ys))


def _with_handle(addr_as, skey, placed: str, base: int, src=None):
    """Fill in a structure's handle, once the placed string exists.

    WHAT A HANDLE IS: exactly the five things `Anchors.record` needs, and nothing
    else --

        placed      the PLACED fragment string, so the coordinates cannot drift
                    from what was emitted (they ARE what was emitted)
        base        the renumbering offset place_with applied
        smiles      \\ from recipe["anchor_src"], because only the MCP layer knows
        index_map   / them: they are made while the SMILES is converted
        base_id     the pre-place id of the first atom
        s           the structures key, for reports and for joining back

    `src` is `recipe["anchor_src"][skey]` -- the (smiles, index_map, base_id)
    triple the tool layer captured at build time. Absent (an engine-level recipe
    written by hand) leaves smiles None, and address.py then refuses an atom
    address with a message naming what was missing rather than resolving to
    nothing.

    NOTHING IS PARSED HERE. This is references into a dict; `node_points()` is not
    called and no RDKit work happens -- the object address costs a rect the
    packer already had, and the atom parse waits until somebody asks for an
    atom. address.py never looks inside a handle, so this dict can grow a field
    without that module changing.
    """
    if not isinstance(addr_as, Addressed):
        return addr_as                  # a NotAddressed has nothing to fill
    smiles, imap, base_id = (src or (None, None, 1))
    return addr_as._replace(handle={"s": skey, "placed": placed, "base": base,
                                    "smiles": smiles, "index_map": imap,
                                    "base_id": base_id})


def _file(addr: "ADDR.Registry", addr_as, rect: Rect, where: str):
    """Hand one placement to the registry, or record that it is not addressable.

    The single line every placement wrapper calls, sitting BESIDE the append that
    draws the thing -- which is the whole of the one-operation invariant (see the
    note above Addressed). `addr_as` must be one of the two declarations; anything
    else is a programming error and says so loudly rather than defaulting to "not
    addressed", because a silent default is the hole.
    """
    if isinstance(addr_as, Addressed):
        return addr.register(addr_as.aid, addr_as.kind, rect,
                             handle=addr_as.handle, ink=addr_as.ink, where=where)
    if isinstance(addr_as, NotAddressed):
        return None
    raise TypeError(
        "%s: addr_as must be Addressed(...) or NotAddressed(why), got %r. Every "
        "placement declares which it is."
        % (where, addr_as))

block_box = LL.block_box
fold_lines = LL.fold_lines
fold_box = LL.fold_box
ASCENT_EM = getattr(LL, "ASCENT_EM", 0.72)

# Feeder feel. Every number here was MEASURED off two hand-drawn reference figures
# (one "split" feeder pair, one "single" tangent feeder), which is why there are no
# free "tune it till it looks right" knobs left: the shape is a circle tangent to
# the ring, and a circle has nothing to tune. RATIOS of R and DEGREES here, no
# absolute lengths except the label gaps, which are multiples of the bond length.
FEEDER = {
    # r as a fraction of the ring radius, up to the cap in feeder_radius() -- both
    # halves of that rule have a reference drawing behind them.
    "r_in":         0.230,   # split, the reagent arc     (split ref: 21.90 / 95)
    "r_out":        0.200,   # split, the by-product arc  (split ref: 18.99 / 95)
    "r_single":     0.185,   # single                     (single ref: 17.57 / 95)
    # how far each arc carries, in degrees of ITS OWN circle
    "sweep_in":      76.0,
    "sweep_out":    124.0,
    "sweep_single": 150.0,
    # free end -> label box centre, x C.BOND_LENGTH. A tail wants more room than a
    # head: at a head the arrowhead triangle has already filled the gap. Measured
    # 1.55/1.48 at the tails and 0.69/1.15 at the heads.
    "tail_gap":     1.50,
    "head_gap":     1.00,
    "aside_gap":    1.60,    # tangency point -> aside box, along the ring normal
    "fan":          0.12,    # extra `at` per additional in/out pair, dimensionless
}

# The clamp window for feeder_radius(), in bond lengths, and the ratio it was
# measured against.
FEEDER_R_NOMINAL = 0.20
FEEDER_R_BOND = (1.0, 1.8)

# How far a box may GIVE WAY along its own tangent when its measured spot is taken,
# and how finely that is searched. The budget itself is DERIVED FROM THE BOX -- see
# give_way_cap, which is where the two numbers below are spent.
#
# NOT from the reference drawings, deliberately kept out of FEEDER for that reason:
# every number in there is a measurement off those drawings, and these are a repair
# budget, measured against failures instead.
#
# Past the cap the box stays exactly where the measurement puts it and the step
# warns. Sliding further would trade a visible overlap for a label quietly attached
# to the wrong part of the figure, which is the worse failure: an overlap is obvious
# to whoever opens the file, a mis-bound label is not.
GIVE_WAY_BOND = 2.0        # floor, in bond lengths: the smallest useful budget
GIVE_WAY_STEP_BOND = 0.25  # how finely the offset is searched, in bond lengths
GIVE_WAY_SPAN = 1.0        # ...and the budget for a box bigger than that floor


def feeder_radius(R: float, bond: float, ratio: float = FEEDER_R_NOMINAL) -> float:
    """Radius of one feeder arc. A FUNCTION OF ITS OWN ON PURPOSE.

    SETTLED by two reference drawings, one either side of the knee. A single drawing
    cannot decide it when its ring is 6.6 bond lengths across, because there
    `0.20 * R` and `1.3 * bond` are the same number: scale by R and a big ring grows
    fat feeders with it, peg to the bond length and a feeder stays a constant-sized
    decoration that only tracks the type size, and such a figure cannot tell the two
    apart.

    The second drawing (a big hand-drawn ring, R = 16.18 bond) does, and the answer
    is that BOTH readings are right over their own range. `0.20 * R` meets the
    1.8-bond cap at R = 9 bond, and the two references straddle it: at 6.6 bond the
    measured r is 1.32 = 0.20 * 6.6, uncapped, so a small ring really does scale;
    at 16.18 bond it is 1.80, not the 3.24 that scaling would give, so a big one
    really does stop. So: proportional, then capped -- which is what a clamp IS.

    Its own function, and the only place the choice lives, so that the answer stays
    a one-function question if another reference drawing ever bears on it.
    The window is rescaled by `ratio` / FEEDER_R_NOMINAL so that capping a crowded
    figure shrinks all three arcs together instead of collapsing them onto one
    radius: the caps are 2.07 / 1.80 / 1.665 bond for in / out / single. The LOWER
    bound has no drawing behind it -- it would take a ring under about 5 bond lengths
    to reach -- and stays as a guard, not a finding.
    """
    lo, hi = FEEDER_R_BOND
    k = ratio / FEEDER_R_NOMINAL
    return min(max(ratio * R, lo * bond * k), hi * bond * k)


class FeederMotif(NamedTuple):
    """One step's in/out/aside: what to draw, what it claimed, how wide it sits.

    `up` and `down` are how far the motif reaches either side of its own tangency
    point, in degrees at the ring centre, measured along the ring's direction of
    travel; `angular_extent` is their sum, the total width, which is what the
    report shows. TWO numbers rather than one because a feeder is not symmetric
    about the point where it meets the ring -- a reagent arc reaches back upstream,
    a by-product arc forward downstream -- and a gap wide enough for the total is
    still too narrow when the whole of it has to fit on one side.

    Boxes that landed in the free middle of the cycle are left out of both: they cost
    disk, not circumference (split's by-product label belongs inside the ring),
    and being near the centre they subtend an enormous angle there, so
    charging them to the angle budget would let them dominate the one thing they do
    not sit on. What they do have to clear is the centre note, which is _disk_clash's
    job and is asked of every box.

    `angular_extent` is None when the motif holds a box whose circumscribed circle
    COVERS THE RING CENTRE. Such a box subtends every angle there is, so it has no
    angular width to report -- 360 would be a saturation, not a measurement, and
    reporting it next to real extents invites it to be compared with them. A box
    like that goes under the radial check (_disk_clash) instead, and
    `extent_defined` says which regime this motif fell into.

    The BUDGET is unaffected: `up` and `down` still carry the 180 that such a box
    forces, because an unpayable demand is what makes the solver grow the ring, and
    that is the right answer to a label lying across the middle of the cycle. What
    changes here is only what gets published.

    The packer consumes all of this (see solve_radius).
    """
    prims: list
    reservations: list
    angular_extent: float | None
    up: float
    down: float
    extent_defined: bool = True


ORIGIN = (0.0, 0.0)


def _ring_frame(t0: float, sweep: float, R: float, at: float, centre=ORIGIN):
    """(P, N) at arc parameter `at` in [0,1]: the point on the ring, and the unit
    outward normal there.

    `centre` DEFAULTS TO THE ORIGIN, and every geometry helper in this module that
    takes a centre does the same. That is not politeness: a single-ring figure
    passes nothing, so the arithmetic is bit-for-bit the single-ring arithmetic,
    which a byte-identity regression check verifies. The normal is
    centre-independent.
    """
    th = math.radians(t0 + sweep * at)
    c, s = math.cos(th), math.sin(th)
    return (centre[0] + c * R, centre[1] + s * R), (c, s)


def _feeder_circle(t0: float, sweep: float, R: float, at: float, r: float,
                   outside: bool, centre=ORIGIN):
    """The circle of radius `r` TANGENT TO THE RING at arc parameter `at`.

    `at` runs along the arc AS DRAWN -- `t0` and `sweep` arrive already trimmed clear
    of both species -- so its default 0.5 is the midpoint of the visible arrow, which
    is the point the eye reads as the middle of the step. That is the house
    convention and it is deliberate: do not "fix" this to measure from one species
    centre to the next, which would put 0.5 somewhere that is not the middle of any
    ink. Off-centre values are for the exceptions (pulling an in/out pair onto two
    tangency points, steering a motif out of a crowded quarter).

    This one line -- put the feeder's centre at distance R +/- r along the radius
    through the tangency point -- is the whole model. Everything a feeder does
    follows from it: the arc touches the ring at exactly one point, so it reads as
    joining the cycle rather than being laid across it, and it needs no free
    parameters, because a circle through a point with a given tangent has none.

    Returns (Cf, alpha_T, dirn):

      Cf       the feeder circle's centre.
      alpha_T  the tangency point's angle AS SEEN FROM Cf. An externally tangent
               circle touches on its far side, hence the 180.
      dirn     the sign that walks alpha away from the tangency point on the
               UPSTREAM side of the ring's travel for an external circle, and on
               the DOWNSTREAM side for an internal one. The two tangency modes
               reverse the rolling sense -- a point at alpha_T + d sits at
               T - r*d*(travel direction) outside and at T + r*d*(travel) inside --
               so ONE sign serves both an arriving feeder (external, wants its free
               end upstream) and a departing one (internal, wants it downstream).
               Measured that way round on the split reference drawing; flip this
               to flip both at once.
    """
    th = math.radians(t0 + sweep * at)
    d = R + (r if outside else -r)
    Cf = (centre[0] + math.cos(th) * d, centre[1] + math.sin(th) * d)
    alpha_T = math.degrees(th) + (180.0 if outside else 0.0)
    return Cf, alpha_T, (1.0 if sweep >= 0 else -1.0)


def _circle_pt(Cf, r: float, alpha: float):
    a = math.radians(alpha)
    return Cf[0] + math.cos(a) * r, Cf[1] + math.sin(a) * r


def _circle_tangent(alpha: float, away: float):
    """Unit tangent at angle `alpha`, pointing AWAY from the arc.

    `away` is -sign(arc sweep) at the tail end (back the way the arc came) and
    +sign(arc sweep) at the head (on the way it was going). Labels sit along this
    line, so a feeder reads as one continuous flow of ink -- word, arc, ring --
    rather than an arc with a caption parked near it. Checked against the single
    reference drawing, whose tail label sits 0.1 deg off this direction.
    """
    a = math.radians(alpha)
    return -math.sin(a) * away, math.cos(a) * away


def give_way(box_at, alloc: Allocator, allow, cap: float, step: float,
             pad: float = 0.0, kind: str = "caption", owner: str = "?",
             moved: list | None = None):
    """Smallest offset in [0, cap] at which `box_at(e)` is clear of everything already
    reserved -- or None, meaning it never got clear.

    `pad` is how much CLEARANCE is wanted, not just separation. Two labels whose
    boxes miss each other by 1.7 units are disjoint, so nothing that measures
    interpenetration can see them -- and yet they print as one block, because a text
    box is smaller than the ink it stands for. The caller passes the block leading,
    so two separate labels never sit closer than two lines of one label would.

    Takes a FUNCTION OF THE OFFSET and returns an OFFSET, which is the whole point.
    The caller builds the box it commits by calling the same `box_at` with the answer,
    so the rectangle that was tested, the rectangle that gets reserved and the
    rectangle the type is written into are one piece of arithmetic evaluated once.
    A routine that returned a rectangle would invite the caller to re-derive the
    text's anchor from it, and "measured here, drawn there" is exactly the fork the
    single folding path closes.

    Stepped rather than solved. What is being avoided is a specific box, at a
    specific place, and the clearance test is the allocator's own occupancy map --
    not a distance with a closed form. Stepping also makes the answer a round number
    of bond lengths, so a figure does not shift by a hair when a neighbour moves by
    one.

    CLEAR MEANS TOUCHING NOTHING, tol=0.0, not "within the reporting tolerance". Left
    at the policy value this stops one step early and settles for a 4.5-unit
    overlap that the allocator then declines to report -- exactly the kind of silent
    failure the allocator exists to prevent. The tolerance is there to absorb
    prediction error after the fact; aiming at it spends it in advance. Kind pairs
    that overlap BY DESIGN are still exempt (a caption may cross the arrow it
    annotates), because that exemption is about meaning, not about slack.

    `kind` is what is being moved, and it decides only which of those by-design
    exemptions apply: a caption may cross an arrow, a structure may not, so a
    structural feeder searching for room has to clear the ring arcs as well (its own
    arc is on `allow`, being tangent to them by construction).

    `moved` is where a SUCCESSFUL displacement is recorded, and it exists for the
    reason free_near's docstring gives: reporting only FAILURES (`no_room`) makes a
    label that quietly slid 28 units to get clear look identical to one that never
    had to move -- the figure right and the account of it missing. Callers with
    nothing to report it to pass nothing and get a throwaway list -- the radius
    trial run does, since its whole allocator is discarded a few lines later.
    """
    def cands():
        e = 0.0
        while e <= cap + 1e-9:
            yield (e, box_at(e),
                   ("out along its own tangent", e,
                    "%d of %d steps" % (round(e / step), round(cap / step))))
            e += step
    # A caller with nowhere to report to still must not silently drop the
    # record: `moved=None` becomes a throwaway list HERE, visibly, rather than
    # free_near quietly accepting None.
    got, _i = alloc.free_near(cands(), kind, owner,
                              moved if moved is not None else [],
                              ignore=allow, pad=pad)
    return got


def give_way_cap(w: float, h: float, ux: float, uy: float) -> float:
    """How far a w x h box may give way along (ux, uy), in page units.

    A FUNCTION OF THE BOX, and that is the whole point. The floor is the flat
    GIVE_WAY_BOND, calibrated on feeder LABELS: the worst measured clash (a folded
    by-product label against its own folded reagent label) clears at 1.40 bond
    lengths and the next at 1.00, so 2.0 bond -- 28.8 units -- covers the known
    cases with 40% to spare. A STRUCTURAL feeder makes that number useless rather
    than merely tight: a molecule's box is 60 to 100 units across, so a structure
    asked to give way would exhaust the whole budget inside its own width, warn, and
    stay put. A repair that cannot repair anything is worse than none, because it
    reports itself as having tried.

    So the budget is the box's own EXTENT along the direction it is being pushed
    (2 * _support, i.e. exactly the length the box occupies on that line), with the
    label constant kept as a floor. The window is pinned at both ends, as fold_bond's
    is, and it is narrow:

      floor    2.0 bond = 28.8 units, the label measurement above. A small box must
               still get the budget that was measured for it; `span` overtakes the
               floor at a box 28.8 units long, which is about a two-line caption.
               MEASURED, on the structural clashes to hand: two reagents fanned onto
               one step need 36.0 units of tangent to come apart (phenylacetylene and
               benzoic acid), and a structure at each end of a "single" arc needs
               39.6. Both are past the flat 28.8, and with the span term switched off
               both fail exactly as described above -- the first leaves an 11.2-unit
               molecule-on-molecule overlap and a "nowhere to go" warning about a
               budget that was never the right size for the thing it would not move.
      ceiling  span exactly, ONE box length, geometric rather than measured -- which
               is the honest way round, because what this end guards against is not an
               overlap that can be measured but a species that has come loose from its
               arc. Displaced by its own extent a box finishes edge-to-edge with the
               footprint it was measured at: it still touches the spot the geometry
               chose for it, so a reader still binds it to that arc. One box length
               further and nothing does.
               It is also, and not by coincidence, the WORST CASE THAT CAN ARISE
               between two boxes of the same size: separating two boxes along u costs
               their two extents along u less the offset already between their
               centres, so two equal boxes exactly on top of each other cost span and
               nothing costs more. Hence the measured offsets sit where they do -- 96%
               and 99% of budget for two synthetic test pairs, 70% and 38% for the
               real molecules above.
               A box clashing with something LARGER than itself can still run out, and
               then it warns and stays put, which is the visible failure rather than
               the silent one.

    So this is not the label constant with a coat of paint: 4.2 bond lengths of
    budget for a 60-unit fragment and 6.9 for a 100-unit one, and still exactly 2.0
    for every text label -- which is why no label-only figure moves.
    """
    return max(GIVE_WAY_BOND * C.BOND_LENGTH,
               GIVE_WAY_SPAN * 2.0 * _support(w, h, ux, uy))


def _box_disk(rc, centre=ORIGIN):
    """A reserved box as (how far its centre is from the ring centre, its own
    circumscribed radius).

    EVERY packing question in this module is asked of that pair rather than of the
    four corners, and the reason is that the packer measures a feeder before it knows
    where on the circle the feeder will sit. The boxes are axis-aligned, so their
    corners are not: rotate a motif and the corner angles change, which made a
    demand measured in one frame wrong in the frame that got drawn -- measured, on a
    three-species ring, as a downstream reach of 19.6 degrees that became 42.7 once
    the ring was repacked, and a species overlap that followed from exactly that gap.
    A disk containing the box contains it at every rotation, so the answer stops
    depending on the frame. Conservative for a wide label, which is the safe
    direction, and it is the convention the rest of the module already uses: a
    species is its `reach` and the centre note is its `note_reach`, both this same
    half-diagonal.
    """
    return (math.hypot(rc.cx - centre[0], rc.cy - centre[1]),
            math.hypot(rc.w / 2.0, rc.h / 2.0))


def _reach_span(rects, disk: float, th_t: float, travel: float,
                centre=ORIGIN, interior_ok: bool = False):
    """(up, down, defined) for a motif's reserved boxes -- see FeederMotif.

    `defined` is False when some box covers the ring centre, which makes the motif's
    angular width a saturation rather than a measurement. The two angles still come
    back saturated, because the packer's job is to grow the ring until such a box
    stops lying across the middle; `defined` is what stops that 360 being published
    as if it were a number anyone could compare.

    Angles are measured RELATIVE TO th_t, the motif's own tangency angle, and signed
    by `travel` so positive is downstream. Relative because that is what the packer
    needs (how much arc must lie either side of the crossing point), and because it
    sidesteps the 0/360 wrap that breaks every max-minus-min on a circle.

    `disk` is the radius of the FREE INTERIOR: the ring's radius less the largest
    species reach. A box wholly inside it cannot touch a species -- every one of
    those lies within that reach of the circle -- so it spends no circumference and
    is skipped here; it still answers to _disk_clash for the centre note. Anything
    reaching into the species annulus pays for its angle like everything else there,
    however far inside its centre happens to sit: a by-product label 9 bond lengths
    wide can be centred well inside the ring and still stick out past the species on
    either side of it.
    """
    up = down = 0.0
    defined = True
    for rc in rects:
        d, rho = _box_disk(rc, centre)
        if d + rho < disk:
            continue
        # 180 when the box covers the ring centre: it then subtends every angle
        # there is, which is unpayable and makes the solver grow the ring, which
        # is the right answer to a caption lying across the middle of the cycle.
        #
        # With `interior_ok` the middle of the ring is a place things are SUPPOSED
        # to be -- a nested ring, a species at the centre -- so a box lying across
        # it is no longer a reason to grow the circle. The flag arrives per RING
        # from cfg, never as a change to any global default: see
        # `interior_allocatable` in CYCLE_DEFAULTS for why that distinction is the
        # whole point.
        if d <= rho:
            if not interior_ok:
                defined = False
            else:
                # It spends no circumference: it is not in the annulus the species
                # occupy, it is in the middle, which is now allocatable. The centre
                # note still answers to _disk_clash.
                continue
        half = 180.0 if d <= rho else math.degrees(math.asin(rho / d))
        rel = ((math.degrees(math.atan2(rc.cy - centre[1], rc.cx - centre[0]))
                - th_t + 540.0) % 360.0 - 180.0) * travel
        up, down = max(up, half - rel), max(down, half + rel)
    return up, down, defined


CYCLE_DEFAULTS = {
    "node_gap": 40.0,      # clear space demanded between adjacent node boxes
    "arc_margin": 12.0,    # how far an arc backs off from a node before starting
    "cap_gap": 14.0,       # step caption clear of the ring
    "label_gap": 6.0,      # node label clear of its structure
    "cap_size": 9,
    "num_size": 11,
    "center_size": 10,
    "text_size": 11,       # a {"text": ...} node's own type size...
    "text_bold": True,     # ...and weight. Matches num_size above, so swapping a
                           # faked {"s": "[Cu]", "label": "Cu(I)"} for a real text
                           # node changes nothing but the loss of the stray atom.
    "min_radius": 95.0,
    "min_sweep": 12.0,     # never emit an arc shorter than this, it reads as a dot
    "direction": "cw",     # "cw" | "ccw"
    "label_pos": "auto",   # "below" (house style) | "outward" | "auto"
    "feeder_style": "split",   # every step's default; a step may override it
    "center_gap": 26.0,    # centre note must stay this clear of the species
    "fold_bond": 10.0,     # widest a text line may render, IN BOND LENGTHS: see _fold
    "arc_style": {},       # default arc_arrow kwargs for every step
    "start_angle": 270.0,  # 270 puts the first species at the TOP (CDXML y grows down)
    "arc_samples": 48,
    "margin": 30.0,
    # Whether the ring's middle may hold things. FALSE by default, so no existing
    # figure changes and NOTHING GLOBAL MOVES: this is a cycle cfg key, not an
    # allocator.py default -- see the note in _graph_from_tables on why that
    # distinction is the whole point. Set per figure by a caller who knows the
    # middle is in use; desugar also turns it on per RING when it can see a nested
    # ring or a pinned interior node.
    "interior_allocatable": False,
    # The radius of a `crossing: "bridge"` hump.
    #
    # COMPARED TO WHAT: the arrowhead already on the same line. HeadSize 1000 is
    # 10 layout units (every head number in cdxml_core is /100), so a hump of
    # radius 5 spans the same 10 units the head does -- the two marks on one line
    # read at one scale, and the number is not a fresh guess.
    "bridge_hop_r": 5.0,
    # The colortable INDEX a muted path is drawn in, or
    # None for "nothing is muted". An INDEX, not an [r, g, b] -- rgb has to be
    # appended to the document's colortable before it has one, and the tool layer
    # is what holds the header (cdxml_core.color_index says so in its own error).
    # None by default, so every call site can pass it through unbranched and an
    # unmuted figure emits exactly the bytes it always did.
    "muted_color": None,
}

_CW_SIGN: float | None = None


def cw_sign() -> float:
    """+1 if a POSITIVE sweep advances clockwise on screen, else -1.

    arc_arrow's docstring says a positive sweep is visually counter-clockwise,
    but arc_points advances the parameter through `cy + sin(t) * ry` in a
    coordinate system where y grows DOWN, which makes it clockwise. Rather than
    pick a side and risk drawing every cycle backwards, sample the primitive and
    read the answer off the geometry. Cached; call reset_cw_sign() in a test.
    """
    global _CW_SIGN
    if _CW_SIGN is None:
        try:
            pts = C.arc_points(C.arc_arrow(1, 0.0, 0.0, 10.0, 10.0, 0.0, 90.0), n=8)
            # start at angle 0 -> right of centre. A clockwise quarter turn ends
            # up BELOW centre, i.e. at a larger y once y grows down.
            _CW_SIGN = 1.0 if (pts and pts[-1][1] > pts[0][1]) else -1.0
        except Exception:
            _CW_SIGN = 1.0
    return _CW_SIGN


def reset_cw_sign() -> None:
    global _CW_SIGN
    _CW_SIGN = None


def _cfg(user: dict | None) -> dict:
    out = dict(LL.CFG)
    out.update(CYCLE_DEFAULTS)
    out.update(user or {})
    return out


# ---- what each level of a recipe accepts ---------------------------------
#
# WHITELISTS, not blacklists. A recipe is hand-written, or written from a docstring
# by something that half-remembers an older schema, and a key this engine does not
# read is nearly always a typo or a leftover -- `{"bananna": 7}` or a key that
# has since gone. Ignoring it silently is the worst outcome available: the figure
# comes out subtly not-what-was-asked-for and there is nothing to point at. So an
# unknown key raises and names its own path, and adding a key to the schema means
# adding it here on purpose.
NODE_KEYS = frozenset({"s", "text", "label", "id", "as",
                       # the graph form's node keys
                       "kind", "at", "pairs",
                       # `frame` -- the square-bracket form only
                       "frame",
                       # `content` -- a slot holding a SEQUENCE
                       "content"})
# `content`: what ONE ITEM of that sequence takes. A species written or drawn,
# optionally with its own name under it -- deliberately NOT the node vocabulary:
# `at` / `kind` / `pairs` / `frame` / `as` are properties of the SLOT, and an
# item is a thing inside the slot, not a slot of its own.
CONTENT_KEYS = frozenset({"s", "text", "label"})
# ...and the connectives an item may be instead of a dict. ONE, because "+" is the
# only one the schema shows and the only one the known uses need (`2a + 2b`).
# A second one is a line here plus a size, not a design.
CONNECTIVES = ("+",)
# Edges and rings, the NEW input form. The engine validates these because it is
# what consumes them; the MCP layer's own whitelists sit above.
EDGE_KEYS = frozenset({"as", "from", "to", "via", "ring", "name", "feeders",
                       "style", "reverse",
                       # `muted` is NOT here: it lives on the path row -- see
                       # PATH_KEYS.
                       "bulge", "crossing"})
RING_KEYS = frozenset({"as", "nodes", "centre", "radius", "start_angle",
                       "direction", "label",
                       # PLACED BY TOUCHING ANOTHER RING. `centre` says where;
                       # this says "wherever that one is, touching it, in this
                       # direction" -- which is how published dual cycles are
                       # mostly drawn (tangent, with the SET / PT / HAT badge ON
                       # the touch point).
                       # The two are alternatives, not a pair; see the check below.
                       "tangent_to", "at"})
# The FOURTH table. A path is a NAMED SET OF EDGES: the host of a group-level name
# (a "TM" over a pair of cross-ring edges, an "ROH" over a pair of edges from one
# source, "path A" / "path B") and the unit of muting. Nodes, edges and rings are
# the first three tables; this is the fourth, and it exists because both of those
# are properties of a SET and neither has anywhere else to live -- writing a group
# name on one member means a second member can write it again and the label draws
# twice, which no rule could catch.
PATH_KEYS = frozenset({"as", "edges", "name", "muted"})
# An edge's `crossing`. THREE values, and the third one names its counterpart:
# `{"crossed": "<edge as>"}` says "the crossing you see between me and that edge is
# deliberate -- draw nothing, and do not report it". Published cycles need it: a
# deliberate coupling mark crosses plainly, with no hump, and without this the only
# way to draw one would be to leave a collision in the report.
CROSSING = ("plain", "bridge", "crossed")
FEEDERS_KEYS = frozenset({"in", "out", "style", "at", "aside"})
VIA = ("arc", "line")           # "polyline" is in the schema but not built
STEP_KEYS = frozenset({"above", "below", "inside", "reverse", "style", "at",
                       "in", "out", "aside", "feeder_style", "as"})
# A feeder entry is a species like a node is, so it takes a node's vocabulary plus
# the two placement keys only a feeder has. NOT "id": a node reads one (a text node
# needs a name for its reservations) and a feeder does not -- its boxes are named by
# side and index, "in[0]", which is what the recipe calls it too. The MCP layer takes
# an "id" and spends it on the structures dict before the engine ever sees the entry.
FEEDER_KEYS = frozenset({"text", "s", "label", "at", "angle", "as"})

# TWO NAMESPACES, AND THEY DO NOT MEET. The principle, stated once:
#
#   `s` / the structures key answers  "does this molecule need building again?"
#                                     -> CHEMICAL IDENTITY. Shared on purpose:
#                                     three nodes may name one key, and the tool
#                                     documents a feeder and a node sharing one.
#   `as` answers                      "which of the five ascorbic acids on this
#                                     page do you mean?" -> INSTANCE.
#
# Construction may be shared; an address may not. So `as` is a SEPARATE optional
# key rather than a reinterpretation of `id` or `s`, and neither namespace gives
# anything up to the other.
#
# Omitted -> the positional fallback (`node#0`, `edge#2`, `edge#0.in[1]`). A
# recipe with no `as` at all keeps the positional names and draws the same
# bytes.
_AS_RESERVED = "#"


# A node's `frame`. Only `shape: "bracket"` is built.
FRAME_KEYS = frozenset({"shape", "lip"})
FRAME_SHAPES = ("bracket",)
# The lip is the little return at each end of the stroke. 60 is what ChemDraw
# writes: every hand-drawn bracket measured on published figures carries
# LipSize="60".
BRACKET_LIP = 60
# How far outside the node's own box each stroke stands. COMPARED TO WHAT: the
# same clearance a connection keeps off a box -- address.py's ATTACH_TRIM,
# 12.0 -- so a bracket and an arrowhead stop the same distance from the species.
BRACKET_PAD = 12.0


def _frame(value, where: str):
    """A node's `frame`, validated. -> a dict or None.

    Square brackets round a species -- e.g. `[Ar-Pd(IV)]`, or a bracketed
    intermediate. The spelling is copied from hand-drawn ChemDraw files:
    `<graphic GraphicType="Bracket" BracketType="Square" LipSize="60"/>`, one
    VERTICAL STROKE per side (x1 == x2), the two written in opposite directions.

    The schema also lists `fill` and `label` on a frame. Neither is built, and
    writing either one raises rather than being dropped -- a box the caller asked
    for and did not get is a silent failure.
    """
    if value is None:
        return None
    if not isinstance(value, dict):
        raise ValueError('%s: a frame is a dict, e.g. {"shape": "bracket"}, got %s'
                         % (where, type(value).__name__))
    check_keys(value, FRAME_KEYS, where)
    shape = value.get("shape")
    if shape not in FRAME_SHAPES:
        raise ValueError(
            "%s.shape: one of %s, got %r. Only a bracket frame is implemented; "
            "`fill` and `label` are not."
            % (where, " / ".join(map(repr, FRAME_SHAPES)), shape))
    lip = value.get("lip", BRACKET_LIP)
    if isinstance(lip, bool) or not isinstance(lip, (int, float)) or lip <= 0:
        raise ValueError("%s.lip: a positive number, got %r" % (where, lip))
    return {"shape": shape, "lip": int(lip)}


def bracket_pair(gid: int, rect, lip: int) -> list:
    """The two `<graphic>` strokes that bracket `rect`. Hand-drawn spelling.

    A bracket is one vertical stroke, so a pair is two objects. THE DIRECTION IS
    WHAT PICKS THE SIDE: in a hand-drawn ChemDraw file the left stroke runs
    bottom -> top and the right one top -> bottom, and nothing else in the element
    says left or right. A test file writing both the same way was used to check
    this, and confirmed that the two graphics alone draw, with no
    <bracketedgroup>.
    """
    x0, x1 = rect.x - BRACKET_PAD, rect.right + BRACKET_PAD
    yt, yb = rect.y - BRACKET_PAD, rect.bottom + BRACKET_PAD
    return [
        '<graphic id="%d" BoundingBox="%.2f %.2f %.2f %.2f" GraphicType="Bracket"'
        ' BracketType="Square" LipSize="%d"/>' % (gid, x0, yb, x0, yt, lip),
        '<graphic id="%d" BoundingBox="%.2f %.2f %.2f %.2f" GraphicType="Bracket"'
        ' BracketType="Square" LipSize="%d"/>' % (gid + 1, x1, yt, x1, yb, lip),
    ]


def _crossing(value, where: str):
    """An edge's `crossing`, validated. -> (value, counterpart or None).

    Absent means "plain". `bridge` finds its own crossings and draws a hump at each,
    so it names nothing. `crossed` draws NOTHING, so it has to name the edge it
    deliberately crosses -- otherwise "exempt this crossing" would have to mean
    "exempt this edge against everything", which is exactly the blanket exemption
    this scoping avoids.
    """
    if isinstance(value, dict):
        unknown = set(value) - {"crossed"}
        if unknown:
            raise ValueError(
                "%s.crossing: %s is not a crossing key. The dict form is "
                '{"crossed": "<the edge this one deliberately crosses>"}.'
                % (where, ", ".join(sorted(map(repr, unknown)))))
        tgt = str(value.get("crossed") or "")
        if not tgt:
            raise ValueError(
                '%s.crossing: {"crossed": ...} needs the name of the edge this one '
                "deliberately crosses." % where)
        return "crossed", tgt
    v = "plain" if value is None else str(value)
    if v == "crossed":
        # THE VALUE EXISTS, the bare string just cannot carry the counterpart.
        # Saying so beats "one of 'plain' / 'bridge' / 'crossed'", which would send
        # the caller looking for a typo in a value they spelled correctly.
        raise ValueError(
            '%s.crossing: "crossed" has to name the edge it crosses, because the '
            "exemption is scoped to that one edge and to nothing else. Write "
            '{"crossed": "<edge as>"}.' % where)
    if v not in CROSSING:
        raise ValueError(
            "%s.crossing: one of %s, got %r. 'bridge' draws the hop over whatever "
            "this edge crosses and 'crossed' draws nothing and silences that one "
            "crossing; both are declared and never inferred, because the same "
            "crossing can be a deliberate coupling marker or a typesetting accident "
            "and the geometry is identical."
            % (where, " / ".join(map(repr, CROSSING)), value))
    return v, None


def _bulge(value, where: str):
    """An edge's `bulge`: how far the edge bows, as a fraction of its own
    length. None means a straight line. Refuses 0 rather than treating it as
    straight -- writing `"bulge": 0` means the caller wanted a bow and will get a
    line, and `cdxml_core.arc_chord` raises on it anyway."""
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("%s.bulge: a number (fraction of the chord), got %r"
                         % (where, value))
    v = float(value)
    if v == 0.0:
        raise ValueError(
            "%s.bulge: 0 is not a bow. Leave 'bulge' out for a straight edge."
            % where)
    if abs(v) > 2.0:
        # A sagitta of 2x the chord is already a >270-degree arc; past that the
        # bow wraps around and reads as a different edge. Named rather than
        # clamped: clamping would draw something the caller did not ask for.
        raise ValueError(
            "%s.bulge: %r bows by %.1fx the edge's own length, which wraps past "
            "three quarters of a circle and stops reading as the same edge. Parallel "
            "paths in published cycles typically use 0.2-0.4." % (where, value, abs(v)))
    return v


def circle_through(p, m, q):
    """(cx, cy, r) of the circle through three points, or None if they are collinear.

    A marker's `pairs` are consumed by drawing each pair as ONE continuous curve
    through the marker, and three points define exactly one circle -- so "the curve
    that leaves P, passes through the marker and arrives at Q" is not a shape
    anyone has to choose. That is what makes a wrong pairing VISIBLE: pair P with a
    different Q and the circle is a different circle.
    """
    (px, py), (mx, my), (qx, qy) = p, m, q
    d = 2.0 * (px * (my - qy) + mx * (qy - py) + qx * (py - my))
    if abs(d) < 1e-9:
        return None                      # collinear: no finite circle
    p2, m2, q2 = px * px + py * py, mx * mx + my * my, qx * qx + qy * qy
    cx = (p2 * (my - qy) + m2 * (qy - py) + q2 * (py - my)) / d
    cy = (p2 * (qx - mx) + m2 * (px - qx) + q2 * (mx - px)) / d
    return cx, cy, math.hypot(px - cx, py - cy)


def arc_through(p, m, q):
    """(cx, cy, r, t0, sweep) for the arc P -> Q that passes through M.

    The sweep is the one that reaches M before Q, which is what picks the correct
    half of the circle -- the other half is the same circle and the wrong gesture.
    Returns None when the three points are collinear.
    """
    c = circle_through(p, m, q)
    if c is None:
        return None
    cx, cy, r = c
    ang = lambda z: math.degrees(math.atan2(z[1] - cy, z[0] - cx))
    t0 = ang(p)
    dm = (ang(m) - t0) % 360.0
    dq = (ang(q) - t0) % 360.0
    sweep = dq if dm < dq else dq - 360.0
    return cx, cy, r, t0, sweep


def _polyline_at(pts, s: float):
    """The point `s` units along a polyline, measured by arc length."""
    if s <= 0.0:
        return pts[0]
    run = 0.0
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        d = math.hypot(x1 - x0, y1 - y0)
        if run + d >= s:
            f = (s - run) / d if d else 0.0
            return (x0 + (x1 - x0) * f, y0 + (y1 - y0) * f)
        run += d
    return pts[-1]


def _polyline_len(pts) -> float:
    return sum(math.hypot(b[0] - a[0], b[1] - a[1]) for a, b in zip(pts, pts[1:]))


def _seg_cross(p, q, r, t):
    """Where segments p-q and r-t cross, or None. Proper crossings only.

    Endpoints excluded (the open interval), because two edges that MEET at a shared
    node touch there by construction and a hop drawn on a node would be nonsense.
    """
    (px, py), (qx, qy), (rx, ry), (tx, ty) = p, q, r, t
    dx1, dy1, dx2, dy2 = qx - px, qy - py, tx - rx, ty - ry
    den = dx1 * dy2 - dy1 * dx2
    if abs(den) < 1e-12:
        return None                       # parallel or degenerate
    a = ((rx - px) * dy2 - (ry - py) * dx2) / den
    b = ((rx - px) * dy1 - (ry - py) * dx1) / den
    if not (1e-9 < a < 1 - 1e-9 and 1e-9 < b < 1 - 1e-9):
        return None
    return (px + dx1 * a, py + dy1 * a)


def bridge_hops(pts, others, hop_r: float):
    """Where a polyline has to hop over the others, as arc-length positions.

    An edge's `crossing: "bridge"`. WHERE the hop goes is pure geometry and is
    computed; WHETHER to draw one must be declared and is never inferred -- the
    same crossing can be a deliberate coupling mark or a typesetting accident:
    opposite meaning, identical geometry.

    Hops closer together than two radii are merged, and hops within one radius of
    either end are dropped: a hump half-off the end of the line reads as a kink.

    `others` is {edge name: polyline}; the third return value names the edges that
    were actually crossed, because a hop TOUCHES the line it hops -- at a right-
    angle crossing the arc's apex lands on it however the bow is chosen -- and the
    ledger needs that contact exempted against those edges and no others.
    """
    hits = []
    run = 0.0
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        seg = math.hypot(x1 - x0, y1 - y0)
        for name, op in others.items():
            for oa, ob in zip(op, op[1:]):
                c = _seg_cross((x0, y0), (x1, y1), oa, ob)
                if c is not None:
                    hits.append((run + math.hypot(c[0] - x0, c[1] - y0), name))
        run += seg
    total = run
    hits = sorted(h for h in hits if hop_r < h[0] < total - hop_r)
    merged, crossed = [], []
    for h, name in hits:
        if name not in crossed:
            crossed.append(name)
        if merged and h - merged[-1] < 2.0 * hop_r:
            continue
        merged.append(h)
    return merged, total, crossed


def _end_addr(v, where: str, by_as: dict, gnodes, eid_of: dict):
    """One edge endpoint, validated. -> a node id (str) or a parametric dict.

    `from` / `to` take address-layer addresses, and specifically `{"of": "e3",
    "at": 0.5}` lets one edge terminate on another edge.

    A NODE endpoint has a node index and everything downstream can use it. A
    PARAMETRIC one does not: it is a point on another edge's ink, so it belongs to
    no ring, has no reach and cannot be a ring arc. Every pass that indexes by node
    skips it, and `edge_shape` calls the edge OFF-RING -- which is exactly the
    rule "at least one end belongs to no ring".
    """
    if isinstance(v, str) or v is None:
        aid = str(v)
        if aid not in by_as:
            raise ValueError(
                "%s: %r is not a node in this figure. The nodes are %s. To "
                "terminate on an EDGE instead, name it parametrically: "
                '{"of": "<edge as>", "at": 0.5}.'
                % (where, aid, ", ".join(repr(g["as"]) for g in gnodes)))
        return aid
    if not isinstance(v, dict):
        raise ValueError(
            "%s: an endpoint is a node id, or {\"of\": <edge>, \"at\": 0..1} for a "
            "point on an edge; got %s." % (where, type(v).__name__))
    unknown = set(v) - {"of", "at"}
    if unknown:
        # `atom` is a legal address and the address layer resolves it, but
        # nothing has wired it as an EDGE ENDPOINT -- say that, rather than let
        # a key the caller believes in fall through.
        raise ValueError(
            "%s: %s is not usable on an endpoint. An endpoint takes 'of' + 'at' "
            "(a point on another edge's ink); an atom refinement on an endpoint is "
            "not implemented." % (where, ", ".join(sorted(map(repr, unknown)))))
    if "at" not in v:
        raise ValueError(
            '%s: {"of": ...} on an endpoint needs "at" (0..1 along that edge\'s '
            "ink). Without it, name the node directly." % where)
    tgt = str(v.get("of") or "")
    if tgt not in eid_of:
        raise ValueError(
            "%s.of: %r is not an edge in this figure. The edges are %s."
            % (where, tgt, ", ".join(map(repr, eid_of)) or "none"))
    t = v["at"]
    if isinstance(t, bool) or not isinstance(t, (int, float)):
        raise ValueError("%s.at: a number 0..1 along the edge's ink, got %r"
                         % (where, t))
    if not (0.0 <= float(t) <= 1.0):
        raise ValueError("%s.at: 0..1 along the edge's ink, got %g" % (where, t))
    return {"of": tgt, "at": float(t)}


def _as_id(value, where: str, fallback: str) -> str:
    """Validate an `as` and return it, or the positional fallback when absent.

    Refuses '#' so the two families of name stay disjoint: the positional ids are
    built with it, and a caller writing `"as": "node#1"` on a different node would
    otherwise collide with node 1's own fallback and get a duplicate-id error about
    a name they never chose.
    """
    if value is None:
        return fallback
    if not isinstance(value, str) or not value.strip():
        raise ValueError(
            "%s: 'as' is the name you want to address this object by, a non-empty "
            "string -- got %r. Leave it out and the engine names it %r."
            % (where, value, fallback))
    aid = value.strip()
    if _AS_RESERVED in aid:
        raise ValueError(
            "%s: 'as' may not contain %r -- that character belongs to the names "
            "the engine generates when you do not give one (%r), and reusing it "
            "would collide with an object you did not name."
            % (where, _AS_RESERVED, fallback))
    return aid

_BEZIER_GONE = (
    "feeders are drawn as an arc tangent to the ring, which has no shape to "
    "tune; remove this key")

# Keys that MEANT something once, and what to reach for instead. Worth the words:
# whoever hits one of these has a working older-schema recipe in front of them, and
# "unknown key" alone would send them looking for a typo that is not there.
_V1_GONE = {
    "join": (
        'both of its values were tangent arcs -- "cross" and "tangent" were never '
        'two different shapes -- so the choice moved up a level and changed meaning. '
        'Put "feeder_style": "split" (two arcs sharing one tangency point, the '
        'default) or "single" (one arc straddling it) on the STEP'),
    "head": (
        'which end carries the arrowhead now follows from the step\'s '
        '"feeder_style": "split" draws the reagent arc headless and gives the '
        'by-product arc the only Full head, "single" heads the product end. Both '
        'were measured off the reference drawings, so there is nothing left to '
        'choose per entry. To restyle the RING arc instead, use the step\'s "style"'),
    "cross_depth": _BEZIER_GONE,
    "arm": _BEZIER_GONE,
    "bulge": _BEZIER_GONE,
    "eps": _BEZIER_GONE,
    "label_gap": (
        _BEZIER_GONE + ". The node-label gap that shares the name is a cfg key, "
        "not a recipe key"),
}

_NODE_GONE = dict(_V1_GONE, smiles=(
    "a node names a fragment this engine was handed, so use 's' -- a key into "
    "recipe['structures']. The MCP layer is what turns a SMILES into one"))

# Same message, same reason: a feeder entry draws a structure, and it names that
# structure exactly as a node does.
_FEEDER_GONE = dict(_V1_GONE, smiles=_NODE_GONE["smiles"].replace("a node", "a feeder"))


def check_keys(d: dict, allowed: frozenset, where: str,
               gone: dict = _V1_GONE) -> None:
    """Raise on the first key `where` does not read. Path-first, as every other
    message in this engine is, so the complaint is greppable against the recipe.

    Public because the MCP wrapper has its own node vocabulary ("smiles" where the
    engine takes "s") and rebuilds each node before this engine sees it -- so an
    unknown key there would be dropped on the floor before any check here could
    catch it. One checker, one message shape, both entry points.
    """
    for k in d:
        if k in allowed:
            continue
        hint = gone.get(k)
        raise ValueError(
            "%s: unknown key %r%s. This level takes: %s"
            % (where, k, " -- %s" % hint if hint else "",
               ", ".join(sorted(allowed))))


def _check_steps(steps, cfg: dict) -> None:
    """Every step and every feeder entry, before ANYTHING is measured or drawn.

    Up front rather than as each step is reached: a step past the end of the ring
    (an open arc draws n-1 of them) is never visited, and a recipe whose last step
    is malformed should still say so rather than quietly dropping it.

    Takes `cfg` for one question only -- whether a step that named no feeder_style
    is going to be "split" -- and that question has to be asked here because the
    answer decides whether a STRUCTURAL by-product is drawable at all.
    """
    for i, st in enumerate(steps):
        if not st:
            continue
        if not isinstance(st, dict):
            raise ValueError("steps[%d]: a step is a dict, got %s"
                             % (i, type(st).__name__))
        check_keys(st, STEP_KEYS, "steps[%d]" % i, _V1_GONE)
        style = str(st.get("feeder_style") or cfg["feeder_style"])
        for tag in ("in", "out"):
            for k, it in enumerate(st.get(tag) or []):
                where = "steps[%d].%s[%d]" % (i, tag, k)
                if not isinstance(it, dict):
                    raise ValueError(
                        "%s: a feeder entry is a dict like {\"text\": \"H~2~O\"}, "
                        "got %s" % (where, type(it).__name__))
                check_keys(it, FEEDER_KEYS, where, _FEEDER_GONE)
                has_t, has_s = bool(it.get("text")), it.get("s") is not None
                if has_t and has_s:
                    raise ValueError(
                        "%s: 'text' and 's' are mutually exclusive -- a feeder "
                        "entry is either a structure or a piece of type, not both"
                        % where)
                if it.get("label") and not has_s:
                    # 'label' belongs to a DRAWN species -- it is the caption under
                    # the molecule. On a text entry there is nothing to caption and
                    # nothing would be drawn with it, and a key that goes nowhere is
                    # the failure this whitelist exists to prevent.
                    raise ValueError(
                        "%s: 'label' captions a DRAWN structure, so it needs an 's'. "
                        "A text entry is already its own label -- put the extra words "
                        "in the 'text' (\"PhC≡CH (#2a#)\"), or give the entry an 's' "
                        "and let the label sit under the molecule." % where)
                if has_s and tag == "out" and style == "split":
                    raise ValueError(
                        "%s: a structural by-product cannot use feeder_style "
                        "'split'. Split's by-product arc is tangent INSIDE the ring "
                        "-- that is what makes the pair read as one S through it -- "
                        "so its free end, and this structure with it, would land in "
                        "the middle of the cycle, across the centre note and the "
                        "species opposite. Every reference drawing puts a structure "
                        "OUTSIDE the ring. Use \"feeder_style\": \"single\", whose "
                        "one arc keeps both ends outside; or write this by-product "
                        "as 'text', which is small enough to live in there. A "
                        "structural 'in' is fine either way -- that arc is "
                        "externally tangent." % where)


# Only these reach arc_arrow. Anything else in a step's "style" is dropped here --
# but it does not GET here: preflight.preflight_cycle refuses an unknown key in
# `arc_style`, `steps[].style` and `edges[].style` before this module is entered,
# so the filter below is the second line of defence rather than the policy.
#
# THAT IS DELIBERATE. Dropping the key so that "a recipe written against a newer
# primitive still draws" would be forward compatibility that buys a recipe nothing,
# because the key it carried would be silently ignored and the figure would come
# back with the engine default and no note. `preflight.ARC_KEYS` restates this tuple
# (that module imports no engine) and a test asserts the two are equal.
_ARC_KEYS = ("head", "tail", "head_type", "head_size", "head_width",
             "head_center", "line")


def _arc_kw(cap: dict, cfg: dict) -> dict:
    """arc_arrow kwargs for one step: engine default, then the step's override.

    head="Full" + tail="Full" is the reversible/equilibrium step; head=None is a
    plain headless arc; "HalfLeft"/"HalfRight" give the single-electron half head;
    line="dashed" for a dashed shaft.
    """
    out = {k: v for k, v in (cfg.get("arc_style") or {}).items() if k in _ARC_KEYS}
    for k, v in (cap.get("style") or {}).items():
        if k in _ARC_KEYS:
            out[k] = v
    return out


def _support(w: float, h: float, ux: float, uy: float) -> float:
    """Half-extent of a w x h box along the unit direction (ux, uy).

    Exactly how far the box reaches that way, so a caption can be pushed just
    clear of the ring instead of by a guessed constant.
    """
    return (abs(ux) * w + abs(uy) * h) / 2.0


def _fold(lines, size: float, cfg: dict):
    """(folded lines, (w, h)) for one text block. EVERY piece of type this engine
    writes comes through here, and takes BOTH halves of the answer from one call.

    That is the whole design: without folding, a 41-character reagent label is
    reserved as one 171-unit line and drawn as one 171-unit line, so the engine is
    self-consistent and still wrong -- two such labels on one step overlap by 57
    units and the allocator's depth metric, which reports the SHORTER way out of an
    overlap, scores the 1.8-unit vertical sliver and stays quiet. Folding is what
    makes a long label a tall box instead of a flat one, which is both the readable
    drawing and the honest reservation.

    cfg["fold_bond"] is the cap, IN BOND LENGTHS, and it is a FIXED constant on
    purpose -- not the room this particular label happens to have. Deriving it from
    available space closes a loop: fold width -> box height -> feeder reach -> ring
    radius -> available space -> fold width, and the radius solver is a bisection
    over exactly that quantity, so the loop would make its demand non-monotone and
    the bisection would chatter between a folded and an unfolded answer.

    WHY 10 (= 144 units at the ACS 14.4-unit bond). The window is pinned at both
    ends by measurement, and it is narrow:

      floor   the widest line in any of the reference regression figures is 8.41
              bond lengths (a centre note, "ArX + HNR~2~ -> ArNR~2~ + HX",
              121.0 units). Fold below that and a figure that was already right
              changes, so the cap must clear 8.41 -- 10 clears it by 19%.
      ceiling the narrower of the two overlapping labels in the failing case
              (a by-product label, 10.70 bond / 154.1 units) has to FOLD, or the
              overlap survives. So the cap must sit under 10.70.

    (8.41, 10.70) admits 10 and very little else; 12 would leave that case's
    11.90-bond reagent label on one line and change nothing about the overlap. Two
    independent sanity checks land in the same place: 144 units is about 32
    characters at cap_size 9, the usual line length inside published artwork, and
    it is 3/4 of the diameter of the smallest ring this engine will draw
    (min_radius 95), so a caption line never grows wider than its own cycle. The
    hand-drawn references are no help here -- their longest caption is
    "PhC≡CH (2a)" at 3.7 bond -- which is why the regression figures are the
    evidence.
    """
    return fold_box(lines, size, cfg["fold_bond"] * C.BOND_LENGTH)


def _content_layout(items, S: dict, cfg: dict, where: str):
    """A node's `content`: ONE SLOT HOLDING A SEQUENCE, laid out inline.

    Published cycles want this in two forms and it is the same mechanism both
    times -- a metal centre written `L~n~Pd^II^` with a real allyl chain drawn
    beside it, and a pair of products `2a` + `2b`. Neither is "a node that is
    secretly two nodes": both are ONE species slot on the ring, spending one share
    of the angle budget, with more than one thing in it.

    LAYOUT_LINEAR'S ROW, reused rather than reinvented (by design -- it is the
    same shape): items advance left to right, a species advances by its own
    width and a "+" by its width plus `plus_gap` each side. That inner loop is
    `layout_linear.render_linear`'s, and `cfg` already carries its numbers --
    `_cfg` starts from `LL.CFG` -- so a "+" between two species in a ring slot is
    the same mark at the same spacing as the "+" between two products in a scheme.

    A TEXT ITEM IS MEASURED AS A TEXT NODE IS, not as linear's text cell: `_fold`
    at `text_size`, bold. Otherwise `{"text": "Cu(I)"}` would be one size as a node
    and another size inside a `content`, which is the kind of difference nobody
    finds by reading. What IS taken from linear is the horizontal breathing room,
    `text_cell_pad`, whose stated reason -- "so type does not touch a bond" -- is
    exactly the situation here: `L~n~Pd^II^` sits immediately left of a drawn chain.

    THE COMMON AXIS IS THE SPECIES AXIS, not the box axis. Every species is centred
    on one horizontal line and a label hangs BELOW its own item, so a pair of
    products sits on one line with their numbers under them rather than being pushed
    up and down by which of them happens to be named.

    Returns (parts, (w, h), axis); `axis` is measured DOWN from the box top.
    """
    if not isinstance(items, list) or not items:
        raise ValueError(
            "%s: 'content' is a NON-EMPTY LIST -- the sequence in one slot, e.g. "
            "[{\"text\": \"L~n~Pd^II^\"}, {\"s\": \"allyl\"}] or "
            "[{\"s\": \"2a\"}, \"+\", {\"s\": \"2b\"}]. Got %s."
            % (where, "an empty list" if isinstance(items, list)
               else type(items).__name__))
    parts = []
    for j, it in enumerate(items):
        w2 = "%s[%d]" % (where, j)
        if isinstance(it, str):
            # A BARE STRING IS A CONNECTIVE. Not a shorthand for {"text": ...}:
            # a species and a "+" ride at different sizes and take different
            # spacing, so guessing which one a loose string meant would be
            # guessing about the drawing.
            if it not in CONNECTIVES:
                raise ValueError(
                    "%s: a bare string in 'content' is a CONNECTIVE, and %s is the only one "
                    "available. Got %r -- a species written out goes in a dict: {\"text\": "
                    "%r}."
                    % (w2, " / ".join(repr(c) for c in CONNECTIVES), it, it))
            pw, ph = block_box([it], cfg["plus_size"])
            parts.append({"what": "plus", "text": it, "lines": [it],
                          "size": cfg["plus_size"], "s": None, "frag": None,
                          "w": pw, "h": ph,
                          "lbl": None, "llines": None, "lw": 0.0, "lh": 0.0})
            continue
        if not isinstance(it, dict):
            raise ValueError(
                "%s: an item of 'content' is a dict -- {\"text\": ...} or "
                "{\"s\": ...} -- or the bare string \"+\". Got %s."
                % (w2, type(it).__name__))
        check_keys(it, CONTENT_KEYS, w2, _NODE_GONE)
        has_t, has_s = bool(it.get("text")), it.get("s") is not None
        if has_t and has_s:
            raise ValueError(
                "%s: 'text' and 's' are mutually exclusive -- an item is either a "
                "structure or a piece of type, not both. Two things in one slot is "
                "what the SEQUENCE is for: write them as two items." % w2)
        if not (has_t or has_s):
            raise ValueError(
                "%s: needs either 's' (a key into structures) or 'text'" % w2)
        if has_s:
            bb = S[it["s"]][1]
            p = {"what": "s", "s": it["s"], "frag": S[it["s"]][0],
                 "text": None, "lines": None,
                 "size": None, "w": bb[2] - bb[0], "h": bb[3] - bb[1]}
        else:
            lines, (tw, th) = _fold([it["text"]], cfg["text_size"], cfg)
            p = {"what": "text", "s": None, "frag": None,
                 "text": it["text"], "lines": lines,
                 "size": cfg["text_size"],
                 "w": tw + 2 * cfg["text_cell_pad"], "h": th}
        # EXPLICIT ONLY, on both kinds of item -- a structure NODE defaults its
        # label to its structures key, and doing that here would print "allyl"
        # under the chain hanging off a written metal centre, which is the first
        # of the two uses. The pair of products names both of its items on
        # purpose; the metal centre names neither.
        p["lbl"] = it.get("label") or None
        if p["lbl"]:
            p["llines"], (p["lw"], p["lh"]) = _fold([p["lbl"]], cfg["num_size"], cfg)
        else:
            p["llines"], p["lw"], p["lh"] = None, 0.0, 0.0
        parts.append(p)

    x = 0.0
    for p in parts:
        if p["what"] == "plus":
            x += cfg["plus_gap"]
            p["dx"] = x
            x += p["w"] + cfg["plus_gap"]
        else:
            p["dx"] = x
            x += p["w"]
    # A NAME WIDER THAN THE THING IT NAMES still has to be inside the box the ring
    # reserved for this slot, or an arc backs off from the species and lands on the
    # number. The ADVANCE stays the species width -- items must stay adjacent -- so
    # the overhang is taken out of the box's edges instead.
    left = min(p["dx"] + p["w"] / 2.0 - max(p["w"], p["lw"]) / 2.0 for p in parts)
    right = max(p["dx"] + p["w"] / 2.0 + max(p["w"], p["lw"]) / 2.0 for p in parts)
    for p in parts:
        p["dx"] -= left
    axis = max(p["h"] for p in parts) / 2.0
    below = max((p["h"] / 2.0 + cfg["label_gap"] + p["lh"]) if p["lbl"]
                else p["h"] / 2.0 for p in parts)
    return parts, (right - left, axis + below), axis


def _inferred_kind(m: dict) -> str:
    """`structure` or `text` for a node that did not say.

    A CONTENT SLOT IS A STRUCTURE SLOT IF ANYTHING IN IT IS DRAWN. That is the
    whole rule, and it is the one both known uses want: a written metal centre
    with a real chain beside it behaves like a structure -- an edge may name an
    atom of the chain -- while `2a` + "+" + `2b` written out has no atoms and is
    type. `marker` is never inferred, here or below: it is the caller's own word,
    because a SET box and a written species are both type on the page.
    """
    if m["content"] is not None:
        return "structure" if any(p["what"] == "s" for p in m["content"]) else "text"
    return "text" if m["frag"] is None else "structure"


def measure(recipe: dict, cfg: dict):
    """One dict per node: its box, its label's box, and how far it reaches.

    Two kinds of node, mutually exclusive. A STRUCTURE node names an entry in
    recipe["structures"] through "s" and is measured by that fragment's bbox. A
    TEXT node carries "text" instead and IS that text -- nothing is drawn but the
    type, so its box is the text box.

    Nothing downstream learns the difference. The ring only ever pushes
    rectangles around: the radius solver, the arc trimming and the caption placer
    all work off "w"/"h" and the "reach" derived from them, never off atoms. So a
    text node that reports an honest box needs no special case anywhere else.
    """
    S = recipe.get("structures") or {}
    meas = []
    for i, nd in enumerate(recipe["nodes"]):
        if not isinstance(nd, dict):
            raise ValueError("nodes[%d]: a node is a dict like {\"s\": \"cat\"} or "
                             "{\"text\": \"Cu(I)\"}, got %s"
                             % (i, type(nd).__name__))
        check_keys(nd, NODE_KEYS, "nodes[%d]" % i, _NODE_GONE)
        has_t, has_s = bool(nd.get("text")), nd.get("s") is not None
        has_c = nd.get("content") is not None
        if has_c and (has_t or has_s):
            raise ValueError(
                "nodes[%d]: 'content' IS the slot's contents, so it replaces 's' / "
                "'text' rather than sitting beside them. Put the %s in the list: "
                "{\"content\": [{\"%s\": ...}, ...]}."
                % (i, "structure" if has_s else "type", "s" if has_s else "text"))
        if has_t and has_s:
            raise ValueError(
                "nodes[%d]: 'text' and 's' are mutually exclusive -- a node is "
                "either a structure or a piece of type, not both" % i)
        if not (has_t or has_s or has_c):
            raise ValueError(
                "nodes[%d]: needs either 's' (a key into structures) or 'text', or "
                "'content' for a slot holding several things" % i)

        sid = str(nd["s"] if has_s else nd.get("id", "n%d" % (i + 1)))
        # `tlines` / `llines` are the lines the drawing loop WRITES, and the boxes
        # beside them are what it reserves -- one _fold call each, so a folded
        # species name cannot be measured on two lines and drawn on one.
        tlines, content, axis = None, None, 0.0
        if has_c:
            # `w`/`h` are the COMPOSITE box, which is the only thing the rest
            # of the engine ever sees: the radius solver, the arc trimming and the
            # caption placer all work off "w"/"h" and the "reach" derived from
            # them. A slot with three things in it spends one species' share of the
            # angle budget, sized by that box -- which is what makes it one node.
            txt, frag = None, None
            content, (w, h), axis = _content_layout(
                nd["content"], S, cfg, "nodes[%d].content" % i)
            # No default, for the text node's reason: the items carry their own
            # names, and the slot's name is a separate claim the caller has to make.
            lbl = nd.get("label") or None
        elif has_t:
            txt, frag = nd["text"], None
            tlines, (w, h) = _fold([txt], cfg["text_size"], cfg)
            # The type already reads as the species' name, so defaulting a label
            # on would print it twice -- the exact duplication this node kind
            # exists to remove. Only an explicit label earns a second line.
            lbl = nd.get("label") or None
        else:
            txt, frag = None, S[nd["s"]][0]
            bb = S[nd["s"]][1]
            w, h = bb[2] - bb[0], bb[3] - bb[1]
            # An explicit "" means no label; only a missing one falls back to
            # the structures key, which a recipe may use as the compound number.
            lbl = nd.get("label")
            if lbl is None:
                lbl = str(nd["s"])
        if lbl:
            llines, (lw, lh) = _fold([lbl], cfg["num_size"], cfg)
        else:
            llines, lw, lh = None, 0.0, 0.0
        meas.append({"nd": nd, "sid": sid, "frag": frag, "text": txt,
                     "tlines": tlines, "w": w, "h": h,
                     "content": content, "axis": axis,
                     "lbl": lbl, "llines": llines, "lw": lw, "lh": lh,
                     "reach": math.hypot(w / 2.0, h / 2.0)})
    return meas


# ---- the angle budget ----------------------------------------------------
#
# 360 degrees is a budget. Node i spends the angle its own box subtends at the
# centre; the step from node i to node i+1 spends what its arc and its feeder need.
# Both shrink as R grows -- they are lengths divided by a radius -- so there is a
# smallest R at which the bill comes to 360, and solve_radius finds it.

RING_AS = "ring#0"          # the one ring nodes + steps can express; see desugar()


def _inside_ring(pt, gr) -> bool:
    """Is `pt` inside ring `gr`'s free middle?

    "Inside" means inside the circle, not merely inside the bounding box -- a ring
    is round. The radius may still be "auto" at this point (the packer has not run),
    so `min_radius` stands in: it is the smallest circle this engine will draw, so a
    point inside THAT is inside any radius the solver can choose. Conservative in
    the safe direction -- it can miss a nesting the solver later creates, never
    invent one.
    """
    r = gr.get("radius")
    r = float(r) if isinstance(r, (int, float)) else CYCLE_DEFAULTS["min_radius"]
    return math.hypot(pt[0] - gr["centre"][0], pt[1] - gr["centre"][1]) < r


def _direction(v, where: str) -> str:
    """'cw' or 'ccw'. Validated where it enters, not where it is used."""
    if v not in ("cw", "ccw"):
        raise ValueError("%s: direction is 'cw' or 'ccw', got %r" % (where, v))
    return v


def _xy(v, where: str) -> tuple:
    """A [x, y] pair from a recipe. LAYOUT frame -- negatives are normal."""
    if not isinstance(v, (list, tuple)) or len(v) != 2:
        raise ValueError("%s: coordinates are [x, y], got %r" % (where, v))
    try:
        return (float(v[0]), float(v[1]))
    except (TypeError, ValueError):
        raise ValueError("%s: coordinates hold a non-number: %r" % (where, v))


def edge_shape(e: dict, rings, node_rings: dict) -> str:
    """Which of the edge shapes this edge IS. Derived, never declared.

    These fall out of the data, and a test asserts the kind the engine judged, not
    the picture. So the judgement lives in one function:

        arc        both ends on ONE ring and ADJACENT in its cyclic order
        chord      both ends on one ring, NOT adjacent
        crossring  the ends share NO ring, but each is on one
        offring    at least one end is on no ring at all

    The order matters. A SHARED node (one id listed by two rings) belongs to
    several rings, so "same ring" is an INTERSECTION question, not an equality one --
    and that is exactly what makes a shared-node double ring work without a new
    shape: an edge from the shared node to a member of either ring still finds a
    common ring and is that ring's arc.

    ADJACENCY IS CYCLIC AND DIRECTED-BLIND: on a 4-ring, 3 -> 0 is adjacent. An
    open ring's missing gap is simply an edge nobody wrote, not a special case.
    A 2-ring has no non-adjacent pair, hence the `n > 2` guard.
    """
    a, b = e["_from"], e["_to"]
    if a is None or b is None:
        # A PARAMETRIC END belongs to no ring, so the off-ring rule ("at least
        # one end belongs to no ring") is satisfied by construction and there is
        # nothing to intersect.
        return "offring"
    ra, rb = set(node_rings.get(a, ())), set(node_rings.get(b, ()))
    if not ra or not rb:
        return "offring"
    both = ra & rb
    if not both:
        return "crossring"
    by_as = {g["as"]: g for g in rings}
    # A shared node can put two rings in `both`; prefer one where the pair is
    # ADJACENT, because that is the ring whose arc this edge is.
    for r_as in sorted(both):
        gr = by_as[r_as]
        pos, n = gr["_pos"], len(gr["_idx"])
        pa, pb = pos.get(a), pos.get(b)
        if pa is None or pb is None:
            continue
        if n > 2 and ((pa + 1) % n == pb or (pb + 1) % n == pa):
            e["_ring_pick"] = r_as
            return "arc"
        if n == 2:
            e["_ring_pick"] = r_as
            return "arc"
    e["_ring_pick"] = sorted(both)[0]
    return "chord"


def _graph_from_tables(recipe: dict, cfg: dict, gnodes, by_as) -> dict:
    """The NEW input form: `edges` / `rings` are given, so validate and resolve.

    The old form's desugar INVENTS the tables; this one is handed them and has to
    turn addresses into the indices the geometry uses. Everything else -- the shape
    classification, the payload view, the ring row -- is shared, so there is one
    definition of what an edge means regardless of how the caller wrote it.
    """
    # ---- RINGS: ANY n --------------------------------------------------------
    #
    # As many as the caller writes -- not just two; published cycles include
    # three-ring figures -- each with its own centre, radius, start_angle and
    # direction.
    rings_in = list(recipe.get("rings") or [])
    if not rings_in:
        # No `rings`: the obvious one -- every node the caller did NOT pin with an
        # `at`, in the order they were written. That is exactly the ring the old
        # form has implicitly, so "edges but no rings" is not a special mode.
        _free = [g["as"] for g in gnodes if g["at"] is None]
        if len(_free) < 2:
            raise ValueError(
                "with no 'rings' the ring is every node that has no 'at', and %d "
                "of %d node(s) qualify -- a cycle needs at least two species on "
                "the circle. Either drop an 'at', or write the ring you mean in "
                "'rings'." % (len(_free), len(gnodes)))
        rings_in = [{"nodes": _free}]

    grings = []
    ring_of: dict = {}          # node index -> [ring `as`, ...] in listing order
    for ri, raw in enumerate(rings_in):
        where = "rings[%d]" % ri
        if not isinstance(raw, dict):
            raise ValueError("%s: a ring is a dict with 'nodes', got %s"
                             % (where, type(raw).__name__))
        rg = dict(raw)
        check_keys(rg, RING_KEYS, where)
        member_as = [str(v) for v in (rg.get("nodes") or [])]
        if not member_as:
            raise ValueError(
                "%s: a ring needs 'nodes' -- the addresses of its members in "
                "cyclic order." % where)
        seen: set = set()
        for a in member_as:
            if a not in by_as:
                raise ValueError(
                    "%s.nodes: %r is not a node in this figure. The nodes are %s."
                    % (where, a, ", ".join(repr(g["as"]) for g in gnodes)))
            if a in seen:
                raise ValueError(
                    "%s.nodes: %r appears twice. A ring lists each member once, in cyclic "
                    "order. (Two RINGS may share a member -- that is a shared-node double "
                    "ring -- but one ring may not list the same node twice.)" % (where, a))
            seen.add(a)
        for g in gnodes:
            if g["at"] is not None and g["as"] in seen:
                raise ValueError(
                    "nodes: %r has an explicit 'at' AND is listed in %s.nodes. A node with "
                    "'at' is placed by its coordinates, one without is placed by its ring -- "
                    "it cannot be both. Drop one."
                    % (g["as"], where))
        r_as = _as_id(rg.get("as"), where, "ring#%d" % ri)
        idx = [by_as[a] for a in member_as]
        for ix in idx:
            ring_of.setdefault(ix, []).append(r_as)
        lab = rg.get("label")
        if lab is not None:
            if not isinstance(lab, dict):
                raise ValueError("%s.label: a dict with 'text' and 'position', "
                                 "got %s" % (where, type(lab).__name__))
            _p2 = lab.get("position", "center")
            if _p2 not in ("center", "above", "below"):
                raise ValueError(
                    "%s.label.position: one of 'center' / 'above' / 'below', got "
                    "%r." % (where, _p2))
        # Two rules for `tangent_to`, both refused rather than resolved: `centre`
        # and `tangent_to` are two answers to one question, and a direction is not
        # guessable -- 0 (to the right) would be a silent choice of layout, which
        # is the thing this key exists to state.
        if rg.get("tangent_to") is not None:
            if rg.get("centre") is not None:
                raise ValueError(
                    "%s: 'centre' and 'tangent_to' both say where this ring goes. "
                    "Give one -- 'centre' to place it yourself, 'tangent_to' to "
                    "hang it off ring %r." % (where, rg["tangent_to"]))
            if rg.get("at") is None:
                raise ValueError(
                    "%s: 'tangent_to' needs 'at' -- which way from ring %r this "
                    "one sits, in degrees (the CDXML frame has y DOWNWARD, so 0 is "
                    "right, 90 is BELOW, 180 left, 270 above)." 
                    % (where, rg["tangent_to"]))
        grings.append({
            "as": r_as,
            "nodes": member_as,
            # A NON-ZERO CENTRE is allowed: with several rings, each one needs
            # somewhere to go.
            "centre": (_xy(rg["centre"], where + ".centre")
                       if rg.get("centre") is not None else ORIGIN),
            # WHETHER THE CALLER CHOSE IT, which ORIGIN cannot say on its own:
            # a ring that was given no centre and one placed at (0, 0) by hand
            # arrive here identical, and the packer has to tell them apart to
            # know whether it may move this ring. See _spread_rings.
            "_centre_given": rg.get("centre") is not None,
            # Resolved once the radii are known (nothing can place a
            # tangent ring before it knows how big it is), down where the rings are
            # spread -- see the `tangent_to` block there.
            "tangent_to": rg.get("tangent_to"),
            "tangent_at": (float(rg["at"]) if rg.get("at") is not None else None),
            # WHETHER THE CALLER PINNED THE ROTATION, for the same reason
            # `_centre_given` exists: a tangent ring with no stated start_angle is
            # rotated to meet its neighbour ARC-FIRST, and one that states it is
            # left alone. A default that arrives here as a number is not a choice.
            "_start_given": rg.get("start_angle") is not None,
            "radius": rg.get("radius", recipe.get("radius")
                             or cfg.get("radius") or "auto"),
            "start_angle": float(rg.get("start_angle", cfg["start_angle"])),
            "direction": _direction(rg.get("direction", cfg["direction"]),
                                    where + ".direction"),
            "label": lab,
            "_idx": idx,
            "_pos": {ix: pp for pp, ix in enumerate(idx)},
        })

    # A node's placing ring is the FIRST that lists it (a shared node); every
    # ring that lists it is kept, because "is this shared" is a question the shape
    # classifier asks.
    for g in gnodes:
        rs = tuple(ring_of.get(g["_i"], ()))
        g["_rings"] = rs
        g["_ring"] = rs[0] if rs else None
        # A node has to be placeable. There are exactly two ways -- a ring places
        # it, or its own `at` does -- so a node with neither has no position at all.
        # Refused here, where the caller can see which node it is: without this the
        # drawing loop would crash on a None ring, which says nothing useful.
        if not rs and g["at"] is None:
            raise ValueError(
                "nodes: %r is in no ring and has no 'at', so nothing says where it goes. "
                "A node is placed by a ring that lists it, or by its own coordinates. "
                "Add it to a ring's 'nodes', give it an 'at', or drop it." % g["as"])

    rg = dict(rings_in[0])
    ring_idx = grings[0]["_idx"]
    ring_pos = grings[0]["_pos"]
    n_ring = len(ring_idx)

    gedges = []
    # EVERY EDGE'S NAME, BEFORE THE LOOP. An edge may terminate on an edge written
    # later in the table, so the set of edge ids has to exist before the first
    # endpoint is validated. Same ids the loop then assigns -- one rule.
    _eids = [_as_id((ed or {}).get("as"), "edges[%d]" % k, "edge#%d" % k)
             for k, ed in enumerate(recipe["edges"])]
    _eid_of = {}
    for k, i2 in enumerate(_eids):
        _eid_of.setdefault(i2, k)
    _param_of = {}          # edge index -> {"from"/"to": target edge id}

    for k, ed in enumerate(recipe["edges"]):
        where = "edges[%d]" % k
        check_keys(dict(ed), EDGE_KEYS, where)
        fr, to = _end_addr(ed["from"], where + ".from", by_as, gnodes, _eid_of), \
            _end_addr(ed["to"], where + ".to", by_as, gnodes, _eid_of)
        for end, v in (("from", fr), ("to", to)):
            if isinstance(v, dict):
                _param_of.setdefault(k, {})[end] = v["of"]
        if fr == to:
            raise ValueError("%s: 'from' and 'to' are both %r; an edge needs two "
                             "different ends." % (where, fr))
        feed = dict(ed.get("feeders") or {})
        if feed:
            check_keys(feed, FEEDERS_KEYS, where + ".feeders")
        e = {
            "as": _eids[k],
            "from": fr, "to": to,
            "ring": ed.get("ring", (rg.get("as") or RING_AS)),
            "name": dict(ed["name"]) if ed.get("name") else None,
            "feeders": feed or None,
            "style": ed.get("style"),
            "reverse": bool(ed.get("reverse")),
            # `crossing` says whether a layout crossing is DELIBERATE (the
            # engine cannot tell a coupling marker from a typesetting accident --
            # opposite meaning, identical geometry), so it is declared, never
            # inferred.
            "crossing": _crossing(ed.get("crossing"), where)[0],
            # The edge `crossed` names, or None. Kept beside `crossing` rather than
            # folded into it so every existing `== "bridge"` read still works.
            "crossed_with": _crossing(ed.get("crossing"), where)[1],
            # `bulge` bows the edge; it is how one of several parallel edges gets
            # out of the way of its twins. Explicit, not derived -- there is no
            # automatic fan.
            "bulge": _bulge(ed.get("bulge"), where),
            # DERIVED from path membership below, never written on the edge. One
            # rule, one owner: muting is a property of a SET (the path table).
            "muted": False,
            # None when that end is a POINT ON ANOTHER EDGE rather than a node.
            # Every pass that indexes nodes checks for it; `edge_shape` reads it as
            # "on no ring", which makes the edge off-ring by the off-ring rule.
            "_from": by_as.get(fr) if isinstance(fr, str) else None,
            "_to": by_as.get(to) if isinstance(to, str) else None,
        }
        e["_shape"] = edge_shape(e, grings, {g["_i"]: g["_rings"] for g in gnodes})
        # WHICH ring's arc it is: the classifier picked one when the ends share
        # several (a shared node). An explicit `ring` on the edge wins, because the
        # caller said it.
        if e["_shape"] == "arc":
            e["ring"] = ed.get("ring") or e.get("_ring_pick") or grings[0]["as"]
        else:
            # A chord, a cross-ring edge or an off-ring edge is NOT any ring's arc,
            # so it does not claim one unless the caller said so. Leaving the
            # default in place would make a cross-ring edge report `ring: r1`, which
            # would then feed the shared-step test a ring the edge does not belong
            # to.
            e["ring"] = ed.get("ring")
        e.pop("_ring_pick", None)
        # `via` DEFAULTS FROM THE SHAPE and may be overridden. An arc of the ring
        # is drawn as an arc; a chord or an off-ring edge has no ring arc to lie
        # on, so it is a straight line. The caller may say otherwise.
        via = ed.get("via") or ("arc" if e["_shape"] == "arc" else "line")
        if via not in VIA:
            raise ValueError(
                "%s: 'via' is %s, got %r. 'polyline' is not supported."
                % (where, " / ".join(map(repr, VIA)), via))
        if via == "arc" and e["_shape"] != "arc":
            raise ValueError(
                "%s: via 'arc' means THIS RING's arc, and %r -> %r is a %s -- "
                "there is no ring arc between them. Use via 'line'."
                % (where, fr, to, e["_shape"]))
        e["via"] = via
        # A HOP IS BUILT ON A STRAIGHT RUN. The pieces are cut by arc length along
        # the edge's own polyline, and reconstructing them along a bow or around a
        # ring arc is a different piece of geometry -- refused rather than drawn
        # approximately, since an approximate hop sits beside the crossing instead
        # of on it and still looks deliberate.
        if e["crossing"] == "bridge":
            if via != "line":
                raise ValueError(
                    "%s: 'crossing: \"bridge\"' is built for a straight edge "
                    "(via 'line'); this one is via %r. A hop along a ring arc is "
                    "not implemented." % (where, via))
            if e["bulge"] is not None:
                raise ValueError(
                    "%s: 'crossing: \"bridge\"' and 'bulge' together are not "
                    "implemented -- the hop is cut out of a straight run. Bow the "
                    "edge OR hop it, not both." % where)
        # NOT SILENTLY IGNORED. Four keys below are meaningful on a ring arc and
        # meaningless on a straight line; refusing costs a caller nothing, and
        # ignoring leaves a caller who never finds out why the figure is wrong.
        # Same format as `crossing` / `muted` / `content` / `graph` / `frame`: say
        # what it is and that it is not supported.
        if via == "line":
            if (e["name"] or {}).get("inside"):
                raise ValueError(
                    "%s: 'inside' puts a caption in the MIDDLE OF THE RING, which "
                    "a straight edge has no notion of -- it is drawn between two "
                    "points, not around a circle. This edge is a %s, so it is a "
                    "line. Use 'above' / 'below', which on a line mean the two "
                    "sides of it." % (where, e["_shape"]))
            if e["feeders"]:
                raise ValueError(
                    "%s: 'feeders' is an arc TANGENT to the ring with the reagent at its "
                    "free end, so it only exists on a ring arc; this edge is a %s and is "
                    "drawn as a line. Feeders on a straight edge are not supported -- put "
                    "the species on its own node and draw a second edge to it." % (where, e["_shape"]))
        gedges.append(e)


    # ---- ONE LEVEL OF PARAMETRIC ENDPOINT ------------------------------------
    #
    # An edge terminating on an edge needs that edge's INK, and ink exists only
    # once the target has been drawn and filed. Targets are drawn first (see the
    # loop's sort key), which works for one level; a chain would need the order
    # solved. Refused rather than ordered by luck -- the failure would be an
    # address that resolves to nothing, blamed on the wrong edge.
    for _k, _ends in _param_of.items():
        for _end, _tgt in _ends.items():
            _ti = _eid_of[_tgt]
            # SELF FIRST: an edge pointing at itself is also "param-ended", so the
            # chain check below would catch it and blame the wrong thing.
            if _ti == _k:
                raise ValueError(
                    "edges[%d].%s: an edge cannot terminate on ITSELF (%r)."
                    % (_k, _end, _tgt))
            if _ti in _param_of:
                raise ValueError(
                    "edges[%d].%s: %r itself terminates on another edge, and a "
                    "chain of parametric endpoints is not implemented -- its ink "
                    "does not exist until it is drawn. Point at an edge whose two "
                    "ends are nodes." % (_k, _end, _tgt))
    for _e in gedges:
        _e["_param"] = bool(isinstance(_e["from"], dict)
                            or isinstance(_e["to"], dict))

    # ---- `crossed`: whom it exempts, and in both directions -------------------
    #
    # ONE DECLARATION COVERS THE PAIR. The collision is reported against whichever
    # of the two is placed second, and the caller cannot know the placement order,
    # so writing it on one edge has to be enough. Recorded symmetrically here
    # rather than requiring both edges to say it: two keys that have to agree is a
    # fragile shape.
    _xset: dict = {}
    for _e in gedges:
        _t = _e.get("crossed_with")
        if _t is None:
            continue
        if _t not in _eid_of:
            raise ValueError(
                "edges[%s].crossing.crossed: %r is not an edge in this figure. "
                "The edges are %s."
                % (_e["as"], _t, ", ".join(map(repr, _eid_of)) or "none"))
        if _t == _e["as"]:
            raise ValueError(
                "edges[%s].crossing.crossed: an edge cannot deliberately cross "
                "ITSELF (%r)." % (_e["as"], _t))
        _xset.setdefault(_e["as"], set()).add(_t)
        _xset.setdefault(_t, set()).add(_e["as"])
    for _e in gedges:
        _e["_crossed_set"] = tuple(sorted(_xset.get(_e["as"], ())))
        # TOUCHES, not crosses -- a separate set because `_crossed_set` carries an
        # ASSERTION (there is a crossing here; the engine refuses the declaration if
        # the inks do not actually cross) and two tangent arcs meet at one point
        # without crossing. Only `tangent_to` fills this in; no caller key sets it.
        _e["_touch_set"] = ()

    # ---- the remaining two shapes, which are properties of the SET -----------
    #
    # `arc` / `chord` / `crossring` / `offring` are properties of ONE edge, so
    # `edge_shape` decides them. These two are not:
    #
    #   PARALLEL PATH   several edges between the SAME unordered pair. No single
    #                   edge can tell; it takes counting the set.
    #   SHARED STEP     a marker whose incident edges belong to MORE THAN ONE ring
    #                   -- one marker node carrying edges from two different
    #                   rings. A property of the node, read off its edges.
    #
    # Recorded, not acted on here: a parallel path is bowed away from its twin by
    # the caller's `bulge`. This says WHAT it is, which is what a test asserts --
    # the kind the engine judged, not the picture.
    # KEYED ON THE ADDRESS, not on the node index: a parametric end has no index,
    # so two edges landing on two DIFFERENT edges would both key as None and be
    # called parallel when they share nothing.
    def _endkey(v):
        return v if isinstance(v, str) else ("@", v["of"], v["at"])
    _pairs_seen: dict = {}
    for _e in gedges:
        _pairs_seen.setdefault(
            frozenset((_endkey(_e["from"]), _endkey(_e["to"]))), []).append(_e)
    for _grp in _pairs_seen.values():
        for _e in _grp:
            _e["_parallel"] = len(_grp) if len(_grp) > 1 else 0
    # THE RINGS OF THE FAR ENDS, not of the edges. A marker is pinned, so it is on
    # no ring, so its own edges are `offring` and carry no `ring` at all -- asking
    # them always answers False. What makes a marker a SHARED STEP is that the
    # things it connects belong to different rings.
    _far: dict = {}
    for _e in gedges:
        for _end, _other in ((_e["_from"], _e["_to"]), (_e["_to"], _e["_from"])):
            if _end is None or _other is None:
                continue            # a parametric end is on no node and no ring
            _far.setdefault(_end, set()).update(gnodes[_other]["_rings"])
    for _g in gnodes:
        _g["_shared_step"] = bool(_g["kind"] == "marker"
                                  and len(_far.get(_g["_i"], set())) > 1)

    # ---- MARKER NODES: the edge pairing must be DECLARED ----------------------
    #
    # This is among the things the engine CANNOT infer: an e- circle joined by four
    # edges, paired crosswise two by two -- pair them wrongly and the arrows point
    # the wrong way, while the drawing still looks right. That last clause is the
    # whole reason this raises instead of guessing: a wrong pairing produces a
    # figure that looks plausible and says something false, which is the one
    # outcome declared-not-inferred rules exist to prevent.
    #
    # Two incident edges need no declaration -- there is only one pairing. Three or
    # more do, and the count is what decides, not the kind.
    _inc: dict = {}
    for _k, _e in enumerate(gedges):
        for _end in (_e["_from"], _e["_to"]):
            if _end is None:
                continue            # a point on an edge is not incident to a node
            _inc.setdefault(_end, []).append(_e["as"])
    for _g in gnodes:
        if _g["kind"] != "marker":
            if _g.get("pairs"):
                raise ValueError(
                    "nodes[%d]: 'pairs' belongs on a 'marker' node; %r is a %s."
                    % (_g["_i"], _g["as"], _g["kind"]))
            continue
        _mine = _inc.get(_g["_i"], [])
        _pairs = _g.get("pairs")
        if _pairs is None:
            if len(_mine) > 2:
                raise ValueError(
                    "nodes[%d]: marker %r has %d edges (%s) and no 'pairs'. Which edge "
                    "continues into which CANNOT be inferred: getting it wrong reverses an "
                    "arrow and still draws something that looks right. Declare it, e.g. "
                    "\"pairs\": [[%r, %r]]."
                    % (_g["_i"], _g["as"], len(_mine),
                       ", ".join(map(repr, _mine)), _mine[0], _mine[1]))
            _g["_pairs"] = ()
            continue
        if not isinstance(_pairs, (list, tuple)):
            raise ValueError(
                "nodes[%d].pairs: a list of two-edge pairs, e.g. "
                "[[\"in1\", \"out1\"], [\"in2\", \"out2\"]]; got %s"
                % (_g["_i"], type(_pairs).__name__))
        _used: dict = {}
        _out = []
        for _pi, _pr in enumerate(_pairs):
            _w = "nodes[%d].pairs[%d]" % (_g["_i"], _pi)
            if not isinstance(_pr, (list, tuple)) or len(_pr) != 2:
                raise ValueError("%s: a pair is exactly two edge names, got %r"
                                 % (_w, _pr))
            _a2, _b2 = str(_pr[0]), str(_pr[1])
            if _a2 == _b2:
                raise ValueError("%s: an edge cannot be paired with itself (%r)"
                                 % (_w, _a2))
            for _nm in (_a2, _b2):
                if _nm not in _mine:
                    raise ValueError(
                        "%s: %r is not an edge at marker %r. Its edges are %s."
                        % (_w, _nm, _g["as"], ", ".join(map(repr, _mine)) or "none"))
                if _nm in _used:
                    raise ValueError(
                        "%s: %r is already paired with %r. Each edge continues "
                        "into exactly one other." % (_w, _nm, _used[_nm]))
            _used[_a2], _used[_b2] = _b2, _a2
            _out.append((_a2, _b2))
        _unpaired = [x for x in _mine if x not in _used]
        if _unpaired and len(_mine) > 2:
            raise ValueError(
                "nodes[%d].pairs: marker %r leaves %s unpaired. With more than two edges "
                "every one has to say what it continues into, or the ones left over "
                "would be drawn with no indication of where they lead."
                % (_g["_i"], _g["as"], ", ".join(map(repr, _unpaired))))
        _g["_pairs"] = tuple(_out)

    # The old `center` argument still works in the new form: it becomes ring 0's
    # label at position "center", which is what it desugars to in the old form too.
    # An explicit `label` on that ring wins -- the caller said it twice and the
    # specific one is the one they wrote on the ring.
    note = recipe.get("center")
    if note and grings[0].get("label") is None:
        grings[0]["label"] = {"text": note, "position": "center"}

    # WHICH RINGS HAVE AN ALLOCATABLE INTERIOR.
    #
    # Per ring, never a global default -- that is the whole point. A ring gets it
    # when something is SUPPOSED to be in its middle: another ring nested inside
    # it, or a node pinned inside its circle. Only then does a box lying across
    # the centre stop being a reason to grow the circle.
    #
    # `cfg["interior_allocatable"]` forces it on for every ring, for a caller who
    # knows better; it defaults to False, so no existing figure changes.
    _forced = bool(cfg.get("interior_allocatable"))
    # ...AND RINGS THAT ARE ABOUT TO BE SPREAD ARE NOT NESTED. Every ring with no
    # centre of its own sits at ORIGIN here, so two of them read as each being
    # inside the other and BOTH would claim an allocatable interior -- on the
    # strength of a coincidence the packer is about to undo. _spread_rings gives
    # them their own centres once the radii are known; until then, treat them as
    # what they will be, which is side by side.
    # A `tangent_to` ring is also about to be given a centre, so it is no more
    # "nested in" its neighbour than a spread one is -- same coincidence at ORIGIN,
    # same undoing a few hundred lines down.
    _will_spread = len(grings) > 1 and (
        not any(g["_centre_given"] for g in grings)
        or any(g["tangent_to"] for g in grings))
    for gr in grings:
        inner = [] if _will_spread else [
            o for o in grings
            if o is not gr and _inside_ring(o["centre"], gr)]
        pinned_in = [g for g in gnodes
                     if g["at"] is not None and _inside_ring(g["at"], gr)]
        gr["_interior_ok"] = bool(_forced or inner or pinned_in)
        gr["_nested_in"] = [] if _will_spread else [
            o["as"] for o in grings
            if o is not gr and _inside_ring(gr["centre"], o)]

    # ---- PATHS: the fourth table ---------------------------------------------
    #
    # A path is a NAMED SET OF EDGES. It exists because two features are
    # properties of a set and of nothing smaller:
    #
    #   the group-level NAME -- a "TM", an "ROH", "path A" / "path B". Written
    #     on a member instead, a second member could write it again and the label
    #     would draw twice with no rule able to catch it.
    #   MUTING       -- greying out as a group: edges, their labels and attached
    #     decorations all go grey together. That is state spanning several
    #     objects, not a color set on each edge separately. Declared once here;
    #     the engine propagates it to each member's arrow, caption and feeders.
    #
    # Muting ONE edge is a one-edge path. That is the cost of having one owner for
    # the rule, and it is the cheaper side of the trade.
    by_edge = {e["as"]: e for e in gedges}
    gpaths = []
    seen_path: dict = {}
    for pi, raw in enumerate(recipe.get("paths") or []):
        where = "paths[%d]" % pi
        if not isinstance(raw, dict):
            raise ValueError("%s: a path is a dict with 'edges', got %s"
                             % (where, type(raw).__name__))
        check_keys(raw, PATH_KEYS, where)
        mem = [str(v) for v in (raw.get("edges") or [])]
        if not mem:
            raise ValueError(
                "%s: a path needs 'edges' -- the names of the edges it groups. An "
                "empty path has nothing to name and nothing to mute." % where)
        for m in mem:
            if m not in by_edge:
                raise ValueError(
                    "%s: %r is not an edge in this figure. The edges are %s."
                    % (where, m, ", ".join(repr(e["as"]) for e in gedges)))
        if len(set(mem)) != len(mem):
            dup = [m for m in dict.fromkeys(mem) if mem.count(m) > 1]
            raise ValueError("%s.edges: %s listed twice in the same path."
                             % (where, ", ".join(map(repr, dup))))
        pid = _as_id(raw.get("as"), where, "path#%d" % pi)
        if pid in seen_path:
            raise ValueError("%s: the path name %r is already used by %s."
                             % (where, pid, seen_path[pid]))
        seen_path[pid] = where
        gp = {"as": pid, "edges": mem,
              "name": dict(raw["name"]) if raw.get("name") else None,
              "muted": bool(raw.get("muted")),
              "_i": pi}
        if gp["name"]:
            check_keys(gp["name"], STEP_KEYS, where + ".name")
        gpaths.append(gp)
        # OR, not assignment: an edge in two paths is muted if EITHER mutes it.
        # Nothing here can un-mute, so the order the paths are written in cannot
        # change the answer.
        if gp["muted"]:
            for m in mem:
                by_edge[m]["muted"] = True

    return {"nodes": gnodes, "edges": gedges, "rings": grings,
            "paths": gpaths, "_by_as": by_as}


def _step_view(edge: dict) -> dict:
    """An edge's payload, back in the dict shape `feeder_motif` / `_arc_kw` take.

    DERIVED FROM THE TABLE, not carried across from the original step. That is the
    difference between the edge table driving the drawing and the edge table being
    a parallel structure nobody reads: a mutation that blanks `name` / `feeders` /
    `style` / `reverse` would leave the regression figures GREEN if the consumers
    still held the original dict, and turns them red because they read this.

    Nothing is lost in the round trip, and that is checkable rather than hopeful:
    STEP_KEYS is closed -- above / below / inside / reverse / style / at / in / out
    / aside / feeder_style / as -- and every one of those has a home in the edge
    row (`name`, `feeders`, `style`, `reverse`, `as`). A test asserts the round
    trip over a step that uses all eleven.

    `at` lives under `feeders` in the table (it is where the feeder meets the arc,
    which is a feeder property) and comes back to the top level here, because the
    step CAPTION also consults it to decide which way to slide off the crossing
    point. One value, one place it is stored, two readers.
    """
    st = dict(edge.get("name") or {})
    feed = dict(edge.get("feeders") or {})
    if "style" in feed:
        st["feeder_style"] = feed.pop("style")
    st.update(feed)                     # in / out / aside / at
    if edge.get("style") is not None:
        st["style"] = edge["style"]
    if edge.get("reverse"):
        st["reverse"] = True
    st["as"] = edge["as"]
    return st


def desugar(recipe: dict, cfg: dict, meas) -> dict:
    """`nodes` + `steps` -> the three tables (nodes, edges, rings). INTERNAL.

    WHAT IS WRONG WITH THE OLD SHAPE, in one sentence: `steps[i]` connects
    `nodes[i]` to `nodes[i+1]` and NOTHING SAYS SO -- the topology lives in the
    list order, so there is no field that could ever say "this arrow goes from 4
    to 5". This function is where that stops being true: after it, the edge table
    states each edge's two endpoints, and the drawing loop reads them instead of
    recomputing `(i + 1) % n`.

    THIS CHANGES NO PUBLIC INTERFACE. `draw_catalytic_cycle`'s signature,
    parameters and documentation are untouched; this is a translation that happens
    inside the engine, and figures written in the old form must come out byte for
    byte identical. So every value here is carried across unchanged and in the
    same order -- the tables are a different SHAPE for the same numbers, not a
    recalculation of them.

    ONE RING. The old form yields exactly one `rings` entry. The `graph` field of a
    subgraph node is in the schema and is not produced or consumed here: the schema
    is recursive, and the implementation may support one level.

    FIELDS PREFIXED `_` ARE PRIVATE and are not part of the table schema. They
    carry what the existing geometry already had in hand -- the resolved endpoint
    INDICES, and the original step dict -- so that the representation can change
    without touching one arithmetic expression. `_from` / `_to` are what the loops
    use; `from` / `to` are the addresses that say the same thing, and a test
    asserts the two agree for every edge, so the address fields are load-bearing
    rather than decorative.
    """
    nodes_in = recipe.get("nodes") or []
    steps = list(recipe.get("steps") or [])
    n = len(nodes_in)

    # ---- nodes -----------------------------------------------------------
    # Shared by BOTH input forms: a node row looks the same whichever table the
    # edges came from, which is what lets `nodes` be one parameter.
    # `as` resolved HERE, once, so the address name has a single definition and the
    # drawing loop cannot disagree with the edge table about what a node is called.
    gnodes = []
    for i, m in enumerate(meas):
        nd = m["nd"]
        gnodes.append({
            "as": _as_id(nd.get("as"), "nodes[%d]" % i, "node#%d" % i),
            # The node vocabulary is a SUBSET of the address layer's: a node is
            # only ever structure / text / marker / subgraph. The old nodes + steps
            # form produces the two it can express.
            # `marker` is the caller's own word and cannot be inferred -- a SET box
            # and a written species are both type on the page, and only the caller
            # knows which is which. That is why the marker pairing rule can key off
            # it.
            "kind": (nd.get("kind") or _inferred_kind(m)),
            "pairs": nd.get("pairs"),
            # The node's `frame`. Only the SQUARE BRACKET form is built; see _frame.
            "frame": _frame(nd.get("frame"), "nodes[%d].frame" % i),
            "s": nd.get("s"),
            "text": nd.get("text"),
            "label": m["lbl"],
            # Explicit coordinates, or None and a ring places it.
            "at": (None if nd.get("at") is None
                   else _xy(nd["at"], "nodes[%d].at" % i)),
            "_i": i,
            # Which ring places this node, and every ring that lists it. A SHARED
            # node (the same id in two rings' `nodes`) is placed by the FIRST
            # one -- decided here, once, so the drawing loop does not re-derive it.
            "_ring": None,
            "_rings": (),
        })
    # NO duplicate check here. The address Registry already refuses a repeat and
    # its message names BOTH origins, which a check at this point cannot --
    # nothing has been registered yet, so there is no first origin to name. A
    # duplicate therefore reaches the node loop and raises there, with the better
    # message. One rule, one owner.
    by_as = {g["as"]: g["_i"] for g in gnodes}

    # ---- the NEW form: edges + rings come from the recipe ------------------
    if recipe.get("edges") is not None:
        return _graph_from_tables(recipe, cfg, gnodes, by_as)

    # ---- edges -----------------------------------------------------------
    # The arc count, unchanged: one per step, closing the ring (n) or leaving it
    # open (n - 1); a step past the end is validated and never drawn.
    n_arcs = n if len(steps) >= n else max(0, min(len(steps), n - 1)) or (n - 1)
    if not steps:
        n_arcs = n

    gedges = []
    for k in range(n_arcs):
        st = (steps[k] if k < len(steps) else {}) or {}
        j = (k + 1) % n
        cap = {key: st[key] for key in ("above", "below", "inside") if key in st}
        feed = {key: st[key] for key in ("in", "out", "aside", "at")
                if key in st}
        if st.get("feeder_style"):
            feed["style"] = st["feeder_style"]
        gedges.append({
            "as": _as_id(st.get("as"), "steps[%d]" % k, "edge#%d" % k),
            # THE POINT OF THE WHOLE EXERCISE: the endpoints are stated.
            "from": gnodes[k]["as"],
            "to": gnodes[j]["as"],
            "via": "arc",                   # only the new form can say "line"
            "ring": RING_AS,                # the old form has exactly one
            "name": cap or None,
            "feeders": feed or None,
            # Carried because the old shape has both and the translation may not
            # drop a public key. `style` is the per-step arc_arrow override,
            # `reverse` swaps the heads against the ring.
            "style": st.get("style"),
            "reverse": bool(st.get("reverse")),
            # The "must be declared" pair, at their documented defaults. Nothing
            # in the old shape can set them, so this branch only ever writes these.
            "crossing": "plain",
            "crossed_with": None,       # the old form has no crossing syntax
            "muted": False,
            "_from": k,
            "_to": j,
            # The old form can only ever produce ring arcs: every node is on the
            # ring and every step joins two neighbours. Set explicitly so that
            # render_cycle reads ONE field whichever form the caller used.
            "_shape": "arc",
        })

    # ---- rings -----------------------------------------------------------
    for g in gnodes:
        g["_ring"], g["_rings"] = RING_AS, (RING_AS,)
    note = recipe.get("center")
    ring = {
        "as": RING_AS,
        "nodes": [g["as"] for g in gnodes],
        # Every node, in order: the old form's ring IS the whole node list.
        "_idx": list(range(n)),
        "_pos": {i: i for i in range(n)},
        "centre": (0.0, 0.0),           # LAYOUT origin; multi-ring layout moves rings
        "radius": recipe.get("radius") or cfg.get("radius") or "auto",
        "start_angle": cfg["start_angle"],
        "direction": cfg["direction"],
        # The centre note becomes the ring's label at position "center". The old
        # form has only one thing in the middle, so contention there (label vs
        # note vs nested ring) cannot arise.
        "label": ({"text": note, "position": "center"} if note else None),
        # The old form has one ring and nothing pinned, so its middle holds at most
        # the centre note -- the regime the packer was built for.
        "_interior_ok": bool(cfg.get("interior_allocatable")),
        "_nested_in": [],
    }
    for _e in gedges:
        _e["_param"] = False        # the old form has no parametric endpoint
        _e["_crossed_set"] = ()
        _e["_touch_set"] = ()           # see the note on the new-form branch
    # `paths` is always empty from the OLD form: it has no syntax for one, and an
    # empty list is what keeps every downstream loop from needing a `.get`.
    return {"nodes": gnodes, "edges": gedges, "rings": [ring], "paths": [],
            "_by_as": by_as}


_BUDGET = 360.0

# The identity that protects the regression figures -- n * _node_sep(legacy_radius)
# == 360 -- is exact in real arithmetic and misses by up to ~6e-14 degrees in IEEE
# double (measured: 6 of 100 synthetic n/reach pairs land just the wrong side of
# it). 1e-9 degrees is 2e-9 page units at R = 100: eleven orders of magnitude below
# the two-decimal grid every CDXML coordinate is written on, and the difference
# between "a ring with no feeder is byte-identical BY CONSTRUCTION" and
# "byte-identical as long as the arithmetic keeps being kind".
_BUDGET_TOL = 1e-9


# ---- WHERE ON THE CIRCLE A STEP'S ARC SITS -------------------------------
#
# Two rings drawn TANGENT (the usual way published figures draw a dual cycle, with
# the SET / PT / HAT badge printed ON the touch point) put their closest approach
# at ONE point. A node sitting at that point lands its label on the other ring --
# measured 14.3 units of block-on-block on a tangent pair. The fix is not a nudge:
# the ring has to be ROTATED so the tangency direction falls in the MIDDLE OF AN
# ARC, and which arc is a chemistry question (the shared step is the one that
# belongs there), so the caller names the step and this pair of functions does the
# arithmetic.
#
# MEASURED, not reasoned: 6 test drawings (n = 3 / 4 / 6, start_angle = 0 / 90 / 270,
# cw and ccw) read the drawn arcs' own Center3D / Tail3D / AngularSize back out.
# Every midpoint landed on the formula below to 0.0 degrees. Angles are in the
# CDXML frame -- y grows DOWNWARD, so 0 is right and 90 is BELOW, and `cw` means
# increasing angle, which is what `_arc_point` already does.

def arc_mid_angle(start_angle: float, step: int, n: int,
                  direction: str = "cw") -> float:
    """Angle of the MIDDLE of step `step`'s arc on an n-node ring."""
    sgn = 1.0 if direction == "cw" else -1.0
    return (start_angle + sgn * (360.0 / n) * (step + 0.5)) % 360.0


def start_angle_for(theta: float, step: int, n: int,
                    direction: str = "cw") -> float:
    """The `start_angle` that puts step `step`'s arc midpoint at `theta`.

    The inverse of arc_mid_angle, and the reason it exists: `theta` is where the
    other ring touches, so this is how a ring gets rotated to meet it arc-first
    instead of node-first.
    """
    sgn = 1.0 if direction == "cw" else -1.0
    return (theta - sgn * (360.0 / n) * (step + 0.5)) % 360.0


def rotate_to_face(start_angle: float, theta: float, n: int,
                   direction: str = "cw") -> float:
    """`start_angle` rotated by the SMALLEST amount that puts an arc midpoint at
    `theta` -- i.e. the ring meets that direction arc-first instead of node-first.

    Smallest, not "put step 0 there", because rotating a ring moves every species
    on it: the figure's reading is the caller's, and this only has to buy the one
    thing tangency needs. The move is never more than half a step (180/n degrees),
    and which step ends up facing `theta` falls out of where the ring already was.

    A ring whose arcs already straddle `theta` comes back unchanged.
    """
    adv = 360.0 / n
    best = None
    for k in range(n):
        d = ((theta - arc_mid_angle(start_angle, k, n, direction) + 180.0)
             % 360.0) - 180.0
        if best is None or abs(d) < abs(best):
            best = d
    return (start_angle + best) % 360.0


def legacy_radius(meas, cfg, note_reach: float = 0.0) -> float:
    """Radius that satisfies BOTH constraints on the ring.

    Neighbours must not touch each other, and the centre note must not reach the
    species. Without the second, a wide net-reaction line in the middle creeps up
    on the left and right species.

    This is the pre-packer formula, and it stays -- PUBLIC, and named for what it
    is -- because it is the solver's floor.
    solve_radius never returns less than this, so no figure that was already right
    can be made smaller, and every figure with no feeder in it gets this number
    back unchanged. That is the whole regression guarantee, and it is one `<=`
    away: at exactly this radius the n gaps demand exactly 360 degrees.
    """
    n = len(meas)
    reach = max(m["reach"] for m in meas)
    by_neighbours = 2.0 * reach + cfg["node_gap"]
    if n >= 2:
        by_neighbours /= (2.0 * math.sin(math.pi / n))
    by_note = note_reach + reach + cfg["center_gap"] if note_reach else 0.0
    return max(cfg["min_radius"], by_neighbours, by_note)


def _trim(reach: float, R: float, cfg: dict) -> float:
    """Half-angle of a species box plus its arc margin: how far an arc backs off."""
    return math.degrees(math.atan2(reach + cfg["arc_margin"], R))


def _node_sep(meas, cfg, R: float) -> float:
    """Centre-to-centre angle two neighbouring species need at radius R.

    legacy_radius asks for a CHORD of 2*reach + node_gap between neighbouring
    centres; this is that same demand read as an angle, 2*asin(half the chord / R).
    Hence, at R = legacy_radius, exactly 2*pi/n -- n of these come to exactly 360
    and the bisection returns its own lower bound untouched.

    ONE reach for every node, the largest, exactly as legacy_radius does it. Using
    each node's own reach would be more honest and would ask for LESS, which is the
    one direction this must not go: it would shrink a ring that was already right.
    """
    h = max(m["reach"] for m in meas) + cfg["node_gap"] / 2.0
    return 2.0 * math.degrees(math.asin(min(1.0, h / R)))


def _arc_of(ang_i: float, adv_i: float, reach_i: float, reach_j: float,
            R: float, cfg: dict, s: float):
    """(t0, sweep) of the arc carrying node i to node i+1, in degrees.

    `adv_i` is the signed angle from node i to node i+1 -- 360/n on an evenly
    divided ring, whatever the packer allotted otherwise. ONE definition, called
    by the drawing loop and by the packer's trial run both: a budget measured in a
    different frame from the one the ink lands in measures one thing and draws
    another.
    """
    di, dj = _trim(reach_i, R, cfg), _trim(reach_j, R, cfg)
    t0 = ang_i + s * di
    sweep = adv_i - s * (di + dj)
    if abs(sweep) < cfg["min_sweep"]:
        sweep = s * cfg["min_sweep"]
    return t0, sweep


def _arc_need(motif: FeederMotif, at: float, cfg: dict) -> float:
    """How long this step's arc must be for its feeder to fit beside the species.

    The tangency sits at parameter `at` along the arc, so the upstream reach has
    only `at` of the arc to live in and the downstream reach the remaining 1 - at:
    hence the division, not up + down. An arc also has a shortest readable length
    whether or not anything hangs off it, which is the min_sweep floor -- inert on
    the reference regression figures (on all of them the neighbour spacing is the
    binding constraint by 3 degrees or more), and biting only where an arc would
    otherwise come out too short to read.

    `at` is held off the ends: an author may pin a tangency at 0.02, and 50x its
    upstream reach is a demand no radius can meet, whereas the small overlap left
    by clamping is something the allocator reports honestly.
    """
    a = min(0.95, max(0.05, at))
    return max(cfg["min_sweep"], motif.up / a, motif.down / (1.0 - a))


def _gap_demands(meas, cfg, R: float, arc_need: dict) -> list[float]:
    """Angle each gap between neighbouring species demands at radius R.

    Two claims on the same gap, whichever is larger: the neighbours must not touch
    (_node_sep, the legacy formula), and the arc between them plus its feeder must
    fit between their boxes (the trim angles are where the arc already backs off
    from each species, so adding the feeder's reach to them is asking for room
    beside the species rather than through them).
    """
    n = len(meas)
    sep = _node_sep(meas, cfg, R)
    out = []
    for i in range(n):
        j = (i + 1) % n
        need = max(cfg["min_sweep"], arc_need.get(i, 0.0))
        out.append(max(sep, _trim(meas[i]["reach"], R, cfg)
                       + _trim(meas[j]["reach"], R, cfg) + need))
    return out


def _disk_clash(rects, note_reach: float, centre=ORIGIN) -> float:
    """How far the innermost feeder box intrudes on the centre note, or 0.

    The other half of the packing problem: split's by-product label is SUPPOSED to
    land inside the cycle, where it spends no ring circumference at all but can
    walk straight into the net-reaction note. Growing R moves it outwards -- its
    arc is tangent to the ring, so the whole motif rides out with the circle --
    while the note keeps its size at the centre. Monotone, so it rides along inside
    the same bisection.

    Asked of EVERY box, not only the ones counted as inside the ring: a wide label
    can be charged to the angle budget for poking out through the species and still
    have its far end lying across the middle. Measured on a five-species ring: a
    134-unit by-product label reaching from radius 2 to radius 137, correctly billed
    for its angle, 4.5 units into the note.

    Radial, and against the note's circumscribed radius rather than its box, for the
    same reason _box_disk exists: the budget is solved before the angles are handed
    out, so an answer that depends on where round the ring the feeder sits clears the
    note in the frame it was measured in and walks back into it in the frame that
    gets drawn (measured: 6.5 units on a two-species ring). This is the same
    rotation-independent question legacy_radius already asks of the species.
    """
    if not note_reach:
        return 0.0
    return max((note_reach - (d - rho)
                for d, rho in (_box_disk(rc, centre) for rc in rects)),
               default=0.0)


def solve_radius(meas, cfg, note_reach: float = 0.0, probe=None,
                 trace: dict | None = None) -> float:
    """Smallest radius whose 360 degrees pay for everything. BISECTION.

    `probe(R)` returns ({step index: arc angle needed}, [in-ring boxes]) for the
    feeders at radius R -- None when the ring has none, which is the common case and
    returns legacy_radius immediately.

    Bisection, not fixed-point iteration: `need` is monotone decreasing in R, on
    which bisection converges unconditionally, while iterating R := f(R) can
    oscillate whenever a feeder's own radius scales with R (it does, up to its
    cap, and the packer must be correct either way).

    Not converging RAISES. A radius quietly rounded off to whatever the loop had
    reached would produce a figure that is subtly wrong with nothing to point at --
    the same failure mode as a silently dropped recipe key, and the message carries
    the demand at both ends so the caller can see what was unpayable.
    """
    lo = legacy_radius(meas, cfg, note_reach)
    if probe is None:
        if trace is not None:
            trace.update(solver="legacy", tries=0, need=None)
        return lo

    def demand(R: float):
        arc_need, boxes = probe(R)
        return (sum(_gap_demands(meas, cfg, R, arc_need)),
                _disk_clash(boxes, note_reach))

    def fits(R: float) -> bool:
        need, clash = demand(R)
        return need <= _BUDGET + _BUDGET_TOL and clash <= 0.0

    tries = 1
    if fits(lo):
        # No feeder asked for more than the ring was already going to give it.
        # Return the SAME float legacy_radius produced, not a bisected
        # approximation of it: this is the byte-for-byte path.
        if trace is not None:
            trace.update(solver="floor", tries=tries, need=round(demand(lo)[0], 2))
        return lo

    hi = lo
    for _ in range(20):
        hi *= 1.5
        tries += 1
        if fits(hi):
            break
    else:
        need_lo, clash_lo = demand(lo)
        need_hi, clash_hi = demand(hi)
        raise ValueError(
            "the ring cannot be made big enough for its feeders: at R = %.1f the "
            "%d gaps demand %.1f of 360 degrees (in-ring overlap %.1f), and at R = "
            "%.1f -- %.0fx bigger -- they still demand %.1f (overlap %.1f). "
            "Something is asking for space that does not shrink with the radius: "
            "check for an enormous species box, a caption pinned into a feeder's "
            "'at', or a centre note wider than the cycle."
            % (lo, len(meas), need_lo, clash_lo, hi, hi / lo, need_hi, clash_hi))

    for _ in range(60):
        if hi - lo < 0.01:
            break
        mid = (lo + hi) / 2.0
        tries += 1
        if fits(mid):
            hi = mid
        else:
            lo = mid
    else:
        raise ValueError(
            "the radius solver did not converge: 60 halvings left the answer "
            "somewhere in [%.4f, %.4f], which cannot happen for a monotone "
            "demand -- suspect a feeder whose reach GROWS with the radius."
            % (lo, hi))
    if trace is not None:
        trace.update(solver="bisect", tries=tries, need=round(demand(hi)[0], 2))
    return hi    # the feasible end. lo is known NOT to fit.


def allocate_gaps(demands, warn=None) -> list[float]:
    """Hand out the 360 degrees: every gap its demand, then the slack shared EVENLY.

    Evenly, not in proportion to demand. Proportional sharing gives the widest step
    even more room and the ring visibly leans towards it; an equal share keeps the
    species as close to evenly spaced as their demands allow. This is the polar form
    of layout_linear's adaptive gutter.
    """
    n = len(demands)
    slack = _BUDGET - sum(demands)
    if slack < 0.0:
        # Only reachable with an author-pinned radius -- otherwise solve_radius
        # would have grown the ring until this was payable.
        if warn is not None:
            warn("the ring's %d steps need %.0f degrees but there are only 360: "
                 "the species are spaced evenly instead and something will "
                 "overlap. Drop 'radius' and let the packer pick one."
                 % (n, sum(demands)))
        return [_BUDGET / n] * n
    share = slack / n
    return [d + share for d in demands]


# ---- feeder motif: reagent in / by-product out, tangent to the ring -------

_MOTIF_KEYS = ("in", "out", "aside")


def _has_motif(step: dict) -> bool:
    """One predicate, asked by the packer and by the drawing loop. Two spellings of
    "does this step have a feeder" would eventually disagree, and the figure would
    be drawn on a budget that never counted it."""
    return any(step.get(k) for k in _MOTIF_KEYS)


def feeder_motif(step: dict, si: int, t0: float, sweep: float, R: float,
                 disk: float, cfg: dict, alloc: Allocator, owner: str, ring_rids,
                 ids: dict, addr: "ADDR.Registry", warn=None,
                 structures: dict | None = None,
                 anchor_src: dict | None = None,
                 edge_as: str | None = None,
                 centre=ORIGIN, interior_ok: bool = False,
                 decorations: dict | None = None,
                 behind: list | None = None,
                 moved: list | None = None) -> FeederMotif:
    """Build one step's `in` / `out` / `aside` on the ring of radius `R`.

    `disk` is the free interior radius -- R less the largest species reach -- and is
    used only to sort the reserved boxes into the ones that cost ring circumference
    and the ones that cost middle-of-the-cycle (see _reach_span). It changes nothing
    about what is drawn.

    Every primitive RESERVES BEFORE IT IS EMITTED. That ordering is the whole
    point: drawing first and accounting afterwards ships overlaps, and a motif that
    meets the ring is exactly the shape of thing that goes wrong that way.

    MODULE LEVEL, and handed its allocator, because the packer runs this very
    function against a throwaway Allocator to find out how much room a feeder wants
    before the radius is settled (render_cycle.probe_feeders). Measured and drawn by
    one routine: a separate "just for the budget" copy of the geometry is a fork
    that lets the ink and the accounting drift out of step.

    `ring_rids` are this step's own arc pieces. A feeder TOUCHES the arc it feeds by
    construction, and in "split" the two feeder arcs touch each other there too, so
    those go on the allow-list. Everything else -- a neighbouring step's arc, a
    structure, the centre note -- is still a real collision and still gets reported.

    `ids` carries the running <t>, <arrow> and fragment ids in and out (keys "t",
    "a" and "f"), so a probe can be handed its own counters and leave the drawing's
    numbering alone.

    `addr` is the address registry (see the one-operation invariant) and is
    REQUIRED, for the same reason `alloc` is: the packer's probe gets a THROWAWAY
    one, exactly as it gets a throwaway Allocator, so the registering wrappers
    below never have to ask whether this run is for real. A conditional there is
    the hole this whole invariant is about -- and a probe runs this function once
    per bisection step, so registering into the live registry would raise on the
    second try.

    `structures` is recipe["structures"], the same {key: (fragment, bbox)} the nodes
    are built from: a feeder entry naming one through "s" DRAWS THAT MOLECULE at the
    free end of its arc instead of a piece of type. Same dict, same convention, so a
    reagent may be a node in one figure and a feeder in the next without being
    rebuilt. Structures go OUTSIDE the ring only -- see _check_steps, which refuses a
    structural by-product on the split style rather than quietly restyling the step.

    `warn` is where a label that could not be pushed clear says so. The PACKER'S
    PROBE PASSES NONE, on purpose: it runs this function a dozen times or more while
    the bisection closes in, and a repair that is reported once per trial radius
    would bury the report in near-duplicates of a complaint about a ring that was
    never drawn.

    That the probe pushes labels at all is what keeps the budget honest, and it works
    because a motif's own boxes go into whichever allocator it was handed -- so the
    in/out clash this repair exists for is visible to the throwaway allocator too,
    and gets repaired there identically before _reach_span measures the result. What
    the probe cannot see is the rest of the page (no species, no neighbouring motif),
    so a label that only gives way at drawing time moves further out than the budget
    was told. That direction is the tangent, mostly away from the ring, and the
    allocator reports the outcome either way.
    """
    prims: list = []
    res: list = []
    # A give_way that finds room MOVED something, and that has to be sayable.
    moved = moved if moved is not None else []
    allow = list(ring_rids)
    step_at = float(step.get("at", 0.5))
    S = structures or {}
    ids.setdefault("f", 0)      # fragment id base; a probe may arrive without one

    def no_room(what: str, cap: float) -> None:
        """Stays at the measured spot, and says so. The overlap itself is the
        allocator's to report -- this names the step and the way out. The cap is
        quoted because it is the box's own size (give_way_cap), so "it could not
        get clear in 2 bond lengths" would be a number from the wrong figure."""
        if warn is not None:
            warn("steps[%d]: '%s' has nowhere to go: sliding it up to %.1f bond "
                 "lengths along its own tangent never got it clear, so it stays "
                 "where the geometry puts it and overlaps. Give this entry its "
                 "own 'at', or shorten it."
                 % (si, what[:60], cap / C.BOND_LENGTH))

    def claim_text(lines, cx_, cy_, w_, h_, sub, *, addr_as, push=None):
        """Reserve a text box, register it, and write the type into it.

        `addr_as` is REQUIRED (the one-operation invariant): a written species IS
        an addressable object, an `aside` and a drawn species' number are not,
        and this signature is where that gets decided rather than at some later
        line that can be omitted.

        `push`, when given, is the unit direction the box may GIVE WAY along if its
        measured spot is already taken -- for a feeder label that is the tangent it
        was set out along in the first place, so giving way means "further out the
        way you were already going", never sideways into a new relationship with the
        figure. In "split" the reagent's tangent and the by-product's are about 132
        degrees apart (measured), so the two of them backing off along their own
        tangents pulls the pair apart, which is why a fraction of a bond length is
        usually enough.

        `at_offset` is the single piece of geometry: the rect the search probes, the
        rect that gets reserved and the anchor the type is written from all come out
        of one call to it. At offset 0 the arithmetic is bit-for-bit the unpushed
        placement.
        """
        o = "%s.%s" % (owner, sub)
        ux, uy = push if push else (0.0, 0.0)

        def at_offset(e: float):
            px, py = cx_ + ux * e, cy_ + uy * e
            return (Rect(px - w_ / 2.0, py - h_ / 2.0, w_, h_),
                    (px, py - h_ / 2.0 + ASCENT_EM * cfg["cap_size"]))

        e = 0.0
        if ux or uy:
            cap = give_way_cap(w_, h_, ux, uy)
            got = give_way(lambda x: at_offset(x)[0], alloc, allow, cap,
                           GIVE_WAY_STEP_BOND * C.BOND_LENGTH,
                           LL.block_leading(lines, cfg["cap_size"]),
                           owner=o, moved=moved)
            if got is None:
                no_room(" ".join(lines), cap)
            elif got:
                e = got
        rect, anchor = at_offset(e)
        r = alloc.reserve(rect, o, "caption")
        alloc.check_emit(r.rid, o)
        res.append(r)
        _file(addr, addr_as, rect, o)
        prims.append(C.styled_text(
            ids["t"], anchor[0], anchor[1],
            lines, size=cfg["cap_size"], align="center"))
        ids["t"] += 1

    # A DRAWN species that had to give way leaves its displacement here, keyed by
    # `sub`. The arc that delivers it is emitted BEFORE the species is placed --
    # it has to be, because the species' spot is measured off the arc's own free
    # end -- so the arc cannot know at drawing time that the box will move. This
    # dict is how it finds out afterwards; see `reach`.
    _gave: dict = {}

    def claim_struct(it, sub, ex, ey, ux, uy, gap_key, *, addr_as):
        """Reserve, register and draw a MOLECULE at the free end of an arc.

        `addr_as` is REQUIRED (the one-operation invariant). This is the SECOND
        `C.place_with` in the module, and the first one is the ring node -- exactly
        the kind of pair where two placement sites exist and only one records.
        Neither can be written here without a declaration.

        The structure and its label are one thing to place and two things to reserve:
        `at_offset` returns the box that contains both, which is what the give-way
        search tests, and the two boxes cut out of it at the chosen offset are what
        get reserved and drawn. Their union is inside the tested box by construction
        (the label is centred under the structure and no wider than the combined
        width), so nothing can appear anywhere the search did not test -- while the
        allocator still gets the honest pair of rectangles rather than one block with
        the gap between them declared as occupied.

        THE STRUCTURE IS NOT ROTATED. Only the position travels along the tangent;
        the molecule stays upright, because a rotated structure is not acceptable in
        a chemistry figure whatever it does for the packing. `angle` (below) turns the
        DIRECTION the box is set out along, never the box's contents.
        """
        o = "%s.%s" % (owner, sub)
        frag, bb = S[it["s"]]
        w_, h_ = bb[2] - bb[0], bb[3] - bb[1]
        lbl = it.get("label") or str(it["s"])
        llines, lw, lh = None, 0.0, 0.0
        if lbl:
            llines, (lw, lh) = _fold([lbl], cfg["num_size"], cfg)
        lgap = cfg["label_gap"] if llines else 0.0
        W, H = max(w_, lw), h_ + lgap + lh

        # THE GAP IS MEASURED TO THE NEAR EDGE, not to the centre.
        #
        # FEEDER's tail_gap / head_gap were measured on the reference LABELS and to
        # the anchor in the middle of them, which is harmless for a 16-unit word and
        # useless for a molecule: at 1.50 and 1.00 bond the arrow's free end lands
        # INSIDE the box of every fragment tried -- phenylacetylene (50 x 36),
        # benzoic acid (25 x 50), the pinacol boronate (97 x 60) -- i.e. the arrow
        # stops in the middle of the structure it is delivering. So a structure asks
        # instead for the clearance this engine already keeps between an arrow and a
        # drawn species, cfg["arc_margin"], measured from the edge of the box along
        # the direction it is set out: 37.0 / 37.2 / 60.4 units for those three,
        # against the 21.6 and 14.4 the centre rule gave.
        #
        # THE EQUIVALENT TEXT CONSTANT, back-derived, is why text is left alone. Read
        # the four reference labels the same way -- free end to the near edge of the
        # box, along the placement direction -- and they come to -0.14, +0.03, -0.30
        # and +0.50 bond lengths: a label STARTS WHERE THE ARC ENDS, straddling the
        # free end as often as not, because a text box is a loose flat abstraction
        # whose near corner is off to one side of the word. Any positive edge
        # clearance would therefore move all four, and every text feeder drawn with
        # them; a clearance of about zero is not a rule worth writing. The
        # centre measurement already says the same thing for a small box and says it
        # in the units it was measured in, so it stays, and it stays as a FLOOR here
        # too -- a fragment small enough (a lone metal atom) not to reach the edge
        # rule should sit where a piece of type would have.
        #
        # The two rules disagreeing is the point, not a wart: a text box under-states
        # its ink and reads as beginning at the arc's end, a fragment's bbox is tight
        # around ink an arrow must not be drawn into. The allocator draws the same
        # line between the two kinds, for the same reason (Policy.text_kinds).
        d = max(FEEDER[gap_key] * C.BOND_LENGTH,
                cfg["arc_margin"] + _support(W, H, ux, uy))

        def at_offset(e: float):
            px, py = ex + ux * (d + e), ey + uy * (d + e)
            return Rect(px - W / 2.0, py - H / 2.0, W, H)

        cap = give_way_cap(W, H, ux, uy)
        # pad=0: a fragment's bbox comes off its atom coordinates, so unlike a text
        # box it does not under-state its own ink and needs no clearance added to it.
        got = give_way(at_offset, alloc, allow, cap,
                       GIVE_WAY_STEP_BOND * C.BOND_LENGTH, kind="structure",
                       owner=o, moved=moved)
        if got is None:
            no_room(lbl or it["s"], cap)
        box = at_offset(got or 0.0)
        if got:
            # The arc ends at (ex, ey) and this box now sits `got` further out along
            # (ux, uy). Recorded rather than fixed up here: this function does not
            # know which arc delivered it, nor which end of that arc carries the head.
            _gave[sub] = (float(got), ex, ey, ux, uy)

        sx, sy = box.cx - w_ / 2.0, box.y
        r = alloc.reserve(Rect(sx, sy, w_, h_), o, "structure")
        alloc.check_emit(r.rid, o)
        res.append(r)
        # The decoration goes to the caller's `behind` list, not into `prims`.
        #
        # Per-motif ordering -- decoration immediately before its own structure,
        # on the argument that a feeder species stands alone outside the ring so
        # nothing else is under its glow -- does not work, and the reason is worth
        # keeping: `prims` ends up in `feeders`, which render_cycle joins LAST,
        # while stamp_z only runs when `behind` is non-empty. So the feeder's glow
        # would never be stamped at all, and an unstamped page graphic is the one
        # case measured to come out ON TOP. Document order only binds if something
        # stamps it.
        _ffx, _fdx = C.place_with(frag, (decorations or {}).get(it.get("s")),
                                  ids["f"], sx, sy)
        if _fdx and behind is not None:
            behind.append(_fdx)
        # The one-operation invariant: registered HERE, in the same breath as the
        # append that draws it and after place_with, so the handle can carry the
        # PLACED string rather than the pre-place one. Nothing between these two
        # lines can fail.
        _file(addr, _with_handle(addr_as, it.get("s"), _ffx, ids["f"],
                                 (anchor_src or {}).get(it.get("s"))),
              Rect(sx, sy, w_, h_), o)
        prims.append(_ffx)
        ids["f"] += 1000
        if llines:
            o_l = o + ".label"
            ly = sy + h_ + lgap
            rl = alloc.reserve(Rect(box.cx - lw / 2.0, ly, lw, lh), o_l, "label")
            alloc.check_emit(rl.rid, o_l)
            res.append(rl)
            prims.append(C.styled_text(
                ids["t"], box.cx, ly + ASCENT_EM * cfg["num_size"],
                llines, size=cfg["num_size"], bold=True, align="center"))
            ids["t"] += 1

    def claim_curve(xml, sub, *, addr_as):
        """One arc of a feeder motif. `addr_as` is REQUIRED (the one-operation
        invariant).

        Every caller passes NotAddressed, and that is a deliberate decision rather
        than an oversight: a feeder arc is one stroke of a motif that the STEP owns,
        it carries no id in any recipe, and the `edge` table's rows mean the ring's
        own steps. Making it an object would invent an id no caller
        wrote and no caller can point at.
        """
        o = "%s.%s" % (owner, sub)
        pts = C.arc_points(xml, n=cfg["arc_samples"])
        if not pts:
            # unreadable as an arc: fall back to reserving the chord, so the
            # extent stays honest instead of silently missing
            pts = [(0.0, 0.0)]
        rs = alloc.reserve_polyline(pts, o, "arrow", allow=allow)
        if rs:
            alloc.check_emit(rs[-1].rid, o)
            res.extend(rs)
            # the next curve in this motif may cross this one
            allow.extend(r.rid for r in rs)
        # Filed even though every caller declares NotAddressed, so the declaration
        # is not decorative: the argument reaches the registry, and the day a
        # feeder arc becomes addressable this line already works. Requiring an
        # argument and then ignoring it is the same hole one level up -- caught by
        # a test that asserts each wrapper reaches _file.
        _file(addr, addr_as, _ink_rect(pts), o)
        prims.append(xml)

    def side(items, tag: str) -> list:
        """Validate one side's entries. Nothing is drawn until both sides have
        passed, so a bad `out` cannot leave half a motif reserved.

        Keys were already whitelisted by _check_steps before any of this ran, along
        with the two questions that need the step's style to answer (text-vs-structure
        exclusivity, and a structure on split's inward arc); what is left here is what
        a well-named key can still say wrong.
        """
        out = []
        for k, it in enumerate(items or []):
            where = "steps[%d].%s[%d]" % (si, tag, k)
            if it.get("s") is not None:
                if it["s"] not in S:
                    raise ValueError(
                        "%s: 's' is %r, which is not a key in recipe['structures'] "
                        "(it has: %s). A feeder names a fragment this engine was "
                        "handed, exactly as a node does."
                        % (where, it["s"], ", ".join(sorted(map(str, S))) or "nothing"))
            elif not it.get("text"):
                raise ValueError("%s: needs a 'text' or an 's'" % where)
            a_t = float(it.get("at", step_at))
            if not 0.0 < a_t < 1.0:
                raise ValueError(
                    "%s: 'at' must be strictly between 0 and 1, got %r"
                    % (where, a_t))
            # PINNED means this entry named its own tangency point and must not
            # be fanned off it -- that is how an author pulls an in/out pair
            # apart into two separate meetings with the ring.
            out.append((it, a_t, it.get("at") is not None))
        return out

    style = str(step.get("feeder_style") or cfg["feeder_style"])
    if style not in ("split", "single"):
        raise ValueError(
            "steps[%d]: 'feeder_style' is 'split' or 'single', got %r"
            % (si, style))
    ins, outs = side(step.get("in"), "in"), side(step.get("out"), "out")
    # `single` takes ONE `out` and AT MOST one `in`.
    #
    # `split` refuses a structural `out`, so demanding exactly one of each here
    # would block "the product peels off this step and there is no co-reactant"
    # on both sides. That shape recurs independently in published cycles, which
    # is what makes it a need rather than a preference.
    #
    # With no `in` the arc is unchanged and its tangency comes from the `out`
    # instead: one arc still straddles the tangency point, and only the tail label
    # is absent. A figure that DOES give an `in` takes exactly the one-in-one-out
    # path, tangency included, so its bytes do not change.
    if style == "single" and not (len(ins) <= 1 and len(outs) == 1):
        raise ValueError(
            "steps[%d]: feeder_style 'single' is ONE arc whose head is the "
            "product, so it needs exactly one 'out' and at most one 'in' (got "
            "%d and %d). Use 'split', which draws an arc per entry."
            % (si, len(ins), len(outs)))

    def emit(Cf, r, a_tail, a_head, head, sub):
        """One tangent arc, reserved then drawn.

        The SAME primitive the ring is drawn with, deliberately: a feeder is a
        circular arc, so arc_arrow already writes it exactly, arc_points already
        reads it back for the reservation, and ChemDraw gets one kind of arc to
        round-trip instead of two.
        """
        _aid = ids["a"]
        xml = C.arc_arrow(_aid, Cf[0], Cf[1], r, r, a_tail, a_head - a_tail,
                          head=head)
        ids["a"] += 1
        claim_curve(xml, sub, addr_as=NotAddressed(
            "a stroke of the step's own motif, not an object: no recipe gives it an "
            "id"))
        # The handle `reach` needs, and nothing else: WHERE in `prims` this arc
        # landed (claim_curve appends it last), and enough to write it again
        # without its head. Every existing caller ignores the return value.
        return (len(prims) - 1, _aid, Cf, r, a_tail, a_head, head)

    def reach(sub, handle, head_here: bool):
        """Carry the ink on to a structure that gave way from this arc's free end.

        WHY NOT GROW THE ARC. `feeder_radius` is a clamp measured off two reference
        drawings and `_feeder_circle` has no free parameters by construction -- a
        feeder is a fixed-size decoration tangent to the ring, not a rubber band. So
        the arc keeps its shape and carries on in the direction it was already
        pointing, which is exactly the direction the box was pushed: `give_way`
        moves a box ALONG ITS OWN TANGENT at this end (see claim_struct), so a
        straight segment of that same length continues the same line and leaves the
        species the `arc_margin` clearance it would have had if it had never moved.

        WHY THE HEAD MOVES WITH IT. An arrowhead left in the middle of a stroke that
        carries on to a molecule reads as two arrows, the first pointing at nothing.
        So when the head sat at the end that gave way, the arc is written again
        HEADLESS under its own id and the head goes on the new segment.

        A FIGURE WHERE NOBODY GAVE WAY IS UNTOUCHED: `_gave` is empty, this returns
        immediately, and no id is spent. That is what keeps the regression figures
        byte-exact.

        Measured on an N-SCF3 phthalimide and its anion on one `single` arc:
        without this, the by-product is pushed 46.8 units while the arc stops 46.8
        units short of it -- two molecules floating beside a cycle with no ink
        between them, reported as "the figure is correct".
        """
        g = _gave.get(sub)
        if g is None or handle is None:
            return
        got, ex, ey, ux, uy = g
        idx, aid_used, Cf, r, a_tail, a_head, head = handle
        x1, y1 = ex + ux * got, ey + uy * got
        if head_here and head:
            prims[idx] = C.arc_arrow(aid_used, Cf[0], Cf[1], r, r, a_tail,
                                     a_head - a_tail, head=None)
            xml = C.arrow_xy(ids["a"], ex, ey, x1, y1)
        else:
            xml = C._line(ids["a"], ex, ey, x1, y1)
        ids["a"] += 1
        o = "%s.%s.reach" % (owner, sub)
        rs = alloc.reserve_polyline([(ex, ey), (x1, y1)], o, "arrow", allow=allow)
        if rs:
            alloc.check_emit(rs[-1].rid, o)
            res.extend(rs)
            allow.extend(_r.rid for _r in rs)
        _file(addr, NotAddressed(
            "the tail the arc grew to reach its own species: one stroke of the "
            "motif, with no id in any recipe"),
            _ink_rect([(ex, ey), (x1, y1)]), o)
        prims.append(xml)
        if warn is not None:
            warn("steps[%d]: '%s' gave way %.1f units, so its feeder arc was "
                 "carried on by a straight segment of the same length to reach "
                 "it%s." % (si, sub, got,
                            " (the arrowhead moved to the new end)"
                            if head_here and head else ""))

    def label(Cf, r, alpha, away, it, gap_key, sub):
        """What sits at the free end `alpha`: a word, or a molecule, set out along
        the TANGENT there.

        Out along the tangent, not out along the radius: that puts the species on
        the continuation of the arc's own line, so reagent, arc and ring read as
        one stroke of flow. The single reference drawing's tail label sits 0.1 deg
        off this direction.

        `angle`, when the entry carries one, ROTATES THAT DIRECTION -- the manual
        nudge, for the cases an author can see and the geometry cannot (a published
        figure may set a reagent at a slant because its author chose to). It turns
        where the box goes and nothing about the box: there is no automatic
        rotation and no automatic dodging here, on purpose.

        FOLDED here, at the one place a label's box is worked out, so the lines
        that go to claim_text are the lines that were measured. A reagent name is
        the longest piece of type on a cycle -- these are what overlapped -- and a
        folded one is a taller box at the same centre, which the gap below already
        measures through _support and the packer already charges for through
        _box_disk. Nothing else in the motif needs to know.
        """
        ex, ey = _circle_pt(Cf, r, alpha)
        ux, uy = _circle_tangent(alpha, away)
        turn = float(it.get("angle") or 0.0)
        if turn:
            c_, s_ = math.cos(math.radians(turn)), math.sin(math.radians(turn))
            ux, uy = ux * c_ - uy * s_, ux * s_ + uy * c_
        # Named as a node is: the molecule takes the entry's own name and its caption
        # takes ".label". A TEXT entry keeps the original name -- the type is the
        # whole of that entry, and ".label" is what layout reports and tests already
        # call it -- so ".label" always means "type belonging to this entry" either
        # way.
        # THE ENTRY'S OWN `as` IF IT HAS ONE, else composed off whatever this step
        # is called: `<edge name>.in[0]`. So naming the step renames its feeders
        # with it, and naming one entry overrides just that one.
        #
        # It cannot be the entry's "id": that is a key into `structures`, and the
        # tool's own docs promise a reagent drawn as a feeder and as a node "means
        # the same molecule built once" -- one key, two instances, which the
        # address layer forbids as one address.
        #
        # It cannot be `owner` either: that is "step:<key i>-><key j>", so a ring
        # whose nodes share a structures key gives SEVERAL steps the same owner.
        # Deriving from it would raise `duplicate id 'step:s1->s1.in[0]'` on a
        # three-node ring of one repeated species -- a test builds exactly that.
        # The edge's name comes from the edge TABLE, which is the single
        # definition of it; recomputing it from `step` here would be a second
        # one that merely happens to agree. `si` remains the fallback for a caller
        # that builds a motif without a table -- the packer's trial run does
        # exactly that.
        _base = edge_as or _as_id(step.get("as"), "steps[%d]" % si,
                                  "edge#%d" % si)
        _aid = _as_id(it.get("as"), "steps[%d].%s" % (si, sub),
                      "%s.%s" % (_base, sub))
        if it.get("s") is not None:
            claim_struct(it, sub, ex, ey, ux, uy, gap_key,
                         addr_as=Addressed(_aid, "structure"))
            return
        lines, (lw, lh) = _fold([it["text"]], cfg["cap_size"], cfg)
        d = FEEDER[gap_key] * C.BOND_LENGTH
        # the gap direction IS the give-way direction: if this spot is taken the
        # label carries on out along the same line rather than looking for a new one
        claim_text(lines, ex + ux * d, ey + uy * d, lw, lh,
                   sub + ".label", push=(ux, uy),
                   # A WRITTEN species is a `text` object, exactly as a written node
                   # is: text nodes are the objects published figures most want as
                   # an arrow endpoint. A reagent written at an arc's end is the
                   # same kind of thing.
                   addr_as=Addressed(_aid, "text"))

    if style == "split":
        # in[k] and out[k] SHARE a tangency point, so the pair is one unbroken
        # S of ink crossing the ring: the reagent arc arrives at the tangency
        # HEADLESS and the by-product arc sets off from that same point carrying
        # the only arrowhead -- one head between them, which is what makes the
        # pair read as one gesture. Successive pairs fan along the arc, or two
        # reagents entering at one point would stack on top of each other.
        for k in range(max(len(ins), len(outs))):
            fanned = min(0.95, max(0.05, step_at + k * FEEDER["fan"]))
            if k < len(ins):
                it, a_t, pinned = ins[k]
                a_t = a_t if pinned else fanned
                r = feeder_radius(R, C.BOND_LENGTH, FEEDER["r_in"])
                Cf, aT, dirn = _feeder_circle(t0, sweep, R, a_t, r, True, centre)
                a_free = aT + dirn * FEEDER["sweep_in"]     # upstream
                _h = emit(Cf, r, a_free, aT, None, "in[%d].curve" % k)
                label(Cf, r, a_free, dirn, it, "tail_gap", "in[%d]" % k)
                # The free end is the arc's tail here and carries no head.
                reach("in[%d]" % k, _h, head_here=False)
            if k < len(outs):
                it, a_t, pinned = outs[k]
                a_t = a_t if pinned else fanned
                r = feeder_radius(R, C.BOND_LENGTH, FEEDER["r_out"])
                Cf, aT, dirn = _feeder_circle(t0, sweep, R, a_t, r, False, centre)
                a_free = aT + dirn * FEEDER["sweep_out"]    # downstream
                _h = emit(Cf, r, aT, a_free, "Full", "out[%d].curve" % k)
                label(Cf, r, a_free, dirn, it, "head_gap", "out[%d]" % k)
                # A structure cannot be split's `out` (it is refused, because that
                # arc dives inside the ring), so this call is for the day that
                # changes -- and it costs nothing while `_gave` stays empty.
                reach("out[%d]" % k, _h, head_here=True)
    else:
        # One externally tangent arc STRADDLING the tangency point -- it touches
        # the ring at the middle of its sweep and both ends stay outside, which
        # is why both labels can sit out there clear of the cycle. The arc is the
        # `in` and the `out` at once, so its tangency comes from the `in`.
        it_out, _a_out, _p = outs[0]
        # The tangency comes from the `in` when there is one -- unchanged, so a
        # figure with both entries is byte-for-byte the same -- and from the
        # `out` when there is not.
        it_in, a_t = (ins[0][0], ins[0][1]) if ins else (None, _a_out)
        r = feeder_radius(R, C.BOND_LENGTH, FEEDER["r_single"])
        Cf, aT, dirn = _feeder_circle(t0, sweep, R, a_t, r, True, centre)
        half = FEEDER["sweep_single"] / 2.0
        # WITH NO `in`, THE ARC IS HALF OF ONE. The straddling arc's two ends both
        # sit 14.4 units OUTSIDE the ring (measured) -- which is right when each end
        # carries a species, because the pair reads as "this goes in, that comes
        # out". With no reagent the tail end carries nothing, and a bare arc end
        # hanging in the air outside the circle reads as a stray tail with nothing
        # attached. So the tail starts AT the tangency, on the ring, and the arc
        # sweeps out to the product -- which is what "the product peels off this
        # step" draws.
        a_tail = (aT + dirn * half) if it_in is not None else aT
        a_head = aT - dirn * half
        _h = emit(Cf, r, a_tail, a_head, "Full", "single.curve")
        if it_in is not None:
            label(Cf, r, a_tail, dirn, it_in, "tail_gap", "in[0]")
        label(Cf, r, a_head, -dirn, it_out, "head_gap", "out[0]")
        # This is THE case: two drawn species on one arc, the second of which has
        # to give way from the first. The head sits at the `out` end.
        reach("in[0]", _h, head_here=False)
        reach("out[0]", _h, head_here=True)

    aside = step.get("aside")
    if aside:
        lines = [aside] if isinstance(aside, str) else list(aside)
        (px, py), (nx, ny) = _ring_frame(t0, sweep, R, step_at, centre)
        lines, (aw, ah) = _fold(lines, cfg["cap_size"], cfg)
        # split's by-product arc dives INSIDE the ring at the tangency, so the
        # additive goes outside with the reagent; single keeps both arcs outside
        # and leaves the inside free for it. Measured that way round on both
        # references. Outside it also has to clear the ring arc itself, hence
        # the half-extent; inside there is nothing to clear.
        s_out = 1.0 if style == "split" else -1.0
        d = FEEDER["aside_gap"] * C.BOND_LENGTH
        if s_out > 0:
            d += _support(aw, ah, nx, ny)
        claim_text(lines, px + nx * s_out * d, py + ny * s_out * d,
                   aw, ah, "aside",
                   # The tool's own words: "Type only -- it is a reagent note, not
                   # a species." A note about a step, so the step is the object.
                   addr_as=NotAddressed(
                       "a reagent note on the step, not a species; the tool "
                       "documents it as type only"))

    up, down, defined = _reach_span([r.rect for r in res], disk,
                                    t0 + sweep * step_at,
                                    1.0 if sweep >= 0 else -1.0,
                                    centre, interior_ok)
    return FeederMotif(prims, res, (up + down) if defined else None,
                       up, down, defined)


def render_cycle(recipe: dict, header_prefix: str, cfg: dict | None = None,
                 page_mode: str = "fit", alloc: Allocator | None = None,
                 policy: Policy | None = None,
                 report: dict | None = None,
                 shift: list | None = None,
                 addresses: "ADDR.Registry | None" = None) -> str:
    """Render the ring. Returns the CDXML string, as the other engines do.

    `addresses` is an address registry to fill in. Passing one is optional;
    REGISTERING IS NOT -- when the caller does not supply a registry this function
    makes its own and fills it either way, then publishes it on
    `report["addresses"]`. That asymmetry is the invariant: there is no mode in
    which a placement goes unregistered, only a mode in which nobody kept the
    answer. Cost is a dict of rects the packer had already computed, which is why
    it can be unconditional.

    A REGISTRY IS PER FIGURE. Handing the same one to two calls raises on the
    second -- the names are per figure, so the second render re-registers `node#0`.
    Pass a fresh one, or pass none. See Registry's own docstring; it matters for
    multi-ring figures, where sharing one looks obvious.

    `report`, if given, is filled in with anything the caller needs that is not
    part of the drawing -- "warnings", "feeder_extents" and "ring" (what the packer
    decided). Passing a dict rather than changing the return type keeps every
    existing caller (and the regression tests) working unchanged.
    """
    cfg = _cfg(cfg)
    nodes = recipe.get("nodes") or []
    if not nodes:
        raise ValueError("recipe['nodes'] is empty")
    alloc = alloc if alloc is not None else Allocator(policy)

    meas = measure(recipe, cfg)          # checks the nodes' keys as it goes
    steps = list(recipe.get("steps") or [])
    _check_steps(steps, cfg)
    n = len(meas)

    # ---- nodes + steps -> the three tables ------------------------------------
    #
    # HERE, before anything is measured against the ring, because the centre note
    # is the ring's own label and it sizes the ring. From this point the
    # function reads the tables: the note below, the edge loop's endpoints, the
    # per-step payload. Every arithmetic expression is unchanged -- a byte-identity
    # regression check is what proves that.
    G = desugar(recipe, cfg, meas)
    EDGES = G["edges"]
    RINGS = G["rings"]
    # The fourth table (paths). Always present -- the old form emits [] --
    # so nothing downstream needs a `.get`.
    PATHS = G["paths"]
    RING = RINGS[0]
    # THE PACKER SEES THE RING'S MEMBERS, NOT EVERY NODE (a node with an explicit
    # `at` is placed by its coordinates, one without is placed by its ring). A node
    # the caller pinned spends no circumference and must not enter the angle budget
    # or the radius it forces.
    #
    # On the old form the ring IS every node in order, so `RIDX` is [0..n-1] and
    # `RMEAS` is `meas` -- every expression below is then bit-for-bit the
    # single-ring arithmetic, which the byte-identity check verifies.
    RIDX = RING["_idx"]              # node index per ring position
    RPOS = RING["_pos"]              # node index -> ring position
    RMEAS = [meas[i] for i in RIDX]
    n_ring = len(RIDX)
    if n_ring < 2:
        raise ValueError(
            "the ring has %d member(s). A cycle needs at least two species on the "
            "circle; nodes pinned with an explicit 'at' are not on it."
            % n_ring)

    def _ring_note(rg):
        """(lines, w, h, reach) for one ring's CENTRE label, folded ONCE.

        Folded here and the same lines and box are used to size the ring and to
        write the note -- the number that sets the radius has to be the number that
        lands on the page. A label at `above` / `below` is not a centre note: it
        does not compete for the middle, so it contributes no `note_reach`.
        """
        lab = rg.get("label")
        if not lab or lab.get("position", "center") != "center":
            return [], 0.0, 0.0, 0.0
        t = lab.get("text")
        lines = ([t] if isinstance(t, str) else list(t)) if t else []
        if not lines:
            return [], 0.0, 0.0, 0.0
        lines, (w, h) = _fold(lines, cfg["center_size"], cfg)
        return lines, w, h, math.hypot(w / 2.0, h / 2.0)

    note_lines, _nw, _nh, note_reach = _ring_note(RING)
    note = note_lines[0] if note_lines else None

    parts, texts, arcs, feeders = [], [], [], []
    # Halo / dot / locants for the nodes, joined FIRST -- see layout_linear.
    # A feeder's decoration is layered inside its own motif instead; see the note
    # at the C.place_with call in feeder_motif for why that is enough there.
    D = recipe.get("decorations") or {}
    behind: list = []
    warnings: list[str] = []
    extents: dict[str, float | None] = {}
    undefined: list[str] = []
    tid, base, aid = 92000, 1000, 93000
    # Unconditional, per the signature's note: no "if the caller asked for it".
    addr = addresses if addresses is not None else ADDR.Registry()
    # {structures key: (smiles, index_map, base_id)}, put there by the MCP layer at
    # build time -- the two halves of an atom lookup that only IT knows. Read here
    # and stored in the handle UNPARSED: that keeps node_points() and RDKit out of
    # every figure that never asks for an atom. Same recipe key layout_linear uses.
    ASRC = recipe.get("anchor_src") or {}

    # ---- THE PACKER, ONCE PER RING --------------------------------------------
    #
    # The existing angle-packing solver is kept whole and becomes a subroutine
    # run once per ring: `legacy_radius` / `solve_radius` / `allocate_gaps` /
    # `_gap_demands` are untouched, and the loop below calls them per ring
    # instead of once. Each ring carries its own centre, radius, start_angle and
    # direction, so each gets its own answer.
    #
    # A ONE-RING FIGURE TAKES THE IDENTICAL PATH. `RSTATE[ring0]` holds the same
    # numbers the single-ring code computed, produced by the same expressions in the
    # same order -- which the byte-identity regression check verifies.
    def _solve_ring(rg):
        """Everything geometric about one ring. Returns a dict, never mutates."""
        idx = rg["_idx"]
        rmeas = [meas[i] for i in idx]
        nr = len(idx)
        if nr < 2:
            raise ValueError(
                "rings[%r]: %d member(s). A cycle needs at least two species on "
                "the circle; nodes pinned with an explicit 'at' are not on it."
                % (rg["as"], nr))
        centre = tuple(rg.get("centre") or ORIGIN)
        interior_ok = bool(rg.get("_interior_ok"))
        sgn_r = cw_sign() * (1.0 if rg["direction"] == "cw" else -1.0)
        step = 360.0 / nr
        adv_r = sgn_r * step
        s_r = 1.0 if adv_r > 0 else -1.0
        # The species occupy an annulus of this half-width about the ring, so R
        # minus it is the free middle -- what a feeder box has to be wholly inside
        # to count as spending disk rather than circumference.
        rmax = max(m["reach"] for m in rmeas)
        nlines, nw, nh, nreach = _ring_note(rg)

        # EDGE INDICES, not node indices. `e["_from"]` is a NODE index -- equal to
        # the edge index only when the edges are written in ring order, which
        # every old-form recipe is. Measured: on a three-ring written
        # `C->A, A->B, B->C`, building `fed` from node indices DRAWS the feeder on
        # C->A and solves the radius as if it were on B->C. Silent mis-sizing,
        # possible as soon as edge order is free.
        fed_r = [k for k, e in enumerate(EDGES)
                 if e["ring"] == rg["as"] and e["via"] == "arc"
                 and _has_motif(_step_view(e))]

        def probe(R_try):
            """What every feeder on THIS ring wants at radius R_try -- measured by
            BUILDING each one against a throwaway allocator and reading back what it
            reserved. `probe_feeders`' method, per ring.

            The reservations ARE the measurement: one geometry, not a second
            "just for the budget" copy of the same shape.
            """
            need, boxes = {}, []
            for k in fed_r:
                e = EDGES[k]
                cap = _step_view(e)
                pf, pt = rg["_pos"][e["_from"]], e["_to"]
                t0_, sweep_ = _arc_of(rg["start_angle"] + adv_r * pf, adv_r,
                                      meas[e["_from"]]["reach"],
                                      meas[pt]["reach"], R_try, cfg, s_r)
                m = feeder_motif(cap, k, t0_, sweep_, R_try, R_try - rmax,
                                 cfg, Allocator(), "probe:%d" % k, (),
                                 {"t": 0, "a": 0, "f": 0}, ADDR.Registry(),
                                 structures=recipe.get("structures"),
                                 edge_as=e["as"], centre=centre,
                                 interior_ok=interior_ok)
                need[k] = _arc_need(m, float(cap.get("at", 0.5)), cfg)
                boxes.extend(r.rect for r in m.reservations)
            return need, boxes

        rad = rg.get("radius")
        pinned_r = rad if isinstance(rad, (int, float)) else None
        trace_r: dict = {}
        lo_r = legacy_radius(rmeas, cfg, nreach)
        if pinned_r:
            R_r = float(pinned_r)
            trace_r.update(solver="pinned", tries=0, need=None)
        else:
            R_r = solve_radius(rmeas, cfg, nreach,
                               probe if fed_r else None, trace_r)

        if fed_r:
            gaps_r = allocate_gaps(_gap_demands(rmeas, cfg, R_r, probe(R_r)[0]),
                                   warnings.append)
            ang_r, acc = [], 0.0
            for g in gaps_r:
                ang_r.append(rg["start_angle"] + sgn_r * acc)
                acc += g
            advs_r = [sgn_r * g for g in gaps_r]
        else:
            # THE RED LINE. No feeder anywhere -> the pre-packer ring, to the bit.
            # Written out rather than routed through allocate_gaps because equal
            # demands come to 360/n only in real arithmetic, and because 270 + adv*3
            # is not 270+adv+adv+adv in floating point -- differences far below the
            # two-decimal grid, but the guarantee being tested is byte equality.
            gaps_r = [step] * nr
            ang_r = [rg["start_angle"] + adv_r * i for i in range(nr)]
            advs_r = [adv_r] * nr
        return {
            "rg": rg, "idx": idx, "pos": rg["_pos"], "meas": rmeas, "n": nr,
            "centre": centre, "interior_ok": interior_ok,
            "sgn": sgn_r, "s": s_r, "adv": adv_r, "step_deg": step,
            "reach_max": rmax, "R": R_r, "lo": lo_r, "ang": ang_r,
            "advs": advs_r, "gaps": gaps_r, "trace": trace_r,
            "note_lines": nlines, "nw": nw, "nh": nh, "note_reach": nreach,
            "fed": fed_r,
        }

    RSTATE = {rg["as"]: _solve_ring(rg) for rg in RINGS}

    # ---- SEVERAL RINGS AND NOBODY SAID WHERE: put them side by side -----------
    #
    # A ring with no `centre` of its own gets ORIGIN, which is right for one ring
    # and is wrong for two: a photoredox / nickel dual cycle would come back with
    # the Ni species drawn at the IDENTICAL coordinates as the Ir ones -- Ir(II)
    # and Ni(II) both at (235, 180) -- and 113 declared collisions. Giving each
    # ring a `centre` by hand produces two rings side by side with nothing
    # overlapping at all; this supplies that as the default.
    #
    # WHY IT HAS TO HAPPEN HERE and not where the centres are normalised: a ring's
    # own extent is R + reach_max, and R is what the packer has just solved. Before
    # this line there is no radius to space them by.
    #
    # ONLY WHEN NOBODY CHOSE. One explicit centre means the caller is placing the
    # rings themselves -- nested, stacked, wherever -- and a figure that moves
    # under them is worse than one that overlaps where they put it.
    #
    # SPACING IS cfg["center_gap"], not a new number: its own comment reads
    # "centre note must stay this clear of the species", which is the same
    # question asked of two rings' species instead of a note and a ring's.
    # ---- ...OR A RING HUNG OFF ANOTHER ONE, TOUCHING IT ----------------------
    #
    # WHY THIS AND NOT THE SPACING ABOVE. The default puts rings clear of each
    # other by their OUTER extent (R + reach_max + center_gap), which is the right
    # answer to "do not overlap" and the wrong answer to "how do published figures
    # draw a dual cycle". Measured on published figures: most two-ring cycles draw
    # the circles TANGENT and print the SET / PT / HAT badge on the touch point --
    # the tangency is where the electron transfer is, so it is
    # chemistry, not spacing. The default spread leaves the circles 37 units apart
    # even with single-letter species, and no `center_gap` can express "touching",
    # because the number that would is R1 + R2 and this loop is the only place that
    # knows both.
    #
    # ROTATION COMES WITH IT. Tangency puts the closest approach at ONE point; a
    # species sitting there lands its label on the other ring (measured: 14.3 and
    # 18.0 units of block-on-block on two tangent pairs). So each
    # ring whose rotation the caller did NOT pin is turned by rotate_to_face to
    # meet the other one arc-first. A stated `start_angle` is left alone -- see
    # `_start_given`.
    _tan = [rg for rg in RINGS if rg.get("tangent_to")]
    if _tan:
        _known = {rg["as"] for rg in RINGS if not rg.get("tangent_to")}
        for rg in _tan:
            if rg["tangent_to"] not in {r["as"] for r in RINGS}:
                raise ValueError(
                    "rings[%r].tangent_to: no ring called %r (this figure has %s)."
                    % (rg["as"], rg["tangent_to"],
                       ", ".join(sorted(r["as"] for r in RINGS))))
            if rg["tangent_to"] == rg["as"]:
                raise ValueError(
                    "rings[%r].tangent_to: a ring cannot be tangent to ITSELF."
                    % rg["as"])
        _left = list(_tan)
        while _left:
            _moved = False
            for rg in list(_left):
                if rg["tangent_to"] not in _known:
                    continue
                _me, _ref = RSTATE[rg["as"]], RSTATE[rg["tangent_to"]]
                _th = math.radians(rg["tangent_at"])
                _d = _me["R"] + _ref["R"]
                _me["centre"] = (_ref["centre"][0] + _d * math.cos(_th),
                                 _ref["centre"][1] + _d * math.sin(_th))
                rg["centre"] = _me["centre"]
                # Turn each ring to face the other ARC-FIRST. `ang` holds the slot
                # angles and every later pass reads it (never start_angle again), so
                # a rotation is one delta added to that list.
                for _st, _rg, _face in ((_me, rg, rg["tangent_at"] + 180.0),
                                        (_ref, _ref["rg"], rg["tangent_at"])):
                    if _rg.get("_start_given"):
                        continue
                    _n = len(_st["ang"])
                    _new = rotate_to_face(_rg["start_angle"], _face, _n,
                                          _rg["direction"])
                    _dl = _new - _rg["start_angle"]
                    if _dl:
                        _st["ang"] = [a + _dl for a in _st["ang"]]
                        _rg["start_angle"] = _new
                # THE TWO ARCS THAT MEET THERE TOUCH BY CONSTRUCTION, so they are
                # put on each other's allow-list rather than reported. Same rule as
                # the feeder's ("A feeder TOUCHES the arc it feeds by construction")
                # and narrowed the same way: ONE arc per ring -- the one whose span
                # crosses the tangency direction -- not the whole ring. The wide
                # version would hide a second ring crossing this one somewhere else,
                # which is the mistake a too-wide feeder exemption would make.
                #
                # It is the ENGINE's contact here, not the caller's, which is why it
                # is exempted and not left for `crossing: {"crossed": ...}`: the
                # caller never wrote the crossing and cannot be asked to declare it.
                _pair = []
                for _st, _rgd, _face in ((_me, rg, rg["tangent_at"] + 180.0),
                                         (_ref, _ref["rg"], rg["tangent_at"])):
                    _pick = None
                    for _e in EDGES:
                        if _e.get("_shape") != "arc" or _e.get("ring") != _rgd["as"]:
                            continue
                        _pp = _st["pos"].get(_e["_from"])
                        if _pp is None:
                            continue
                        _mid = _st["ang"][_pp] + _st["advs"][_pp] / 2.0
                        _dd = abs(((_face - _mid + 180.0) % 360.0) - 180.0)
                        if _pick is None or _dd < _pick[0]:
                            _pick = (_dd, _e)
                    if _pick is not None:
                        _pair.append(_pick[1])
                if len(_pair) == 2:
                    for _one, _other in (_pair, _pair[::-1]):
                        _one["_touch_set"] = tuple(
                            set(_one["_touch_set"]) | {_other["as"]})
                _known.add(rg["as"])
                _left.remove(rg)
                _moved = True
            if not _moved:
                # Every survivor is waiting on another survivor: a ring of
                # `tangent_to`s with no ring standing still to hang it off.
                raise ValueError(
                    "rings: 'tangent_to' goes round in a circle (%s) -- one ring "
                    "has to be placed on its own, by `centre` or by being the one "
                    "nobody hangs off."
                    % ", ".join("%s->%s" % (r["as"], r["tangent_to"])
                                for r in _left))
    elif len(RINGS) > 1 and not any(rg.get("_centre_given") for rg in RINGS):
        _x = 0.0
        for rg in RINGS:
            _st = RSTATE[rg["as"]]
            _half = _st["R"] + _st["reach_max"]
            _st["centre"] = (_x + _half, 0.0)
            rg["centre"] = _st["centre"]
            _x += 2.0 * _half + cfg["center_gap"]

    _S0 = RSTATE[RING["as"]]
    # Single-ring aliases: every consumer below that is not per-ring reads ring 0
    # through these.
    sgn, step_deg, adv, s = _S0["sgn"], _S0["step_deg"], _S0["adv"], _S0["s"]
    reach_max = _S0["reach_max"]
    R, lo, ang, advs, gaps, trace = (_S0["R"], _S0["lo"], _S0["ang"],
                                     _S0["advs"], _S0["gaps"], _S0["trace"])

    def put(xml_list, xml, rect, owner, kind, *, addr_as, addr_where=None):
        """Reserve, register, draw -- ONE operation (the one-operation invariant).

        `addr_as` is a REQUIRED keyword and must be `Addressed(...)` or
        `NotAddressed(why)`. Required rather than defaulted, because a default is
        precisely the shape of the hole where registration is a separate optional
        line and one of two placement paths simply lacks it. Here Python refuses
        the call.

        `addr_where` is the recipe path for the record's provenance, defaulting to
        the allocator owner. Given separately where the owner has lost something
        the author would want back -- a node's owner is `node:<structures key>`,
        which two nodes can share, while `nodes[3]` cannot.
        """
        res = alloc.reserve(rect, owner, kind)
        alloc.check_emit(res.rid, owner)
        _file(addr, addr_as, rect, addr_where or owner)
        xml_list.append(xml)
        return res

    def put_edge(xml, pts, owner, eid, *, check, allow=()):
        """One EDGE's ink -- a ring arc or a straight line: reserve, register, draw.

        Its own wrapper because it is the ONE placement in this function that does
        not go through `put` -- i.e. exactly where a bypass would live. The
        check_emit runs only when `check` is set (the readable-arc branch).

        No `addr_as` here, and that is not a weaker rule but a stronger one: this
        wrapper places exactly one kind of thing, so there is nothing to decide and
        registration is unconditional. THE INK IS THE POLYLINE IT RESERVED -- one
        `pts`, used for both -- so an edge's `at` can never be measured along a
        curve different from the one the ledger accounted. Same discipline as
        "The reservations ARE the measurement".
        """
        rs = alloc.reserve_polyline(pts, owner, "arrow", allow=allow)
        if check:
            alloc.check_emit(alloc.reservations()[-1].rid, owner)
        _file(addr, Addressed(eid, "edge",
                              ink=tuple((float(x), float(y)) for x, y in pts)),
              _ink_rect(pts), owner)
        arcs.append(xml)
        return rs

    def free(rect, kind="caption"):
        return not alloc.probe(rect, kind=kind)

    def at(t_deg, r, centre=ORIGIN):
        t = math.radians(t_deg)
        return centre[0] + math.cos(t) * r, centre[1] + math.sin(t) * r

    def _line_edge_caption(e, cap, pa, pb):
        """`above` / `below` on a straight edge: beside its midpoint.

        PERPENDICULAR to the line, which is what "above" means on something that is
        not horizontal -- `above` goes to one side and `below` to the other, and on
        a horizontal edge that reduces to the usual up / down. Tried at increasing
        offsets like the arc captions are, then placed anyway and reported: a
        silent overlap is the one outcome that must not happen.

        This is deliberately the simple version -- no sliding ALONG the line (an
        arc caption does that to dodge a feeder, and a line edge has no feeder).
        """
        nonlocal tid
        for key, side in (("above", -1.0), ("below", +1.0)):
            lines = cap.get(key)
            if not lines:
                continue
            lines, (cw_, ch_) = _fold(lines, cfg["cap_size"], cfg)
            mx, my = (pa[0] + pb[0]) / 2.0, (pa[1] + pb[1]) / 2.0
            dx, dy = pb[0] - pa[0], pb[1] - pa[1]
            d = math.hypot(dx, dy) or 1.0
            nx, ny = -dy / d, dx / d          # unit normal
            px = py = None
            for extra in (0.0, ch_, ch_ * 2.0, ch_ * 3.0):
                gap = cfg["cap_gap"] + _support(cw_, ch_, nx, ny) + extra
                tx, ty = mx + nx * side * gap, my + ny * side * gap
                if free(Rect(tx - cw_ / 2.0, ty - ch_ / 2.0, cw_, ch_)):
                    px, py = tx, ty
                    break
            if px is None:
                gap = cfg["cap_gap"] + _support(cw_, ch_, nx, ny)
                px, py = mx + nx * side * gap, my + ny * side * gap
                warnings.append(
                    "edges[%s]: nothing beside the line is clear for its '%s' "
                    "caption; it is drawn at the conventional spot and overlaps. "
                    "Move an endpoint, or shorten the caption." % (e["as"], key))
            put(texts, C.styled_text(tid, px, py - ch_ / 2.0
                                     + ASCENT_EM * cfg["cap_size"],
                                     lines, size=cfg["cap_size"], align="center",
                                     color=(cfg["muted_color"] if e["muted"]
                                            else None)),
                Rect(px - cw_ / 2.0, py - ch_ / 2.0, cw_, ch_),
                "edge:%s.caption" % e["as"], "caption",
                addr_as=NotAddressed(
                    "the edge's own caption; address the edge"))
            tid += 1

    # ---- nodes: on the ring by angle, or wherever an explicit `at` says -------
    #
    # "Inside the ring / outside the ring" is NOT a property of a node. A node
    # with `at` sits at those coordinates; a node without one is placed by its
    # ring. A species at the centre of the circle and an off-ring chain are the
    # same thing in the model, at different coordinates -- so there is one loop,
    # not two.
    #
    # Old form: every node is a ring member, so every iteration takes the `at(ang
    # [...])` branch.
    for i, m in enumerate(meas):
        _pinned = G["nodes"][i]["at"]
        # A node is placed by ITS OWN ring (there may be several), or by its
        # explicit `at`. A SHARED node -- one listed by two rings -- is
        # placed by the FIRST ring that lists it; `_ring` records which, so the
        # choice is made once in desugar and not re-derived here.
        _st = RSTATE.get(G["nodes"][i]["_ring"]) if _pinned is None else None
        if _pinned is None:
            cx, cy = at(_st["ang"][_st["pos"][i]], _st["R"], _st["centre"])
        else:
            cx, cy = _pinned
        sx, sy = cx - m["w"] / 2.0, cy - m["h"] / 2.0
        sid = m["sid"]
        # THE ADDRESS NAME IS `as`, OR POSITIONAL WHEN THERE IS NO `as`.
        #
        # It cannot be `sid`, and that is measured rather than assumed: `sid` is the
        # key into `structures`, and two ring nodes are allowed to share one --
        # nodes[{"smiles": X, "id": "A"}, {"smiles": X, "id": "A"}] renders with
        # zero findings and both reserve under `node:A`. Two instances are two
        # addresses, so a shared key cannot be an address. See `_as_id` for the two
        # namespaces and why they do not meet.
        _naid = G["nodes"][i]["as"]
        _nwhere = "nodes[%d] (s=%s)" % (i, sid)
        if m["content"] is not None:
            # A `content` slot holds a SEQUENCE. ONE reservation, ONE
            # address record, ONE blob of XML -- deliberately, because that is
            # exactly what a one-thing slot produces, so nothing downstream learns
            # that this node is different. Reserving the pieces separately would
            # leave the gaps between them unreserved (an arc could thread the
            # `plus_gap`) and would put an intra-node exemption in every report.
            _blob, _handle = [], None
            _nstruct = sum(1 for _p in m["content"] if _p["what"] == "s")
            for _p in m["content"]:
                _px = sx + _p["dx"]
                _ptop = sy + m["axis"] - _p["h"] / 2.0
                if _p["what"] == "s":
                    _pfx, _pdx = C.place_with(_p["frag"], D.get(_p["s"]),
                                              base, _px, _ptop)
                    _blob.append(_pfx)
                    if _pdx:
                        behind.append(_pdx)
                    if _nstruct == 1:
                        # THE SLOT'S ATOMS ARE THAT STRUCTURE'S ATOMS when exactly
                        # one thing is drawn -- so `{"content": [{"text": ...},
                        # {"s": X}]}` stays as atom-addressable as `{"s": X}` was.
                        # With two drawn species there is no single answer and the
                        # handle stays None; address.py then refuses an atom
                        # address rather than picking one.
                        _handle = (_p["s"], _pfx, base)
                    base += 1000
                else:
                    _blob.append(C.styled_text(
                        tid, _px + _p["w"] / 2.0, _ptop + ASCENT_EM * _p["size"],
                        _p["lines"], size=_p["size"],
                        bold=(cfg["text_bold"] if _p["what"] == "text" else False),
                        align="center"))
                    tid += 1
                if _p["lbl"]:
                    _blob.append(C.styled_text(
                        tid, _px + _p["w"] / 2.0,
                        _ptop + _p["h"] + cfg["label_gap"]
                        + ASCENT_EM * cfg["num_size"],
                        _p["llines"], size=cfg["num_size"], bold=True,
                        align="center"))
                    tid += 1
            _cas = Addressed(_naid, _inferred_kind(m))
            if _handle is not None:
                _cas = _with_handle(_cas, _handle[0], _handle[1], _handle[2],
                                    ASRC.get(_handle[0]))
            put(parts, "\n".join(_blob),
                Rect(sx, sy, m["w"], m["h"]), "node:%s" % sid, "structure",
                addr_as=_cas, addr_where=_nwhere)
        elif m["frag"] is None:
            # Text node: the type IS the species. Reserved as "structure", not
            # "label", so arcs and captions back off from it by exactly the same
            # rules they use for a drawn fragment. The <t> anchor is the first
            # line's BASELINE, hence + ASCENT_EM * size off the box top.
            put(texts, C.styled_text(tid, cx, sy + ASCENT_EM * cfg["text_size"],
                                     m["tlines"], size=cfg["text_size"],
                                     bold=cfg["text_bold"], align="center"),
                Rect(sx, sy, m["w"], m["h"]), "node:%s" % sid, "structure",
                # A WRITTEN species: the kind published figures most want as an
                # arrow endpoint.
                addr_as=Addressed(_naid, "text"), addr_where=_nwhere)
            tid += 1
        else:
            _fx, _dx = C.place_with(m["frag"], D.get(sid), base, sx, sy)
            put(parts, _fx,
                Rect(sx, sy, m["w"], m["h"]), "node:%s" % sid, "structure",
                # The FIRST of the module's two place_with sites; the feeder's
                # claim_struct is the other. Both declare, so neither can place
                # without registering.
                addr_as=_with_handle(Addressed(_naid, "structure"), sid, _fx,
                                     base, ASRC.get(sid)),
                addr_where=_nwhere)
            if _dx:
                behind.append(_dx)
            base += 1000

        # THE NODE'S BRACKETS (the node's `frame`), round the box just placed.
        #
        # Emitted here rather than by `decorate`, because a bracket is a property
        # of the SLOT, not of the molecule -- a published figure brackets a written
        # intermediate as readily as a drawn one, and a text node has no fragment
        # for a decoration to live inside. Reserved: the strokes stand
        # BRACKET_PAD outside the species and an arc has to keep off them.
        _fr = G["nodes"][i].get("frame")
        if _fr:
            _brect = Rect(sx - BRACKET_PAD, sy - BRACKET_PAD,
                          m["w"] + 2 * BRACKET_PAD, m["h"] + 2 * BRACKET_PAD)
            for _bi, _bx in enumerate(bracket_pair(tid, _brect, _fr["lip"])):
                _stroke = Rect(_brect.x if _bi == 0 else _brect.right,
                               _brect.y, 0.0, _brect.h)
                put(parts, _bx, _stroke,
                    "node:%s.bracket%d" % (sid, _bi), "structure",
                    addr_as=NotAddressed(
                        "the node's own bracket; address the node"),
                    addr_where="%s.frame" % _nwhere)
                tid += 1

        # A text node without an explicit label has nothing more to draw.
        if not m["lbl"]:
            continue

        # House style puts the label BELOW its structure, as the linear, grid and
        # radial engines all do. On a ring "below" sometimes falls inside the
        # circle and lands on an arc or the centre note, so try it, then fall back
        # to radially outward, which is always clear of the interior.
        # A PINNED node has no ring angle, so "radially outward" has no meaning for
        # it -- its label goes below, which is the house style anyway. Straight down
        # is also the direction the outward maths would give for an angle of 90.
        _lang = _st["ang"][_st["pos"][i]] if _pinned is None else 0.0
        ux, uy = ((math.cos(math.radians(_lang)),
                   math.sin(math.radians(_lang))) if _pinned is None
                  else (0.0, 1.0))
        below = (cx, cy + m["h"] / 2.0 + cfg["label_gap"] + m["lh"] / 2.0)
        d0 = (_support(m["w"], m["h"], ux, uy) + cfg["label_gap"]
              + _support(m["lw"], m["lh"], ux, uy))
        outward = [(cx + ux * (d0 + e), cy + uy * (d0 + e))
                   for e in (0.0, m["lh"], m["lh"] * 2.0)]
        mode = cfg["label_pos"] if _pinned is None else "below"
        if mode == "below":
            cands = [below]
        elif mode == "outward":
            cands = outward
        else:
            cands = [below] + outward
        lx = ly = None
        for tx, ty in cands:
            if free(Rect(tx - m["lw"] / 2.0, ty - m["lh"] / 2.0,
                         m["lw"], m["lh"]), "label"):
                lx, ly = tx, ty
                break
        if lx is None:
            lx, ly = cands[0]
        put(texts, C.styled_text(tid, lx, ly - m["lh"] / 2.0 + ASCENT_EM * cfg["num_size"],
                                 m["llines"], size=cfg["num_size"], bold=True,
                                 align="center"),
            _LBOX := Rect(lx - m["lw"] / 2.0, ly - m["lh"] / 2.0,
                          m["lw"], m["lh"]),
            "node:%s.label" % sid, "label",
            # A mark ON the node, not an object: it is the node's own name, so
            # anything wanting to point at it points at the node.
            addr_as=NotAddressed("the node's name, drawn beside it; address the "
                                 "node itself"))
        # The ATTACHMENT BOX is the structure UNION its label. The label is not
        # an object (nothing can point at it) but it is ink that belongs to the
        # node, and an edge leaving on the label's side would otherwise be trimmed
        # clear of the STRUCTURE and land inside the LABEL -- measured: structure
        # box ended at y -64.7, label occupied -58.7..-50.5, and the downward
        # attachment point came out at -52.7.
        #
        # Grown here rather than registered here, because the label's position is
        # only known once it has been searched for, while the record must be filed
        # in the same breath as the placement (the one-operation invariant).
        addr.grow(_naid, _LBOX, _nwhere)
        tid += 1

    # ---- SHARED NODES: is the caller's geometry actually consistent? ----------
    #
    # A shared-node double ring puts one id in two rings' `nodes`. It is placed
    # by the FIRST ring, so it sits on that circle by construction -- and it lies on
    # the SECOND ring's circle only if the caller chose centres and radii that make
    # it. Nothing here can fix that: those coordinates come from the caller.
    #
    # But it must not be silent. Measured on a two-ring figure: without this check
    # the shared node sits 95.0 from its own ring's centre (right) and 315.9 from
    # the other ring's, whose radius is 95 -- so that ring draws two arcs to a point
    # they come nowhere near, and `warnings` is EMPTY. A silently wrong figure is
    # the failure to avoid, so the offset is reported with the number and with
    # the constraint that would satisfy it.
    for _g in G["nodes"]:
        if len(_g["_rings"]) < 2 or _g["at"] is not None:
            continue
        _own = RSTATE[_g["_ring"]]
        _p = at(_own["ang"][_own["pos"][_g["_i"]]], _own["R"], _own["centre"])
        for _other in _g["_rings"][1:]:
            _os = RSTATE[_other]
            _d = math.hypot(_p[0] - _os["centre"][0], _p[1] - _os["centre"][1])
            if abs(_d - _os["R"]) > cfg["arc_margin"]:
                warnings.append(
                    "node %r is shared by rings %r and %r, and it does not lie on "
                    "%r's circle: it sits %.1f from that centre, which wants %.1f "
                    "(off by %+.1f). It is placed by %r, the first ring that lists "
                    "it, so %r draws its arcs to a point they do not reach. Put "
                    "%r's centre %.1f from (%.1f, %.1f), or pin that ring's "
                    "'radius'."
                    % (_g["as"], _g["_ring"], _other, _other, _d, _os["R"],
                       _d - _os["R"], _g["_ring"], _other, _other, _os["R"],
                       _p[0], _p[1]))

    # ---- ring arcs: one per EDGE, trimmed clear of both species ---------------
    #
    # `for e in EDGES` rather than `for i in range(n_arcs)`, and `j` comes off
    # `e["_to"]` rather than `(i + 1) % n`. On a single ring the two agree by
    # construction -- which is exactly why this loop can iterate edges and stay
    # byte-identical -- and an edge that is NOT a ring arc (a chord, a cross-ring
    # edge) is expressible without touching this loop.
    # ARCS FIRST, THEN LINES. A stable sort, so a figure whose edges are all arcs
    # (every old-form figure) keeps its exact order and its exact bytes.
    #
    # The reason for the order: an OFF-RING edge crosses the ring circle BY
    # CONSTRUCTION -- it starts on the circle and ends outside it -- so the
    # allocator would report an overlap on every such edge, and a feature that
    # reports an ERROR every time it is used is a feature whose reports get turned
    # off. The ring's own arcs therefore go on that edge's allow-list, which is
    # exactly what a feeder already does with the arc it feeds ("A feeder TOUCHES
    # the arc it feeds by construction"). Doing that needs the arcs' rids, hence
    # arcs first. A CHORD gets no such exemption: it runs through the interior,
    # where hitting the centre note is a real finding.
    _ring_rids: list = []
    # Arcs first, then plain lines, then BRIDGES. A bridging edge has to know where
    # every other edge runs in order to hop over it, so it is drawn once the others
    # have filed their ink. Two bridges crossing each other: the second hops the
    # first, which is the same answer a human would draw.
    _drawn_ink: dict = {}
    _edge_rids: dict = {}
    # node index -> rids of the ring arcs INCIDENT to that node. Filled by the arc
    # branch, which runs first; read by the line branch to scope its exemption.
    _node_arc_rids: dict = {}
    # ---- MARKER PAIRS become ONE curve each ----------------------------------
    #
    # If `pairs` were only validated and stored, with NOTHING reading it, a correct
    # and a deliberately wrong pairing would produce byte-identical files and the
    # declaration could not be checked by any means -- not by a rule ("one in, one
    # out" passes a pairing that sends the electron straight back where it came
    # from) and not by eye. Consuming it here is what makes it checkable: the two
    # edges of a pair are drawn as two halves of ONE CIRCLE through the marker, and
    # pairing P with a different Q gives a different circle.
    _by_e = {e["as"]: e for e in EDGES}
    _pair_jobs, _paired = [], set()
    for _g in G["nodes"]:
        for _a2, _b2 in (_g.get("_pairs") or ()):
            _ea, _eb = _by_e[_a2], _by_e[_b2]
            _ins = [x for x in (_ea, _eb) if x["_to"] == _g["_i"]]
            _outs = [x for x in (_ea, _eb) if x["_from"] == _g["_i"]]
            if len(_ins) != 1 or len(_outs) != 1:
                raise ValueError(
                    "nodes[%d].pairs: %r and %r both point the same way at marker "
                    "%r. A pair is ONE gesture -- something arrives and continues "
                    "out -- so one edge has to end at the marker and the other to "
                    "start there. (This does NOT catch a wrong pairing: pairing an "
                    "arrival with the departure back to where it came from is also "
                    "one in and one out. Only the drawn curve shows that.)"
                    % (_g["_i"], _a2, _b2, _g["as"]))
            _pair_jobs.append((_g, _ins[0], _outs[0]))
            _paired.update((_a2, _b2))
    for e in sorted(EDGES, key=lambda e: (e["via"] != "arc", e["_param"],
                                          e["crossing"] == "bridge")):
        if e["as"] in _paired:
            continue                 # drawn below, as half of its pair's curve
        i, j = e["_from"], e["_to"]
        cap0 = _step_view(e)

        if e["via"] == "line":
            # A CHORD or an OFF-RING edge: a straight arrow between the two nodes'
            # attachment points. Not a ring arc, so none of the arc machinery below
            # applies -- and nothing here is new rendering: `C.curve` already
            # carries every head variant and `graphics.straight_head` already turns
            # two points into the collinear cubic that draws a line (this is
            # wiring, not new writing).
            #
            # The endpoints come from the ADDRESS LAYER, which was built for exactly
            # this: `endpoints` attaches each end toward the other and trims by
            # ATTACH_TRIM, symmetrically, so flipping the edge moves neither end.
            # Its known limitation applies -- the reference point is the box
            # centre, so a very wide box attaches from the middle of its edge.
            _where = "edges[%s]" % e["as"]
            _pa, _pb = addr.endpoints(e["from"], e["to"], where=_where)
            if e["bulge"] is not None:
                # A BOWED EDGE LEAVES ITS NODE ALONG THE BOW, not along the chord.
                #
                # Not cosmetic. Two parallel edges between the same pair get the
                # same two attachment points from `endpoints` -- it aims each end
                # at the other's box centre, which is the same ray for both -- so
                # they start and end on the same pixel and the ledger reports a
                # real overlap at each end, on every correct use of the feature.
                # Attaching along the departure TANGENT separates them, and it is
                # also where the eye expects an arrow to leave a box.
                #
                # The tangent is read off the drawn curve rather than derived, via
                # the same canonical decode the reservation uses -- deriving it
                # here would be arc_chord's sign convention implemented a second
                # time. One pass only: moving the ends turns the tangent by a
                # second-order amount, and a fixed point would buy nothing a
                # reader could see.
                _sp = _GX.bowed_points(_pa, _pb, e["bulge"], 32)
                _far = 1000.0
                for _end, _p, _q in (("from", _sp[0], _sp[1]),
                                     ("to", _sp[-1], _sp[-2])):
                    _dx, _dy = _q[0] - _p[0], _q[1] - _p[1]
                    _n = math.hypot(_dx, _dy) or 1.0
                    _t = (_p[0] + _dx / _n * _far, _p[1] + _dy / _n * _far)
                    if _end == "from":
                        _pa = addr.attach(e["from"], _t, _where + ".from")
                    else:
                        _pb = addr.attach(e["to"], _t, _where + ".to")
            _kw = _arc_kw(cap0, cfg)
            _head = _kw.get("head", "Full")
            if cap0.get("reverse"):
                _pa, _pb = _pb, _pa
            # GREY IF THE EDGE IS ON A MUTED PATH. One value, read once, and handed
            # to the arrow, its caption and its feeders alike -- that is the whole
            # of "edges, their labels and attached decorations all go grey
            # together": the caller mutes the SET and the engine carries it to
            # every object the set owns.
            _col = cfg["muted_color"] if e["muted"] else None
            if e["bulge"] is None:
                # A NATIVE <arrow>, not a <curve>. A curve carrying an arrowhead
                # draws the head alone with no line -- measured in ChemDraw, five
                # encodings on one test page; only Head3D / Tail3D drew. See
                # graphics.straight_arrow.
                xml = _GX.straight_arrow(aid, _pa, _pb,
                                         head=(_head or "None"),
                                         head_type="Solid", color=_col)
                # The ink is generated from the SAME two endpoints the arrow was
                # built from -- an <arrow> has no CurvePoints to read back, so this
                # is the straight-line counterpart of arc_points, and one source
                # means the reservation cannot describe a different line from the
                # drawn one.
                _pts = _GX.straight_points(_pa, _pb, cfg["arc_samples"] // 2)
            else:
                # BOWED (the edge's `bulge`): how one of several parallel
                # edges gets clear of its twins. Still one native <arrow>, just an
                # elliptical one -- same object, same head, so a bowed edge and a
                # straight one are the same kind of thing to every later pass.
                xml = _GX.bowed_arrow(aid, _pa, _pb, e["bulge"],
                                      head=(_head or "None"),
                                      head_type="Solid", color=_col)
                _pts = _GX.bowed_points(_pa, _pb, e["bulge"],
                                        cfg["arc_samples"] // 2)
            aid += 1
            owner = "edge:%s" % e["as"]
            if e["crossing"] == "bridge":
                # THE HOP (`crossing: "bridge"`).
                #
                # Drawn in PIECES, because the notation needs a gap in the line and
                # nothing here can erase: segment, hump, segment, ... and only the
                # last piece carries the arrowhead. The hump is a semicircle --
                # bulge 0.5 is exactly half a circle, since sagitta = radius and
                # chord = 2 x radius.
                _hops, _tot, _crossed = bridge_hops(
                    _pts, {k: v for k, v in _drawn_ink.items() if k != e["as"]},
                    cfg["bridge_hop_r"])
                if not _hops:
                    # REFUSED, not drawn plain. `bridge` says "the crossing you
                    # see here is deliberate"; if there is no crossing, the caller
                    # is describing a figure that is not the one they wrote, and
                    # drawing a plain line would hide that forever.
                    raise ValueError(
                        "edges[%s]: 'crossing: \"bridge\"' draws the hop over "
                        "whatever this edge crosses, and it crosses nothing -- no "
                        "other edge in the figure meets it between its two ends. "
                        "Either the hop is not needed, or the edge that should "
                        "cross it is somewhere else." % e["as"])
                _r = cfg["bridge_hop_r"]
                _parts, _ink, _s0 = [], [], 0.0
                _n2 = max(2, cfg["arc_samples"] // 8)
                for _h in _hops:
                    _pA = _polyline_at(_pts, _h - _r)
                    _pB = _polyline_at(_pts, _h + _r)
                    _parts.append(("seg", _polyline_at(_pts, _s0), _pA))
                    _parts.append(("hop", _pA, _pB))
                    _s0 = _h + _r
                _parts.append(("seg", _polyline_at(_pts, _s0), _pts[-1]))
                _xmls = []
                for _k, (_kind, _u, _v) in enumerate(_parts):
                    _last = _k == len(_parts) - 1
                    if _kind == "hop":
                        _xmls.append(_GX.bowed_arrow(aid, _u, _v, 0.5, head=None,
                                                     color=_col))
                        _ink += _GX.bowed_points(_u, _v, 0.5, _n2)
                    else:
                        _xmls.append(_GX.straight_arrow(
                            aid, _u, _v,
                            head=((_head or "None") if _last else "None"),
                            head_type="Solid", color=_col))
                        _ink += _GX.straight_points(_u, _v, _n2)
                    aid += 1
                # ONE object to the ledger and to the registry, because it is one
                # EDGE -- the pieces are how the hop is drawn, not what it is.
                xml = "\n".join(_xmls)
                _pts = _ink
            # A CROSS-RING edge crosses a ring circle by construction too -- it
            # leaves one and arrives at another -- so it gets the same exemption an
            # off-ring edge gets, and for the same reason a feeder gets it against
            # the arc it feeds. A CHORD still does not: it runs through the
            # interior, where hitting the centre note or a nested ring is a real
            # finding. Measured: without this every cross-ring edge reported two to
            # four collisions, and a feature that reports an ERROR every time it is
            # used is a feature whose reports get switched off.
            # EXEMPT AGAINST ITS OWN TWO ENDS' ARCS, not against every arc on the
            # page. What is by construction is the contact AT the endpoints: the
            # edge's ink starts just outside a node's box and so do the two arcs
            # incident to that node, so they share that neighbourhood whatever the
            # figure. An arc somewhere else is a real crossing.
            #
            # The wide version (every ring arc) was measured hiding one: on a
            # Sonogashira dual cycle, B->D crossed the right ring's C->D arc and
            # reported nothing. Narrowing does NOT fix that particular case -- C
            # and D are its own endpoints, so C->D is incident to them either way.
            # It stops the unrelated arcs.
            _allow = []
            if e["_shape"] in ("offring", "crossring"):
                for _nd in (e["_from"], e["_to"]):
                    if _nd is not None:
                        _allow += _node_arc_rids.get(_nd, [])
                for _v in (e["from"], e["to"]):
                    # A PARAMETRIC end lies ON its target's ink by definition.
                    if isinstance(_v, dict):
                        _allow += _edge_rids.get(_v["of"], [])
            # `crossed`: the crossing with the NAMED edge is deliberate, so it is
            # exempt -- and only against that edge, whatever shape this one is. A
            # chord that deliberately crosses one thing still reports everything
            # else it lands on.
            for _cw in e["_crossed_set"]:
                _allow += _edge_rids.get(_cw, [])
            if e["crossing"] == "bridge":
                # EXEMPT AGAINST THE EDGES IT HOPS, and against nothing else. The
                # hump's ends sit on this edge and its apex is one radius to the
                # side; where the two lines meet at a right angle that apex lands
                # on the crossed line whatever side the bow takes, so the contact
                # is the notation rather than a fault. Naming the crossed edges
                # keeps every other collision this edge could have still audible.
                for _cn in _crossed:
                    _allow += _edge_rids.get(_cn, [])
            _res = put_edge(xml, _pts, owner, e["as"], check=True, allow=_allow)
            _edge_rids[e["as"]] = [r.rid for r in _res]
            _drawn_ink[e["as"]] = _pts
            _line_edge_caption(e, cap0, _pa, _pb)
            continue

        # THE ARC BELONGS TO ITS OWN RING. Its frame -- start angle, allotted
        # advance, radius, sweep sign -- and its CENTRE all come off
        # `RSTATE[e["ring"]]`, so a second ring somewhere else on the page draws in
        # its own frame without this loop changing.
        _st = RSTATE[e["ring"]]
        _p = _st["pos"][i]
        t0, sweep = _arc_of(_st["ang"][_p], _st["advs"][_p], meas[i]["reach"],
                            meas[j]["reach"], _st["R"], cfg, _st["s"])

        kw = _arc_kw(cap0, cfg)
        if cap0.get("reverse"):
            # draw this step against the ring: swap the heads rather than the
            # geometry, so the arc still follows the same piece of the circle
            kw["head"], kw["tail"] = kw.get("tail"), kw.get("head", "Full")
        xml = C.arc_arrow(aid, _st["centre"][0], _st["centre"][1],
                          _st["R"], _st["R"], t0, sweep, **kw)
        # A RING ARC can be on a muted path too (a figure may grey a stretch of
        # the circle, not just an off-ring branch). `recolor` is a no-op when the
        # colour is None, so an unmuted arc is the same string it always was.
        if e["muted"]:
            xml = _GX.recolor(xml, cfg["muted_color"])
        aid += 1
        pts = C.arc_points(xml, n=cfg["arc_samples"])
        owner = "step:%s->%s" % (meas[i]["sid"], meas[j]["sid"])
        # `as` on the STEP names the edge; otherwise positional, and for the same
        # reason the nodes are: `owner` is built from two structures keys, so two
        # nodes sharing a key give two steps the same owner. Step index is unique
        # by construction.
        # Off the edge table: ONE definition of what this edge is called, shared
        # with the address the table's `from` / `to` refer to.
        eid = e["as"]
        # A LINE CAN HOP OVER A RING ARC too, so an arc files its ink beside the
        # lines'. Arcs are drawn first, which is what makes that possible.
        _drawn_ink[eid] = pts
        if pts:
            # An ARC can be the one that declares a deliberate crossing too (two
            # rings' arcs can cross). Arcs are placed first, so symmetry usually
            # covers it from the other side; carried here as well so the answer
            # does not depend on which of the two was written first.
            _arc_allow = []
            for _cw in e["_crossed_set"]:
                _arc_allow += _edge_rids.get(_cw, [])
            # The other ring's arc at a `tangent_to` touch point. ARCS ONLY: this
            # set is filled in by the tangency resolver and by nothing else, and the
            # two entries it ever holds are both ring arcs.
            for _tw in e["_touch_set"]:
                _arc_allow += _edge_rids.get(_tw, [])
            ring_res = put_edge(xml, pts, owner, eid, check=True,
                                allow=_arc_allow)
            _ring_rids.extend(r.rid for r in ring_res)
            _edge_rids[eid] = [r.rid for r in ring_res]
            for _nd in (i, j):
                _node_arc_rids.setdefault(_nd, []).extend(
                    r.rid for r in ring_res)
        else:
            # not an arc that can be read back: reserve the chord so the extent is
            # still honest rather than silently missing
            a0, a1 = (at(t0, _st["R"], _st["centre"]),
                      at(t0 + sweep, _st["R"], _st["centre"]))
            ring_res = put_edge(xml, [a0, a1], owner, eid, check=False)

        # The motif goes in BEFORE the caption, so the caption is the one that has
        # to give way -- the feeder's position carries meaning (it marks where the
        # reagent enters), a caption's does not.
        motif = None
        if _has_motif(cap0):
            # "f" carries on from the species loop's own fragment base, so a feeder
            # structure cannot collide ids with a node's -- place() renumbers by it.
            ids = {"t": tid, "a": aid, "f": base}
            # `_moved` rather than a module-level list: a cycle's own report
            # channel is `warnings`, which render_cycle already returns and the
            # tool already prints. One more list would be a second place to
            # forget to read.
            _moved: list = []
            motif = feeder_motif(cap0, i, t0, sweep, _st["R"],
                                 _st["R"] - _st["reach_max"], cfg,
                                 alloc, owner, [r.rid for r in ring_res], ids,
                                 addr,
                                 warnings.append, recipe.get("structures"),
                                 anchor_src=ASRC, edge_as=e["as"],
                                 centre=_st["centre"],
                                 interior_ok=_st["interior_ok"],
                                 decorations=D, behind=behind, moved=_moved)
            for _m in _moved:
                warnings.append(
                    "%s had to give way %s to get clear (%.1f units, %s). The "
                    "figure is correct, but that piece now sits further from the "
                    "step it belongs to than the geometry chose -- usually "
                    "because two feeders share one arc. Give one an 'at' of its "
                    "own, or shorten a label."
                    % (_m[0], _m[1], _m[2], _m[3]))
            tid, aid, base = ids["t"], ids["a"], ids["f"]
            feeders.extend(motif.prims)
            # None, never 360: a motif with a box across the ring centre has no
            # angular width to quote, and a saturated number sitting in the same
            # dict as measured ones is the kind of thing that gets averaged.
            extents[owner] = (round(motif.angular_extent, 2)
                              if motif.extent_defined else None)
            if not motif.extent_defined:
                undefined.append(owner)

        cap = cap0
        lines = cap.get("above") or cap.get("below")
        if not lines:
            continue
        lines, (cw_, ch_) = _fold(lines, cfg["cap_size"], cfg)
        inward = bool(cap.get("inside"))
        # Where along the arc to try. With no feeder that is the midpoint and
        # nothing else, which is why a figure with no feeder stays byte-identical.
        # With a feeder, slide along the arc AWAY from the crossing point first,
        # and only then push radially out.
        params = [0.5]
        if motif is not None:
            away = 1.0 if float(cap0.get("at", 0.5)) <= 0.5 else -1.0
            params += [min(0.88, max(0.12, 0.5 + away * d))
                       for d in (0.14, 0.24, 0.34, 0.44)]
        # Slide ALONG the arc before flying outward: a caption that keeps its
        # normal distance from the ring and shifts sideways still reads as
        # belonging to this step, while one pushed three line-heights out does
        # not. Hence the arc sweep is the INNER loop -- with the radial push
        # outside it, the very first candidate (still sitting on the feeder) would
        # find clear air one line further out and never move sideways at all.
        px = py = None
        for extra in (0.0, ch_, ch_ * 2.0, ch_ * 3.0):
            for prm in params:
                tm = t0 + sweep * prm
                ux, uy = math.cos(math.radians(tm)), math.sin(math.radians(tm))
                if inward:
                    ux, uy = -ux, -uy
                d0 = cfg["cap_gap"] + _support(cw_, ch_, ux, uy)
                tx, ty = (_st["centre"][0] + ux * (_st["R"] + d0 + extra),
                          _st["centre"][1] + uy * (_st["R"] + d0 + extra))
                if free(Rect(tx - cw_ / 2.0, ty - ch_ / 2.0, cw_, ch_)):
                    px, py = tx, ty
                    break
            if px is not None:
                break
        if px is None:
            # Nowhere on this arc is clear. Place it at the midpoint anyway so the
            # figure is still complete, but SAY SO -- a silent overlap is the one
            # outcome that must not happen.
            tm = t0 + sweep / 2.0
            ux, uy = math.cos(math.radians(tm)), math.sin(math.radians(tm))
            if inward:
                ux, uy = -ux, -uy
            d0 = cfg["cap_gap"] + _support(cw_, ch_, ux, uy)
            px, py = (_st["centre"][0] + ux * (_st["R"] + d0),
                      _st["centre"][1] + uy * (_st["R"] + d0))
            if motif is not None:
                warnings.append(
                    "steps[%d]: the arc is too short to fit both the feeder and '%s'; the "
                    "caption overlaps. Give the feeder a different 'at', or shorten the "
                    "caption."
                    % (i, " / ".join(lines)[:60]))
        put(texts, C.styled_text(tid, px, py - ch_ / 2.0 + ASCENT_EM * cfg["cap_size"],
                                 lines, size=cfg["cap_size"], align="center"),
            Rect(px - cw_ / 2.0, py - ch_ / 2.0, cw_, ch_),
            owner + ".caption", "caption",
            # A mark ON the step. The step is the object; this is its wording.
            addr_as=NotAddressed("the step's own caption; address the step"))
        tid += 1

    # ---- `crossed` must actually cross ---------------------------------------
    #
    # Same rule as `bridge`: the declaration asserts "the crossing you see here is
    # deliberate". If there is no crossing, the caller is describing a figure that
    # is not the one they wrote, and silently accepting it would hide that for good.
    # Checked once every edge has filed its ink, so an arc counts as a counterpart
    # too -- deliberately crossing a ring arc is exactly a case this has to
    # cover.
    for e in EDGES:
        for _cw in e["_crossed_set"]:
            if e["as"] > _cw:
                continue                 # one message per pair, not two
            _ia, _ib = _drawn_ink.get(e["as"]), _drawn_ink.get(_cw)
            if not _ia or not _ib:
                continue                 # neither drew; a later pass will say why
            _hits = sum(1 for _pa, _pb in zip(_ia, _ia[1:])
                        for _qa, _qb in zip(_ib, _ib[1:])
                        if _seg_cross(_pa, _pb, _qa, _qb))
            if not _hits:
                raise ValueError(
                    "edges[%s]: 'crossing: {\"crossed\": %r}' says the crossing "
                    "between these two is deliberate, and they do not cross -- "
                    "their drawn inks meet nowhere between their ends. Either the "
                    "declaration is not needed, or the edge that should cross this "
                    "one is somewhere else." % (e["as"], _cw))

    # ---- a marker's PAIRED edges, as two halves of one circle ----------------
    #
    # Three points define exactly one circle, so "the curve that leaves P, passes
    # through the marker and arrives at Q" is not a shape anyone chose. Each half is
    # filed as its OWN edge -- the one-operation invariant still holds, both ends
    # stay addressable -- and what makes the pair read as one gesture is that both
    # halves lie on the same circle. Only the second half carries the arrowhead: one
    # electron, one destination, one head.
    for _mk, _ein, _eout in _pair_jobs:
        _mrec = addr.resolve(_mk["as"], "nodes[%d]" % _mk["_i"]).record
        _M = (_mrec.rect.cx, _mrec.rect.cy)
        _P = addr.attach(_ein["from"], _M, "edges[%s].from" % _ein["as"])
        _Q = addr.attach(_eout["to"], _M, "edges[%s].to" % _eout["as"])
        _cap_i, _cap_o = _step_view(_ein), _step_view(_eout)
        _hd = _arc_kw(_cap_o, cfg).get("head", "Full") or "None"
        _sol = arc_through(_P, _M, _Q)
        _halves = []
        if _sol is None:
            # COLLINEAR: P, the marker and Q are on one straight line, so there is
            # no finite circle. Two straight arrows is the same gesture -- a
            # straight line IS the curve through three collinear points -- and it
            # is still one gesture rather than a fallback to something else.
            _halves = [(_GX.straight_arrow(aid, _P, _M, head="None",
                                           color=(cfg["muted_color"]
                                                  if _ein["muted"] else None)),
                        _GX.straight_points(_P, _M, cfg["arc_samples"] // 2),
                        _ein, _P, _M),
                       (_GX.straight_arrow(aid + 1, _M, _Q, head=_hd,
                                           head_type="Solid",
                                           color=(cfg["muted_color"]
                                                  if _eout["muted"] else None)),
                        _GX.straight_points(_M, _Q, cfg["arc_samples"] // 2),
                        _eout, _M, _Q)]
            aid += 2
        else:
            _cx3, _cy3, _r3, _t3, _sw3 = _sol
            # SPLIT AT THE MARKER'S CENTRE, with no gap.
            #
            # A gap the width of the marker's box would give the label clear
            # space, but measured on a two-pair figure: BOTH gestures are
            # near-horizontal where they pass through the marker -- they have to
            # be, since both pass through the same point -- so their two ends on
            # the right of the box come out 0.4 units apart, and the gap hides the
            # very crossing it exists to show. So the halves meet AT the marker;
            # the label is drawn after the arcs and sits on top, which is how
            # published figures draw it.
            _tM = math.degrees(math.atan2(_M[1] - _cy3, _M[0] - _cx3))
            _dm = (_tM - _t3) % 360.0
            if _sw3 < 0:
                _dm -= 360.0
            _spans = [(_t3, _dm, _ein, "None"),
                      (_t3 + _dm, _sw3 - _dm, _eout, _hd)]
            for _ts, _sp, _ee, _hh in _spans:
                if abs(_sp) < 1e-6:
                    # The marker's box swallowed this half whole. Nothing to draw,
                    # and a zero-length arrow is exactly the degenerate object to
                    # avoid.
                    warnings.append(
                        "edges[%s]: the marker %r sits so close to this end that "
                        "its half of the paired curve has no length; it is not "
                        "drawn. Move the marker further from %r."
                        % (_ee["as"], _mk["as"],
                           _ee["from"] if _ee is _ein else _ee["to"]))
                    continue
                _x = C.arc_arrow(aid, _cx3, _cy3, _r3, _r3, _ts, _sp,
                                 head=(None if _hh == "None" else _hh),
                                 head_type="Solid")
                if _ee["muted"]:
                    _x = _GX.recolor(_x, cfg["muted_color"])
                _pp = C.arc_points(_x, n=cfg["arc_samples"])
                _halves.append((_x, _pp, _ee,
                                _pp[0] if _pp else _P, _pp[-1] if _pp else _Q))
                aid += 1
        # EVERY HALF AT THIS MARKER MEETS EVERY OTHER THERE, by construction: all
        # four terminate at the same point, which is what a marker IS. Exempt
        # against the marker's own halves and nothing else, so a gesture that runs
        # along something else is still reported.
        # ...and against the MARKER'S OWN box. The curve passes through the
        # marker's centre -- that is the definition of the gesture -- so it runs
        # under the label, and the label is drawn after the arcs and covers it.
        # Looked up by owner rather than remembered at placement time: one source
        # for "which reservation is this node's", and it cannot go stale.
        _mk_owner = "node:%s" % meas[_mk["_i"]]["sid"]
        _mine_rids = _mk.setdefault(
            "_half_rids",
            [r.rid for r in alloc.reservations() if r.owner == _mk_owner])
        for _x, _pp, _ee, _u, _v in _halves:
            _own = "edge:%s" % _ee["as"]
            _rs = put_edge(_x, _pp, _own, _ee["as"], check=True,
                           allow=list(_ring_rids) + list(_mine_rids))
            _mine_rids.extend(r.rid for r in _rs)
            _edge_rids[_ee["as"]] = [r.rid for r in _rs]
            _drawn_ink[_ee["as"]] = _pp
            _line_edge_caption(_ee, _step_view(_ee), _u, _v)

    # ---- each path's name ----------------------------------------------------
    #
    # A GROUP-LEVEL label: a "TM" over a pair of cross-ring edges, an "ROH" over a
    # pair of edges from one source, "path A" / "path B". Drawn ONCE for the set,
    # which is the whole reason the path row exists -- written on a member
    # instead, a second member could write it again and nothing could catch the
    # duplicate.
    #
    # The members are reached through the ADDRESS LAYER's group form, `{"of":
    # [...]}` -- the form the address layer has for exactly this and refuses
    # as an endpoint. So the ink the label is positioned against is the ink the
    # ledger accounted, not a second copy computed here.
    for _gp in PATHS:
        _nm = _gp.get("name")
        if not _nm:
            continue
        _recs = addr.resolve({"of": _gp["edges"]},
                             "paths[%d].name" % _gp["_i"]).records
        _ink = [pt for r in _recs for pt in r.ink]
        if not _ink:
            # Every member is an object with no ink of its own. Nothing to sit
            # beside, and guessing a spot would put the name somewhere the caller
            # cannot predict.
            raise ValueError(
                "paths[%d]: %r has a name but none of its edges drew any ink, so "
                "there is nothing for the name to sit beside."
                % (_gp["_i"], _gp["as"]))
        # The set's overall travel: first member's first point to last member's
        # last, in the order the caller listed them. `above` / `below` are
        # perpendicular to THAT, the same way they are on a single line edge.
        _p0, _p1 = _recs[0].ink[0], _recs[-1].ink[-1]
        _cx = sum(q[0] for q in _ink) / len(_ink)
        _cy = sum(q[1] for q in _ink) / len(_ink)
        _dx, _dy = _p1[0] - _p0[0], _p1[1] - _p0[1]
        _dn = math.hypot(_dx, _dy)
        if _dn < 1e-9:
            # A closed group -- it ends where it starts, so it has no travel
            # direction. Fall back to the page's up/down, which is what "above"
            # means when nothing else defines it.
            _ux, _uy = 0.0, -1.0
        else:
            _ux, _uy = -_dy / _dn, _dx / _dn
        for _key, _side in (("above", -1.0), ("below", +1.0)):
            _ls = _nm.get(_key)
            if not _ls:
                continue
            _ls, (_pw, _ph) = _fold(_ls, cfg["cap_size"], cfg)
            _px = _py = None
            for _extra in (0.0, _ph, _ph * 2.0, _ph * 3.0, _ph * 4.0):
                _g = cfg["cap_gap"] + _support(_pw, _ph, _ux, _uy) + _extra
                _tx, _ty = _cx + _ux * _side * _g, _cy + _uy * _side * _g
                if free(Rect(_tx - _pw / 2.0, _ty - _ph / 2.0, _pw, _ph)):
                    _px, _py = _tx, _ty
                    break
            if _px is None:
                _g = cfg["cap_gap"] + _support(_pw, _ph, _ux, _uy)
                _px, _py = _cx + _ux * _side * _g, _cy + _uy * _side * _g
                warnings.append(
                    "paths[%d]: nothing beside %r is clear for its '%s' name; it "
                    "is drawn at the conventional spot and overlaps."
                    % (_gp["_i"], _gp["as"], _key))
            put(texts,
                C.styled_text(tid, _px, _py - _ph / 2.0 + ASCENT_EM * cfg["cap_size"],
                              _ls, size=cfg["cap_size"], align="center",
                              color=(cfg["muted_color"] if _gp["muted"] else None)),
                Rect(_px - _pw / 2.0, _py - _ph / 2.0, _pw, _ph),
                "path:%s.name" % _gp["as"], "caption",
                # THE PATH is the object; this is its wording. Same rule as a
                # step's caption -- address the path, not its name.
                addr_as=NotAddressed("the path's own name; address the path"),
                addr_where="paths[%d].name" % _gp["_i"])
            tid += 1

    # ---- each ring's label ---------------------------------------------------
    #
    # ONE PER RING, at that ring's own centre -- and the label's `position` decides
    # whether it goes in the middle at all. `above` / `below` exist because
    # published figures need them: when two circles are tangent the centre is
    # squeezed out, and the label goes to the top. With more than one ring on a
    # page the contention is real.
    for _rgs in RSTATE.values():
        _rg = _rgs["rg"]
        _lab2 = _rg.get("label")
        if not _lab2 or not _lab2.get("text"):
            continue
        _pos2 = _lab2.get("position", "center")
        _cxr, _cyr = _rgs["centre"]
        if _pos2 == "center":
            _l2, _w2, _h2 = _rgs["note_lines"], _rgs["nw"], _rgs["nh"]
            _lx2, _ly2 = _cxr, _cyr
        else:
            _t2 = _lab2["text"]
            _l2, (_w2, _h2) = _fold([_t2] if isinstance(_t2, str) else list(_t2),
                                    cfg["center_size"], cfg)
            # OUTSIDE the circle, clear of the species annulus, on the side asked
            # for. `cap_gap` is the same clearance a step caption keeps.
            _off = (_rgs["R"] + _rgs["reach_max"] + cfg["cap_gap"] + _h2 / 2.0)
            _lx2 = _cxr
            _ly2 = _cyr - _off if _pos2 == "above" else _cyr + _off
        # The address id is the RING's own name for a centre label -- a figure has
        # at most one per ring -- and `center` stays the name of the single-ring
        # one, so no old figure's address set moves.
        _aid2 = "center" if _rg["as"] == RING["as"] else ("%s.label" % _rg["as"])
        put(texts, C.styled_text(tid, _lx2,
                                 _ly2 - _h2 / 2.0 + ASCENT_EM * cfg["center_size"],
                                 _l2, size=cfg["center_size"], align="center"),
            Rect(_lx2 - _w2 / 2.0, _ly2 - _h2 / 2.0, _w2, _h2),
            ("center.note" if _rg["as"] == RING["as"]
             else "ring:%s.label" % _rg["as"]), "caption",
            # A ring's row covers "ring label, centre note" -- so a ring's label
            # IS an object.
            addr_as=Addressed(_aid2, "ring"), addr_where=_aid2)
        tid += 1

    if report is not None:
        # ONE channel for "something is wrong with this figure". With `warnings`
        # and the allocator's findings as two places to look, the one every caller
        # actually reads can come back [] while the other holds an ERROR --
        # measured on a figure whose two feeder labels were declared on top of
        # each other and whose "warnings" was empty. An instrument with two dials
        # that disagree is worse than one dial.
        #
        # ERROR only: NOTE and WARN are prediction drift, which this module absorbs
        # by design, and promoting them would train the reader to skip the list.
        # The whole allocator is read, not just this engine's own reservations -- if
        # a caller hands in a page with something already on it, an overlap with
        # that is still this figure's problem.
        errs = sorted((f for f in alloc.findings if f.severity >= Severity.ERROR),
                      key=lambda f: -f.amount)
        for f in errs[:3]:
            warnings.append("allocator/%s: %s" % (f.code, f.message))
        if len(errs) > 3:
            warnings.append(
                "allocator: %d further finding(s) -- read layout_report for all %d"
                % (len(errs) - 3, len(errs)))
        report["warnings"] = warnings
        report["feeder_extents"] = extents
        # Which motifs came back None above, and why -- so "extent is null" is
        # answerable without reading this module. These are governed by the radial
        # check against the centre note (_disk_clash), not by the angle budget.
        report["feeder_extent_undefined"] = undefined
        # What the packer decided, so a caller can see WHY the ring came out this
        # size -- and so a test can assert that a ring with no feeder got the
        # legacy radius back rather than a number that merely looks like it.
        # `ring` stays RING 0, so every existing caller and test reads what it
        # always did; `rings` adds one entry per ring.
        report["ring"] = dict(trace, radius=R, legacy_radius=lo,
                              gaps=[round(g, 2) for g in gaps])
        report["rings"] = {
            k: dict(v["trace"], radius=v["R"], legacy_radius=v["lo"],
                    gaps=[round(g, 2) for g in v["gaps"]],
                    centre=list(v["centre"]), nodes=list(v["rg"]["nodes"]),
                    # HOW FAR THE SPECIES STICK OUT past the circle. Without it
                    # `centre` and `radius` do not add up to anything a caller can
                    # check: two rings are placed clear of each other by their
                    # OUTER extent, which is radius + this, and a spacing that
                    # forgot it would read as perfectly reasonable numbers.
                    reach=round(v["reach_max"], 2),
                    # AS RESOLVED, not as asked: `tangent_to` computes the centre
                    # and may rotate the ring, and a caller that cannot read back
                    # where its species ended up has to measure the ink to find out.
                    start_angle=round(v["rg"]["start_angle"], 2),
                    interior_allocatable=v["interior_ok"])
            for k, v in RSTATE.items()}
        # The path table, as the engine resolved it -- which edges each
        # path holds, whether it is muted, and what it is called. Published for the
        # same reason `rings` is: a caller cannot otherwise see that the group it
        # wrote is the group the engine drew.
        report["paths"] = {
            gp["as"]: {"edges": list(gp["edges"]), "muted": gp["muted"],
                       "named": bool(gp["name"])}
            for gp in PATHS}
        # The address registry, published whether or not the caller passed one in
        # -- see the signature. LAYOUT frame, so these rects are the engine's own
        # and may be negative; nothing here has been shifted to PAGE.
        report["addresses"] = addr

    bb = alloc.bbox()
    # feeders last: appending an empty list is a no-op, which is what keeps a
    # figure with no feeder byte-identical to the pre-feeder output.
    # `behind` is the same trick at the front, and for the same reason.
    body = "\n".join(behind + arcs + parts + texts + feeders)
    if behind:
        body, _z = C.stamp_z(body)     # only when there is something to layer
    return C.finalize(header_prefix, body,
                      [bb.x, bb.y, bb.right, bb.bottom],
                      page_mode=page_mode, margin=cfg["margin"], shift=shift)
