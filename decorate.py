"""decorate.py — a primitive, drawn at an anchor.

A survey of published figures listed twenty decoration features. Almost all of
them are one sentence with the nouns swapped: A PRIMITIVE, DRAWN AT AN ANCHOR.
Glow, dot, position number, region colour -- not four features, two ideas and a
short list.

WHAT COUNTS AS AN ANCHOR

A substructure query, resolved by smiles_to_frag.match_node_ids. Not coordinates
(dead the moment a substituent moves) and not atom indices (dead the moment the
SMILES is rewritten). That distinction is the whole point: it is what lets one
spec run over sixteen scope cells.

WHAT COUNTS AS A PRIMITIVE

Two, not three:

    a DISC   -- cdxml_core.disc(). A dot on an atom and a glow behind a group are
                the same object at two radii: both, drawn by hand in ChemDraw,
                come back as orbitals.
    TEXT     -- styled_text, which the pipeline has had all along.

plus region COLOUR, which is not a primitive at all: beautify.color_region
takes atom ids, and the anchor is only what hands it some.

and then a third, later:

    a CURVE  -- cdxml_core.curve(). The squiggle, and so far only the squiggle.

DELIBERATELY NOT GENERAL. Five styles, TWO resolvers, no plug-in table. This is
the part with the highest over-design risk, and the defence is: cap the list, and
only when a new style arrives discuss whether there is an abstraction here.

For the fifth style the answer was AGAINST abstracting. What squiggle would have
had to share is a resolver that returns atoms, and it needs one that returns a
bond; and a home next to the structure, when it needs one inside it. The two
things a plug-in table would have factored out are exactly the two it does
differently, so the table would have had two entries and one of them would have
been an exception. Ask again at the sixth.

THE STRUCTURE IS ALWAYS THE FRONT LAYER

Two ways to get this wrong. On the ordering: a dot marking an atom drawn in
FRONT of it is a coloured disc sitting on top of the thing it is supposed to
point at. On the mechanism: "emit it first so it is behind" assumes document
order decides, and it does not -- ChemDraw sorts by Z (see cdxml_core.stamp_z).

So there is no front pile. apply() returns decorations back-to-front -- glows,
then dots and numbers -- and the caller puts the structure after them and stamps
Z over the lot, which is what puts the molecule on top of its own annotations.
"""

from __future__ import annotations

import itertools
import math
import re

import beautify as B
import cdxml_core as C
import graphics as GX
import smiles_to_frag as S

# How far a glow reaches beyond the group it covers, in bond lengths.
#
# LOOSELY held on purpose, and the reason is a finding rather than a shrug:
# 0.25, 0.5 and 0.8 shown side by side on the same molecule cannot be told
# apart, so the largest was picked. A shaded orbital fades out at its edge, so the
# radius written here is nowhere near where the colour visibly stops -- several
# units of it disappear into the gradient. Anything in [0.25, 0.8] looks the
# same; 0.8 is the end of the range that was actually looked at, and a bigger
# glow is the safer error when the alternative is one that crops the structure.
#
# So: do NOT tighten this to a decimal on the strength of one figure. It is not
# that kind of constant, unlike FLOAT_DEPTH or LOCANT_OFFSET, which were measured
# off a drawing and do change the picture.
HALO_PAD = 0.8

# The squiggle, in fractions of C.BOND_LENGTH. PROVISIONAL -- these are a
# starting guess and have not yet been measured against a drawn figure,
# unlike LOCANT_OFFSET or FLOAT_DEPTH. Tune them against the picture, not against
# each other.
#
# NOT the "wavy" bond style below, and the two are one letter apart in intent:
# this draws a separate wave ACROSS an intact bond (the retrosynthetic cut mark),
# while "wavy" makes the BOND ITSELF a wave (the scope-table attachment point).
# Different objects, different CDXML, different meaning to a chemist.
SQUIGGLE_AMP = 0.25            # how far it swings, along the bond
SQUIGGLE_WAVELENGTH = 0.20     # one full back-and-forth, along the spine
SQUIGGLE_OVERHANG = 0.30       # how far past each end of the bond it reaches

# --- per-bond appearance ------------------------------------------------------
# One CDXML attribute per style, set on the bond an anchor names. All four are
# properties of the bond ELEMENT, so they live inside the <fragment> and reach
# every figure tool for free -- the same reason "fill" does and "halo" does not.
#
# "wavy" IS NOT "squiggle", and the names are close enough to need a written
# warning:
#
#     wavy      the BOND ITSELF is drawn as a wave. The attachment-point mark a
#               scope table puts where the shared fragment was cut away -- about
#               44 of them in one published table.
#     squiggle  a wave drawn ACROSS an intact bond, at right angles to it: the
#               retrosynthetic "cut here" mark. A separate <curve>, not a bond
#               property, which is why it is handled further down.
#
# Evidence, since it differs per value:
#   Bold / Hash / Dash / Wavy   Display, DTD-enumerated AND used in ChemDraw's own
#                               documents (83 / 36 / 8 / 6 occurrences).
#   dative                      Order, and the ONLY value here with no DTD entry
#                               and no file precedent anywhere -- verified instead
#                               by round trip: written, opened, saved by
#                               ChemDraw 26, still there. The
#                               three values on the same spec page -- ionic,
#                               hydrogen, threecenter -- were dropped by that same
#                               save, so they are NOT offered here.
#
# Hash vs Dash is deliberately left to the caller. Both are flat, non-wedge bonds
# (the tapered stereo hash is WedgedHashBegin/End, which the wedge writer already
# emits); ChemDraw's own samples use them in different documents with no
# discernible convention, so claiming one means "stereochemistry" and the other
# "secondary interaction" would be inventing a rule from two teaching files.
BOND_STYLES = {
    "bold":   {"Display": "Bold"},
    "hash":   {"Display": "Hash"},
    "dash":   {"Display": "Dash"},
    "wavy":   {"Display": "Wavy"},
    "dative": {"Order": "dative"},
}

# dative is the only one whose direction is chemistry rather than drawing.
_ORIENTED = ("dative",)

# --- one annotation primitive, five variants ---------------------------------
# A cubic Bezier between two points, with an arrowhead and a caption that are
# each optional. Five things a figure asks for, and one emitter; which variant
# you get is which keys you wrote:
#
#   on          text   offset   arrow   -> what it draws
#   1 atom      "α"    --       --         the letter beside the atom
#   1 atom      "95%D" [dx,dy]  full       leader line out to the caption
#   a bond      "..."  [dx,dy]  --         a break-this-bond leader
#   2 atoms     --     --       full       an electron-pushing arrow
#   2 atoms     --     --       half       ...with ONE barb: the fishhook a
#                                          homolysis draws. Not a new mechanism:
#                                          ArrowheadHead="HalfLeft", which
#                                          ChemDraw's own documents use 30 times.
#   2 atoms     --     --       --         a bare arc, e.g. the 25 grey lines
#                                          one figure draws to mark where a minor
#                                          regioisomer connects
#
# `offset` DECIDES THE LEADER, by being present -- not by being large enough.
# A magnitude test would need a threshold, a threshold needs something to be
# measured against, and there is nothing here to measure against. Presence is a
# question the caller has already answered by typing the key or not.
#
# NOT ACROSS STRUCTURES. Anchors resolve through one fragment's index_map, so
# both ends of an arrow live in one cell. The piperidine-attacks-pyridinium arrow
# spans two, and that needs to know where a packer put each cell -- the spatial
# ledger the engines do not share. Left as a limit rather than worked around: a
# caller cannot supply those coordinates either, because the engine chose them.
ANNOTATE_ARROWS = {"full": "Full", "half": "HalfLeft", "none": None}
ANNOTATE_BULGE = 0.35        # sideways swell, as a fraction of the chord
# MINIMUM SAG, in fractions of C.BOND_LENGTH -- a FLOOR, not a replacement.
# `bulge` is a fraction of the CHORD, which is right until the chord is short:
# an arrow from a bond's midpoint onto one of that bond's own atoms -- a pi bond
# collapsing onto its oxygen, a sigma bond onto the atom that leaves, i.e. half
# of every mechanism arrow drawn -- has a chord of HALF A BOND. At 0.35 that is
# a 2.5-unit sag, which renders as an unreadable hook, while the same 0.35 on a
# full-bond chord is fine. So the two ends' distance cannot be the only
# reference; the drawing's own unit has to put a floor under it.
# WHAT THIS MOVES, stated rather than wished away: every annotate arrow whose
# sag was under the floor changes shape. At the default bulge 0.35 and a floor of
# 0.60 that is every chord under 0.60*14.4/0.35 = 24.7 units, i.e. under 1.71
# bond lengths. Long arrows are untouched. The test suite passes either way
# because nothing in it asserts on an annotate arc's geometry -- green here means
# uncovered, not unchanged.
# WHERE 0.60 COMES FROM. Not a measurement -- there was no hand-drawn
# electron-pushing arrow at hand to read one off. So it was CHOSEN, off a
# six-rung ladder rendered at 0.30 / 0.45 / 0.60 / 0.75 / 0.90 / 1.20 on the
# same five arrows: 0.60 is the rung where a bond-onto-its-own-atom arrow first
# reads clearly without any arc swelling over an atom label. Re-render such a
# ladder before changing it.
ANNOTATE_SAG_FLOOR = 0.60

# ARROW ENDS STOP SHORT OF WHAT THEY POINT AT. An imported mechanism
# (mech-USPTO-31K) names the atom or the bond and says nothing about clearance,
# so the clearance has to be added here. In fractions of C.BOND_LENGTH:
#
#   "label"   a drawn atom label -- any heteroatom, an H, a charged carbon.
#             MEASURED, not chosen: _label_overhang puts an atom label's half
#             extent at 0.34-0.46 bond (O 0.38, N 0.36, Br 0.46, S 0.34), so
#             anything under ~0.46 leaves the arrowhead inside the glyph. 0.55
#             is the widest of those plus a small gap.
#   "bond"    the arrow starts on or lands on a particular bond. No glyph to
#             clear -- only enough to read as "beside the bond, not on it".
#   "carbon"  a plain carbon vertex, which is a corner and nothing else.
#   "centre"  a ring centre ("attacks the ring"): the head BELONGS on the
#             centroid, so nothing is pulled back.
#   "ion"     a CHARGED atom, and it gets LESS room than a neutral label, not
#             more. The rule is chemistry rather than typography: when the
#             nucleophile is an anion it is the CHARGE that reacts, so the
#             arrow has to read as leaving the minus sign. At the neutral label's
#             0.55 bond the tail floats away from it and looks like it starts
#             from nowhere. PROVISIONAL, to be tuned off a rendered ladder.
ANNOTATE_INSET = {"label": 0.55, "ion": 0.35, "bond": 0.25, "carbon": 0.20,
                  "centre": 0.0}
# ... and a floor on what survives. Two ends that each pull back can meet in the
# middle: the commonest mechanism arrow of the lot -- a bond onto one of that
# bond's own atoms -- has a chord of HALF a bond, which 0.25 + 0.55 would
# overrun. Both setbacks then scale down together, keeping at least this
# fraction of the chord. GEOMETRY, NOT TUNING: when the two targets are 0.5 bond
# apart and the label alone is 0.4 bond wide, no pair of on-chord endpoints can
# be outside both labels AND between the targets. On those the bow is what
# separates the ends; the clearance is whatever is left.
ANNOTATE_INSET_KEEP = 0.35


def endpoint_kind(element: int = 6, charge: int = 0, natoms: int = 1,
                  node_type: str | None = None) -> str:
    """Which ANNOTATE_INSET bucket one end of an arrow falls in.

    `natoms` is how many atoms the anchor named: one is an atom, two is the bond
    between them, three or more is a centroid (a ring). A CHARGED carbon counts
    as "label" and not as "carbon" -- ChemDraw draws it with its symbol and its
    hydrogen count, so it takes label-sized room even though it is a carbon.
    """
    if natoms >= 3:
        return "centre"
    if natoms == 2:
        return "bond"
    # A NodeType means ChemDraw draws TEXT there and not a vertex: a contracted
    # group (NodeType="Fragment"), an R placeholder ("GenericNickname"), an
    # attachment point. MEASURED, not assumed -- with abbreviate=True a Cbz
    # collapses to one NodeType="Fragment" node carrying a <t>, and that node has
    # no Element attribute at all, so classifying on element alone would file
    # the widest label in the figure as a bare carbon and give it the NARROWEST
    # clearance: measured, 2.88 units where the label bucket wants 7.92.
    if node_type:
        return "label"
    if charge:
        return "ion"
    if element != 6:
        return "label"
    return "carbon"


def endpoint_backs(kind_a: str, kind_b: str):
    """(setback at a, setback at b) in CDXML units, for C.arc_chord_trimmed."""
    return (ANNOTATE_INSET.get(kind_a, 0.0) * C.BOND_LENGTH,
            ANNOTATE_INSET.get(kind_b, 0.0) * C.BOND_LENGTH)


ANNOTATE_TEXT_SIZE = 9
ANNOTATE_TEXT_PAD = 4.0      # a hugging caption's clearance from its atom

# Which structures got a DELOCALIZED depiction this batch. Per-batch, cleared by
# the tool layer, and read by _chem_report -- the same shape as
# smiles_to_frag.TEMPLATE_MISSES. It exists so the round-trip mismatch such a
# figure PROVOKES can be explained instead of suppressed: the double bonds really
# are gone, the file really is not the molecule, and a checker that quietly
# stopped saying so would be the worse bug (see the comment on the inchi-mismatch
# warning in chemwand_mcp).
DELOCALIZED: list = []


STYLES = (("color", "fill", "halo", "dot", "locants", "squiggle", "annotate",
           "delocalized", "lonepair")
          + tuple(BOND_STYLES))


def _encloses_a_ring(bonds) -> bool:
    """True when these (id, endA, endB) bonds contain a cycle. Union-find.

    Not `len(bonds) >= len(atoms)`, which is the tempting one-liner and is wrong
    on a disconnected match: a benzene plus a loose ethyl is 5 atoms and 4 bonds,
    fails that test, and does contain a ring.
    """
    parent: dict = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for _bid, a, b in bonds:
        ra, rb = find(a), find(b)
        if ra == rb:
            return True                  # this bond closes a loop
        parent[ra] = rb
    return False


def apply_setback(a, b, back):
    """Move an arrow's two ends back along the chord, toward each other.

    ENDS FIRST, ARC SECOND -- the right causality. The two
    endpoints are what the figure is about; the curve is derived from them. So
    the clearance moves the ENDPOINTS and the arc is built through wherever they
    landed. Trimming the finished arc instead lets the
    curve decide where its own ends are, which is that causality reversed.

    Both setbacks scale down together when the chord cannot afford them; see
    ANNOTATE_INSET_KEEP for why that is the common case and not the corner one.
    """
    back0, back1 = back
    dx, dy = b[0] - a[0], b[1] - a[1]
    L = math.hypot(dx, dy)
    if not L or not (back0 or back1):
        return a, b
    room = L * (1.0 - ANNOTATE_INSET_KEEP)
    if back0 + back1 > room:
        k = room / (back0 + back1)
        back0, back1 = back0 * k, back1 * k
    ux, uy = dx / L, dy / L
    return ((a[0] + ux * back0, a[1] + uy * back0),
            (b[0] - ux * back1, b[1] - uy * back1))


# AUTO-ROUTING an annotate arc. The arc has exactly TWO degrees of freedom --
# which side it bows to and how deep -- so routing is a search over those two and
# nothing else. No new geometry: every candidate is an arc the emitter already
# draws. The bar is "avoid overlap where possible; where it cannot be avoided,
# stay readable" -- so this minimises overlap, it does not promise to abolish it.
ANNOTATE_ROUTE_DEPTHS = (0.30, 0.45, 0.60, 0.85, 1.20)
ANNOTATE_ROUTE_SAMPLES = 24
# Crossing a BOND is what makes an arrow unreadable, so it dominates. Grazing an
# atom's drawn form is worse than clear air and better than a crossing, hence a
# separate, smaller weight rather than one number for both.
ANNOTATE_ROUTE_CROSS = 100.0
ANNOTATE_ROUTE_GRAZE = 8.0
# RUNNING ALONGSIDE a bond is ugly even when it never crosses one, and crossings
# alone cannot tell apart the two candidates for one measured arenium arrow:
# both crossed exactly once, the tie went to the shallower, and the shallower was
# the one threading between the N and H labels while the other sat in the ring's
# empty middle. So proximity to a bond is scored too, over this clearance.
ANNOTATE_ROUTE_NEAR = 0.25          # fractions of C.BOND_LENGTH
ANNOTATE_ROUTE_NEAR_COST = 6.0
ANNOTATE_ROUTE_DEPTH_COST = 1.0     # tie-break toward the shallower arc


def _point_seg_dist(p, r, s) -> float:
    """Shortest distance from point p to segment r-s."""
    dx, dy = s[0] - r[0], s[1] - r[1]
    L2 = dx * dx + dy * dy
    if L2 == 0.0:
        return math.hypot(p[0] - r[0], p[1] - r[1])
    t = ((p[0] - r[0]) * dx + (p[1] - r[1]) * dy) / L2
    t = 0.0 if t < 0.0 else (1.0 if t > 1.0 else t)
    return math.hypot(p[0] - (r[0] + t * dx), p[1] - (r[1] + t * dy))


def _seg_cross(p, q, r, s) -> bool:
    """Do segments p-q and r-s properly cross?"""
    def side(o, u, v):
        return (u[0] - o[0]) * (v[1] - o[1]) - (u[1] - o[1]) * (v[0] - o[0])
    d1, d2 = side(p, q, r), side(p, q, s)
    d3, d4 = side(r, s, p), side(r, s, q)
    return ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0))


def arc_polyline(a, b, bulge, n: int = ANNOTATE_ROUTE_SAMPLES) -> list:
    """The arc arc_chord would draw, sampled into n+1 points."""
    cx, cy, r, t0, sweep = C.arc_circle(a[0], a[1], b[0], b[1], bulge)
    out = []
    for k in range(n + 1):
        t = math.radians(t0 + sweep * float(k) / n)
        out.append((cx + r * math.cos(t), cy + r * math.sin(t)))
    return out


def route_bulge(a, b, segments, keepouts, depths=ANNOTATE_ROUTE_DEPTHS,
                exempt=()):
    """(bulge, score) for the arc from a to b that stays clearest of a structure.

    `segments` are (i, j, p0, p1) from C.bond_segments; `keepouts` a list of
    ((x, y), radius) for the atoms' drawn forms. Ties go to the shallower arc, so
    a clean shallow bow always beats an equally clean deep one.

    `exempt` is the set of node ids the arrow is ANCHORED to. A bond with both
    ends in it is skipped: an arrow starting at a bond's midpoint touches that
    bond by construction, and without the exemption every candidate scores a
    crossing and the search has nothing to choose between.

    THE TWO ENDS ARE EXCLUDED from scoring. An arrow is supposed to begin and end
    beside the structure -- charging it for that would push every arrow out to
    infinity. Only the middle of the run has to be clear.
    """
    best, best_score = None, None
    skip = 2
    for depth in depths:
        for sgn in (1.0, -1.0):
            bulge = sgn * depth
            try:
                # NEGATED, exactly as _push_arrow negates it on the way into
                # arc_chord. Scoring the un-negated bulge scores the mirror image
                # of the arc that gets drawn -- measured: the router returned a
                # route with zero crossings and the figure came out crossing one.
                pts = arc_polyline(a, b, -bulge)
            except ValueError:
                continue
            score = ANNOTATE_ROUTE_DEPTH_COST * depth
            for k in range(skip, len(pts) - 1 - skip):
                p, q = pts[k], pts[k + 1]
                for i, j, r, t in segments:
                    if i in exempt and j in exempt:
                        continue
                    if _seg_cross(p, q, r, t):
                        score += ANNOTATE_ROUTE_CROSS
            near = ANNOTATE_ROUTE_NEAR * C.BOND_LENGTH
            for p in pts[skip:len(pts) - skip]:
                for (cx, cy), rad in keepouts:
                    if not rad:
                        continue
                    d = math.hypot(p[0] - cx, p[1] - cy)
                    if d < rad:
                        score += ANNOTATE_ROUTE_GRAZE * (1.0 - d / rad)
                for i, j, r, t in segments:
                    if i in exempt and j in exempt:
                        continue
                    d = _point_seg_dist(p, r, t)
                    if d < near:
                        score += ANNOTATE_ROUTE_NEAR_COST * (1.0 - d / near)
            if best_score is None or score < best_score:
                best, best_score = bulge, score
    return best, best_score


def _push_arrow(gid, a, b, bulge, head, colour, back=None) -> str:
    """One electron-pushing arrow, as a NATIVE <arrow> rather than a <curve>.

    WHY NOT A <curve>. `C.curve(..., ArrowheadHead=...)` -- a `<curve>` carrying
    an arrowhead -- draws the HEAD ALONE, with no line. Measured three ways:
      - four plausible curve encodings rendered: no line;
      - CurveType 0 and 8 both draw only the head;
      - a file with eleven curves with arrowheads has been through ChemDraw and
        none has a SupersededBy twin, i.e. ChemDraw never understood them (one
        twin appears per arrow it DOES understand).
    In an opened figure the arrows came back as loose red triangles.

    STILL INSIDE THE <fragment>, which the DTD does not allow (line 143 lists
    n | b | t | graphic | curve | ... and not arrow) -- a render test asked
    whether ChemDraw tolerates it anyway and it does, it draws. That
    tolerance is what keeps this a one-line swap: an annotation still rides along
    with C.place() in all six engines, and nothing has to be threaded through
    place_with the way a halo is.

    THE BULGE IS NEGATED. `bow()` pushes along +normal for a positive bulge and
    `arc_chord` bows the other way for the same number (measured: bow puts its
    controls at y=+35, arc_chord's arc reaches y=-35). Negating keeps every
    existing spec bowing the side it always meant.

    A circular arc, not the cubic `bow()` builds -- ChemDraw's arc arrow has no
    Bezier form; the rendered shape was checked and accepted.
    """
    # FLOOR THE SAG, then TRIM THE ENDS. ANNOTATE_SAG_FLOOR for the first,
    # C.arc_chord_trimmed for the second: `bulge` alone is a fraction of the
    # chord, so a half-bond chord bows by nothing, and an untrimmed arrow lands
    # on top of the very label it points at. Raise the bulge -- never lower it --
    # until the sag clears the floor, keeping the side it was bowing to.
    L = math.hypot(b[0] - a[0], b[1] - a[1])
    floor = ANNOTATE_SAG_FLOOR * C.BOND_LENGTH
    if L and bulge and abs(bulge) * L < floor:
        bulge = math.copysign(floor / L, bulge)
    back0, back1 = back or (0.0, 0.0)
    if L and (back0 or back1):
        a, b = apply_setback(a, b, (back0, back1))
        L = math.hypot(b[0] - a[0], b[1] - a[1])
        back0 = back1 = 0.0
    if L and not bulge:
        # BULGE 0 MEANS DRAW IT STRAIGHT, and a circular arc cannot express that
        # -- a zero sagitta is an infinite radius, and arc_circle rightly refuses.
        # It is a real input and not a degenerate one: a caller may flatten the
        # bow on purpose to drive an arrow through a caption, and routing that
        # through the arc would turn a drawn arrow into a raise.
        # On a straight run, setting back along the chord IS setting back along
        # the path, so the same two numbers apply with no trigonometry.
        room = L * (1.0 - ANNOTATE_INSET_KEEP)
        if back0 + back1 > room:
            k = room / (back0 + back1)
            back0, back1 = back0 * k, back1 * k
        ux, uy = (b[0] - a[0]) / L, (b[1] - a[1]) / L
        return GX.straight_arrow(
            gid, (a[0] + ux * back0, a[1] + uy * back0),
            (b[0] - ux * back1, b[1] - uy * back1),
            head=(head or None), head_type="Solid",
            triple=(C.ARC_HEAD_SIZE, C.ARC_HEAD_CENTER, C.ARC_HEAD_WIDTH),
            color=colour)
    xml = C.arc_chord(gid, a[0], a[1], b[0], b[1], bulge=-bulge,
                      head=(head or None),
                      head_type="Solid" if head else None)
    if colour is None:
        return xml
    return B._set_attr_on_open_tag(xml, "color", C.color_index(colour))


def bow(p0, p1, bulge):
    """Cubic Bezier control points for a chord from p0 to p1, bowed sideways.

    Four points, which is what C.curve wants. The two controls sit a third and
    two thirds along and are pushed out along the chord's normal, so the curve
    leaves and arrives roughly along it -- what a pushing arrow looks like when a
    chemist draws one, rather than a circular arc through both atoms.
    """
    (x0, y0), (x1, y1) = p0, p1
    dx, dy = x1 - x0, y1 - y0
    d = math.hypot(dx, dy) or 1.0
    nx, ny = -dy / d, dx / d          # unit normal
    off = bulge * d
    return [(x0, y0),
            (x0 + dx / 3.0 + nx * off, y0 + dy / 3.0 + ny * off),
            (x0 + 2 * dx / 3.0 + nx * off, y0 + 2 * dy / 3.0 + ny * off),
            (x1, y1)]


def _centre(points):
    return (sum(p[0] for p in points) / len(points),
            sum(p[1] for p in points) / len(points))


def _unit(dx, dy):
    d = math.hypot(dx, dy)
    return (dx / d, dy / d) if d else (0.0, -1.0)


# BEZIER PUSHING ARROWS. Four numbers decide the shape: the tangent direction at
# each end and the handle length at each end. The two DIRECTIONS are chemistry
# and are derived, not tuned -- an arrow leaves a lone pair along that lone pair,
# leaves a bond perpendicular to it, and arrives anti to the bond that is
# breaking. Only the LENGTHS are taste, and they are the one thing measured.
#
# WHERE THE TABLE COMES FROM. 1460 curved arrows read straight out of the vector
# PDF of a reference set of published mechanism figures (PyMuPDF, so these are
# the authors' own control points, not pixels). Of those, 531 are a
# SINGLE cubic -- the same two-knot shape ChemDraw's two-step pen tool draws, so
# the only apples-to-apples subset. Their medians, by chord length in bond
# lengths, with n:
#
#     chord      n     out    in    |angle out|
#     0.7-1.0    61   1.09  1.14     75 deg
#     1.0-1.5    76   0.98  0.92     68
#     1.5-2.5   105   0.86  0.78     51
#     2.5-4.0   140   0.75  0.72     45
#     4.0-7.0   109   0.57  0.57     44
#     7.0+       27   0.37  0.36     37
#
# THE FRACTION IS NOT A CONSTANT -- it runs 1.09 down to 0.37 -- so "handle =
# 0.60 * chord" is the wrong parameterisation and this has to be a table.
#
# THE SHAPE IS THE REFERENCE SET'S, THE LEVEL WAS CHOSEN. 0.60 / 0.50 was picked
# off a three-rung ladder on a 2.51-bond chord, where the reference set's own
# median is 0.75 / 0.72; that is just below the comparable group's lower quartile
# (0.64 / 0.61), i.e. inside the distribution on the gentle side. So the table is
# scaled by that ratio rather than replaced by a single number.
BEZ_HANDLE_BOOK = ((1.0, 1.09, 1.14), (1.5, 0.98, 0.92), (2.5, 0.86, 0.78),
                   (4.0, 0.75, 0.72), (7.0, 0.57, 0.57))
BEZ_HANDLE_LONG = (0.37, 0.36)              # past the last row
BEZ_HANDLE_SCALE = (0.60 / 0.75, 0.50 / 0.72)
# HEAD CLEARANCE, chosen off a four-rung ladder (0.20 / 0.35 / 0.50 / 0.65) with
# the handles fixed at the chosen level. Measured to the atom CENTRE; the useful
# number is the daylight to the nearest drawn ink, which at 0.65 is 0.56 bond.
# 0.20 is too tight -- the ink gap there is 0.17 bond, and the arrowhead is a
# solid triangle.
BEZ_HEAD_CLEARANCE = 0.65
# HOW FAR OFF THE CHORD A TANGENT POINTS -- a TARGET, not a ceiling. Read as a
# ceiling (clamp only when the derived direction points further round than this)
# it leaves nothing pushing the angle UP, and the result was measured: all four
# buckets of one test figure came out BELOW the reference set's lower quartile,
# one of them at 5 degrees where the reference arrows run 44-69. That is an
# overcorrection from fixing "very bold": the arrows stop swinging and go flat
# instead.
#
# So the split is: CHEMISTRY DECIDES THE SIDE, THE REFERENCE SET DECIDES HOW FAR. A
# derived direction (a lone pair's, a bond's outward normal) keeps the side it
# came out on and has its magnitude set to these medians. A direction an author
# NAMED -- `anti` for an Sn2 backside attack -- is a hard chemical constraint and
# is left alone to the degree.
#
# The numbers are the (out, in) medians of the reference set's 531 single-cubic
# arrows, by chord in bond lengths, same source as the handle table above.
BEZ_TANGENT_BOOK = ((1.0, 75.0, 85.0), (1.5, 68.0, 76.0), (2.5, 51.0, 49.0),
                    (4.0, 45.0, 42.0), (7.0, 44.0, 38.0))
BEZ_TANGENT_LONG = (37.0, 33.0)


def bez_tangent_aim(chord: float):
    """(out, in) TARGET in degrees off the chord, for a chord of `chord` units."""
    bl = chord / C.BOND_LENGTH
    for hi, a, b in BEZ_TANGENT_BOOK:
        if bl < hi:
            return a, b
    return BEZ_TANGENT_LONG


def bez_handles(chord: float):
    """(out, in) handle LENGTHS in drawing units for a chord of `chord` units."""
    bl = chord / C.BOND_LENGTH
    fo, fi = BEZ_HANDLE_LONG
    for hi, a, b in BEZ_HANDLE_BOOK:
        if bl < hi:
            fo, fi = a, b
            break
    return (fo * BEZ_HANDLE_SCALE[0] * chord,
            fi * BEZ_HANDLE_SCALE[1] * chord)


def lone_pair_dirs(i, pts, nbrs, n: int):
    """`n` angles (radians) for the electrons on node `i`: the widest GAPS between
    its bonds, widest first.

    Chemistry, not taste: a lone pair sits where there is no bond. An ether
    oxygen has two bonds and two gaps and gets one pair in each; a carbonyl
    oxygen has one bond, so the pairs spread over the open side; a phosphine's
    phosphorus has three bonds and one gap. Measured support for doing it this
    way rather than off the free valence: against the free-valence direction the
    lone pairs in published figures deviate by a median of 78 degrees, because a
    single bisector cannot describe two pairs on opposite sides.
    """
    angs = sorted(math.atan2(pts[j][1] - pts[i][1], pts[j][0] - pts[i][0])
                  for j in nbrs.get(i, ()) if j in pts)
    if not angs:                                   # a bare ion: spread evenly
        return [-math.pi / 2 + 2 * math.pi * k / max(n, 1) for k in range(n)]
    if len(angs) == 1:                             # one bond: fan the open side
        # 120 degrees apart for a pair, which is what a carbonyl or an alkoxide
        # oxygen is drawn with; 55 reads as pinched.
        base = angs[0] + math.pi
        step = min(math.radians(120.0), 2 * math.pi / (n + 1))
        return [base + step * (k - (n - 1) / 2.0) for k in range(n)]
    gaps = []
    for k in range(len(angs)):
        a0 = angs[k]
        w = (angs[(k + 1) % len(angs)] - a0) % (2 * math.pi)
        gaps.append((w, a0 + w / 2.0))
    gaps.sort(key=lambda g: -g[0])
    if n <= len(gaps):
        return [g[1] for g in gaps[:n]]
    w, mid = gaps[0]                               # more pairs than gaps
    step = min(math.radians(55.0), w / (n + 1))
    return [mid + step * (k - (n - 1) / 2.0) for k in range(n)]


def _bonded(frag: str, i: int, j: int) -> bool:
    """Is there a drawn bond between these two nodes right now?"""
    for a, b, _p0, _p1 in C.bond_segments(frag):
        if {a, b} == {i, j}:
            return True
    return False


def _norm(v):
    L = math.hypot(v[0], v[1])
    return (v[0] / L, v[1] / L) if L else (1.0, 0.0)


def _aim_off_chord(d, chord_dir, target):
    """Rotate `d` to exactly `target` degrees off the chord, keeping its SIDE.

    Which side is the chemistry -- a lone pair points where the atom has no bond,
    a bond's normal points away from the structure -- and that survives. Only the
    magnitude is taken from the reference set.
    """
    d, c = _norm(d), _norm(chord_dir)
    cross = c[0] * d[1] - c[1] * d[0]
    t = math.radians(target) * (1.0 if cross >= 0 else -1.0)
    return (c[0] * math.cos(t) - c[1] * math.sin(t),
            c[0] * math.sin(t) + c[1] * math.cos(t))


def _end_dir(ids, pts, nbrs, toward=None, centre=None):
    """The tangent at one end, before clamping. Chemistry, three cases.

    `toward` is a node id an author named explicitly: this end's tangent points at
    it. That is the Sn2 case -- the head arrives at the carbon travelling toward
    where the leaving group sits, which is the backside, and no rule about the
    drawing could have guessed it.
    """
    here = _centre([pts[i] for i in ids])
    if toward is not None and toward in pts:
        return (pts[toward][0] - here[0], pts[toward][1] - here[1])
    if len(ids) >= 2:
        p0, p1 = pts[ids[0]], pts[ids[1]]
        n = (-(p1[1] - p0[1]), p1[0] - p0[0])          # normal to the bond
        if centre is not None:
            if (here[0] - centre[0]) * n[0] + (here[1] - centre[1]) * n[1] < 0:
                n = (-n[0], -n[1])                      # the outward side
        return n
    i = ids[0]
    sx = sy = 0.0
    for j in nbrs.get(i, ()):
        if j in pts:
            u = _norm((pts[j][0] - pts[i][0], pts[j][1] - pts[i][1]))
            sx += u[0]
            sy += u[1]
    return (-sx, -sy) if (sx or sy) else (0.0, -1.0)


def _push_bezier(gid, pts, nbrs, tail_ids, head_ids, kinds, head, colour,
                 toward_tail=None, toward_head=None, centre=None):
    """The whole arrow, from the two ends' node ids. Returns the <curve> string."""
    a0 = _centre([pts[i] for i in tail_ids])
    b0 = _centre([pts[i] for i in head_ids])
    # NORMALISED HERE, at the one place both come from. _end_dir returns a plain
    # difference of two atom positions for a NAMED direction, i.e. about one bond
    # long, and the setback multiplies whatever it is given: unnormalised, an Sn2
    # arrow's head landed 11 bond lengths away -- measured as a chord that jumped
    # from 1.22 to 11.44 bond.
    raw_out = _norm(_end_dir(tail_ids, pts, nbrs, toward_tail, centre))
    raw_in = _norm(_end_dir(head_ids, pts, nbrs, toward_head, centre))
    if toward_head is None and len(head_ids) == 1:
        raw_in = (-raw_in[0], -raw_in[1])   # a free valence points OUT, the head arrives IN
    # SETBACK ALONG THE TANGENT, not along the chord. The causality is the same
    # -- ends first, curve second -- but the end knows its own direction, so
    # backing off along it keeps the arrow on its approach line instead of
    # sliding it toward the other end. That is also why ANNOTATE_INSET_KEEP is
    # not used here: setting back along the chord, a bond-to-its-own-atom arrow
    # loses 65% of its chord and comes back as a 326-degree loop.
    ins = ANNOTATE_INSET.get(kinds[0], 0.2) * C.BOND_LENGTH
    out = (BEZ_HEAD_CLEARANCE if len(head_ids) == 1
           else ANNOTATE_INSET["bond"]) * C.BOND_LENGTH
    # THE HOOK, and it is its own case because it is half of every mechanism drawn:
    # the arrow starts on a bond and lands on one of THAT BOND'S OWN atoms -- a
    # pi bond collapsing onto its oxygen, a sigma bond onto the atom that leaves.
    # Deriving the head's direction from the atom's free valence puts the two
    # tangents nearly opposed, and then each end's setback pushes past the other:
    # measured, the chord collapsed to 0.15 bond and the tangent came out 173
    # degrees off it. Iterating the cap made it worse, because the feedback is the
    # problem, not the starting point.
    #
    # What published mechanisms draw here is a small hook sitting OUTSIDE the
    # bond: out along the bond's outward normal, over, and back in onto the
    # atom. So both tangents
    # are that one normal, +n leaving and -n arriving, and both ends are displaced
    # along +n -- which leaves the chord running along the bond at its full half
    # bond instead of eating it.
    shared = set(tail_ids) & set(head_ids)
    bond_end = tail_ids if len(tail_ids) >= 2 else (
        head_ids if len(head_ids) >= 2 else None)
    if shared and bond_end is not None:
        p0, p1 = pts[bond_end[0]], pts[bond_end[1]]
        n = _norm((-(p1[1] - p0[1]), p1[0] - p0[0]))
        if centre is not None:
            mid = ((p0[0] + p1[0]) / 2.0, (p0[1] + p1[1]) / 2.0)
            if (mid[0] - centre[0]) * n[0] + (mid[1] - centre[1]) * n[1] < 0:
                n = (-n[0], -n[1])
        a = (a0[0] + n[0] * ins, a0[1] + n[1] * ins)
        b = (b0[0] + n[0] * out, b0[1] + n[1] * out)
        # Both ends move along the SAME normal, so the chord is known in closed
        # form and the aim needs no second pass.
        chord_dir = (b[0] - a[0], b[1] - a[1])
        aim_out, aim_in = bez_tangent_aim(math.hypot(*chord_dir))
        d_out = _aim_off_chord(n, chord_dir, aim_out)
        d_in = _aim_off_chord((-n[0], -n[1]), chord_dir, aim_in)
    else:
        chord_dir = (b0[0] - a0[0], b0[1] - a0[1])
        aim_out, aim_in = bez_tangent_aim(math.hypot(*chord_dir))
        # A NAMED direction is chemistry to the degree; only a derived one is
        # re-aimed. `anti` says the head arrives opposite the leaving group, and
        # rotating that to the reference median would be drawing a different
        # mechanism.
        d_out = raw_out if toward_tail is not None else _aim_off_chord(
            raw_out, chord_dir, aim_out)
        d_in = raw_in if toward_head is not None else _aim_off_chord(
            raw_in, chord_dir, aim_in)
        a = (a0[0] + d_out[0] * ins, a0[1] + d_out[1] * ins)
        b = (b0[0] - d_in[0] * out, b0[1] - d_in[1] * out)
    l0, l1 = bez_handles(math.hypot(b[0] - a[0], b[1] - a[1]))
    return C.bezier_arrow(gid, a, d_out, l0, b, d_in, l1,
                          head=head, color=colour)


# WHICH COPY, decided ONCE for the whole structure. `node_ids_by_atom_map` says
# exactly which atom a map number means; it cannot say which of four equivalent
# bromines the packer put within reach. Chemically they are one atom, so any is
# correct and the readable one is right -- but two arrows of the SAME step often
# name the same bromine from two separate specs (P attacks Br; the Br-C pair
# leaves), and optimising each spec on its own would let them choose DIFFERENT
# bromines and draw a mechanism nobody meant. So the substitution is one choice
# shared by every map-addressed spec, scored on the total chord.
MAP_SUB_COMBOS = 256      # past this, keep the atoms as written rather than search


_MAPPED_ATOM = re.compile(r"\[([^\]]*?):(\d+)\]")


def _refuse_accidental_radicals(smiles: str) -> None:
    """A BRACKET ATOM HAS NO IMPLICIT HYDROGENS, and writing a map number forces
    brackets. So `[C:2]1(CCCCC1)[O:3]` is not a CH of a ring -- it is a carbon
    with three bonds, no hydrogen and therefore ONE RADICAL ELECTRON, which the
    engine then dutifully stamps a dot on -- which reads as a mistake in the
    chemistry, because it is one.
    RDKit is right and the SMILES was wrong, so this refuses it by name rather
    than drawing it. Write `[CH:2]`; write `[CH0:2]` if a radical is meant.
    """
    mol = S.parse_subject(smiles)
    if mol is None:
        return
    spelled = {int(n): txt for txt, n in _MAPPED_ATOM.findall(smiles or "")}
    for a in mol.GetAtoms():
        n = a.GetAtomMapNum()
        if not n or not a.GetNumRadicalElectrons():
            continue
        if "H" in spelled.get(n, ""):
            continue                       # said so explicitly, leave it alone
        raise ValueError(
            "atom :%d (%s) came out with %d radical electron(s): a bracket atom "
            "carries no implicit hydrogens, and writing :%d forces brackets. "
            "Write [%sH:%d] if you meant a hydrogen there, or [%sH0:%d] if the "
            "radical is deliberate."
            % (n, a.GetSymbol(), a.GetNumRadicalElectrons(), n,
               a.GetSymbol(), n, a.GetSymbol(), n))


def _map_substitution(frag, smiles, index_map, base_id, specs, where):
    """{map number -> node id} for this structure, or None if nothing needs it."""
    live = [sp for sp in specs
            if sp.get("style") == "annotate" and not sp.get("on")]
    # A maps-addressed spec of any style needs the same table; only the ARROWS
    # get a vote on which copy, because only they have two ends to measure.
    if not live and not any(sp.get("maps") for sp in specs):
        return None
    written = S.node_ids_by_atom_map(smiles, index_map, base_id=base_id,
                                     where=where)
    _refuse_accidental_radicals(smiles)
    used = sorted({int(v) for sp in live
                   for v in list(sp.get("tail") or [1]) + list(sp.get("head") or [2])})
    missing = [n for n in used if n not in written]
    if missing:
        raise ValueError(
            "annotate names atom map(s) %s, which the subject %r does not carry; "
            "it maps %s" % (missing, smiles, sorted(written) or "nothing"))
    classes = S.atom_map_classes(smiles, index_map, base_id=base_id)
    free = [n for n in used if len(classes.get(n, ())) > 1]
    if not free:
        return written
    combos = 1
    for n in free:
        combos *= len(classes[n])
    if combos > MAP_SUB_COMBOS:
        return written
    pts = C.node_points(frag)

    def cost(cand):
        total = 0.0
        for sp in live:
            t = [cand[int(v)] for v in (sp.get("tail") or [1])]
            h = [cand[int(v)] for v in (sp.get("head") or [2])]
            if not all(i in pts for i in t + h):
                continue
            a, b = _centre([pts[i] for i in t]), _centre([pts[i] for i in h])
            total += math.hypot(b[0] - a[0], b[1] - a[1])
        return total

    best, best_cost = written, None
    for pick in itertools.product(*[classes[n] for n in free]):
        cand = dict(written)
        cand.update(dict(zip(free, pick)))
        if len(set(cand.values())) != len(cand):
            continue                      # two map numbers collapsed onto one atom
        c = cost(cand)
        if best_cost is None or c < best_cost - 1e-9:
            best, best_cost = cand, c
    return best


def apply(frag: str, smiles: str, index_map: dict, specs, *,
          first_id: int | None = None, base_id: int = 1,
          where: str | None = None):
    """Resolve every spec against `frag` and return (decorations, frag).

    IDS COME FROM THE FRAGMENT'S OWN RANGE, which is `next_free_id`'s rule and the
    reason it exists: "taken from the fragment's own range, an id stays unique for
    exactly as long as the atom ids do". Two fixed pools (say 7000 for graphics,
    7500 for text) that reset on every call would be wrong twice over:

      * N decorated structures in ONE call would each start at 7000, so they
        would all carry the same numbers.
      * 7000 IS THE SEVENTH CELL'S BASE. Every engine hands structures ids with
        `base += 1000`, so cell 7's atoms are 7001, 7002, ... Measured on ten
        structures with six dots on the first: `<graphic id="7001">` and
        `<n id="7001">` in one document. Moving the pool higher only moves the
        threshold -- the number of structures in a figure has no upper bound,
        so any fixed pool is eventually inside somebody's lane.

    ONE COUNTER, not two. Separate pools for graphics and text, 500 apart, would
    be a second fixed gap to be outgrown -- one decoration set emits
    one graphic per SMARTS match, and that count is unbounded too.

    `first_id` is for a caller that already knows the range; None means ask the
    fragment. No caller passes it at present, and place_with renumbers what
    comes back by the same `base` as the fragment, so the two stay together.

    `decorations` is one string, already in back-to-front order: glows first,
    then dots and position numbers. It all belongs BEHIND the structure -- emit
    it, then the fragment, then stamp_z over both -- and BEHIND is correctness, not
    polish: the glow is opaque, measured in ChemDraw (a green halo over a benzene
    hid the ring). EMPTY for the styles that only
    rewrite the fragment ("color", "fill"). Those two needed nothing from a
    caller; the rest need a `decorations` sink and a place_with at the other end,
    which every engine has.

    `frag` comes back modified by "color" (tints atoms and bonds) and by "fill"
    (adds a native ColoredMolecularArea inside the fragment).

    `base_id` must match the one fragment() was built with, or the ids an anchor
    resolves to will not be the ids in the string. `where` labels this structure
    in the miss report.

    Each spec: {"on": <SMARTS>, "style": ..., "color": ..., "labels": [...]}
    """
    behind, front = [], []
    # One past this fragment's own largest id. `frag` is still at its build-time
    # numbering here (place() has not run), so everything emitted below numbers
    # itself just above the atoms and travels with them.
    nid = C.next_free_id(frag) if first_id is None else int(first_id)
    _sub = _map_substitution(frag, smiles, index_map, base_id, specs or [], where)
    for spec in (specs or []):
        style = spec.get("style", "color")
        if style not in STYLES:
            raise ValueError(f"unknown decoration style {style!r}; "
                             f"known: {', '.join(STYLES)}")
        on = spec.get("on")
        # NO `on` IS A REAL SPEC, for exactly one style. An annotate arrow may
        # name its two ends through the SUBJECT's own :n rather than through a
        # query -- see node_ids_by_atom_map for why a query cannot say WHICH copy.
        # Every other style still needs an anchor, and a spec without one is not
        # dropped here without a word: a figure missing an arrow is a different
        # mechanism, so it raises instead.
        if not on and not spec.get("maps") and not (
                style == "annotate" and (spec.get("tail") or spec.get("head"))):
            raise ValueError(
                "decoration %r has no anchor. Give 'on' (a SMARTS), or 'maps' "
                "(atom map numbers written in the subject SMILES); an annotate "
                "arrow may instead give 'tail'/'head' as map numbers."
                % (spec,))
        colour = spec.get("color", "red")

        if style == "annotate":
            # Resolved here rather than through match_bond_ids: that one is about
            # a DRAWN bond and reports bond-shaped failures ("both ends collapsed
            # into one abbreviation"), while two anchors here are just two atoms
            # and are routinely not bonded at all -- an electron-pushing arrow
            # goes from a lone pair to an empty orbital.
            pts = C.node_points(frag)
            # ENDS BY MAP NUMBER, not by match order. `grp[0], grp[-1]` below is
            # an accident of how the query happened to be written -- it can put a
            # documented ":1 tail, :2 head" on the last ring atom instead.
            # An arrow has a tail and a head, and half of every mechanism arrow
            # ever drawn has a BOND at one end (a pi bond opening, a lone pair
            # closing onto one), which two atom slots cannot say at all.
            #   "tail"/"head": lists of :n map numbers. ONE number is that atom;
            #   TWO is the midpoint of those two, i.e. the bond between them.
            # Default [1] -> [2]: exactly what the docstring promises.
            tail_n = [int(v) for v in (spec.get("tail") or [1])]
            head_n = [int(v) for v in (spec.get("head") or [2])]
            if on:
                ids = S.match_node_ids(smiles, on, index_map, base_id=base_id,
                                       where=where)
                mapped = S.match_node_ids_by_map(smiles, on, index_map,
                                                 base_id=base_id, where=where)
            else:
                # THE SUBJECT'S OWN :n. One subject, so exactly one "match":
                # nothing to uniquify, and nothing for `nth` to choose between,
                # which is the whole reason this path exists.
                # _sub already resolved -- and already checked that every map
                # number named here exists -- for the whole structure at once.
                mapped = [_sub]
                ids = [[_sub[n] for n in tail_n + head_n]]
            raw = [[pts[i] for i in grp if i in pts] for grp in ids]
            # WHAT EACH END HAS TO CLEAR. props gives the element and the charge
            # of the node an end landed on, which is the whole of what decides
            # the setback: a drawn label needs the most room, a bond less, a bare
            # carbon vertex least. See ANNOTATE_INSET.
            props = C.node_props(frag)
            ends = []
            for k in range(len(raw)):
                m = mapped[k] if k < len(mapped) else None
                if not m or not all(n in m and m[n] in pts
                                    for n in tail_n + head_n):
                    ends.append(None)
                    continue

                def _kind(ns, _m=m):
                    if len(ns) != 1:
                        return endpoint_kind(natoms=len(ns))
                    pr = props.get(_m[ns[0]], {})
                    return endpoint_kind(int(pr.get("Element", 6)),
                                         int(pr.get("Charge") or 0),
                                         node_type=pr.get("NodeType"))

                # The fourth item is the node ids this arrow is ANCHORED to, for
                # route_bulge's exemption: a bond it starts or ends on is one it
                # is allowed to touch.
                # A BOND THAT DOES NOT EXIST YET HAS NO MIDPOINT WORTH POINTING
                # AT. Naming the sink as the bond that FORMS is the right way to
                # say it -- replay needs both ends to do the arithmetic -- but
                # the two atoms are typically a whole molecule apart at this
                # moment, so their midpoint is empty space. Measured on one
                # published mechanism: three arrows out of fourteen landed on an
                # unrelated substituent (an R4, a ring oxygen, a nitrogen) rather
                # than on the atom being attacked. Published mechanisms draw this
                # arrow AT the atom, so when the named bond is not drawn, the head
                # goes to whichever of its two atoms is not the source.
                _h = list(head_n)
                if len(_h) == 2 and not _bonded(frag, m[_h[0]], m[_h[1]]):
                    _out = [n for n in _h if m[n] not in {m[t] for t in tail_n}]
                    if len(_out) == 1:
                        _h = _out
                ends.append((_centre([pts[m[n]] for n in tail_n]),
                             _centre([pts[m[n]] for n in _h]),
                             endpoint_backs(_kind(tail_n), _kind(_h)),
                             {m[n] for n in tail_n + _h},
                             tuple(m[n] for n in tail_n),
                             tuple(m[n] for n in _h),
                             (_kind(tail_n), _kind(_h))))
            keep = [(g, e) for g, e in zip(raw, ends) if g]
            # ONE MATCH, NAMED. Annotating every match is right for a scope table
            # -- one spec, sixteen cells, each gets its own copy. It is wrong for
            # a mechanism on a SYMMETRIC substrate: `[cH:1][cH:2]` on benzene
            # matches six times, so the one pi bond attacking the nitronium comes
            # out as six arrows. SMARTS cannot say "whichever one", because graph
            # matching is symmetry-blind by construction.
            # So: "nth" (0-based) draws only that match. On a symmetric ring every
            # match is chemically the same arrow, so 0 is as good as any; bump it
            # when the layout puts that copy somewhere ugly. Out of range RAISES
            # rather than drawing nothing -- a silently missing arrow in a
            # mechanism is a different mechanism.
            nth = spec.get("nth")
            if nth is not None:
                nth = int(nth)
                if not -len(keep) <= nth < len(keep):
                    raise ValueError(
                        "annotate nth=%d but the anchor %r matched %d time(s)"
                        % (nth, on, len(keep)))
                keep = [keep[nth]]
            here = [g for g, _ in keep]
            ends = [e for _, e in keep]
            if not here:
                continue
            text = spec.get("text")
            off = spec.get("offset")
            head = ANNOTATE_ARROWS.get(str(spec.get("arrow", "none")).lower(),
                                       "MISSING")
            if head == "MISSING":
                raise ValueError(
                    "annotate arrow=%r; one of %s. 'half' is the fishhook a "
                    "single-electron arrow draws."
                    % (spec.get("arrow"), ", ".join(sorted(ANNOTATE_ARROWS))))
            # "bulge": "auto" hands the side and the depth to route_bulge, which
            # picks whichever arc stays clearest of the drawn structure. OPT-IN,
            # not the default: a spec that passes no bulge gets ANNOTATE_BULGE,
            # and silently re-routing every such arrow is a change nobody asked
            # for.
            _b = spec.get("bulge", ANNOTATE_BULGE)
            auto = isinstance(_b, str) and _b.strip().lower() == "auto"
            bulge = ANNOTATE_BULGE if auto else float(_b)
            if auto:
                _segs = C.bond_segments(frag)
                _kos = [(pts[i], r) for i, r in C.keepout_radii(frag).items()
                        if i in pts]
            size = int(spec.get("size", ANNOTATE_TEXT_SIZE))
            _nbrs: dict = {}
            for _i, _j, _p0, _p1 in C.bond_segments(frag):
                _nbrs.setdefault(_i, []).append(_j)
                _nbrs.setdefault(_j, []).append(_i)
            _frag_centre = _centre(list(pts.values())) if pts else None

            def _toward(v, _s=_sub):
                """A map number an author gave for a tangent -> its node id."""
                if v is None or _s is None:
                    return None
                n = int(v[0] if isinstance(v, (list, tuple)) else v)
                if n not in _s:
                    raise ValueError(
                        "tail_dir/head_dir names atom map :%d, which the subject "
                        "does not carry" % n)
                return _s[n]

            for grp, end in zip(here, ends):
                # one anchor -> its atom (or a bond's midpoint); two -> both ends
                # `backs` stays None on the two legacy paths: a SMARTS with no
                # :n maps at all, and a captioned leader. Neither knows which
                # atom is which end, so neither can say what to clear, and
                # leaving them alone keeps their output unchanged.
                backs = None
                own = frozenset()
                tids = hids = ()
                kinds = ("carbon", "carbon")
                if end is not None and text is None:
                    a, b, backs, own, tids, hids, kinds = end
                elif len(grp) >= 2 and text is None:
                    a, b = grp[0], grp[-1]
                else:
                    a = b = _centre(grp)
                if text is not None and off:
                    # leader: from the anchor out to where the caption sits
                    b = (a[0] + float(off[0]), a[1] + float(off[1]))
                parts = []
                if a != b:
                    # The setback runs HERE and not inside _push_arrow when the
                    # bulge is being chosen, because the router needs the FINAL
                    # endpoints: pick the side and depth from where the arrow
                    # actually starts and stops, not from where its targets sit.
                    bulge_here, back_here = bulge, backs
                    if auto:
                        if backs:
                            a, b = apply_setback(a, b, backs)
                        back_here = None
                        bulge_here = route_bulge(a, b, _segs, _kos,
                                                 exempt=own)[0] or bulge
                    if head and text is None and tids and hids:
                        # EVERY PUSHING ARROW IS A BEZIER. An arc is of little
                        # use for arrow pushing, and the measurements agree: an
                        # arc has one degree of freedom, so its tangents fall out
                        # of the two endpoints, and 25% of the reference set's
                        # 1460 push arrows carry an inflection no conic arc can
                        # have.
                        #
                        # STILL AN ARC, and deliberately, are the two things that
                        # are not pushing arrows: the leader out to a caption
                        # (`text` with an `offset`) and the HEADLESS arc that
                        # marks where a minor regioisomer connects (`arrow:
                        # "none"`, so `head` is None). Both were calibrated as
                        # arcs and neither moves electrons.
                        #
                        # The remaining fallback is a SMARTS with no :n maps at
                        # all: it yields two points and no node ids, and the
                        # tangents need to know which node is which end.
                        parts.append(_push_bezier(
                            C.next_free_id(frag), pts, _nbrs, tids, hids,
                            kinds, head, colour,
                            toward_tail=_toward(spec.get("tail_dir")),
                            toward_head=_toward(spec.get("head_dir")),
                            centre=_frag_centre))
                    else:
                        parts.append(_push_arrow(C.next_free_id(frag), a, b,
                                                 bulge_here, head, colour,
                                                 back=back_here))
                if text is not None:
                    tx, ty = (b if off else (a[0] + ANNOTATE_TEXT_PAD, a[1]))
                    parts.append(C.styled_text(
                        C.next_free_id(frag) + len(parts), tx, ty, [str(text)],
                        size=size, align="center", color=C.color_index(colour)))
                # INSIDE the fragment, both of them. The DTD's content model for
                # <fragment> lists t and curve alongside n and b, so an annotation
                # is part of the molecule object and C.place() carries it, so this
                # style needed nothing from the engines. (halo / dot / locants are
                # NOT inside the fragment and need exactly that; they ride
                # along too, via cdxml_core.place_with.) Same move
                # squiggle made, and the same last-</fragment> rule: a superatom
                # nests one,
                # and inserting at the first would bury the mark in the label.
                k = frag.rindex("</fragment>")
                frag = frag[:k] + "".join(parts) + frag[k:]
            continue

        if style == "delocalized":
            # A CHARGE SPREAD OVER SEVERAL ATOMS, drawn the way a hand-drawn
            # reference structure draws it: the double bonds come OUT of the
            # delocalized part and a dashed arc plus one circled symbol go in. See C.delocalized_marks
            # for where every number came from.
            #
            # The anchor names the DELOCALIZED ATOMS, not the whole ring -- for an
            # arenium that is the five carbons and not the sp3 one, which is also
            # what lets the arc find its own ends (it runs the way that avoids the
            # largest angular gap, and that gap is the one across the sp3 carbon).
            chg = spec.get("charge", 0)
            if isinstance(chg, str):
                chg = {"+": 1, "-": -1, "": 0}.get(chg.strip(), 0)
            chg = int(chg)
            pts = C.node_points(frag)
            for grp in S.match_node_ids(smiles, on, index_map, base_id=base_id,
                                        where=where):
                here = [pts[i] for i in grp if i in pts]
                if len(here) < 3:
                    raise ValueError(
                        "delocalized anchor %r matched %d atom(s); a delocalized "
                        "system needs at least three, and the arc has no ends to "
                        "find with fewer" % (on, len(here)))
                frag = C.singleize_bonds(frag, grp)
                # AND take the charge off the atoms. The hand-drawn reference
                # carries no Charge on any node: the charge is the floating
                # symbol, full stop. Left on, _charge_symbols has already bound a Plus to the
                # cationic carbon and the figure states the charge twice.
                frag = C.clear_atom_charges(frag, grp)
                DELOCALIZED.append(where or on)
                aid = C.next_free_id(frag)
                # "center": "ring" (the DEFAULT) | "system".
                #
                # The rule -- "the positive charge sits at the exact centre of the
                # six-membered ring" -- beats IUPAC
                # GR-5.4's "centered within the delocalized system" as the
                # default here for a plain reason: the two differ by 0.2 bond on
                # an arenium, and centring on the five carbons puts the symbol
                # visibly off the ring's middle, which is what a reader checks.
                # "system" stays available and is what a NON-CYCLIC delocalized
                # system (an allyl cation) has to use, since it has no ring.
                #
                # The arc gets the same centre. Centring the arc on the five
                # carbons while the ring is drawn around a different point makes
                # the dashes hug one side of the ring.
                centre = None
                apothem = 0.0
                if str(spec.get("center", "ring")).lower() == "ring":
                    ring = S.enclosing_ring_node_ids(smiles, grp, index_map,
                                                     base_id=base_id)
                    if not ring:
                        raise ValueError(
                            'delocalized center="ring" but the anchor %r is not '
                            "contained in any one ring, so there is no ring "
                            'centre to use -- pass center="system"' % (on,))
                    rp = [pts[i] for i in ring if i in pts]
                    centre = (sum(q[0] for q in rp) / len(rp),
                              sum(q[1] for q in rp) / len(rp))
                    # The APOTHEM: how far the ring's own bond lines run from that
                    # centre. Scaling this rather than the atom distance is what
                    # turns "does the arc clear the ring" from luck into a number.
                    rset = set(ring)
                    apothem = min(
                        (_point_seg_dist(centre, p0, p1)
                         for i, j, p0, p1 in C.bond_segments(frag)
                         if i in rset and j in rset), default=0.0)
                marks = C.delocalized_marks(
                    here, aid, aid + 1, charge=chg, centre=centre,
                    apothem=apothem, radius_frac=spec.get("radius"),
                    color=(spec.get("color") if "color" in spec else None))
                k = frag.rindex("</fragment>")
                frag = frag[:k] + "".join(marks) + frag[k:]
            continue

        if style in BOND_STYLES:
            # Handled up here with squiggle, and for the first of its two reasons:
            # it names a BOND, so match_bond_ids is the resolver, and falling
            # through to the atom resolver as well would file a second miss report
            # for the same anchor.
            #
            # Unlike squiggle it adds no object at all -- it rewrites an attribute
            # on a <b> that fragment() already wrote. Nothing to place, nothing to
            # reserve space for, and it survives place() because the bond does.
            pairs = list(S.match_bond_ids(smiles, on, index_map,
                                          base_id=base_id, where=where,
                                          style=style))
            if pairs:
                frag, gone = B.set_bond_attrs(frag, pairs, BOND_STYLES[style],
                                              orient=style in _ORIENTED)
                for _bid, attr, was, now in gone:
                    S.BOND_STYLE_OVERWROTE.append(
                        (S._subject(smiles, where), on, style,
                         "%s %s -> %s" % (attr, was, now)))
            continue

        if style == "squiggle":
            # THE ODD ONE OUT, twice over, which is why it is handled before the
            # atom resolver rather than as another arm of the ladder below.
            #
            # It names a BOND, so it needs match_bond_ids -- calling the atom
            # resolver as well would file a second miss report for the same
            # anchor.
            #
            # And it lives INSIDE the fragment, as "fill" does. The CDX spec lists
            # Fragment as a legal container for <curve>, and translate() moves
            # CurvePoints (hand-drawn metallocene arcs need it) -- so place()
            # carries a squiggle along with the structure, and this style reaches
            # every figure tool for free. halo / dot / locants are separate
            # objects and do not; they reach every tool through
            # cdxml_core.place_with, which moves them by the structure's own
            # shift. Being inside the fragment is still the cheaper of the two.
            amp = float(spec.get("amp", SQUIGGLE_AMP)) * C.BOND_LENGTH
            wl = float(spec.get("wavelength", SQUIGGLE_WAVELENGTH)) * C.BOND_LENGTH
            over = float(spec.get("overhang", SQUIGGLE_OVERHANG)) * C.BOND_LENGTH
            if amp <= 0 or wl <= 0:
                raise ValueError(
                    "squiggle 'amp' and 'wavelength' are fractions of a bond "
                    "length and must be positive; got amp=%r wavelength=%r"
                    % (spec.get("amp", SQUIGGLE_AMP),
                       spec.get("wavelength", SQUIGGLE_WAVELENGTH)))
            pts = C.node_points(frag)
            fbb = C.bbox_of(frag)
            for ida, idb in S.match_bond_ids(smiles, on, index_map,
                                             base_id=base_id, where=where):
                if ida not in pts or idb not in pts:
                    continue
                (ax, ay), (bx, by) = pts[ida], pts[idb]
                cpts = C.squiggle_points(ax, ay, bx, by, amp=amp,
                                         wavelength=wl, overhang=over)
                # LAST </fragment>, the same as "fill": a superatom nests one, and
                # the first would put the mark inside the abbreviation.
                k = frag.rindex("</fragment>")
                frag = (frag[:k]
                        + C.curve(C.next_free_id(frag), cpts, color=colour)
                        + frag[k:])
                if fbb:
                    ink = C.curve_polyline(cpts)
                    out = max(fbb[0] - min(p[0] for p in ink),
                              fbb[1] - min(p[1] for p in ink),
                              max(p[0] for p in ink) - fbb[2],
                              max(p[1] for p in ink) - fbb[3])
                    if out > 0:
                        S.SQUIGGLE_OVERFLOW.append(
                            (S._subject(smiles, where), on, out))
            continue

        if on:
            matches = S.match_node_ids(smiles, on, index_map, base_id=base_id,
                                       where=where)
        else:
            # "maps": the subject's own :n, for every style that takes an atom
            # LIST -- color, fill, halo, dot, locants, lonepair. Same address the
            # arrows use, and deliberately through the SAME substitution: colour
            # the bond that broke and you want the bromine the arrow actually
            # chose, not the one that happened to be written first.
            mp = _sub if _sub is not None else S.node_ids_by_atom_map(
                smiles, index_map, base_id=base_id, where=where)
            want = [int(v) for v in spec["maps"]]
            gone = [n for n in want if n not in mp]
            if gone:
                raise ValueError(
                    "%s names atom map(s) %s, which the subject %r does not "
                    "carry; it maps %s"
                    % (style, gone, smiles, sorted(mp) or "nothing"))
            matches = [[mp[n] for n in want]]

        if style == "color" and len(matches) > 1:
            # ONE REGION, not one per match. color_region paints a bond only when
            # BOTH its ends are in the list it is handed, so a pattern whose
            # matches are single atoms, painted match by match, colours NO bond at
            # all: on="*" over a 12-atom molecule tinted 2 atoms and 0 bonds, and
            # on="[#6]" tinted nothing whatsoever -- a carbon has no label, so its
            # colour shows only through its bonds. Neither says a word, because
            # the "matched nothing" note is false here: "*" matched twelve times.
            #
            # "the atoms and the bonds between them" is what the docstring
            # promises and what a chemist asking for a tinted skeleton means, and
            # `them` is every matched atom. Unioning is also what makes a common
            # published style expressible at all -- every structure drawn in one
            # colour with the reactive site in another.
            #
            # ONLY for "color". The others are per-match by nature: two halos are
            # two rings, two dots are two dots, and locants consume a match's
            # ORDER.
            _seen, _one = set(), []
            for _ids in matches:
                for _i in _ids:
                    if _i not in _seen:
                        _seen.add(_i)
                        _one.append(_i)
            matches = [_one]

        # Recorded here, AFTER the color merge: by this point color always has
        # exactly one region, which is the behaviour it promises and should not be
        # reported. The remaining styles are drawn once per match, and dot /
        # locants / lonepair once more per atom in each match -- two multipliers,
        # so both counts are recorded.
        #
        # Both conditions are required; chosen after measuring 110 published
        # figures:
        #   matches>1 alone -> 20 reports / 8 figures, 13 of them [#0] placed on
        #                      every R position, which is DELIBERATE, so reporting
        #                      it is noise;
        #   plus atoms>matches -> 4 reports / 3 figures, including one real bug
        #                      (locant "10" drawn on two different carbons: 2
        #                      occurrences in the file, against 1 each for 19/17/8).
        # The matches>1 condition cannot be dropped: one match naming a whole
        # group is the normal case for halo / fill (a coumarin core, 10 atoms, one
        # match), and dropping it would report every halo.
        if on and len(matches) > 1:
            _n_at = sum(len(_m) for _m in matches)
            if _n_at > len(matches):
                _row = (S._subject(smiles, where), on, style, len(matches), _n_at)
                if _row not in S.ANCHOR_MULTI:
                    S.ANCHOR_MULTI.append(_row)

        pts = C.node_points(frag)

        for ids in matches:
            here = [pts[i] for i in ids if i in pts]
            if not here:
                continue

            if style == "lonepair":
                # NO "radical" SIBLING, deliberately. The engine already stamps
                # an Electron for every atom RDKit gives radical electrons to
                # (smiles_to_frag, GetNumRadicalElectrons), so a style for it
                # would draw the dot TWICE (only one of the two graphics carries
                # IgnoreWarnings). Lone pairs are the real gap: nothing else emits
                # them, and 53% of a reference set of published mechanism
                # figures draw them.
                # INSIDE the fragment, like annotate and for the same reason:
                # C.place() carries a fragment's contents, so this works in all
                # six engines with nothing asked of them. ChemDraw's own samples
                # write these at BOTH page level and fragment level, so depth 1
                # has precedent (Carbon Hybridization.cdxml: 28 of them).
                n_e = int(spec.get("n", 1))
                if n_e < 1:
                    raise ValueError("%s n=%d; it must be at least 1"
                                     % (style, n_e))
                _nb: dict = {}
                for _i, _j, _p0, _p1 in C.bond_segments(frag):
                    _nb.setdefault(_i, []).append(_j)
                    _nb.setdefault(_j, []).append(_i)
                marks = []
                for i in ids:
                    if i not in pts:
                        continue
                    marks += C.electron_symbols(
                        C.next_free_id(frag) + len(marks), pts[i][0], pts[i][1],
                        lone_pair_dirs(i, pts, _nb, n_e), i,
                        kind="lonepair")
                if marks:
                    k = frag.rindex("</fragment>")
                    frag = frag[:k] + "".join(marks) + frag[k:]

            elif style == "color":
                frag = B.color_region(frag, ids, C.color_index(colour))

            elif style == "fill":
                # NATIVE, and it goes INSIDE the fragment: a ColoredMolecularArea
                # is a property of the molecule, not a shape laid over it, so it
                # survives editing and moves with the structure. Same tag a ring
                # highlight drawn by hand in ChemDraw carries.
                bonds = B.region_bonds(frag, ids)
                # a bond with no id cannot be named; fragment() gives every bond
                # one, a hand-drawn library entry is whatever ChemDraw wrote
                bids = [b for b, _a, _c in bonds if b is not None]
                if not _encloses_a_ring(bonds) or not bids:
                    # PROVISIONAL, and the reason is honest ignorance: nobody has
                    # yet drawn a non-cyclic highlight by hand, so what ChemDraw
                    # renders for an area bounded by a chain is unknown. Tinting
                    # is the outcome that cannot look wrong -- and the caller is
                    # told, because asking for a fill and getting a tint is
                    # visible on paper and would otherwise be a mystery.
                    #
                    # Once per ANCHOR, not once per match: one query landing on
                    # three chains is one thing to fix, and "3 anchors" would be a
                    # false count of how many specs need attention.
                    _subj = S._subject(smiles, where)
                    if (_subj, on) not in S.FILL_NOT_RING:
                        S.FILL_NOT_RING.append((_subj, on))
                    frag = B.color_region(frag, ids, C.color_index(colour))
                    continue
                area = C.colored_area(nid, C.color_index(colour), bids)
                nid += 1
                # LAST </fragment> -- a superatom nests one, and inserting at the
                # first would put the highlight inside the abbreviation
                k = frag.rindex("</fragment>")
                frag = frag[:k] + area + frag[k:]

            elif style == "dot":
                # one disc per atom, centred ON it -- three of four hand-drawn
                # examples sit exactly on an atom, the fourth is the odd one out.
                # It goes BEHIND the atom: a marker that hides what it marks is no use.
                for x, y in here:
                    front.append(C.disc(nid, x, y,
                                        C.ATOM_DOT_R * C.BOND_LENGTH, colour))
                    nid += 1

            elif style == "halo":
                # An ELLIPSE fitted to the group's own proportions, not a circle
                # round its longest reach. A coumarin core is twice as wide as it
                # is tall; a circle big enough to hold it is half empty, and the
                # empty half makes it read as a highlighter blot rather than as
                # a glow behind a molecule.
                xs = [x for x, _y in here]
                ys = [y for _x, y in here]
                cx, cy = (min(xs) + max(xs)) / 2.0, (min(ys) + max(ys)) / 2.0
                # No separate floor for a one-atom group: the pad alone is 0.8 of
                # a bond, which is already a visible glow, and a max() against a
                # smaller floor could never fire.
                pad = HALO_PAD * C.BOND_LENGTH
                rx = (max(xs) - min(xs)) / 2.0 + pad
                ry = (max(ys) - min(ys)) / 2.0 + pad
                # furthest back of all: a glow under the dots as well as under
                # the structure, so two overlapping specs still read
                behind.append(C.disc(nid, cx, cy, rx, colour, shaded=True, ry=ry))
                nid += 1

            elif style == "locants":
                labels = spec.get("labels") or []
                cx, cy = _centre(here)
                for k, (x, y) in enumerate(here):
                    if k >= len(labels) or labels[k] in (None, ""):
                        continue
                    # Radially outward from the ring centre through the atom --
                    # which is how five of six hand-placed numbers sit; the sixth
                    # was pushed sideways by a substituent, and steering around
                    # neighbours is a judgement this does not try to make.
                    ux, uy = _unit(x - cx, y - cy)
                    px = x + ux * C.LOCANT_OFFSET
                    py = y + uy * C.LOCANT_OFFSET
                    w, h = C.M.box([str(labels[k])], C.LOCANT_SIZE)
                    front.append(C.styled_text(
                        nid, px, py - h / 2.0 + 0.716 * C.LOCANT_SIZE,
                        [str(labels[k])], size=C.LOCANT_SIZE,
                        color=None if colour == "black" else colour,
                        align="center"))
                    nid += 1

    # glows, then dots and numbers, then (the caller's job) the structure
    return "\n".join([p for p in ("\n".join(behind), "\n".join(front)) if p]), frag
