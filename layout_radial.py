"""
layout_radial.py — radial / divergent figure engine (render once).

Follows the earlier JavaScript implementation: one center structure with up to 8
products on a 3x3 grid, native angled arrows from the center (or from another
product for secondary steps) trimmed to sit just outside both boxes AND clear the
product labels, bold "id, yield" + superscript-footnote labels, and an optional
rounded box around the center (drawn with Lines, since the Rectangle graphic is
not exported to SVG).

Fragments are pre-built (frag, bbox) via smiles_to_fragment, and pre-colored via
beautify if desired - same convention as the linear engine.

Recipe:
  {
    "center": {"frag": (frag, bbox),
               # set UNDER the structure and centred, the way a product's id and
               # yield are; one string, or a list for the two- and three-line
               # form published figures use ("2a" / "70%" / "from 2.0 mmol of 1a").
               "label": "BF~4~^-^" | ["2a", "70%"]?,
               "box": True?},          # encloses the structure AND that label
    "products": {                 # any subset of the 8 grid slots
      #  TL  T  TR
      #  L   .  R
      #  BL  B  BR
      "TR": {"frag": (frag, bbox), "id": "9a", "yield": "94%", "foot": "a",
             "direction": "in" | "out"?,   # DEFAULT "out"; see below
             "above": ["Pd(OAc)~2~", "K~2~CO~3~"]?,   # annotates the arrow that
             "below": ["DMF, 80 degC"]?},             # POINTS AT this product
      ...
    },
    "arrowsFrom": {"R": "TR"}?,    # product R's arrow starts at product TR (secondary step)
    "footnotes": ["^a^ 0.2 mmol scale.", ...]?,   # under the whole figure
    "trunk": True?,   # same-side slots hang off one vertical line (see below)
  }

WHICH WAY EACH SPOKE POINTS
---------------------------
`direction` per slot, "out" (the default) or "in". "out" is the hub feeding a
product; "in" is a species feeding the hub -- the reaction's ENTRY, which outgoing
spokes alone cannot draw: the hub would have outgoing arrows and no incoming one.

PER SPOKE, NOT PER FIGURE. The figures that motivate it are "one in, several out",
so a figure-level switch covers neither of them.

IT IS AN ARGUMENT SWAP, NOT A FLIPPED ARROWHEAD. draw_arrow's two trims are by
ROLE -- src_margin at the `from` end, tgt_margin at the `to` end -- so an inbound
spoke is drawn by handing it (slot, hub), and every margin lands on the right end
by construction. The cheaper-looking alternative, keeping the call and swapping
Head3D / Tail3D at the emitter, would leave src_margin on the HUB and tgt_margin on
the slot: the margins for an arrow pointing the other way.

NOTHING ELSE IS ROLE-DEPENDENT, and that is measured rather than assumed. s_hh and
t_hh take their below-overhang from the sign of uy, which flips with the direction
and hands the two terms to each other; caption_box canonicalises both its normal and
its reading direction, so a caption block is invariant under the swap (checked on 16
combinations of slot x side x cap_rotate, identical to the last bit). The two
margins are the whole of the asymmetry: reversing one spoke of an eight-spoke figure
changes exactly TWO elements in the document -- that spoke's arrow and that spoke's
own caption -- each by 2.0 units toward the hub, which is src_margin minus
tgt_margin. The shaft's LENGTH does not change.

AN INBOUND SLOT IS NEVER TRUNKED, the same rule an arrowsFrom slot already has and
for the same reason: a trunk is the HUB's fan-out, and a spoke that does not leave
the hub has nothing to hang off it.

THE TRUNK, AND WHY BOTH SIDES ARE THE **FORK** FORM
---------------------------------------------------
With "trunk" on, the slots down one side (TL/L/BL, or TR/R/BR) stop radiating
separately: one stub leaves the hub, meets a vertical line, and each product
hangs off that line on an ordinary horizontal arrow. Drawn by rowspan.emit_bus,
which already does exactly this for a linear fork.

The obvious reading is that the RIGHT side is emit_bus's fork (spine to the right
of the cell) and the LEFT side its merge (spine to the left). The x is right and
the ARROWHEAD IS NOT. Measured on emit_bus directly: with fork=False the stub
comes back as

    <arrow Tail3D="-84 0" Head3D="-60 0">      head ON the hub

because a merge means several things converging INTO one, so emit_bus gives the
stub the head and leaves the branches headless. A radial figure diverges -- the
hub is the source -- so that arrowhead points backwards, and the branches are the
things that need heads.

So both sides use fork=True, and the left one is mirrored by giving it a NEGATIVE
lead: bx = x1 + lead lands hub-ward of the left edge, and the stub stays the plain
line the fork branch draws. The geometry is then a true mirror -- verified by
reading the emitted elements on both sides -- and every arrowhead sits on a
product, where a divergent figure puts it.

WHICH SIDE "above" MEANS ON A DIAGONAL
--------------------------------------
A caption goes on the shaft's normal, at `cap_gap` clear of it -- but the normal
of a direction is a PAIR, and picking the one by a fixed handedness (rotate u by
-90 degrees) makes "above" flip to the underside as soon as the arrow points
left. The L slot's arrow does point left, so that is not a corner case.

So the normal is oriented canonically: the one pointing into the upper half of
the page (ny < 0; ties, i.e. a vertical shaft, break to nx > 0). "above" is then
the upper side of every shaft, whichever way it runs, and for a rightward arrow
that is the same side linear puts it on -- the degenerate case agrees.

For a vertical spoke (T, B) there is no upper side, so above/below become right
and left. That is already this codebase's convention: linear's txt_gap_side
exists precisely because a vertical arrow's captions sit beside the shaft.

CAPTION TEXT ON A DIAGONAL MAY BE SET ALONG THE SHAFT
-----------------------------------------------------
`cap_rotate` (off by default), and DIAGONALS ONLY. Horizontal caption text is the
default, but published radial figures also set the conditions along the spoke.

WHY ONLY THE DIAGONALS, and it is the same half-plane rule twice over rather than a
restriction bolted on. Type has to read left-to-right, so the shaft direction is
canonicalised to ux > 0; the normal is canonicalised to ny < 0. For n = (uy, -ux)
those two conditions are THE SAME -- ny < 0 exactly when ux > 0 -- so on every
diagonal and horizontal slot "above" keeps meaning the upper side of the page and the
u/-u invariance is untouched. They can only disagree where ux == 0, which is exactly
T and B, and there rotating would move "above" from the right of the shaft to the
left and silently relocate every existing vertical caption. So T and B stay
horizontal, and a trunked side does too because the trunk has already turned its
diagonals into horizontal arrows.

WHAT IT COSTS, so nobody turns it on expecting it to be free. Rotating collapses the
caption's displacement along the normal from cap_gap + _support(w, h, ..) to
cap_gap + h/2, which moves the block toward the hub; all 44 measured row-pitch points
move, by up to 131.7 units, and the pitch has to be re-solved (it is -- see
diagonal_caption_row and test_radial_diagonal_caption). It also opens two
obstacle classes the horizontal geometry never had: a caption against ANOTHER
caption, and a caption against another spoke's shaft.

WHAT ChemDraw DOES WITH IT, measured on a hand-drawn sample rather than assumed:
RotationAngle is degrees x 65536, normalised to [0, 360) (a -45 rotation comes back
as 315, so cdxml_core reduces mod 360 before writing); the rotation is about `p`, the
first baseline's centre, which is where styled_text already points; and the box
ChemDraw then reports is the AABB of the rotated INK, which is SMALLER than the
rotated advance box this engine reserves. See caption_shape for why that direction is
the safe one.
"""

from __future__ import annotations

import math
import re

import cdxml_core as C

POS = {"TL": (-1, -1), "T": (0, -1), "TR": (1, -1), "L": (-1, 0), "R": (1, 0),
       "BL": (-1, 1), "B": (0, 1), "BR": (1, 1)}
ASCENT_EM = 0.716

# THE ARROWHEAD'S OWN LENGTH, read off the emitter rather than restated here. arrow_xy
# hard-codes HeadSize="1000" (cdxml_core.py:1032), and HeadSize is in HUNDREDTHS of a
# page unit -- the convention arrows.py records for ArrowShaftSpacing, where ChemDraw's
# own 300 was measured as 3.00 units. Two independent checks agree with that reading
# and with no other:
#
#   ArrowheadWidth = 0.250 * HeadSize = 250, i.e. 2.50 units of half-width, which is
#   exactly the max single-side reach measured for a plain head at HeadSize 1000 --
#   2.50, taken from ChemDraw's own recomputed BoundingBox.
#
#   ArrowheadCenterSize = 0.875 * HeadSize = 875, i.e. 8.75 units, putting the barb's
#   notch 1.25 units behind a 10.00-unit tip. That is the shallow concave back of the
#   ACS arrowhead; under any other scaling the pair is nonsense.
#
# WHAT IS NOT DIRECTLY MEASURED, and it is the soft spot: the ALONG-AXIS length cannot
# be recovered from any BoundingBox in this repo, because on that side the box is
# bounded by the SHAFT. What a hand-drawn ChemDraw export does pin (an arrow at
# HeadSize 1000, ArrowheadHead="Full") is that the head overhangs the
# Head3D point by +0.00 -- the glyph is drawn entirely BEHIND the head point. So this
# length is spent out of the shaft, which is what makes it the right thing to hold a
# shaft length against.
_HEAD_SIZE = float(re.search(r'HeadSize="(\d+)"',
                             C.arrow_xy(0, 0.0, 0.0, 1.0, 0.0)).group(1))
HEAD_LEN = _HEAD_SIZE / 100.0             # 10.00 units at the emitter's HeadSize 1000

# The shortest shaft this engine draws without complaining: the arrowhead itself, plus
# half a bond length of actual line behind it. C.BOND_LENGTH (14.4) is the engine's
# bond -- NOT allocator.BOND, which is 30.0 and is that module's tolerance scale.
MIN_ARROW_LEN = HEAD_LEN + 0.5 * C.BOND_LENGTH        # 10.00 + 7.20 = 17.20
# MINIMUM means the smallest ACCEPTABLE span, so a shaft exactly that long is fine and
# only a shorter one is reported. The tolerance is there to absorb IEEE noise, nothing
# more: auto_col_gap solves for span == MIN_ARROW_LEN exactly and lands 1e-14 under it,
# which would otherwise make the solver trip its own guard. Nothing observable is being
# waived -- arrow_xy writes coordinates at one decimal place, so a shortfall under 0.05
# units cannot reach the file at all, and this is eleven orders of magnitude below that.
_SPAN_TOL = 1e-6

CFG = dict(
    # THE FALLBACK, NOT WHAT THE TOOL USES. draw_radial_figure solves both gaps from
    # the figure (auto_col_gap below, auto_row_gap for a trunk) and passes numbers
    # down, so these two are what a DIRECT caller of render_radial gets -- the demo
    # scripts and this suite. They are deliberately left at fixed default values so
    # that a test which pins nothing still lays out the figure it was written against;
    # the solving lives one layer up, where the recipe's trunk flag and the caller's
    # own pinned values are both in hand.
    col_gap=300.0, row_gap=285.0,
    label_gap=14.0,
    num_size=10, foot_size=7,
    # The CENTRE's label sets one point larger than a product's, which is the
    # hierarchy published figures use. It is a named value rather than a bare 11 in
    # the put_label call because auto_col_gap has to measure the very label the
    # renderer draws, and two
    # spellings of one size is how a solved gap ends up disagreeing with the ink.
    cen_size=11,
    src_margin=24.0, tgt_margin=22.0,
    box_pad_x=34.0, box_pad_t=20.0, box_pad_b=18.0,
    margin=40.0,
    cap_size=9,          # arrow caption type. 9 in linear and 9 in cycle, so a
                         # reader moving between the three engines sees one size.
    cap_gap=4.0,         # caption box edge -> shaft, measured PERPENDICULAR to the
                         # shaft. This is linear's txt_gap_side, not its
                         # txt_gap_above/below pair: that pair is asymmetric because
                         # it was tuned for type sitting over and under a HORIZONTAL
                         # line, and linear's own comment says it is "the wrong axis
                         # once the arrow turns vertical". Every radial spoke is off
                         # the horizontal, so the symmetric side value is the one
                         # that applies, and it is already measured.
    footnote_size=8,     # the footnote block under the whole figure. The three
    footnote_gap=22.0,   # names and the three values are layout_grid.GRID_DEFAULTS'
    footnote_lead=3.0,   # and layout_linear.FOOTNOTE_DEFAULTS', so a footnote under a
                         # radial figure, under a linear scheme and under a scope
                         # table are ONE object at one size in one place.
    #
    # THE GAP IS 22.0, NOT table_layout's 10.0. There notes_gap is 2 x row_pad: the
    # distance from a data table's bottom RULE down to the notes tucked under it,
    # which is a caption hanging off a table. This block is not that. It is running
    # prose under a whole FIGURE, and the other two engines that set one put it 22.0
    # below -- so the numbers to share are theirs, not table_layout's.
    #
    # footnote_width: WHICH FIGURE the prose runs under, not how wide to make it.
    #
    # A footnote is running prose under THIS figure, so there is no second honest
    # value for a caller to choose -- but a radial can BE a part of a larger figure:
    # a composer draws it beside a cycle and the thing the note runs under is the
    # composed figure, not this call's bbox. Measured on a 608-wide composition,
    # folding to this call's own width gave 186.64 -- the right definition applied to
    # the wrong figure.
    #
    # So None means "the figure this call drew", and a number means "the figure this
    # call is part of". Not a second honest value; the same one, pointed at the right
    # figure. layout_grid and layout_linear carry this key for the same reason.
    footnote_width=None,
    # Widest a caption line may render before it is re-broken at word boundaries,
    # IN BOND LENGTHS -- the same 10.0 the cycle engine folds at, and the same
    # quantity, but pinned on RADIAL's evidence rather than inherited. The window
    # here is measured, two-sided, and wider than cycle's:
    #
    #   floor    6.18 bond (89.0 units), the widest caption line in the
    #            reference radial figures -- "morpholine, Pd~2~(dba)~3~". Fold below
    #            that and a figure that was already right changes. 10 clears it
    #            by 62%.
    #   ceiling  13.86 bond (199.6 units) is the widest single line that still
    #            draws clean on four diagonal spokes; at 15.88 bond (228.6) the
    #            TR and TL captions collide. So the cap must sit under 13.86,
    #            and 10 sits 28% under it.
    #
    # Cycle's own reasoning does NOT carry over wholesale: it argues the cap must
    # be a fixed constant because deriving it from available space closes a loop
    # through its radius bisection. Radial has no such loop -- the grid is fixed
    # at col_gap/row_gap whatever the captions do -- so the ceiling above could be
    # measured directly, which is why it is a measurement here and an argument
    # there. What does carry over is the typographic check: 144 units is about 32
    # characters at cap_size 9, the usual line length in published artwork, and
    # radial sets cap_size 9 too.
    #
    # It is radial's OWN key rather than a read of cycle's cfg on purpose: cycle's
    # floor and ceiling are properties of its ring, and retuning them there must
    # not silently move type in a radial figure.
    cap_fold_bond=10.0,
    # Set the conditions ALONG a diagonal shaft instead of horizontally beside it,
    # which is what published radial figures do. OFF by default and it has to be:
    # rotating a caption collapses its normal displacement from
    # cap_gap + _support(w, h, ..) to cap_gap + h/2, so the block moves toward the
    # hub and the solved row pitch moves on all 44 measured points, by up to 131.7.
    # There is no version of this that leaves an existing figure where it was, so
    # the only way an existing figure cannot move is for the feature not to run
    # unless it is asked for -- the same argument the trunk makes.
    #
    # DIAGONALS ONLY. A horizontal or vertical shaft is untouched even with this on:
    # for those the reading direction and the page normal disagree about which side
    # "above" is (they agree on all six diagonal and horizontal slots and conflict
    # on exactly T and B), and the published figures that motivate the feature set
    # type along diagonals. See the module docstring.
    cap_rotate=False,
)


from allocator import (ARROW_INK_WIDTH as _ARROW_INK_WIDTH,
                       Allocator as _Allocator, OBox as _OBox, Rect as _Rect,
                       Severity as _Severity, run_shapes as _run_shapes,
                       sat_depth as _sat_depth)
# Share the linear engine's analytical text model so all three engines set type
# identically -- and so radial labels are measured rather than guessed at.
from layout_linear import block_box as _block_box
# THE one folding routine, not a second one. fold_box hands back the folded lines
# and the box they actually occupy in a single call, so the rectangle reserved for
# a caption cannot disagree with the lines drawn into it -- and a line that already
# fits is passed through untouched, which is what keeps a hand-split caption
# byte-identical to the same caption drawn without folding.
from layout_linear import fold_box as _fold_box
# The clearance auto_row_gap adds over the collision floor. Read from linear rather
# than restated, because it is the same relation linear defined it for: "clearance
# between adjacent rows' ACTUAL content extents".
from layout_linear import CFG as _LL_CFG
LL_ROW_CLEAR = _LL_CFG["row_clear"]
# The fan's spine and its stub, from the engine that already draws them. Verified
# callable from outside layout_linear: `plans` and `shifts` are unreferenced in
# emit_bus (checked at the bytecode level -- no LOAD_FAST, no cellvars, no nested
# code objects), so None is passed for both. It reads a span's
# {owner, last, x0, x1, cy, lead, fork} and a row_base indexed by row number.
from rowspan import emit_bus as _emit_bus, fork_stub_y as _fork_stub_y

# HOW THICK THE TRUNK'S LINE IS, measured off emit_bus rather than restated. It
# reserves its spine as (bx - 1, ...) to (bx + 1, ...), so the line reaches 1.0 either
# side -- but that 1 lives in rowspan and a copy here would be free to drift from it.
# One test span with a known geometry, and the box that comes back IS the answer.
def _trunk_ink_half() -> float:
    _els, _boxes, _g = _emit_bus(
        [dict(owner=0, last=1, x0=0.0, x1=0.0, cy=0.0, lead=1.0, fork=True)],
        None, None, [-10.0, 10.0], {"gap": 1.0}, 0, ink=True)
    spine = max(_boxes, key=lambda b: b[3] - b[1])      # the tall one
    return (spine[2] - spine[0]) / 2.0


TRUNK_INK_HALF = _trunk_ink_half()

# Which grid slots can share one vertical trunk. T and B sit at x = 0, on no side
# at all, so they can never hang on a vertical line to either hand and keep their
# own spoke whatever the trunk setting is.
TRUNK_SIDES = {"left": ("TL", "L", "BL"), "right": ("TR", "R", "BR")}
# ONE definition of the support formula, not a second copy of it. `_support` is
# where "how far must this box move to clear that line" was worked out (for the
# cycle's ring); a radial spoke asks the identical question of a straight shaft.
from layout_cycle import _support


def _reach(hw, hh, ux, uy):
    return min(abs(hw / ux) if ux != 0 else 1e9, abs(hh / uy) if uy != 0 else 1e9)


def label_text(d: dict) -> str:
    """The compound label under one product: bold "id, yield" plus a footnote mark.

    A function rather than an inline f-string because auto_row_gap has to measure
    the very label the renderer draws. Two spellings of one string is how a
    reservation ends up disagreeing with the ink inside it.
    """
    raw = d.get("id")
    num = "" if raw is None else str(raw).strip()
    yld = d.get("yield", "")
    foot = d.get("foot", "")
    # No id, no bold: #...# around an empty id left a bold ", 88%" with
    # nothing before the comma. Only the yield is left to write.
    if not num:
        return (yld or "") + (f"^[{foot}]^" if foot else "")
    return (f"#{num}{', ' + yld if yld else ''}#"
            + (f"^[{foot}]^" if foot else ""))


def auto_row_gap(products: dict, cfg: dict | None = None) -> float:
    """Row pitch derived from the products, for a TRUNK figure. Opt-in.

    WHAT SETS THE FLOOR, measured rather than argued. Sweeping row_gap down in
    2.5-unit steps over five product profiles (h from 12.5 to 72.9) and both spoke
    styles, the first collision is always the same pair -- a product's own compound
    LABEL against the STRUCTURE in the row below it -- and it arrives at

        row_gap  =  product height + label_gap + label height

    predicted against measured, in units:  87.5/87.5, 94.6/95.0, 95.6/97.5,
    51.5/52.5. Every residual is positive and inside one sweep step, so the
    expression is the floor and never over-states it. It is the same in both spoke
    styles, because a label hangs below its product whichever way the arrow came in.

    THE MARGIN IS layout_linear's row_clear (28.0), not a new number: that constant
    is defined as "clearance between adjacent rows' ACTUAL content extents", which
    is exactly the relation here. So the answer is the floor plus one row_clear.

    WHY TRUNK ONLY. With horizontal arrows the caption lives in the corridor
    between the trunk and the product, so it is not in the vertical gap at all and
    the floor above is independent of how many caption lines there are -- measured
    52.5/95.0/87.5/97.5 for 0, 1 and 3 lines alike. With DIAGONAL spokes the
    caption is pushed along the shaft normal INTO that gap, and the floor then
    moves in a way none of these quantities predicts: at three caption lines it
    went to 87.5 for a product whose formula floor is 35.2, and to 150.0 for one
    whose floor is 87.5. Two profiles, opposite directions, so there is no margin
    that fixes it -- which is why this refuses to guess for diagonal figures rather
    than returning a number that is wrong for some of them.
    """
    cfg = {**CFG, **(cfg or {})}
    h_max = max((e["frag"][1][3] - e["frag"][1][1]) for e in products.values())
    lab_h = max(_block_box([label_text(e)], cfg["num_size"])[1]
                for e in products.values())
    return h_max + cfg["label_gap"] + lab_h + LL_ROW_CLEAR


# WHICH WAY ONE SPOKE POINTS, and THE only place the string is read. Three callers
# ask -- the renderer, trunked_slots and _spoke_ends -- and three spellings of a
# default is how one of them draws what the other two did not solve for.
#
# RAISES ON ANYTHING ELSE rather than falling back to "out". A misspelt
# "direction": "input" that quietly meant "out" draws the arrow the wrong way round
# with nothing on the page to say so -- the same silent-wrong-chemistry failure the
# MIN_ARROW_LEN guard exists to catch, and the same reason render_radial's cfg gate
# raises rather than coercing.
def _inbound(entry: dict) -> bool:
    d = entry.get("direction", "out")
    if d not in ("in", "out"):
        raise ValueError(
            'direction must be "in" or "out", got %r. "out" (the default) points '
            'from the hub AT this slot; "in" points from this slot AT the hub.'
            % (d,))
    return d == "in"


def trunked_slots(products: dict, trunk: bool,
                  arrows_from: dict | None = None) -> dict:
    """{side: [slots]} for the sides that actually get a trunk. ONE definition.

    render_radial decides this to draw the bus; auto_col_gap has to ask the same
    question, because a trunked slot's arrow is HORIZONTAL and a horizontal arrow is
    the one whose caption is squeezed between the hub and the product. Two copies of
    "is this slot trunked" is how a solved gap ends up sizing a spoke the renderer
    does not draw.

    TWO OR MORE PER SIDE, never a slot fed by another product, and never an INBOUND
    slot: all three rules are the renderer's, and they are here rather than there so
    there is one of each. An inbound spoke does not leave the hub, and a trunk is the
    hub's fan-out -- there is nothing for it to hang off.
    """
    af = arrows_from or {}
    out: dict = {}
    if not trunk:
        return out
    for side, slots in TRUNK_SIDES.items():
        members = [t for t in slots if t in products and t not in af
                   and not _inbound(products[t])]
        if len(members) >= 2:
            out[side] = sorted(members, key=lambda t: POS[t][1])
    return out


def caption_width(entry: dict, cfg: dict | None = None) -> float:
    """The widest caption block on one product's arrow, FOLDED as the renderer folds.

    Both sides, because "above" and "below" are each centred on the shaft midpoint
    independently, so the wider of the two is what has to fit. fold_box is the same
    call draw_arrow makes -- a second spelling of the fold is a reservation that
    disagrees with the ink in it.
    """
    cfg = {**CFG, **(cfg or {})}
    w = 0.0
    for side in ("above", "below"):
        lines = entry.get(side)
        if not lines:
            continue
        lines = [lines] if isinstance(lines, str) else list(lines)
        _l, (w_, _h) = _fold_box(lines, cfg["cap_size"],
                                 cfg["cap_fold_bond"] * C.BOND_LENGTH)
        w = max(w, w_)
    return w


def caption_box(tail: tuple, head: tuple, lines, side: str,
                cfg: dict | None = None) -> tuple:
    """(folded lines, (cx, cy, w, h), rot) for one caption block on one shaft.

    THE renderer's own placement, not a second spelling of it. draw_arrow calls this
    and so does diagonal_caption_row, because a row pitch solved against a block
    somewhere other than where the block is drawn does not solve anything -- and this
    placement is subtle enough that two copies would drift. It takes the shaft's two
    TRIMMED ends, which is what the caption is actually centred on: the untrimmed
    hub-to-product line is not where the arrow is.

    (cx, cy) is the block's CENTRE in both modes, and `rot` is its angle in degrees,
    0.0 when the block is set horizontally. The caller turns those into either a Rect
    or an OBox; nothing else needs to know which mode it is in.

    THE NORMAL IS ORIENTED, not taken by a fixed handedness -- see the module
    docstring. "above" is the upper of the two normals, so it means the same side of
    the page whichever way round the spoke was written.

    ROTATION (cap_rotate, off by default, diagonals only) changes two terms and
    nothing else:

      THE READING DIRECTION. Type has to run left-to-right, so u is canonicalised to
      ux > 0 -- and that is the SAME half-plane rule already applied to the normal,
      because for n = (uy, -ux) the condition ny < 0 holds exactly when ux > 0. So
      "above" keeps meaning the upper side of the page on every diagonal, the u/-u
      invariance survives untouched, and a left-hand spoke's type reads back toward
      the hub, which is what published figures do. The two rules disagree only where
      ux == 0, i.e. the T and B spokes, and those are excluded.

      THE DISPLACEMENT. A shaft-ALIGNED box's support along the normal is exactly
      h/2, so the _support term -- which is the half-extent of an AXIS-ALIGNED box
      along a slanted direction, and the reason a wide diagonal caption is pushed so
      far out -- is replaced rather than adjusted. It cannot be replaced inside
      _support itself: layout_cycle imports that for its ring feeders and must not
      move.
    """
    cfg = {**CFG, **(cfg or {})}
    x0, y0 = tail
    x1, y1 = head
    dx, dy = x1 - x0, y1 - y0
    L = math.hypot(dx, dy) or 1.0
    ux, uy = dx / L, dy / L
    nx, ny = (uy, -ux)
    if ny > 0 or (ny == 0 and nx < 0):
        nx, ny = -nx, -ny
    lines = [lines] if isinstance(lines, str) else list(lines)
    # FOLD FIRST, then measure -- and take both from one call, so the box below is
    # the box of the lines that get drawn. Each of the caller's lines is folded on
    # its own, so a deliberate break survives and only an over-wide line is broken
    # further.
    lines, (w_, h_) = _fold_box(lines, cfg["cap_size"],
                                cfg["cap_fold_bond"] * C.BOND_LENGTH)
    sgn = 1.0 if side == "above" else -1.0
    sx, sy = nx * sgn, ny * sgn
    # DIAGONAL means neither component is zero, which is the geometric question and
    # not the slot's name. A trunked spoke has been turned into a horizontal arrow by
    # the time it gets here, so it is correctly left unrotated with nothing extra
    # plumbed through to say so.
    rot = 0.0
    if cfg["cap_rotate"] and ux != 0.0 and uy != 0.0:
        rux, ruy = (ux, uy) if ux > 0 else (-ux, -uy)
        rot = math.degrees(math.atan2(ruy, rux))
        d_ = cfg["cap_gap"] + h_ / 2.0
    else:
        # _support is how far the box's own corner reaches this way, so the block
        # clears the SHAFT LINE by exactly cap_gap and no guessing is involved. On a
        # diagonal that pushes a wide block a long way out -- that is the true
        # geometry of an axis-aligned box beside a slanted line, not slack to be
        # tuned away. It is also the term that makes a diagonal caption claim ROW
        # space: sy is non-zero there, and _support carries the block's HEIGHT into
        # the displacement.
        d_ = cfg["cap_gap"] + _support(w_, h_, sx, sy)
    return (lines,
            ((x0 + x1) / 2.0 + sx * d_, (y0 + y1) / 2.0 + sy * d_, w_, h_),
            rot)


def caption_shape(cx: float, cy: float, w: float, h: float, rot: float,
                  size: float) -> tuple:
    """(reservation Rect, OBox or None, styled_text's `p`) for one caption block.

    ONE PLACE, because the renderer reserves the block and the solver tests it and a
    third spelling of "where is this block" is how the two come apart. rot == 0.0
    returns shape None and exactly the Rect and `p` that shipped before rotation
    existed -- so an unrotated figure is byte-identical by construction rather than
    by a guard.

    WHERE `p` GOES. styled_text puts the FIRST BASELINE at `p`, ASCENT_EM * size
    below the block's top, so `p` sits e = (0, ascent - h/2) from the block's centre.
    ChemDraw rotates the block about `p` -- that is the sample's finding, not an
    assumption -- so the centre lands at p + Rot(-e) and `p` therefore has to be
    centre + Rot(e) for the block to end up centred on (cx, cy).

    WHY THE ADVANCE BOX IS SAFE TO RESERVE, measured on the hand-drawn sample. The
    rotated BoundingBox ChemDraw writes is the AABB of the rotated INK, and ink is
    smaller than the advance box by two terms: RAGGEDNESS (narrower lines do not fill
    the block's width -- 9.55 units of w+h on the real 3-line NaBH4 caption, 0 on a
    single line) and BEARINGS (1.3-1.6 units, ink narrower than the advance width).
    Total 11.11 on that caption, about 9%.

    So this OVER-reserves, and the containment is EXACT rather than approximate: in
    the unrotated frame the ink is a subset of the advance box, rotation about `p` is
    rigid, and a subset stays a subset under a rigid transform. Nothing has to be
    measured to know the ink fits inside what is reserved -- which is why the 9% is
    spent on purpose and not chased. Modelling the raggedness is a documented
    follow-up and cheap when wanted (fold_box already hands back the folded lines, so
    the per-line widths are free): one decision later, not a re-investigation.
    """
    ex, ey = 0.0, ASCENT_EM * size - h / 2.0
    if not rot:
        return (_Rect(cx - w / 2.0, cy - h / 2.0, w, h), None,
                (cx + ex, cy + ey))
    th = math.radians(rot)
    c, s = math.cos(th), math.sin(th)
    shape = _OBox(cx, cy, c, s, w / 2.0, h / 2.0)
    return shape.aabb(), shape, (cx + ex * c - ey * s, cy + ex * s + ey * c)


def center_half_width(center: dict, cfg: dict | None = None,
                      trunk: bool = False) -> float:
    """How far the hub's ink reaches sideways -- the same number draw_arrow trims to.

    PUBLIC, and shared with auto_col_gap rather than restated there, because a solved
    gap that measures a different hub from the one the arrow trims to is a gap that
    does not solve anything.

    Three cases, and they are the renderer's own, above:
      * boxed   -- the box encloses the wider of structure and label, plus box_pad_x.
      * bare    -- the spoke trims to the STRUCTURE and a wider label is left to
                   overhang, which is this engine's convention for products too.
      * trunked -- a trunk is a LINE DRAWN DOWN THE PAGE, so it clears the hub's whole
                   FOOTPRINT including a label wider than the structure (trunk_hw).
    """
    cfg = {**CFG, **(cfg or {})}
    bb = center["frag"][1]
    cw = bb[2] - bb[0]
    lab = center.get("label")
    lab_w = 0.0
    if lab:
        lines = [lab] if isinstance(lab, str) else [str(t) for t in lab]
        lab_w = _block_box(lines, cfg["cen_size"])[0]
    half = (max(cw / 2, lab_w / 2) + cfg["box_pad_x"]) if center.get("box") else cw / 2
    return max(half, lab_w / 2.0) if trunk else half


def center_half_height(center: dict, cfg: dict | None = None) -> tuple:
    """(how far the hub reaches ABOVE its origin, how much hangs BELOW that).

    The vertical twin of center_half_width, and it mirrors the pair render_radial
    builds as `center_half[1]` / `center_below`. It is a second spelling of those
    lines and there is a test that holds the two together against the renderer's own
    reservations -- the alternative was refactoring the renderer to call this, but its
    version is computed mid-layout from cen_lab_bot and the refactor would risk moving
    a figure to buy tidiness.
    """
    cfg = {**CFG, **(cfg or {})}
    bb = center["frag"][1]
    ch = bb[3] - bb[1]
    lab = center.get("label")
    lab_h = 0.0
    if lab:
        lines = [lab] if isinstance(lab, str) else [str(t) for t in lab]
        lab_h = _block_box(lines, cfg["cen_size"])[1]
    if center.get("box"):
        y0 = -ch / 2 - cfg["box_pad_t"]
        y1 = (ch / 2 + cfg["label_gap"] + lab_h if lab_h else ch / 2) + cfg["box_pad_b"]
        return -y0, y1 - (-y0)
    return ch / 2, (cfg["label_gap"] + lab_h if lab_h else 0.0)


def spoke_axis(slot: str, arrows_from: dict | None = None,
               on_trunk: frozenset | set = frozenset()) -> str:
    """"horizontal" | "vertical" | "diagonal" -- which way this product's arrow runs.

    Not the same question as where the slot sits. A trunked slot is fed by an ORDINARY
    HORIZONTAL arrow off the trunk however far up the page it is, and a slot listed in
    arrowsFrom is fed from another product rather than from the hub.
    """
    if slot in on_trunk:
        return "horizontal"
    af = arrows_from or {}
    if slot in af and af[slot] in POS:
        dx = POS[slot][0] - POS[af[slot]][0]
        dy = POS[slot][1] - POS[af[slot]][1]
    else:
        dx, dy = POS[slot]
    if dy == 0:
        return "horizontal"
    return "vertical" if dx == 0 else "diagonal"


def vertical_spoke_row(center: dict, products: dict, cfg: dict | None = None,
                       trunk: bool = False, arrows_from: dict | None = None) -> float:
    """Row pitch the VERTICAL spokes need. The col formula on the other axis -- almost.

    Same constraint: the arrow's two trims are measured inward from opposite ends and
    below some pitch they cross. But row_gap is NOT symmetric between up and down,
    because draw_arrow is not:

        s_hh = from_half[1] + (from_below if uy > 0 else 0.0)
        t_hh = to_half[1]   + (to_below   if uy < 0 else 0.0)

    So an UPWARD spoke pays the product's own compound label hanging below it, and a
    DOWNWARD spoke pays the hub's overhang instead. Writing those out:

        up    row >= c_half            + src + tgt + h/2 + label_gap + lab_h + min_arrow
        down  row >= c_half + c_below  + src + tgt + h/2                     + min_arrow

    A NAIVE MIRROR OF THE COLUMN FORMULA IS WRONG, and wrong by a constant: it misses
    22.70 going up (label_gap 14.0 + a measured 8.70 label) and over-states by 2.00
    going down (box_pad_b 18 is 2 less than box_pad_t 20, so a boxed hub's c_below is
    NEGATIVE). Measured across six profiles and both directions: the naive form is off
    by exactly -22.70 / +2.00 every time, and the direction-aware pair above is exact
    to 0.0000 on all twelve.

    Only PURE VERTICAL spokes are counted. A diagonal one is covered by col_gap (the
    diagonal span is the horizontal span divided by |ux| <= 1, so it is never the
    tighter of the two -- measured: at the solved col_gap a diagonal spoke never
    reports arrow-too-short at ANY row_gap down to 20). A horizontal one does not
    involve row_gap at all.
    """
    cfg = {**CFG, **(cfg or {})}
    on_trunk = {t for m in trunked_slots(products, trunk, arrows_from).values()
                for t in m}
    c_half, c_below = center_half_height(center, cfg)
    need = 0.0
    for slot, e in products.items():
        if spoke_axis(slot, arrows_from, on_trunk) != "vertical":
            continue
        src = (arrows_from or {}).get(slot)
        if src and src in products:
            sb = products[src]["frag"][1]
            s_half = (sb[3] - sb[1]) / 2.0
            s_below = cfg["label_gap"] + _block_box(
                [label_text(products[src])], cfg["num_size"])[1]
        else:
            s_half, s_below = c_half, c_below
        bb = e["frag"][1]
        lab_h = _block_box([label_text(e)], cfg["num_size"])[1]
        up = POS[slot][1] < (POS[src][1] if src and src in POS else 0)
        need = max(need,
                   s_half + (0.0 if up else s_below)
                   + cfg["src_margin"] + cfg["tgt_margin"]
                   + (bb[3] - bb[1]) / 2.0
                   + (cfg["label_gap"] + lab_h if up else 0.0)
                   + MIN_ARROW_LEN)
    return need


# HOW FINE THE LADDER IS that diagonal_caption_row walks before it bisects, in page
# units. It is a step size and not a tolerance: the answer is bisected to machine
# precision once the last dirty rung is found, and this only has to be small enough
# not to STEP OVER a clash band entirely.
#
# PINNED AGAINST THE NARROWEST BAND, MEASURED, rather than derived: the obvious
# derivation does not hold, so the arithmetic is set out here rather than left to
# contradict the measurement.
#
# THE OBVIOUS DERIVATION. A band exists while the caption overlaps an obstacle
# in x and closes when the block clears that obstacle's edge in y; the thinnest
# obstacle on a radial page is a one-line compound label, 8.70 units tall (measured,
# "9a, 94%"), and the caption's vertical position moves at roughly half the row pitch
# -- so about 17 units, which would make 1.0 seventeen times finer than needed.
#
# WHY IT DOES NOT BOUND ANYTHING. That argument assumes the obstacle is a compound
# label. With cap_rotate on, a caption also has to clear ANOTHER CAPTION and another
# spoke's SHAFT (see _caption_clash), and a shaft is ARROW_INK_WIDTH = 4.00 thick,
# less than half the 8.70 the bound was taken from. The derivation is sound for
# compound labels and says nothing about these.
#
# WHAT IS MEASURED. All 44 points swept at 1.0 to find the profiles with more
# than one band, then those refined at 0.05:
#
#   unrotated   2 of 44 have a second band; narrowest is 26.30 units
#               (3-slot benzene n=1, 59.75-86.05)
#   rotated     3 of 44 have a second band; narrowest is  7.40 units
#               (4-corner tolyl-ester n=2..4, 112.35-119.75)
#
# So the margin is 7.4x, not 17x. 1.0 still lands seven rungs inside the narrowest
# band, which is what this constant has to do, and the figure to beat if anyone
# widens the step is 7.40 -- not 17, and not the 26.30 the unrotated geometry offers.
#
# CLEAN GAPS NEED NO MARGIN and are deliberately not part of this bound. The narrowest
# measured is 6.15 units rotated (the same figure, 106.20-112.35) and 10.55 unrotated;
# stepping over a gap merges two bands and returns the top of the upper one, which is
# the conservative answer. Only stepping over a DIRTY band loses, and that is the
# direction bounded above.
#
# Also measured directly: at a 4.0 step the unrotated naphthalene 65-94 band is
# stepped clean over and the answer comes back 39.85 units low; at 1.0 every band in
# the six profiles is caught, in both modes.
_ROW_LADDER = 1.0
# The highest row pitch worth walking to. A radial figure is drawn on paper; 1200
# units is over three letter pages tall, so a figure needing more than this has a
# problem no row pitch is going to fix. Reaching the top means "no clean row exists",
# which is reported as such rather than silently answered with the top rung.
_ROW_CEILING = 1200.0


def diagonal_caption_row(center: dict, products: dict, col_gap: float,
                         cfg: dict | None = None, trunk: bool = False,
                         arrows_from: dict | None = None) -> float:
    """Row pitch a DIAGONAL spoke's caption needs. The term the plain formula misses.

    WHAT IT IS. On a diagonal the caption block is displaced along the shaft NORMAL,
    and that normal has a vertical component -- for a falling diagonal it is col/L --
    so the block is carried into the very gap being solved. The displacement is
    cap_gap + _support(w, h, ...), and _support carries the block's HEIGHT, which is
    why a bare or one-line caption does not break auto_row_gap's formula and a
    three-line one does. On a horizontal or vertical spoke the normal is axis-aligned
    and this term is empty, which is why nothing but the diagonals is walked.

    WHAT IT CLEARS, and these are the only two obstacles measured to bind:
      * the PRODUCT (or its compound label) in the row the block is pushed into --
        a BL spoke's "above" caption against the L slot, and its mirrors.
      * the HUB's own footprint, when no product sits between: its box if it is
        boxed, otherwise its structure and the label under it. This is what binds on
        a four-corner figure, which has nothing in the middle row at all.

    NOT A CLOSED FORM, and this is the honest description rather than a shortfall.
    row_gap appears inside L = hypot(col, row), inside both min() branches of _reach,
    and inside _support. Substituting the branch that binds gives

        (row + K)(col^2 + row^2) + L * (linear in row) - col * row * w - col^2 * h = 0

    and clearing the radical makes it a sextic. So the CLEARANCE is closed form -- it
    is written out below, exactly as draw_arrow places the block -- and only the
    inversion is numeric. Nothing is rendered and the allocator is not consulted.

    THE CLASH SET IS A UNION OF BANDS, NOT A HALF-LINE. As the pitch grows the spoke
    steepens and the caption SWEEPS ACROSS an obstacle, so a figure can be clean,
    then dirty, then clean again -- measured on naphthalene at a solved col_gap:
    clean over 54-64, dirty again over 65-93, clean from 94. The answer is therefore
    the top of the LAST band, walked down from _ROW_CEILING, and not the first clean
    pitch found from below. Taking the first clean one instead reports 53.27 where
    the true answer is 93.12.

    SOLVED TO CONTACT, not to contact plus a clearance. auto_row_gap adds one
    row_clear over its floor and it was natural to do the same here; measured over
    280 configurations, contact plus _SPAN_TOL leaves two figures worse than a
    flat 285 and contact plus row_clear leaves four -- the extra 28 units push a
    figure INTO the next band up. So the measurement chose contact, and the
    _SPAN_TOL is there for the same reason auto_col_gap carries one: the allocator
    holds a caption against a label at ZERO tolerance, so landing exactly on the
    boundary lets a float residue read as a collision.

    Returns 0.0 when no diagonal spoke carries a caption, and inf when no pitch under
    _ROW_CEILING is clean -- see row_gap_is_sound for the family that does that.
    """
    cfg = {**CFG, **(cfg or {})}
    on_trunk = {t for m in trunked_slots(products, trunk, arrows_from).values()
                for t in m}
    live = [(s, e) for s, e in products.items()
            if spoke_axis(s, arrows_from, on_trunk) == "diagonal"
            and (e.get("above") or e.get("below"))]
    if not live:
        return 0.0
    rungs = int((_ROW_CEILING - _ROW_LADDER) / _ROW_LADDER)
    ladder = [_ROW_LADDER * (i + 1) for i in range(rungs)]
    dirty = [r for r in ladder if _caption_clash(center, products, live, col_gap, r, cfg)]
    if not dirty:
        return 0.0
    last = max(dirty)
    if last >= ladder[-1]:
        return float("inf")
    lo, hi = last, last + _ROW_LADDER
    # 60 halvings takes a 1.0-unit bracket to 1e-18, which is below the float noise
    # in the clearance itself; the loop is bounded rather than tolerance-driven so it
    # cannot fail to terminate on a clearance that is flat in row.
    for _ in range(60):
        mid = (lo + hi) / 2.0
        if _caption_clash(center, products, live, col_gap, mid, cfg):
            lo = mid
        else:
            hi = mid
    return hi + _SPAN_TOL


def _caption_clash(center: dict, products: dict, live: list, col_gap: float,
                   row_gap: float, cfg: dict) -> bool:
    """Does any diagonal caption block land on something, at this row pitch?

    The clearance diagonal_caption_row inverts, written once. The overlap test is
    allocator.sat_depth -- the SAME function the allocator uses when the renderer
    reserves this block, so the solver and the report cannot disagree about what
    counts as touching. For two axis-aligned boxes it is min(dx, dy), i.e. exactly
    the per-axis test; for a ROTATED caption it is
    the oriented answer, which is why this needed no separate rotated branch.

    An overlap on both axes is also the reason a band closes: the x overlap is what
    switches on and off as the spoke steepens.

    TWO MORE OBSTACLE CLASSES ONCE THE BLOCK ROTATES, and they are not a refinement
    -- they are the difference between a solved pitch that is clean and one that
    ships a collision. _row_obstacles holds the hub and the products, and its
    docstring's claim that those are "the only two obstacles measured to bind" was
    measured on HORIZONTAL blocks. A rotated block swings along its own shaft, and
    then:

      ANOTHER CAPTION. Two adjacent spokes' captions rotate toward each other and
      meet. Measured: 3-slot tolyl-ester with four caption lines reports
      'center->L.above' on 'center->TL.above', 0.5 units deep over a 72-unit
      overlap, at the pitch the obstacle-only solver returned.

      ANOTHER SPOKE'S SHAFT. The block reaches across a neighbouring arrow.
      Measured on the same figure: 'center->BL.above' on 'center->L[0]', 2.8 deep.

    Both are invisible to a solver that only knows the hub and the products, and
    without them the answer is 5.35 units under on that figure -- the one point of
    the 44 that does not otherwise come back 0.00000.

    ITS OWN SHAFT IS EXEMPT, and by the same rule the allocator uses rather than a
    second spelling of it: related() treats 'arrow:center->BL.above' and
    'arrow:center->BL[0]' as one object because they share component 0, so a caption
    sitting on the arrow it annotates is not a finding there and must not be one
    here. Compared by SLOT, which is what that component reduces to.
    """
    obstacles = _row_obstacles(center, products, col_gap, row_gap, cfg)
    on_trunk: set = set()
    # Every caption that exists at this pitch, and every shaft. Built once for the
    # whole test rather than per block, because caption-vs-caption needs both of
    # them placed and a block must not be compared against itself.
    caps: list = []
    for slot, e in live:
        t, hd = _spoke_ends(center, e, slot, col_gap, row_gap, cfg)
        for side in ("above", "below"):
            if not e.get(side):
                continue
            _l, (cx, cy, w_, h_), _rot = caption_box(t, hd, e[side], side, cfg)
            _bb, _shape, _p = caption_shape(cx, cy, w_, h_, _rot, cfg["cap_size"])
            caps.append((slot, _shape if _shape is not None else _bb))
    # ASKED FOR, not restated. allocator.run_shapes is the very function
    # reserve_polyline reserves through, so the shaft this clears is the shaft that
    # is reserved -- including the half-width it reaches PAST each endpoint, which an
    # OBox of the bare segment does not have. Restating it as a bare strip was
    # measured to under-answer by 1.95 units on 3-slot tolyl-ester with four caption
    # lines, because the block caught exactly that end cap.
    shafts: list = []
    # built only when it will be used; see the gate below
    for slot, e in products.items():
        t, hd = _spoke_ends(center, e, slot, col_gap, row_gap, cfg)
        for _r, _s in _run_shapes([t, hd], _ARROW_INK_WIDTH):
            shafts.append((slot, _s if _s is not None else _r))
    # GATED ON ROTATION, and this is not caution -- it is what the unrotated
    # measurement says. The 44-point table came back 0.00000 on every point with the
    # hub and the products alone, so unrotated those really are the only two classes
    # that bind, exactly as _row_obstacles claims. Adding the two rotated ones there
    # would clear obstacles that do not bind and move a shipped figure for nothing;
    # measured, it moves 4-corner naphthalene with four caption lines. So an
    # unrotated figure keeps the solver it was measured with.
    rotated = bool(cfg["cap_rotate"])
    for i, (slot, mine) in enumerate(caps):
        for bx0, by0, bx1, by1 in obstacles:
            if _sat_depth(mine, _Rect(bx0, by0, bx1 - bx0, by1 - by0)) > 0.0:
                return True
        if not rotated:
            continue
        for j, (other_slot, other) in enumerate(caps):
            if j <= i:                      # each pair once; never against itself
                continue
            if _sat_depth(mine, other) > 0.0:
                return True
        for other_slot, shaft in shafts:
            if other_slot == slot:          # its own arrow -- related(), not a clash
                continue
            if _sat_depth(mine, shaft) > 0.0:
                return True
    return False


def _row_obstacles(center: dict, products: dict, col_gap: float, row_gap: float,
                   cfg: dict) -> list:
    """Everything a displaced caption can land on: the hub's footprint, and every
    product with the compound label under it.

    A caption IS checked against its own product -- allocator.related() compares the
    head of the owner path, and "arrow:center->BL" is not "prod:BL", so the renderer
    checks that pair too and this has to agree with it.
    """
    out = []
    bb = center["frag"][1]
    cw, ch = bb[2] - bb[0], bb[3] - bb[1]
    lab = center.get("label")
    cen_lab_w = cen_lab_h = 0.0
    if lab:
        lines = [lab] if isinstance(lab, str) else [str(t) for t in lab]
        cen_lab_w, cen_lab_h = _block_box(lines, cfg["cen_size"])
    cen_lab_bot = (ch / 2 + cfg["label_gap"] + cen_lab_h) if cen_lab_h else ch / 2
    if center.get("box"):
        # The BOX is the hub's footprint once there is one -- it encloses the
        # structure and the label both, so the two inside it are not obstacles in
        # their own right.
        hw = max(cw / 2, cen_lab_w / 2)
        out.append((-hw - cfg["box_pad_x"], -ch / 2 - cfg["box_pad_t"],
                    hw + cfg["box_pad_x"], cen_lab_bot + cfg["box_pad_b"]))
    else:
        out.append((-cw / 2, -ch / 2, cw / 2, ch / 2))
        if cen_lab_h:
            out.append((-cen_lab_w / 2, ch / 2 + cfg["label_gap"],
                        cen_lab_w / 2, cen_lab_bot))
    for slot, e in products.items():
        b = e["frag"][1]
        w, h = b[2] - b[0], b[3] - b[1]
        cx, cy = POS[slot][0] * col_gap, POS[slot][1] * row_gap
        out.append((cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2))
        lw, lh = _block_box([label_text(e)], cfg["num_size"])
        top = cy + h / 2 + cfg["label_gap"]
        out.append((cx - lw / 2, top, cx + lw / 2, top + lh))
    return out


def _spoke_ends(center: dict, entry: dict, slot: str, col_gap: float,
                row_gap: float, cfg: dict) -> tuple:
    """Where one hub spoke's shaft actually starts and stops, both trims applied.

    A SECOND SPELLING of draw_arrow's trim, the same way vertical_spoke_row carries
    one, and for the same reason: the alternative is threading the renderer's
    mid-layout state out to the solver. test_radial_diagonal_caption holds
    it against the renderer's own reservations so the two cannot drift.
    """
    c_half, c_below = center_half_height(center, cfg)
    c_hw = center_half_width(center, cfg)
    b = entry["frag"][1]
    w, h = b[2] - b[0], b[3] - b[1]
    lab_h = _block_box([label_text(entry)], cfg["num_size"])[1]
    dx, dy = POS[slot][0] * col_gap, POS[slot][1] * row_gap
    L = math.hypot(dx, dy) or 1.0
    ux, uy = dx / L, dy / L
    s_hh = c_half + (c_below if uy > 0 else 0.0)
    t_hh = h / 2 + ((cfg["label_gap"] + lab_h) if uy < 0 else 0.0)
    # BY ROLE, so an inbound spoke swaps them: its source is the SLOT and its target
    # is the hub. This is the whole of what the reversal changes here -- s_hh and
    # t_hh above are geometric and hand their terms to each other when uy flips --
    # and getting it wrong costs 2.0 units of caption displacement in a pitch that is
    # bisected against exactly that placement.
    _in = _inbound(entry)
    sR = _reach(c_hw, s_hh, ux, uy) + (cfg["tgt_margin"] if _in else cfg["src_margin"])
    tR = L - (_reach(w / 2, t_hh, ux, uy)
              + (cfg["src_margin"] if _in else cfg["tgt_margin"]))
    return (ux * sR, uy * sR), (ux * tR, uy * tR)


def row_gap_is_sound(products: dict, trunk: bool = False,
                     arrows_from: dict | None = None) -> bool:
    """May row_gap be solved for this figure at all?

    NO IF A CAPTIONED VERTICAL SPOKE STANDS BESIDE A DIAGONAL ONE, and that is
    the whole rule.

    A CAPTIONED DIAGONAL SPOKE ALONE IS NOT REFUSED. There the floor goes to 149.4
    where auto_row_gap's formula says 115.5, but that formula is INCOMPLETE, not
    invalid. The missing term is diagonal_caption_row --
    the block displaced along the shaft normal, whose vertical component carries it
    into the gap -- and with it in hand the floor is predicted to 0.00000 on all 44
    profile x caption-line-count points where that family binds, bisected on the
    renderer at a solved col_gap. So that case is solved rather than refused, and
    it covers the typical radial figure.

    WHAT GENUINELY IS NOT COVERED is a different pair. A caption on a T or B spoke
    sits in the COLUMN corridor, and the SHAFT of a diagonal spoke beside it sweeps
    through that corridor as the pitch changes. Measured on a captioned T with TR
    present: the pair collides over row 112 to 496 and over 140 to 496 on two
    profiles, so there is no pitch in any usable range that clears it, and it is not a
    floor that could be added to anything. Increasing the pitch does not help and
    neither does decreasing it.

    A CAPTIONED VERTICAL SPOKE ON ITS OWN IS FINE, which is why the rule names the
    pair and not the slot: T and B captioned together, with no diagonal present,
    collide only below row 24 -- far under any pitch this solver returns.

    The other caption-against-caption pairs measured -- TL with TR, TL with BL, L with
    TL -- band only at the extremes (under 68, or over 316) and a solved pitch of 50
    to 170 sits clear of them. They are not gated, and the allocator reports them if a
    caller pins a pitch that lands in one.
    """
    on_trunk = {t for m in trunked_slots(products, trunk, arrows_from).values()
                for t in m}
    axis = {s: spoke_axis(s, arrows_from, on_trunk) for s in products}
    if "diagonal" not in axis.values():
        return True
    return not any(axis[s] == "vertical" and (e.get("above") or e.get("below"))
                   for s, e in products.items())


def solve_row_gap(center: dict, products: dict, cfg: dict | None = None,
                  trunk: bool = False, arrows_from: dict | None = None,
                  col_gap: float | None = None) -> float:
    """The row pitch, from ALL THREE constraints that act on this axis.

    THREE DERIVATIONS, deliberately not merged into one expression:
      * auto_row_gap  -- a product's label against the STRUCTURE in the row below it.
        Product against product; nothing to do with any arrow.
      * vertical_spoke_row -- a vertical arrow running out of shaft between its trims.
        Arrow against its own two ends; nothing to do with any neighbour.
      * diagonal_caption_row -- a diagonal spoke's CAPTION, displaced along the shaft
        normal into the gap, against the hub or the row it lands in. Arrow's type
        against somebody else's structure.
    They bind on different figures -- the first on a side with two stacked slots, the
    second on a T or B with a tall product, the third on any captioned diagonal -- so
    the answer is the largest, and keeping them separate is what lets each be checked
    against its own measurement.

    col_gap IS REQUIRED BY THE THIRD ONE and by neither of the others: a diagonal
    spoke's direction is set by both pitches at once, so where its caption lands
    depends on the column. Left None the third term is skipped, which is what a caller
    that has not solved the column yet should get -- a term measured against the wrong
    column would be worse than an absent one. draw_radial_figure solves col first and
    passes it.
    """
    need = max(auto_row_gap(products, cfg) if products else 0.0,
               vertical_spoke_row(center, products, cfg, trunk, arrows_from))
    if col_gap is not None and products:
        need = max(need, diagonal_caption_row(center, products, col_gap, cfg,
                                              trunk, arrows_from))
    return need


def auto_col_gap(center: dict, products: dict, cfg: dict | None = None,
                 trunk: bool = False, arrows_from: dict | None = None) -> float:
    """Column pitch derived from the hub and the products. CLOSED FORM, no search.

    WHAT SETS THE FLOOR. Sweeping
    col_gap down over six product profiles, the first thing to go is the HORIZONTAL
    spoke: the arrow's own trims, measured inward from opposite ends, run out of shaft
    between them. For an L or R slot uy is 0, so both _reach calls collapse to the
    half-widths and the span is exactly

        span  =  col_gap - hub_half - src_margin - tgt_margin - w/2

    Setting that to MIN_ARROW_LEN and solving for col_gap is this function. Predicted
    against measured, in units, over the six profiles in test_radial_spacing:
    every residual is non-negative and inside one 2.5-unit sweep step.

    THE LABEL IS DELIBERATELY NOT IN THIS. The tempting reasoning -- "a label wider
    than its own structure gets clipped by the incoming arrow" -- rests on a benzene
    reading of 182.5 that is not a floor at all: sweeping at 1.0 resolution it is a
    six-unit ERROR BAND (175-180) with clean layout on BOTH sides, and it is the
    allocator's 2.0-unit polyline inflation grazing the label's bottom-left corner,
    not ink. Deflate that inflation and the shaft clears the label by 1.3-2.0 units
    everywhere inside the band. It also flips on ONE LETTER of label width -- "2TL"
    against "2TR" is 1.16 units and makes the band vanish. A constraint that moves
    that far on a letter is a measurement artifact, and encoding it would bake the
    artifact into every figure.

    WHY THE HORIZONTAL FORM COVERS THE DIAGONALS TOO, which is what makes one formula
    enough. For any off-centre slot |ux| <= 1 and hypot(col, row) = col / |ux|, so

        span  =  col/|ux| - hub_half/|ux|' - w/(2|ux|)' - src - tgt
              >= (col - hub_half - w/2) - src - tgt

    -- _reach is a min, so it never exceeds the half-width term, and dividing a
    non-negative quantity by |ux| <= 1 only makes it larger. The horizontal case is
    therefore the tightest of them, and solving it clears every other spoke. Verified
    numerically over all eight slots.

    Slots at x = 0 (T and B) are excluded: col_gap multiplies POS[p][0], which is 0
    for them, so it does not move them at all. A figure holding only those keeps the
    default 300 -- the value is unused, and returning it leaves such a figure
    byte-identical.

    THE CAPTION ALSO CLAIMS SHAFT, and on a HORIZONTAL spoke it is usually the term
    that binds. Its normal there is (0, -1), so the block is centred on the shaft
    MIDPOINT in x and reaches w_cap/2 each way. Writing A = hub_half + src_eff for the
    tail's distance from the hub centre and B = col - w/2 - tgt_margin for the head's,
    the midpoint is at (A + B) / 2 and the two obstacles give

        vs the hub edge      at hub_half:      span >= w_cap - 2 * src_eff
        vs the product edge  at col - w/2:     span >= w_cap - 2 * tgt_margin

    so the binding one is 2 * min(src_eff, tgt_margin) -- with this engine's 24 and 22,
    2 * tgt_margin = 44. DERIVED, not fitted: predicted against measured over seven
    caption widths from 27.0 to 120.6 units, every residual inside 3e-5, and there is
    no tolerance term in it.

    A TRUNK BRANCH TAKES THE SAME 44 even though its src_pad is 0, because its tail
    sits on a thin LINE rather than against a box -- the trunk_lead corridor behind it
    is deliberately empty (see the ink note in render_radial), so the hub-side bound
    never binds and the product side decides again. Measured: identical spans, and the
    src_eff = 0 form over-predicts by exactly 44.

    MORE EXACTLY, A TRUNK BRANCH IS THE SAME FORMULA WITH A NEGATIVE src_eff, not
    the same 44 -- what decides it is what the caption actually hits. Its tail sits
    ON the trunk, and the trunk is a drawn LINE reserving 1.0 unit
    either side (emit_bus, ink=True; read off it as TRUNK_INK_HALF rather than
    restated). So the clearance behind the tail is not 0 but -1.0, and the bound is
    span >= w_cap + 2 * TRUNK_INK_HALF = w_cap + 2. Measured: the L and BL captions on
    a three-slot trunk cross the spine and the stub at exactly that point.

    IT OVER-SIZES ONE CASE ON PURPOSE. The TOPMOST trunked slot's "above" caption sits
    beyond the spine's upper end and could not cross it whatever its width -- measured
    clean where L and BL both collide. Excluding it would mean predicting where the
    spine stops against how tall the caption is, which is a rule that turns over on
    caption height; a figure captioned only on its top slot therefore gets a slightly
    wider grid than it strictly needs, and every residual stays non-negative.

    AN UNROTATED DIAGONAL SPOKE DOES NOT NEED IT. There the caption is displaced along
    the shaft normal by cap_gap + _support(...), which carries it clear of both
    obstacles; measured at every caption width above, on a diagonal slot the floor
    never moves at all. So the term is applied only where the arrow really is
    horizontal -- the L and R slots, and any slot hanging off a trunk.

    A ROTATED ONE NEEDS A DIFFERENT TERM, and it is a different CONSTRAINT rather than
    a corrected coefficient. Set along the shaft, the block is COLLINEAR with it, so
    it no longer overhangs sideways past the two ends and there is nothing for the two
    ends to absorb: the shaft's own span has to CONTAIN the caption's along-axis
    extent. That extent is the advance box's WIDTH -- not a projection of it, because
    a shaft-aligned box's long axis IS the shaft.

        rotated diagonal:    span >= w_cap
        horizontal:          span >= w_cap - 2 * min(src_eff, tgt_margin)

    So the 44 simply is not there. It is not a fudge factor; it is the sideways
    overhang the two margins swallow, and a collinear block has no sideways overhang.

    ABOVE AND BELOW ARE TWO INDEPENDENT BLOCKS, so the binding extent is
    max(w_above, w_below) and never their sum. That is already what caption_width
    returns and already how caption_box and _caption_clash treat them -- one call and
    one shape per side -- so a split caption is correctly cheaper than a merged one on
    both axes.

    IT IS SOLVED ON THE HORIZONTAL FORM AND THAT CLEARS EVERY DIAGONAL, which is what
    lets the column keep being solved BEFORE the row even though a diagonal's angle
    depends on both pitches. Same argument as above, one step further: with
    D = col/|ux| and reach = min(hw/|ux|, hh/|uy|) <= hw/|ux|,

        span_diag  >=  (col - hub_hw - prod_hw)/|ux| - src - tgt
                   >=  (col - hub_hw - prod_hw)      - src - tgt   =  span_horizontal

    because |ux| <= 1 and the bracket is non-negative. The horizontal form therefore
    never OVER-predicts a diagonal span, so a col that satisfies it satisfies every
    diagonal AT EVERY ROW. Measured over six profiles x four line counts x five
    columns x six rows: min(L_diag - span_horizontal) = +0.245, never negative. So
    col-then-row survives and there is no fixed point to iterate.

    AND THE ARGUMENT IS A ONE-SIDED BOUND, NOT A MONOTONICITY CLAIM, which matters
    because both the tempting monotonicity claims are FALSE. Measured: the containment
    set is a union of BANDS in col as well as in row (up to three transitions over a
    60..620 sweep at a small fixed row -- at a nearly-vertical spoke the span is long,
    shrinks as the spoke swings toward horizontal, then grows again), and the slack is
    not monotone in row either (8 of 36 swept cases decrease). Neither breaks anything
    here, and that is the point of solving a SUFFICIENT condition rather than
    inverting the predicate: "span_horizontal >= w_cap" needs no ladder because it is
    not looking for the top of a band. The row axis, which really does invert its
    clearance, needs the ceiling walk. This does not.

    WHY THIS IS NOT BISECTED ON THE RENDERER like the terms above it. A caption that
    overruns its arrow lands in EMPTY PAGE -- measured, the allocator reports nothing
    at all -- so the report is not an oracle for it, exactly as it is not for a
    reversed arrow (see radial_arrow_guard). The predicate is containment,
    w_cap <= span, checked on the renderer's own geometry instead.
    """
    cfg = {**CFG, **(cfg or {})}
    off = [(s, e) for s, e in products.items() if POS[s][0] != 0]
    if not off:
        return cfg["col_gap"]
    base = (center_half_width(center, cfg, trunk) + cfg["src_margin"]
            + cfg["tgt_margin"])
    on_trunk = {t for m in trunked_slots(products, trunk, arrows_from).values()
                for t in m}
    need = 0.0
    for slot, e in off:
        bb = e["frag"][1]
        # PER SLOT, not w_max against cap_max. The widest structure and the widest
        # caption need not be the same product, and pairing them would size the grid
        # for a slot that does not exist.
        span = MIN_ARROW_LEN
        if POS[slot][1] == 0 or slot in on_trunk:
            # WHAT SITS BEHIND THE TAIL decides the backward bound. A hub spoke has
            # src_margin of empty corridor there; a trunk branch has the trunk's own
            # LINE, which the caption may not cross, so its clearance is NEGATIVE.
            src_eff = (-TRUNK_INK_HALF if slot in on_trunk else cfg["src_margin"])
            # CLEAR the obstacle rather than land on it. Solving for equality puts the
            # caption's edge and the trunk's edge on the same coordinate, and the
            # allocator holds a caption against a decoration at ZERO tolerance -- so a
            # 1.4e-14 float residue reads as a collision. The arrow guard could be
            # loosened because it is this engine's; the allocator's check is not, so the
            # solver gives way instead. _SPAN_TOL is eleven orders below the 0.05 the emitter
            # can express, so nothing that reaches the file is being bought.
            span = max(span, caption_width(e, cfg)
                       - 2.0 * min(src_eff, cfg["tgt_margin"]) + _SPAN_TOL)
        elif cfg["cap_rotate"]:
            # A ROTATED DIAGONAL. The block is collinear with the shaft, so the shaft
            # has to contain it outright -- no sideways overhang for the margins to
            # absorb, hence no 44. See the derivation above. `off` has already dropped
            # T and B, and a trunked slot took the branch above, so this is exactly
            # the four corner slots that really do turn.
            span = max(span, caption_width(e, cfg) + _SPAN_TOL)
        need = max(need, base + (bb[2] - bb[0]) / 2.0 + span)
    return need


def render_radial(recipe: dict, header_prefix: str, cfg: dict | None = None,
                  page_mode: str = "fit", alloc=None,
                  shift: list | None = None,
                  footnote_report: dict | None = None) -> str:
    cfg = {**CFG, **(cfg or {})}
    # THE ENGINE'S OWN GATE, and it RAISES rather than reports. cfg goes straight to
    # arithmetic a few lines down -- `POS[p][0] * cfg["col_gap"]` -- and POS holds
    # -1 / 0 / 1, so a string arriving here does not fail, it MULTIPLIES: -1 * "300"
    # is '', 0 * "300" is '', 1 * "300" is '300'. A product would be placed at a
    # string coordinate and the figure would come out silently wrong with nothing to
    # point at, which is the same failure mode as a dropped recipe key.
    #
    # The tool layer coerces (see _gap_value there), but the engine is called directly
    # by every test in this suite and by the demo scripts, so it cannot assume that.
    # "auto" is resolved ABOVE this line, at the tool layer, against the products it
    # needs to measure -- the engine takes numbers.
    for _k in ("col_gap", "row_gap"):
        _v = cfg[_k]
        # bool first: True is an int in Python and would pass as a 1.0-unit gap.
        if isinstance(_v, bool) or not isinstance(_v, (int, float)) or _v <= 0:
            raise ValueError(
                'cfg[%r] must be a positive number, got %r (%s). "auto" is resolved '
                "at the tool layer; render_radial takes numbers only."
                % (_k, _v, type(_v).__name__))
    # Load-bearing, not a shadow: this IS the page-extent bookkeeping.
    alloc = alloc if alloc is not None else _Allocator()
    prod = recipe.get("products", {})

    def whf(entry):
        bb = entry[1]
        return bb[2] - bb[0], bb[3] - bb[1]

    # grid centers
    cen = {"center": (0.0, 0.0)}
    for p in prod:
        cen[p] = (POS[p][0] * cfg["col_gap"], POS[p][1] * cfg["row_gap"])

    parts, texts, arrows, graphics = [], [], [], []
    # Halo / dot / locants, keyed the way this engine keys its products --
    # by SLOT ("TL", "B", ... plus "center"), because that is the only name a
    # radial product has here. See layout_linear for why `behind` is joined first.
    D = recipe.get("decorations") or {}
    behind: list = []
    base, tid, aid, gid = 1000, 91000, 90000, 96000

    def put_box(x0, y0, x1, y1, owner="?", kind="structure", shape=None):
        # `shape` is the ORIENTED ink for a rotated caption; x0..y1 stay its hull, so
        # the page extent is the same question as for any caption. None for everything
        # else, which is every caller but the rotated-caption one.
        res = alloc.reserve(_Rect(x0, y0, x1 - x0, y1 - y0), owner, kind,
                            shape=shape)
        alloc.check_emit(res.rid, owner)
        return res

    def place(entry, cx, cy, owner="?", deco=None):
        nonlocal base
        frag, bb = entry
        w, h = bb[2] - bb[0], bb[3] - bb[1]
        _fx, _dx = C.place_with(frag, deco, base, cx - w / 2, cy - h / 2)
        parts.append(_fx); base += 1000
        if _dx:
            behind.append(_dx)
        put_box(cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2,
                owner=owner, kind="structure")
        return w, h

    def put_label(cx, top, text, size, owner="?"):
        """Centered styled label; top = box top. Returns (box bottom, height, width).

        `text` is one line or a list of them. Every label this engine writes itself
        is a single line, but a CENTRE label is routinely two or three in published
        figures ("2a" / "70%" / "from 2.0 mmol of 1a"), and both styled_text and
        _block_box take a list -- so this does not flatten it either.

        The HEIGHT is returned because the arrow trim needs it: an arrow entering a
        product from below has to stop clear of that product's label, and the only
        honest figure for "how tall is that label" is the one measured here. The
        WIDTH, because a box drawn round the centre has to enclose a label wider
        than the structure it sits under.
        """
        nonlocal tid
        lines = [text] if isinstance(text, str) else [str(t) for t in text]
        asc = ASCENT_EM * size
        parts_local = C.styled_text(tid, cx, top + asc, lines, size=size, align="center")
        texts.append(parts_local); tid += 1
        # Real metrics. An estimate of len(text_str) * 0.6 * size would count the
        # #bold# / ^sup^ markers as if they printed -- so a label with a footnote
        # marker would reserve room for characters that never appear -- and would
        # pin the height at exactly `size`.
        w, h = _block_box(lines, size)
        put_box(cx - w / 2, top, cx + w / 2, top + h, owner=owner, kind="label")
        return top + h, h, w

    # center
    cw, ch = place(recipe["center"]["frag"], 0.0, 0.0, owner="center",
                   deco=D.get("center"))

    # products + labels
    lab_h: dict = {}          # slot -> the MEASURED height of its label block
    for p, d in prod.items():
        cx, cy = cen[p]
        w, h = place(d["frag"], cx, cy, owner="prod:%s" % p, deco=D.get(p))
        _bot, lab_h[p], _lw = put_label(cx, cy + h / 2 + cfg["label_gap"],
                                        label_text(d), cfg["num_size"],
                                        owner="prod:%s.label" % p)

    # ---- the centre's own label: UNDER the structure, horizontally centred ---
    # Not at (cw/2 + 20, -ch/2 - 14) -- two written-in offsets that would put it
    # off the structure's top-right corner. Three things are wrong with that, and
    # none of them is cosmetic:
    #
    #   IT SITS IN THE TRUNK'S LANE. A trunk stands src_margin (24) outside the hub
    #   and such a label is CENTRED 20 outside it, so the two overlap by 4 before
    #   the label has any width at all, and every unit of half-width past that goes
    #   straight in. Measured on two different hubs: a 9-character label overlapped
    #   the trunk 7 wide / 7.2 deep in BOTH of them, a 14-character one 25 / 10.2 --
    #   the same readings from different hub widths, because the depth was set by
    #   the label's own length and by nothing about the hub. (Untrunked figures hide
    #   it: a diagonal spoke leaves a ~25-wide corridor there, which caps the
    #   reading, so the monotone relation does not show.)
    #
    #   IT IS ONE LINE. A published centre label is routinely two or three -- "2a /
    #   70% / from 2.0 mmol of 1a", "3a, 72% / 1.39 g".
    #
    #   IT STRADDLES ITS OWN BOX instead of sitting inside it.
    #
    # Under the structure, at label_gap, centred, is where this engine ALREADY puts
    # a compound's id and yield -- see the product labels above. The centre
    # follows that rule rather than a second one of its own. The position is the one
    # published figures use: the centre label below the structure and inside the
    # box, not at a corner.
    cen_lab_h = cen_lab_w = 0.0
    cen_lab_bot = ch / 2
    if recipe["center"].get("label"):
        cen_lab_bot, cen_lab_h, cen_lab_w = put_label(
            0.0, ch / 2 + cfg["label_gap"], recipe["center"]["label"],
            cfg["cen_size"], owner="center.label")

    # center box (Line-based), padded around the center bbox and the label under it
    box_half = box_below = None
    if recipe["center"].get("box"):
        # max(), because a two-word label is easily wider than a small hub and a box
        # has to hold what is drawn inside it.
        _hw = max(cw / 2, cen_lab_w / 2)
        x0 = -_hw - cfg["box_pad_x"]; y0 = -ch / 2 - cfg["box_pad_t"]
        x1 = _hw + cfg["box_pad_x"]; y1 = cen_lab_bot + cfg["box_pad_b"]
        graphics.append(C.box_lines(gid, x0, y0, x1, y1)); gid += 4
        put_box(x0, y0, x1, y1, owner="center.box", kind="decoration")
        # (half width, extent ABOVE the origin), and then what hangs BELOW it --
        # not (w/2, h/2). The box is not symmetric about the hub's centre once a
        # label grows under it, so an arrow leaving downwards has to
        # clear the bottom edge while one leaving upwards must not be pushed out to
        # meet it. draw_arrow already carries exactly that asymmetry for products
        # (half-height, plus whatever hangs below); the hub uses the same pair.
        box_half = ((x1 - x0) / 2, -y0)
        box_below = y1 - (-y0)

    # arrows: center -> product, or product -> product (secondary step)
    af = recipe.get("arrowsFrom", {})
    center_half = box_half if box_half else (cw / 2, ch / 2)
    center_below = (box_below if box_half is not None else
                    (cfg["label_gap"] + cen_lab_h if cen_lab_h else 0.0))
    # WHERE A TRUNK MAY STAND, which is not the same question as where a spoke's
    # tail starts. A spoke is trimmed to the hub's STRUCTURE and a label wider than
    # its own structure is left to overhang -- that is this engine's convention for
    # products too, and col_gap's docstring names the consequence. A trunk is a LINE
    # DRAWN DOWN THE PAGE, so the same overhang is not a trim to argue about: the
    # line goes through the type. Measured with no box, a 43.1-wide hub and a
    # 144.3-wide centre label: the trunk stands at x 45.6 and the label reaches
    # 72.1, so 26 units of it are struck through.
    #
    # So the trunk clears the hub's whole FOOTPRINT. With a box this is already what
    # center_half says (the box encloses the label plus box_pad_x), and with a label
    # no wider than the structure it is the structure -- so this can only move a
    # figure that would otherwise draw a line through its own label.
    trunk_hw = max(center_half[0], cen_lab_w / 2.0)

    def draw_arrow(from_xy, from_half, from_below, to_xy, to_half, to_below,
                   owner="arrow", cap=None, src_pad=None):
        """Returns the arrow's own reservations, so a trunk can allow-list them.

        `from_below` / `to_below` are how far each end's ink reaches BELOW its own
        half-height: for a product that is the gap plus its compound label, for a
        boxed hub it is the box's lower half. One number rather than a label height
        the callee adds a gap to, because the hub's is not a label -- and a caller
        that measures its own overhang cannot disagree with the trim about it.

        `src_pad` overrides how far the tail sits outside the source. It is
        src_margin for a spoke leaving the hub, and 0.0 for a branch leaving a
        trunk -- a branch tail must sit ON the trunk line, not a margin off it,
        or the line it is supposed to hang from does not touch it.
        """
        nonlocal aid, tid
        dx, dy = to_xy[0] - from_xy[0], to_xy[1] - from_xy[1]
        L = math.hypot(dx, dy) or 1.0
        ux, uy = dx / L, dy / L
        # MEASURED, not a constant. label_gap + a flat label_h of 16.0 is not the
        # height of anything: a real "9a, 94%" block is 8.70 units and one carrying
        # a footnote marker is 9.27, so a flat trim would reserve room for roughly
        # 7 units of label that does not exist and every vertical spoke would stop
        # that far short of the product it points at. The label is measured
        # in put_label; this uses that number.
        s_hh = from_half[1] + (from_below if uy > 0 else 0.0)
        t_hh = to_half[1] + (to_below if uy < 0 else 0.0)
        sR = _reach(from_half[0], s_hh, ux, uy) + (
            cfg["src_margin"] if src_pad is None else src_pad)
        tR = L - (_reach(to_half[0], t_hh, ux, uy) + cfg["tgt_margin"])
        x0, y0 = from_xy[0] + ux * sR, from_xy[1] + uy * sR
        x1, y1 = from_xy[0] + ux * tR, from_xy[1] + uy * tR
        # THE SHAFT HAS TO STILL BE A SHAFT, and this is the one thing on the page
        # that nothing else can check. sR and tR are trims measured INWARD from
        # opposite ends and nothing stops them crossing: at col_gap 147.5 a 92.8-wide
        # product puts tR 19.6 units BEHIND sR, and the arrow is emitted with its head
        # on the hub side -- pointing back at its own source. Measured, on the L and R
        # spokes of the profile sweep.
        #
        # THE ALLOCATOR CANNOT SEE IT. Both endpoints are legal points, the polyline
        # pieces are legal boxes in legal places, and a reversed arrow overlaps
        # nothing, so the figure draws "clean" while stating the reverse chemistry.
        # That is why the check lives here, at the only point where the intended
        # source->target sense and the drawn sense are both in hand.
        #
        # FLAGGED AT MIN_ARROW_LEN RATHER THAN AT ZERO: a 3-unit shaft behind a
        # 10-unit arrowhead is already not a drawing anyone would keep, and catching
        # it one bond-half early means the report names the figure before it
        # degenerates rather than after.
        _span = tR - sR
        if _span < MIN_ARROW_LEN - _SPAN_TOL:
            alloc.note(
                _Severity.ERROR, "arrow-too-short",
                "%r spans %.1f units between its trims, under the %.1f minimum "
                "(arrowhead %.1f + half a bond)%s"
                % (owner, _span, MIN_ARROW_LEN, HEAD_LEN,
                   " -- NEGATIVE, so it is drawn pointing back at its source"
                   if _span < 0 else ""),
                amount=MIN_ARROW_LEN - _span)
        arrows.append(C.arrow_xy(aid, x0, y0, x1, y1)); aid += 1
        # A radial spoke is a diagonal line: its axis-aligned hull is mostly empty
        # space, so reserve the run in short pieces the way connectors do.
        _res = alloc.reserve_polyline([(x0, y0), (x1, y1)], owner, "arrow")

        # ---- captions, on the shaft's normal ---------------------------------
        cap = cap or {}
        if not (cap.get("above") or cap.get("below")):
            return _res
        for side in ("above", "below"):
            lines = cap.get(side)
            if not lines:
                continue
            lines, (tx, ty, w_, h_), _rot = caption_box((x0, y0), (x1, y1),
                                                        lines, side, cfg)
            _bb, _shape, _p = caption_shape(tx, ty, w_, h_, _rot,
                                            cfg["cap_size"])
            texts.append(C.styled_text(tid, _p[0], _p[1], lines,
                                       size=cfg["cap_size"], align="center",
                                       rot=_rot))
            tid += 1
            # Owner extends the ARROW's own name, so allocator.related() sees one
            # object (it compares the head of the dotted path) and a caption
            # grazing the shaft it annotates is not reported -- while it stays
            # strict against every structure, label and other arrow on the page.
            put_box(_bb.x, _bb.y, _bb.right, _bb.bottom,
                    owner="%s.%s" % (owner, side), kind="caption",
                    shape=_shape)
        return _res

    # ---- which slots share a trunk ------------------------------------------
    # OPT-IN, and it has to be: switching a side to a trunk turns its diagonal
    # spokes into horizontal ones, which is a different drawing. There is no
    # version of this that leaves an existing figure byte-identical the way the
    # caption fold did, so the only way an existing figure cannot move is for the
    # feature not to run unless it is asked for.
    #
    # TWO OR MORE PER SIDE. A trunk shared by one arrow shares nothing: emit_bus
    # spans row_base[owner]..row_base[last], those collapse to one y, and it emits
    # a zero-length spine -- BoundingBox="84.00 0.00 84.00 0.00", measured -- laid
    # on top of the stub. So a lone slot keeps its own spoke, which is also the
    # right drawing: one arrow off a hub is one arrow off a hub.
    #
    # A slot fed by another PRODUCT (arrowsFrom) is never trunked. Its arrow does
    # not start at the hub, so it has nothing to do with the hub's trunk.
    trunked = trunked_slots(prod, bool(recipe.get("trunk")), af)

    # The trunk stands where the arrow tails ALREADY stand: src_margin outside the
    # hub, the same quantity that decides where an untrunked spoke leaves the hub.
    # Nothing new is measured, and turning the trunk on therefore does not move a
    # single tail -- it only joins them up.
    trunk_lead = cfg["src_margin"]
    side_sign = {"left": -1.0, "right": 1.0}
    trunk_res: dict = {s: [] for s in trunked}

    def _below(slot):
        """How far a product's ink reaches below its own half-height: the gap, plus
        the compound label hanging in it. Measured in put_label, not re-guessed."""
        return cfg["label_gap"] + lab_h.get(slot, 0.0)

    for p, d in prod.items():
        w, h = whf(d["frag"])
        to_half = (w / 2, h / 2)
        _cap = {k: d[k] for k in ("above", "below") if d.get(k)}
        _side = next((s for s, m in trunked.items() if p in m), None)
        if p in af:
            if _inbound(d):
                raise ValueError(
                    'products[%r] sets direction="in" and is also listed in '
                    "arrowsFrom. `direction` says which way the spoke BETWEEN THIS "
                    "SLOT AND THE HUB points, and an arrowsFrom slot has no such "
                    "spoke -- its arrow starts at product %r. Drop one of the two."
                    % (p, af[p]))
            src = af[p]
            sw, sh = whf(prod[src]["frag"])
            draw_arrow(cen[src], (sw / 2, sh / 2), _below(src),
                       cen[p], to_half, _below(p),
                       owner="arrow:%s->%s" % (src, p), cap=_cap)
        elif _side:
            # A branch off the trunk is an ORDINARY HORIZONTAL arrow, which is the
            # whole reason the trunk is worth drawing: on a horizontal shaft the
            # caption's "above" degenerates to plain screen-up, and the tail needs
            # no diagonal reasoning at all. It leaves a POINT on the trunk, so the
            # source half-size is (0, 0) and src_pad is 0.
            _tx = side_sign[_side] * (trunk_hw + trunk_lead)
            trunk_res[_side] += draw_arrow(
                (_tx, cen[p][1]), (0.0, 0.0), 0.0,
                cen[p], to_half, _below(p),
                owner="arrow:trunk-%s->%s" % (_side, p), cap=_cap, src_pad=0.0)
        elif _inbound(d):
            # SWAP THE ARGUMENTS, do not flip the arrowhead -- the module docstring
            # says why. draw_arrow trims by ROLE, so handing it (slot, hub) puts
            # src_margin on the source (the slot) and tgt_margin on the target (the
            # hub), which is what an arrow pointing AT the hub means. The hub's own
            # overhang is now the TARGET's, and the slot's is the source's, which is
            # the same swap one level down.
            draw_arrow(cen[p], to_half, _below(p),
                       cen["center"], center_half, center_below,
                       owner="arrow:%s->center" % p, cap=_cap)
        else:
            # The centre's overhang is real: its optional label hangs under
            # it exactly as a product's does, and a downward spoke has to start
            # clear of that (or of the box drawn round both). 0.0 only for a bare
            # centre.
            draw_arrow(cen["center"], center_half, center_below,
                       cen[p], to_half, _below(p),
                       owner="arrow:center->%s" % p, cap=_cap)

    # ---- the trunk itself: spine + stub, from rowspan.emit_bus ---------------
    for side, members in trunked.items():
        sgn = side_sign[side]
        shaft_ys = [cen[s][1] for s in members]
        # fork_stub_y is the rule: odd n is colinear with the middle shaft, even n
        # takes the middle PAIR's midpoint. Every vertically symmetric side -- and
        # on a three-row grid that is every side except two -- comes back 0, the
        # hub's own centre line.
        #
        # THE TWO IT CANNOT ANSWER HERE are {top, middle} and {middle, bottom}:
        # it returns +-142.5, which is colinear with no arrow AND misses the hub
        # entirely (the hub's half-height is nearer 40), so the stub would start in
        # mid-air beside the hub's corner -- the exact "nothing touches it" defect
        # emit_bus exists to prevent. fork_stub_y is not wrong; its even-n rule is
        # written for a spanning cell that align_y MOVES to meet the stub, and
        # radial's hub is pinned at the origin and cannot be moved. So when its
        # answer is off the hub, the stub takes the hub's own centre -- which for
        # both those cases is colinear with the middle branch anyway, the very
        # property fork_stub_y is reaching for.
        cy = _fork_stub_y(shaft_ys)
        if abs(cy) > center_half[1]:
            cy = 0.0
        # BOTH SIDES ARE THE FORK FORM, and the left one is not the merge form --
        # see the module docstring. Mirrored by a negative lead, so bx lands inside
        # the hub-ward side and the stub stays a plain line.
        sp = dict(owner=0, last=len(members) - 1,
                  x0=sgn * trunk_hw, x1=sgn * trunk_hw,
                  cy=cy, lead=sgn * trunk_lead, fork=True)
        # ink=True: the SPINE and the STUB as themselves, not the L-shaped hull
        # they span. The hull's empty corner is this figure's corridor -- the
        # trunk_lead (24) units between the hub's edge and the trunk -- so a hull
        # reservation declares that corridor occupied and reports anything placed
        # in it as sitting on a line 24 units away. Measured on a published radial
        # figure at four label lengths: the hull came back 26.0 x 320.7 and the
        # centre label was reported 9.6 deep in it, where 9.6 is the LABEL's own
        # height rather than any distance to the trunk -- so cutting the label
        # from 33 characters to 9 moved the reading by 0.0, and only pushing the
        # trunk further out (adding the centre box) changed it. That is the
        # signature of a reservation that is not the shape of its ink.
        #
        # THE PAGE EXTENT DOES NOT MOVE. The spine already spans the full row
        # range and the stub the full corridor width, so the union of the two is
        # the hull to the unit; only the interior of the corner is given back.
        _els, _boxes, gid = _emit_bus([sp], None, None, shaft_ys,
                                      {"gap": trunk_lead}, gid, ink=True)
        graphics.extend(_els)
        # The branch arrows are MEANT to touch the trunk -- that is what hanging
        # off it means -- so they are allowed rather than reported. Without this
        # each branch tail scores ~3 units into the trunk box against a 2.0
        # decoration/arrow tolerance and the report fills with by-design contacts.
        # The trunk's OWN pieces join the list as they are reserved: the spine and
        # the stub cross at the hub's centre line by construction, and one drawing
        # touching itself is not a finding either.
        _allow = [r.rid for r in trunk_res[side]]
        for _b in _boxes:
            _r = alloc.reserve(_Rect(_b[0], _b[1], _b[2] - _b[0], _b[3] - _b[1]),
                               "trunk:%s" % side, "decoration", allow=_allow)
            alloc.check_emit(_r.rid, "trunk:%s" % side)
            _allow.append(_r.rid)

    # footnotes, under everything. Placed LAST because "under everything" is a
    # question only the allocator can answer, and it can only answer it once every
    # structure, label, arrow and caption has been reserved.
    notes = recipe.get("footnotes") or []
    if isinstance(notes, str):
        notes = [notes]
    if notes:
        _cb = alloc.bbox()
        ny = _cb.bottom + cfg["footnote_gap"]
        # THE FIGURE THE PROSE RUNS UNDER -- this call's own width unless the caller
        # is composing and says which figure that really is. See CFG's footnote_width.
        # _cb is read ONCE here, before the first line is emitted -- put_box reserves,
        # and re-reading the bbox inside the loop would let each footnote widen the
        # cap and push the next one down.
        _ncap = cfg.get("footnote_width") or _cb.w
        _nlines, _nover, _nj = [], [], 0
        for i, note in enumerate(notes):
            # FOLD FIRST, then measure -- and take both from one call, which is the
            # same invariant the arrow captions in this file keep, for the same reason:
            # the box reserved below has to be the box of the lines that actually get
            # drawn. Each note is folded on its own, so a deliberate break survives and
            # only an over-wide note is broken further.
            _lines, (_bw, _) = _fold_box([str(note)], cfg["footnote_size"], _ncap)
            _nlines.append(len(_lines))
            # block_box's width IS the widest of the folded lines, so ONE comparison
            # says whether this note still has something too wide in it at all.
            _wide = _bw > _ncap + 0.01
            for k, line in enumerate(_lines):
                w_, h_ = _block_box([line], cfg["footnote_size"])
                # A WORD WIDER THAN THE FIGURE CANNOT BE FOLDED -- fold_lines gives it
                # a line of its own rather than hyphenating it -- so the only honest
                # thing left is to SAY SO. Silence here is a footnote running off the
                # side of the figure while the reply reports a clean drawing.
                if _wide and w_ > _ncap + 0.01:
                    _nover.append({"note": i, "line": k, "w": round(w_, 2),
                                   "cap": round(_ncap, 2), "text": line})
                # LEFT, and the x is the figure's own left edge. table_layout found that
                # centring footnotes and back-solving for a left edge leaves a ragged
                # margin, because every note is a different width.
                texts.append(C.styled_text(
                    tid, _cb.x, ny + ASCENT_EM * cfg["footnote_size"], [line],
                    size=cfg["footnote_size"], align="left"))
                tid += 1
                # PER FOLDED LINE, not per note: a folded note is several boxes, and
                # something landing between two of its lines has to surface as a
                # collision. layout_grid numbers its footnote reservations the same
                # way, so the two engines' reports read alike.
                put_box(_cb.x, ny, _cb.x + w_, ny + h_,
                        owner="footnote[%d]" % _nj, kind="caption")
                _nj += 1
                ny += h_ + cfg["footnote_lead"]
        # THE REPLY'S ACCOUNT OF THE FOLD, same rule as the gaps: a caller cannot see
        # a fold in the file, and `over` non-empty is the only way an unfoldable word
        # reaches them. Absent when nobody passed a dict, so the engine's own tests
        # and the demo scripts are unaffected.
        if footnote_report is not None:
            footnote_report.update({"width": round(_ncap, 2), "lines": _nlines,
                                 "over": _nover})

    # finalize
    body = "\n".join(behind + graphics + parts + arrows + texts)
    if behind:
        body, _z = C.stamp_z(body)     # only when there is something to layer
    # One account. alloc.bbox() is the page extent -- there is no second list that
    # could disagree with it.
    _bb = alloc.bbox()
    content_bbox = (_bb.x, _bb.y, _bb.right, _bb.bottom)
    return C.finalize(header_prefix, body, content_bbox,
                      page_mode=page_mode, margin=cfg["margin"], shift=shift)
