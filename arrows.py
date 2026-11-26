"""
arrows.py — route x style arrows for the desktop CDXML pipeline.

Two orthogonal dimensions:
  route (path shape): straight | elbow_hv | elbow_vh | via([[x,y],...])
  style (look):       kind + heads + headType + line + noGo + curve + sizes

Verified ChemDraw arrow attributes:
  ArrowheadType : Solid | Hollow | Angle
  ArrowheadHead / ArrowheadTail : Full | HalfLeft | HalfRight | (none)
  ArrowShaftSpacing : double shaft (retro / equilibrium)
  NoGo : Cross | Hash        AngularSize : curved (deg)      LineType : Dashed | Bold | Wavy

kind: 'normal' | 'retro'(=> open-V double shaft) | 'equilibrium'(<=> two half arrows)
      | 'harpoon'(single half arrow) | 'fishhook'(curled harpoon)
      | 'noGoCross' | 'noGoHash' | 'resonance'(<-> = heads:'both')
heads: 'end'(default) | 'both' | 'none' | 'start'   headType: 'solid'|'hollow'|'angle'
line: 'solid'|'dashed'|'bold'|'wavy'   barb: 'HalfLeft'|'HalfRight'   curve: degrees
      NOTE 'wavy' is accepted but ChemDraw stores and then IGNORES LineType="Wavy" on an
      arrow (verified: survives a round trip, renders identically to plain). 'bold' does
      work; for thicker still, set LineWidth or BoldWidth per arrow.

Endpoint edge-trimming (straight/elbow/via) uses each node's center + half-size
(hw, hh) plus a margin so the arrow stops just outside the box.
"""

from __future__ import annotations

import math
import re

import cdxml_core as C

# ONE DEFINITION, aliased -- the same arrangement the key sets have with preflight.
# `apply_style` here and `C.arc_arrow` both write `LineType` on an <arrow>, so they
# share one transcription of the mapping. It lives in cdxml_core because the import
# runs that way.
line_type = C.line_type


def _set_attr(s: str, attr: str, val) -> str:
    pat = re.compile(r'\b' + attr + r'="[^"]*"')
    if pat.search(s):
        return pat.sub(f'{attr}="{val}"', s)
    return s.replace('FillType="None"', f'FillType="None" {attr}="{val}"', 1)


def set_attrs(s: str, attrs: dict) -> str:
    """Write arbitrary attributes onto one arrow element. Public because the round-trip
    harness needs to write SPECULATIVE attributes -- ones the style vocabulary does not
    cover yet -- and let ChemDraw adjudicate: if the attribute survives an open+save it
    is real and supported, if it is dropped it is not. Values are not validated."""
    for k, v in attrs.items():
        s = _set_attr(s, k, v)
    return s


def apply_style(s: str, style: dict | None = None, is_first: bool = True, is_last: bool = True) -> str:
    """Apply a style dict to one arrow-segment string. is_first/is_last mark which
    physical end(s) of the whole route this segment touches (single seg: both True)."""
    style = dict(style or {})
    kind = style.get("kind", "normal")
    if kind == "resonance":
        kind, style = "normal", {**style, "heads": "both"}
    heads = style.get("heads", "end")
    type_map = {"solid": "Solid", "hollow": "Hollow", "angle": "Angle"}
    head_type = style.get("headType") or ("angle" if kind == "retro" else "solid")
    s = _set_attr(s, "ArrowheadType", type_map.get(head_type, "Solid"))

    if kind == "equilibrium":
        if not is_last:
            s = _set_attr(s, "ArrowheadHead", "(none)")
        else:
            # BOTH barbs are HalfLeft: "left" is relative to each end's own direction
            # of travel, so HalfLeft/HalfLeft puts them on opposite visual sides (the
            # correct equilibrium glyph). HalfLeft/HalfRight puts both barbs on the
            # SAME side. Measured from a hand-drawn ChemDraw export of an
            # equilibrium arrow.
            s = _set_attr(s, "ArrowheadHead", "HalfLeft")
            s = _set_attr(s, "ArrowheadTail", "HalfLeft")
        # ArrowShaftSpacing is in HUNDREDTHS of a unit, like HeadSize. ChemDraw uses
        # 300 (= 3.00 units). A value of 4 would mean 0.04 units, i.e. the two
        # shafts drawn on top of each other.
        s = _set_attr(s, "ArrowShaftSpacing", style.get("shaftSpacing", EQUIL_SHAFT_SPACING))
        s = _set_attr(s, "HeadSize", style.get("headSize", 1500))
        s = _set_attr(s, "ArrowheadCenterSize", 1313)
        s = _set_attr(s, "ArrowheadWidth", style.get("width", 375))
    elif kind in ("harpoon", "fishhook"):
        s = _set_attr(s, "ArrowheadHead", (style.get("barb", "HalfLeft") if is_last else "(none)"))
        if kind == "fishhook" and is_last:
            s = _set_attr(s, "AngularSize", style.get("curve", 90))
    else:
        s = _set_attr(s, "ArrowheadHead",
                      "Full" if (is_last and heads not in ("none", "start")) else "(none)")
        if is_first and heads in ("both", "start"):
            s = _set_attr(s, "ArrowheadTail", "Full")
        if kind == "retro":
            # same double-shaft mechanism, so the same hundredths unit. NOT measured:
            # the calibration export had no retro arrow. Verify before trusting.
            s = _set_attr(s, "ArrowShaftSpacing", style.get("shaftSpacing", RETRO_SHAFT_SPACING))

    nogo = style.get("noGo") or ("Cross" if kind == "noGoCross" else "Hash" if kind == "noGoHash" else None)
    if nogo and is_last:
        s = _set_attr(s, "NoGo", "Hash" if re.search("hash", str(nogo), re.I) else "Cross")
    if style.get("curve"):
        s = _set_attr(s, "AngularSize", style["curve"])
    line = style.get("line")
    if line and line != "solid":
        s = _set_attr(s, "LineType", line_type(line))
    if style.get("headSize"):
        s = _set_attr(s, "HeadSize", style["headSize"])
    if style.get("width"):
        s = _set_attr(s, "ArrowheadWidth", style["width"])

    # ---- attributes the round-trip loop verified as real (they survive an open+save in
    # ChemDraw AND change what is drawn). Each was a guess until it was measured. ----
    if style.get("color") is not None:
        s = _set_attr(s, "color", C.color_index(style["color"]))
    if style.get("dipole"):
        s = _set_attr(s, "Dipole", "yes")
    if style.get("hashSpacing") is not None:
        s = _set_attr(s, "HashSpacing", style["hashSpacing"])
    if style.get("lineWidth") is not None:
        s = _set_attr(s, "LineWidth", style["lineWidth"])
    if style.get("boldWidth") is not None:
        s = _set_attr(s, "BoldWidth", style["boldWidth"])
    if style.get("equilibriumRatio") is not None:
        # asymmetric equilibrium: one shaft drawn shorter than the other, which is how a
        # scheme shows WHICH WAY the equilibrium lies. Only the value 50 has been observed
        # (visibly asymmetric); the scale is not calibrated yet.
        s = _set_attr(s, "ArrowEquilibriumRatio", style["equilibriumRatio"])
    return s


def styled_arrow(aid: int, x0: float, y0: float, x1: float, y1: float, style: dict | None = None) -> str:
    """A single-segment styled arrow between two explicit points."""
    return apply_style(C.arrow_xy(aid, x0, y0, x1, y1), style, is_first=True, is_last=True)


# ---- route builders: return a list of segment strings (styled), head on last seg ----
def _reach(hw, hh, ux, uy):
    return min(abs(hw / ux) if ux != 0 else 1e9, abs(hh / uy) if uy != 0 else 1e9)


def _seg(aid, x0, y0, x1, y1, style, is_first, is_last):
    return apply_style(C.arrow_xy(aid, x0, y0, x1, y1), style, is_first=is_first, is_last=is_last)


def route(spec: dict, aid: int = 90001, margin: float = 16.0):
    """Build a (possibly multi-segment) styled route.
    spec: {kind:'straight'|'elbow_hv'|'elbow_vh'|'via',
           f:(cx,cy), fh:(hw,hh), t:(cx,cy), th:(hw,hh), via:[[x,y],...], style:{...}}
    Returns {'segs':[...], 'head':(hx0,hy0,hx1,hy1), 'pts':[(x,y), ...]}
    head = final segment endpoints; `pts` is the WHOLE run, corner points included,
    already trimmed to the node edges.

    `pts` exists so a caller can reserve what was actually drawn instead of
    re-deriving it. Reserving the bounding box of the endpoints for every
    routed kind is wrong: the bounding box of an L is mostly the empty
    corner -- the connector then "collides" with things it never touches.
    Handing back the run keeps the geometry computed once, here, which is the same
    reason `head` is returned rather than recomputed by the caller."""
    kind = spec.get("kind", "straight")
    style = spec.get("style")
    (fx, fy), (fhw, fhh) = spec["f"], spec["fh"]
    (tx, ty), (thw, thh) = spec["t"], spec["th"]

    if kind == "straight":
        dx, dy = tx - fx, ty - fy
        L = math.hypot(dx, dy) or 1.0
        ux, uy = dx / L, dy / L
        sR = _reach(fhw, fhh, ux, uy) + margin
        tR = L - (_reach(thw, thh, ux, uy) + margin)
        x0, y0, x1, y1 = fx + ux * sR, fy + uy * sR, fx + ux * tR, fy + uy * tR
        return {"segs": [_seg(aid, x0, y0, x1, y1, style, True, True)],
                "head": (x0, y0, x1, y1), "pts": [(x0, y0), (x1, y1)]}

    if kind == "elbow_hv":                    # horizontal from source, then vertical into target
        sgx = (1 if tx >= fx else -1)
        sgy = (1 if ty >= fy else -1)
        sx, sy = fx + sgx * (fhw + margin), fy
        cx, cy = tx, fy
        ex, ey = tx, ty - sgy * (thh + margin)
        return {"segs": [_seg(aid, sx, sy, cx, cy, style, True, False),
                         _seg(aid + 1, cx, cy, ex, ey, style, False, True)],
                "head": (cx, cy, ex, ey),
                "pts": [(sx, sy), (cx, cy), (ex, ey)]}

    if kind == "elbow_vh":                    # vertical from source, then horizontal into target
        sgx = (1 if tx >= fx else -1)
        sgy = (1 if ty >= fy else -1)
        sx, sy = fx, fy + sgy * (fhh + margin)
        cx, cy = fx, ty
        ex, ey = tx - sgx * (thw + margin), ty
        return {"segs": [_seg(aid, sx, sy, cx, cy, style, True, False),
                         _seg(aid + 1, cx, cy, ex, ey, style, False, True)],
                "head": (cx, cy, ex, ey),
                "pts": [(sx, sy), (cx, cy), (ex, ey)]}

    if kind == "via":
        pts = [[fx, fy]] + [list(p) for p in spec.get("via", [])] + [[tx, ty]]
        # trim first/last to node edges
        ax, ay = pts[0]; bx, by = pts[1]
        dx, dy = bx - ax, by - ay; L = math.hypot(dx, dy) or 1.0
        r = _reach(fhw, fhh, dx / L, dy / L) + margin
        pts[0] = [ax + dx / L * r, ay + dy / L * r]
        n = len(pts); ax, ay = pts[n - 1]; bx, by = pts[n - 2]
        dx, dy = ax - bx, ay - by; L = math.hypot(dx, dy) or 1.0
        r = _reach(thw, thh, dx / L, dy / L) + margin
        pts[n - 1] = [ax - dx / L * r, ay - dy / L * r]
        segs = []
        for i in range(len(pts) - 1):
            segs.append(_seg(aid + i, pts[i][0], pts[i][1], pts[i + 1][0], pts[i + 1][1],
                             style, i == 0, i == len(pts) - 2))
        n = len(pts)
        return {"segs": segs,
                "head": (pts[n - 2][0], pts[n - 2][1], pts[n - 1][0], pts[n - 1][1]),
                "pts": [tuple(q) for q in pts]}

    # fallback: straight
    return route({**spec, "kind": "straight"}, aid, margin)


# ---------------------------------------------------------------------------
# Composite arrows — one LOGICAL arrow may emit SEVERAL CDXML elements.
#
# apply_style() above handles "one <arrow>, attributes overlaid".
# emit_arrow() sits on top of it and is allowed to return a LIST, which is what
# every multi-element pattern needs: multi-step twins, attached arcs.
#
# NOMINAL_HALF is the vertical half-band cdxml_core.arrow() already writes into
# its BoundingBox. ink_extra() reports ink BEYOND that band, so a plain arrow
# returns all zeros and existing figures are byte-identical.
# ---------------------------------------------------------------------------

# ---- MEASURED from a ChemDraw round trip (a sheet of test arrows, ACS template,
# BondLength 14.4). Each number is the MAX single-side reach of the rendered ink from the
# shaft axis. Max reach, not half-height, because this extent model is symmetric: a
# one-sided barb (harpoon) measures half 1.40 but still reaches 2.50 on the barb side, so
# using half would under-reserve exactly where the ink is. Over-reserving costs whitespace;
# under-reserving collides.
#
# Reported relative to a PLAIN arrow rather than to the nominal band, because the caption
# gaps in layout_linear were tuned against a plain arrow -- "how much taller than plain is
# this style" is the quantity the layout actually needs, and it keeps plain output
# byte-identical.
NOMINAL_HALF = 2.00        # what cdxml_core.arrow() writes into its BoundingBox
PLAIN_REACH = 2.50         # measured reach of a plain arrow: the datum
SHAFT_REACH = 0.30         # headless arrow -- i.e. the ENTIRE nominal band is arrowhead

# reach from the axis, per style. dashed / bold / wavy / reverse / resonance / harpoon all
# measured exactly PLAIN_REACH, so they are absent here and contribute no extra.
#
# Two of those deserve an explanation, because "same extent" means different things:
#
#   bold  DOES work -- the shaft is visibly thicker. It just cannot change the extent: the
#         document's BoldWidth is 2.01 against a LineWidth of 0.60, so the shaft goes from
#         +/-0.30 to +/-1.005 and stays well inside the arrowhead's 2.50 reach. So extra=0
#         is right -- but not because bold is a no-op. For anything thicker than
#         bold, set LineWidth directly (LineWidth=2
#         measured a reach of 5.39) or BoldWidth, both of which the DTD allows per arrow.
#
#   wavy  really is a no-op on arrows. LineType="Wavy" SURVIVES the round trip -- ChemDraw
#         stores it -- but renders identically to plain, confirmed visually. There is no
#         wavy arrow available this way; the palette's wavy entries are separate presets.
REACH = {
    "noGoCross": 3.30, "noGoHash": 3.30,   # + the cross/hash overhangs each end
    "equilibrium": 4.15,                    # at its own default HeadSize 1500
    "retro": 5.55,                          # NOT 4.15 by analogy -- measured 1.4 more
    "hollowHead": 3.10,                     # hollow and angle heads measure the same
}
OVERHANG = {"noGoCross": 0.30, "noGoHash": 0.30, "equilibrium": 0.90,
            "retro": 0.60, "hollowHead": 0.60}

EQUIL_SHAFT_SPACING = 300  # hundredths of a unit => 3.00 units
RETRO_SHAFT_SPACING = 300  # same unit; survived the round trip and gave a distinct reach


# LineWidth thickens the shaft AND scales the head. Only two points measured:
# the document default 0.60 -> reach 2.50, and LineWidth 2.0 -> reach 5.39. Straight line
# through them, which is a 2-POINT FIT and wants a third point before it is trusted far
# outside that span.
LW_DEFAULT, LW_REACH_DEFAULT = 0.60, 2.50
LW_SLOPE = (5.39 - LW_REACH_DEFAULT) / (2.0 - LW_DEFAULT)


def line_reach(line_width: float) -> float:
    """Ink reach of a shaft drawn at `line_width`, from the two measured points."""
    return LW_REACH_DEFAULT + LW_SLOPE * (float(line_width) - LW_DEFAULT)


def head_reach(head_size: float) -> float:
    """Ink reach of a solid arrowhead, from HeadSize. Measured at 800 / 1000 / 1500 and
    exactly linear: 1.70 / 2.00 / 2.75 half-height, i.e. 2.20 / 2.50 / 3.25 reach. So the
    three-entry HEAD_PRESETS table is not needed for extents -- this formula covers it."""
    return head_size / 1000.0 * 1.5 + 1.0


def _repeat_offsets(spec: dict | None) -> list[tuple[float, float]]:
    """Per-copy (dx, dy) for a multi-step twin arrow, CENTERED on the axis so the
    group still reads as sitting on the row baseline (rather than hanging below
    it). spec['repeat'] = {'n':2, 'dx':10, 'dy':10}. No repeat -> [(0.0, 0.0)]."""
    rp = (spec or {}).get("repeat")
    if not rp:
        return [(0.0, 0.0)]
    n = max(1, int(rp.get("n", 2)))
    if n == 1:
        return [(0.0, 0.0)]
    dx, dy = float(rp.get("dx", 10.0)), float(rp.get("dy", 10.0))
    mid = (n - 1) / 2.0
    return [((i - mid) * dx, (i - mid) * dy) for i in range(n)]


def _resolve_arcs(spec: dict | None, cx: float, y: float, length: float,
                  aid: int = 1) -> tuple[list[str], list]:
    """Build every arc attached to this arrow. Returns (elements, bboxes).

    ONE resolver serves both emit_arrow and ink_extra: the bbox is read back out of the
    element cdxml_core already wrote, so what gets measured is exactly what gets drawn.
    Deriving the extent separately is how the same geometry ends up implemented twice.

    Each spec["arcs"] entry:
        at    0..1  position along the shaft where the arc touches it. FRACTIONAL, so it
                    survives the arrow being re-lengthened or its caption re-wrapped.
        off   [dx, dy]  the free end, relative to that anchor point.
        dir   "in"  (default) the arc ARRIVES at the shaft -- a reagent curving into the
                    flow. "out" it DEPARTS -- a by-product peeling off.
        bulge optional. Given, the arc is a plain bow through both points (arc_chord) and
                    the landing angle falls out. Omitted, the arc meets the shaft
                    TANGENTIALLY (arc_tangent), which is what reads as merging into it.
        head/tail   arrowhead ends, default head only.
    """
    out, boxes = [], []
    left = cx - length / 2.0
    for a in (spec or {}).get("arcs") or []:
        ax = left + float(a.get("at", 0.5)) * length
        dx, dy = [float(v) for v in a.get("off", (-40.0, -40.0))]
        if abs(dy) < 1e-6:
            raise ValueError("arcs[].off needs a non-zero dy: a free end level with the "
                             "shaft has no finite arc through it, use a straight arrow")
        fx, fy = ax + dx, y + dy
        kw = {"head": a.get("head", "Full"), "tail": a.get("tail")}
        if a.get("line"):
            kw["line"] = a["line"]
        going_out = a.get("dir", "in") == "out"
        # the tangent points along the direction of TRAVEL in x, so the arc never has to
        # double back on itself to satisfy the constraint
        t = (1.0 if dx >= 0 else -1.0) if going_out else (1.0 if dx <= 0 else -1.0)
        if a.get("bulge") is not None:
            p = ((ax, y), (fx, fy)) if going_out else ((fx, fy), (ax, y))
            el = C.arc_chord(aid, p[0][0], p[0][1], p[1][0], p[1][1],
                             float(a["bulge"]), **kw)
        elif going_out:
            el = C.arc_tangent(aid, ax, y, fx, fy, tangent=(t, 0.0), at="tail", **kw)
        else:
            el = C.arc_tangent(aid, fx, fy, ax, y, tangent=(t, 0.0), at="head", **kw)
        out.append(el)
        boxes.append(C.bbox_of(el))
        aid += 1
    return out, boxes


def ink_extra(spec: dict | None, length: float = 0.0) -> dict:
    """Ink reaching BEYOND cdxml_core.arrow()'s nominal +/-NOMINAL_HALF band.

    Returns TWO sets of numbers, because they answer different questions:

      up / down / left / right      the SHAFT's own ink. Captions clear this, and only
                                    this -- an attached arc can reach 46 units above the
                                    shaft, and pushing the caption above the arc's apex
                                    would be absurd when the two are at different x.
      box_up / box_down / ...       the TOTAL footprint including attached arcs. The page
                                    bbox, the cell advance and the row-spacing pre-pass
                                    use these, so an arc still cannot overlap a neighbour
                                    or the row above.

    All values >= 0."""
    spec = spec or {}
    st = spec.get("style") or {}
    kind = st.get("kind", "normal")

    # reach of ONE copy from its own axis, then converted to "extra beyond plain"
    reach, over = PLAIN_REACH, 0.0
    key = kind if kind in REACH else ("hollowHead" if st.get("headType") in
                                      ("hollow", "angle") else None)
    if st.get("noGo") and key is None:
        key = "noGoCross"
    if key:
        reach = max(reach, REACH[key])
        over = max(over, OVERHANG.get(key, 0.0))
    if st.get("headSize"):
        reach = max(reach, head_reach(float(st["headSize"])))
    if st.get("lineWidth") is not None:
        reach = max(reach, line_reach(st["lineWidth"]))
    half = NOMINAL_HALF + max(0.0, reach - PLAIN_REACH)

    offs = _repeat_offsets(spec)
    xs = [o[0] for o in offs]
    ext = {
        "up":    max(0.0, max(half - o[1] for o in offs) - NOMINAL_HALF),
        "down":  max(0.0, max(half + o[1] for o in offs) - NOMINAL_HALF),
        "left":  max(0.0, -min(xs)) + over,
        "right": max(0.0,  max(xs)) + over,
    }
    for k in ("up", "down", "left", "right"):
        ext["box_" + k] = ext[k]

    # attached arcs enlarge the FOOTPRINT only. Measured from the elements themselves, so
    # the reservation can never disagree with what gets drawn.
    if (spec or {}).get("arcs") and length:
        _, boxes = _resolve_arcs(spec, 0.0, 0.0, length)
        if boxes:
            span = length / 2.0
            ext["box_up"] = max(ext["box_up"], -min(b[1] for b in boxes) - NOMINAL_HALF)
            ext["box_down"] = max(ext["box_down"], max(b[3] for b in boxes) - NOMINAL_HALF)
            ext["box_left"] = max(ext["box_left"], -span - min(b[0] for b in boxes))
            ext["box_right"] = max(ext["box_right"], max(b[2] for b in boxes) - span)
    return {k: max(0.0, v) for k, v in ext.items()}


def emit_arrow(aid: int, cx: float, y: float, length: float,
               spec: dict | None = None) -> tuple[list[str], int]:
    """Every CDXML element for one logical arrow. Returns (elements, next_aid).
    `spec` is the recipe's arrow dict (style / repeat / later arcs+badge), NOT a
    bare style dict — pass spec.get('style') to apply_style yourself if you only
    have a style."""
    spec = spec or {}
    style = spec.get("style")
    if spec.get("reverse"):
        # RIGHT-TO-LEFT arrow: the shaft is still laid out left-to-right by the packer,
        # so reversing means moving the head to the tail end. heads="start" does exactly
        # that (ArrowheadHead none + ArrowheadTail Full) - this just gives it a name, so
        # a caller doesn't have to know the heads vocabulary to point an arrow backwards.
        style = {**(style or {}), "heads": "start"}
    out = []
    for dx, dy in _repeat_offsets(spec):
        s = C.arrow(aid, cx + dx, y + dy, length)
        if style:
            s = apply_style(s, style)
        out.append(s)
        aid += 1
    arcs, _ = _resolve_arcs(spec, cx, y, length, aid)
    out.extend(arcs)
    return out, aid + len(arcs)
