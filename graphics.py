"""Page-level graphics a caller places by hand: boxes, rules, blocks, dots.

The sibling of decorate.py, and the difference between them is the whole design:

    decorate   lives INSIDE a <fragment>, anchored by SMARTS to atoms and bonds.
               It moves when the structure moves, and reaches every figure tool
               for free because the fragment carries it.
    graphics   lives on the PAGE at coordinates. Nothing anchors it; nothing
               moves it. Every tool that wants it has to emit it.

WHY THE CALLER GIVES COORDINATES, which is not where this wants to end up. The
figure's own layers -- a section title block, a panel frame, a slot for a ligand
-- all need to know where a section or a cell ENDED UP, and that answer lives in
a spatial ledger the engines do not share (the six of them were grown one at a
time). Building the positioning half here would mean writing it
once per engine, which is the shape of the problem rather than a fix for it. So
this ships the half that has one home -- the objects, their styles, and a
coordinate interface -- and automatic placement waits for the shared ledger.

    Measured, not assumed: ChemDraw's own documents have no layout container
    either. `<group>` carries BoundingBox / id / Integral / Z and nothing else,
    and in Citric Acid Cycle 18 rectangles enclose zero fragments. There is
    nothing to copy.

THE COORDINATE FRAME is the figure's own content box, (0, 0) at its top-left,
y downwards -- NOT the page, whose origin depends on page_mode and the margin.
The engine offsets these by its content box before finalize, so a caller draws
against the box the return value reported on the previous call. Two passes, and
that is honest: nothing here can know how wide a grid came out.
"""

from __future__ import annotations

import cdxml_core as C

KINDS = ("rect", "line", "dot", "text", "arc", "poly")

# One spec's allowed keys, per kind. Rejected rather than ignored: a graphic that
# silently does not appear is indistinguishable from one nobody asked for, which
# is the failure this module exists to remove.
_KEYS = {
    "rect": frozenset({"kind", "at", "style", "color", "radius", "shadow"}),
    "line": frozenset({"kind", "at", "color", "width", "dash", "bold",
                       # A headed line is not a
                       # GraphicType="Line" at all -- see _straight_arrow.
                       "head", "head_type"}),
    "dot": frozenset({"kind", "at", "r", "color", "shaded"}),
    # FREE TEXT AT A COORDINATE. Among published table-of-contents graphics,
    # 38 of 105 introduction figures carry a selling-point block, a banner title
    # or a grey caption slab, and without this every one of them has to be
    # smuggled into some engine's own text slot. That costs something: bullets
    # that ride in as FOOTNOTES wrap on the footnote rule and can break
    # "* Aryl radical-mediated intermolecular HAT" across two lines -- one
    # selling point cut in half, no warning.
    #
    # NOTHING IS RENDERED HERE THAT IS NOT ALREADY RENDERED ELSEWHERE.
    # `C.styled_text` carries nested markup, colour, bold/italic, justification
    # and rotation, and `slot()` twenty lines down calls it for its labels. What
    # this adds is a way to ASK for one that is not a label of something else.
    "text": frozenset({"kind", "at", "text", "size", "color",
                       "bold", "italic", "align", "rot",
                       # A TEXT BLOCK'S OWN BACKGROUND SLAB. 30 of those 105
                       # introduction figures set their captions on a pale
                       # rectangle, and drawing it as a separate `rect` cannot
                       # work once the text is PLACED rather than given
                       # coordinates: the slab has to be emitted BEFORE the text
                       # to sit behind it (document order is what Z is stamped
                       # from), but its size is only known AFTER the text has
                       # been measured and placed. One spec breaks the cycle.
                       "fill", "fill_pad"}),
    # AN ELLIPTICAL ARC ARROW -- the ring that circles a catalyst badge (21 of
    # the 105 draw one). `GraphicType="Arc"` occurs 196 times in a reference set
    # of ChemDraw files, so this is a spelling that was read off real files rather
    # than reasoned out.
    "arc": frozenset({"kind", "at", "r", "ry", "t0", "sweep",
                      "head", "tail", "head_type", "color"}),
    # AN ARBITRARY POLYGON, OUTLINE OR FILLED. The shape published figures keep
    # asking for: a fat hollow block arrow, a lightning bolt, a pie sector, a
    # tilted band, an LED. Without it every one of them has to be assembled by
    # the caller out of 7..10 separate `line` graphics, and the filled half of
    # them cannot be drawn at all.
    #
    # `points` is the caller's corners. `fill` and `closed` pick the CurveType
    # bitmask (1 = Closed, 128 = Filled, both = 129); see the comment on the
    # emit branch for the part that had to be MEASURED.
    "poly": frozenset({"kind", "points", "fill", "closed", "color", "width"}),
}

DEFAULT_TEXT_SIZE = 10          # == label_size, what a titled thing already uses
TEXT_ALIGNS = ("left", "center", "right")

# A STRAIGHT ARROW IS AN `<arrow>`, NOT A `<curve>`. Measured, not reasoned.
#
# A `<curve>` with four collinear Bezier control points looks like the obvious
# spelling, since `C.curve` already supports every arrowhead variant, and its XML
# is correct by every check available without ChemDraw. Opened in ChemDraw, a
# `<curve>` carrying an arrowhead draws THE ARROWHEAD ALONE, with no line -- and
# that is true of the bowed shape too, so it is not about the straight case or
# the zero-area bounding box.
#
# A test file put five encodings of the same two points on one page. Only the
# native `<arrow>` with Head3D / Tail3D drew a line at all:
#
#   1  <curve> straight, head, no HeadSize    arrowhead only, no line
#   2  <curve> straight, head, HeadSize=800   arrowhead only, no line
#   3  <curve> BOWED, head                    arrowhead only, no line
#   4  <curve> straight, DIAGONAL             arrowhead only, no line
#   5  native <arrow>, Head3D/Tail3D          DRAWS
#
# The arrow subsystem is the UNDOCUMENTED part of the CDXML spec, and for
# undocumented objects a hand-checked file is the source of truth. The curve route
# relies on a documented-object argument and so skips that check.
HEAD_TRIPLE = C.HEAD_PRESETS["straight"]     # (1000, 875, 250), measured


def straight_arrow(aid: int, p0, p1, *, head: str = "Full",
                   head_type: str = "Solid", triple=None, color=None) -> str:
    """One `<arrow>` from p0 to p1, in any direction.

    Written here rather than in cdxml_core because that module is not to be
    touched. The attribute set is `cdxml_core.arrow()`'s -- which is horizontal-only
    -- with Head3D / Tail3D / Center3D generalised to two arbitrary points; the
    head triple comes from `C.HEAD_PRESETS`, which was measured off a hand-drawn
    export, so no number here is invented. Default is the "straight" preset,
    because that is what ChemDraw itself uses for a plain straight arrow.
    """
    # "NO HEAD" ON AN <arrow> IS SPELLED "(none)", not "None".
    #
    # Read off real ChemDraw files, which are the source of truth for a spelling:
    # across 545 hand-drawn files ChemDraw writes ArrowheadHead="(none)" 82 times
    # on an <arrow> and "None" not once. ("None" IS in cdxml_core.CURVE_HEADS, but
    # that is the <curve> element's enumeration -- a different element.) With
    # "None", the pieces of a `crossing: "bridge"` line come out with an
    # arrowhead each.
    if head in (None, "None", "none"):
        head = "(none)"
    (x0, y0), (x1, y1) = (float(p0[0]), float(p0[1])), (float(p1[0]), float(p1[1]))
    hs, hc, hw = triple or HEAD_TRIPLE
    cx, cy = (x0 + x1) / 2.0, (y0 + y1) / 2.0
    pad = hw / 100.0
    bx0, by0 = min(x0, x1) - pad, min(y0, y1) - pad
    bx1, by1 = max(x0, x1) + pad, max(y0, y1) + pad
    col = "" if color is None else ' color="%d"' % C.color_index(color)
    return (f'<arrow id="{aid}"{col}'
            f' BoundingBox="{bx0:.2f} {by0:.2f} {bx1:.2f} {by1:.2f}"'
            f' FillType="None" ArrowheadHead="{head}" ArrowheadType="{head_type}"'
            f' HeadSize="{hs}" ArrowheadCenterSize="{hc}" ArrowheadWidth="{hw}"'
            f' Head3D="{x1:.2f} {y1:.2f} 0" Tail3D="{x0:.2f} {y0:.2f} 0"'
            f' Center3D="{cx:.2f} {cy:.2f} 0"/>')


def straight_points(p0, p1, n: int = 24):
    """The segment p0 -> p1 sampled into n+1 points, for a ledger reservation.

    The counterpart of `C.arc_points` for a straight arrow: an `<arrow>` has no
    CurvePoints to read back, so the ink is generated from the same two endpoints
    the object was built from -- one source, so the reservation cannot describe a
    different line from the one drawn.
    """
    (x0, y0), (x1, y1) = p0, p1
    return [(x0 + (x1 - x0) * k / n, y0 + (y1 - y0) * k / n)
            for k in range(n + 1)]


DEFAULT_LINE_WIDTH = 1.0
DEFAULT_DOT_R = 3.0


def _at(spec, n, where):
    at = spec.get("at")
    if not isinstance(at, (list, tuple)) or len(at) != n:
        raise ValueError(
            "%s: 'at' is %d numbers%s, got %r"
            % (where, n, " (x0 y0 x1 y1)" if n == 4 else " (x y)", at))
    try:
        return [float(v) for v in at]
    except (TypeError, ValueError):
        raise ValueError("%s: 'at' holds a non-number: %r" % (where, at))


def _line_type(spec):
    """dash / bold -> the LineType word ChemDraw writes, or None for plain.

    Two booleans rather than a `line_type` string because the combination is real
    -- ChemDraw writes "Dashed Bold" -- and asking a caller to know that spelling
    is asking them to know the DTD is wrong about it.
    """
    dash, bold = bool(spec.get("dash")), bool(spec.get("bold"))
    if dash and bold:
        return "Dashed Bold"
    if dash:
        return "Dashed"
    if bold:
        return "Bold"
    return None


def emit(specs, gid: int = 98000, origin=(0.0, 0.0), tid: int = 99500):
    """(xml, [x0, y0, x1, y1]) for a caller's graphics list, or ("", None).

    `tid` numbers the `<t>` objects a "text" graphic emits, and it is a SEPARATE
    range from `gid` because a `<t>` and a `<graphic>` sharing an id is a
    duplicate id in the finished file -- which compose's write() guard reports
    and ChemDraw itself resolves unpredictably. 99500 sits clear of slot()'s own
    tid block (98800 + 40 per slot).

    `origin` is where the caller's (0, 0) lands in the engine's own layout space
    -- the content box's top-left corner. Every spec is shifted by it here, once,
    so no tool has to remember to.

    The returned box is the extent of what was drawn, for the caller to union
    into the figure's content box. Without that a frame drawn round the edge of a
    figure would fall outside the page: finalize sizes the page from the box it
    is handed, and it is handed the ENGINE's, which knows nothing about these.

    COLOURS MUST ALREADY BE INDICES. An [r, g, b] belongs in the document's
    colortable, which is written into the header before any body exists, so the
    tool layer runs _extend_palette over these specs first and rewrites `color`
    in place -- the same order of operations a highlight colour goes through.
    Passing raw rgb here raises out of color_index rather than guessing an index
    against a palette that never received the colour.
    """
    if not specs:
        return "", None
    ox, oy = origin
    out, xs, ys = [], [], []

    for i, spec in enumerate(specs):
        where = "graphics[%d]" % i
        if not isinstance(spec, dict):
            raise ValueError("%s: each graphic is a dict, got %s"
                             % (where, type(spec).__name__))
        kind = spec.get("kind")
        if kind not in KINDS:
            raise ValueError("%s: 'kind' is one of %s, got %r"
                             % (where, ", ".join(KINDS), kind))
        unknown = set(spec) - _KEYS[kind]
        if unknown:
            raise ValueError(
                "%s: %s is not a key of a %r graphic. Known: %s"
                % (where, ", ".join(sorted(unknown)), kind,
                   ", ".join(sorted(_KEYS[kind]))))

        if kind == "rect":
            x0, y0, x1, y1 = _at(spec, 4, where)
            style = spec.get("style", "Plain")
            if style not in C.RECTANGLE_TYPES:
                raise ValueError(
                    "%s: style %r has no precedent in any ChemDraw file. Known: %s"
                    % (where, style, ", ".join(C.RECTANGLE_TYPES)))
            out.append(C.rect(gid, x0 + ox, y0 + oy, x1 + ox, y1 + oy,
                              kind=style, color=spec.get("color"),
                              radius=int(spec.get("radius", C.CORNER_RADIUS)),
                              shadow=int(spec.get("shadow", C.SHADOW_SIZE))))
            gid += 1
            xs += [x0 + ox, x1 + ox]
            ys += [y0 + oy, y1 + oy]

        elif kind == "line":
            x0, y0, x1, y1 = _at(spec, 4, where)
            head = spec.get("head")
            if head:
                # A HEADED line is a native <arrow>: a Line graphic cannot carry an
                # arrowhead, and a <curve> that carries one draws the head with no
                # line (measured -- see straight_arrow). Only when asked, so a
                # headless `line` emits unchanged bytes, which the scope and cycle
                # goldens pin character by character.
                if head not in C.CURVE_HEADS:
                    raise ValueError(
                        "%s: 'head' is one of %s, got %r"
                        % (where, ", ".join(C.CURVE_HEADS), head))
                ht = spec.get("head_type", "Solid")
                if ht not in C.CURVE_HEAD_TYPES:
                    raise ValueError(
                        "%s: 'head_type' is one of %s, got %r"
                        % (where, ", ".join(C.CURVE_HEAD_TYPES), ht))
                if spec.get("dash") or spec.get("bold") or "width" in spec:
                    # Said rather than silently dropped: an <arrow> carries no
                    # LineWidth or LineType, so these would vanish.
                    raise ValueError(
                        "%s: 'head' turns this line into an arrow object, which "
                        "carries no 'width' / 'dash' / 'bold'. Drop the head, or "
                        "drop those." % where)
                out.append(straight_arrow(gid, (x0 + ox, y0 + oy),
                                          (x1 + ox, y1 + oy),
                                          head=head, head_type=ht,
                                          color=spec.get("color")))
            else:
                out.append(C._line(gid, x0 + ox, y0 + oy, x1 + ox, y1 + oy,
                                   lw=float(spec.get("width", DEFAULT_LINE_WIDTH)),
                                   color=spec.get("color"),
                                   line_type=_line_type(spec)))
            gid += 1
            xs += [x0 + ox, x1 + ox]
            ys += [y0 + oy, y1 + oy]

        elif kind == "dot":
            x, y = _at(spec, 2, where)
            r = float(spec.get("r", DEFAULT_DOT_R))
            # disc(), not a second circle emitter. A hand-drawn reference draws its
            # endpoint dots as GraphicType="Oval" and ChemDraw writes both, but a
            # filled circle already has one home here and it is the one measured
            # against a hand-drawn file.
            out.append(C.disc(gid, x + ox, y + oy, r,
                              spec.get("color", "black"),
                              shaded=bool(spec.get("shaded"))))
            gid += 1
            xs += [x + ox - r, x + ox + r]
            ys += [y + oy - r, y + oy + r]

        elif kind == "text":
            # deferred: the text metrics live in the layout module, and importing
            # it at module scope would make graphics depend on a layout engine.
            from layout_linear import block_box, ASCENT_EM
            x, y = _at(spec, 2, where)
            lines = spec.get("text")
            if isinstance(lines, str):
                lines = [lines]
            if (not isinstance(lines, (list, tuple)) or not lines
                    or not all(isinstance(t, str) for t in lines)):
                raise ValueError(
                    "%s: 'text' is a string or a list of strings (one per line), "
                    "got %r" % (where, spec.get("text")))
            align = spec.get("align", "left")
            if align not in TEXT_ALIGNS:
                raise ValueError("%s: 'align' is one of %s, got %r"
                                 % (where, ", ".join(TEXT_ALIGNS), align))
            size = float(spec.get("size", DEFAULT_TEXT_SIZE))
            # 9 AND NOT 9.0 IN THE FILE. A reference set of ChemDraw files writes
            # size="9" 26103 times and a decimal only where the size really is
            # fractional (size="7.5", 22 times) -- so an integral size goes in as
            # an int, and this module never produces a spelling no ChemDraw file
            # uses.
            size_out = int(size) if float(size).is_integer() else size
            bold, italic = bool(spec.get("bold")), bool(spec.get("italic"))
            w, h = block_box(lines, size, bold=bold)
            # `at` IS THE TOP-LEFT OF THE BLOCK, whatever `align` says, so it
            # reads the same as rect / line / dot and a caller can stack a text
            # block under a rect by adding the rect's height. `align` justifies
            # the lines INSIDE that block; styled_text anchors Left at the left
            # edge, Center at the middle and Right at the right, so the anchor
            # is derived here rather than pushed onto the caller.
            anchor = x + {"left": 0.0, "center": w / 2.0, "right": w}[align]
            # + ASCENT_EM * size: styled_text's `p` is the FIRST LINE'S BASELINE,
            # not the top of the block (slot() does the same conversion for its
            # labels). Without it every text graphic would ride one ascender high.
            fill = spec.get("fill")
            if fill is not None:
                fp = float(spec.get("fill_pad", 5.0))
                out.append(C.rect(gid, x + ox - fp, y + oy - fp,
                                  x + ox + w + fp, y + oy + h + fp,
                                  kind="Filled", color=fill))
                gid += 1
                xs += [x + ox - fp, x + ox + w + fp]
                ys += [y + oy - fp, y + oy + h + fp]
            out.append(C.styled_text(
                tid, anchor + ox, y + oy + ASCENT_EM * size, list(lines),
                size=size_out, color=spec.get("color"), bold=bold, italic=italic,
                align=align, rot=float(spec.get("rot", 0.0))))
            tid += 1
            xs += [x + ox, x + ox + w]
            ys += [y + oy, y + oy + h]

        # NAMED RATHER THAN LEFT AS THE CATCH-ALL. An `else: # arc` is fine
        # while arc is the last kind and stops being fine the moment another one
        # is added -- `kind` is checked against KINDS at the top of the loop, so
        # nothing can reach here unnamed anyway.
        elif kind == "arc":
            cx, cy = _at(spec, 2, where)
            if "r" not in spec:
                raise ValueError("%s: an arc needs 'r' (its radius)" % where)
            rx = float(spec["r"])
            ry = float(spec.get("ry", rx))
            head = spec.get("head", "Full")
            if head is not None and head not in C.CURVE_HEADS:
                raise ValueError("%s: 'head' is one of %s, got %r"
                                 % (where, ", ".join(C.CURVE_HEADS), head))
            tail = spec.get("tail")
            if tail is not None and tail not in C.CURVE_HEADS:
                raise ValueError("%s: 'tail' is one of %s, got %r"
                                 % (where, ", ".join(C.CURVE_HEADS), tail))
            ht = spec.get("head_type", "Solid")
            if ht not in C.CURVE_HEAD_TYPES:
                raise ValueError("%s: 'head_type' is one of %s, got %r"
                                 % (where, ", ".join(C.CURVE_HEAD_TYPES), ht))
            el = C.arc_arrow(gid, cx + ox, cy + oy, rx, ry,
                             float(spec.get("t0", 0.0)),
                             float(spec.get("sweep", 270.0)),
                             head=(None if head == "None" else head),
                             tail=(None if tail == "None" else tail),
                             head_type=ht)
            # arc_arrow predates colour and cdxml_core is not to be touched for
            # it; recolor() twenty lines down is where that attribute is set.
            out.append(recolor(el, spec.get("color")))
            gid += 1
            # MEASURED THROUGH arc_points, not from cx +/- r. A 90-degree sweep
            # occupies one quadrant, and boxing it as the whole circle would push
            # the figure's content box out by up to 2r on two sides -- which
            # finalize turns into blank page. arc_points is the canonical reader
            # for exactly this geometry and it takes the string back.
            pts = C.arc_points(el)
            if pts:
                xs += [p[0] for p in pts]
                ys += [p[1] for p in pts]
            else:                       # sweep 0: nothing drawn, but still placed
                xs.append(cx + ox)
                ys.append(cy + oy)

        elif kind == "poly":
            pts = spec.get("points")
            if not (isinstance(pts, (list, tuple)) and len(pts) >= 2
                    and all(isinstance(p, (list, tuple)) and len(p) == 2
                            for p in pts)):
                raise ValueError(
                    "%s: a 'poly' takes 'points' = [(x, y), ...] with at least "
                    "two corners; got %r" % (where, pts))
            closed = bool(spec.get("closed", True))
            filled = bool(spec.get("fill", False))
            if filled and not closed:
                raise ValueError(
                    "%s: a filled polygon has to be closed. CurveType 128 without 1 fills "
                    "the shape but leaves the outline open, so the last edge is drawn "
                    "straight across. Drop 'closed': False, or drop 'fill'."
                    % where)
            # CurveType is a BITMASK: 1 = Closed, 128 = Filled. 129 is a filled
            # polygon; 1 alone is the outline that otherwise takes seven `line`
            # graphics; 0 is an open polyline.
            ctype = (1 if closed else 0) | (128 if filled else 0)
            # THE CONTROL POINTS GO ON THE ENDPOINTS, and this is the part that
            # had to be measured rather than reasoned out. CurvePoints is a chain
            # of cubic Beziers, so the obvious spelling for a straight edge is to
            # put the two controls a third and two thirds of the way along it --
            # each segment is then dead straight and the tangents at a corner
            # differ, which ought to be a corner. It is NOT: ChemDraw draws that
            # with every corner ROUNDED (measured in a ChemDraw test file).
            # Putting both controls ON the endpoints gives edges that are straight
            # AND corners that are sharp, identical to the seven-line outline. A
            # tripled knot works too; this is the simpler statement, and
            # it generalises -- move a control off the endpoint and that one edge
            # bows, which is what a rounded banner will want later.
            ring = list(pts) + ([pts[0]] if closed else [])
            ctrl = [(ring[0][0] + ox, ring[0][1] + oy)]
            for _a, _b in zip(ring, ring[1:]):
                ctrl += [(_a[0] + ox, _a[1] + oy), (_b[0] + ox, _b[1] + oy),
                         (_b[0] + ox, _b[1] + oy)]
            out.append(C.curve(gid, ctrl, color=spec.get("color"),
                               curve_type=ctype,
                               line_width=spec.get("width")))
            gid += 1
            xs += [p[0] + ox for p in pts]
            ys += [p[1] + oy for p in pts]

    return "\n".join(out), [min(xs), min(ys), max(xs), max(ys)]


# --- a side slot -- structures parked beside the reaction, in a frame ---------
# The ligand box a scheme puts next to its arrow (L8, Xantphos, 4CzIPN), the two
# separate boxes a paper draws for (R)-CPA and L43, the shared frame holding
# ChCl / EG / Urea each under its own name -- a recurring item in published
# schemes.
#
# WHAT MAKES THIS THE DOWNGRADED VERSION, and where the line is: the caller says
# WHERE. Working out where a slot should go means knowing where the scheme ended
# and how much room is left beside it, which is the shared spatial ledger the six
# engines do not have -- a separate problem, and not one to solve twice.
#
# THE LINE IS ENFORCED BY THE SIGNATURE, not by a comment. slot() takes its own
# corner and its own content and NOTHING ELSE: no page, no figure box, no sibling
# positions. The tempting extension is to let a slot "just check" what is beside
# it, and that one line is the whole shared-ledger problem arriving unannounced.
# It cannot be written without changing this signature, which is the point.
DEFAULT_SLOT_GAP = 14.0
DEFAULT_SLOT_PAD = 10.0
DEFAULT_SLOT_LABEL_SIZE = 9
DEFAULT_SLOT_LABEL_GAP = 4.0

ARRANGEMENTS = ("row", "column", "grid")


def slot(at, items, *, frame=None, fill=None, arrange="row", cols=None,
         gap=DEFAULT_SLOT_GAP, pad=DEFAULT_SLOT_PAD,
         label_size=DEFAULT_SLOT_LABEL_SIZE, label_gap=DEFAULT_SLOT_LABEL_GAP,
         label_color=None, gid=98500, tid=98600, base=980000):
    """(xml, [x0, y0, x1, y1]) for one framed group of structures at `at`.

    `at` is the slot's top-left corner, `items` is [(frag, bbox, label_or_None)]
    ALREADY BUILT -- this module does not know what a SMILES is, and keeping it
    that way is why the tool layer builds them: they go through _register, so a
    slot inherits the duplicate-id check and the empty-id fallback for free.

    Each structure's name goes UNDERNEATH it, where a product grid puts a
    compound number, so a reader learns one rule rather than two. `label_color`
    for green names, as a deep-eutectic-solvent figure uses.

    frame=None DRAWS NO FRAME, and that is a real case rather than a degenerate
    one -- a slot can be separated by a vertical dashed rule instead, which
    is a `line` graphic beside a frameless slot. A white frame and no frame look
    identical on paper, so the test for this asserts there is no rect object.
    """
    if arrange not in ARRANGEMENTS:
        raise ValueError("slot arrange is one of %s, got %r"
                         % (", ".join(ARRANGEMENTS), arrange))
    if not items:
        raise ValueError("a slot with no structures in it has nothing to frame")
    if frame is not None and frame not in C.RECTANGLE_TYPES:
        raise ValueError("slot frame %r has no precedent. Known: %s"
                         % (frame, ", ".join(C.RECTANGLE_TYPES)))

    # measure: each cell is its structure, plus its name under it
    cells = []
    for item in items:
        # 3-tuple or 4-tuple: the 4th is the structure's decoration, already
        # built. Accepting both keeps every 3-tuple caller working unchanged.
        frag, bb, label = item[0], item[1], item[2]
        deco = item[3] if len(item) > 3 else None
        w, h = bb[2] - bb[0], bb[3] - bb[1]
        lw = lh = 0.0
        if label:
            from layout_linear import block_box       # deferred: text metrics only
            lw, lh = block_box([str(label)], label_size)
            h += label_gap + lh
            w = max(w, lw)
        cells.append([frag, bb, label, w, h, deco])

    n = len(cells)
    ncols = n if arrange == "row" else 1 if arrange == "column" else max(1, int(cols or 2))
    rows = [cells[i:i + ncols] for i in range(0, n, ncols)]
    col_w = [max((r[j][3] for r in rows if len(r) > j), default=0.0)
             for j in range(ncols)]
    row_h = [max(c[4] for c in r) for r in rows]

    x0, y0 = float(at[0]) + pad, float(at[1]) + pad
    out, deco_out, cx_end, cy = [], [], x0, y0
    for r, cells_in_row in enumerate(rows):
        cx = x0
        for j, (frag, bb, label, w, h, deco) in enumerate(cells_in_row):
            sw, sh = bb[2] - bb[0], bb[3] - bb[1]
            sx = cx + (col_w[j] - sw) / 2.0             # centred in its column
            _fx, _dx = C.place_with(frag, deco, base, sx, cy)
            out.append(_fx)
            if _dx:
                deco_out.append(_dx)
            base += 1000
            if label:
                from layout_linear import block_box, ASCENT_EM
                lw, lh = block_box([str(label)], label_size)
                ly = cy + sh + label_gap
                out.append(C.styled_text(
                    tid, cx + col_w[j] / 2.0, ly + ASCENT_EM * label_size,
                    [str(label)], size=label_size, align="center",
                    color=label_color))
                tid += 1
            cx += col_w[j] + gap
        cx_end = max(cx_end, cx - gap)
        cy += row_h[r] + gap
    inner = [x0, y0, cx_end, cy - gap]

    box = [inner[0] - pad, inner[1] - pad, inner[2] + pad, inner[3] + pad]
    head = []
    if fill is not None:
        head.append(C.rect(gid, *box, kind="Filled", color=fill)); gid += 1
    if frame is not None:
        head.append(C.rect(gid, *box, kind=frame)); gid += 1
    # fill, frame, then the halos, then the structures. The slot's own fill is
    # BEHIND a halo for the same reason it is behind a structure -- and the halo
    # is behind the structure because it is opaque (measured on ChemDraw
    # output). _add_graphics re-stamps Z over the whole
    # document afterwards, which is what makes this document order bind.
    return "\n".join(head + deco_out + out), box


def union(box, other):
    """The smallest box covering both. `other` may be None."""
    if other is None:
        return box
    if box is None:
        return list(other)
    return [min(box[0], other[0]), min(box[1], other[1]),
            max(box[2], other[2]), max(box[3], other[3])]


# --------------------------------------------------------------------------
# A bowed edge, and colour on an arc.
# --------------------------------------------------------------------------

def recolor(xml: str, color) -> str:
    """Return `xml` with color="<index>" set on its OPENING tag.

    `cdxml_core.arc_arrow` takes no colour -- it predates anything that needed one
    -- and that module is not to be touched, so the attribute is set here instead.
    The spelling is `cdxml_core.color_index`'s and `straight_arrow`'s: a colortable
    INDEX, never an [r, g, b], because rgb has to be appended to the document's
    table before it has an index and only the tool layer holds the header.

    `color=None` is a no-op, so a call site can pass whatever the config holds
    without branching -- which is what keeps an unmuted figure byte-identical.
    """
    if color is None:
        return xml
    idx = C.color_index(color)
    end = xml.find(">")
    head, rest = xml[:end], xml[end:]
    if ' color="' in head:
        import re as _re
        head = _re.sub(r'\bcolor="[^"]*"', 'color="%d"' % idx, head)
    else:
        head = head.replace("<arrow", '<arrow color="%d"' % idx, 1) \
            if head.startswith("<arrow") else head + ' color="%d"' % idx
    return head + rest


def bowed_arrow(aid: int, p0, p1, bulge: float, *, head="Full",
                head_type="Solid", color=None, **kw) -> str:
    """One `<arrow>` bowing from p0 to p1 by `bulge` x the chord length.

    Thin wrapper on `cdxml_core.arc_chord`, which already solves the circle through
    two points with a given sagitta -- the form a layout wants, since it names the
    two points the edge has to touch and lets the radius fall out. Here only to add
    colour and to keep the head preset the same one a straight edge uses, so a
    bowed parallel edge and its straight twin carry identical arrowheads.
    """
    hs, hc, hw = HEAD_TRIPLE
    xml = C.arc_chord(aid, float(p0[0]), float(p0[1]), float(p1[0]), float(p1[1]),
                      bulge=float(bulge), head=(head or None), head_type=head_type,
                      head_size=hs, head_center=hc, head_width=hw, **kw)
    return recolor(xml, color)


def bowed_points(p0, p1, bulge: float, n: int = 24):
    """The bowed edge sampled into n+1 points, for a ledger reservation.

    Decoded from the element `bowed_arrow` writes, via `cdxml_core.arc_points` --
    THE canonical reader for that geometry. Deriving the polyline independently
    here would be the same convention implemented twice, and the two would drift.
    """
    return C.arc_points(bowed_arrow(0, p0, p1, bulge, head=None), n)
