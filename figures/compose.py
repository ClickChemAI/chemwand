"""compose.py - multi-engine composition prototype. **Not one line of engine code
changed.**

Each engine runs one complete call as usual (headless, producing a complete CDXML); the
composer works only at the finished-document level: split <page>, renumber, translate,
union, merge tables, put back in the positive quadrant.

This is exactly panels.py's nesting model moved to the document level: each engine
measures its own bbox in its own ledger, and the composer learns only "this block takes
this much room". Engines cannot see each other's ink, so the composer adds its own union
collision check (report["collisions"]).

    c = Composer()
    c.cycle("ring", nodes=[...], steps=[...])
    c.linear("top",  items=[...])
    c.above("top", "ring", gap=26)          # stack by the actual extent of the ink
    c.join("tail", to="ring", at="II")      # shared anchor: one species is drawn only once
    rep = c.write("~/chemwand_figures/x151.cdxml")
"""
import math
import os
import re
import sys
import contextlib
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import cdxml_core as CORE       # noqa: E402
import cdxdoc as CD             # noqa: E402
import graphics as GX           # noqa: E402
import beautify as BT           # noqa: E402

# id range for cross-part furniture: above any id a single part uses, and below
# cross-part arrows (4.9M), and 4.9M itself is below part 1's offset of 5M. The three do
# not overlap.
FURNITURE_GID = 4_800_000
FURNITURE_TID = 4_850_000

# id range given to each part. The largest auto id in an engine is _add_graphics's
# 970000+si*20000; even 20 slots is only 1.37M, so a 5M stride leaves an order of
# magnitude of headroom.
ID_STRIDE = 5_000_000

# renumber() covers 6 scalars + 4 id lists; it does not cover these, and <scheme>/<step>
# rely entirely on them. Miss them and reaction steps point at another part's objects,
# while the file still opens fine.
_EXTRA_ID_ATTRS = ("ReactionStepReactants", "ReactionStepProducts", "ReactionStepArrows",
                   "ReactionStepObjectsAboveArrow", "ReactionStepObjectsBelowArrow",
                   "ReactionStepAtomMap", "BondOrdering", "BracketedObjectIDs",
                   "GraphicID", "SupersededBy")
# AtomID is deliberately not in here: it is the chemical atom number (1,2,3...), not an
# object id.

_COLOR_ATTRS = ("color", "bgcolor", "bcolor", "facecolor")


def _renumber_extra(body: str, base: int) -> str:
    for attr in _EXTRA_ID_ATTRS:
        body = re.sub(r'\b' + attr + r'="([\d\s]+)"',
                      lambda m: '%s="%s"' % (attr, " ".join(
                          str(base + int(t)) if int(t) else "0" for t in m.group(1).split())),
                      body)
    return body


def _shift_z(body: str, base: int) -> str:
    return re.sub(r'\bZ="(\d+)"', lambda m: 'Z="%d"' % (base + int(m.group(1))), body)


def _table(header: str, tag: str) -> list:
    m = re.search(r'<%s>(.*?)</%s>' % (tag, tag), header, re.S)
    return re.findall(r'<\w+\b[^>]*/>', m.group(1)) if m else []


# ---------------------------------------------------------------- small geometry helpers

def _ovl(a, b):
    """Overlap (dx, dy) of two bboxes; non-positive means no overlap."""
    return (min(a[2], b[2]) - max(a[0], b[0]), min(a[3], b[3]) - max(a[1], b[1]))


def _seg_seg(p, q) -> bool:
    def d(ax, ay, bx, by, cx, cy):
        return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
    d1 = d(p[0], p[1], p[2], p[3], q[0], q[1])
    d2 = d(p[0], p[1], p[2], p[3], q[2], q[3])
    d3 = d(q[0], q[1], q[2], q[3], p[0], p[1])
    d4 = d(q[0], q[1], q[2], q[3], p[2], p[3])
    return ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0))


def _seg_box(s, b) -> bool:
    if max(s[0], s[2]) < b[0] or min(s[0], s[2]) > b[2]:
        return False
    if max(s[1], s[3]) < b[1] or min(s[1], s[3]) > b[3]:
        return False
    if b[0] <= s[0] <= b[2] and b[1] <= s[1] <= b[3]:
        return True
    if b[0] <= s[2] <= b[2] and b[1] <= s[3] <= b[3]:
        return True
    e = [(b[0], b[1], b[2], b[1]), (b[2], b[1], b[2], b[3]),
         (b[2], b[3], b[0], b[3]), (b[0], b[3], b[0], b[1])]
    return any(_seg_seg(s, t) for t in e)


_FIND_RE = re.compile(r"^\s+(graphics\[|slots\[|\S+ covers )", re.M)


def _findings(report: dict) -> list:
    """The lines of a part report that are real "findings" (covers / out of bounds
    etc.), without the long explanation.

    This family is graphics **inside a part** covering a structure ("graphics[5] (rect)
    covers node:IntB by 7.0 pt"). The composer's covered() deliberately ignores the
    inside of a part - layer order within one source is that engine's own business - so
    this channel cannot be removed; remove it and nobody reports #010's electrode.
    """
    out = []
    for n in (report.get("notes") or []):
        for line in str(n).splitlines():
            if _FIND_RE.match(line) or " covers " in line:
                out.append(line.strip())
    return out


# The engine emits 28 kinds of NOTE, and _findings let not a single headline through
# (measured: 90 figure scripts, 712 parts, 1310 notice lines; the filter passed 201
# lines, **all** of them the graphics[N] covers cell:... family above). So what survived
# is exactly the family that should stay silent, and not one line of the family that
# really should be surfaced got through - 11 of them were "no anchor matched at all, the
# decoration was simply not drawn".
#
# This table picks one family: **you asked for it and it did not happen**, or **it
# quietly did something else**. The rest are not forwarded, each for its own reason:
#   "carry placeholder labels"    395 lines / 72 figures - pure information; every
#                                                  Markush figure has it
#   "reaches outside the page"    103 lines / 38 figures - the part's own page box;
#                                                  build() reframes by the union at the end,
#                                                  so it is discarded
#   "graphic(s) land on top"       75 lines / 27 figures - covered() answers this
#                                                  more precisely by draw order
#   "came through undeformed"       2 lines        - information, and good news
#
# The picking relies on wording, so test_every_forwarded_phrase_is_still_in_the_engine
# pins each one to the engine source: change the wording and it goes red first, and can
# say which one.
_DRIFT_PHRASES = (
    "highlight anchor(s) matched nothing",              # 11 lines / 5 figures (this round)
    "each match names a GROUP rather than one atom",    # 4 lines / 3 figures: one spec drew several markers
    "'fill' anchor(s) matched something with no ring in it",   # 2 lines: silently downgraded to color
    "did not match the depiction template",             # 3 lines: template did not match, each drawn on its own
    "came out of the constrained layout with uneven bonds",    # 2 lines: template was discarded
    "carry no group in the abbreviation table",         # 4 lines: abbreviate was on but had no effect
    "is itself a Markush structure",                    # 1 line: the template itself carries R groups
    "were drawn in FULL anyway",                        # abbreviate accepted but not honoured
    "were asked to float and could not",                # floating accepted but not honoured
    "matched no placeholder in the SMILES",             # rlabels key mistyped
    "had nowhere to land and were NOT written",         # AND/OR stereo marks not written to the file
    "did not load, so every molecule was computed",     # hand-drawn depiction library did not load
    "named an id another cell already held",            # duplicate id: that SMILES was not drawn
    "could not be built and are NOT in the file",       # whole entry missing
    "span(s) were not reached",                         # column span had no effect
    "cannot be shrunk with them",                       # reagent was shrunk, its decorations did not follow
)


def _drift(report: dict) -> list:
    """"You asked for it and it did not happen" family - what the engine said and
    _findings always threw away."""
    out = []
    for n in (report.get("notes") or []):
        for line in str(n).splitlines():
            line = line.strip()
            if any(ph in line for ph in _DRIFT_PHRASES):
                out.append(line[6:].strip() if line.startswith("NOTE:") else line)
    return out


# ---------------------------------------------------------------- parts

# ALIGN has two vocabularies, and callers should not have to remember which layer they
# are talking to. panels.py says start/center/end; here we speak per axis, left/right or
# top/bottom, because "left" is meaningful when stacking vertically and meaningless when
# laying out horizontally. Both are accepted: start/end are translated per axis.
# Existing spellings untouched.
#
# On a typo, list the words valid **on this axis**. Listing all four regardless of axis
# would throw away the information "you used the wrong axis".
_ALIGN_ALIAS = {"y": {"start": "left", "end": "right"},
                "x": {"start": "top", "end": "bottom"}}


def _pick(align, table, axis):
    key = _ALIGN_ALIAS[axis].get(align, align)
    if key not in table:
        raise ValueError(
            "align %r not recognised. Valid along this axis: %s (panels' start/center/end also accepted)"
            % (align, " / ".join(sorted(table))))
    return table[key]


_AT_DEG = {"right": 0.0, "below": 90.0, "left": 180.0, "above": 270.0}


def _one_ring(part, i):
    """Ring i of the part; if not found, say what was measured."""
    rr = part.rings()
    if len(rr) <= i:
        raise ValueError(
            "part %r has no ring %d: measured %d (radius %s). Only cycle parts have rings, "
            "and the arcs must still be there - if join / drop discarded the arcs they cannot be measured."
            % (part.name, i, len(rr),
               ", ".join("%.1f" % x["r"] for x in rr) or "none"))
    return rr[i]


class Part:
    def __init__(self, name, doc, report):
        self.name, self.doc, self.report = name, doc, report
        self.dx = self.dy = 0.0
        self.drop = set()          # indices of top-level objects to delete whole

    def kept(self):
        return [o for i, o in enumerate(self.doc.objs) if i not in self.drop]

    def ink_bbox(self):
        bs = [o.ink_bbox() for o in self.kept() if o.ink_bbox()]
        return [min(b[0] for b in bs) + self.dx, min(b[1] for b in bs) + self.dy,
                max(b[2] for b in bs) + self.dx, max(b[3] for b in bs) + self.dy]

    def rings(self):
        """The circle of each ring in this part: [{"c": (cx, cy), "r": R, "arcs": n}],
        with dx/dy applied.

        Taken from the **ink actually drawn**, not the coordinates the engine promised:
        ring arcs are ChemDraw <arrow>s, and the centre and radius are written in
        Center3D / MajorAxisEnd3D (the engine itself draws from these two numbers). The
        centre in the report is in the engine's layout space, off from the finished
        CDXML by finalize's translation, and the composer's placement needs the latter.

        The radius is not guessed: the report's rings[*].radius (ring.radius for a
        single-ring figure) already solved R; here only arc groups whose radius matches
        that number count. Feed arcs and pendant curved arrows have radii 11~26, rings
        start at 95, so they separate anyway; but "they separate" is a coincidence of
        this corpus, "equals the number the engine solved" is not.
        """
        want = []
        rep = self.report or {}
        rr = rep.get("rings")
        if isinstance(rr, dict):
            want = [float(v["radius"]) for v in rr.values()
                    if isinstance(v.get("radius"), (int, float))]
        elif isinstance((rep.get("ring") or {}).get("radius"), (int, float)):
            want = [float(rep["ring"]["radius"])]
        if not want:
            return []
        groups: dict = {}
        for m in re.finditer(r"<arrow\b[^>]*?>", self.body(), re.S):
            t = m.group(0)
            if "AngularSize" not in t or "MajorAxisEnd3D" not in t:
                continue
            a = dict(re.findall(r'(\w+)="([^"]*)"', t))
            cx, cy = [float(v) for v in a["Center3D"].split()[:2]]
            ax, ay = [float(v) for v in a["MajorAxisEnd3D"].split()[:2]]
            r = math.hypot(ax - cx, ay - cy)
            if not any(abs(r - w) < 0.5 for w in want):
                continue
            k = (round(cx, 1), round(cy, 1), round(r, 1))
            groups[k] = groups.get(k, 0) + 1
        # A ring has at least two arcs (one per edge). Two independently drawn arrows
        # will not have centres within 0.1 unit of each other, so this also blocks a
        # lone arc whose radius happens to equal R.
        out = [{"c": (cx + self.dx, cy + self.dy), "r": r, "arcs": n}
               for (cx, cy, r), n in groups.items() if n >= 2]
        out.sort(key=lambda d: (-d["r"], d["c"]))
        return out

    def body(self):
        spans = CD.top_objects(self.doc.body)
        return "\n".join(self.doc.body[s:e] for i, (_, s, e) in enumerate(spans)
                         if i not in self.drop)


# ---------------------------------------------------------------- composer

class Composer:
    def __init__(self, margin: float = 30.0, workdir: str | None = None):
        self.parts: list[Part] = []
        self.margin = margin
        self.workdir = workdir or os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "parts")
        self.notes: list[str] = []
        self.conflicts: list[dict] = []
        self.links: list[dict] = []    # cross-part arrows, see arrow()
        # Cross-part furniture: full-width title bars, frames around a group of parts,
        # section backgrounds, rules between sections. Their anchor is **the union of
        # several parts**, a number only this layer knows - on the engine side the
        # figure in `of: "figure"` is the ink of that one call, not the whole figure.
        self.furniture: list[dict] = []

    # -- run engines --------------------------------------------------------
    # The name is this layer's only address. Placement, `of`, `span`, join, arrow - each
    # accepts only a name, and lookup by name **finds only the first**. So two things
    # with one name is not "the later overrides the earlier", it is the later one **can
    # never be addressed again**: in #054 three bulbs and three donors were all called
    # b1/b2/b3, every placement landed on the bulbs, the three donors stayed at the
    # origin, and write() said nothing.
    #
    # Reject rather than note, because there is no correct fallback: picking the first
    # is arbitrary, picking the last just as arbitrary, and the one the caller wanted
    # may be neither. Renaming takes a second.
    #
    # Parts and furniture share this one namespace - `of=` accepts both, `_union` checks
    # furniture first then parts, so a piece of furniture named after a part **hides**
    # that part. One rule covers both.
    def _claim(self, name, what):
        """Claim a name; reject if already claimed. `what` is only used to make the
        error clear."""
        name = str(name)
        parts = [p.name for p in self.parts]
        furn = [f["name"] for f in self.furniture]
        held = "part" if name in parts else ("furniture" if name in furn else None)
        if held:
            raise ValueError(
                "name %r is already taken by a %s, cannot also give it to %s - this layer finds things by name and only finds"
                " the first; the later one could not be placed. Existing: parts %s; furniture %s"
                % (name, held, what, ", ".join(parts) or "none", ", ".join(furn) or "none"))
        return name

    def _run(self, tool, name, kw):
        import chemwand_mcp as M
        import json
        # Claim the name before calling the engine: with a duplicate name the drawing is
        # wasted, and the error should come before that layout run.
        name = self._claim(name, "part")
        fn = getattr(M, tool)
        with M.headless(self.workdir):
            rep = fn(filename="_part_" + name, **kw)
        rep = json.loads(rep) if isinstance(rep, str) else rep
        path = rep.get("path")
        if not path or not os.path.exists(path):
            raise RuntimeError("%s(%s) wrote no file: %s" % (tool, name, str(rep)[:400]))
        p = Part(name, CD.Doc(open(path).read(), name), rep)
        self.parts.append(p)
        return p

    def cycle(self, name, **kw):
        kw.setdefault("center", [])
        return self._run("draw_catalytic_cycle", name, kw)

    def linear(self, name, **kw):
        kw.setdefault("footnotes", [])
        kw.setdefault("page", "fit")
        return self._run("draw_linear_scheme", name, kw)

    def radial(self, name, **kw):
        return self._run("draw_radial_figure", name, kw)

    def table(self, name, **kw):
        return self._run("draw_scope_table", name, kw)

    def molecule(self, name, **kw):
        return self._run("draw_molecule", name, kw)

    def snaking(self, name, **kw):
        return self._run("draw_snaking_scheme", name, kw)

    def __getitem__(self, name) -> Part:
        for p in self.parts:
            if p.name == name:
                return p
        raise KeyError("no part named %r (have %s)"
                       % (name, ", ".join(p.name for p in self.parts)))

    # -- placement: by the actual extent of the ink -------------------------
    def _stack(self, who, ref, gap, axis, sign, align):
        """ref can be one name or a list of names - for a list, place against their
        **union**. A top bar has to sit above two rings; picking just one of them as
        reference skews it."""
        a = self[who]
        names = [ref] if isinstance(ref, str) else list(ref)
        boxes = [self[n].ink_bbox() for n in names]
        bb = [min(x[0] for x in boxes), min(x[1] for x in boxes),
              max(x[2] for x in boxes), max(x[3] for x in boxes)]
        ab = a.ink_bbox()
        if axis == "y":
            a.dy += (bb[1] - gap - ab[3]) if sign < 0 else (bb[3] + gap - ab[1])
            k = _pick(align, {"center": ((bb[0] + bb[2]) - (ab[0] + ab[2])) / 2,
                              "left": bb[0] - ab[0], "right": bb[2] - ab[2]}, "y")
            a.dx += k
            return a
        else:
            a.dx += (bb[0] - gap - ab[2]) if sign < 0 else (bb[2] + gap - ab[0])
            k = _pick(align, {"center": ((bb[1] + bb[3]) - (ab[1] + ab[3])) / 2,
                              "top": bb[1] - ab[1], "bottom": bb[3] - ab[3]}, "x")
            a.dy += k
        return a

    def above(self, who, ref, gap=26.0, align="center"):
        return self._stack(who, ref, gap, "y", -1, align)

    def below(self, who, ref, gap=26.0, align="center"):
        return self._stack(who, ref, gap, "y", +1, align)

    def left_of(self, who, ref, gap=26.0, align="center"):
        return self._stack(who, ref, gap, "x", -1, align)

    def right_of(self, who, ref, gap=26.0, align="center"):
        return self._stack(who, ref, gap, "x", +1, align)

    def nudge(self, who, dx=0.0, dy=0.0):
        p = self[who]
        p.dx += dx
        p.dy += dy
        return p

    # -- placement: by circle, not by bounding box --------------------------
    #
    # Why a separate path rather than a gap for below / right_of: those place by the
    # **ink bounding box**, and a ring part's bounding box is "ring + every molecule and
    # every label hung on it". When the boxes just touch, the two circles are still
    # apart by however far the pendants stick out. Measured the 17 multi-ring figures
    # already drawn: **not one is tangent**; the tightest, 050, is off by 0.3 bond
    # lengths, 016 by 2.2 (31.7 units). The engine's default placement has the same flaw
    # (layout_cycle's `_x += 2*(R+reach_max) + center_gap`), so even two rings with
    # single-letter nodes get pushed 37 units apart.
    #
    # And of the originals of these 17 figures in the corpus, 14 are **tangent**, with
    # an SET / PT / HAT badge sitting right on the tangent point - the tangent point is
    # where the electron transfer is drawn. So tangency is not a look preference, it is
    # the position of the chemical logic.
    #
    # Tangency has no threshold to tune: centre distance == R1 + R2, one equation. `gap`
    # is kept not for fine-tuning but so the two originals that already intersect (118 /
    # 120, a small Ir ring passing through a large Cu ring, with the shared electron
    # transfer step inside the lens) can be written as a negative number.
    def tangent(self, who, ref, at="right", gap=0.0, ring=0, ref_ring=0):
        """Move who so it is exactly tangent to ref. Returns (tangent point, measured
        gap).

        at  direction of who's centre from ref's centre. Give degrees (CDXML y points
            down, so 0=right, 90=below, 180=left, 270=above), or right / below / left /
            above.
        gap extra distance between the two circles. 0 = truly tangent; negative =
            deliberately intersecting.
        """
        a, b = self[who], self[ref]
        ca, cb = _one_ring(a, ring), _one_ring(b, ref_ring)
        if isinstance(at, str):
            if at not in _AT_DEG:
                raise ValueError("at %r not recognised; use %s, or give degrees directly"
                                 % (at, " / ".join(_AT_DEG)))
            deg = _AT_DEG[at]
        else:
            deg = float(at)
        th = math.radians(deg)
        d = ca["r"] + cb["r"] + gap
        a.dx += cb["c"][0] + d * math.cos(th) - ca["c"][0]
        a.dy += cb["c"][1] + d * math.sin(th) - ca["c"][1]
        # Tangent point: from ref's centre walk ref's own radius in that direction. With
        # gap=0 it lies on both circles, and that is where the badge belongs.
        pt = (cb["c"][0] + cb["r"] * math.cos(th),
              cb["c"][1] + cb["r"] * math.sin(th))
        got = self.ring_gap(who, ref, ring=ring, ref_ring=ref_ring)
        self.notes.append("%s tangent to %s: direction %.0f°, tangent point (%.1f, %.1f), "
                          "measured gap %+.2f" % (who, ref, deg, pt[0], pt[1], got))
        return pt, got

    def ring_gap(self, who, ref, ring=0, ref_ring=0):
        """Measured gap between the rings of two parts: centre distance - (R1+R2).

        0 = tangent, positive = apart, negative = intersecting. Use it to re-check after
        placement - tangency must be a check that runs, not "I tuned it".
        """
        ca, cb = _one_ring(self[who], ring), _one_ring(self[ref], ref_ring)
        d = math.hypot(ca["c"][0] - cb["c"][0], ca["c"][1] - cb["c"][1])
        return d - (ca["r"] + cb["r"])


    # -- placement: shared anchors -------------------------------------------
    def anchor_of(self, part, key, by):
        """(anchor coordinates, [objects to drop]) - by='species' uses the structure box centre, by='label' the label box centre."""
        p = self[part]
        if by == "label":
            lab = p.doc.label(key)
            return ((lab.bbox[0] + lab.bbox[2]) / 2, (lab.bbox[1] + lab.bbox[3]) / 2), [lab]
        f, lab = p.doc.species(key)
        if f is None:
            raise ValueError("%s: no structure found under label %r (for text-only species use by='label')"
                             % (part, key))
        return ((f.bbox[0] + f.bbox[2]) / 2, (f.bbox[1] + f.bbox[3]) / 2), [f, lab]

    def residual(self, who, to, at, my=None, by="species", by_to=None):
        """When both sides know this species, how far apart their positions still are.

        Each engine lays out on its own, so a species' relative position in the two
        outputs is set independently. Anchoring one fixes the difference at the other -
        nobody can cancel both at once; it is **the constraint conflict between the two
        engines itself**, not a bug. What the composer can do is measure it and report
        it.
        """
        a, b = self[who], self[to]
        ca, _ = self.anchor_of(who, my or at, by)
        cb, _ = self.anchor_of(to, at, by_to or by)
        return (round(cb[0] + b.dx - (ca[0] + a.dx), 1),
                round(cb[1] + b.dy - (ca[1] + a.dy), 1))

    def residual_shapes(self, who, to, drop=False):
        """After two parts are placed, list how far apart every shared structure still
        is (moves nothing).

        The `#098` kind: two rings share one **node** (text) and also one **species**
        (structure). Only one can be the anchor; the other's difference has to be
        reported. With drop=True also delete the duplicate structures in who - matched
        or not, drawing one molecule twice is always wrong.
        """
        a, b = self[who], self[to]
        sa, sb = a.doc.shapes(), b.doc.shapes()
        out = []
        for k in sa:
            if k not in sb or len(sa[k]) != 1 or len(sb[k]) != 1:
                continue
            r = (round(sb[k][0].bbox[0] + b.dx - sa[k][0].bbox[0] - a.dx, 1),
                 round(sb[k][0].bbox[1] + b.dy - sa[k][0].bbox[1] - a.dy, 1))
            out.append(r)
            self.conflicts.append({"part": who, "vs": to, "key": "shape", "gap": list(r)})
            self.notes.append("  ⚠ %s/%s also share a structure, off by (%.1f, %.1f)" % (who, to, r[0], r[1]))
            if drop:
                f = sa[k][0]
                lab = a.doc.label_of(f)
                for i, o in enumerate(a.doc.objs):
                    if o is f or (lab is not None and o is lab):
                        a.drop.add(i)
        return out

    def join_shape(self, who, to, on=0, drop=True):
        """Align the structures in who and to that **draw the same molecule**, without
        relying on labels.

        Find all species the two share, align on the `on`-th one, and measure and report
        the residuals of the rest - a chain attached to rings at both ends (#026's
        imidoyl->...->adduct) necessarily has one end that does not match; that is the
        constraint conflict after each engine laid out on its own, not something with a
        "right" answer to pick. With drop=True the duplicate structures in who are
        deleted whole, so shared species are drawn only once.
        """
        a, b = self[who], self[to]
        sa, sb = a.doc.shapes(), b.doc.shapes()
        common = [k for k in sa if k in sb and len(sa[k]) == 1 and len(sb[k]) == 1]
        if not common:
            raise ValueError("%s and %s do not draw the same molecule" % (who, to))
        common.sort(key=lambda k: sa[k][0].bbox[1])       # top to bottom, easy to point at
        pick = common[on]
        fa, fb = sa[pick][0], sb[pick][0]
        a.dx += (fb.bbox[0] + fb.bbox[2]) / 2 + b.dx - (fa.bbox[0] + fa.bbox[2]) / 2
        a.dy += (fb.bbox[1] + fb.bbox[3]) / 2 + b.dy - (fa.bbox[1] + fa.bbox[3]) / 2
        self.notes.append("%s aligned to %s by shared structure (%d shared, used #%d)"
                          % (who, to, len(common), on))
        for i, k in enumerate(common):
            if i != on:
                r = (round(sb[k][0].bbox[0] + b.dx - sa[k][0].bbox[0] - a.dx, 1),
                     round(sb[k][0].bbox[1] + b.dy - sa[k][0].bbox[1] - a.dy, 1))
                self.notes.append("  ⚠ shared structure #%d does not match: off by (%.1f, %.1f)" % (i, r[0], r[1]))
                self.conflicts.append({"part": who, "vs": to, "key": "shape#%d" % i,
                                       "gap": list(r)})
        if drop:
            for k in common:
                f = sa[k][0]
                lab = a.doc.label_of(f)        # the label must go too, or the number stays behind
                for i, o in enumerate(a.doc.objs):
                    if o is f or (lab is not None and o is lab):
                        a.drop.add(i)
        return a

    def join(self, who, to, at, my=None, drop=True, by="species",
             by_to=None, also=None):
        """Align the species called `my` (default: same name as `at`) in who onto the
        species called `at` in to.

        What is aligned is the **structure's bbox centre**, not the label position - the
        two engines do not necessarily give labels the same line height and
        justification; the structure itself is what both sides recognise as the same
        thing. With drop=True the duplicate structure and label in who are deleted
        whole, so this species is drawn only once in the finished figure, and no
        cross-engine arrow is needed.
        """
        a, b = self[who], self[to]
        my = my or at
        ca, mine = self.anchor_of(who, my, by)
        cb, theirs = self.anchor_of(to, at, by_to or by)
        a.dx += cb[0] - ca[0] + b.dx
        a.dy += cb[1] - ca[1] + b.dy
        self.notes.append("anchored %s.%s → %s.%s, offset (%.1f, %.1f)"
                          % (who, my, to, at, a.dx, a.dy))
        if drop:
            for i, o in enumerate(a.doc.objs):
                if o in mine:
                    a.drop.add(i)
            self.notes.append("  in %s, dropped the duplicate %r, keeping only the one in %s" % (who, my, to))
        if also:
            k, kb = (also, also) if isinstance(also, str) else also
            r = self.residual(who, to, kb, my=k, by=by, by_to=by_to)
            self.notes.append("  ⚠ second anchor %r does not match: off by (%.1f, %.1f). "
                              "Each engine laid out on its own; nobody can cancel this difference at the same time" % (k, r[0], r[1]))
            self.conflicts.append({"part": who, "vs": to, "key": k, "gap": list(r)})
        return a

    def arrow_between(self, part, a, b):
        """The arrow in part connecting species a and b.

        Found geometrically: the one spanning between the two structures and overlapping
        them vertically. Does not use the engine's `placement` output, because linear's
        placement reports only the arrow tail/head, **no geometry for any species**
        (layout_linear.py:1000-1030), so names cannot be matched.
        """
        p = self[part]
        fa, _ = p.doc.species(a)
        fb, _ = p.doc.species(b)
        if fa is None or fb is None:
            raise ValueError("%s: no structure found for %r or %r" % (part, a, b))
        lo, hi = sorted((fa.bbox[2], fb.bbox[0])) if fa.bbox[2] <= fb.bbox[0] \
            else sorted((fb.bbox[2], fa.bbox[0]))
        ymid = (fa.bbox[1] + fa.bbox[3] + fb.bbox[1] + fb.bbox[3]) / 4
        hits = [o for o in p.kept() if o.tag == "arrow" and o.ink_bbox()
                and lo - 4 <= (o.ink_bbox()[0] + o.ink_bbox()[2]) / 2 <= hi + 4
                and o.ink_bbox()[1] - 30 <= ymid <= o.ink_bbox()[3] + 30]
        if len(hits) != 1:
            raise ValueError("%s: %r→%r: found %d arrows between them" % (part, a, b, len(hits)))
        return hits[0]

    def under_arrow(self, who, ref, a, b, gap=0.0, align="center"):
        """Place the whole of who under the a→b arrow in ref, centred horizontally.

        This is the placement when there is **no shared species**: a catalytic cycle and
        the main chain share not one molecule; the only relation between them is "this
        cycle explains this step", and that is a geometric relation. gap can be negative
        - in the original the ring's top arc touches the arrow.
        """
        w, r = self[who], self[ref]
        ar = self.arrow_between(ref, a, b)
        ab = ar.ink_bbox()
        wb = w.ink_bbox()
        cx = (ab[0] + ab[2]) / 2 + r.dx
        w.dx += cx - (wb[0] + wb[2]) / 2 if align == "center" else 0.0
        w.dy += (ab[3] + r.dy + gap) - wb[1]
        self.notes.append("%s placed below %s's %s→%s arrow (gap %.0f)" % (who, ref, a, b, gap))
        return w

    # -- cross-part arrows: the one arrow only the parent layer can draw ------

    def arrow(self, frm, to, *, bulge=0.0, gap=10.0, head="Full", tail=None,
              color=None, line=None):
        """An arrow with one end in this block and the other in that one.

        Engines cannot draw it, and this is no omission: each engine knows only what
        **it** laid out; for an arrow whose ends belong to two calls, only the parent
        layer, after placement, knows both endpoints. #173's original has exactly this
        shape - `•CCl₃` comes out of the cycle and goes straight to the 8-3→8-4 step of
        the main chain; today it can only be drawn once in the ring and once in the
        chain, and the reader has to recognise it is the same species.

        Endpoints are written as (part name, designator):
          "8-3"            the structure that label designates
          ("8-2", "8-3")   the arrow between these two structures (a reagent landing
                           on a step is this)

        Resolution is deferred to build(): placement may still change (under_arrow,
        nudge both come later), so endpoints must be computed from the **final**
        positions, not those at declaration time.
        """
        self.links.append(dict(frm=frm, to=to, bulge=float(bulge), gap=float(gap),
                               head=head, tail=tail, color=color, line=line))
        return self

    def _link_box(self, ref):
        """(part name, designator) -> this block's ink box in composed coordinates."""
        name, what = ref
        p = self[name]
        if isinstance(what, (tuple, list)):
            o = self.arrow_between(name, *what)
        else:
            # species() returns (structure, label). **For a species written as text the
            # structure is None** - the ring's `•CCl₃ + Br⁻` is like that; it has only a
            # <t>. That is a legal endpoint too: an arrow pointing at "the species
            # written out" is no different from pointing at a drawn one.
            frag, lab = p.doc.species(str(what))
            o = frag or lab
        if o is None:
            raise ValueError("%s: %r not found" % (name, what))
        b = o.ink_bbox() or o.bbox
        return [b[0] + p.dx, b[1] + p.dy, b[2] + p.dx, b[3] + p.dy]

    @staticmethod
    def _edge(box, toward, gap):
        """Walk from box's centre toward toward, land on the box edge, then back off by
        gap.

        Endpoints take the edge, not the centre: an arrow growing out of a structure's
        centre would sit on its own bonds.
        """
        cx, cy = (box[0] + box[2]) / 2, (box[1] + box[3]) / 2
        dx, dy = toward[0] - cx, toward[1] - cy
        n = math.hypot(dx, dy) or 1.0
        hw, hh = (box[2] - box[0]) / 2, (box[3] - box[1]) / 2
        # scale to the rectangle edge: take whichever direction hits first
        t = min(hw / abs(dx) if dx else 1e18, hh / abs(dy) if dy else 1e18)
        return (cx + dx * t + dx / n * gap, cy + dy * t + dy / n * gap)

    def _build_links(self, aid0=4_900_000):
        """Resolve the declared cross-part arrows into (xml, polyline). Called in
        build().

        Id range 4.9M: above any id a single part uses (the largest is _add_graphics's
        970000+si*20000), and below part 1's offset of 5M.
        """
        out = []
        for i, lk in enumerate(self.links):
            a, b = self._link_box(lk["frm"]), self._link_box(lk["to"])
            ca = ((a[0] + a[2]) / 2, (a[1] + a[3]) / 2)
            cb = ((b[0] + b[2]) / 2, (b[1] + b[3]) / 2)
            p0 = self._edge(a, cb, lk["gap"])
            p1 = self._edge(b, ca, lk["gap"])
            if lk["bulge"]:
                ch = (p1[0] - p0[0], p1[1] - p0[1])
                L = math.hypot(*ch) or 1.0
                nrm = (-ch[1] / L, ch[0] / L)
                m = ((p0[0] + p1[0]) / 2 + nrm[0] * lk["bulge"] * L,
                     (p0[1] + p1[1]) / 2 + nrm[1] * lk["bulge"] * L)
                xml = CORE.bezier_arrow(aid0 + i, p0, (m[0] - p0[0], m[1] - p0[1]),
                                        L * 0.35, p1, (p1[0] - m[0], p1[1] - m[1]),
                                        L * 0.35, head=lk["head"], tail=lk["tail"],
                                        color=lk["color"], line=lk["line"])
                pts = CORE.bezier_points_from(
                    [(p0[0] - (m[0] - p0[0]), p0[1] - (m[1] - p0[1])), p0, m,
                     (p1[0] - (p1[0] - m[0]), p1[1] - (p1[1] - m[1])), p1,
                     (p1[0] + (p1[0] - m[0]), p1[1] + (p1[1] - m[1]))])
            else:
                xml = CORE.arrow_xy(aid0 + i, p0[0], p0[1], p1[0], p1[1])
                pts = [p0, p1]
            out.append((xml, pts))
        return out

    # -- cross-part arrows: the one arrow only the parent layer can draw ------

    def arrow(self, frm, to, *, bulge=0.0, gap=10.0, head="Full", tail=None,
              color=None, line=None):
        """An arrow with one end in this block and the other in that one.

        Engines cannot draw it, and this is no omission: each engine knows only what
        **it** laid out; for an arrow whose ends belong to two calls, only the parent
        layer, after placement, knows both endpoints. #173's original has exactly this
        shape - `•CCl₃` comes out of the cycle and goes straight to the 8-3→8-4 step of
        the main chain; today it can only be drawn once in the ring and once in the
        chain, and the reader has to recognise it is the same species.

        Endpoints are written as (part name, designator):
          "8-3"            the structure that label designates
          ("8-2", "8-3")   the arrow between these two structures (a reagent landing
                           on a step is this)

        Resolution is deferred to build(): placement may still change (under_arrow,
        nudge both come later), so endpoints must be computed from the **final**
        positions, not those at declaration time.
        """
        self.links.append(dict(frm=frm, to=to, bulge=float(bulge), gap=float(gap),
                               head=head, tail=tail, color=color, line=line))
        return self

    def _link_box(self, ref):
        """(part name, designator) -> this block's ink box in composed coordinates."""
        name, what = ref
        p = self[name]
        if isinstance(what, (tuple, list)):
            o = self.arrow_between(name, *what)
        else:
            # species() returns (structure, label). **For a species written as text the
            # structure is None** - the ring's `•CCl₃ + Br⁻` is like that; it has only a
            # <t>. That is a legal endpoint too: an arrow pointing at "the species
            # written out" is no different from pointing at a drawn one.
            frag, lab = p.doc.species(str(what))
            o = frag or lab
        if o is None:
            raise ValueError("%s: %r not found" % (name, what))
        b = o.ink_bbox() or o.bbox
        return [b[0] + p.dx, b[1] + p.dy, b[2] + p.dx, b[3] + p.dy]

    @staticmethod
    def _edge(box, toward, gap):
        """Walk from box's centre toward toward, land on the box edge, then back off by
        gap.

        Endpoints take the edge, not the centre: an arrow growing out of a structure's
        centre would sit on its own bonds.
        """
        cx, cy = (box[0] + box[2]) / 2, (box[1] + box[3]) / 2
        dx, dy = toward[0] - cx, toward[1] - cy
        n = math.hypot(dx, dy) or 1.0
        hw, hh = (box[2] - box[0]) / 2, (box[3] - box[1]) / 2
        # scale to the rectangle edge: take whichever direction hits first
        t = min(hw / abs(dx) if dx else 1e18, hh / abs(dy) if dy else 1e18)
        return (cx + dx * t + dx / n * gap, cy + dy * t + dy / n * gap)

    def _build_links(self, aid0=4_900_000):
        """Resolve the declared cross-part arrows into (xml, polyline). Called in
        build().

        Id range 4.9M: above any id a single part uses (the largest is _add_graphics's
        970000+si*20000), and below part 1's offset of 5M.
        """
        out = []
        for i, lk in enumerate(self.links):
            a, b = self._link_box(lk["frm"]), self._link_box(lk["to"])
            ca = ((a[0] + a[2]) / 2, (a[1] + a[3]) / 2)
            cb = ((b[0] + b[2]) / 2, (b[1] + b[3]) / 2)
            p0 = self._edge(a, cb, lk["gap"])
            p1 = self._edge(b, ca, lk["gap"])
            if lk["bulge"]:
                ch = (p1[0] - p0[0], p1[1] - p0[1])
                L = math.hypot(*ch) or 1.0
                nrm = (-ch[1] / L, ch[0] / L)
                m = ((p0[0] + p1[0]) / 2 + nrm[0] * lk["bulge"] * L,
                     (p0[1] + p1[1]) / 2 + nrm[1] * lk["bulge"] * L)
                xml = CORE.bezier_arrow(aid0 + i, p0, (m[0] - p0[0], m[1] - p0[1]),
                                        L * 0.35, p1, (p1[0] - m[0], p1[1] - m[1]),
                                        L * 0.35, head=lk["head"], tail=lk["tail"],
                                        color=lk["color"], line=lk["line"])
                pts = CORE.bezier_points_from(
                    [(p0[0] - (m[0] - p0[0]), p0[1] - (m[1] - p0[1])), p0, m,
                     (p1[0] - (p1[0] - m[0]), p1[1] - (p1[1] - m[1])), p1,
                     (p1[0] + (p1[0] - m[0]), p1[1] + (p1[1] - m[1]))])
            else:
                xml = CORE.arrow_xy(aid0 + i, p0[0], p0[1], p1[0], p1[1])
                pts = [p0, p1]
            out.append((xml, pts))
        return out

    # -- collisions: the only place that knows where two engines both took space --
    def collisions(self, tol: float = 1.0, pad: float = 0.0):
        """With pad>0 each object is inflated by pad before testing, to catch "not
        overlapping but too close".

        Overlap and whitespace are two problems - the allocator splits them the same way
        (overlap uses depth, whitespace uses `probe(pad=...)`). In #190 the ring's 5CzBN
        node label was only 12 units from the chain's arrow, no overlap anywhere, yet it
        read as part of the chain."""
        objs = []
        # The parent layer's own arrows must go through the ruler too. They are the only
        # ink in this figure that belongs to no engine; skipping them means putting an
        # unchecked thing on a figure that "exists specifically to check overlaps".
        class _Link:                      # only to supply the tag/text fields
            tag = "arrow"
            def __init__(self, i): self.i = i
            def text(self): return "parent arrow %d" % self.i
        for i, (_xml, pts) in enumerate(self._build_links()):
            xs = [q[0] for q in pts]; ys = [q[1] for q in pts]
            objs.append(("(parent)", _Link(i),
                         [min(xs) - pad, min(ys) - pad, max(xs) + pad, max(ys) + pad],
                         [] if pad else [(pts[k][0], pts[k][1], pts[k + 1][0], pts[k + 1][1])
                                         for k in range(len(pts) - 1)]))
        for p in self.parts:
            for o in p.kept():
                bb = o.ink_bbox()
                if not bb:
                    continue
                # A <t> that is only whitespace draws nothing. They come from the
                # label:" " trick - the engine requires every species to have a label,
                # and to print none you give a space, leaving in the file an object with
                # a box but no ink. Testing it for collisions is a pure false positive.
                if o.tag == "t" and not o.text().strip():
                    continue
                objs.append((p.name, o,
                             [bb[0] + p.dx - pad, bb[1] + p.dy - pad,
                              bb[2] + p.dx + pad, bb[3] + p.dy + pad],
                             [(s[0] + p.dx, s[1] + p.dy, s[2] + p.dx, s[3] + p.dy)
                              for s in (o.ink or [])] if not pad else []))
        out = []
        for i in range(len(objs)):
            for j in range(i + 1, len(objs)):
                pa, oa, ba, sa = objs[i]
                pb, ob, bb, sb = objs[j]
                if pa == pb:
                    continue                     # overlaps inside one engine are that engine's own business
                dx, dy = _ovl(ba, bb)
                if dx <= 0 or dy <= 0:
                    continue
                if sa and sb:
                    hit = any(_seg_seg(s, t) for s in sa for t in sb)
                    kind, sev = "line×line", 0.0
                elif sa or sb:
                    segs, box = (sa, bb) if sa else (sb, ba)
                    hit = any(_seg_box(s, box) for s in segs)
                    kind, sev = "line×block", 0.0
                else:
                    hit = min(dx, dy) > tol
                    kind, sev = "block×block", min(dx, dy)
                if hit:
                    out.append({"a": "%s/%s %s" % (pa, oa.tag, oa.text()[:18]),
                                "b": "%s/%s %s" % (pb, ob.tag, ob.text()[:18]),
                                "kind": kind, "depth": round(sev, 1),
                                "where": [round(max(ba[0], bb[0])), round(max(ba[1], bb[1]))]})
        return sorted(out, key=lambda r: -r["depth"])

    # -- who is on top: draw order + opacity ---------------------------------
    #
    # Measured 11 RectangleTypes (scratchpad/opaque_probe, text emitted first, rectangle
    # after): those with Filled or Shaded in the name swallow the text beneath entirely;
    # Shadow / Shadow Bold / RoundEdge Shadow have only a border, the inside is
    # transparent and the text still reads. Discs (Orbital's OvalType) follow the same
    # rule; polygons look at bit 128 of CurveType (1=Closed, 128=Filled).
    _OPAQUE = re.compile(r"Filled|Shaded")

    @staticmethod
    def _opaque(o) -> bool:
        """Whether this object completely covers what was drawn before it."""
        if Composer._OPAQUE.search(o.a.get("RectangleType")
                                   or o.a.get("OvalType") or ""):
            return True
        if o.tag == "curve":
            try:
                return bool(int(o.a.get("CurveType", "0")) & 128)
            except ValueError:
                return False
        return False

    # An opaque object hides **its own shape**, not its bounding box. Of 26 first-pass
    # hits in the corpus, 20 were blobs (polygons); by bounding box, text outside a
    # blob's corners was also judged "covered" - #090's "a" sits in the top-left corner,
    # still some way from the grey blob.
    _COVER_GRID = 9          # take 9x9=81 points on the victim box, resolution about 1.2%

    @staticmethod
    def _region(o):
        """(bounding box, "is the point inside?" test) of an opaque object. None means the shape cannot be read.

        The bounding box does not read BoundingBox: ChemDraw's Orbital / Symbol write it
        "right side first, zero height" (the comment on hand-drawn files in cdxml_core
        records this); used as a rectangle it would give an inverted box, so that object
        could never be judged overlapping - a silent miss. Center3D / MajorAxisEnd3D /
        MinorAxisEnd3D are geometric quantities, and both rectangles and discs carry
        them.
        """
        a = o.a
        if o.tag == "curve" and o.ink:
            pts = [(s[0], s[1]) for s in o.ink] + [(o.ink[-1][2], o.ink[-1][3])]
            xs = [q[0] for q in pts]; ys = [q[1] for q in pts]

            def inside(x, y, _p=pts):
                hit = False
                for k in range(len(_p) - 1):
                    x0, y0 = _p[k]; x1, y1 = _p[k + 1]
                    if (y0 > y) != (y1 > y) and \
                            x < x0 + (y - y0) * (x1 - x0) / (y1 - y0):
                        hit = not hit
                return hit
            return [min(xs), min(ys), max(xs), max(ys)], inside
        if "Center3D" in a and "MajorAxisEnd3D" in a and "MinorAxisEnd3D" in a:
            cx, cy = CD._nums(a["Center3D"])[:2]
            rx = abs(CD._nums(a["MajorAxisEnd3D"])[0] - cx)
            ry = abs(CD._nums(a["MinorAxisEnd3D"])[1] - cy)
            box = [cx - rx, cy - ry, cx + rx, cy + ry]
            if a.get("GraphicType") == "Orbital":
                def inside(x, y):
                    return ((x - cx) / max(rx, 1e-9)) ** 2 + \
                           ((y - cy) / max(ry, 1e-9)) ** 2 <= 1.0
            else:
                def inside(x, y):
                    return box[0] <= x <= box[2] and box[1] <= y <= box[3]
            return box, inside
        return None, None

    @classmethod
    def _hidden_fraction(cls, victim, box, inside) -> float:
        """How much of the victim box falls inside that shape.

        Sample rather than intersect bounding boxes: a structure's box is mostly white,
        and blobs are not square, so both ends would skew. Lay a grid on the victim box
        and ask point by point "is it inside": for text (the box is solid ink) this is
        accurate; for structures it still runs high - a molecule's box has empty space
        in it anyway. Better high: missing one invisible label costs more than one extra
        report.
        """
        n = cls._COVER_GRID
        w, h = victim[2] - victim[0], victim[3] - victim[1]
        if w <= 0 or h <= 0:
            return 1.0 if inside((victim[0] + victim[2]) / 2.0,
                                 (victim[1] + victim[3]) / 2.0) else 0.0
        hits = 0
        for i in range(n):
            x = victim[0] + (i + 0.5) * w / n
            for j in range(n):
                if inside(x, victim[1] + (j + 0.5) * h / n):
                    hits += 1
        return hits / float(n * n)

    def _furniture_objs(self):
        """Each top-level object of the furniture + which piece of furniture it belongs
        to, in the parts' shared coordinates.

        _build_furniture returns only one union box per piece of furniture, while a text
        furniture with fill emits two objects, "solid rectangle + text" - by the union,
        that rectangle would judge its own text covered. So the xml is split again here,
        with the same top_objects / Obj used for reading parts.

        The colour table is passed empty: this pass measures geometry only; whether
        colour indices start at 2 or 9 changes no coordinate.
        """
        if not self.furniture:
            return []
        xmls, _boxes = self._build_furniture([])
        out = []
        for f, x in zip(self.furniture, xmls):
            for t, s, e in CD.top_objects(x):
                out.append((f["name"], CD.Obj(t, x[s:e])))
        return out

    # Default 0.10, picked by Ge from a ladder run over 90 figure scripts:
    #     >= 1%  38 hits/16 figures    >= 10%  23 hits/11 figures    >= 50%  6 hits/5 figures
    # Before choosing the number, "compared with what" was settled - with the victim's
    # own area (see the hid section). The value sits at 10% because **the lowest
    # verified rung is here**: #090 panels a/b/d titles report 22% / 11% / 10%, all
    # three on the real figure are titles bitten through the middle by a grey blob,
    # obviously wrong at a glance. The rung below (1%~7%, 14 hits, mostly decorative
    # spheres grazing a label edge) was not checked one by one on real figures, so it is
    # not in the default. To see it, pass covered(hid=0.01) explicitly.
    COVERED_HID = 0.10

    def covered(self, hid: float | None = None):
        """Who is covered by a solid block drawn **later**, and what percentage of
        itself is covered.

        A different question from collisions(). collisions() asks "did two things
        collide"; this asks
        "of the two that collided, which is on top". A whole class of overlaps in figures is **deliberate** - protein blobs under
        molecules, background bars behind text, badges on frames - count them as
        collisions and one figure can report dozens, drowning the real accidents. What
        separates them is no marker but **draw order**: the block drawn first is
        background; the block drawn later is the accident.

        Draw order is fixed in build(), so "who is on top" needs no guessing: furniture
        (in order added) -> parts
        (in declaration order, each one's Z stacked on top in turn) -> parent arrows. A figure without furniture does not rewrite Z; then
        parts keep their own Z as before, and the order is still declaration order.

        `hid` is compared with **the victim's own area**, not with the covering block,
        and not in absolute pt: a 400x60 background bar covering a 30x10 label has an
        overlap depth of only 10 pt, which sounds small, yet that label is 100% gone.
        collisions() reports depth because it asks "how deep did it go in".

        Only two kinds of victim count: text (<t>) and structures (<fragment>). A buried
        line or arrow is not reported here - the engine's own _reconcile already reports
        it inside a part, and cross-part lines are rare anyway.
        """
        hid = self.COVERED_HID if hid is None else hid
        objs = []                       # strictly in draw order: (source, object, box, dx, dy)
        for src, o in self._furniture_objs():
            bb = o.ink_bbox()
            if bb:
                objs.append((src, o, bb, 0.0, 0.0))
        for p in self.parts:
            for o in p.kept():
                bb = o.ink_bbox()
                if bb:
                    objs.append((p.name, o, [bb[0] + p.dx, bb[1] + p.dy,
                                             bb[2] + p.dx, bb[3] + p.dy],
                                 p.dx, p.dy))
        out = []
        for j, (sb, ob, bb, ox, oy) in enumerate(objs):
            if not self._opaque(ob):
                continue
            box, inside = self._region(ob)
            if inside is None:
                continue                         # if the shape cannot be read, do not guess
            if ox or oy:                         # an object inside a part: its attributes are still in part coordinates
                box = [box[0] + ox, box[1] + oy, box[2] + ox, box[3] + oy]
                inside = (lambda f, a, b: lambda x, y: f(x - a, y - b))(inside, ox, oy)
            for sa, oa, ba, _ax, _ay in objs[:j]:      # only look at what was drawn **before** it
                if sa == sb:
                    continue                     # layer order within one source is that engine's own business
                if oa.tag not in ("t", "fragment"):
                    continue
                if oa.tag == "t" and not oa.text().strip():
                    continue                     # the label:" " placeholder trick, no ink
                dx, dy = _ovl(ba, box)
                if dx <= 0 or dy <= 0:
                    continue                     # coarse bounding-box screen first, saves most of the sampling
                frac = self._hidden_fraction(ba, box, inside)
                if frac < hid:
                    continue
                out.append({"a": "%s/%s %s" % (sa, oa.tag, oa.text()[:18]),
                            "by": "%s/%s %s" % (
                                sb, ob.tag, ob.a.get("RectangleType")
                                or ob.a.get("OvalType") or "poly"),
                            "hid": round(frac, 2),
                            "where": [round(max(ba[0], box[0])),
                                      round(max(ba[1], box[1]))]})
        return sorted(out, key=lambda r: -r["hid"])

    # -- output -------------------------------------------------------------
    def build(self):
        if not self.parts:
            raise ValueError("no parts at all")
        head = self.parts[0].doc.header
        fonts = {re.search(r'id="(\d+)"', f).group(1): f
                 for f in _table(head, "fonttable")}
        colors = _table(head, "colortable")

        bodies, zbase, idbase = [], 0, 0
        for k, p in enumerate(self.parts):
            body = p.body()
            if k:
                body = _renumber_extra(CORE.renumber(body, idbase), idbase)
            # The next part's id range starts above the largest id this part actually
            # uses. A fixed STRIDE × k assumes each part fits in one range; an
            # already-merged file takes several ranges, and merging it again lands the
            # next part in its second range, colliding ids (3 merged figures merged with
            # 1 more: 25 duplicate ids). Ordinary parts fit in one range, so the result
            # is still STRIDE × k and the output is unchanged.
            top = max((int(i) for i in re.findall(r'\bid="(\d+)"', body)), default=0)
            idbase = max(idbase, (top // ID_STRIDE + 1) * ID_STRIDE)
            if k:
                # font table: union by id; same id with a different font must be
                # reported
                for f in _table(p.doc.header, "fonttable"):
                    fid = re.search(r'id="(\d+)"', f).group(1)
                    if fid in fonts and fonts[fid] != f:
                        self.notes.append("font id=%s is defined differently in two parts; the first one was kept" % fid)
                    fonts.setdefault(fid, f)
                # colour table: append entries beyond the common prefix, and rewrite
                # this part's indices to match
                mine = _table(p.doc.header, "colortable")
                remap = {}
                for idx, c in enumerate(mine):
                    if idx < len(colors) and colors[idx] == c:
                        continue
                    if c in colors:
                        remap[idx + 2] = colors.index(c) + 2
                    else:
                        colors.append(c)
                        remap[idx + 2] = len(colors) + 1
                if remap:
                    for attr in _COLOR_ATTRS:
                        body = re.sub(
                            r'\b' + attr + r'="(\d+)"',
                            lambda m: '%s="%d"' % (attr, remap.get(int(m.group(1)),
                                                                   int(m.group(1)))),
                            body)
                    self.notes.append("%s colour indices remapped %s" % (p.name, remap))
            # Z: every engine renumbers from 1, so plain concatenation scrambles layer
            # order. Stack in declaration order.
            if p.doc.max_z():
                body = _shift_z(body, zbase)
                zbase += p.doc.max_z()
            bodies.append(CORE._translate_body(body, p.dx, p.dy))

        # Furniture is built **after the tables are merged and before the header is
        # rewritten**: its [r,g,b] must be appended to colors in place, and colors is
        # written back to the header a few lines below. Emitted before all parts,
        # together with the Z restamp later, it really ends up behind the ink.
        furn_xml, furn_boxes = self._build_furniture(colors)
        bodies = furn_xml + bodies

        head = re.sub(r'<fonttable>.*?</fonttable>',
                      "<fonttable>\n%s\n</fonttable>" % "\n".join(
                          fonts[k] for k in sorted(fonts, key=int)), head, flags=re.S)
        head = re.sub(r'<colortable>.*?</colortable>',
                      "<colortable>\n%s\n</colortable>" % "\n".join(colors),
                      head, flags=re.S)

        # The parent layer's own arrows. Emitted after all parts, so they sit on top -
        # they connect two parts, and covered by either one they could not be read.
        links = self._build_links()
        bodies += [xml for xml, _pts in links]

        bs = [p.ink_bbox() for p in self.parts]
        cb = [min(b[0] for b in bs), min(b[1] for b in bs),
              max(b[2] for b in bs), max(b[3] for b in bs)]
        for _xml, _pts in links:      # the page box must hold them, or ChemDraw clips them
            cb[0] = min(cb[0], min(q[0] for q in _pts))
            cb[1] = min(cb[1], min(q[1] for q in _pts))
            cb[2] = max(cb[2], max(q[0] for q in _pts))
            cb[3] = max(cb[3], max(q[1] for q in _pts))
        for _b in furn_boxes:         # likewise: full-width bars are wider than the ink, frames bigger than it
            if _b:
                cb[0], cb[1] = min(cb[0], _b[0]), min(cb[1], _b[1])
                cb[2], cb[3] = max(cb[2], _b[2]), max(cb[3], _b[3])
        # One positive-quadrant translation + one page box. Each engine already
        # finalized once, and after stitching there are negative coordinates again, and
        # when ChemDraw opens the file it pins the page origin at (0,0) and clips the
        # negative part.
        ox, oy = self.margin - cb[0], self.margin - cb[1]
        body = CORE._translate_body("\n".join(bodies), ox, oy)
        if furn_xml:
            # DOCUMENT ORDER IS NOT DRAWING ORDER - cdxml_core.stamp_z proved this with
            # three discs covered by a halo. Parts come in each carrying its own Z, and
            # _shift_z stacks them; furniture has to sit behind all ink, and "emitted
            # early" alone cannot do that. So strip once over the whole body and restamp
            # once in document order, the same thing _add_graphics does inside a single
            # engine, and likewise **only when there is furniture** - a figure without
            # furniture is not touched by a single byte, otherwise every existing
            # figure's golden would have to be re-recorded.
            body = re.sub(r'\s+Z="\d+"', "", body)
            body, _znext = CORE.stamp_z(body)
        pw, ph = cb[2] - cb[0] + 2 * self.margin, cb[3] - cb[1] + 2 * self.margin
        return CORE.page(head, body, [0, 0, round(pw), round(ph)]), (ox, oy), (pw, ph)

    # -- cross-part furniture ------------------------------------------------
    #
    # Why this layer is indispensable: the engine-side `of: "figure"` means **the ink
    # drawn by that one call**. The dark-grey title bar across the top of 126 has to
    # span the whole figure; it belongs to no single engine call, and it is not an arrow
    # either - apart from `_build_links`, the composer had nothing else it drew itself.
    #
    # Same vocabulary as the engine-side `of`/`side`
    # (above/below/left/right/center/around), because to the reader they are the same
    # thing. **The only difference is written here**: the composer layer has no ledger,
    # so it does not search or avoid; it lands wherever the anchor computation puts it.
    # Real collisions are reported by collisions() - that is
    # "the only place that knows where two engines both took space" anyway.

    _FURN_SIDES = ("above", "below", "left", "right", "center", "around")

    def graphic(self, kind, of="all", side="above", *, align="center",
                gap=8.0, pad=6.0, box=None, span=None, name=None, **spec):
        """One piece of cross-part furniture. `of` / `span` accept a part name, a
        furniture name, a list of names, or "all".

        `of` positions it, `span` sets the size on the other axis (defaults to `of`).
        They are separate because 126's full-width bar has to "sit above top" but be "as
        wide as the whole figure" - one anchor cannot answer both.

        `box` applies only to placed rect / line; None on an axis = take `span`'s size
        on that axis.
        """
        if side not in self._FURN_SIDES:
            raise ValueError("side must be one of %s, got %r"
                             % (" / ".join(self._FURN_SIDES), side))
        # The automatic name furn%d goes through this check too: after a caller
        # hand-writes "furn3", the 4th furniture's automatic name would collide with it,
        # and _build_furniture's resolved is a dict, so a collision is a silent
        # overwrite.
        name = self._claim(name or "furn%d" % len(self.furniture), "furniture")
        self.furniture.append(dict(
            kind=kind, of=of, side=side, align=align, gap=float(gap),
            pad=float(pad), box=box, span=span, name=name, spec=spec))
        return self.furniture[-1]["name"]

    def _union(self, of, resolved):
        names = ["all"] if of is None else ([of] if isinstance(of, str) else list(of))
        boxes = []
        for n in names:
            if n == "all":
                boxes += [p.ink_bbox() for p in self.parts]
            elif n in resolved:
                boxes.append(resolved[n])
            else:
                boxes.append(self[n].ink_bbox())
        if not boxes:
            raise ValueError("of=%r resolved to nothing" % (of,))
        return [min(b[0] for b in boxes), min(b[1] for b in boxes),
                max(b[2] for b in boxes), max(b[3] for b in boxes)]

    def _build_furniture(self, colors):
        """(list of xml, list of boxes). In the parts' shared coordinates; build() translates everything at the end.

        `colors` is the merged colour table, appended in place: furniture [r,g,b] takes
        the same path as on the engine side
        (_add_graphics's _extend_palette), except there is no tool layer here, so it is wired up directly.
        """
        if not self.furniture:
            return [], []
        gid, tid = FURNITURE_GID, FURNITURE_TID
        out, boxes, resolved = [], [], {}

        def idx_of(rgb):
            ent = BT._RGB.format(*rgb)
            if ent not in colors:
                colors.append(ent)
            return colors.index(ent) + 2

        for f in self.furniture:
            anchor = self._union(f["of"], resolved)
            sp = self._union(f["span"], resolved) if f["span"] is not None else anchor
            side, pad, gap = f["side"], f["pad"], f["gap"]
            clean = dict(f["spec"], kind=f["kind"])
            for k in ("color", "fill"):
                if isinstance(clean.get(k), (list, tuple)) and len(clean[k]) == 3:
                    clean[k] = idx_of(list(clean[k]))

            if f["kind"] in ("rect", "line"):
                if side == "around":
                    w, h = anchor[2] - anchor[0] + 2 * pad, anchor[3] - anchor[1] + 2 * pad
                else:
                    bx = f["box"] or [None, None]
                    w = (sp[2] - sp[0]) if bx[0] is None else float(bx[0])
                    h = (sp[3] - sp[1]) if bx[1] is None else float(bx[1])
                bx0 = by0 = 0.0
            else:
                clean["at"] = [0.0, 0.0]
                _probe, pb = GX.emit([dict(clean)], origin=(0.0, 0.0))
                w, h = pb[2] - pb[0], pb[3] - pb[1]
                bx0, by0 = pb[0], pb[1]

            def across(lo, hi, size):
                return {"start": lo, "center": lo + (hi - lo - size) / 2.0,
                        "end": hi - size}[f["align"]]

            if side == "around":
                x, y = anchor[0] - pad, anchor[1] - pad
            elif side == "center":
                x = (anchor[0] + anchor[2] - w) / 2.0
                y = (anchor[1] + anchor[3] - h) / 2.0
            elif side == "above":
                x, y = across(sp[0], sp[2], w), anchor[1] - gap - h
            elif side == "below":
                x, y = across(sp[0], sp[2], w), anchor[3] + gap
            elif side == "left":
                x, y = anchor[0] - gap - w, across(sp[1], sp[3], h)
            else:                                            # right
                x, y = anchor[2] + gap, across(sp[1], sp[3], h)

            clean["at"] = ([x, y, x + w, y + h] if f["kind"] in ("rect", "line")
                           else [x - bx0, y - by0])
            xml, bb = GX.emit([clean], gid=gid, tid=tid, origin=(0.0, 0.0))
            gid += 20
            tid += 20
            out.append(xml)
            boxes.append(bb)
            resolved[f["name"]] = [x, y, x + w, y + h]
        return out, boxes

    def size(self):
        """Width and height of the whole figure's ink under the current placement (not
        written yet, pure measurement).

        To hand footnotes up to the parent layer you first need to know "how wide the
        parent is", and that waits until all parts are placed - so it takes two passes:
        place once to measure the width, then pass the width back to the part carrying
        the footnotes and redraw it.
        """
        bs = [p.ink_bbox() for p in self.parts]
        return (max(b[2] for b in bs) - min(b[0] for b in bs),
                max(b[3] for b in bs) - min(b[1] for b in bs))

    def write(self, path):
        # FIGDIR overrides the output directory. Every figure script writes under
        # ~/chemwand_figures, which is where people look at figures; regression tests
        # write to their own temp dir and should not flood that one.
        path = os.path.expanduser(path)
        _d = os.environ.get("FIGDIR")
        if _d:
            path = os.path.join(os.path.expanduser(_d), os.path.basename(path))
        doc, shift, size = self.build()
        # Check id uniqueness before going out. The composer gives every part a
        # non-overlapping id range, so a collision here means **some part is duplicated
        # internally** - measured: draw_catalytic_cycle drawing N decorated structures
        # in one call emits N <graphic>s all with id 7000, because decoration numbers
        # restart per structure. That is the engine's business, but the finished file is
        # mine, so it must be said. Count from **after** the <page ...> tag. The page
        # itself is id="1", and part 0 gets no range offset (stride × 0), so any engine
        # that starts with atom 1 - draw_molecule does - would trigger a false alarm
        # from this guard. A false alarm is worse than a miss: it makes real duplicates
        # look like noise.
        _body = doc[doc.index("<page"):]
        _ids = re.findall(r'\bid="(\d+)"', _body[_body.index(">") + 1:])
        _dup = sorted((i for i, n in Counter(_ids).items() if n > 1), key=int)
        if _dup:
            self.notes.append("⚠ duplicate ids %s in the finished figure - the id ranges are separate, so some part"
                              " is duplicated internally (decoration ids restart per structure)" % _dup)
        open(path, "w").write(doc)
        col = self.collisions()
        near = [r for r in self.collisions(pad=5.0)
                if not any(r["a"] == q["a"] and r["b"] == q["b"] for q in col)]
        # Served up in a separate column from collisions, not merged. Merged, the reader
        # would have to sort out which ones are
        # "deliberately underneath", and that is exactly what covered() has already sorted for them.
        cov = self.covered()
        return {"path": path, "near": near, "covered": cov,
                "dup_ids": _dup, "page": [round(size[0]), round(size[1])],
                "parts": [{"name": p.name, "objs": len(p.kept()),
                           "at": [round(p.dx, 1), round(p.dy, 1)],
                           "warnings": p.report.get("warnings") or [],
                           # The part's **own** notes must be served up too. The
                           # composer only looks across parts; problems inside a part
                           # are visible only to the engine's own instruments - #010's
                           # electrode sits on Int-B, and the engine's own words are
                           # "graphics[5] (rect) covers node:IntB (structure) by 7.0
                           # pt", but I forwarded only warnings, and so swallowed the
                           # one sentence that reported it.
                           "findings": _findings(p.report),
                           # Separate column from findings: findings is "a graphic
                           # covers a structure", drift is "the thing you asked for
                           # never happened". Two things, two actions.
                           "drift": _drift(p.report)}
                          for p in self.parts],
                "collisions": col, "conflicts": self.conflicts, "notes": self.notes}
