"""cdxdoc.py - read-only split of one CDXML into top-level objects, plus each object's
true "ink".

Foundation of the composer. Three jobs:
  1) Exact split into top-level objects (depth scan, not regex guessing), because a
     whole object must be removable in one piece
  2) Read geometry: fragment/t use BoundingBox; arrow(arc)/curve use a polyline,
     because an arc's bbox is its bounding rectangle, and using it for collisions
     would mark the whole inside of the ring as "occupied"
  3) Find named anchors: label text printed on the figure -> the structure beneath it
"""
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import cdxml_core as CORE     # noqa: E402
import text_metrics as TM     # noqa: E402

# ---------------------------------------------------------------- top-level object split

_ATTR = re.compile(r'\s+([\w:]+)="([^"]*)"')


def _end_of_tag(s: str, i: int) -> int:
    """i points at '<'; return the index of that open tag's '>'. Skips '>' inside
    quotes."""
    q = None
    j = i + 1
    while j < len(s):
        c = s[j]
        if q:
            if c == q:
                q = None
        elif c in '"\'':
            q = c
        elif c == '>':
            return j
        j += 1
    raise ValueError("unterminated tag at %d" % i)


def top_objects(body: str):
    """[(tag, start, end)], end is the exclusive right end. Lists only depth-0
    objects."""
    out, depth, i, start, tag0 = [], 0, 0, None, None
    while True:
        i = body.find('<', i)
        if i < 0:
            break
        if body.startswith('<?', i) or body.startswith('<!', i):
            i = body.find('>', i) + 1
            continue
        gt = _end_of_tag(body, i)
        closing = body[i + 1] == '/'
        selfclose = body[gt - 1] == '/'
        name = re.match(r'</?([\w:]+)', body[i:gt + 1]).group(1)
        if closing:
            depth -= 1
            if depth == 0:
                out.append((tag0, start, gt + 1))
                start = tag0 = None
        elif selfclose:
            if depth == 0:
                out.append((name, i, gt + 1))
            # a self-closing tag does not change depth
        else:
            if depth == 0:
                start, tag0 = i, name
            depth += 1
        i = gt + 1
    return out


def attrs(xml: str) -> dict:
    gt = _end_of_tag(xml, 0)
    return dict(_ATTR.findall(xml[:gt + 1]))


def _nums(v):
    return [float(t) for t in v.split()]


# ---------------------------------------------------------------- geometry

def _segs(p):
    return [(p[i][0], p[i][1], p[i + 1][0], p[i + 1][1]) for i in range(len(p) - 1)]


_RUN = re.compile(r'<s\b([^>]*)>(.*?)</s>', re.S)
_UNESC = (("&lt;", "<"), ("&gt;", ">"), ("&amp;", "&"))


def text_envelope(xml: str, a: dict):
    """A <t> written by the engine has **no BoundingBox** (CORE.text writes only p= and
    Justification).

    Two consequences, both real: the page-box union misses all text, so the page is too
    tight and ChemDraw wraps text at the page edge (the tBuOOH split into two lines in
    the stitching experiment); collision checks also miss text as a whole class. So we
    measure it here.

    Calibrated against 13 t objects that ChemDraw itself wrote back (cycle_010.cdxml;
    Name= had already become the file name, i.e. ChemDraw saved it, so those
    BoundingBoxes are its own ground truth):
      · p.x is the box centre (for Center, |Δ| < 0.5)
      · p.y is the baseline of the last line; the box extends further down by the
        descender
      · the repo's text_metrics width is always 0.29~1.45 wider than ChemDraw (right
        direction: enough, never short)
      · with sub/superscripts the height is 2.0~3.6 less than ChemDraw, so add
        0.45*size when they are present
    This is an **outer envelope**, not an exact box: better slightly too big.
    """
    lines, cur, has_script, size = [], [], False, 9.0
    for at, t in _RUN.findall(xml):
        d = dict(re.findall(r'(\w+)="([^"]*)"', at))
        if t == "\n":
            lines.append(cur); cur = []; continue
        f, sz = int(d.get("face", 0)), float(d.get("size", 9))
        size = max(size, sz)
        script = bool(f & 32 or f & 64)
        has_script |= script
        for k, v in _UNESC:
            t = t.replace(k, v)
        cur.append((t, sz * TM.SUBSCRIPT_W if script else sz, bool(f & 1), bool(f & 2)))
    lines.append(cur)
    if not any(lines):
        return None
    w = max(sum(TM.width(t, e, bold=b, italic=i, trim=False) for t, e, b, i in L) or 0.0
            for L in lines)
    up = TM.height(["".join(t for t, *_ in L) for L in lines], size,
                   bold=any(r[2] for L in lines for r in L))
    if has_script:
        up += 0.45 * size
    dn = 0.25 * size
    px, py = [float(v) for v in a["p"].split()[:2]]
    just = a.get("Justification", "Center")
    x0 = px - w / 2 if just == "Center" else (px - w if just == "Right" else px)
    return [x0, py - up, x0 + w, py + dn]


_NODE = re.compile(r'<n\b([^>]*)>(.*?)</n>', re.S)


def frag_label_boxes(xml: str):
    """Envelope of the **long atom labels** in a structure, estimated symmetrically.

    The `<fragment>` BoundingBox **does not cover atom labels**: #190's rlabel is
    `TT^+^BF~4~^−^`, label p.x=31.59 while the fragment's left edge is 30.0, and the
    first T was clipped on export.

    Why symmetric rather than direction-aware: atom labels use `LabelJustification` (not
    the top-level caption's `Justification`), and ChemDraw's `p` for a GenericNickname
    is not the text origin - measured: a `LabelJustification="Left"` label is drawn to
    the **left** of p. I have not calibrated this against ChemDraw ground truth, so no
    guessing the direction: leave half a width on each side. For the page box this only
    adds whitespace; for collisions it means "rather over-report".
    """
    out = []
    for m in re.finditer(r'<t\b([^>]*)>(.*?)</t>', xml, re.S):
        ta = dict(re.findall(r'(\w+)="([^"]*)"', m.group(1)))
        if "p" not in ta:
            continue
        ta.pop("Justification", None)          # always use the centred, symmetric envelope
        env = text_envelope(m.group(0), ta)
        if env:
            out.append(env)
    return out


class Obj:
    __slots__ = ("tag", "xml", "a", "bbox", "ink")

    def __init__(self, tag, xml):
        self.tag, self.xml = tag, xml
        self.a = attrs(xml)
        bb = self.a.get("BoundingBox")
        self.bbox = _nums(bb) if bb else None
        if self.bbox is None and tag == "t" and "p" in self.a:
            self.bbox = text_envelope(xml, self.a)
        if tag == "fragment" and self.bbox:
            for e in frag_label_boxes(xml):
                self.bbox = [min(self.bbox[0], e[0]), min(self.bbox[1], e[1]),
                             max(self.bbox[2], e[2]), max(self.bbox[3], e[3])]
        self.ink = self._ink()

    def _ink(self):
        """Shape for collision checks: list of polyline segments [(x0,y0,x1,y1)], or
        None to use the bbox."""
        a = self.a
        if self.tag == "arrow":
            # CORE.arc_points / CORE.curve_polyline are the canonical readers for these
            # two geometries; their docstrings say plainly "every renderer must call
            # this; do not re-derive it on your own".
            p = CORE.arc_points(a)
            if p:
                return _segs(p)
            if "Head3D" in a and "Tail3D" in a:
                h, t = _nums(a["Head3D"])[:2], _nums(a["Tail3D"])[:2]
                return [(t[0], t[1], h[0], h[1])]
        if self.tag == "curve" and "CurvePoints" in a:
            v = _nums(a["CurvePoints"])
            return _segs(CORE.curve_polyline([(v[i], v[i + 1]) for i in range(0, len(v), 2)]))
        return None

    def ink_bbox(self):
        if self.ink:
            xs = [v for s in self.ink for v in (s[0], s[2])]
            ys = [v for s in self.ink for v in (s[1], s[3])]
            return [min(xs), min(ys), max(xs), max(ys)]
        return self.bbox

    def shape(self):
        """Translation-invariant structure fingerprint: recognises "these two parts draw
        the same molecule".

        Much more robust than matching by label - shared species often have no printed
        label at all (#026's imine / adduct have none), and the same SMILES in two parts
        is drawn by the same layout code, with identical coordinates; the only
        difference is a translation.

        The fingerprint takes the **whole fragment**, not a few chosen attributes.
        Choosing breaks: an imine and an imidoyl radical have identical atom counts,
        coordinates, bond orders, even the attributes on <n>; the difference is the
        radical marker, and that is not in <n>. Hash the whole thing and the difference
        cannot escape, wherever it is.
        """
        if self.tag != "fragment" or "BoundingBox" not in self.a:
            return None
        bb = _nums(self.a["BoundingBox"])
        x = CORE._translate_body(self.xml, -bb[0], -bb[1])
        x = re.sub(r'\b(id|B|E|BeginAttach|EndAttach|object)="\d+"', "", x)
        x = re.sub(r'\b(Attachments|CrossingBonds|BondCircularOrdering|BasisObjects)'
                   r'="[\d\s]+"', "", x)
        x = re.sub(r'\bZ="\d+"', "", x)
        return re.sub(r"\s+", " ", x)

    def text(self) -> str:
        return "".join(re.findall(r'<s\b[^>]*>(.*?)</s>', self.xml, re.S)).strip()


class Doc:
    """A complete CDXML produced by an engine."""

    def __init__(self, src: str, name: str = "?"):
        self.name = name
        m = re.search(r'(<page\b[^>]*BoundingBox="([-\d. ]+)"[^>]*>)\n?(.*)\n?</page>',
                      src, re.S)
        if not m:
            raise ValueError("%s: no <page> found" % name)
        self.header = src[:m.start(1)]
        self.page_bbox = _nums(m.group(2))
        self.body = m.group(3)
        self.objs = [Obj(t, self.body[s:e]) for t, s, e in top_objects(self.body)]

    # -- geometry ----------------------------------------------------------
    def ink_bbox(self):
        bs = [o.ink_bbox() for o in self.objs if o.ink_bbox()]
        return [min(b[0] for b in bs), min(b[1] for b in bs),
                max(b[2] for b in bs), max(b[3] for b in bs)]

    def max_id(self) -> int:
        return max([int(v) for v in re.findall(r'\bid="(\d+)"', self.body)] or [0])

    def max_z(self) -> int:
        return max([int(v) for v in re.findall(r'\bZ="(\d+)"', self.body)] or [0])

    # -- anchors -----------------------------------------------------------
    def label(self, txt: str):
        """The <t> that prints txt. Raises if several share the name - a duplicate name
        in a composition is an ambiguity."""
        hits = [o for o in self.objs if o.tag == "t" and o.text() == txt]
        if len(hits) != 1:
            raise KeyError("%s: label %r matched %d (%s)"
                           % (self.name, txt, len(hits),
                              ", ".join(sorted({o.text() for o in self.objs
                                                if o.tag == "t"})[:14])))
        return hits[0]

    def shapes(self):
        out = {}
        for o in self.objs:
            k = o.shape()
            if k:
                out.setdefault(k, []).append(o)
        return out

    def label_of(self, frag):
        """The label directly below the structure (x overlaps, closest). Inverse of
        `species()`.

        When dropping a duplicate structure its label must go too: drop only the
        structure and the number stays behind, so one species ends up with two numbers
        in the finished figure.
        """
        if not frag.bbox:
            return None
        cx = (frag.bbox[0] + frag.bbox[2]) / 2
        best, bo = 1e18, None
        for o in self.objs:
            if o.tag != "t" or not o.bbox or not o.text().strip():
                continue
            if not (frag.bbox[0] - 10 <= (o.bbox[0] + o.bbox[2]) / 2 <= frag.bbox[2] + 10):
                continue
            d = o.bbox[1] - frag.bbox[3]        # label top - structure bottom
            if -4 <= d < best:
                best, bo = d, o
        return bo

    def species(self, txt: str):
        """The structure that label txt designates: the fragment directly above the
        label, x overlapping, closest."""
        lab = self.label(txt)
        lx = (lab.bbox[0] + lab.bbox[2]) / 2 if lab.bbox else _nums(lab.a["p"])[0]
        ly = lab.bbox[1] if lab.bbox else _nums(lab.a["p"])[1]
        best, bo = 1e18, None
        for o in self.objs:
            if o.tag != "fragment" or not o.bbox:
                continue
            if not (o.bbox[0] - 10 <= lx <= o.bbox[2] + 10):
                continue
            d = ly - o.bbox[3]          # label top - structure bottom
            if -4 <= d < best:
                best, bo = d, o
        return bo, lab
