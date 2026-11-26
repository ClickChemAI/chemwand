"""cdxdoc.py — 只读地把一份 CDXML 拆成顶层对象，并算出每个对象真正的"墨".

组合器的地基。三件事：
  1) 精确切分顶层对象（深度扫描，不是正则猜）——因为要能整段删掉一个对象
  2) 读几何：fragment/t 用 BoundingBox；arrow(弧)/curve 用折线，因为弧的 bbox
     是外接矩形，拿它判碰撞会把整个环的内部都算成"占用"
  3) 找命名锚点：印在图上的标签文字 -> 它底下那个结构
"""
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import cdxml_core as CORE     # noqa: E402
import text_metrics as TM     # noqa: E402

# ---------------------------------------------------------------- 顶层对象切分

_ATTR = re.compile(r'\s+([\w:]+)="([^"]*)"')


def _end_of_tag(s: str, i: int) -> int:
    """i 指向 '<'，返回该开标签 '>' 的下标。跳过引号内的 '>'。"""
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
    """[(tag, start, end)]，end 为开区间右端。只列深度 0 的对象。"""
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
            # 自闭合不改变深度
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


# ---------------------------------------------------------------- 几何

def _segs(p):
    return [(p[i][0], p[i][1], p[i + 1][0], p[i + 1][1]) for i in range(len(p) - 1)]


_RUN = re.compile(r'<s\b([^>]*)>(.*?)</s>', re.S)
_UNESC = (("&lt;", "<"), ("&gt;", ">"), ("&amp;", "&"))


def text_envelope(xml: str, a: dict):
    """引擎写出来的 <t> **没有 BoundingBox**（CORE.text 只写 p= 和 Justification）。

    后果有两条，都很实：页框求并时漏掉全部文字，于是页太紧、ChemDraw 在页边把
    文字折行（就是拼接实验里 tBuOOH 竖成两行那件事）；碰撞检查也会整类漏掉文字。
    所以这里自己量。

    量法对着 ChemDraw 自己回写的 13 个 t 标定过（cycle_010.cdxml，Name= 已经
    变成文件名，是 ChemDraw 存回去的，那些 BoundingBox 是它自己的真值）：
      · p.x 就是盒子中心（Center 时，|Δ| < 0.5）
      · p.y 是最后一行的基线，盒子往下还有降部
      · 仓库的 text_metrics 宽度一律比 ChemDraw 宽 0.29~1.45（方向是对的，够用不会短）
      · 高度在有上下标时比 ChemDraw 少 2.0~3.6，所以有上下标就补 0.45*size
    这是个**外包络**，不是精确框：宁可略大。
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
    """结构里那些**长原子标签**的包络，对称估。

    `<fragment>` 的 BoundingBox **不覆盖原子标签**：#190 的 rlabel 是
    `TT^+^BF~4~^−^`，标签 p.x=31.59 而 fragment 左边界是 30.0，导出来第一个 T 被裁掉。

    为什么对称估而不是按方向算：原子标签用的是 `LabelJustification`（不是顶层
    caption 的 `Justification`），而且 ChemDraw 对 GenericNickname 的 `p` 不是
    文字原点——实测 `LabelJustification="Left"` 的标签画在 p 的**左边**。
    这一条我没有拿 ChemDraw 的真值标定过，所以不猜方向，两边各留半个宽度。
    页框方向上这只是多留白；碰撞方向上是"宁可多报"。
    """
    out = []
    for m in re.finditer(r'<t\b([^>]*)>(.*?)</t>', xml, re.S):
        ta = dict(re.findall(r'(\w+)="([^"]*)"', m.group(1)))
        if "p" not in ta:
            continue
        ta.pop("Justification", None)          # 一律按居中的对称包络算
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
        """判碰撞用的形状：折线段列表 [(x0,y0,x1,y1)]，或 None 表示用 bbox。"""
        a = self.a
        if self.tag == "arrow":
            # CORE.arc_points / CORE.curve_polyline 是这两种几何的规范读入口，
            # 它们的 docstring 明写"每个渲染器都必须调它，不要各自再推一遍"。
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
        """平移无关的结构指纹：认"这两个部件画的是同一个分子"。

        比按标签认稳得多——共享的物种常常根本没有印出来的标签（#026 的
        imine / adduct 都没有），而两个部件里的同一个 SMILES 是同一份 layout
        代码画的，连坐标都一样，差的只是一个平移。

        指纹取**整段 fragment**，不是挑几个属性。挑着取会出事：亚胺和亚胺酰
        自由基的原子数、坐标、键级、连 <n> 上的属性都一模一样，差别在自由基那个
        标记上，而它不在 <n> 里。整段一起算，差别在哪儿都跑不掉。
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
    """一份引擎产出的完整 CDXML。"""

    def __init__(self, src: str, name: str = "?"):
        self.name = name
        m = re.search(r'(<page\b[^>]*BoundingBox="([-\d. ]+)"[^>]*>)\n?(.*)\n?</page>',
                      src, re.S)
        if not m:
            raise ValueError("%s: 找不到 <page>" % name)
        self.header = src[:m.start(1)]
        self.page_bbox = _nums(m.group(2))
        self.body = m.group(3)
        self.objs = [Obj(t, self.body[s:e]) for t, s, e in top_objects(self.body)]

    # -- 几何 --------------------------------------------------------------
    def ink_bbox(self):
        bs = [o.ink_bbox() for o in self.objs if o.ink_bbox()]
        return [min(b[0] for b in bs), min(b[1] for b in bs),
                max(b[2] for b in bs), max(b[3] for b in bs)]

    def max_id(self) -> int:
        return max([int(v) for v in re.findall(r'\bid="(\d+)"', self.body)] or [0])

    def max_z(self) -> int:
        return max([int(v) for v in re.findall(r'\bZ="(\d+)"', self.body)] or [0])

    # -- 锚点 --------------------------------------------------------------
    def label(self, txt: str):
        """印着 txt 的那个 <t>。同名多个则抛——组合里的重名就是歧义。"""
        hits = [o for o in self.objs if o.tag == "t" and o.text() == txt]
        if len(hits) != 1:
            raise KeyError("%s: 标签 %r 命中 %d 个（%s）"
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
        """结构正下方那个标签（x 有重叠、离得最近）。`species()` 的反向。

        丢重复结构的时候必须连标签一起丢：只丢结构的话，编号会留在原地，
        同一个物种在成品里出现两个号。
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
            d = o.bbox[1] - frag.bbox[3]        # 标签顶 - 结构底
            if -4 <= d < best:
                best, bo = d, o
        return bo

    def species(self, txt: str):
        """标签 txt 指认的那个结构：标签正上方、x 有重叠、离得最近的 fragment。"""
        lab = self.label(txt)
        lx = (lab.bbox[0] + lab.bbox[2]) / 2 if lab.bbox else _nums(lab.a["p"])[0]
        ly = lab.bbox[1] if lab.bbox else _nums(lab.a["p"])[1]
        best, bo = 1e18, None
        for o in self.objs:
            if o.tag != "fragment" or not o.bbox:
                continue
            if not (o.bbox[0] - 10 <= lx <= o.bbox[2] + 10):
                continue
            d = ly - o.bbox[3]          # 标签顶 - 结构底
            if -4 <= d < best:
                best, bo = d, o
        return bo, lab
