"""compose.py — 多引擎组合原型。**一行引擎代码都没改。**

每个引擎照旧跑完整的一次调用（headless，产出一份完整 CDXML），组合器只在
成品文档这一层做事：切开 <page>、重编号、平移、并集、合表、重新落到正象限。

这正是 panels.py 的嵌套模型搬到文档层：每个引擎在自己的账本里量自己的 bbox，
组合器只学到一条"这块占了这么大"。引擎之间看不见对方的墨，所以组合器自己
补一道并集碰撞检查（report["collisions"]）。

    c = Composer()
    c.cycle("ring", nodes=[...], steps=[...])
    c.linear("top",  items=[...])
    c.above("top", "ring", gap=26)          # 按墨的实际范围叠
    c.join("tail", to="ring", at="II")      # 共享锚点：同一个物种只画一次
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

# 跨部件家具的号段：在任何单个部件用得到的号之上，又在跨部件箭头（4.9M）之下，
# 而 4.9M 本身在第 1 个部件的偏移 5M 之下。三者互不重叠。
FURNITURE_GID = 4_800_000
FURNITURE_TID = 4_850_000

# 每个部件分到的 id 段。引擎里最大的自动 id 是 _add_graphics 的 970000+si*20000，
# 20 个 slot 也才 1.37M，所以 5M 的步长留了一个量级的余量。
ID_STRIDE = 5_000_000

# renumber() 覆盖 6 个标量 + 4 个 id 列表，这几个它不覆盖，而 <scheme>/<step> 全靠它们。
# 漏掉的话反应步会指向别的部件的物件，文件照样打得开。
_EXTRA_ID_ATTRS = ("ReactionStepReactants", "ReactionStepProducts", "ReactionStepArrows",
                   "ReactionStepObjectsAboveArrow", "ReactionStepObjectsBelowArrow",
                   "ReactionStepAtomMap", "BondOrdering", "BracketedObjectIDs",
                   "GraphicID", "SupersededBy")
# AtomID 故意不在里面：它是化学上的原子编号（1,2,3…），不是物件 id。

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


# ---------------------------------------------------------------- 几何小工具

def _ovl(a, b):
    """两个 bbox 的重叠 (dx, dy)，非正表示不重叠。"""
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
    """部件报告里真正是"发现"的那几行（covers / 越界之类），不要那一大段说明。

    这一族是**部件内部**的图元压住结构（"graphics[5] (rect) covers node:IntB by
    7.0 pt"）。组合器的 covered() 故意不看部件内部——同来源的层序是那个引擎自己
    的账——所以这条通道不能撤，撤了 #010 那个电极就没人报了。
    """
    out = []
    for n in (report.get("notes") or []):
        for line in str(n).splitlines():
            if _FIND_RE.match(line) or " covers " in line:
                out.append(line.strip())
    return out


# 引擎会发 28 种 NOTE，_findings 一条抬头都放不过去（量过：90 个图脚本、712 个
# 部件、1310 行提醒，过滤器放过去 201 行，**全部**是上面那一族 graphics[N] covers
# cell:…）。也就是说活下来的恰恰是最该沉默的那一族，而真正该端上来的那一族一条都
# 没过——其中 11 条是"锚点一个都没匹配上，装饰直接没画"。
#
# 这张表挑的是一族：**你要了，它没发生**，或者**它悄悄换了别的做法**。剩下的不转
# 发，各有各的理由：
#   "carry placeholder labels"    395 条 / 72 张 —— 纯信息，Markush 图都有
#   "reaches outside the page"    103 条 / 38 张 —— 部件自己的页框，build() 最后
#                                                  会按并集重新框一次，丢掉了
#   "graphic(s) land on top"       75 条 / 27 张 —— covered() 按画序答得更准
#   "came through undeformed"       2 条        —— 信息，是好消息
#
# 挑法靠措辞，所以 test_every_forwarded_phrase_is_still_in_the_engine 把每一条钉在
# 引擎源码上：改了措辞先红，并且说得出是哪一条。
_DRIFT_PHRASES = (
    "highlight anchor(s) matched nothing",              # 11 条 / 5 张（这一轮）
    "each match names a GROUP rather than one atom",    # 4 条 / 3 张：一条 spec 画了好几个标记
    "'fill' anchor(s) matched something with no ring in it",   # 2 条：静默降级成 color
    "did not match the depiction template",             # 3 条：模板没匹配上，各画各的
    "came out of the constrained layout with uneven bonds",    # 2 条：模板被丢掉
    "carry no group in the abbreviation table",         # 4 条：abbreviate 开了没生效
    "is itself a Markush structure",                    # 1 条：模板本身带 R 基
    "were drawn in FULL anyway",                        # abbreviate 收下了没照做
    "were asked to float and could not",                # floating 收下了没照做
    "matched no placeholder in the SMILES",             # rlabels 的键打错了
    "had nowhere to land and were NOT written",         # AND/OR 立体标记没写进文件
    "did not load, so every molecule was computed",     # 手画图库没加载
    "named an id another cell already held",            # 同名 id：那个 SMILES 没画
    "could not be built and are NOT in the file",       # 整条目缺失
    "span(s) were not reached",                         # 跨列没生效
    "cannot be shrunk with them",                       # 试剂缩了，它身上的装饰没跟着
)


def _drift(report: dict) -> list:
    """"你要了，它没发生"那一族——引擎说了，而 _findings 一直扔掉的那些。"""
    out = []
    for n in (report.get("notes") or []):
        for line in str(n).splitlines():
            line = line.strip()
            if any(ph in line for ph in _DRIFT_PHRASES):
                out.append(line[6:].strip() if line.startswith("NOTE:") else line)
    return out


# ---------------------------------------------------------------- 部件

# ALIGN 有两套词，而调用者不该记得自己在跟哪一层说话。panels.py 说
# start/center/end；这里按轴说 left/right 或 top/bottom，因为"左"在竖着叠时
# 有意义、在横着排时没有。两套都收：start/end 按轴翻译。原有拼写一个字没动。
#
# 写错时列出的是**这条轴上**可用的词。不分轴地把四个都列出来，等于把
# "你用错轴了"这条信息扔掉。
_ALIGN_ALIAS = {"y": {"start": "left", "end": "right"},
                "x": {"start": "top", "end": "bottom"}}


def _pick(align, table, axis):
    key = _ALIGN_ALIAS[axis].get(align, align)
    if key not in table:
        raise ValueError(
            "align %r 不认识。沿这条轴可用：%s（也收 panels 那套 start/center/end）"
            % (align, " / ".join(sorted(table))))
    return table[key]


_AT_DEG = {"right": 0.0, "below": 90.0, "left": 180.0, "above": 270.0}


def _one_ring(part, i):
    """部件里第 i 个环，找不到就把量到了什么一起说出来。"""
    rr = part.rings()
    if len(rr) <= i:
        raise ValueError(
            "部件 %r 没有第 %d 个环：量到 %d 个（半径 %s）。只有 cycle 部件有环，"
            "而且弧得还在——join / drop 把弧丢掉了就量不到。"
            % (part.name, i, len(rr),
               ", ".join("%.1f" % x["r"] for x in rr) or "无"))
    return rr[i]


class Part:
    def __init__(self, name, doc, report):
        self.name, self.doc, self.report = name, doc, report
        self.dx = self.dy = 0.0
        self.drop = set()          # 要整段删掉的顶层对象下标

    def kept(self):
        return [o for i, o in enumerate(self.doc.objs) if i not in self.drop]

    def ink_bbox(self):
        bs = [o.ink_bbox() for o in self.kept() if o.ink_bbox()]
        return [min(b[0] for b in bs) + self.dx, min(b[1] for b in bs) + self.dy,
                max(b[2] for b in bs) + self.dx, max(b[3] for b in bs) + self.dy]

    def rings(self):
        """这个部件里每个环的圆：[{"c": (cx, cy), "r": R, "arcs": n}]，带上 dx/dy。

        取自**画出来的墨**，不是引擎承诺的坐标：环弧是 ChemDraw 的 <arrow>，圆心
        和半径就写在 Center3D / MajorAxisEnd3D 上（引擎自己也是拿这两个数画的）。
        报告里的 centre 是引擎布局空间的，跟成品 CDXML 差一个 finalize 的平移，
        而组合器摆位要的是后者。

        半径不猜：报告 rings[*].radius（单环图是 ring.radius）已经解出了 R，这里
        只认半径对得上那个数的弧组。喂料弧、挂件的弯箭头半径是 11~26，环是 95 起，
        本来就分得开；但"分得开"是这一份语料的巧合，"等于引擎解出的那个数"不是。
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
        # 一个环至少两段弧（一条边一段）。两根各自画的箭头不会把圆心撞到 0.1 单位
        # 以内，所以这一条同时挡掉了半径碰巧等于 R 的单根弧。
        out = [{"c": (cx + self.dx, cy + self.dy), "r": r, "arcs": n}
               for (cx, cy, r), n in groups.items() if n >= 2]
        out.sort(key=lambda d: (-d["r"], d["c"]))
        return out

    def body(self):
        spans = CD.top_objects(self.doc.body)
        return "\n".join(self.doc.body[s:e] for i, (_, s, e) in enumerate(spans)
                         if i not in self.drop)


# ---------------------------------------------------------------- 组合器

class Composer:
    def __init__(self, margin: float = 30.0, workdir: str | None = None):
        self.parts: list[Part] = []
        self.margin = margin
        self.workdir = workdir or os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "parts")
        self.notes: list[str] = []
        self.conflicts: list[dict] = []
        self.links: list[dict] = []    # 跨部件的箭头，见 arrow()
        # 跨部件的家具：通栏标题条、框住一组部件的外框、分区底色、段间横线。
        # 它们的锚点是**几个部件的并集**，那个数只有这一层知道——引擎那边
        # `of: "figure"` 里的 figure 是它自己那一次调用画出来的墨，不是整张图。
        self.furniture: list[dict] = []

    # -- 跑引擎 -------------------------------------------------------------
    # 名字是这一层唯一的地址。摆位、`of`、`span`、join、arrow —— 每一个都只收名
    # 字，而按名字找**只找得到第一个**。所以两个东西同名不是"后一个覆盖前一个"，
    # 是后一个从此**再也指不到**：#054 里三只灯泡和三个供体都叫 b1/b2/b3，摆位全
    # 落在灯泡上，三个供体留在原点，write() 一声不吭。
    #
    # 拒掉而不是报一句，因为没有可退的正确行为：挑第一个是随意的，挑最后一个一样
    # 随意，而调用者要的那个，两种挑法都可能不是。改名是一秒钟的事。
    #
    # 部件和家具共用这一个命名空间——`of=` 两边都收，`_union` 先查家具再查部件，
    # 所以一件家具叫了部件的名字就把那个部件**遮掉**了。一条规则管两边。
    def _claim(self, name, what):
        """占下一个名字，占过了就拒。`what` 只用来把错误说清楚。"""
        name = str(name)
        parts = [p.name for p in self.parts]
        furn = [f["name"] for f in self.furniture]
        held = "部件" if name in parts else ("家具" if name in furn else None)
        if held:
            raise ValueError(
                "名字 %r 已经被一个%s占了，不能再给%s —— 这一层按名字找东西只找得到"
                "第一个，后一个会摆不动。现有：部件 %s；家具 %s"
                % (name, held, what, ", ".join(parts) or "无", ", ".join(furn) or "无"))
        return name

    def _run(self, tool, name, kw):
        import chemwand_mcp as M
        import json
        # 先占名字再调引擎：重名的话画出来也是白画，而且报错要早于那一趟布局。
        name = self._claim(name, "部件")
        fn = getattr(M, tool)
        with M.headless(self.workdir):
            rep = fn(filename="_part_" + name, **kw)
        rep = json.loads(rep) if isinstance(rep, str) else rep
        path = rep.get("path")
        if not path or not os.path.exists(path):
            raise RuntimeError("%s(%s) 没写出文件：%s" % (tool, name, str(rep)[:400]))
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
        raise KeyError("没有名为 %r 的部件（有 %s）"
                       % (name, ", ".join(p.name for p in self.parts)))

    # -- 摆位：按墨的实际范围 ------------------------------------------------
    def _stack(self, who, ref, gap, axis, sign, align):
        """ref 可以是一个名字，也可以是一串名字 —— 一串的话按它们的**并集**摆。
        顶栏要压在两个环上面，只挑其中一个当参照会歪。"""
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

    # -- 摆位：按圆，不按外框 ------------------------------------------------
    #
    # 为什么单开一条路，而不是给 below / right_of 调一个 gap：那几个按**墨的外框**
    # 摆，而一个环部件的外框是"环 + 挂在环上的每个分子和每个标号"。外框刚好贴上的
    # 时候，两个圆还隔着挂件伸出去的那一截。量过已经画好的 17 张多环图：**一张都没
    # 相切**，最紧的 050 差 0.3 根键长，016 差 2.2 根（31.7 单位）。引擎那边的默认
    # 摆法是同一个毛病（layout_cycle 的 `_x += 2*(R+reach_max) + center_gap`），
    # 连单字母节点的两个环都会被推开 37 单位。
    #
    # 而语料里这 17 张的原图，14 张是**相切**的，切点上正压着 SET / PT / HAT 徽章
    # —— 切点就是电子转移画在哪儿。所以相切不是外观偏好，是化学逻辑的位置。
    #
    # 相切没有阈值可调：两圆心距 == R1 + R2，一个等式。`gap` 留着不是为了微调，是
    # 为了让原图本来就相交的那两张（118 / 120，小 Ir 环穿过大 Cu 环，共用的那步电
    # 子转移就在透镜里）能写成负数。
    def tangent(self, who, ref, at="right", gap=0.0, ring=0, ref_ring=0):
        """把 who 挪到跟 ref 正好相切。返回 (切点, 实测缺口)。

        at  who 的圆心落在 ref 圆心的哪个方向。给度数（CDXML 的 y 朝下，所以
            0=右、90=下、180=左、270=上），或者给 right / below / left / above。
        gap 两圆之间额外的距离。0 = 真相切；负 = 故意相交。
        """
        a, b = self[who], self[ref]
        ca, cb = _one_ring(a, ring), _one_ring(b, ref_ring)
        if isinstance(at, str):
            if at not in _AT_DEG:
                raise ValueError("at %r 不认识，可用 %s，或者直接给度数"
                                 % (at, " / ".join(_AT_DEG)))
            deg = _AT_DEG[at]
        else:
            deg = float(at)
        th = math.radians(deg)
        d = ca["r"] + cb["r"] + gap
        a.dx += cb["c"][0] + d * math.cos(th) - ca["c"][0]
        a.dy += cb["c"][1] + d * math.sin(th) - ca["c"][1]
        # 切点：从 ref 的圆心往那个方向走 ref 自己的半径。gap=0 时它同时在两个圆上，
        # 徽章就该压在这儿。
        pt = (cb["c"][0] + cb["r"] * math.cos(th),
              cb["c"][1] + cb["r"] * math.sin(th))
        got = self.ring_gap(who, ref, ring=ring, ref_ring=ref_ring)
        self.notes.append("%s 相切于 %s：方向 %.0f°，切点 (%.1f, %.1f)，"
                          "实测缺口 %+.2f" % (who, ref, deg, pt[0], pt[1], got))
        return pt, got

    def ring_gap(self, who, ref, ring=0, ref_ring=0):
        """两个部件的环之间量出来的缺口：圆心距 - (R1+R2)。

        0 = 相切，正 = 散开，负 = 相交。摆完位拿它复核——相切必须是一条能跑的
        检查，不是"我调过了"。
        """
        ca, cb = _one_ring(self[who], ring), _one_ring(self[ref], ref_ring)
        d = math.hypot(ca["c"][0] - cb["c"][0], ca["c"][1] - cb["c"][1])
        return d - (ca["r"] + cb["r"])


    # -- 摆位：共享锚点 ------------------------------------------------------
    def anchor_of(self, part, key, by):
        """(锚点坐标, [要丢掉的对象]) —— by='species' 用结构框中心，by='label' 用标签框中心。"""
        p = self[part]
        if by == "label":
            lab = p.doc.label(key)
            return ((lab.bbox[0] + lab.bbox[2]) / 2, (lab.bbox[1] + lab.bbox[3]) / 2), [lab]
        f, lab = p.doc.species(key)
        if f is None:
            raise ValueError("%s: 标签 %r 底下找不到结构（text-only 的物种请用 by='label'）"
                             % (part, key))
        return ((f.bbox[0] + f.bbox[2]) / 2, (f.bbox[1] + f.bbox[3]) / 2), [f, lab]

    def residual(self, who, to, at, my=None, by="species", by_to=None):
        """两边都认这个物种时，它们各自落点还差多少。

        两个引擎各自排完版，同一个物种在两份产出里的相对位置是各自定的。锚定一处
        就把另一处的差额固定下来了——这个差额没人能同时消掉，它是**两个引擎的
        约束冲突本身**，不是 bug。组合器能做的是把它量出来报给人看。
        """
        a, b = self[who], self[to]
        ca, _ = self.anchor_of(who, my or at, by)
        cb, _ = self.anchor_of(to, at, by_to or by)
        return (round(cb[0] + b.dx - (ca[0] + a.dx), 1),
                round(cb[1] + b.dy - (ca[1] + a.dy), 1))

    def residual_shapes(self, who, to, drop=False):
        """两个部件已经摆好之后，列出所有共有结构还差多少（不移动任何东西）。

        `#098` 那种：两个环共用一个**节点**（文字），又共用一个**物种**（结构）。
        锚点只能挑一个，另一个的差额得报出来。drop=True 时顺手把 who 里重复的
        结构删掉——不管对不对得上，同一个分子画两遍总是错的。
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
            self.notes.append("  ⚠ %s/%s 还共用一个结构，差 (%.1f, %.1f)" % (who, to, r[0], r[1]))
            if drop:
                f = sa[k][0]
                lab = a.doc.label_of(f)
                for i, o in enumerate(a.doc.objs):
                    if o is f or (lab is not None and o is lab):
                        a.drop.add(i)
        return out

    def join_shape(self, who, to, on=0, drop=True):
        """把 who 和 to 里**画的是同一个分子**的结构对上，不靠标签。

        找出两边共有的全部物种，用第 `on` 个对齐，其余的把残差量出来报给人看——
        一条两端都接在环上的链（#026 的 imidoyl→…→adduct）必然有一头对不上，
        那是两个引擎各自排完版之后的约束冲突，不是可以选一个"对"的答案。
        drop=True 时 who 里那些重复的结构整段删掉，共享的物种只画一次。
        """
        a, b = self[who], self[to]
        sa, sb = a.doc.shapes(), b.doc.shapes()
        common = [k for k in sa if k in sb and len(sa[k]) == 1 and len(sb[k]) == 1]
        if not common:
            raise ValueError("%s 和 %s 没有画同一个分子" % (who, to))
        common.sort(key=lambda k: sa[k][0].bbox[1])       # 从上到下，好指认
        pick = common[on]
        fa, fb = sa[pick][0], sb[pick][0]
        a.dx += (fb.bbox[0] + fb.bbox[2]) / 2 + b.dx - (fa.bbox[0] + fa.bbox[2]) / 2
        a.dy += (fb.bbox[1] + fb.bbox[3]) / 2 + b.dy - (fa.bbox[1] + fa.bbox[3]) / 2
        self.notes.append("%s 按共有结构对齐到 %s（共有 %d 个，用第 %d 个）"
                          % (who, to, len(common), on))
        for i, k in enumerate(common):
            if i != on:
                r = (round(sb[k][0].bbox[0] + b.dx - sa[k][0].bbox[0] - a.dx, 1),
                     round(sb[k][0].bbox[1] + b.dy - sa[k][0].bbox[1] - a.dy, 1))
                self.notes.append("  ⚠ 第 %d 个共有结构对不上：差 (%.1f, %.1f)" % (i, r[0], r[1]))
                self.conflicts.append({"part": who, "vs": to, "key": "shape#%d" % i,
                                       "gap": list(r)})
        if drop:
            for k in common:
                f = sa[k][0]
                lab = a.doc.label_of(f)        # 标签必须一起丢，否则编号留在原地
                for i, o in enumerate(a.doc.objs):
                    if o is f or (lab is not None and o is lab):
                        a.drop.add(i)
        return a

    def join(self, who, to, at, my=None, drop=True, by="species",
             by_to=None, also=None):
        """把 who 里叫 `my`（默认同名 `at`）的物种，对齐到 to 里叫 `at` 的物种上。

        对齐的是**结构的 bbox 中心**，不是标签的位置——两个引擎给标签的行高和
        对齐方式不一定一样，结构本身才是两边都认得的同一个东西。
        drop=True 时把 who 里那一份重复的结构和标签整段删掉，于是这个物种在
        成品里只画一次，跨引擎的箭头也就不用画了。
        """
        a, b = self[who], self[to]
        my = my or at
        ca, mine = self.anchor_of(who, my, by)
        cb, theirs = self.anchor_of(to, at, by_to or by)
        a.dx += cb[0] - ca[0] + b.dx
        a.dy += cb[1] - ca[1] + b.dy
        self.notes.append("锚定 %s.%s → %s.%s，位移 (%.1f, %.1f)"
                          % (who, my, to, at, a.dx, a.dy))
        if drop:
            for i, o in enumerate(a.doc.objs):
                if o in mine:
                    a.drop.add(i)
            self.notes.append("  丢掉 %s 里重复的 %r，只保留 %s 那一份" % (who, my, to))
        if also:
            k, kb = (also, also) if isinstance(also, str) else also
            r = self.residual(who, to, kb, my=k, by=by, by_to=by_to)
            self.notes.append("  ⚠ 第二个锚点 %r 对不上：差 (%.1f, %.1f)。"
                              "两个引擎各自排完版，这个差没人能同时消掉" % (k, r[0], r[1]))
            self.conflicts.append({"part": who, "vs": to, "key": k, "gap": list(r)})
        return a

    def arrow_between(self, part, a, b):
        """part 里连接物种 a 和 b 的那根箭头。

        找法是几何的：横跨在两个结构之间、竖直方向与它们重叠的那根。不走引擎的
        `placement` 出参，因为 linear 的 placement 只报箭头的 tail/head，**一个物种的
        几何都没有**（layout_linear.py:1000-1030），对不上名字。
        """
        p = self[part]
        fa, _ = p.doc.species(a)
        fb, _ = p.doc.species(b)
        if fa is None or fb is None:
            raise ValueError("%s: %r 或 %r 找不到结构" % (part, a, b))
        lo, hi = sorted((fa.bbox[2], fb.bbox[0])) if fa.bbox[2] <= fb.bbox[0] \
            else sorted((fb.bbox[2], fa.bbox[0]))
        ymid = (fa.bbox[1] + fa.bbox[3] + fb.bbox[1] + fb.bbox[3]) / 4
        hits = [o for o in p.kept() if o.tag == "arrow" and o.ink_bbox()
                and lo - 4 <= (o.ink_bbox()[0] + o.ink_bbox()[2]) / 2 <= hi + 4
                and o.ink_bbox()[1] - 30 <= ymid <= o.ink_bbox()[3] + 30]
        if len(hits) != 1:
            raise ValueError("%s: %r→%r 之间找到 %d 根箭头" % (part, a, b, len(hits)))
        return hits[0]

    def under_arrow(self, who, ref, a, b, gap=0.0, align="center"):
        """把 who 整块摆到 ref 里 a→b 那根箭头的下面，横向对中。

        这是**没有共享物种**时的摆法：催化剂循环和主链一个分子都不共用，
        它们之间只有"这个循环解释的是这一步"这层关系，而那是几何关系。
        gap 可以为负——原图里环的上弧就是顶到箭头上的。
        """
        w, r = self[who], self[ref]
        ar = self.arrow_between(ref, a, b)
        ab = ar.ink_bbox()
        wb = w.ink_bbox()
        cx = (ab[0] + ab[2]) / 2 + r.dx
        w.dx += cx - (wb[0] + wb[2]) / 2 if align == "center" else 0.0
        w.dy += (ab[3] + r.dy + gap) - wb[1]
        self.notes.append("%s 摆在 %s 的 %s→%s 箭头下方（间距 %.0f）" % (who, ref, a, b, gap))
        return w

    # -- 跨部件的箭头：父层唯一能画的那一根 -----------------------------------

    def arrow(self, frm, to, *, bulge=0.0, gap=10.0, head="Full", tail=None,
              color=None, line=None):
        """一头在这一块、另一头在那一块的箭头。

        引擎画不了它，而且这不是遗漏：每个引擎只知道**自己**排了什么，两端分属
        两次调用的那根箭头，两个端点的位置只有摆完位的父层同时知道。#173 的原图
        就是这个形状 —— `•CCl₃` 从循环里出来，直接接到主链 8-3→8-4 那一步上；
        今天只能在环里和链里各画一份，读者得自己认出那是同一个物种。

        端点写成 (部件名, 指认)：
          "8-3"            那个标签指认的结构
          ("8-2", "8-3")   这两个结构之间的那根箭头（reagent 落在一步上就是这个）

        解析推迟到 build()：摆位还可能变（under_arrow、nudge 都在后面），
        端点得按**最终**位置算，而不是声明那一刻的。
        """
        self.links.append(dict(frm=frm, to=to, bulge=float(bulge), gap=float(gap),
                               head=head, tail=tail, color=color, line=line))
        return self

    def _link_box(self, ref):
        """(部件名, 指认) -> 这一块在拼合坐标里的墨框。"""
        name, what = ref
        p = self[name]
        if isinstance(what, (tuple, list)):
            o = self.arrow_between(name, *what)
        else:
            # species() 返回 (结构, 标签)。**写成文字的物种结构是 None** —— 环上
            # 的 `•CCl₃ + Br⁻` 就是这样，它只有一个 <t>。那也是一个合法的端点：
            # 一根箭头指向"写出来的那个物种"和指向画出来的没有区别。
            frag, lab = p.doc.species(str(what))
            o = frag or lab
        if o is None:
            raise ValueError("%s 里找不到 %r" % (name, what))
        b = o.ink_bbox() or o.bbox
        return [b[0] + p.dx, b[1] + p.dy, b[2] + p.dx, b[3] + p.dy]

    @staticmethod
    def _edge(box, toward, gap):
        """从 box 的中心朝 toward 走，落在框边上再退开 gap。

        端点取边而不是中心：箭头从结构中心长出来会压在它自己的键上。
        """
        cx, cy = (box[0] + box[2]) / 2, (box[1] + box[3]) / 2
        dx, dy = toward[0] - cx, toward[1] - cy
        n = math.hypot(dx, dy) or 1.0
        hw, hh = (box[2] - box[0]) / 2, (box[3] - box[1]) / 2
        # 缩放到矩形边界：取两个方向上先撞到的那一个
        t = min(hw / abs(dx) if dx else 1e18, hh / abs(dy) if dy else 1e18)
        return (cx + dx * t + dx / n * gap, cy + dy * t + dy / n * gap)

    def _build_links(self, aid0=4_900_000):
        """把声明过的跨部件箭头解析成 (xml, 折线)。在 build() 里调用。

        号段 4.9M：在任何单个部件用得到的号之上（最大是 _add_graphics 的
        970000+si*20000），又在第 1 个部件的偏移 5M 之下。
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

    # -- 跨部件的箭头：父层唯一能画的那一根 -----------------------------------

    def arrow(self, frm, to, *, bulge=0.0, gap=10.0, head="Full", tail=None,
              color=None, line=None):
        """一头在这一块、另一头在那一块的箭头。

        引擎画不了它，而且这不是遗漏：每个引擎只知道**自己**排了什么，两端分属
        两次调用的那根箭头，两个端点的位置只有摆完位的父层同时知道。#173 的原图
        就是这个形状 —— `•CCl₃` 从循环里出来，直接接到主链 8-3→8-4 那一步上；
        今天只能在环里和链里各画一份，读者得自己认出那是同一个物种。

        端点写成 (部件名, 指认)：
          "8-3"            那个标签指认的结构
          ("8-2", "8-3")   这两个结构之间的那根箭头（reagent 落在一步上就是这个）

        解析推迟到 build()：摆位还可能变（under_arrow、nudge 都在后面），
        端点得按**最终**位置算，而不是声明那一刻的。
        """
        self.links.append(dict(frm=frm, to=to, bulge=float(bulge), gap=float(gap),
                               head=head, tail=tail, color=color, line=line))
        return self

    def _link_box(self, ref):
        """(部件名, 指认) -> 这一块在拼合坐标里的墨框。"""
        name, what = ref
        p = self[name]
        if isinstance(what, (tuple, list)):
            o = self.arrow_between(name, *what)
        else:
            # species() 返回 (结构, 标签)。**写成文字的物种结构是 None** —— 环上
            # 的 `•CCl₃ + Br⁻` 就是这样，它只有一个 <t>。那也是一个合法的端点：
            # 一根箭头指向"写出来的那个物种"和指向画出来的没有区别。
            frag, lab = p.doc.species(str(what))
            o = frag or lab
        if o is None:
            raise ValueError("%s 里找不到 %r" % (name, what))
        b = o.ink_bbox() or o.bbox
        return [b[0] + p.dx, b[1] + p.dy, b[2] + p.dx, b[3] + p.dy]

    @staticmethod
    def _edge(box, toward, gap):
        """从 box 的中心朝 toward 走，落在框边上再退开 gap。

        端点取边而不是中心：箭头从结构中心长出来会压在它自己的键上。
        """
        cx, cy = (box[0] + box[2]) / 2, (box[1] + box[3]) / 2
        dx, dy = toward[0] - cx, toward[1] - cy
        n = math.hypot(dx, dy) or 1.0
        hw, hh = (box[2] - box[0]) / 2, (box[3] - box[1]) / 2
        # 缩放到矩形边界：取两个方向上先撞到的那一个
        t = min(hw / abs(dx) if dx else 1e18, hh / abs(dy) if dy else 1e18)
        return (cx + dx * t + dx / n * gap, cy + dy * t + dy / n * gap)

    def _build_links(self, aid0=4_900_000):
        """把声明过的跨部件箭头解析成 (xml, 折线)。在 build() 里调用。

        号段 4.9M：在任何单个部件用得到的号之上（最大是 _add_graphics 的
        970000+si*20000），又在第 1 个部件的偏移 5M 之下。
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

    # -- 碰撞：父层唯一知道两个引擎同时占了哪里的地方 -------------------------
    def collisions(self, tol: float = 1.0, pad: float = 0.0):
        """pad>0 时把每个物件先胀 pad 再判，用来抓"没重叠但贴太近"。

        重叠和留白是两个问题——allocator 自己也是这么分的（重叠用 depth，
        留白用 `probe(pad=...)`）。#190 里环的 5CzBN 节点标签离链条的箭头只有
        12 个单位，一处重叠都没有，但读起来像是链条的一部分。"""
        objs = []
        # 父层自己画的箭头也得过尺子。它是这张图里唯一不属于任何引擎的墨，
        # 漏掉它等于在一个"专为查重叠而存在"的图上放一个不查的东西。
        class _Link:                      # 只为凑齐 tag/text 两个字段
            tag = "arrow"
            def __init__(self, i): self.i = i
            def text(self): return "父层箭头 %d" % self.i
        for i, (_xml, pts) in enumerate(self._build_links()):
            xs = [q[0] for q in pts]; ys = [q[1] for q in pts]
            objs.append(("(父层)", _Link(i),
                         [min(xs) - pad, min(ys) - pad, max(xs) + pad, max(ys) + pad],
                         [] if pad else [(pts[k][0], pts[k][1], pts[k + 1][0], pts[k + 1][1])
                                         for k in range(len(pts) - 1)]))
        for p in self.parts:
            for o in p.kept():
                bb = o.ink_bbox()
                if not bb:
                    continue
                # 只有空白的 <t> 什么都不画。它们来自 label:" " 那个技巧——
                # 引擎要求每个物种有标签，不想印就得给一个空格，于是文件里留下
                # 一个有框但没有墨的对象。拿它判碰撞是纯粹的假阳性。
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
                    continue                     # 引擎内部的重叠是引擎自己的账
                dx, dy = _ovl(ba, bb)
                if dx <= 0 or dy <= 0:
                    continue
                if sa and sb:
                    hit = any(_seg_seg(s, t) for s in sa for t in sb)
                    kind, sev = "线×线", 0.0
                elif sa or sb:
                    segs, box = (sa, bb) if sa else (sb, ba)
                    hit = any(_seg_box(s, box) for s in segs)
                    kind, sev = "线×块", 0.0
                else:
                    hit = min(dx, dy) > tol
                    kind, sev = "块×块", min(dx, dy)
                if hit:
                    out.append({"a": "%s/%s %s" % (pa, oa.tag, oa.text()[:18]),
                                "b": "%s/%s %s" % (pb, ob.tag, ob.text()[:18]),
                                "kind": kind, "depth": round(sev, 1),
                                "where": [round(max(ba[0], bb[0])), round(max(ba[1], bb[1]))]})
        return sorted(out, key=lambda r: -r["depth"])

    # -- 谁在上面：画序 + 不透明 ---------------------------------------------
    #
    # 实测 11 种 RectangleType（scratchpad/opaque_probe，文字先发、矩形后发）：
    # 名字里带 Filled 或 Shaded 的把底下的字整个吃掉；Shadow / Shadow Bold /
    # RoundEdge Shadow 只有边框，肚子是透的，字照样读得出来。圆盘（Orbital 的
    # OvalType）同一条规则，多边形看 CurveType 的 128 位（1=Closed, 128=Filled）。
    _OPAQUE = re.compile(r"Filled|Shaded")

    @staticmethod
    def _opaque(o) -> bool:
        """这个对象会不会把它之前画的东西整个盖掉。"""
        if Composer._OPAQUE.search(o.a.get("RectangleType")
                                   or o.a.get("OvalType") or ""):
            return True
        if o.tag == "curve":
            try:
                return bool(int(o.a.get("CurveType", "0")) & 128)
            except ValueError:
                return False
        return False

    # 一个不透明物件遮住的是它**自己的形状**，不是它的外框。语料里 26 条初判
    # 有 20 条是团块（多边形），按外框算，团块四角外面的字也被判成"盖住了"——
    # #090 那个 "a" 就在左上角，离灰团还有一截。
    _COVER_GRID = 9          # 受害者框上取 9x9=81 个点，分辨率约 1.2%

    @staticmethod
    def _region(o):
        """不透明物件的 (外框, 判"点在里面吗")。None 表示读不出形状。

        外框不读 BoundingBox：ChemDraw 的 Orbital / Symbol 把它写成"右边先、
        高度为零"（cdxml_core 里那段手画文件的注释记着这件事），当成矩形用会
        得到一个倒着的框，于是那个物件永远判不出重叠——静默漏报。Center3D /
        MajorAxisEnd3D / MinorAxisEnd3D 是几何量，矩形和圆盘都带。
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
        """受害者框里有多大一块落在那个形状内。

        取样不取外框交集：结构的框大半是白的，团块又不是方的，两头都会算歪。
        在受害者框上铺格子、逐点问"在不在里面"，对字（框就是一块实心的墨）是
        准的，对结构仍然偏大——一个分子的框里本来就有空地。宁可偏大：漏报一条
        看不见的标签，比多报一条要贵。
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
        """家具的每个顶层对象 + 它属于哪件家具，在部件的共享坐标系里。

        _build_furniture 一件家具只还一个并集框，而一件带 fill 的 text 家具发的是
        "实心矩形 + 文字"两个对象——按并集算，那块矩形会把它自己的文字判成被盖住。
        所以这里把 xml 再拆一遍，用的还是读部件那套 top_objects / Obj。

        色表传空的：这一趟只量几何，颜色索引从 2 还是从 9 开始不影响任何一个坐标。
        """
        if not self.furniture:
            return []
        xmls, _boxes = self._build_furniture([])
        out = []
        for f, x in zip(self.furniture, xmls):
            for t, s, e in CD.top_objects(x):
                out.append((f["name"], CD.Obj(t, x[s:e])))
        return out

    # 默认 0.10，Ge 从 90 个图脚本跑出来的阶梯上挑的：
    #     >= 1%  38 条/16 张    >= 10%  23 条/11 张    >= 50%  6 条/5 张
    # 定这个数之前先定的是"跟什么比"——跟受害者自己的面积比（见 hid 那段）。取值
    # 落在 10% 是因为**验过的最低那一档就在这里**：#090 面板 a/b/d 的标题分别报
    # 22% / 11% / 10%，三条在真图上全是被灰团从中间啃掉的标题，肉眼一看就不对。
    # 再往下那一档（1%~7%，14 条，多是装饰球擦到标签边）没有逐条看过真图，所以
    # 不放进默认值里。要看就显式传 covered(hid=0.01)。
    COVERED_HID = 0.10

    def covered(self, hid: float | None = None):
        """谁被**后**画的实心块盖住了，盖掉了它自己的百分之几。

        和 collisions() 是两个问题。collisions() 问"两样东西撞上了吗"，这里问
        "撞上的那两样，哪个在上面"。图里有一整类重叠是**故意的**——蛋白团块垫在
        分子底下、底色条垫在文字后面、徽章压在框上——把它们算成碰撞，一张图能报
        出几十处，真事故就淹在里面了。分它们的不是什么标记，是**画的先后**：先画
        的那块是背景，后画的那块才是事故。

        画序在 build() 里是定死的，所以"谁在上面"不用猜：家具（按加入顺序）→ 部件
        （按声明顺序，各自的 Z 依次往上摞）→ 父层箭头。没家具的图不重盖 Z，那时
        部件之间照旧按各自的 Z，顺序还是声明顺序。

        `hid` 跟**受害者自己的面积**比，不跟盖它的那块比，也不用绝对 pt：一条
        400x60 的底色条盖掉一个 30x10 的标签，重叠深度只有 10 pt 听着不多，而那个
        标签是 100% 没了。collisions() 报 depth 是因为它问的是"扎进去多深"。

        只认两种受害者：文字（<t>）和结构（<fragment>）。一根被埋掉的线或箭头这里
        不报——引擎自己的 _reconcile 在部件内部已经报它，而跨部件的线本来就少。
        """
        hid = self.COVERED_HID if hid is None else hid
        objs = []                       # 严格按画序：(来源, 对象, 框, dx, dy)
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
                continue                         # 形状读不出来就不猜
            if ox or oy:                         # 部件里的物件，属性还在部件坐标系
                box = [box[0] + ox, box[1] + oy, box[2] + ox, box[3] + oy]
                inside = (lambda f, a, b: lambda x, y: f(x - a, y - b))(inside, ox, oy)
            for sa, oa, ba, _ax, _ay in objs[:j]:      # 只看它**之前**画的
                if sa == sb:
                    continue                     # 同来源内部的层序是那个引擎的账
                if oa.tag not in ("t", "fragment"):
                    continue
                if oa.tag == "t" and not oa.text().strip():
                    continue                     # label:" " 那个占位技巧，没有墨
                dx, dy = _ovl(ba, box)
                if dx <= 0 or dy <= 0:
                    continue                     # 外框先粗筛，省掉绝大多数取样
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

    # -- 出片 ---------------------------------------------------------------
    def build(self):
        if not self.parts:
            raise ValueError("一个部件都没有")
        head = self.parts[0].doc.header
        fonts = {re.search(r'id="(\d+)"', f).group(1): f
                 for f in _table(head, "fonttable")}
        colors = _table(head, "colortable")

        bodies, zbase, idbase = [], 0, 0
        for k, p in enumerate(self.parts):
            body = p.body()
            if k:
                body = _renumber_extra(CORE.renumber(body, idbase), idbase)
            # 下一个部件的号段从本部件实际用到的最大 id 之上开始。固定的 STRIDE × k
            # 假设每个部件都在一段之内；一个已经合并过的文件占了好几段，再合一次时
            # 下一个部件就落进它的第二段，撞号（3 张合成的再并 1 张，25 个重号）。
            # 普通部件都在一段之内，算出来仍是 STRIDE × k，输出不变。
            top = max((int(i) for i in re.findall(r'\bid="(\d+)"', body)), default=0)
            idbase = max(idbase, (top // ID_STRIDE + 1) * ID_STRIDE)
            if k:
                # 字体表：按 id 求并，同 id 不同字体要说出来
                for f in _table(p.doc.header, "fonttable"):
                    fid = re.search(r'id="(\d+)"', f).group(1)
                    if fid in fonts and fonts[fid] != f:
                        self.notes.append("字体 id=%s 两个部件定义不同，用了第一个" % fid)
                    fonts.setdefault(fid, f)
                # 颜色表：公共前缀之外的追加，并把本部件的索引改写过去
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
                    self.notes.append("%s 的颜色索引重映射 %s" % (p.name, remap))
            # Z：各引擎都从 1 开始重编，直接拼会让层序错乱。按声明顺序往上摞。
            if p.doc.max_z():
                body = _shift_z(body, zbase)
                zbase += p.doc.max_z()
            bodies.append(CORE._translate_body(body, p.dx, p.dy))

        # 家具在**合表之后、重写表头之前**建：它的 [r,g,b] 要就地追加进 colors，
        # 而 colors 下面几行就要写回表头去。发在所有部件之前，配合后面的 Z 重盖，
        # 才真的在墨的后面。
        furn_xml, furn_boxes = self._build_furniture(colors)
        bodies = furn_xml + bodies

        head = re.sub(r'<fonttable>.*?</fonttable>',
                      "<fonttable>\n%s\n</fonttable>" % "\n".join(
                          fonts[k] for k in sorted(fonts, key=int)), head, flags=re.S)
        head = re.sub(r'<colortable>.*?</colortable>',
                      "<colortable>\n%s\n</colortable>" % "\n".join(colors),
                      head, flags=re.S)

        # 父层自己画的那几根箭头。放在所有部件之后发，所以它压在上面 ——
        # 它连接的是两个部件，被任何一个盖住都读不出来。
        links = self._build_links()
        bodies += [xml for xml, _pts in links]

        bs = [p.ink_bbox() for p in self.parts]
        cb = [min(b[0] for b in bs), min(b[1] for b in bs),
              max(b[2] for b in bs), max(b[3] for b in bs)]
        for _xml, _pts in links:      # 页框要装得下它，否则 ChemDraw 会裁掉
            cb[0] = min(cb[0], min(q[0] for q in _pts))
            cb[1] = min(cb[1], min(q[1] for q in _pts))
            cb[2] = max(cb[2], max(q[0] for q in _pts))
            cb[3] = max(cb[3], max(q[1] for q in _pts))
        for _b in furn_boxes:         # 同理：通栏条比墨宽，框比墨大
            if _b:
                cb[0], cb[1] = min(cb[0], _b[0]), min(cb[1], _b[1])
                cb[2], cb[3] = max(cb[2], _b[2]), max(cb[3], _b[3])
        # 一次正象限平移 + 一个页框。引擎各自 finalize 过一遍，拼完又有负坐标，
        # 而 ChemDraw 打开文件时把页原点钉在 (0,0) 并裁掉负的那部分。
        ox, oy = self.margin - cb[0], self.margin - cb[1]
        body = CORE._translate_body("\n".join(bodies), ox, oy)
        if furn_xml:
            # DOCUMENT ORDER IS NOT DRAWING ORDER —— cdxml_core.stamp_z 用三个
            # 被光晕盖住的圆盘证过这件事。部件各自带着自己的 Z 进来，_shift_z
            # 把它们摞起来；家具要压在所有墨的后面，光靠"发得早"办不到。
            # 所以整篇 strip 一次、按文档顺序重盖一次，和 _add_graphics 在单个
            # 引擎里做的是同一件事，同样**只在有家具时做**——没家具的图一个字节
            # 都不动，否则每张现有图的 golden 都要重录。
            body = re.sub(r'\s+Z="\d+"', "", body)
            body, _znext = CORE.stamp_z(body)
        pw, ph = cb[2] - cb[0] + 2 * self.margin, cb[3] - cb[1] + 2 * self.margin
        return CORE.page(head, body, [0, 0, round(pw), round(ph)]), (ox, oy), (pw, ph)

    # -- 跨部件家具 ---------------------------------------------------------
    #
    # 为什么这一层非有不可：引擎那边的 `of: "figure"` 指的是**那一次调用画出来
    # 的墨**。126 顶上那条深灰标题条要横跨整张图，它既不属于任何一次引擎调用，
    # 也不是箭头——`_build_links` 之外，组合器原来没有第二样自己会画的东西。
    #
    # 和引擎那边 `of`/`side` 用同一套词汇（above/below/left/right/center/around），
    # 因为它们对读的人是同一件事。**唯一的差别写在这里**：组合层没有账本，所以
    # 不搜索、不避让，按锚点算出来就落在那儿。真撞上了由 collisions() 报——那本来
    # 就是"父层唯一知道两个引擎同时占了哪里的地方"。

    _FURN_SIDES = ("above", "below", "left", "right", "center", "around")

    def graphic(self, kind, of="all", side="above", *, align="center",
                gap=8.0, pad=6.0, box=None, span=None, name=None, **spec):
        """一件跨部件家具。`of` / `span` 收部件名、家具名、名字列表，或 "all"。

        `of` 定位，`span` 定另一根轴的尺寸（不给就等于 `of`）。两者分开是因为
        126 的通栏条要"压在 top 上面"但"和整张图一样宽"——一个锚点答不了。

        `box` 只给放置的 rect / line；某一轴写 None = 取 `span` 在那一轴的尺寸。
        """
        if side not in self._FURN_SIDES:
            raise ValueError("side 是 %s 之一，收到 %r"
                             % ("、".join(self._FURN_SIDES), side))
        # 自动名 furn%d 也过这一关：调用者手写一个 "furn3" 之后，第 4 件家具的自动
        # 名就会撞上它，而 _build_furniture 的 resolved 是个 dict，撞上就是静默覆盖。
        name = self._claim(name or "furn%d" % len(self.furniture), "家具")
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
            raise ValueError("of=%r 没解析出任何东西" % (of,))
        return [min(b[0] for b in boxes), min(b[1] for b in boxes),
                max(b[2] for b in boxes), max(b[3] for b in boxes)]

    def _build_furniture(self, colors):
        """(xml 列表, 盒子列表)。在部件的共享坐标系里，build() 最后统一平移。

        `colors` 是合并后的色表，就地追加：家具的 [r,g,b] 和引擎那边走同一条路
        （_add_graphics 的 _extend_palette），只是这里没有 tool layer，所以自己接。
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
        """当前摆位下，整张图的墨有多宽多高（还没落盘，纯量）。

        脚注要上交父层就得先知道"父层有多宽"，而那要等所有部件都摆好——
        所以是两趟：先摆一趟量出宽度，再把宽度回传给那个背脚注的部件重画一次。
        """
        bs = [p.ink_bbox() for p in self.parts]
        return (max(b[2] for b in bs) - min(b[0] for b in bs),
                max(b[3] for b in bs) - min(b[1] for b in bs))

    def write(self, path):
        # FIGDIR 改写落盘目录。每个图脚本都写在 ~/chemwand_figures 下，那是
        # 人看图的地方；回归测试要写到自己的临时目录去，不该往那里灌文件。
        path = os.path.expanduser(path)
        _d = os.environ.get("FIGDIR")
        if _d:
            path = os.path.join(os.path.expanduser(_d), os.path.basename(path))
        doc, shift, size = self.build()
        # 出门前查一遍 id 唯一。组合器给每个部件分了互不重叠的号段，所以这里
        # 再撞上就说明**某个部件自己内部**就重了 —— 实测 draw_catalytic_cycle
        # 一次画 N 个带装饰的结构，会发出 N 个 id 都等于 7000 的 <graphic>，
        # 装饰号是按结构重新开始数的。那是引擎的账，但成品文件是我出的，得说出来。
        # 从 <page ...> 这一标签**之后**开始数。页面自己就是 id="1"，而第 0 个部件
        # 不加号段偏移（stride × 0），所以任何以 1 号原子开头的引擎——draw_molecule
        # 就是——会被这条守卫误报一次。误报比漏报更糟：它会让真的重号看起来像噪音。
        _body = doc[doc.index("<page"):]
        _ids = re.findall(r'\bid="(\d+)"', _body[_body.index(">") + 1:])
        _dup = sorted((i for i, n in Counter(_ids).items() if n > 1), key=int)
        if _dup:
            self.notes.append("⚠ 成品里有重复 id %s —— 号段是分开的，所以是某个部件"
                              "内部自己重了（装饰 id 按结构重数）" % _dup)
        open(path, "w").write(doc)
        col = self.collisions()
        near = [r for r in self.collisions(pad=5.0)
                if not any(r["a"] == q["a"] and r["b"] == q["b"] for q in col)]
        # 和 collisions 分成两栏端上来，不合并。合起来读的人就得自己分辨哪几条是
        # "故意垫在下面"，而那正是 covered() 已经替他分完的那件事。
        cov = self.covered()
        return {"path": path, "near": near, "covered": cov,
                "dup_ids": _dup, "page": [round(size[0]), round(size[1])],
                "parts": [{"name": p.name, "objs": len(p.kept()),
                           "at": [round(p.dx, 1), round(p.dy, 1)],
                           "warnings": p.report.get("warnings") or [],
                           # 部件**自己**的 notes 也要端上来。组合器只看跨部件，
                           # 部件内部的问题只有引擎自己的仪表看得见——#010 的电极
                           # 压在 Int-B 上，引擎原话是 "graphics[5] (rect) covers
                           # node:IntB (structure) by 7.0 pt"，而我只转了 warnings，
                           # 就把唯一报出来的那句吞了。
                           "findings": _findings(p.report),
                           # 和 findings 分两栏：findings 是"图元压住了结构"，
                           # drift 是"你要的这件事根本没发生"。两件事，两个动作。
                           "drift": _drift(p.report)}
                          for p in self.parts],
                "collisions": col, "conflicts": self.conflicts, "notes": self.notes}
