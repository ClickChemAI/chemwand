"""
beautify.py — coloring core for the desktop CDXML pipeline.

Per-atom / per-bond coloring (more accurate than whole-group coloring), plus
colortable extension. Operates on fragment or full
CDXML strings. Colors are palette indices; add them with extend_colortable first.

CDXML color facts: indices 0=black, 1=white are
implicit; the <colortable> lists user colors starting at index 2 (a fresh doc ships
idx2=white, idx3=black), so appended colors start at index 4. Carbons carry no atom
label, so a carbon is colored via its BONDS; heteroatoms via the node color.

NOTE: exact on-canvas color rendering is best confirmed once in ChemDraw (no
renderer here) - the structure of the emitted attributes matches ChemDraw's schema.
"""

from __future__ import annotations

import re

_RGB = "<color r=\"{0}\" g=\"{1}\" b=\"{2}\"/>"


def extend_colortable(cdxml: str, colors: list) -> tuple[str, list[int]]:
    """Append colors [[r,g,b], ...] (0..1 floats) to a full CDXML's <colortable>.
    Returns (cdxml, idx) where idx[i] is the ChemDraw color index of colors[i]."""
    m = re.search(r"<colortable>(.*?)</colortable>", cdxml, re.S)
    existing = len(re.findall(r"<color\b", m.group(1))) if m else 0
    start = 2 + existing
    add = "\n".join(_RGB.format(*c) for c in colors)
    out = cdxml.replace("</colortable>", add + "\n</colortable>", 1)
    return out, [start + i for i in range(len(colors))]


def _set_attr_on_open_tag(tag: str, attr: str, value) -> str:
    """Add or replace attr="value" within the OPENING tag (up to the first '>')."""
    end = tag.find(">")
    head, rest = tag[:end], tag[end:]
    if re.search(rf'\b{attr}="[^"]*"', head):
        head = re.sub(rf'\b{attr}="[^"]*"', f'{attr}="{value}"', head)
    else:
        head = re.sub(r"^(<\w+)", rf'\1 {attr}="{value}"', head, count=1)
    return head + rest


def color_nodes(frag: str, atom_ids, color_idx: int) -> str:
    """Set label color on the given atom nodes (by id). Heteroatoms show a colored
    symbol; bare carbons have no label (color them via color_bonds instead)."""
    ids = {str(i) for i in atom_ids}

    def repl(m):
        tag = m.group(0)
        nid = re.search(r'\bid="(\d+)"', tag)
        if not nid or nid.group(1) not in ids:
            return tag
        tag = _set_attr_on_open_tag(tag, "color", color_idx)
        # if the node carries a <t><s> label, color that run too
        tag = re.sub(r"<s\b([^>]*)>",
                     lambda sm: f"<s{_merge_color(sm.group(1), color_idx)}>", tag)
        return tag

    return _NODE_RE.sub(repl, frag)


def _merge_color(attrs: str, color_idx: int) -> str:
    if 'color="' in attrs:
        return re.sub(r'\bcolor="\d+"', f'color="{color_idx}"', attrs)
    return attrs + f' color="{color_idx}"'


def color_bonds(frag: str, color_idx: int, endpoints=None) -> str:
    """Color bonds. endpoints: a set of frozenset({idA, idB}) selecting specific
    bonds by their end-atom ids; None colors every bond in the fragment."""
    def repl(m):
        bond = m.group(0)
        if endpoints is not None:
            B = re.search(r'\bB="(\d+)"', bond)
            E = re.search(r'\bE="(\d+)"', bond)
            if not (B and E) or frozenset({B.group(1), E.group(1)}) not in endpoints:
                return bond
        if 'color="' in bond:
            return re.sub(r'\bcolor="\d+"', f'color="{color_idx}"', bond)
        return re.sub(r"\s*/>$", f' color="{color_idx}"/>', bond)

    return re.sub(r"<b\b[^>]*?/>", repl, frag)


def set_bond_attrs(frag: str, pairs, attrs: dict, *, orient: bool = False):
    """Set attributes on the bonds joining each ORDERED (idA, idB) pair.

    `pairs` is ordered, unlike color_bonds' frozensets, and that is the whole
    reason this exists. A dative bond has a DIRECTION -- the CDX spec defines
    Order="dative" as running from the begin atom to the end atom -- so which of
    the two ids sits in B= and which in E= is the difference between N->B and
    B->N, i.e. between the right arrow and a chemically wrong one.

    orient=True rewrites B=/E= so the bond runs idA -> idB. That is honouring the
    caller's SMARTS, not overriding it: the order in `pairs` came from the atom
    maps they wrote, while the order already in the file came from RDKit's bond
    direction, which is an artefact of parsing and means nothing chemically.
    Swapping the two ends of an ordinary bond changes nothing else -- region_bonds
    reads them as a set, and ColoredMolecularArea names bond ids, not endpoints.

    Bonds are matched on the id PAIR rather than on a bond id because that is what
    an anchor resolves to; a bond whose pair is not in the fragment is skipped
    silently here and reported by the caller, which knows the anchor.

    Returns (frag, replaced), where `replaced` is every attribute that already had
    a value and now has a different one: (bond id, attr, old, new). CDXML gives one
    `Display` and one `Order` per bond, so restyling can only ever be destructive
    when the bond was already saying something -- a wedge, or a double bond -- and
    the loss is invisible afterwards. The caller decides what to do about it; this
    function's job is to make sure nobody has to notice it on paper.
    """
    want = {frozenset({str(a), str(b)}): (str(a), str(b)) for a, b in pairs}
    replaced = []

    def repl(m):
        bond = m.group(0)
        B = re.search(r'\bB="(\d+)"', bond)
        E = re.search(r'\bE="(\d+)"', bond)
        if not (B and E):
            return bond
        key = frozenset({B.group(1), E.group(1)})
        if key not in want:
            return bond
        bid = re.search(r'\bid="(\d+)"', bond)
        for attr, value in attrs.items():
            was = re.search(rf'\b{attr}="([^"]*)"', bond)
            if was and was.group(1) != str(value):
                replaced.append((bid.group(1) if bid else None,
                                 attr, was.group(1), str(value)))
            bond = _set_attr_on_open_tag(bond, attr, value)
        if orient:
            a, b = want[key]
            bond = _set_attr_on_open_tag(bond, "B", a)
            bond = _set_attr_on_open_tag(bond, "E", b)
        return bond

    return re.sub(r"<b\b[^>]*?/>", repl, frag), replaced


def region_bonds(frag: str, atom_ids) -> list[tuple[str | None, str, str]]:
    """The bonds INSIDE a set of atoms: (bond id or None, endA, endB) for every
    bond whose BOTH endpoints are in the set.

    One definition of "the bonds of this region", used twice: color_region tints
    them, and decorate's "fill" names them in a ColoredMolecularArea. Two
    implementations of this predicate would drift the first time either changed,
    and they would drift silently -- a tint and a fill covering different bonds of
    the same matched ring is not something anyone would spot on paper.

    The id is None for a bond that carries none. Every bond cdxml_core.fragment()
    writes has had one since bonds needed to be nameable, but a fragment from the
    hand-drawn depiction library is whatever ChemDraw wrote.
    """
    ids = {str(i) for i in atom_ids}
    out = []
    for m in re.finditer(r"<b\b[^>]*?/>", frag):
        tag = m.group(0)
        B = re.search(r'\bB="(\d+)"', tag)
        E = re.search(r'\bE="(\d+)"', tag)
        if not (B and E) or B.group(1) not in ids or E.group(1) not in ids:
            continue
        bid = re.search(r'\bid="(\d+)"', tag)
        out.append((bid.group(1) if bid else None, B.group(1), E.group(1)))
    return out


# ONE node regex, used by color_region and color_nodes: a node is either
# self-closing or carries children (a <t> label, an Electron, ...).
_NODE_RE = re.compile(r"<n\b(?:[^>]*?/>|[^>]*?>.*?</n>)", re.S)


def color_region(frag: str, atom_ids, color_idx: int) -> str:
    """Color a connected region: the LABELLED atoms in the set AND every bond whose
    BOTH endpoints are in the set. Bare ring carbons show via those internal bonds.

    "Labelled" cannot mean `Element="`: that is every real atom that prints a
    symbol -- and NO placeholder. So an R group, or the `[B]` a paper uses as
    shorthand for a boron unit, could not be tinted by any spelling: not
    on="[#0]", not maps=[1]. Nothing would be raised, the label would just stay
    black.
    A placeholder node carries no Element (it is a GenericNickname or an
    ElementList) but it does carry the <t> that prints it, and color_nodes
    already tints a <t>'s runs -- so the test is whether the node DRAWS
    something, which is the only thing colouring an atom can mean.
    """
    ids = {str(i) for i in atom_ids}
    het = []
    for tag in _NODE_RE.findall(frag):
        nid = re.search(r'\bid="(\d+)"', tag)
        if nid and nid.group(1) in ids and ('Element="' in tag or "<t" in tag):
            het.append(nid.group(1))
    frag = color_nodes(frag, het, color_idx)
    eps = {frozenset({a, b}) for _bid, a, b in region_bonds(frag, atom_ids)}
    return color_bonds(frag, color_idx, endpoints=eps)
