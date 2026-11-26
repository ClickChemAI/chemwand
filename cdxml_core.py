"""
cdxml_core.py — desktop CDXML generation core (offline, no web renderer).

Port of the ChemDraw+ Web library's PART 1 build helpers (the `B.*` family)
plus a fragment writer, so molecule fragments and reaction schemes can be
assembled as CDXML *strings* and opened ONCE in desktop ChemDraw.

Design note (render-once):
  This module does zero rendering and needs no chemistry toolkit. It works on
  coordinates + bonds that are already computed. The SMILES -> 2D coordinates
  step is separate (RDKit on the Mac) and feeds Atom/Bond lists into here.

Coordinate system matches ChemDraw CDXML: BondLength = 14.4 units, y grows down.
Fragments are normalized to a (0,0) top-left origin so bbox is [0,0,w,h] and
`place()` can drop them at an absolute page position (mirrors the web engine's
`dim(id) -> [0,0,w,h]` assumption).
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass

import text_metrics as M

BOND_LENGTH = 14.4
SUPERATOM_LABEL_SIZE = 10   # pt size of a contracted superatom <t> label (see _superatom_node)
ATOM_LABEL_SIZE = 10        # matches LabelSize in the ACS header below
ATOM_LABEL_MARGIN = 1.59    # matches MarginWidth in the same header: the space
                            # ChemDraw keeps around a label so bonds do not
                            # touch the glyphs. It is part of what gets drawn,
                            # so it belongs in the box.

# How far BELOW a node's point the label's <t> baseline sits. Measured, not
# guessed: every labelled node in a hand-drawn Markush reference file (a ChemDraw
# export -- two R, one X, one Y, one R^1) puts its <t> at exactly node.y + 3.58, with no
# spread at all. _superatom_node still writes 3.60, the value it was built with;
# the 0.02 is inside the line width and correcting it would rewrite every golden
# that contains an abbreviation, so that is a separate, deliberate regen.
LABEL_BASELINE_DROP = 3.58

# How far into a ring a FLOATING substituent bond reaches, as a fraction of the
# ring's INRADIUS -- centre to the middle of an edge, not centre to a vertex.
# The inradius is the right yardstick because the bond crosses an EDGE, and that
# edge is also what a reader looks at to decide whether it went through.
#
# Window, measured off the three floating bonds drawn by hand in a ChemDraw
# reference file: 0.698, 0.729, 0.729 (the last two are one object
# copied, so this is really two independent judgements, 0.698 and 0.729). 0.71 is
# between them. In the same units the bond crosses the edge and carries on
# 0.27-0.30 of an inradius past it, and all three bonds are exactly 14.40 long --
# an ordinary bond, moved, which is what makes this cosmetic rather than special.
FLOAT_DEPTH = 0.71

# ----------------------------------------------------------------------
# Decoration constants, all measured off a reference file in which the
# primitives were drawn by hand in ChemDraw (see disc() for why they are orbitals).
# ----------------------------------------------------------------------
# A dot marking one atom, as a fraction of the bond length. Window: the four
# hand-drawn discs came out at 0.577, 0.600, 0.600 and 0.689 -- the last one is
# the odd one, sitting a full bond away from any atom rather than on one, so the
# atom-marking window is really (0.577, 0.600) and 0.59 is inside it.
ATOM_DOT_R = 0.59

# A ring position number. Two measurements, both from the same six digits:
#   SIZE 6 pt, against LabelSize 10 for the atoms -- a locant is deliberately
#   smaller than anything it labels.
#   OFFSET from its vertex, radially outward, to the digit's ink centre:
#   4.43 4.43 4.49 5.02 5.04 5.49. The tight end is where a substituent crowded
#   the digit sideways (the "1" also came out 13.8 degrees off radial), so the
#   window is honest rather than tight, and 4.8 sits in the middle of it.
LOCANT_SIZE = 6
LOCANT_OFFSET = 4.8

# Ordinary valences, for working out how many hydrogens ChemDraw will draw on a
# heteroatom. Charge shifts it the obvious way: N+ carries four bonds, O- one.
_VALENCE = {"B": 3, "C": 4, "N": 3, "O": 2, "P": 3, "S": 2, "Se": 2,
            "F": 1, "Cl": 1, "Br": 1, "I": 1, "Si": 4}


def _implicit_h(element: str, charge: int, bond_order_sum: int) -> int:
    """Hydrogens ChemDraw will render on this atom. 0 for anything unusual."""
    v = _VALENCE.get(element)
    if v is None:
        return 0
    return max(0, v + charge - bond_order_sum)


def _label_overhang(element: str, nh: int):
    """(along, half_sym, half_h) for a rendered atom label, in drawing units.

    `along` is how far the label reaches on the side the hydrogens hang off;
    `half_sym` how far the bare symbol reaches on the other side. The symbol
    sits ON the atom point, so both are measured from there.
    """
    half_sym = M.width(element, ATOM_LABEL_SIZE) / 2.0
    h = 0.0
    if nh >= 1:
        h += M.width("H", ATOM_LABEL_SIZE)
    if nh >= 2:
        # SCALED by SUBSCRIPT_W. SUBSCRIPT_W comes from a caption's manual
        # subscript (face=32) while an atom label is set with LabelFace="96",
        # so the two could in principle differ. Measured against ChemDraw's own
        # 652 sample files, the two agree:
        # face=96 gives k = 0.74888 over 3377 runs, face=32 gives 0.7525,
        # face=64 gives 0.7524 -- all within 0.4%. One constant covers both.
        #
        # Measured as "H2" minus "H" rather than as "2" alone. The shrink rule
        # is POSITIONAL -- a digit shrinks only when something precedes it --
        # so measuring the digit in isolation would report it at full size and
        # silently undo the fix. In the rendered label it follows the H.
        h += (M.width("H" + str(nh), ATOM_LABEL_SIZE, formula=True)
              - M.width("H", ATOM_LABEL_SIZE))
    half_h = M.height1(element + ("H" if nh else ""), ATOM_LABEL_SIZE) / 2.0
    m = ATOM_LABEL_MARGIN
    return half_sym + h + m, half_sym + m, half_h + m

# ----------------------------------------------------------------------
# Header (colortable / fonttable / <CDXML> attributes)
# ----------------------------------------------------------------------
# Two ways to get a header:
#   * header_prefix_from(real_cdxml)  -> RELIABLE. Borrow the exact prefix from a
#     ChemDraw-exported CDXML (its own colortable + fonttable). Same trick the web
#     library uses via B.page(headerFrom, ...). Guarantees the doc opens and pins
#     the real fonttable (which the metrics model must match).
#   * minimal_header_prefix()         -> BEST EFFORT synthetic header for when no
#     real CDXML is on hand. Openability in desktop ChemDraw is UNVERIFIED until
#     tested on the Mac. font id 21 is declared as Arial here as a starting guess.
_MIN_HEADER = (
    '<?xml version="1.0" encoding="UTF-8" ?>\n'
    '<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >\n'
    '<CDXML CreationProgram="ChemWand 1.0" BondLength="14.400000"'
    ' LabelFont="21" CaptionFont="21" LabelSize="10" CaptionSize="10" LabelFace="96"'
    ' HashSpacing="2.500000" MarginWidth="1.600000" LineWidth="0.600000"'
    ' BoldWidth="2.000000" BondSpacing="18.000000" ChainAngle="120.000000"'
    ' LabelJustification="Auto" CaptionJustification="Left" FractionalWidths="yes"'
    ' InterpretChemically="yes" ShowAtomQuery="yes" ShowAtomStereo="no"'
    ' ShowBondQuery="yes" ShowBondStereo="no" HideImplicitHydrogens="no">\n'
    # The SAME eight entries as the ACS header, so COLORS is valid under either one. With
    # only white and black here, color="6" was silently out of range and rendered black --
    # a header that cannot express the standard palette is a trap, not a saving.
    '<colortable>\n'
    '<color r="1" g="1" b="1"/>\n'   # 2 white
    '<color r="0" g="0" b="0"/>\n'   # 3 black
    '<color r="1" g="0" b="0"/>\n'   # 4 red
    '<color r="1" g="1" b="0"/>\n'   # 5 yellow
    '<color r="0" g="1" b="0"/>\n'   # 6 green
    '<color r="0" g="1" b="1"/>\n'   # 7 cyan
    '<color r="0" g="0" b="1"/>\n'   # 8 blue
    '<color r="1" g="0" b="1"/>\n'   # 9 magenta
    '</colortable>\n'
    '<fonttable>\n'
    '<font id="21" charset="iso-8859-1" name="Arial"/>\n'
    '</fonttable>\n'
)


def minimal_header_prefix() -> str:
    """Synthetic best-effort header. Prefer header_prefix_from() with a real CDXML."""
    return _MIN_HEADER


def header_prefix_from(cdxml: str) -> str:
    """Borrow everything before <page> from a real ChemDraw CDXML string."""
    i = cdxml.index("<page")
    return cdxml[:i]


# ----------------------------------------------------------------------
# Fragment writer: atoms + bonds -> <fragment>
# ----------------------------------------------------------------------
# EVERY element, 1..118, because a hand-typed subset is a silent wrong answer.
#
# A hand-typed table is easy to get wrong: one that JUMPED from Ba(56) straight
# to W(74) would miss Re, Os, Ir and the whole lanthanide block plus Hf and Ta --
# and a reader that defaulted a miss to 6 would turn them into CARBON. An OsO4
# ester would draw as a carbon heterocycle with an unlabelled vertex where the
# metal should be, and the only complaint anywhere would be the InChI round trip
# saying "the drawing is not the molecule that was sent", which says nothing
# about why.
#
# The same silent-carbon default would hit a dummy atom too: its GetSymbol() is
# "*", which is not an element, so a Markush placeholder would come out as a
# METHYL (see smiles_to_frag, which routes dummies to placeholder() before they
# reach here).
#
# Symbols and numbers are GENERATED from rdkit.Chem.GetPeriodicTable(), not typed
# by hand -- that is the same library the symbols arrive from, so the two cannot
# disagree. cdxml_core does not import rdkit; the literal is written out instead.
_ELEMENTS = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8,
    "F": 9, "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15,
    "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22,
    "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29,
    "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
    "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43,
    "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
    "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57,
    "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64,
    "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71,
    "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78,
    "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85,
    "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90, "Pa": 91, "U": 92,
    "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99,
    "Fm": 100, "Md": 101, "No": 102, "Lr": 103, "Rf": 104, "Db": 105,
    "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110, "Rg": 111,
    "Cn": 112, "Nh": 113, "Fl": 114, "Mc": 115, "Lv": 116, "Ts": 117,
    "Og": 118,
}


def _atomic_number(symbol: str) -> int:
    """Atomic number for an element symbol. RAISES rather than guessing.

    Never `_ELEMENTS.get(sym, 6)`. A miss here is always a bug -- nicknames go
    through Placeholder and dummies through placeholder(), so whatever reaches
    this function is meant to be an element --
    and "silently a carbon" is the one outcome that cannot be noticed on paper.
    """
    try:
        return _ELEMENTS[symbol]
    except KeyError:
        raise ValueError(
            "element symbol %r is not one of the 118 this writer knows. A "
            "nickname belongs in a Placeholder, not in Atom.element; if this "
            "really is an element, the table above is generated from RDKit's "
            "periodic table and should already have it." % (symbol,)) from None


@dataclass
class Atom:
    x: float
    y: float
    element: str = "C"          # element symbol; carbons render label-free
    charge: int = 0             # formal charge; emitted as Charge= when nonzero
    # Enhanced stereochemistry -- the MDL ABS / AND / OR groups, which are what
    # says "drawn as one diastereomer, actually racemic with the relative
    # configuration defined". They ride ALONGSIDE the wedges, which stay exactly
    # as they are: the drawing is already faithful, and without this the embedded
    # data would contradict it. Nothing here writes AS -- ChemDraw computes the
    # CIP label itself, and this pipeline never writes one.
    stereo_type: str | None = None    # "Absolute" | "And" | "Or"
    stereo_group: int | None = None   # And/Or REQUIRE one; Absolute has no use for it
    # THE ATOM'S ELECTRON STATE. Both default to None and emit NOTHING when unset,
    # so every molecule without a radical is byte-for-byte unchanged.
    #
    #   radical      DTD line 239: Radical (None|Singlet|Doublet|Triplet).
    #   n_hydrogens  DTD line 237: NumHydrogens. Common in real ChemDraw files.
    #
    # WHY BOTH, when either one alone fixes the round trip (measured: a thiazol-
    # 2-ylidene written with neither reads back as CC1=C(C)N(C)CS1 -- the carbene
    # silently became a CH2 -- and NumHydrogens="0" OR Radical="Triplet" each restore
    # [C]): they say different things. NumHydrogens pins the COUNT, which is what
    # stops ChemDraw filling the free valence with hydrogen; Radical carries the
    # STATE, which is the thing a reader and a dot both need. A file with one and
    # not the other is a file that has to be re-derived to be understood.
    radical: str | None = None        # "Singlet" | "Doublet" | "Triplet"
    n_hydrogens: int | None = None    # explicit H count; pins the free valence
    # THE MASS NUMBER. Same additive standing as the two above: None emits nothing,
    # so every figure without an isotope is byte-for-byte unchanged.
    #
    #   isotope  DTD line 225, CDATA on <n>.
    #
    # LAYER 1 ONLY, and unusually so -- `Isotope="` appears in NONE of the 1173
    # CDXML files ChemDraw ships nor in this pipeline's output, so there is no layer-2
    # confirmation to be had and no point hunting for one. The DTD is enough here
    # because a mass number is single-valued: the case where the DTD misleads is
    # space-separated combination values (LineType, RectangleType), and this is not
    # one. What the DTD cannot say is whether ChemDraw DRAWS the superscript; that
    # takes opening a file, the same standing the radical dot has.
    isotope: int | None = None


@dataclass
class Bond:
    i: int                      # 0-based index into the atom list
    j: int
    order: int = 1              # 1/2/3 — always Kekulé, never 1.5 (that renders a ring circle)
    display: str | None = None  # e.g. "WedgeBegin", "WedgedHashBegin"
    # CDXML spends ONE `Order` attribute on two different ideas: how many bonds
    # there are (1/2/3) and what KIND of bond it is ("dative", and per the spec
    # also "ionic"/"hydrogen"/"threecenter", which ChemDraw drops on a round trip).
    # The data model keeps them apart -- multiplicity stays an int that arithmetic
    # and every reader can rely on, and the kind rides here and wins at write time.
    # Widening `order` to int|str instead would push the two-type check into every
    # place that reads it, to buy nothing.
    order_override: str | None = None    # "dative" -- VERIFIED by round trip


@dataclass
class Superatom:
    """A collapsed group shown as a text label (e.g. "NO2") that ChemDraw can
    Expand/Contract. Emitted as an <n NodeType="Fragment"> node whose nested
    <fragment> holds the real expanded atoms plus an ExternalConnectionPoint,
    followed by a <t> label. Put instances directly in the `fragment()` atom list
    alongside normal Atom objects."""
    x: float          # collapsed-node position (== the group's attachment atom)
    y: float
    label: str        # contracted text, e.g. "NO2"
    sub: list         # list[Atom]: the expanded group atoms
    sub_bonds: list   # list[Bond]: bonds among `sub` (0-based indices into `sub`)
    attach: int       # index in `sub` of the atom bonded to the parent (via the ECP)
    ecp_x: float      # parent attachment-atom position (ExternalConnectionPoint)
    ecp_y: float


@dataclass
class Placeholder:
    """A node whose label is NOT an element: R, R^1^, X, Ar, Nu, [N], [Pd].

    NOT a Superatom. A superatom is a real group that happens to be contracted --
    it carries its own expansion and one ExternalConnectionPoint, and ChemDraw can
    expand it back. A placeholder has no expansion (that is what makes it a
    placeholder) and no fixed number of bonds.

    Emitted the way ChemDraw itself writes one, read off a hand-drawn Markush file:

        <n NodeType="GenericNickname" GenericNickname="R" NumHydrogens="0"
           NeedsClean="yes" p="x y"><t ...><s ...>R</s></t></n>

    TWO bonds may end on the same placeholder and nothing extra is needed -- the X
    in that file is a chain linker with a bond on each side, and its markup is
    byte-identical to the terminal R's. ChemDraw trims every incident bond to the
    label box on its own, so the three topologies seen in real files (pendant,
    chain linker, ring vertex) are ONE case here, not three.

    Never reachable through Element=: ChemDraw reads a bare "Y" typed into a ring
    as yttrium, which is how the reference file ended up with an organoyttrium in
    it. The node type has to say generic, out loud.
    """
    x: float
    y: float
    label: str                       # drawn text; ^n^ superscript, ~n~ subscript
    nickname: str = ""               # GenericNickname=; "" -> the label, unmarked
    elements: tuple = ()             # non-empty -> an ElementList query node


@dataclass
class MultiAttachment:
    """A node whose bond attaches to a SET of atoms rather than to one.

        <n NodeType="MultiAttachment" Attachments="<ids>" p="x y"/>

    ChemDraw's own object, and the honest reading of a floating substituent --
    "bonded somewhere in this ring" is exactly what the attribute says. It is also
    what ChemDraw wrote for the metallocene in depictions/ferrocene, so the markup
    is not a guess.

    NOT the default for a floating bond, judged by eye: ChemDraw draws a
    visible asterisk at the node, and a printed figure cannot carry one. Kept
    because it is the correct object if that marker ever turns out to be
    suppressible.

    `attachments` are indices into the same atom list, resolved to node ids at
    emit time -- the ids do not exist until fragment() hands them out.
    """
    x: float
    y: float
    attachments: tuple = ()


def _multi_attachment_node(cp: "MultiAttachment", nid: int, ox: float, oy: float,
                           base_id: int) -> str:
    ids = " ".join(str(base_id + k) for k in cp.attachments)
    return (f'<n id="{nid}" p="{cp.x - ox:.2f} {cp.y - oy:.2f}"'
            f' NodeType="MultiAttachment" Attachments="{ids}"/>')


_LABEL_MARK = re.compile(r"(~[^~]+~|\^[^\^]+\^|\*[^*]+\*)")


def _label_plain(text: str) -> str:
    """The label with its markers removed -- what to MEASURE.

    Measuring the stripped string is exact here rather than approximate, because a
    ChemDraw atom-label superscript is set at FULL size (see _label_runs), so
    "R^1^" really does occupy as much width as "R1".
    """
    return _LABEL_MARK.sub(lambda m: m.group(0)[1:-1], text)


def parse_element_list(label: str) -> tuple:
    """Atomic numbers for a bracketed element label, or () if it is not one.

    "[N]" -> (7,);  "[N,O]" -> (7, 8);  "[Pd]" -> (46,);  "Ar" -> ().
    A bracketed element is a QUERY, not a nickname -- ChemDraw writes
    NodeType="ElementList" for it (verified: the [N] and [Pd] nodes in a
    hand-drawn Markush file carry ElementList="7" and ElementList="46"). Anything
    with a symbol outside _ELEMENTS falls back to (), i.e. to a generic nickname,
    so "[R]" or "[Nu]" stay placeholders instead of failing.
    """
    s = _label_plain(label).strip()
    if not (s.startswith("[") and s.endswith("]")):
        return ()
    syms = [p.strip() for p in s[1:-1].split(",")]
    if not syms or any(p not in _ELEMENTS for p in syms):
        return ()
    return tuple(_ELEMENTS[p] for p in syms)


def placeholder(x: float, y: float, label: str, nickname: str | None = None) -> Placeholder:
    """Build a Placeholder, picking the node type from the label itself."""
    return Placeholder(x, y, label, nickname or "", parse_element_list(label))


_SUPER_GLYPH = {
    "0": "\u2070", "1": "\u00b9", "2": "\u00b2", "3": "\u00b3", "4": "\u2074",
    "5": "\u2075", "6": "\u2076", "7": "\u2077", "8": "\u2078", "9": "\u2079",
    "+": "\u207a", "-": "\u207b", "\u2212": "\u207b",
    "(": "\u207d", ")": "\u207e", "n": "\u207f", "i": "\u2071",
}


def _fold_superscripts(text: str) -> str:
    """Turn ^...^ into single Unicode superscript glyphs where glyphs exist.

    ChemDraw puts a face="64" RUN on the side of the label AWAY from the bond. It
    is the right rule for "OMe" -> "MeO", and the wrong one for "R^1^", which comes
    out "1R" on every left-attached position. No CDXML attribute turns it off
    (checked: dropping GenericNickname, LabelJustification="Right", dropping the
    justification, NodeType="Unspecified" -- all five files render the same), and
    emitting the runs in the other order does not help either, because ChemDraw is
    placing the superscript run, not keeping the order it was given.

    A Unicode superscript is not a run, so that rule has nothing to move: "R\u00b9"
    renders the same on both sides, and matches how the journals print it.

    Only characters with a real superscript glyph fold. "^II^", "^a^", "^\u2022-^" have
    none, so they stay face="64" runs -- still reversed when left-attached, but a
    missing-glyph box would be worse.
    """
    def sub(m):
        body = m.group(0)[1:-1]
        if body and all(ch in _SUPER_GLYPH for ch in body):
            return "".join(_SUPER_GLYPH[ch] for ch in body)
        return m.group(0)
    return re.sub(r"\^[^\^]+\^", sub, text)


def _label_runs(text: str, size: int = ATOM_LABEL_SIZE) -> str:
    """Styled <s> runs for an ATOM label. Deliberately not styled_text().

    Two things an atom label does differently from a caption, both measured off
    the R^1 node in a hand-drawn Markush file:

      * a manual superscript is face="64" at the SAME size. styled_text() drops a
        caption's superscript two points; ChemDraw does not do that here. Most
        superscripts never reach that path -- see _fold_superscripts.

      * a label with NO markup keeps face="96", ChemDraw's auto-formula mode --
        which is what turns a plain "CF3" into CF-subscript-3 for free. It is also
        the trap: auto-formula would render "R1" as R-subscript-1. Write R^1^ when
        you mean R-superscript-1, and the marker turns auto-formula off, exactly
        as it does in ChemDraw's own file. A folded "R\u00b9" is safe under
        auto-formula: it has no ASCII digit left to subscript.
    """
    text = _fold_superscripts(text)
    parts = [p for p in _LABEL_MARK.split(text) if p != ""]
    if len(parts) == 1 and not _LABEL_MARK.fullmatch(parts[0]):
        return f'<s font="70" size="{size}" color="0" face="96">{_esc(text)}</s>'
    out = []
    for p in parts:
        face = 0
        if p.startswith("~") and p.endswith("~"):
            face, p = 32, p[1:-1]
        elif p.startswith("^") and p.endswith("^"):
            face, p = 64, p[1:-1]
        elif p.startswith("*") and p.endswith("*"):
            face, p = 2, p[1:-1]
        fa = f' face="{face}"' if face else ""
        out.append(f'<s font="70" size="{size}" color="0"{fa}>{_esc(p)}</s>')
    return "".join(out)


def _placeholder_node(ph: "Placeholder", nid: int, ox: float, oy: float) -> str:
    """Emit one placeholder <n>. No nested fragment, no connection point."""
    x, y = ph.x - ox, ph.y - oy
    # formula=True: this label is emitted with face="96" (see _label_run), so
    # ChemDraw shrinks its digits and charge signs when it draws it.
    w, _h = M.box([_label_plain(ph.label)], ATOM_LABEL_SIZE, formula=True)
    if ph.elements:
        head = (f'<n id="{nid}" p="{x:.2f} {y:.2f}" NodeType="ElementList"'
                f' NeedsClean="yes"'
                f' ElementList="{" ".join(str(z) for z in ph.elements)}">')
    else:
        # NumHydrogens="0" or ChemDraw adds implicit H to a generic node; NeedsClean
        # asks it to re-lay the label on open, which is also what absorbs the
        # half-pixel of side bearing the anchor below does not model.
        head = (f'<n id="{nid}" p="{x:.2f} {y:.2f}" NodeType="GenericNickname"'
                f' NumHydrogens="0" NeedsClean="yes"'
                f' GenericNickname="{_esc(ph.nickname or _label_plain(ph.label))}">')
    return (head
            + f'<t p="{x - w / 2.0:.2f} {y + LABEL_BASELINE_DROP:.2f}"'
              ' LabelJustification="Left">'
            + _label_runs(ph.label) + "</t></n>")


def regular_polygon(n: int, side: float = BOND_LENGTH, cx: float = 0.0, cy: float = 0.0,
                    start_deg: float = 0.0) -> list[tuple[float, float]]:
    """Vertices of a regular n-gon whose edge length == side (radius = side for n=6)."""
    r = side / (2 * math.sin(math.pi / n))
    pts = []
    for k in range(n):
        a = math.radians(start_deg + k * 360.0 / n)
        pts.append((cx + r * math.cos(a), cy + r * math.sin(a)))
    return pts


def _superatom_node(sa: "Superatom", nid: int, ox: float, oy: float, alloc) -> str:
    """Emit one <n NodeType="Fragment"> superatom: nested expansion + ECP + label.

    `alloc()` hands out fresh ids (kept above the main atom-id range so nothing
    collides). Coordinates are shifted by (ox, oy) to match the parent fragment.
    """
    fid = alloc()
    sub_ids = [alloc() for _ in sa.sub]
    ecp_id = alloc()
    out = [f'<n id="{nid}" p="{sa.x - ox:.2f} {sa.y - oy:.2f}" NodeType="Fragment">'
           f'<fragment id="{fid}">']
    for a, sid in zip(sa.sub, sub_ids):
        anum = _atomic_number(a.element)
        chg = f' Charge="{a.charge}"' if a.charge else ""
        nh = ' NumHydrogens="0"' if a.charge else ""   # pin 0 H on charged group atoms
        out.append(f'<n id="{sid}" p="{a.x - ox:.2f} {a.y - oy:.2f}" Element="{anum}"{chg}{nh}/>')
    out.append(f'<n id="{ecp_id}" p="{sa.ecp_x - ox:.2f} {sa.ecp_y - oy:.2f}"'
               ' NodeType="ExternalConnectionPoint"/>')
    for b in sa.sub_bonds:
        attrs = f'B="{sub_ids[b.i]}" E="{sub_ids[b.j]}"'
        if b.order != 1:
            attrs += f' Order="{b.order}"'
        out.append(f'<b id="{alloc()}" {attrs}/>')
    out.append(f'<b id="{alloc()}" B="{ecp_id}" E="{sub_ids[sa.attach]}"/>')
    out.append(f'</fragment><t p="{sa.x - ox:.2f} {sa.y - oy + 3.60:.2f}"'
               ' LabelJustification="Left">'
               f'<s font="70" size="{SUPERATOM_LABEL_SIZE}" color="0" face="96">{_esc(sa.label)}</s></t></n>')
    return "".join(out)


_ENHANCED_STEREO = ("Absolute", "And", "Or")


RADICALS = ("Singlet", "Doublet", "Triplet")

# How many drawn electrons each state shows, and which symbol shows them.
# Doublet is one unpaired electron; Triplet is two, drawn as two separate dots
# rather than a pair, because a pair is what a Singlet has. Measured in a ChemDraw
# test file: Radical alone draws NOTHING, `<graphic GraphicType="Symbol">` draws,
# and an Orbital draws a dot that is not attached to anything.
_RADICAL_SYMBOL = {"Doublet": ("Electron", 1),
                   "Triplet": ("Electron", 2),
                   "Singlet": ("LonePair", 1)}
# Where the dot sits: this far from the atom's centre, along its FREE VALENCE --
# the direction with no bond in it. Compared against the bond this pipeline draws:
# BOND_LENGTH is 14.4, so 0.45 of it clears a one- or two-character label and
# still reads as belonging to that atom rather than to the next one.
ELECTRON_OFFSET = 0.45 * BOND_LENGTH
ELECTRON_GAP = 0.30 * BOND_LENGTH      # between the two dots of a triplet

# A CHARGE ON A CARBON, drawn as a SYMBOL and not spelled into the atom label.
# Showing the CH is unnecessary and looks unprofessional -- an arenium's cationic
# ring carbon is a bare vertex with a + beside it, never the text "CH+". Left
# alone, ChemDraw does spell it: this emitter writes <n p="..." Charge="1"/> with
# no label at all and ChemDraw renders "CH+" from the Charge attribute on open. A
# Symbol graphic bound to the atom with <represent attribute="Charge"> is what
# stops that.
#
# EVERY NUMBER HERE IS MEASURED off a reference file drawn by hand with ChemDraw's
# own charge tool (Tools row 1, button 2) -- three ways in one file, and this is
# its first structure:
#   node          <n id="7" Charge="1" AS="N"/>        -- no label, no NumHydrogens
#   graphic       BoundingBox="166.36 59.66 158.86 59.66"  Z="28"
#                 GraphicType="Symbol" SymbolType="Plus"
#                 <represent attribute="Charge" object="7"/>
# so the box is written RIGHT EDGE FIRST, is 7.50 wide and ZERO high, its y is the
# atom's y, and its centre sits 0.527 bond from the atom. A sane
# top-left/bottom-right rectangle instead is NOT drawn by ChemDraw.
# The file's second structure is the same thing with SymbolType="CirclePlus" (the
# circled form), so switching is a one-word change.
# THREE FAMILIES, all bound the same way (<represent attribute="Charge">) and all
# drawn by hand in a ChemDraw reference file. Besides the plain and circled pairs
# there is RadicalCation / RadicalAnion -- that glyph carries the charge AND the
# odd electron in one mark, which is why a radical ion must NOT also get a
# separate dot.
#
# WIDTHS, off that file: the plain and circled pairs are 0.521 bond, the radical
# ions HALF that. The box is written (anchor, anchor - width): the glyph sits at
# the FIRST coordinate, the same trap the lone pairs did not have.
#
# OFFSETS, from the same file, and they differ by family and by whether the atom has a
# drawn label to clear:
#
#                       labelled atom     bare carbon vertex
#     plain                  0.51               0.27
#     circled                0.66               0.40      <- the circle is bigger
#     radical ion            0.56               0.36
#
# Cross-check against a reference set of published mechanism figures: 977 circled
# charges measured there have a circle of radius 0.290 bond, i.e. 0.58 bond across,
# against ChemDraw's 0.52 -- inside that set's own 10th-90th spread of 0.41-0.69.
CHARGE_SYMBOL_WIDTH = 7.50             # 0.521 bond; Plus, Minus, CirclePlus, CircleMinus
RADICAL_ION_WIDTH = 0.260 * BOND_LENGTH
CHARGE_SYMBOL_OFFSET = {
    "plain":   (0.51 * BOND_LENGTH, 0.27 * BOND_LENGTH),
    "circled": (0.66 * BOND_LENGTH, 0.40 * BOND_LENGTH),
    "radical": (0.56 * BOND_LENGTH, 0.36 * BOND_LENGTH),
}
_CHARGE_SYMBOL = {1: "Plus", -1: "Minus"}
_CHARGE_SYMBOL_CIRCLED = {1: "CirclePlus", -1: "CircleMinus"}
_CHARGE_SYMBOL_RADICAL = {1: "RadicalCation", -1: "RadicalAnion"}
# CIRCLED IS THE DEFAULT. Published mechanism figures circle their charges (977
# in the reference set), and the circled form is the house choice for both the
# delocalized arene cation and an atom-bound charge. The plain form draws a mark
# that published mechanism figures do not use.
CHARGE_CIRCLED_DEFAULT = True


def _free_valence_dir(k, atoms, bonds):
    """A unit vector pointing where the atom has no bond, for atom index `k`.

    The negated sum of the directions to its neighbours: for a terminal atom that
    is straight out along the bond, for a two-bond atom it bisects the outside of
    the angle. Falls back to straight up for an isolated atom, which has no
    direction to prefer.
    """
    ax, ay = atoms[k].x, atoms[k].y
    sx = sy = 0.0
    for b in bonds:
        o = b.j if b.i == k else (b.i if b.j == k else None)
        if o is None or o >= len(atoms):
            continue
        dx, dy = atoms[o].x - ax, atoms[o].y - ay
        n = math.hypot(dx, dy)
        if n:
            sx += dx / n
            sy += dy / n
    n = math.hypot(sx, sy)
    if n < 1e-9:
        return 0.0, -1.0
    return -sx / n, -sy / n


def _electron_symbols(atoms, bonds, base_id, ox, oy, alloc) -> list:
    """`<graphic GraphicType="Symbol">` for every atom carrying a radical state.

    THREE THINGS ARE EMITTED TOGETHER and each does a different job:
      Radical= on the node   the chemistry, and what makes the molecule survive
                             the round trip. Draws nothing by itself.
      <graphic Symbol>       the dot. This is the only one that draws.
      <represent>            binds the dot to the atom, so ChemDraw knows the two
                             belong together rather than a dot happening to sit
                             there. NOT required to draw (a test file drew
                             without it) -- emitted because a dot that is not
                             bound to its atom is the same silent-drift failure
                             as an arrow stranded at the origin.

    INSIDE the <fragment>, which is what makes this free in all six engines:
    `place()` renumbers ids (including `object`) and `_translate_body` moves a
    BoundingBox, so the dot travels with its structure without any engine knowing
    it exists. UNVERIFIED: whether ChemDraw KEEPS it nested when it re-saves.
    """
    out = []
    for k, a in enumerate(atoms):
        r = getattr(a, "radical", None)
        if not r:
            continue
        if charge_is_radical_ion(a):
            # the RadicalCation / RadicalAnion glyph already carries this electron
            continue
        sym, count = _RADICAL_SYMBOL[r]
        ux, uy = _free_valence_dir(k, atoms, bonds)
        cx, cy = a.x - ox + ux * ELECTRON_OFFSET, a.y - oy + uy * ELECTRON_OFFSET
        # two dots straddle the free-valence direction, perpendicular to it
        px, py = -uy, ux
        for e in range(count):
            t = 0.0 if count == 1 else (e - (count - 1) / 2.0) * ELECTRON_GAP
            x, y = cx + px * t, cy + py * t
            out.append(
                '<graphic id="%d" BoundingBox="%.2f %.2f %.2f %.2f"'
                ' GraphicType="Symbol" SymbolType="%s">'
                '<represent attribute="Radical" object="%d"/></graphic>'
                % (alloc(), x, y, x, y, sym, base_id + k))
    return out


# HOW FAR THE SYMBOL MAY BE PULLED IN when the roomy offset puts it nearer some
# other atom. The hand-drawn value is the first rung; the rest are fallbacks.
CHARGE_SYMBOL_OFFSET_LADDER = (1.00, 0.83, 0.68, 0.55)


def _gap_dirs(k, atoms, bonds):
    """Every gap bisector at atom k, widest first."""
    ax, ay = atoms[k].x, atoms[k].y
    angs = []
    for b in bonds:
        o = b.j if b.i == k else (b.i if b.j == k else None)
        if o is None or o >= len(atoms):
            continue
        dx, dy = atoms[o].x - ax, atoms[o].y - ay
        if math.hypot(dx, dy):
            angs.append(math.atan2(dy, dx))
    if not angs:
        return [(0.0, -1.0)]
    angs.sort()
    if len(angs) == 1:
        a = angs[0] + math.pi
        return [(math.cos(a), math.sin(a))]
    gaps = []
    for i in range(len(angs)):
        w = (angs[(i + 1) % len(angs)] - angs[i]) % (2 * math.pi)
        gaps.append((w, angs[i] + w / 2.0))
    gaps.sort(key=lambda g: -g[0])
    return [(math.cos(a), math.sin(a)) for _w, a in gaps]


def _place_symbol(k, atoms, bonds, off):
    """(x, y) for atom k's symbol: the first gap and offset that leaves it NEARER
    ITS OWN ATOM than any other.

    THE HARD CONSTRAINT, and it is not cosmetic. ChemDraw re-associates a charge
    symbol with the atom it sits by when it re-reads the file, whatever the
    <represent> says. Measured over three round trips: a symbol 0.34 bond from a
    neighbour moved the charge onto that neighbour, and on an ammonium next to a
    bromide it swapped the two charges outright. The nearest OTHER atom includes
    the other molecules of the same cell, because a cell's whole SMILES is one
    fragment and the packer can set two ions close together.
    """
    ax, ay = atoms[k].x, atoms[k].y
    others = [(a.x, a.y) for j, a in enumerate(atoms) if j != k]
    best = None
    for ux, uy in _gap_dirs(k, atoms, bonds):
        for scale in CHARGE_SYMBOL_OFFSET_LADDER:
            d = off * scale
            px, py = ax + ux * d, ay + uy * d
            near = min((math.hypot(px - qx, py - qy) for qx, qy in others),
                       default=float("inf"))
            if best is None or near - d > best[0]:
                best = (near - d, px, py)
            if near > d * 1.25:              # comfortably its own
                return px, py
    return (best[1], best[2]) if best else (ax, ay - off)


def _widest_gap_dir(k, atoms, bonds):
    """Into the widest ANGULAR GAP between atom k's bonds, not along the summed
    free valence.

    THE FREE VALENCE IS DEGENERATE ON A FOUR-COORDINATE ATOM. _free_valence_dir
    negates the sum of the bond unit vectors; on a phosphonium drawn in 2D the
    four roughly cancel, so the sum is near zero and its direction is whatever
    the rounding left -- measured, the charge symbol can land 0.66 bond from its
    own phosphorus and 0.34 bond from a NEIGHBOUR. ChemDraw then re-associates it
    on load: every phosphonium in a six-row test file came back with the
    phosphorus neutral and the +1 sitting on whichever atom the symbol was
    nearest.

    The gap rule has no degenerate case: an atom with bonds always has a widest
    gap, and for one bond its centre IS the free valence, so nothing with one or
    two bonds moves a byte.
    """
    ax, ay = atoms[k].x, atoms[k].y
    angs = []
    for b in bonds:
        o = b.j if b.i == k else (b.i if b.j == k else None)
        if o is None or o >= len(atoms):
            continue
        dx, dy = atoms[o].x - ax, atoms[o].y - ay
        if math.hypot(dx, dy):
            angs.append(math.atan2(dy, dx))
    if not angs:
        return 0.0, -1.0
    angs.sort()
    if len(angs) == 1:
        a = angs[0] + math.pi
        return math.cos(a), math.sin(a)
    best_w, best_a = -1.0, angs[0] + math.pi
    for i in range(len(angs)):
        w = (angs[(i + 1) % len(angs)] - angs[i]) % (2 * math.pi)
        if w > best_w:
            best_w, best_a = w, angs[i] + w / 2.0
    return math.cos(best_a), math.sin(best_a)


def _charge_symbols(atoms, bonds, base_id, ox, oy, alloc, circled=False) -> list:
    """`<graphic GraphicType="Symbol">` for a charge that would otherwise be text.

    The sibling of _electron_symbols, and only for a CARBON. A heteroatom already
    draws its own symbol, and ChemDraw sets the charge beside it correctly there
    (O with a minus, N with a plus) -- emitting a second one would double it. A
    carbon has no label to hang a charge off, so ChemDraw invents one ("CH+"),
    which is what this replaces.

    Direction is _free_valence_dir, the same rule the radical dots use and the
    same rule ChemDraw itself uses to decide where implicit hydrogens go: out
    into the space the atom has no bond in. The hand-drawn reference places it
    to the left because that is the open side of its ring; the rule generalises
    that rather than copying the coordinate.
    """
    out = []
    for k, a in enumerate(atoms):
        chg = int(getattr(a, "charge", 0) or 0)
        if not chg:
            continue
        sign = 1 if chg > 0 else -1
        if charge_is_radical_ion(a):
            fam, table, width = "radical", _CHARGE_SYMBOL_RADICAL, RADICAL_ION_WIDTH
        elif circled:
            fam, table, width = "circled", _CHARGE_SYMBOL_CIRCLED, CHARGE_SYMBOL_WIDTH
        else:
            fam, table, width = "plain", _CHARGE_SYMBOL, CHARGE_SYMBOL_WIDTH
        sym = table.get(sign)
        if not sym:
            continue
        labelled, bare = CHARGE_SYMBOL_OFFSET[fam]
        off = bare if getattr(a, "element", None) == "C" else labelled
        px, py = _place_symbol(k, atoms, bonds, off)
        cx, cy = px - ox, py - oy
        out.append(
            '<graphic id="%d" BoundingBox="%.2f %.2f %.2f %.2f"'
            ' GraphicType="Symbol" SymbolType="%s">'
            '<represent attribute="Charge" object="%d"/></graphic>'
            % (alloc(), cx, cy, cx - width, cy, sym, base_id + k))
    return out


def charge_is_radical_ion(a) -> bool:
    """A charge AND one unpaired electron, i.e. one RadicalCation/Anion glyph.

    THE ONE PLACE that decides it, because two emitters have to agree: the charge
    symbol draws the combined glyph and the dot emitter must then draw nothing.
    Hand-drawn radical-ion structures carry no separate Electron graphic at all.
    """
    chg = int(getattr(a, "charge", 0) or 0)
    return bool(chg) and getattr(a, "radical", None) == "Doublet"
def _isotope_attr(a) -> str:
    """` Isotope="N"` for a labelled atom, else "".

    Kept apart from _electron_attrs rather than folded into it: a mass number is
    not an electron state, and the two have different evidence behind them (that
    one is DTD + ChemDraw's own files, this one is DTD only -- see Atom.isotope).

    RDKit spells "no label" as 0, not None, so the falsy test is the whole check.
    """
    n = getattr(a, "isotope", None)
    if not n:
        return ""
    return ' Isotope="%d"' % int(n)


def _electron_attrs(a) -> str:
    """` Radical="..." NumHydrogens="N"` for an atom that carries either, else "".

    Order is DTD order (NumHydrogens at line 237, Radical at 239) so the two never
    swap between runs. Returns the empty string for an ordinary atom, which is the
    whole of "additive": nothing without a radical moves a byte.
    """
    out = ""
    nh = getattr(a, "n_hydrogens", None)
    if nh is not None:
        out += f' NumHydrogens="{int(nh)}"'
    r = getattr(a, "radical", None)
    if r is not None:
        if r not in RADICALS:
            # The DTD enumerates four and "None" is the default, i.e. the absence
            # of the attribute -- writing it out would be a fifth spelling of "no".
            raise ValueError(
                "Atom.radical is one of %s (or None for an ordinary atom), got %r"
                % (" / ".join(map(repr, RADICALS)), r))
        out += f' Radical="{r}"'
    return out


def _enhanced_stereo_attrs(a: "Atom") -> str:
    """` EnhancedStereoType=... EnhancedStereoGroupNum=...` for one atom, or "".

    THE evidence levels differ between the three values, and they differ enough to
    write down, because the fix when one of them turns out wrong is different:

      "And"       CD-verified. ChemDraw's own sample document writes it --
                  Samples/b-Cypermethrin (stereochemistry).cdxml, three centres in
                  groups 1 and 2, rendered as "&1 (R)" / "&2 (S)". That is
                  ChemDraw's own bytes, which is stronger than a round trip.
      "Or"        FORMAT DOCUMENTATION ONLY. No file anywhere carries it -- not
                  ChemDraw's 652 samples, not this pipeline's output.
      "Absolute"  FORMAT DOCUMENTATION ONLY, same standing as "Or".

    A specification-only value is a CANDIDATE, not a fact, and there is a
    measured reason: `dative`, `ionic`, `hydrogen` and `threecenter` came off ONE
    page of that same specification, and a round trip through ChemDraw kept
    exactly one of the four. So if an Or or Absolute group comes out wrong, the
    procedure is known and cheap: write a test file, open it, Cmd+S, read it
    back, and see whether the attribute is still there.

    The PAIRING is checked rather than assumed. The specification is explicit that
    Or and And require a group number, so emitting one without it would produce a
    file that says something the format does not define -- the silent-wrong-output
    shape this module exists to prevent. Raising names the atom instead.
    """
    if a.stereo_type is None:
        return ""
    if a.stereo_type not in _ENHANCED_STEREO:
        raise ValueError(
            "stereo_type=%r; EnhancedStereoType is one of %s"
            % (a.stereo_type, " / ".join(_ENHANCED_STEREO)))
    out = ' EnhancedStereoType="%s"' % a.stereo_type
    if a.stereo_type in ("And", "Or"):
        if a.stereo_group is None:
            raise ValueError(
                'EnhancedStereoType="%s" at (%.2f, %.2f) needs a stereo_group: the '
                "CDX spec requires a group number for And and Or, and the number is "
                "what says which centres invert together." % (a.stereo_type, a.x, a.y))
        out += ' EnhancedStereoGroupNum="%d"' % a.stereo_group
    return out


def fragment(atoms, bonds, base_id: int = 1, normalize: bool = True):
    """Build a <fragment> CDXML string. Returns (fragment_str, bbox=[x0,y0,x1,y1]).

    `atoms` may contain Atom, Superatom and Placeholder entries. Atom k gets node
    id (base_id + k); Superatom internals get ids above the main range. If normalize,
    coords are shifted so the bounding box's top-left sits at (0,0) -> [0,0,w,h].
    """
    # bbox-contributing points: a plain atom is its own point; a Superatom also
    # contributes its contracted label's text box (real width/height from
    # text_metrics), placed horizontally on the side away from the parent -- so
    # layout leaves room for a long label like "OTBDPS" instead of packing the next
    # structure against a bbox that only saw the anchor point.
    # Who is bonded to whom, and how heavily -- needed both to count implicit
    # hydrogens and to know which way a label hangs.
    _nbrs = {}
    _order = {}
    for b in bonds:
        _nbrs.setdefault(b.i, []).append(b.j)
        _nbrs.setdefault(b.j, []).append(b.i)
        _order[b.i] = _order.get(b.i, 0) + b.order
        _order[b.j] = _order.get(b.j, 0) + b.order

    def _extent(a, k=None):
        if isinstance(a, MultiAttachment):
            return [(a.x, a.y)]           # no ink; only the bond reaching it is drawn

        if isinstance(a, Placeholder):
            # CENTRED on the point, unlike a superatom. A superatom's node IS its
            # attachment atom, so "OTBDPS" hangs off to one side; a placeholder is a
            # word sitting on a point with no atom inside it, and ChemDraw centres
            # it -- every labelled node in a hand-drawn Markush file has its box
            # centred on the node's y exactly and on its x within a third of a unit,
            # including the ones with a bond coming in from the right.
            lw, lh = M.box([_label_plain(a.label)], ATOM_LABEL_SIZE, formula=True)
            m = ATOM_LABEL_MARGIN
            return [(a.x, a.y),
                    (a.x - lw / 2.0 - m, a.y - lh / 2.0 - m),
                    (a.x + lw / 2.0 + m, a.y + lh / 2.0 + m)]

        if isinstance(a, Superatom):
            lw, lh = M.box([a.label], SUPERATOM_LABEL_SIZE, formula=True)
            if a.x >= a.ecp_x:              # bond from the left -> label extends right
                x0, x1 = a.x, a.x + lw
            else:                           # bond from the right -> label extends left
                x0, x1 = a.x - lw, a.x
            return [(a.x, a.y), (x0, a.y - lh / 2.0), (x1, a.y + lh / 2.0)]

        # A carbon renders label-free (ShowTerminalCarbonLabels=no), so it is
        # just its point. Everything else is drawn as text and that text is what
        # a neighbouring cell gets placed against.
        if a.element == "C" or k is None:
            return [(a.x, a.y)]

        nh = _implicit_h(a.element, a.charge, _order.get(k, 0))
        along, half_sym, half_h = _label_overhang(a.element, nh)

        nb = _nbrs.get(k, [])
        if len(nb) == 1:
            # hydrogens hang off the side away from the single bond; a near
            # vertical bond leaves the choice to ChemDraw, which puts them right
            other = atoms[nb[0]]
            right = (a.x - other.x) >= -1e-9
        else:
            right = True                    # centred symbol, H to the right
        x0 = a.x - (half_sym if right else along)
        x1 = a.x + (along if right else half_sym)
        return [(a.x, a.y), (x0, a.y - half_h), (x1, a.y + half_h)]
    # A FRAGMENT WITH NO ATOMS. `min(exs)` on it raises "min() iterable argument is
    # empty" and that is what an empty SMILES would show the caller -- Python's
    # internals, about a molecule they asked for. RDKit parses "" without
    # complaining (it is a valid, empty molecule), so nothing upstream stops it.
    if not atoms:
        raise ValueError(
            "no atoms to draw: this structure is empty. An empty SMILES string "
            "parses without error but has nothing in it -- pass a real one, "
            "e.g. 'CCO'.")
    ext = [p for k, a in enumerate(atoms) for p in _extent(a, k)]
    exs = [p[0] for p in ext]
    eys = [p[1] for p in ext]
    ox, oy = (min(exs), min(eys)) if normalize else (0.0, 0.0)

    next_id = [base_id + len(atoms)]          # ids for superatom internals live here
    def alloc() -> int:
        v = next_id[0]
        next_id[0] += 1
        return v

    node_lines = []
    for k, a in enumerate(atoms):
        nid = base_id + k
        if isinstance(a, MultiAttachment):
            node_lines.append(_multi_attachment_node(a, nid, ox, oy, base_id))
            continue
        if isinstance(a, Placeholder):
            node_lines.append(_placeholder_node(a, nid, ox, oy))
            continue
        if isinstance(a, Superatom):
            node_lines.append(_superatom_node(a, nid, ox, oy, alloc))
            continue
        x, y = a.x - ox, a.y - oy
        chg = f' Charge="{a.charge}"' if a.charge else ""
        est = _enhanced_stereo_attrs(a)
        # APPENDED, never inserted: an atom with neither field emits the empty
        # string and the line is unchanged.
        el = _electron_attrs(a)
        # Isotope is DTD:225 and NumHydrogens/Radical are 237/239, so it goes
        # BEFORE them -- same rule _electron_attrs follows internally, so the file
        # reads in DTD order however many of the three an atom carries.
        iso = _isotope_attr(a)
        if a.element == "C":
            node_lines.append(f'<n id="{nid}" p="{x:.2f} {y:.2f}"{chg}{est}{iso}{el}/>')
        else:
            anum = _atomic_number(a.element)
            # heteroatom: Element attr is enough — ChemDraw draws the symbol and
            # infers implicit H via InterpretChemically.
            node_lines.append(
                f'<n id="{nid}" p="{x:.2f} {y:.2f}" Element="{anum}"{chg}{est}{iso}{el}/>')

    bond_lines = []
    for b in bonds:
        bi, bj = base_id + b.i, base_id + b.j
        attrs = f'B="{bi}" E="{bj}"'
        if b.order_override:
            # Wins over multiplicity: a dative bond IS order 1 in the connection
            # table, and saying Order="1" as well is not possible -- one attribute.
            attrs += f' Order="{b.order_override}"'
        elif b.order != 1:
            attrs += f' Order="{b.order}"'
        if b.display:
            attrs += f' Display="{b.display}"'
        # id from the SAME allocator the superatom internals use, so a bond id can
        # collide with neither an atom id (those are base_id .. base_id+len(atoms)-1)
        # nor a nested expansion's ids (allocated above them, in the node loop, which
        # has already finished). renumber() shifts "id" along with B/E, so a placed
        # fragment stays internally consistent.
        #
        # Bonds need ids because BasisObjects names BONDS -- it is how a native
        # ChemDraw highlight says which ring it covers (see shaded_band).
        bond_lines.append(f'<b id="{alloc()}" {attrs}/>')

    bbox = [min(exs) - ox, min(eys) - oy, max(exs) - ox, max(eys) - oy]
    dot_lines = (_electron_symbols(atoms, bonds, base_id, ox, oy, alloc)
                 + _charge_symbols(atoms, bonds, base_id, ox, oy, alloc,
                                   circled=CHARGE_CIRCLED_DEFAULT))
    frag = (
        f'<fragment BoundingBox="{bbox[0]:.2f} {bbox[1]:.2f} {bbox[2]:.2f} {bbox[3]:.2f}">\n'
        + "\n".join(node_lines) + "\n" + "\n".join(bond_lines)
        + ("\n" + "\n".join(dot_lines) if dot_lines else "")
        + "\n</fragment>"
    )
    return frag, bbox


# ----------------------------------------------------------------------
# B.* assembly helpers (ported from the web library, PART 1)
# ----------------------------------------------------------------------
_P_RE = re.compile(r'\bp="(-?[\d.]+) (-?[\d.]+)"')
_BB_RE = re.compile(r'\bBoundingBox="(-?[\d.]+) (-?[\d.]+) (-?[\d.]+) (-?[\d.]+)"')
# A <curve> carries its geometry in a flat x y x y ... list and has no p="" at all,
# so a fragment holding one (the arc inside a metallocene ring) leaves its arcs
# behind when moved unless this is translated too.
_CURVE_RE = re.compile(r'\bCurvePoints="([-\d.\s]+)"')
# Attributes that hold a LIST of ids rather than one. BondCircularOrdering uses 0
# for "no bond in this position", and 0 has to stay 0 -- no object has id 0.
#
# BasisObjects is here because a ColoredMolecularArea names the BONDS it covers.
# Leave it out and place() renumbers the bonds while the highlight keeps pointing
# at their old numbers -- the fill silently lands on the wrong bonds,
# or on none. Not reachable through the table path (shaded_band emits at final
# coordinates and is never placed), but a fill lives INSIDE a structure's
# fragment, and structures are placed constantly.
_ID_LIST_ATTRS = ("Attachments", "CrossingBonds", "BondCircularOrdering",
                  "BasisObjects")


def bbox_of(s: str) -> list[float] | None:
    m = _BB_RE.search(s)
    return [float(g) for g in m.groups()] if m else None


def translate(frag: str, dx: float, dy: float) -> str:
    """Move every coordinate in a fragment/object string by (dx, dy).

    Covers p="", BoundingBox="" and CurvePoints="". NOT the 3D endpoints an <arrow>
    carries (Head3D/Tail3D/Center3D/Major-MinorAxisEnd3D) -- arrows are emitted at
    their final position, never moved, so translating one would be a new feature
    rather than a fix.
    """
    s = _P_RE.sub(lambda m: f'p="{float(m.group(1)) + dx:.2f} {float(m.group(2)) + dy:.2f}"', frag)
    s = _BB_RE.sub(
        lambda m: 'BoundingBox="%.2f %.2f %.2f %.2f"' % (
            float(m.group(1)) + dx, float(m.group(2)) + dy,
            float(m.group(3)) + dx, float(m.group(4)) + dy),
        s,
    )

    def _curve(m):
        v = [float(t) for t in m.group(1).split()]
        v[0::2] = [x + dx for x in v[0::2]]
        v[1::2] = [y + dy for y in v[1::2]]
        return 'CurvePoints="%s"' % " ".join("%.2f" % t for t in v)

    return _CURVE_RE.sub(_curve, s)


def renumber(frag: str, base: int) -> str:
    """Shift every id (and every reference to one) by `base`."""
    # `object` is <represent>'s back-reference to the node a Symbol shows. Only
    # Symbol graphics emit `object="`, so every other string comes out unchanged.
    for attr in ("id", "B", "E", "BeginAttach", "EndAttach", "object"):
        frag = re.sub(
            r'\b' + attr + r'="(\d+)"',
            lambda m: f'{attr}="{base + int(m.group(1))}"',
            frag,
        )
    for attr in _ID_LIST_ATTRS:
        frag = re.sub(
            r'\b' + attr + r'="([\d\s]+)"',
            lambda m: '%s="%s"' % (attr, " ".join(
                str(base + int(t)) if int(t) else "0" for t in m.group(1).split())),
            frag,
        )
    return frag


def place(frag: str, base: int, x: float, y: float) -> str:
    """Renumber ids by `base` (avoid collisions), then translate so the object's
    top-left bbox corner sits at (x, y).

    MOVES THE ABSOLUTE 3D ENDPOINTS TOO, which is why this says "object" and not
    "fragment". A `<fragment>` without arrows has none of them, so it comes out
    byte-identical to a plain translate -- verified across 24 figures covering all
    six engines plus stacks, slots, graphics, decorations, annotations,
    connectors, legends, sections and tables, and across the goldens. What it
    buys is placing a `<group>` that contains an ARROW.

    Plain translate() is not enough for that: it leaves Head3D/Tail3D where they
    were, so a placed group comes out with its structures moved and its arrow
    stranded at the origin. Measured, and measured to be worse than it looks: such
    a group survives a ChemDraw round trip still torn, and ChemDraw then reads the
    stranded arrow as a reaction step whose product is the whole group -- the
    drawing lies and the embedded reaction data lies with it.

    The same test answered the other half: ChemDraw KEEPS a `<group>` with an arrow
    inside it, selects it as one object and moves the arrow with it. So a composite
    panel is a real thing to build, and this line is its entry ticket.

    NOT the same as place_with(), which exists for a decoration that is a SIBLING of
    its structure and must move by the structure's shift rather than its own.
    """
    bb = bbox_of(frag)
    s = renumber(frag, base)
    if bb:
        s = _translate_body(s, x - bb[0], y - bb[1])
    return s


def place_with(frag: str, deco: str | None, base: int, x: float, y: float):
    """place() the fragment AND carry its decorations by the same shift.

    Returns (fragment_xml, decoration_xml) -- decoration_xml is "" when there is
    nothing to carry, and then the fragment string is byte-identical to place()'s.
    That is what makes this safe to substitute at every call site: a figure with
    no decorations cannot come out different.

    THE SHIFT IS COMPUTED ONCE AND SPENT TWICE, and that is the whole point.
    Letting the decorations go through place() of their own would shift each piece
    by ITS OWN bbox -- which slides a halo off the group it belongs to, because a
    glow's box is 24 units bigger than the atoms' on every side. The failure is
    nearly invisible on paper: a halo 24 units off centre reads as "marked a bit
    carelessly", not as "wrong". Pinned by
    test_place_with_spends_one_shift_on_both.

    translate_objects for the decorations, _translate_body for the fragment, and
    the difference is not an inconsistency: a disc keeps its geometry in absolute
    Center3D / MajorAxisEnd3D, so moving p= and BoundingBox alone tears it in
    half. Written down in three places (translate's docstring,
    translate_objects', legend.place_legend') because the mistake is easy to make.

    "A fragment has no 3D endpoints, so translate is the whole of its move" is
    FALSE: an electron-pushing annotation is a native <arrow>, so an annotated
    fragment holds Center3D / Tail3D / Head3D / Major- / MinorAxisEnd3D, and plain
    translate leaves them where they were -- the arrow stranded at the origin
    while its structure moves, which is the exact failure place()'s own docstring
    warns about. Measured on a real mechanism figure: the pushing arrows came out
    20 to 211 units from the atoms they belong to. `place()` uses _translate_body
    for this reason; this is the same fix on the other door.

    Byte-safe: _translate_body is translate plus a regex over the five 3D
    attributes, and a fragment carrying none of them comes out unchanged.

    `deco` IS renumbered, by the same `base` as the fragment. Exempting it on the
    grounds that decoration ids come from a high pool (7000+ / 7500+) "which no
    fragment base reaches" would be false: every engine hands structures ids with
    `base += 1000`, so the SEVENTH cell's base is exactly 7000 and its atoms are
    7001, 7002, ... Measured on ten structures with six dots on the first:
    `<graphic id="7001">` and `<n id="7001">` in one document. Raising the pool
    only raises the threshold; the number of structures in a figure has no
    ceiling.

    So decorate.apply numbers from the fragment's own range (next_free_id),
    and they are renumbered together here. A decoration lands inside ITS OWN
    cell's thousand, and stays unique for exactly as long as the atom ids do --
    which is next_free_id's whole argument, applied to the objects beside the
    fragment as well as the ones inside it.
    """
    bb = bbox_of(frag)
    dx, dy = (x - bb[0], y - bb[1]) if bb else (0.0, 0.0)
    frag_xml = renumber(frag, base)
    if bb:
        frag_xml = _translate_body(frag_xml, dx, dy)
    if not deco:
        return frag_xml, ""
    return frag_xml, _translate_body(renumber(deco, base), dx, dy)


def arrow(aid: int, cx: float, y: float, length: float) -> str:
    """Native reaction arrow, head on the right. cx = center, y = baseline."""
    x0, x1 = cx - length / 2, cx + length / 2
    return (
        f'<arrow id="{aid}" BoundingBox="{x0:.2f} {y - 2:.2f} {x1:.2f} {y + 2:.2f}"'
        ' FillType="None" ArrowheadHead="Full" ArrowheadType="Solid" HeadSize="1000"'
        ' ArrowheadCenterSize="875" ArrowheadWidth="250"'
        f' Head3D="{x1:.2f} {y:.2f} 0" Tail3D="{x0:.2f} {y:.2f} 0" Center3D="{cx:.2f} {y:.2f} 0"/>'
    )


_SUB_SPLIT = re.compile(r"(~[^~]+~)")


def _esc(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def text(tid: int, cx: float, top_y: float, lines, size: int = 9, face: int | None = None,
         align: str = "center") -> str:
    """Centered (multi-line) text block. ~x~ marks a subscript run (face 32).
    align: 'center' (default) | 'left' | 'right'."""
    if isinstance(lines, str):
        lines = [lines]

    def runs(line: str) -> str:
        out = []
        for p in (p for p in _SUB_SPLIT.split(line) if p != ""):
            if p.startswith("~") and p.endswith("~"):
                out.append(f'<s font="21" size="{size}" face="32">{_esc(p[1:-1])}</s>')
            else:
                fa = f' face="{face}"' if face else ""
                out.append(f'<s font="21" size="{size}"{fa}>{_esc(p)}</s>')
        return "".join(out)

    inner = ""
    for i, l in enumerate(lines):
        if i > 0:
            inner += f'<s font="21" size="{size}">\n</s>'
        inner += runs(l)
    just = {"left": "Left", "right": "Right", "center": "Center"}.get(align, "Center")
    return (
        f'<t id="{tid}" p="{cx:.2f} {top_y:.2f}" Justification="{just}"'
        f' LineHeight="auto" InterpretChemically="no">{inner}</t>'
    )


def number(tid: int, cx: float, top_y: float, label) -> str:
    """Bold, larger compound-number label."""
    return text(tid, cx, top_y, [str(label)], size=11, face=1)


def page(header_prefix: str, body: str, bbox) -> str:
    """Wrap a body (fragments + arrows + texts) into a full CDXML document."""
    x0, y0, x1, y1 = bbox
    return (
        header_prefix
        + f'<page id="1" BoundingBox="{x0} {y0} {x1} {y1}" Width="{x1 - x0}" Height="{y1 - y0}"'
        ' HeightPages="1" WidthPages="1" DrawingSpace="poster">\n'
        + body
        + "\n</page>\n</CDXML>\n"
    )


# ----------------------------------------------------------------------
# Document finalize: positive-quadrant shift + tight page bbox
# ----------------------------------------------------------------------
# ChemDraw's file-open / Quick Look renderers place the page origin at (0,0) and
# clip negative coordinates (in CDXML, -y is up). The web core tolerated negatives
# because it auto-fit after load; a static file render does not. So the finalize
# step shifts ALL content (atom p=, BoundingBox, and arrow Head3D/Tail3D/Center3D)
# into the positive quadrant and wraps the page tightly around the real content.
# Every ABSOLUTE 3D point an object can carry. MajorAxisEnd3D / MinorAxisEnd3D are
# absolute too (centre + radius vector), so an arc whose centre is shifted without them
# ends up with garbage radii -- and finalize() shifts every document into the positive
# quadrant, so omitting them silently corrupts every arc that goes through the layout
# engine. Verified: a 40x40 arc shifted +200,+300 came back with radii (-160, -260).
_3D_RE = re.compile(r'\b(Head3D|Tail3D|Center3D|MajorAxisEnd3D|MinorAxisEnd3D)='
                    r'"(-?[\d.]+) (-?[\d.]+) (-?[\d.]+)"')


def _translate_body(body: str, dx: float, dy: float) -> str:
    s = translate(body, dx, dy)  # handles p= and BoundingBox
    s = _3D_RE.sub(
        lambda m: f'{m.group(1)}="{float(m.group(2)) + dx:.2f} '
                  f'{float(m.group(3)) + dy:.2f} {m.group(4)}"',
        s,
    )
    return s


# Objects ChemDraw gives a drawing order to. <s> is a run inside a <t>, not an
# object, and does not get one.
_Z_TAGS = ("fragment", "n", "b", "t", "graphic", "arrow", "curve", "group")
_Z_OPEN = re.compile(r"<(" + "|".join(_Z_TAGS) + r")\b(?![^>]*\bZ=)")


def stamp_z(body: str, start: int = 1):
    """Give every object an explicit drawing order. Returns (body, next free Z).

    DOCUMENT ORDER IS NOT DRAWING ORDER. The natural assumption -- emit the thing
    that should be behind first -- is wrong: ChemDraw sorts by Z and ignores
    where the element sits in the file. Verified on a test file where three
    coloured discs emitted BEFORE the fragment still covered it, and again in a
    hand-corrected version of that file, where the structure is written FIRST
    and carries the HIGHEST Z (94, its atoms and bonds 95-121) while the discs
    behind it sit at 39-83.

    Z is flat and global: a fragment's nodes and bonds each get their own, above
    the fragment's. So a structure is only reliably in front if everything inside
    it is numbered above everything outside it, which is why this stamps children
    and not just containers.

    Only called when a figure actually has decorations. Every other figure this
    repo writes carries no Z at all and renders correctly, and stamping those
    would rewrite every golden to say something they never needed to say.
    """
    n = [start]

    def bump(m):
        z = n[0]
        n[0] += 1
        return f"<{m.group(1)} Z=\"{z}\""

    return _Z_OPEN.sub(bump, body), n[0]


def translate_objects(body: str, dx: float, dy: float) -> str:
    """translate() PLUS the absolute 3D endpoints some objects carry.

    Use this, not translate(), whenever the string may hold an arrow or a disc:
    those keep Center3D / MajorAxisEnd3D / MinorAxisEnd3D in absolute page
    coordinates, and moving the p= and BoundingBox while leaving those behind
    tears the object in half. place() cannot be used instead when several pieces
    must move TOGETHER -- it derives its shift from each piece's own bbox.
    """
    return _translate_body(body, dx, dy)


def _content_extent(body: str):
    xs, ys = [], []
    for m in _P_RE.finditer(body):
        xs.append(float(m.group(1))); ys.append(float(m.group(2)))
    for m in _BB_RE.finditer(body):
        g = [float(v) for v in m.groups()]
        xs += [g[0], g[2]]; ys += [g[1], g[3]]
    for m in _3D_RE.finditer(body):
        xs.append(float(m.group(2))); ys.append(float(m.group(3)))
    return (min(xs), min(ys), max(xs), max(ys)) if xs else None


def document(header_prefix: str, body: str, margin: float = 24.0, page_mode: str = "fit") -> str:
    """Finalize a body into a full CDXML doc: shift all content into the positive
    quadrant (top-left at (margin, margin)) and wrap the page BoundingBox tightly.

    NOTE: text extent is approximated by the <t> anchor here (its p=), which
    under-counts a caption's real width/height. The margin absorbs that; the
    text-metrics model is what gives exact text boxes.
    """
    ext = _content_extent(body)
    if ext is None:
        return page(header_prefix, body, [0, 0, 100, 100])
    return finalize(header_prefix, body, ext, page_mode=page_mode, margin=margin)


# ----------------------------------------------------------------------
# Rich styled text (ported from the web engine's B.styledText)
# ----------------------------------------------------------------------
# Block opts style the whole <t>; inline markers style spans WITHIN a line:
#   ~x~ subscript(32)  ^x^ superscript(64,size-2)  *x* italic(2)  #x# bold(1)  _x_ underline(4)
#   {red}x{}  colour ONE SPAN, the rest of the line keeping the block colour
# face bits (combinable): bold1 italic2 underline4 outline8 shadow16 sub32 sup64.
# color is a colortable index (see beautify.extend_colortable). rot in degrees.
#
# THE COLOUR MARKER IS THE ONLY ASYMMETRIC ONE -- `{name}` opens and a bare `{}`
# closes -- because the other five need no argument and this one carries a colour
# name. Said out loud here and in the tool docs, because `{red}...{red}` is what
# a reader of the other five would try.
#
# It exists because a caption routinely changes colour MID-LINE and splitting the
# line into two text blocks is not equivalent: a published scheme can set eighteen
# er values as a red figure and a black remainder ("94" red, ":6 er" black), and two <t>
# objects would have to re-derive the spacing between them and would come apart
# the first time a font changed. The per-run mechanism is already here -- every
# marker above emits its own <s> -- so colour is applied per run as well.
#
# Matching the WHOLE span (open, text, close) rather than the delimiters
# separately is what makes `{{` safe: a doubled brace cannot start a token, so
# escaped text falls through to the literal branch, which un-doubles it. A
# literal brace INSIDE a coloured span is not expressible ([^{}]* excludes it) --
# noted rather than solved, because no real figure examined needed one.
_STYLE_MARK = re.compile(
    r"(~[^~]+~|\^[^\^]+\^|\*[^*]+\*|#[^#]+#|_[^_]+_|\{[^{}]+\}[^{}]*\{\})")


def styled_text(tid: int, cx: float, top_y: float, lines, size: int = 10, font: int = 21,
                color: int | None = None, bold: bool = False, italic: bool = False,
                underline: bool = False, outline: bool = False, shadow: bool = False,
                align: str = "center", rot: float = 0.0, line_height="auto") -> str:
    block_face = ((1 if bold else 0) | (2 if italic else 0) | (4 if underline else 0)
                  | (8 if outline else 0) | (16 if shadow else 0))

    def runs_of(line: str, col=None, face=None, sz=None) -> str:
        """MARKERS NEST, and they have to. One flat pass -- split the line, give
        each piece ONE style, done -- would turn `#{blue}Table 2{}#` -- a bold
        coloured label, which is how a typical figure title opens -- into one with
        the braces PRINTED, because the `#...#` run swallows the colour marker as
        literal text. Same for a marker inside a coloured span: `{teal}see #1#
        and #3#{}` would print the hashes.

        And it would fail by printing the markup, which is the worst way for a
        text feature to fail: nothing raises, nothing is reported, and the figure
        has `{20}` in the title.

        Each branch RECURSES into its own contents carrying the style it just
        added, so the styles compose. Depth is bounded by the markup itself.
        A flat line takes the flat path and comes out byte for byte the same as a
        single-pass split would -- the goldens are what says so.

        One thing does NOT nest: a colour inside a colour. The marker's own regex
        forbids braces in its body, so the inner one is not seen as a marker.
        One colour per span is also all anyone means.
        """
        face = block_face if face is None else face
        sz = size if sz is None else sz
        out = []
        for p in _STYLE_MARK.split(line):
            if p == "":
                continue
            if p.startswith("{") and p.endswith("{}"):
                name, txt = p[1:].split("}", 1)
                # A name from COLORS, or a raw colortable index as digits. Both
                # go through color_index below, which is the one place that knows
                # the difference.
                nc = int(name) if name.strip().isdigit() else name.strip()
                out.append(runs_of(txt[:-2], nc, face, sz))
            elif p.startswith("~") and p.endswith("~"):
                out.append(runs_of(p[1:-1], col, face | 32, sz))
            elif p.startswith("^") and p.endswith("^"):
                out.append(runs_of(p[1:-1], col, face | 64, max(6, sz - 2)))
            elif p.startswith("*") and p.endswith("*"):
                out.append(runs_of(p[1:-1], col, face | 2, sz))
            elif p.startswith("#") and p.endswith("#"):
                out.append(runs_of(p[1:-1], col, face | 1, sz))
            elif p.startswith("_") and p.endswith("_"):
                out.append(runs_of(p[1:-1], col, face | 4, sz))
            else:
                # Not a marker, so it is literal text -- and the one place a
                # doubled brace turns back into a single one. Done here rather
                # than up front so it cannot corrupt a marker being parsed.
                txt = p.replace("{{", "{").replace("}}", "}")
                catt = f' color="{color_index(col)}"' if col is not None else ""
                fc = f' face="{face}"' if face else ""
                out.append(f'<s font="{font}" size="{sz}"{catt}{fc}>{_esc(txt)}</s>')
        return "".join(out)

    arr = lines if isinstance(lines, list) else [lines]
    # `color` is one value for the whole block, or a list parallel to `lines` -- a caption
    # routinely mixes colours by line ("B~2~(pin)~2~" green above "*via C-C coupling*" in
    # black), and splitting that into two text blocks would break the line spacing.
    per_line = color if isinstance(color, (list, tuple)) else None
    if per_line is not None:
        # WHICH ELEMENT, and nothing else. What counts as a legal colour is
        # color_index's question and is asked there -- repeating the type rule
        # here would mean two places to keep in step, and the next palette change
        # would only update one of them. All this adds is the subscript, because
        # "element 0 of above_color" is what the caller has to go and fix and the
        # bare message cannot say it.
        #
        # None is legal and skipped: it means that line keeps the block colour,
        # which is what `below_color: ["green", null]` is for.
        for _i, _v in enumerate(per_line):
            if _v is None:
                continue
            try:
                color_index(_v)
            except ValueError as _e:
                raise ValueError("per-line colour list, element [%d]: %s"
                                 % (_i, _e)) from None
    inner = ""
    for i, l in enumerate(arr):
        if i > 0:
            inner += f'<s font="{font}" size="{size}">\n</s>'
        inner += runs_of(l, (per_line[i] if per_line and i < len(per_line) else
                             (None if per_line else color)))
    just = {"left": "Left", "right": "Right", "center": "Center"}.get(align, "Center")
    # REDUCED MOD 360 BEFORE IT IS WRITTEN, and this is measured rather than tidiness.
    # A hand-drawn test file (one caption, three copies, rotated one way and the
    # other) came back with the ANTICLOCKWISE copy carrying
    # RotationAngle="20643840" -- 315 degrees, not -45. ChemDraw normalises to
    # [0, 360) and a negative value appears nowhere in 278 real ChemDraw files,
    # whose range is 18 to 342. Unreduced, `rot=-45` would emit "-2949120", which
    # is spelling no ChemDraw has ever produced, and unattested spelling does not
    # reach a file.
    #
    # `% 360.0` and not abs() or a clamp: -45 has to become 315, which is the same
    # rotation, and not 45, which is the mirror image.
    #
    # The `if rot` gate is what keeps every unrotated caller byte-identical --
    # 0.0 is falsy, so a figure that never asks for rotation still
    # writes no attribute at all. -360 and 360 also reduce to 0.0 and correctly
    # write nothing, since neither is a rotation.
    _rot = rot % 360.0 if rot else 0.0
    rota = f' RotationAngle="{round(_rot * 65536)}"' if _rot else ""
    return (f'<t id="{tid}" p="{cx:.2f} {top_y:.2f}" Justification="{just}"'
            f' LineHeight="{line_height}" InterpretChemically="no"{rota}>{inner}</t>')


def arrow_xy(aid: int, x0: float, y0: float, x1: float, y1: float) -> str:
    """Reaction arrow in ANY direction: head at (x1,y1), tail at (x0,y0).
    Basis for angled/elbow/via routes (see arrows.py)."""
    mnx, mxx = min(x0, x1), max(x0, x1)
    mny, mxy = min(y0, y1), max(y0, y1)
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    return (f'<arrow id="{aid}" BoundingBox="{mnx:.1f} {mny:.1f} {mxx:.1f} {mxy:.1f}"'
            ' FillType="None" ArrowheadHead="Full" ArrowheadType="Solid" HeadSize="1000"'
            ' ArrowheadCenterSize="875" ArrowheadWidth="250"'
            f' Head3D="{x1:.1f} {y1:.1f} 0" Tail3D="{x0:.1f} {y0:.1f} 0" Center3D="{cx:.1f} {cy:.1f} 0"/>')


# ----------------------------------------------------------------------
# Scale a fragment (shrink an over-wide reagent to fit an arrow)
# ----------------------------------------------------------------------
def scale_fragment(frag: str, factor: float, cx: float | None = None, cy: float | None = None) -> str:
    """Scale all coordinates (and the fragment BoundingBox) by `factor` about
    (cx, cy). Default center is the fragment's bbox center. Shrinks bond lengths
    too - used to fit a reagent structure within its arrow span."""
    bb = bbox_of(frag)
    if cx is None or cy is None:
        cx, cy = ((bb[0] + bb[2]) / 2, (bb[1] + bb[3]) / 2) if bb else (0.0, 0.0)

    def sp(m):
        x = (float(m.group(1)) - cx) * factor + cx
        y = (float(m.group(2)) - cy) * factor + cy
        return f'p="{x:.2f} {y:.2f}"'

    s = _P_RE.sub(sp, frag)

    def sb(m):
        v = [float(x) for x in m.groups()]
        return 'BoundingBox="%.2f %.2f %.2f %.2f"' % (
            (v[0] - cx) * factor + cx, (v[1] - cy) * factor + cy,
            (v[2] - cx) * factor + cx, (v[3] - cy) * factor + cy)

    return _BB_RE.sub(sb, s)


# ----------------------------------------------------------------------
# Line-based box and square brackets (Line renders in SVG; the
# Rectangle graphic survives round-trip but is NOT exported to SVG).
# ----------------------------------------------------------------------
LINE_TYPES = ("Solid", "Dashed", "Bold", "Dashed Bold")


def _line(gid: int, x0: float, y0: float, x1: float, y1: float, lw: float = 1.0,
          color=None, line_type: str | None = None) -> str:
    """One <graphic GraphicType="Line">.

    `color` and `line_type` APPEND, and only when given, so a call that passes
    neither emits unchanged bytes -- which the scope and cycle goldens pin line
    by line.

    A scan that finds neither on any Line in 266 real files does not mean they
    are unavailable: none of those files had been drawn with one. A hand-drawn
    reference file carries `LineType="Dashed"` and a `color`, and
    ChemDraw's own documents carry `Dashed Bold` -- a COMBINATION the DTD's
    enumeration forbids, which is why LINE_TYPES is a tuple here rather than a
    promise that the DTD is complete.
    """
    col = "" if color is None else ' color="%d"' % color_index(color)
    lt = "" if not line_type else ' LineType="%s"' % line_type
    return (f'<graphic id="{gid}" BoundingBox="{x0:.2f} {y0:.2f} {x1:.2f} {y1:.2f}"'
            f'{col} GraphicType="Line" LineWidth="{lw}"{lt}/>')


# RectangleType is a SPACE-SEPARATED FLAG LIST, and the DTD is no help: it
# declares the attribute CDATA and says nothing. These are the thirteen
# combinations ChemDraw actually writes -- nine off a hand-drawn reference
# file, four more mined from the 652 documents ChemDraw ships. A base word,
# optionally with RoundEdge
# and/or Bold:
#
#     Plain  Dashed  Filled  Shaded  Shadow        <- base
#     RoundEdge (needs CornerRadius)  Bold         <- modifiers
#
# "Shaded" is NOT flat: it is the soft-edged variant, the same distinction
# disc()'s `shaded` draws. Scope was pinned to FLAT fill, so prefer "Filled".
RECTANGLE_TYPES = (
    "Plain", "Dashed", "Filled", "Shaded", "Shadow",
    "RoundEdge", "RoundEdge Dashed", "RoundEdge Filled", "RoundEdge Shaded",
    "RoundEdge Shadow", "RoundEdge Dashed Bold", "Shadow Bold", "Shadow Shaded",
)

# What ChemDraw wrote for both, in every rectangle of every file scanned. Neither
# is a proportion of anything -- a 600 corner on a 25 pt tall block and a 600
# corner on a 124 pt one are the same radius on paper.
CORNER_RADIUS = 600
SHADOW_SIZE = 400


def rect(gid: int, x0: float, y0: float, x1: float, y1: float, *,
         kind: str = "Plain", color=None, radius: int = CORNER_RADIUS,
         shadow: int = SHADOW_SIZE) -> str:
    """A native <graphic GraphicType="Rectangle"> -- frame, or filled block.

    A scan of 266 real files finds no Rectangle at all, so a fill primitive can
    also be built out of a four-atom ring coloured to the background
    (shaded_band), which works and is what the table engine uses. The scan only
    says nobody had drawn a box by hand in those files. This is the real object,
    and it does what the ring cannot -- rounded corners.

    Geometry is written the way ChemDraw writes it: a BoundingBox AND the centre
    with two axis ends, the same triple disc() carries. All three move under
    _translate_body (they are in _3D_RE), so a rectangle survives finalize; they
    are NOT touched by the plain translate(), so do not put one through place().

    NOT a replacement for shaded_band. That one is a fragment and rides along with
    C.place(); this one is a page object at fixed coordinates. Which to reach for
    is decided by whether the thing moves, not by which looks tidier.
    """
    if kind not in RECTANGLE_TYPES:
        raise ValueError(
            "RectangleType %r has no precedent. Known: %s. RectangleType is CDATA "
            "in the DTD, so this list is what ChemDraw was OBSERVED to write, not "
            "what it will accept -- a new combination needs one hand-drawn file "
            "before it goes in here." % (kind, ", ".join(RECTANGLE_TYPES)))
    cx, cy = (x0 + x1) / 2.0, (y0 + y1) / 2.0
    col = "" if color is None else ' color="%d"' % color_index(color)
    extra = ""
    if "RoundEdge" in kind:
        extra += ' CornerRadius="%d"' % radius
    if "Shadow" in kind:
        extra += ' ShadowSize="%d"' % shadow
    return (f'<graphic id="{gid}"'
            f' BoundingBox="{x0:.2f} {y0:.2f} {x1:.2f} {y1:.2f}"{col}'
            f' GraphicType="Rectangle" RectangleType="{kind}"{extra}'
            f' Center3D="{cx:.2f} {cy:.2f} 0"'
            f' MajorAxisEnd3D="{max(x0, x1):.2f} {cy:.2f} 0"'
            f' MinorAxisEnd3D="{cx:.2f} {max(y0, y1):.2f} 0"/>')


def disc(gid: int, cx: float, cy: float, rx: float, color, shaded: bool = False,
         ry: float | None = None) -> str:
    """A filled circle. THE decoration primitive: a dot on an atom and a glow
    behind a group are this, at two radii.

    Emitted as an ORBITAL, which is not a guess and not a pun -- it is the object
    ChemDraw itself produced when both of them were drawn by hand in a
    reference file. There is no filled rectangle or filled curve in the real
    files examined to copy instead: across them, FillType is "None" 1053 times
    and nothing else, so an orbital is the one filled shape with a verified
    spelling.

    shaded=True is OvalType "Circle Shaded" -- soft-edged, which is what makes a
    large one read as a glow rather than as a disc laid over the structure.

    Geometry is centre + two axis endpoints, and BoundingBox is written the
    degenerate way ChemDraw writes it (centre to the major-axis end, zero height)
    rather than the honest extent, because that is the form known to round-trip.
    Consequence worth knowing: document()'s page extent under-counts a disc by up
    to r, which the default margin absorbs -- and a halo sits behind a structure
    whose own box is bigger anyway.

    DRAW ORDER IS THE CALLER'S JOB. A glow has to be emitted BEFORE the fragment
    it sits behind; last drawn is on top, the same rule that decides which of two
    crossing bonds gets the break.
    """
    # ry != rx gives an ELLIPSE, and it costs nothing: the two axis endpoints are
    # already separate attributes. The four hand-drawn discs are all circles
    # only because a circle is what the tool draws by default -- and a circle is
    # the wrong shape for most cores, which are wider than they are tall, so a
    # round glow round a coumarin is mostly empty above and below it.
    # UNVERIFIED: OvalType still says "Circle" here, because that is the string
    # ChemDraw wrote and there is no non-circular orbital in any real file to copy.
    # If a stretched one comes back round, this is the attribute to look at.
    ry = rx if ry is None else ry
    oval = "Circle Shaded" if shaded else "Circle Filled"
    kind = "sShaded" if shaded else "sFilled"
    return (f'<graphic id="{gid}" BoundingBox="{cx + rx:.2f} {cy:.2f} {cx:.2f} {cy:.2f}"'
            f' color="{color_index(color)}" GraphicType="Orbital"'
            f' OvalType="{oval}" OrbitalType="{kind}"'
            f' Center3D="{cx:.2f} {cy:.2f} 0"'
            f' MajorAxisEnd3D="{cx + rx:.2f} {cy:.2f} 0"'
            f' MinorAxisEnd3D="{cx:.2f} {cy + ry:.2f} 0"/>')


# ----------------------------------------------------------------------
# <curve>: the Bezier primitive
# ----------------------------------------------------------------------
# Source is the published CDX/CDXML specification, NOT a hand-drawn file, and
# the spec is enough here because Curve is one of the documented objects -- it
# is the ARROW element the documentation drops, which is what the "hand-drawn
# CDXML is the source of truth" rule was really about.
#
# Two lines of that spec decide the whole design:
#
#   "May be contained by: Page, Group, FRAGMENT" -- so a curve written inside the
#   <fragment> is sanctioned, and C.translate() already carries CurvePoints along
#   (that is why _CURVE_RE exists). A decoration that lives there is moved by
#   place() for free, which is what halo / dot / locants cannot do.
#
#   CurvePoints holds CUBIC BEZIER CONTROL POINTS, not polyline vertices. Feeding
#   it sampled points of a sine draws something that merely passes near a sine.
#
# CurveType is a bitmask and 0 is Plain, which is all a squiggle wants.

# One cubic Bezier per HALF-PERIOD of the sine, in a frame where the half-period
# runs 0..1 along the spine and the amplitude is 1:
#
#     P0=(0,0)   P1=(S, A)   P2=(1-S, A)   P3=(1,0)
#
# The two control points are mirror images because a sine half-arch is symmetric
# about its peak, so the fit has exactly two free numbers. They were FITTED, once,
# by minimising the maximum vertical deviation from sin(pi*x) -- not guessed and
# not eyeballed. The solution equioscillates (5 interior extrema, all at exactly
# +/-0.0027), which is what says it is the true minimax and not a local dip.
#
# Worth knowing what these are NOT: A is 1.337, well above the peak the curve
# reaches. A control point is pulled at, not passed through, so the drawn crest
# lands on 1.000 -- do not "correct" A down to 1.0.
SQUIGGLE_CTRL_S = 0.412365      # along the spine, as a fraction of a half-period
SQUIGGLE_CTRL_A = 1.336955      # across it, as a fraction of the amplitude
# What that fit costs, as a fraction of the amplitude. Pinned so a later change to
# the two constants above has to move this one too (a test checks it).
SQUIGGLE_FIT_MAXDEV = 0.00272


# Arrowheads a <curve> may carry. Straight from the DTD, and unusually well
# corroborated for a style vocabulary: ChemDraw's own 652 documents write Full
# 433 times, HalfLeft 30 and HalfRight 4.
#
# **HalfLeft / HalfRight IS THE FISHHOOK.** A single-barbed head is what a
# mechanism draws for a one-electron arrow (homolysis), and it needs no new
# mechanism -- only an attribute value.
CURVE_HEADS = ("Unspecified", "None", "Full", "HalfLeft", "HalfRight")
CURVE_HEAD_TYPES = ("Solid", "Hollow", "Angle")


def curve(cid: int, points, color=None, curve_type: int = 0,
          head: str | None = None, tail: str | None = None,
          head_type: str | None = None, head_size: int | None = None,
          line_width: float | None = None) -> str:
    """A <curve> from a flat list of cubic-Bezier CONTROL points.

    `points` is [(x, y), ...] of length 3n+1 for n chained cubics: a start point,
    then (control, control, end) per segment. Any other length is a caller bug
    that would draw a curve of the wrong degree, so it raises rather than emits.

    The BoundingBox is the extent of the DRAWN curve (sampled), not of the control
    polygon, which sits well outside it -- and it is worth emitting rather than
    omitting: document()'s page extent reads BoundingBox but not CurvePoints, so
    without one a curve that overhangs its structure would fall off the page.

    No LineWidth BY DEFAULT: the bonds do not carry one either, so both inherit the
    document's 0.60 and a squiggle stays the same weight as the bond it crosses.
    `line_width` appends one only when asked, so a call without it emits
    unchanged bytes -- a test compares these character by character. A closed
    polygon drawn as an outline is the caller that needs it: it is furniture, not
    chemistry, and it is drawn next to a rule rather than a bond.
    """
    if len(points) < 4 or len(points) % 3 != 1:
        raise ValueError("CurvePoints holds cubic Bezier control points, so the "
                         "count must be 3n+1; got %d" % len(points))
    drawn = curve_polyline(points)
    xs = [p[0] for p in drawn]
    ys = [p[1] for p in drawn]
    flat = " ".join("%.2f %.2f" % (x, y) for x, y in points)
    col = "" if color is None else ' color="%d"' % color_index(color)
    # Appended, and only when asked, so a call without it emits unchanged bytes
    # -- a test compares these character by character.
    arrow = ""
    for name, val, allowed in (("ArrowheadHead", head, CURVE_HEADS),
                               ("ArrowheadTail", tail, CURVE_HEADS),
                               ("ArrowheadType", head_type, CURVE_HEAD_TYPES)):
        if val is None:
            continue
        if val not in allowed:
            raise ValueError("curve %s=%r; the DTD enumerates %s"
                             % (name, val, " | ".join(allowed)))
        arrow += ' %s="%s"' % (name, val)
    if head_size is not None:
        arrow += ' HeadSize="%d"' % head_size
    lw = "" if line_width is None else ' LineWidth="%.2f"' % float(line_width)
    return (f'<curve id="{cid}"{col} CurveType="{curve_type}"'
            f' BoundingBox="{min(xs):.2f} {min(ys):.2f} {max(xs):.2f} {max(ys):.2f}"'
            f'{lw}{arrow} CurvePoints="{flat}"/>')


def curve_polyline(points, per: int = 24) -> list[tuple[float, float]]:
    """Decode CurvePoints back into a polyline: `per` samples per cubic segment.

    THE canonical reader for what curve() writes, the same way arc_points() is for
    arc_arrow() -- an extent, a test and a reservation all need the drawn curve
    rather than the control polygon, and three separate de Casteljau loops would
    drift apart the first time one of them changed.
    """
    out = []
    for k in range(0, len(points) - 1, 3):
        (x0, y0), (x1, y1), (x2, y2), (x3, y3) = points[k:k + 4]
        for i in range(per + 1):
            if i == 0 and k:
                continue                     # shared with the previous segment
            t = i / per
            m = 1.0 - t
            out.append((m * m * m * x0 + 3 * m * m * t * x1
                        + 3 * m * t * t * x2 + t * t * t * x3,
                        m * m * m * y0 + 3 * m * m * t * y1
                        + 3 * m * t * t * y2 + t * t * t * y3))
    return out


def squiggle_points(x0: float, y0: float, x1: float, y1: float, *,
                    amp: float, wavelength: float, overhang: float):
    """Control points for a squiggle drawn ACROSS the bond (x0,y0)-(x1,y1).

    The bond-cleavage mark: its spine runs PERPENDICULAR to the bond, centred on
    the midpoint, and it waves back and forth ALONG the bond direction.

    `amp` / `wavelength` / `overhang` arrive in page units, already multiplied out
    of their bond-length fractions by the caller.

    HOW LONG. Half-spine = (this bond's own length)/2 + overhang. The crossing
    part is measured off the real bond so a long bond is still crossed end to end;
    the overhang is a fixed distance so every squiggle in a figure sticks out by
    the same amount.

    WAVELENGTH IS SNAPPED to fit a whole number of half-periods on each side of
    the bond. Left exact, the spine would end mid-swing -- which reads as a
    truncated drawing rather than as a mark -- and the end segments would have to
    be split mid-Bezier. Snapping moves it by at most half a half-period and puts
    both ends, and the bond crossing itself, on a zero of the sine.
    """
    dx, dy = x1 - x0, y1 - y0
    blen = math.hypot(dx, dy)
    if blen == 0:
        raise ValueError("a squiggle needs two distinct bond endpoints")
    ux, uy = dx / blen, dy / blen            # along the bond: the wave's swing
    nx, ny = -uy, ux                         # across it: the spine
    mx, my = (x0 + x1) / 2.0, (y0 + y1) / 2.0
    half = blen / 2.0 + overhang
    k = max(1, int(round(2.0 * half / wavelength)))   # half-periods PER SIDE
    hp = half / k                                     # snapped half-period
    S, A = SQUIGGLE_CTRL_S, SQUIGGLE_CTRL_A

    def at(s, d):
        return (mx + nx * s + ux * d, my + ny * s + uy * d)

    pts = [at(-half, 0.0)]
    for j in range(2 * k):
        s = -half + j * hp
        sign = -1.0 if (j + k) % 2 else 1.0            # sin(pi*(j-k+t))
        pts.append(at(s + S * hp, sign * A * amp))
        pts.append(at(s + (1.0 - S) * hp, sign * A * amp))
        pts.append(at(s + hp, 0.0))
    return pts


# Half the ACS header's LineWidth 0.60 -- the amount a hand-drawn file sets
# the fragment BoundingBox outside its own corner atoms.
_BAND_MARGIN = 0.30


def shaded_band(base: int, x0: float, y0: float, x1: float, y1: float,
                color: int) -> str:
    """A filled rectangle with no visible border. THE shading primitive.

    NOT a Rectangle graphic, and not a guess. Across the real CDXML files
    examined, FillType is "None" 1056 times and nothing else -- there is no
    filled rectangle anywhere to copy a spelling from. So a shaded band was drawn
    by hand in ChemDraw, and ChemDraw wrote THIS:

        a four-atom ring at the corners, every atom and bond color="1", which
        is the background -- that is how the outline becomes invisible -- and a
        <ColoredMolecularArea> inside the fragment naming the four BONDS as its
        BasisObjects, with bgcolor an index into the colortable.

    The shading is therefore a property of a molecule, not a shape you draw,
    which is why nothing resembling it turned up while looking for rectangles.

    `color` is a colortable index. The ACS eight are all saturated, so a band
    colour is appended with beautify.extend_colortable first (indices from 10);
    The hand-drawn file carries its grey at 10 and its blue at 11 for exactly
    this reason.

    Consumes TEN ids from `base`: four atoms, four bonds, the fragment, the area.

    UNVERIFIED, and worth one look on paper: whether ChemDraw pads the shaded
    area outwards from the atoms it is given. If it does, a band drawn to a row's
    exact height renders slightly taller than the row.

    DRAW ORDER IS Z, NOT DOCUMENT ORDER (see stamp_z). A band must carry a LOWER
    Z than the text it sits behind; emitting it earlier is not enough.
    """
    n = (base, base + 1, base + 2, base + 3)
    b = (base + 4, base + 5, base + 6, base + 7)
    m = _BAND_MARGIN
    corners = ((x0, y0), (x0, y1), (x1, y1), (x1, y0))
    nodes = "".join(f'<n id="{i}" p="{x:.2f} {y:.2f}" color="1" AS="N"/>'
                    for i, (x, y) in zip(n, corners))
    bonds = "".join(f'<b id="{i}" color="1" B="{p}" E="{q}" BS="N"/>'
                    for i, (p, q) in zip(b, ((n[0], n[1]), (n[1], n[2]),
                                             (n[2], n[3]), (n[3], n[0]))))
    return (f'<fragment id="{base + 8}" BoundingBox="{x0 - m:.2f} {y0 - m:.2f}'
            f' {x1 + m:.2f} {y1 + m:.2f}">{nodes}{bonds}'
            + colored_area(base + 9, color, b) + '</fragment>')


def colored_area(aid: int, color: int, bond_ids) -> str:
    """ChemDraw's native highlight: a coloured area over the BONDS it names.

    THE tag for a filled region, and it is not a shape -- it is a property of a
    molecule, which is why nothing resembling a filled rectangle turns up when you
    go looking for one. Two callers, deliberately the same emitter:

      shaded_band()  hangs one on a throwaway four-bond ring to shade a table row;
      decorate.apply's "fill" hangs one on bonds of the REAL molecule.

    Verified against a file hand-drawn with ChemDraw's own
    Color > Highlight: three of these, inside <fragment>, after the nodes and
    bonds, each naming bond ids and carrying a bgcolor that indexes the document
    colortable (its bgcolor 5 / 10 / 11 against a 10-entry table are yellow /
    orange / pink under the "0 and 1 reserved, entries start at 2" rule).

    `bond_ids` are ids, NOT endpoint pairs -- BasisObjects names objects. They are
    renumbered with everything else because BasisObjects is in _ID_LIST_ATTRS.
    """
    ids = " ".join(str(i) for i in bond_ids)
    return f'<ColoredMolecularArea id="{aid}" bgcolor="{color}" BasisObjects="{ids}"/>'


_NODE_P_RE = re.compile(r'<n\b[^>]*?\bid="(\d+)"[^>]*?\bp="(-?[\d.]+) (-?[\d.]+)"')


def node_points(frag: str) -> dict:
    """{node id -> (x, y)} for a fragment string.

    How a decoration gets from "these atoms" to "this place": the anchor names
    ids, fragment() has already committed them to coordinates, and this reads
    them back rather than recomputing a layout that would not match.
    """
    return {int(m.group(1)): (float(m.group(2)), float(m.group(3)))
            for m in _NODE_P_RE.finditer(frag)}


_NODE_TAG_RE = re.compile(r'<n\b[^>]*>')
_ATTR_RE = re.compile(r'([A-Za-z_][\w]*)="([^"]*)"')


def node_props(frag: str) -> dict:
    """{node id -> {attribute: value}} for a fragment string, opening tags only.

    The sibling of node_points, for a decoration that has to treat one atom
    differently from another rather than only place something at it -- an arrow
    stops further short of an O or an H than of a ring carbon, because what it
    must not touch is the drawn LABEL and a plain carbon has none.

    Values stay strings; the caller knows which ones are numbers. CDXML omits
    Element for carbon (real files never write Element="6"), so a missing
    attribute is not missing information here -- it means 6, and the caller's
    .get default says so.
    """
    out = {}
    for m in _NODE_TAG_RE.finditer(frag):
        a = dict(_ATTR_RE.findall(m.group(0)))
        if "id" in a:
            out[int(a["id"])] = a
    return out


_BOND_RE = re.compile(r'<b\b[^>]*?\bB="(\d+)"[^>]*?\bE="(\d+)"')
_BY_NUMBER = {v: k for k, v in _ELEMENTS.items()}


def bond_segments(frag: str) -> list:
    """[(i, j, (x0, y0), (x1, y1)), ...] -- every drawn bond, with its node ids.

    For anything that has to route AROUND a structure rather than only land on
    it. node_points says where the atoms are; this says what is strung between
    them, which is what an electron-pushing arrow must not be drawn across.

    The ids come back because a router has to EXEMPT the bonds the arrow is
    anchored to: an arrow whose tail is a bond's midpoint necessarily touches
    that bond, and charging it for that leaves no legal route at all.
    """
    pts = node_points(frag)
    out = []
    for m in _BOND_RE.finditer(frag):
        i, j = int(m.group(1)), int(m.group(2))
        if i in pts and j in pts:
            out.append((i, j, pts[i], pts[j]))
    return out


def keepout_radii(frag: str) -> dict:
    """{node id -> radius} -- how much room each atom's DRAWN form takes.

    A labelled atom is a glyph and an arc through it is unreadable; a plain
    carbon is a corner and an arc may pass close by. So the radius is the label's
    own half extent where there is a label (_label_overhang, the same numbers the
    setbacks were derived from) and a small constant where there is not.
    """
    props = node_props(frag)
    out = {}
    for nid, pr in props.items():
        if pr.get("NodeType"):                       # contracted group, R, ...
            out[nid] = 0.55 * BOND_LENGTH
            continue
        anum = int(pr.get("Element", 6))
        chg = int(pr.get("Charge") or 0)
        if anum == 6 and not chg:
            out[nid] = 0.12 * BOND_LENGTH            # a bare vertex
            continue
        sym = _BY_NUMBER.get(anum, "C")
        _along, half_sym, half_h = _label_overhang(sym, 0)
        out[nid] = max(half_sym, half_h)
    return out


# A DELOCALIZED CHARGE, drawn the way a hand-drawn reference draws it: no double
# bonds in the delocalized part, one DASHED ARC spanning it, and one circled
# charge symbol in the middle, belonging to no single atom.
#
# EVERY NUMBER MEASURED off that hand-drawn reference file (its third
# structure). The important find there is that ChemDraw REPLACED the hand-drawn
# <graphic GraphicType="Arc"> with a native <arrow> of its own
# (SupersededBy="82"), and that arrow carries LineType="Dashed" and NO
# ArrowheadHead attribute at all. So the faithful form is
# arc_arrow(..., head=None, line="Dashed"), which emits that attribute set
# byte-for-byte -- not a <graphic>, which is only what the authoring tool wrote.
#   arc centre    the delocalized atoms' CENTROID (the hand-drawn arc is 1.4
#                 units off it)
#   arc radius    0.766 of the mean centroid-to-atom distance (11.03 of 14.40)
#   sweep         237.5 deg, where the five atoms themselves span 234.1, i.e.
#                 about 1.7 deg of margin past the outermost atom at each end
#   symbol        CirclePlus, box 7.50 wide and zero high, and NO <represent> --
#                 which is why the file carries ChemDraw's own
#                 Warning="Chemical symbol is not associated with an atom."
#                 A delocalized charge belongs to no atom, so that warning is
#                 correct rather than something to design away.
# THE ARC'S RADIUS IS SET FROM THE BONDS, not from the atoms. Measured off the
# hand drawing: radius 11.03 where the ring's bond lines run 12.47 from its centre,
# i.e. only 1.44 units -- 0.10 bond -- of daylight, which is why an arc at that
# radius reads as sitting ON the ring. Scaling the APOTHEM (the
# centre-to-bond-line distance) instead makes the clearance a stated quantity:
# radius = DELOC_ARC_RADIUS * apothem leaves (1 - DELOC_ARC_RADIUS) * apothem.
# WHERE 0.85 COMES FROM. Chosen, not measured -- picked by eye off a five-rung
# ladder of rendered values (0.55 / 0.65 / 0.75 / 0.85 / 0.885, i.e. 0.389 /
# 0.303 / 0.217 / 0.130 / 0.100 bond of daylight). The last rung is the value
# read off the hand-drawn reference file, and it was rejected as too tight: at
# that radius the dashed arc overlaps the six-membered ring. Re-render such a
# ladder before changing this.
DELOC_ARC_RADIUS = 0.85
DELOC_ARC_MARGIN = 1.7                 # degrees, each end
DELOC_SYMBOL = {1: "CirclePlus", -1: "CircleMinus"}
_B_TAG_RE = re.compile(r'<b\b[^>]*/>')


def singleize_bonds(frag: str, ids) -> str:
    """Drop Order="2"/"3" from every bond with BOTH ends in `ids`.

    The delocalized part of a ring is drawn with no double bonds at all -- the
    dashed arc is what says where the electrons are. Only 2 and 3 go: Order can
    also carry a KIND ("dative"), which is not a multiplicity and not this
    function's to remove.
    """
    keep = {int(i) for i in ids}

    def fix(m):
        tag = m.group(0)
        b = re.search(r'\bB="(\d+)"', tag)
        e = re.search(r'\bE="(\d+)"', tag)
        if b and e and int(b.group(1)) in keep and int(e.group(1)) in keep:
            return re.sub(r'\sOrder="[23]"', '', tag)
        return tag

    return _B_TAG_RE.sub(fix, frag)


def clear_atom_charges(frag: str, ids) -> str:
    """Take Charge off `ids`, and remove the per-atom charge symbols bound to them.

    For the delocalized depiction: the charge is no longer ON an atom, it is the
    one free-floating symbol in the middle of the arc, and the hand-drawn
    reference structure has no Charge attribute on any node at all. Left in
    place, the atom keeps its own bound Plus from _charge_symbols and the figure
    shows the charge TWICE.
    """
    keep = {int(i) for i in ids}
    out = frag
    for nid in keep:
        # NO \b after the closing quote: both sides of it are non-word
        # characters, so \b never matches there and the substitution would be
        # a no-op. The quote alone already stops id="6" matching id="60".
        out = re.sub(r'(<n id="%d"[^>]*?)\sCharge="-?\d+"' % nid, r'\1', out)
        out = re.sub(
            r'<graphic\b[^>]*>\s*<represent attribute="Charge" object="%d"\s*/>'
            r'\s*</graphic>\s*' % nid, '', out)
    return out


def delocalized_marks(points, aid_arc: int, aid_sym: int | None = None,
                      charge: int = 0, color=None, sym_points=None,
                      apothem: float = 0.0, centre=None,
                      radius_frac: float | None = None) -> list:
    """[arc_xml, symbol_xml] for a charge delocalized over `points`.

    The arc spans every point: its ends sit at the two outermost atoms and it
    runs the way that AVOIDS the largest angular gap, which for an arenium is the
    gap across the sp3 carbon. Nothing has to say which atoms are the ends.
    """
    n = len(points)
    if centre is not None:
        cx, cy = centre
    else:
        cx = sum(p[0] for p in points) / n
        cy = sum(p[1] for p in points) / n
    frac = DELOC_ARC_RADIUS if radius_frac is None else float(radius_frac)
    if apothem:
        rad = frac * apothem
    else:
        # No bonds handed in (a non-cyclic delocalized system): fall back to the
        # atoms' own distance, which is all there is to scale.
        rad = frac * (sum(math.hypot(p[0] - cx, p[1] - cy)
                          for p in points) / n)
    ang = sorted(math.degrees(math.atan2(p[1] - cy, p[0] - cx)) % 360.0
                 for p in points)
    gaps = [((ang[(i + 1) % n] - ang[i]) % 360.0, i) for i in range(n)]
    gap, k = max(gaps)
    t0 = ang[(k + 1) % n] - DELOC_ARC_MARGIN
    sweep = (360.0 - gap) + 2.0 * DELOC_ARC_MARGIN
    out = [arc_arrow(aid_arc, cx, cy, rad, rad, t0, sweep,
                     head=None, line="Dashed")]
    if color is not None:
        out[0] = out[0].replace('<arrow ', '<arrow color="%d" ' % color_index(color))
    sym = DELOC_SYMBOL.get(1 if charge > 0 else (-1 if charge < 0 else 0))
    if sym and aid_sym is not None:
        # The ARC is always centred on the delocalized atoms -- it has to pass
        # through them. The SYMBOL may not be: IUPAC GR-5.4 centres it "within the
        # delocalized system" (sym_points=None, the five carbons of an arenium),
        # and a house style may want the enclosing RING's centre instead. The two
        # are 0.3 bond apart on an arenium, so this is a real choice and not a
        # rounding difference.
        # THE SYMBOL SHARES THE ARC'S CENTRE when one was given. Without this the
        # + sits 2.88 units off the middle even with the box convention right:
        # the arc moves to the ring's centre while the symbol still averages the
        # five delocalized atoms.
        if sym_points:
            sx = sum(q[0] for q in sym_points) / len(sym_points)
            sy = sum(q[1] for q in sym_points) / len(sym_points)
        else:
            sx, sy = cx, cy
        # THE GLYPH SITS AT THE BOX'S FIRST COORDINATE, not at the box's centre.
        # Read off the hand-drawn reference structure: its floating box is
        # "315.17 59.87 307.67 59.87" and the six-ring's centroid is
        # (315.68, 60.10) -- 0.56 units from the FIRST pair and 4.26 from the
        # box's centre. Centring on the box puts the symbol 3.75 units (0.26 bond)
        # to the right of where it belongs -- visibly not in the middle.
        out.append('<graphic id="%d" BoundingBox="%.2f %.2f %.2f %.2f"'
                   ' GraphicType="Symbol" SymbolType="%s"/>'
                   % (aid_sym, sx, sy, sx - CHARGE_SYMBOL_WIDTH, sy, sym))
    return out


def next_free_id(frag: str) -> int:
    """One past the largest id already in a fragment string.

    For an object that has to live INSIDE the fragment. Its id cannot come from a
    decoration's own pool: place() renumbers a placed fragment by adding `base`,
    and two cells whose bases differ by less than that pool's offset would then
    hand out the same number. Taken from the fragment's own range, an id stays
    unique for exactly as long as the atom ids do.
    """
    ids = [int(m.group(1)) for m in re.finditer(r'\bid="(\d+)"', frag)]
    return max(ids) + 1 if ids else 1


def box_lines(gid: int, x0: float, y0: float, x1: float, y1: float, lw: float = 1.0) -> str:
    """A rectangle drawn as four Line graphics (top, right, bottom, left)."""
    return "\n".join([
        _line(gid, x0, y0, x1, y0, lw),
        _line(gid + 1, x1, y0, x1, y1, lw),
        _line(gid + 2, x1, y1, x0, y1, lw),
        _line(gid + 3, x0, y1, x0, y0, lw),
    ])


def bracket_lines(gid: int, x0: float, y0: float, x1: float, y1: float,
                  serif: float = 6.0, side: str = "both", lw: float = 1.0) -> str:
    """Square brackets around a region. side: 'both' | 'left' | 'right'."""
    out = []
    if side in ("both", "left"):
        out += [_line(gid, x0, y0, x0, y1, lw),
                _line(gid + 1, x0, y0, x0 + serif, y0, lw),
                _line(gid + 2, x0, y1, x0 + serif, y1, lw)]
    if side in ("both", "right"):
        out += [_line(gid + 3, x1, y0, x1, y1, lw),
                _line(gid + 4, x1, y0, x1 - serif, y0, lw),
                _line(gid + 5, x1, y1, x1 - serif, y1, lw)]
    return "\n".join(out)


# ----------------------------------------------------------------------
# ACS 1996 document style: Arial, BondLength 14.4, ACS line widths.
# Gives generated figures true ACS styling: font 70=Arial, BondLength 14.4,
# BondSpacing 18, ACS line widths/hash spacing. Palette indices: 2=white 3=black
# 4=red 5=yellow 6=green 7=cyan 8=blue 9=magenta; extend_colortable appends from 10.
# ----------------------------------------------------------------------
_ACS_HEADER = '<?xml version="1.0" encoding="UTF-8" ?>\n<!DOCTYPE CDXML SYSTEM "https://static.chemistry.revvitycloud.com/cdxml/CDXML.dtd" >\n<CDXML\n CreationProgram="ChemWand 1.0"\n Name="ChemWand"\n BoundingBox="122.53 118.62 161.77 190.70"\n WindowPosition="0 0"\n WindowSize="0 0"\n FractionalWidths="yes"\n InterpretChemically="yes"\n ShowAtomQuery="yes"\n ShowAtomStereo="no"\n ShowAtomEnhancedStereo="yes"\n ShowAtomNumber="no"\n ShowResidueID="no"\n ShowBondQuery="yes"\n ShowBondRxn="yes"\n ShowBondStereo="no"\n ShowTerminalCarbonLabels="no"\n ShowNonTerminalCarbonLabels="no"\n HideImplicitHydrogens="no"\n LabelFont="70"\n LabelSize="10"\n LabelFace="96"\n CaptionFont="70"\n CaptionSize="10"\n HashSpacing="2.49"\n MarginWidth="1.59"\n LineWidth="0.60"\n BoldWidth="2.01"\n BondLength="14.40"\n BondSpacing="18"\n ChainAngle="120"\n LabelJustification="Auto"\n CaptionJustification="Left"\n AminoAcidTermini="HOH"\n ShowSequenceTermini="yes"\n ShowSequenceBonds="yes"\n ShowSequenceUnlinkedBranches="no"\n ShowPrimeLabels="no"\n ResidueWrapCount="40"\n PrintMargins="36 36 36 36"\n MacPrintInfo="00030000004800480000000002F30243FFEDFFEF030502530367052803FC00020000004800480000000002F30243000100000064000000010001010100000001270F000100010000000000000000000000000002001901900000000000400000000000000000000100000000000000000000000000000000"\n ChemPropName=""\n ChemPropFormula="Chemical formula: "\n ChemPropExactMass="Exact mass: "\n ChemPropMolWt="Molecular weight: "\n ChemPropMOverZ="m/z: "\n ChemPropAnalysis="Elemental analysis: "\n ChemPropBoilingPt="Boiling point: "\n ChemPropMeltingPt="Melting point: "\n ChemPropCritTemp="Critical temp: "\n ChemPropCritPres="Critical pres: "\n ChemPropCritVol="Critical vol: "\n ChemPropGibbs="Gibbs energy: "\n ChemPropLogP="Log P: "\n ChemPropMR="MR: "\n ChemPropHenry="Henry&apos;s law: "\n ChemPropEForm="Heat of form: "\n ChemProptPSA="tPSA: "\n ChemPropID=""\n ChemPropFragmentLabel=""\n color="0"\n bgcolor="1"\n RxnAutonumberStart="1"\n RxnAutonumberConditions="no"\n RxnAutonumberStyle="Roman"\n RxnAutonumberFormat="(#)"\n MonomerRenderingStyle="graphic"\n><colortable>\n<color r="1" g="1" b="1"/>\n<color r="0" g="0" b="0"/>\n<color r="1" g="0" b="0"/>\n<color r="1" g="1" b="0"/>\n<color r="0" g="1" b="0"/>\n<color r="0" g="1" b="1"/>\n<color r="0" g="0" b="1"/>\n<color r="1" g="0" b="1"/>\n</colortable><fonttable>\n<font id="70" charset="x-mac-roman" name="Arial"/>\n</fonttable>'


# Colour indices for the ACS header's colortable. CDXML reserves 0 and 1 for background
# and foreground, so the eight <color> entries land at 2..9. Verified empirically for red:
# writing color="4" on an arrow rendered red after a ChemDraw round trip. The rest follow
# from that offset plus the header's entry order, and would each want one test file
# to confirm.
COLORS = {"white": 2, "black": 3, "red": 4, "yellow": 5,
          "green": 6, "cyan": 7, "blue": 8, "magenta": 9}


def color_index(c) -> int:
    """Accept a name from COLORS or a raw colortable index."""
    if isinstance(c, str):
        try:
            return COLORS[c.lower()]
        except KeyError:
            # NEITHER MESSAGE HERE ENUMERATES THE NAMES, and that is deliberate.
            # The eight below are the header's fixed entries, but a document may
            # append more, so a list printed from COLORS would be a confident,
            # incomplete answer -- the kind that sends the reader looking for a
            # spelling error in a name that is simply not in this figure's table.
            # Naming WHERE the names come from stays true as the table grows.
            raise ValueError(
                f"unknown colour {c!r}. This position takes a colour NAME from "
                "the document's palette, or a colortable INDEX (an integer)."
            ) from None
    if isinstance(c, (list, tuple)):
        # An [r, g, b] has to be APPENDED to the document's colortable before it
        # has an index, and only the tool layer holds the header to append it to.
        # int() on a list raises TypeError, which says nothing useful about why.
        #
        # Says WHY this function cannot do it rather than which tools can: a list
        # of tool names is a fact about the current wiring, and it goes stale the
        # moment another tool is wired up. The reason does not.
        raise ValueError(
            "an [r, g, b] colour cannot be used here. Pass a colour NAME or a "
            "colortable INDEX; an [r, g, b] belongs on a highlight spec's \"color\", "
            "which is added to the colortable for you.")
    if isinstance(c, float):
        # int() would TRUNCATE, and the truncation is silent and wrong in the
        # worst possible way: 1.0 becomes index 1, which is the RESERVED
        # FOREGROUND, so the text renders in the default colour and nothing says
        # the request was dropped. Measured on a three-line caption written with
        # above_color=[1.0, 0.60, 0.10] -- an RGB triple in a position where a
        # list already means one colour per LINE. The caller asked for orange and
        # got black.
        #
        # ONE PLACE, checked here rather than at each of the eleven call sites,
        # and verified first: instrumenting this function across the whole test
        # suite plus a tool-layer sweep saw only str and int arguments, so no
        # existing caller relied on the truncation.
        raise ValueError(
            f"colour {c!r} is a float. This position takes a colour NAME or a "
            "colortable INDEX (an integer). If you meant an [r, g, b] triple, "
            'that goes on a highlight spec\'s "color" -- the tool layer appends '
            "it to the document colortable and rewrites it to an index -- and "
            "NOT on a caption colour (above_color / below_color / yield_color), "
            "where a list already means one colour per line.")
    return int(c)


def acs_header_prefix() -> str:
    """ACS 1996 document header (the default for generated figures)."""
    return _ACS_HEADER


def native_bracket(gid: int, x0: float, y0: float, x1: float, y1: float,
                   bracket_type: str = "Square", lip: float = 60) -> str:
    """A native ChemDraw bracket PAIR around the region (x0,y0)=top-left .. (x1,y1)=bottom-right.
    Emitted as two <graphic GraphicType="Bracket"> objects (real, selectable brackets — NOT the
    three-line fake, and decorative so no atom-attachment markers). Direction is encoded by the
    y-order of each bracket's degenerate vertical BoundingBox (verified against a ChemDraw export).
    bracket_type: 'Square' | 'Round' | 'Curly'."""
    left = (f'<graphic id="{gid}" BoundingBox="{x0:.2f} {y1:.2f} {x0:.2f} {y0:.2f}"'
            f' GraphicType="Bracket" BracketType="{bracket_type}" LipSize="{lip}"/>')
    right = (f'<graphic id="{gid + 1}" BoundingBox="{x1:.2f} {y0:.2f} {x1:.2f} {y1:.2f}"'
             f' GraphicType="Bracket" BracketType="{bracket_type}" LipSize="{lip}"/>')
    return left + "\n" + right


# ----------------------------------------------------------------------
# Grouping: emit a ChemDraw <group> around already-built objects (Ctrl+G)
# ----------------------------------------------------------------------
def group(gid: int, children: list[str], bbox=None, integral: bool = False) -> str:
    """Wrap already-built child object strings in a ChemDraw <group> — the Ctrl+G
    equivalent, so the members select / move / align as a single object. `children`
    are element strings (e.g. a <fragment> from place() and a <t> from styled_text).
    The group's own BoundingBox (if given) and the nested children's p= / BoundingBox
    are all handled by translate() / _translate_body() at finalize, so grouping is
    transparent to the positive-quadrant shift.

    VALIDATED against ChemDraw. Two files answer it, and between them they cover
    both directions:

      * WHAT CHEMDRAW WRITES. A hand-made file is a real Ctrl+G of a structure
        plus a caption. Its group carries `id`, `BoundingBox` and `Z`
        and NOTHING ELSE — no Integral. Its children keep ABSOLUTE page
        coordinates (the fragment's box sits inside the group's, in the same
        frame), and the group's BoundingBox is exactly the union of its children's,
        to the last decimal.
      * WHAT CHEMDRAW DOES WITH THIS FUNCTION'S OUTPUT. A test file of this
        function's own output was opened in ChemDraw and saved back. All three groups survived AS
        GROUPS with their children still nested.

    So the form below is right. Two things ChemDraw does on the way through, both
    worth knowing and neither a reason to change anything here:

      * IT RECOMPUTES THE BOUNDING BOX. One went in at 246.74 209.30 385.26 238.10
        and came back 227.74 181.20 385.41 238.45. `bbox` is therefore advisory —
        a caller whose box is an approximation (a text-metrics estimate, say) is
        not shipping that approximation to the reader.
      * IT RENUMBERS A GROUP ID THAT COLLIDES WITH A CHILD'S, silently: 8000/8001/
        8002 came back 8000/8003/8006. Silent means a collision is not visible in
        the bytes, so callers still pick ids from a range of their own.

    Z IS NOT SET HERE. ChemDraw adds one on the way through, and in its own file a
    group's members carry their own Z with other objects interleaved between them
    — so the members' Z is what orders the drawing and the group's is only its own
    slot. Figures this pipeline writes carry no Z at all unless stamp_z runs.
    """
    inner = "\n".join(children)
    bb = ' BoundingBox="%.2f %.2f %.2f %.2f"' % tuple(bbox) if bbox is not None else ""
    integral_attr = ' Integral="yes"' if integral else ""
    return f'<group id="{gid}"{bb}{integral_attr}>\n{inner}\n</group>'


# ----------------------------------------------------------------------
# Page placement: fit / single-column / double-column / full US-Letter
# ----------------------------------------------------------------------
# CDXML units are points (1 unit = 1/72 in; BondLength 14.4 = 0.2 in). ACS figure
# column widths: single ~3.33 in (240 pt), double ~7 in (504 pt). Letter = 612x792.
_PAGE_W = {"single": 240.0, "double": 504.0, "letter": 612.0}

# US Letter's HEIGHT, and the bound for all three fixed-width modes rather than
# for 'letter' alone: an ACS column is a column ON this paper, so a figure taller
# than the sheet cannot be placed on one page whatever the journal's type block
# turns out to be. Deliberately conservative -- a figure that overruns the type
# block but not the paper stays quiet -- because the type block is a number this
# file does not have and would be inventing.
_PAGE_H = 792.0

# Figures whose content forced a FIXED-width page wider than the caller asked for.
# Entries are (page_mode, asked_pt, got_pt); the caller clears it per figure and
# reads it after. Out-parameter by module global rather than by argument for the
# reason the `shift` docstring below gives about re-deriving: finalize is the one
# place the page width is known, and it is called from all five layout modules
# plus table_layout and `document`, so threading a list through would be seven
# signatures for one number.
#
# WHY THIS IS NOT SILENT. 'single' means "ACS single column, 3.33 in",
# and the branch below deliberately GROWS the page instead of clipping. Measured:
# one sotorasib-sized structure in page='single' came out 370 pt wide -- 54% over
# the 240 pt column -- with no warning anywhere in the reply. A figure submitted
# at that width is rescaled by the journal, which shrinks every caption with it.
PAGE_OVERFLOW: list = []

# The smallest width difference that survives into the file: `page(...)` below is
# handed round(pw), so anything under half a point writes the same bytes. Read off
# the emitter rather than chosen -- a figure 0.2 pt over is not over.
_PAGE_OVER_EPS = 0.5


def usable_width(page_mode: str = "letter", margin: float = 30.0):
    """Content width available inside a fixed-width page (page width - 2*margin).
    'fit'/'auto' have no fixed page width, so return None (caller should not wrap)."""
    pw = _PAGE_W.get(page_mode)
    return (pw - 2 * margin) if pw else None


def finalize(header_prefix: str, body: str, content_bbox, page_mode: str = "fit",
             margin: float = 30.0, shift: list | None = None) -> str:
    """Place `body` on a page and return the full CDXML.
      page_mode: 'fit'    -> tight box around the content (default)
                 'single' -> ACS single-column width (240 pt), content centered
                 'double' -> ACS double-column width (504 pt), content centered
                 'letter' -> full US-Letter page (612x792), content centered both ways
                 'auto'   -> smallest standard column that fits, else fit-to-content
    The page NEVER clips: any dimension smaller than the content expands to content+margin.

    `shift` is an out-parameter (same idiom as `notes`, `index_map`, `decorations`):
    pass a list and the (ox, oy) applied here is appended to it.

    THE ONE PLACE THAT NUMBER IS KNOWN. Everything an engine reserves is in the
    engine's own layout space, which can start anywhere -- layout_cycle puts the
    ring centre at (0, 0), so its top-left is negative -- and this function is what
    moves it onto the page. So a caller holding an Allocator full of layout-space
    rects and a written document full of page-space ink has no way to compare them
    without this, and the alternative is re-deriving the four branches below from a
    worse vantage point. _add_graphics' docstring already refuses to do that.

    Deliberately NOT a return value: five engines return a document string and one
    of them (`document`) is called from a dozen places that want nothing else."""
    minx, miny, maxx, maxy = content_bbox
    cw, ch = maxx - minx, maxy - miny
    if page_mode == "letter":
        pw = max(_PAGE_W["letter"], cw + 2 * margin)
        ph = max(_PAGE_H, ch + 2 * margin)
        ox, oy = (pw - cw) / 2 - minx, (ph - ch) / 2 - miny
    elif page_mode in ("single", "double"):
        pw = max(_PAGE_W[page_mode], cw + 2 * margin)
        ph = ch + 2 * margin
        ox, oy = (pw - cw) / 2 - minx, margin - miny
    elif page_mode == "auto":
        want = cw + 2 * margin
        pw = 240.0 if want <= 240.0 else (504.0 if want <= 504.0 else want)
        ph = ch + 2 * margin
        ox, oy = (pw - cw) / 2 - minx, margin - miny
    else:  # 'fit'
        pw, ph = cw + 2 * margin, ch + 2 * margin
        ox, oy = margin - minx, margin - miny
    # The three FIXED widths are a promise the caller made to a journal, so having
    # to break one is news. 'auto' is excluded on purpose -- it asks for "the
    # smallest standard column that fits, else fit-to-content", so going wide is
    # the answer rather than a failure -- and so is 'fit', which has no width to
    # exceed.
    #
    # HEIGHT TOO, on the same terms and for a worse case. 'single' and 'double'
    # fix a width and let the height run, which is right -- a column figure is as
    # tall as it needs to be -- but nothing said where "tall" stopped being
    # printable. Measured: an 8-product scope table at page='single' came back
    # 240 x 1114 pt, i.e. 3.3 x 15.5 INCHES, with `ok: true` and no warning. One
    # column was the honest answer to the arithmetic (each cell needs 127 pt and
    # the column has 180); a figure a page and a half long was not something to
    # hand over in silence.
    if page_mode in _PAGE_W:
        if pw - _PAGE_W[page_mode] >= _PAGE_OVER_EPS:
            PAGE_OVERFLOW.append(("width", page_mode, _PAGE_W[page_mode], pw))
        if ph - _PAGE_H >= _PAGE_OVER_EPS:
            PAGE_OVERFLOW.append(("height", page_mode, _PAGE_H, ph))
    if shift is not None:
        shift.append((ox, oy))
    return page(header_prefix, _translate_body(body, ox, oy), [0, 0, round(pw), round(ph)])


# ----------------------------------------------------------------------
# Elliptical-arc arrows (curved / cycle / oval).
#
# Geometry reverse-engineered from a hand-drawn ChemDraw export and verified
# against 5 independent arcs to < 1 degree:
#
#     P(t) = C + cos(t) * Major + sin(t) * Minor
#            Major = MajorAxisEnd3D - Center3D      (ABSOLUTE coords in the file)
#            Minor = MinorAxisEnd3D - Center3D
#     the arc runs t = theta0 .. theta0 + sweep,  and  AngularSize = -sweep
#
# The sign flip is because CDXML y grows DOWN: a POSITIVE `sweep` here is
# visually counter-clockwise. Presence of AngularSize is what makes an <arrow>
# an arc; a straight arrow simply omits it.
#
# Do NOT read an arc's endpoints back from Head3D/Tail3D -- ChemDraw leaves them
# stale on a flipped/duplicated arc. Recompute them from C/Major/Minor/AngularSize.
# ----------------------------------------------------------------------

# Arrowhead size triples (HeadSize, ArrowheadCenterSize, ArrowheadWidth), measured
# from a hand-drawn ChemDraw export -- one per arrow TYPE, not one per palette size
# slot, so this is not yet a calibrated small/medium/large ladder.
#
# The triple is not three free numbers though: across all three measured cases
# (HeadSize 800, 1000, 1500) the other two are locked to it,
#     ArrowheadCenterSize = 0.875 * HeadSize      (1312.5 rounds to 1313)
#     ArrowheadWidth      = 0.250 * HeadSize
# so this table could collapse to three HeadSize values plus those two ratios.
# Left as explicit triples for now: they are what was actually observed.
HEAD_PRESETS = {
    "straight": (1000, 875, 250),    # plain straight arrow, ChemDraw default
    "arc":      (800, 700, 200),     # curved arrow, ChemDraw default (one size down)
    "equil":    (1500, 1313, 375),   # equilibrium arrow, ChemDraw default (one size up)
}
ARC_HEAD_SIZE, ARC_HEAD_CENTER, ARC_HEAD_WIDTH = HEAD_PRESETS["arc"]


def _arc_point(cx, cy, rx, ry, t_deg):
    t = math.radians(t_deg)
    return cx + math.cos(t) * rx, cy + math.sin(t) * ry


def _arc_bbox(cx, cy, rx, ry, t0, sweep, pad):
    """True extent of the DRAWN arc (not the whole ellipse), sampled + padded."""
    n = max(8, int(abs(sweep) / 3) + 1)
    pts = [_arc_point(cx, cy, rx, ry, t0 + sweep * i / n) for i in range(n + 1)]
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    return min(xs) - pad, min(ys) - pad, max(xs) + pad, max(ys) + pad


# THE ARROW SHAFT STYLE, in one place.
#
# `arrows.apply_style` and `arc_arrow` below both write `LineType` on an <arrow>,
# and both read this one dict rather than each holding a copy. Two hand
# transcriptions of one mapping is a drift hazard with nothing to catch it: adding
# a value to one would leave the other quietly behind, and the two paths are far
# enough apart (a straight shaft's attributes vs an arc's) that nobody would think
# to look.
#
# HERE RATHER THAN IN arrows.py because the dependency runs that way -- arrows
# imports this module, so this is the layer both can see. `arrows.line_type` is an
# alias, the same shape the key sets use with preflight.
#
# The `.get(..., "Dashed")` fallback is KEPT, and it is the second line of
# defence rather than the policy: `preflight.check_line_type` refuses an
# unrecognised value before any of this runs. On its own the fallback would turn
# `line="dotted"` into a dashed line in silence.
#
# `wavy` is still mapped. ChemDraw stores LineType="Wavy" and draws it like a plain
# shaft (measured; see arrows.py), so the tool layer refuses the value -- but this
# function is also called directly by scripts and test files, and dropping the
# mapping would change what those write. Refusing at the edge and mapping
# faithfully here are not in conflict.
_LINE_TYPES = {"dashed": "Dashed", "bold": "Bold", "wavy": "Wavy"}


def line_type(line: str) -> str:
    """The `LineType` word ChemDraw writes for one of the shaft-style names."""
    return _LINE_TYPES.get(line, "Dashed")


def arc_arrow(aid: int, cx: float, cy: float, rx: float, ry: float,
              t0: float, sweep: float, *, head: str | None = "Full",
              tail: str | None = None, head_type: str = "Solid",
              head_size: int = ARC_HEAD_SIZE, head_width: int = ARC_HEAD_WIDTH,
              head_center: int = ARC_HEAD_CENTER, line: str | None = None) -> str:
    """Elliptical-arc arrow. t0/sweep in degrees; +sweep is visually CCW.
    head/tail: 'Full' | 'HalfLeft' | 'HalfRight' | None. sweep=+/-359 -> closed cycle."""
    x0, y0 = _arc_point(cx, cy, rx, ry, t0)
    x1, y1 = _arc_point(cx, cy, rx, ry, t0 + sweep)
    pad = (head_width / 100.0) if (head or tail) else 0.6
    bx0, by0, bx1, by1 = _arc_bbox(cx, cy, rx, ry, t0, sweep, pad)
    a = [f'<arrow id="{aid}"',
         f'BoundingBox="{bx0:.2f} {by0:.2f} {bx1:.2f} {by1:.2f}"',
         'FillType="None"',
         f'ArrowheadType="{head_type}"']
    if head:
        a.append(f'ArrowheadHead="{head}"')
    if tail:
        a.append(f'ArrowheadTail="{tail}"')
    a += [f'HeadSize="{head_size}"',
          f'ArrowheadCenterSize="{head_center}"',
          f'ArrowheadWidth="{head_width}"',
          f'AngularSize="{-sweep:.6g}"',
          f'Center3D="{cx:.2f} {cy:.2f} 0"',
          f'MajorAxisEnd3D="{cx + rx:.2f} {cy:.2f} 0"',
          f'MinorAxisEnd3D="{cx:.2f} {cy + ry:.2f} 0"',
          f'Tail3D="{x0:.2f} {y0:.2f} 0"',
          f'Head3D="{x1:.2f} {y1:.2f} 0"']
    if line and line != "solid":
        a.append('LineType="%s"' % line_type(line))
    return " ".join(a) + "/>"


def arc_points(el, n: int = 180) -> list[tuple[float, float]]:
    """Decode an arc <arrow> back into a polyline. THE canonical reader for the geometry
    arc_arrow() writes -- every renderer must call this rather than re-deriving it, or the
    convention ends up implemented twice and drifts.

    `el` is an lxml Element, an attribute dict, OR the raw element STRING that arc_arrow
    returned -- the writer hands back strings, so the reader accepts them.
    Returns [] if the element is not an arc (no AngularSize), so a caller can branch on
    truthiness. Endpoints come from Center3D + Major/MinorAxisEnd3D + AngularSize, NOT
    from Head3D/Tail3D, which ChemDraw leaves stale on a mirrored or duplicated arc."""
    if isinstance(el, str):
        el = dict(re.findall(r'(\w+)="([^"]*)"', el))
    if not (el.get("AngularSize") and el.get("Center3D") and el.get("MajorAxisEnd3D")):
        return []

    def pt(name):
        a, b = el.get(name).split()[:2]
        return float(a), float(b)

    cx, cy = pt("Center3D")
    rx = pt("MajorAxisEnd3D")[0] - cx
    ry = pt("MinorAxisEnd3D")[1] - cy
    if not rx or not ry:
        return []
    sweep = -float(el.get("AngularSize"))
    tx, ty = pt("Tail3D")
    t0 = math.degrees(math.atan2((ty - cy) / ry, (tx - cx) / rx))
    return [(cx + math.cos(math.radians(t0 + sweep * i / n)) * rx,
             cy + math.sin(math.radians(t0 + sweep * i / n)) * ry)
            for i in range(n + 1)]


def arc_chord(aid: int, x0: float, y0: float, x1: float, y1: float,
              bulge: float = 0.3, **kw) -> str:
    """Circular arc from (x0,y0) to (x1,y1) bowing out by `bulge` * chord length.
    bulge > 0 bows to the LEFT of travel (visually up for a left-to-right chord);
    negative bows the other way. This is the form a layout wants: name the two
    points the arc has to touch and let the radius fall out."""
    ccx, ccy, r, t0, sweep = arc_circle(x0, y0, x1, y1, bulge)
    return arc_arrow(aid, ccx, ccy, r, r, t0, sweep, **kw)


def arc_circle(x0: float, y0: float, x1: float, y1: float, bulge: float):
    """(cx, cy, r, t0, sweep) for the circular arc arc_chord would draw.

    Extracted from arc_chord unchanged, because a second caller needs the
    arc as GEOMETRY rather than as XML: to stop an electron-pushing arrow short
    of the atom it points at, you have to trim it by arc length, and that needs
    the radius. One copy of the derivation, two callers.
    """
    dx, dy = x1 - x0, y1 - y0
    L = math.hypot(dx, dy)
    if L == 0 or bulge == 0:
        raise ValueError("arc_chord needs a non-zero chord and bulge")
    h = abs(bulge) * L
    r = (L * L / 4.0 + h * h) / (2.0 * h)          # circle through both ends, sagitta h
    d = r - h                                       # midpoint -> centre
    mx, my = (x0 + x1) / 2.0, (y0 + y1) / 2.0
    sgn = 1.0 if bulge > 0 else -1.0
    # unit normal; +normal is the side the arc bows toward, so the centre is -normal
    nx, ny = (dy / L) * sgn, (-dx / L) * sgn
    ccx, ccy = mx - nx * d, my - ny * d
    t0 = math.degrees(math.atan2(y0 - ccy, x0 - ccx))
    t1 = math.degrees(math.atan2(y1 - ccy, x1 - ccx))
    sweep = (t1 - t0) % 360.0
    if sgn < 0:                                     # pick the half that bows the right way
        sweep -= 360.0
    return ccx, ccy, r, t0, sweep


def cycle_arrow(aid: int, cx: float, cy: float, r: float, *,
                t0: float = 180.0, gap: float = 0.4, heads: int = 1,
                cw: bool = True, **kw) -> list[str]:
    """Closed circular cycle arrow (the 'catalyst turnover' ring). ALWAYS returns a
    list, since the head count decides the element count.

    heads=1 -> ONE near-360 arc. Same thing as taking the palette's 270-degree arc
      and stretching it to 359.x by hand.
    heads=2 -> TWO half arcs, each with its own head: the ring otherwise built by
      hand as "draw 180 degrees, copy, rotate/mirror".

    `gap` degrees are left open per arc so each head reads as a head and not as a
    join. ChemDraw's own exports sit at 359.3-359.6, i.e. a gap of 0.4-0.7."""
    sgn = -1.0 if cw else 1.0
    if heads == 2:
        span = 180.0 - gap
        return [arc_arrow(aid, cx, cy, r, r, t0, span * sgn, **kw),
                arc_arrow(aid + 1, cx, cy, r, r, t0 + 180.0, span * sgn, **kw)]
    return [arc_arrow(aid, cx, cy, r, r, t0, (360.0 - gap) * sgn, **kw)]


def oval(aid: int, cx: float, cy: float, rx: float, ry: float, **kw) -> str:
    """Headless closed ellipse, drawn as a near-360 arc with no arrowhead — the
    outline used to ring a catalyst / condition label."""
    kw.setdefault("head", None)
    kw.setdefault("tail", None)
    return arc_arrow(aid, cx, cy, rx, ry, 0.0, -359.6, **kw)


# ----------------------------------------------------------------------
# BEZIER pushing arrows -- the shape a mechanism actually draws.
#
# An arc has ONE degree of freedom (the bulge); its tangents at both ends fall
# out of the two endpoints and that one number. A mechanism arrow needs FOUR --
# leave along the lone pair, arrive anti to the leaving group -- and 25% of the
# 1460 curved arrows measured in a reference set of published mechanism figures
# carry an INFLECTION, which no conic arc has at all. So the arc cannot be
# stretched into this; it is a different curve family.
#
# ENCODING, read off a hand-drawn two-step pen-tool file and verified by
# regenerating its CurvePoints byte for byte from six numbers:
#
#     <curve CurveType="8" ArrowheadHead="Full" ArrowheadType="Solid"
#            CurvePoints="in knot out  in knot out  ..."/>
#
#   * CurvePoints is 3N point PAIRS = N knots, each written as its incoming
#     handle, itself, its outgoing handle. A smooth knot sits exactly midway
#     between its two handles, which is how the format was identified.
#   * the drawn path is the chain of cubics knot_i -> knot_{i+1} with controls
#     out_i and in_{i+1}.
#   * `CurveType="8"` is what ChemDraw draws the LINE for.
#
# AND THIS IS WHY A curve()-BUILT ARROW DRAWS A HEAD WITH NO LINE.
# curve() writes 3n+1 control points -- 4 for one cubic -- which is not 3N: four
# points read as one triplet plus a leftover, i.e. ONE knot, i.e. ZERO segments.
# Zero length of line, and the arrowhead still gets stamped at that knot. The
# problem is not the attribute; it is the point count.
# ----------------------------------------------------------------------


# LONE PAIRS AND RADICAL DOTS. Bound to their atom the way a charge symbol is:
# a <graphic GraphicType="Symbol"> carrying <represent attribute="Radical">, so
# ChemDraw treats it as that atom's property rather than as loose ink.
#
# GEOMETRY, off one drawn by hand and judged professional-looking (an oxygen in
# a test file) and cross-checked against 37 lone pairs in real ChemDraw files:
#
#   * the BoundingBox is a SEGMENT, and the two dots of the pair lie along it.
#   * it runs PERPENDICULAR to the atom -> pair direction. Median deviation
#     from perpendicular in the real files: 0.1 degree, over all 37.
#   * the offset is to the segment's CENTRE, not to its first coordinate. This
#     is the one place the charge-symbol rule does NOT hold, and getting it wrong
#     makes the pair look wrong rather than absent: it sits off to one side of
#     its own direction instead of centred on it.
#
#           hand-drawn       real-file median (n=37)
#   centre    0.547 bond          0.473 bond
#   length    0.260 bond          0.208 bond
#
# The hand-drawn value is the pinned one -- that file is the source of truth and
# the real-file median is within 15% of it, so there is nothing to reconcile.
#
# NO IgnoreWarnings. The hand-drawn file carries none; setting it is useful only
# for spotting a duplicate radical dot, and it is not part of the convention.
#
# DIRECTION IS NOT THE FREE VALENCE. Measured against it the real files deviate by a
# median of 78 degrees, because an atom with two pairs puts them on two sides and
# a single bisector cannot say that. The caller supplies angles; decorate derives
# them from the GAPS between the atom's bonds.
LONE_PAIR_OFFSET = 0.547 * BOND_LENGTH
LONE_PAIR_LENGTH = 0.260 * BOND_LENGTH
RADICAL_OFFSET = 0.43 * BOND_LENGTH
_ELECTRON_SYMBOL = {"lonepair": "LonePair", "radical": "Electron"}


def electron_symbols(gid: int, x: float, y: float, angles, node_id: int,
                     kind: str = "lonepair") -> list[str]:
    """One <graphic> per angle, bound to `node_id`. `angles` in radians."""
    sym = _ELECTRON_SYMBOL.get(kind)
    if sym is None:
        raise ValueError("electron kind %r; one of %s"
                         % (kind, ", ".join(sorted(_ELECTRON_SYMBOL))))
    off = LONE_PAIR_OFFSET if kind == "lonepair" else RADICAL_OFFSET
    half = (LONE_PAIR_LENGTH if kind == "lonepair" else 0.0) / 2.0
    out = []
    for k, a in enumerate(angles):
        cx = x + math.cos(a) * off
        cy = y + math.sin(a) * off
        # the segment runs across the direction, so the two dots straddle it
        nx, ny = -math.sin(a) * half, math.cos(a) * half
        out.append('<graphic id="%d" BoundingBox="%.2f %.2f %.2f %.2f"'
                   ' GraphicType="Symbol" SymbolType="%s">'
                   '<represent attribute="Radical" object="%d"/></graphic>'
                   % (gid + k, cx - nx, cy - ny, cx + nx, cy + ny, sym, node_id))
    return out


def bezier_arrow(cid: int, p0, dir0, len0: float, p1, dir1, len1: float, *,
                 head: str | None = "Full", tail: str | None = None,
                 head_type: str = "Solid", color=None,
                 line: str | None = None) -> str:
    """One cubic, named the way the pen tool names it.

    `dir0` is the direction the arrow LEAVES p0 travelling; `dir1` the direction
    it is travelling when it ARRIVES at p1. Both are (dx, dy) and need not be
    unit -- they are normalised here, because every caller has them as a
    difference of two atom positions and normalising at each call site is how a
    sign error gets in. `len0` / `len1` are the handle lengths in drawing units.

    Both handles are written symmetrically about their knot, which is what makes
    the knot smooth and is what ChemDraw's own two-step tool produces.
    """
    d0, d1 = _unit(dir0), _unit(dir1)
    h0 = (d0[0] * len0, d0[1] * len0)
    h1 = (d1[0] * len1, d1[1] * len1)
    pts = [(p0[0] - h0[0], p0[1] - h0[1]), (p0[0], p0[1]),
           (p0[0] + h0[0], p0[1] + h0[1]),
           (p1[0] - h1[0], p1[1] - h1[1]), (p1[0], p1[1]),
           (p1[0] + h1[0], p1[1] + h1[1])]
    drawn = bezier_points_from(pts)
    xs = [q[0] for q in drawn]
    ys = [q[1] for q in drawn]
    a = ['<curve id="%d"' % cid]
    if color is not None:
        a.append('color="%d"' % color_index(color))
    a.append('CurveType="8"')
    a.append('BoundingBox="%.2f %.2f %.2f %.2f"'
             % (min(xs), min(ys), max(xs), max(ys)))
    if head:
        a.append('ArrowheadHead="%s"' % head)
    if tail:
        a.append('ArrowheadTail="%s"' % tail)
    a.append('ArrowheadType="%s"' % head_type)
    if line and line != "solid":
        a.append('LineType="%s"' % line_type(line))
    a.append('CurvePoints="%s"' % " ".join("%.2f %.2f" % q for q in pts))
    return " ".join(a) + "/>"


def bezier_points_from(pts, per: int = 24) -> list[tuple[float, float]]:
    """THE canonical reader for what bezier_arrow writes: triplets -> polyline.

    Same role arc_points() plays for arc_arrow(). Every renderer, extent and test
    must come through here rather than re-deriving the triplet convention, or it
    ends up implemented twice and the two drift. Takes the raw point list so it
    can be called before the XML exists.
    """
    knots = [pts[3 * k + 1] for k in range(len(pts) // 3)]
    out = [knots[0]] if knots else []
    for k in range(len(knots) - 1):
        p0, p3 = knots[k], knots[k + 1]
        c1, c2 = pts[3 * k + 2], pts[3 * (k + 1)]
        for i in range(1, per + 1):
            t = i / per
            u = 1.0 - t
            out.append((u * u * u * p0[0] + 3 * u * u * t * c1[0]
                        + 3 * u * t * t * c2[0] + t * t * t * p3[0],
                        u * u * u * p0[1] + 3 * u * u * t * c1[1]
                        + 3 * u * t * t * c2[1] + t * t * t * p3[1]))
    return out


def bezier_points(el, per: int = 24) -> list[tuple[float, float]]:
    """bezier_points_from() for a <curve> element or its XML string."""
    cp = el.get("CurvePoints") if hasattr(el, "get") else (
        re.search(r'CurvePoints="([^"]*)"', el) or [None, ""])[1]
    if hasattr(el, "get"):
        raw = (cp or "").split()
    else:
        m = re.search(r'CurvePoints="([^"]*)"', el)
        raw = (m.group(1) if m else "").split()
    pts = [(float(raw[i]), float(raw[i + 1])) for i in range(0, len(raw) - 1, 2)]
    if len(pts) % 3:
        raise ValueError(
            "CurvePoints holds (in-handle, knot, out-handle) triplets, so the "
            "count must be a multiple of 3; got %d -- four points is the old "
            "3n+1 cubic convention and ChemDraw draws it as a head with no line"
            % len(pts))
    return bezier_points_from(pts, per)


def _unit(v):
    L = math.hypot(v[0], v[1])
    return (v[0] / L, v[1] / L) if L else (1.0, 0.0)


def arc_tangent(aid: int, x0: float, y0: float, x1: float, y1: float, *,
                tangent: tuple[float, float] = (1.0, 0.0), at: str = "head",
                **kw) -> str:
    """Arc from (x0,y0) to (x1,y1) whose TANGENT AT THE HEAD is `tangent`.

    This is the form a reaction figure actually wants: a curved arrow that merges
    smoothly into the shaft instead of hitting it at an angle. arc_chord() only
    constrains "passes through both points with this sagitta", which leaves the
    landing angle to fall out however it likes -- visibly a kink where the curve
    meets a straight arrow.

    `tangent` is the direction of TRAVEL at the constrained end, so (1, 0) means
    "runs left-to-right there", i.e. tangent to a left-to-right shaft. The circle is
    then unique: its centre lies on the normal at that end, equidistant from both
    points.

    at="head" (default): smooth ARRIVAL - a reagent curving down INTO the flow.
    at="tail":           smooth DEPARTURE - a by-product peeling OFF the shaft.

    The other point lying on the tangent line is degenerate (infinite radius), and
    raises rather than emitting a near-straight arc pretending to be a curve."""
    if at not in ("head", "tail"):
        raise ValueError("at must be 'head' or 'tail'")
    tl = math.hypot(*tangent)
    if tl == 0:
        raise ValueError("tangent must be non-zero")
    tx, ty = tangent[0] / tl, tangent[1] / tl
    nx, ny = -ty, tx                       # normal at the constrained end
    # anchor = the constrained end, free = the other one
    ax_, ay_ = (x1, y1) if at == "head" else (x0, y0)
    fx_, fy_ = (x0, y0) if at == "head" else (x1, y1)
    dx, dy = fx_ - ax_, fy_ - ay_
    denom = 2.0 * (dx * nx + dy * ny)
    if abs(denom) < 1e-9:
        raise ValueError(f"the free point is on the tangent line through the {at} "
                         "(infinite radius) - use a straight arrow instead")
    s = (dx * dx + dy * dy) / denom
    ccx, ccy = ax_ + s * nx, ay_ + s * ny
    r = abs(s)
    t0 = math.degrees(math.atan2(y0 - ccy, x0 - ccx))
    t1 = math.degrees(math.atan2(y1 - ccy, x1 - ccx))
    # travel direction at angle T for a POSITIVE sweep is (-sin T, cos T); pick the
    # sweep sign that makes the arrow run ALONG `tangent` at the constrained end
    tr = math.radians(t1 if at == "head" else t0)
    forward = (-math.sin(tr)) * tx + math.cos(tr) * ty
    sweep = (t1 - t0) % 360.0
    if forward <= 0:
        sweep -= 360.0
    return arc_arrow(aid, ccx, ccy, r, r, t0, sweep, **kw)
