"""
depiction_library.py — canned CDXML depictions that OUTRANK the layout engine.

Some chemistry cannot be said in SMILES, so no amount of work on the converter will
ever draw it. Ferrocene is the type case: the iron is not bonded to ten carbons, it
sits BETWEEN two rings, and a chemist draws that as two tilted cyclopentadienyls,
a dashed bond from the metal to each ring's centre, and an arc inside each ring.
There is no SMILES for it and no rdkit object for it -- rdkit will happily accept
`[Fe+2].c1cc[cH-]c1.c1cc[cH-]c1` and lay it out as three unrelated pieces sitting
in a row, which is not wrong so much as not a drawing of ferrocene.

The answer is not to compute it. It is to keep the picture somebody already drew
correctly, and to recognise when it has been asked for. So this module holds a
folder of hand-drawn CDXML files and a table saying which molecule each one is;
`smiles_to_fragment` asks here FIRST and only falls through to rdkit on a miss.

ADDING AN ENTRY (the whole point -- no Python required beyond one line):
  1. Draw it in ChemDraw.
  2. File > Save As > CDXML into depictions/.
  3. Add a row to LIBRARY below: the file name, and every SMILES spelling that
     means this molecule. Spell them generously -- the point is to catch what a
     person would actually type, including spellings that are formally a different
     compound (see ferrocene's four).
Each file is a COMPLETE CDXML document, so it opens in ChemDraw on a double click:
a library entry can be inspected the same way it was made.

WHAT AN ENTRY MUST CONTAIN: exactly one top-level <fragment>. Everything the
drawing needs has to be inside it -- nodes, bonds, the <t> label of a metal, the
<curve> arcs -- because a fragment is the unit that gets moved onto the page as one
object. A <curve> or a MultiAttachment centroid outside the fragment would be left
behind at its original coordinates.

MATCHING, two ways.

  WHOLE MOLECULE, by canonical SMILES (`lookup`). The canned file is handed over
  as the drawing, untouched.

  AS A CORE (`core_lookup`), for a molecule that merely CONTAINS one. The canned
  coordinates become a depiction template: they pin the core and rdkit builds the
  rest around them, so the drawing that comes out is the caller's molecule with
  somebody else's picture inside it. This is what makes the sugars useful --
  whole-molecule matching reached NOTHING in a real carbohydrate figure, because
  those figures draw glycosides and C-glycosides, never a free sugar.

  A substituted FERROCENE is still only recognised and reported, not drawn from
  the entry: a core has to be a thing rdkit can match, and the sandwich's
  MultiAttachment centroids are not atoms it knows about.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

import cdxml_core as C

try:
    from rdkit import Chem
    from rdkit import RDLogger
    from rdkit.Chem import rdmolops
    RDLogger.DisableLog("rdApp.*")      # a spelling that does not parse is data, not a crash
except ImportError as e:  # pragma: no cover - only meaningful on the Mac
    raise SystemExit("rdkit not found; depiction_library needs it to match SMILES") from e


_DIR = Path(__file__).resolve().parent / "depictions"


# --- the library -------------------------------------------------------------
# (file in depictions/, name for the report, SMILES spellings that mean it)
LIBRARY: list[tuple[str, str, tuple[str, ...]]] = [
    (
        "ferrocene.cdxml", "ferrocene",
        (
            # the standard reading: Fe(II) between two cyclopentadienyl anions
            "[Fe+2].c1cc[cH-]c1.c1cc[cH-]c1",
            # neutral Fe written with the same two anions -- charge-unbalanced, and
            # very common, because people write the metal and the ligands separately
            "[Fe].c1cc[cH-]c1.c1cc[cH-]c1",
            # Fe with two neutral cyclopentaDIENES: two hydrogens too many, and what
            # you get from anyone who types the rings without thinking about charge
            "[Fe].C1=CCC=C1.C1=CCC=C1",
            # sigma-bonded, i.e. someone forcing the connectivity SMILES does have
            "C1=C[CH]([Fe][CH]2C=CC=C2)C=C1",
        ),
    ),
]


# --- what happened this batch (same contract as smiles_to_frag's reports) -----
# A canned depiction is a SUBSTITUTION: the drawing that comes out is not the one
# the SMILES describes, it is the one a chemist would have drawn instead. Right --
# and not something to do silently.
CANNED_USED: list[tuple[str, str]] = []          # (smiles as asked, entry name)
CANNED_CORE_ONLY: list[tuple[str, str]] = []     # (smiles, entry name) -- substituted
# The entry's coordinates were used for PART of the molecule (see core_lookup).
# Worth its own list rather than folding into CANNED_USED: the drawing that comes
# out is still the caller's molecule, only its core is somebody else's picture.
CANNED_CORE_USED: list[tuple[str, str]] = []     # (smiles, entry name)


def reset_reports() -> None:
    CANNED_CORE_USED.clear()
    CANNED_USED.clear()
    CANNED_CORE_ONLY.clear()


class Entry:
    """One canned depiction: its fragment text, its size, and what it matches."""

    def __init__(self, filename: str, name: str, spellings: tuple[str, ...]):
        self.filename = filename
        self.name = name
        self.spellings = spellings
        self._frag: str | None = None
        self._bbox: list[float] | None = None
        self._core = False              # False = not worked out yet, None = no core
        self.canonical: set[str] = set()
        self.query: "Chem.Mol | None" = None
        for smi in spellings:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                raise ValueError("depiction_library: %s spelling does not parse: %r"
                                 % (name, smi))
            self.canonical.add(Chem.MolToSmiles(mol))
            if self.query is None:
                self.query = mol            # first spelling is the substructure query

    # -- the file ------------------------------------------------------------
    def _load(self) -> None:
        path = _DIR / self.filename
        text = path.read_text()
        root = ET.fromstring(text)
        page = root.find("page")
        frags = page.findall("fragment") if page is not None else []
        if len(frags) != 1:
            raise ValueError("%s: a library entry needs exactly one top-level "
                             "<fragment>, found %d" % (path.name, len(frags)))
        # Slice the ORIGINAL TEXT rather than re-serialising the parse: every byte
        # of a hand drawing is deliberate, and a round trip through ElementTree
        # rewrites attribute order and quoting for no reason.
        start = text.index("<fragment")
        end = text.rindex("</fragment>") + len("</fragment>")
        frag = text[start:end]
        # ids down to 1..n, so the fragment behaves like a cdxml_core.fragment() and
        # cannot collide with a neighbour once layout adds its per-structure base.
        frag = _compact_ids(frag)
        bbox = _ink_bbox(frags[0])
        # ...and to (0, 0), which is the contract smiles_to_fragment() returns under
        frag = C.translate(frag, -bbox[0], -bbox[1])
        self._frag = frag
        self._bbox = [0.0, 0.0, bbox[2] - bbox[0], bbox[3] - bbox[1]]

    def fragment(self) -> tuple[str, list[float]]:
        """(fragment_str, bbox) -- exactly what smiles_to_fragment() returns."""
        if self._frag is None:
            self._load()
        return self._frag, list(self._bbox)

    # -- the entry as a CORE ---------------------------------------------------
    def core(self):
        """(reference mol carrying the drawing, match pattern), or None.

        The reference is read out of the entry's own CDXML by rdkit, so its
        conformer IS the drawing -- the chair, not a fresh layout of the same
        molecule. The pattern is that reference with the anomeric substituent
        made a wildcard.

        None for an entry with no anomeric-type position, which is every entry
        that is not a sugar. Those keep whole-molecule matching and nothing else.
        """
        if self._core is not False:
            return self._core
        self._core = None
        mols = Chem.MolsFromCDXMLFile(str(_DIR / self.filename))
        if len(mols) == 1:
            ref = _rescaled(mols[0])
            if ref is not None:
                i = _anomeric_exocyclic(ref)
                if i is not None:
                    rw = Chem.RWMol(ref)
                    rw.GetAtomWithIdx(i).SetAtomicNum(0)
                    q = rw.GetMol()
                    Chem.SanitizeMol(q)
                    # A dummy built this way is a CONCRETE atom-zero unless the
                    # query is adjusted: matching a molecule against a molecule,
                    # rdkit asks for atomic number 0 on the other side too, and
                    # every real target missed. makeDummiesQueries is what turns
                    # it into "anything".
                    p = Chem.AdjustQueryParameters.NoAdjustments()
                    p.makeDummiesQueries = True
                    self._core = (ref, rdmolops.AdjustQueryProperties(q, p))
        return self._core


def _compact_ids(frag: str) -> str:
    """Renumber every id in a fragment to 1..n, references included."""
    import re
    seen: dict[int, int] = {}
    for m in re.finditer(r'\bid="(\d+)"', frag):
        v = int(m.group(1))
        if v not in seen:
            seen[v] = len(seen) + 1

    def one(attr, m):
        v = int(m.group(1))
        return '%s="%d"' % (attr, seen.get(v, v))

    for attr in ("id", "B", "E", "BeginAttach", "EndAttach"):
        frag = re.sub(r'\b' + attr + r'="(\d+)"', lambda m, a=attr: one(a, m), frag)
    for attr in C._ID_LIST_ATTRS:
        frag = re.sub(
            r'\b' + attr + r'="([\d\s]+)"',
            lambda m, a=attr: '%s="%s"' % (a, " ".join(
                str(seen.get(int(t), int(t))) if int(t) else "0"
                for t in m.group(1).split())),
            frag)
    return frag


def _ink_bbox(frag_el) -> tuple[float, float, float, float]:
    """The drawing's real extent: atom points, text boxes and curve points.

    Not the file's own BoundingBox attribute -- ChemDraw pads that, and a canned
    entry has to measure the same way cdxml_core.fragment() does or a scope-table
    cell built from one would reserve a different amount of room than a cell built
    from a converted SMILES.
    """
    xs: list[float] = []
    ys: list[float] = []
    for n in frag_el.iter("n"):
        if n.get("p"):
            x, y = (float(v) for v in n.get("p").split())
            xs.append(x)
            ys.append(y)
    for t in frag_el.iter("t"):
        if t.get("BoundingBox"):
            a, b, c, d = (float(v) for v in t.get("BoundingBox").split())
            xs += [a, c]
            ys += [b, d]
    for cv in frag_el.iter("curve"):
        if cv.get("CurvePoints"):
            v = [float(t) for t in cv.get("CurvePoints").split()]
            xs += v[0::2]
            ys += v[1::2]
    if not xs:
        raise ValueError("library entry has no coordinates")
    return min(xs), min(ys), max(xs), max(ys)


# --- using an entry as a CORE, not only as a whole molecule -------------------
# Twenty-eight canned sugars reached NOTHING in a real carbohydrate figure. Run
# over a published figure on C-glycoside anomeric coupling, the library was
# asked 24 times about 18 distinct structures, 16 of them with a pyranose ring,
# and matched zero -- not even "core recognised". Every sugar in that figure is
# a glucopyranose whose ANOMERIC substituent is something other than OH:
#
#   canned beta-D-glucopyranose  OC[C@H]1O[C@@H](O            )[C@H](O)[C@@H](O)[C@@H]1O
#   the figure's glycosyl sulfone OC[C@H]1O[C@@H](S(=O)(=O)[*:1])[C@H](O)[C@@H](O)[C@@H]1O
#
# Ring, ring oxygen and the other four hydroxyls are identical character for
# character. In carbohydrate chemistry that is not an edge case, it is the field.
#
# So an entry also offers a CORE: itself, with the anomeric substituent free to
# be anything. The anomeric position is defined structurally -- the terminal
# heteroatom on a ring atom that also carries a ring heteroatom -- so it finds
# itself on any pyranose or furanose and finds nothing on ferrocene, which keeps
# every non-sugar entry exactly as it was.
#
# WHAT IS NOT LOOSENED: everything else stays a literal match, stereochemistry
# included. Checked on four targets -- a beta-glycoside matches only the beta
# entry, an alpha-glycoside only the alpha, a galactoside only galactose, and a
# fully O-methylated (protected) glucose still matches glucose, because
# substituting an OH's hydrogen leaves the oxygen in place.

_RDKIT_BOND = 1.5          # what rdkit's own 2D layout uses; the reference has
                           # to be on the same ruler or the atoms it does not
                           # pin come out the wrong size


def _rescaled(mol):
    """The entry's own drawing, on rdkit's bond-length ruler."""
    import math
    import statistics
    conf = mol.GetConformer()

    def xy(i):
        p = conf.GetAtomPosition(i)
        return (p.x, p.y)

    lens = [math.dist(xy(b.GetBeginAtomIdx()), xy(b.GetEndAtomIdx()))
            for b in mol.GetBonds()]
    if not lens or statistics.median(lens) <= 0:
        return None
    k = _RDKIT_BOND / statistics.median(lens)
    for i in range(mol.GetNumAtoms()):
        x, y = xy(i)
        conf.SetAtomPosition(i, (x * k, y * k, 0.0))
    return mol


def _anomeric_exocyclic(mol):
    """Index of the terminal heteroatom at an anomeric-type centre, or None."""
    ri = mol.GetRingInfo()
    for a in mol.GetAtoms():
        if a.GetSymbol() != "C" or not ri.NumAtomRings(a.GetIdx()):
            continue
        in_ring, outside = None, None
        for n in a.GetNeighbors():
            if n.GetAtomicNum() in (1, 6):
                continue
            if ri.AreAtomsInSameRing(a.GetIdx(), n.GetIdx()):
                in_ring = n
            elif n.GetDegree() == 1:
                outside = n
        if in_ring is not None and outside is not None:
            return outside.GetIdx()
    return None


_ENTRIES = [Entry(f, n, s) for f, n, s in LIBRARY]


def lookup(smiles: str):
    """The canned entry for this SMILES, or None.

    A hit means the WHOLE molecule is that entry. A molecule that merely contains
    one is recorded in CANNED_CORE_ONLY and returns None, so the caller draws it
    the usual way and the tool can say what was recognised and not used.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    canon = Chem.MolToSmiles(mol)
    for e in _ENTRIES:
        if canon in e.canonical:
            CANNED_USED.append((smiles, e.name))
            return e
    for e in _ENTRIES:
        # useChirality, BECAUSE THE LIBRARY HOLDS 28 STEREOISOMERS OF EACH OTHER.
        # Without it the first skeleton match wins and this note names whichever
        # sugar happens to sit first in LIBRARY. Measured: 6-O-acetyl-D-glucose and
        # maltose both came back "alpha-D-allopyranose", and L-glucose -- whose
        # mirror image the library does not have at all -- came back with a hit.
        # With it on the same three give beta-D-glucopyranose, alpha-D-gluco-
        # pyranose and nothing. Harmless for a library of one achiral entry
        # (ferrocene), wrong as soon as it holds stereoisomers.
        if e.query is not None and mol.GetSubstructMatch(e.query, useChirality=True):
            CANNED_CORE_ONLY.append((smiles, e.name))
            return None
    return None


def core_lookup(smiles: str):
    """(entry, reference mol with the drawing's coordinates, match pattern).

    For a molecule that CONTAINS a canned entry's core. Returns None on a miss,
    and on a hit the caller is expected to hand the reference to rdkit's
    GenerateDepictionMatching2DStructure so the core is drawn from the canned
    coordinates and the rest is built around it.

    COORDGEN MUST BE OFF for that call. Measured on dapagliflozin: with
    rdDepictor.SetPreferCoordGen(True) the chair comes out FLAT -- ring interior
    angles spread 0.3 degrees against the reference's 105.4 -- because flattening
    a ring to its own idea of one is exactly CoordGen's job. With it off the
    spread is 105.4 to the decimal and bonds land in 1.49-1.55.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    for e in _ENTRIES:
        hit = e.core()
        if hit is None:
            continue
        ref, patt = hit
        if mol.HasSubstructMatch(patt, useChirality=True):
            CANNED_CORE_USED.append((smiles, e.name))
            return e, ref, patt
    return None


def names() -> list[str]:
    return [e.name for e in _ENTRIES]


if __name__ == "__main__":
    for e in _ENTRIES:
        frag, bb = e.fragment()
        print("%-12s %-22s %5.1f x %5.1f units   %d nodes, %d bonds, %d curves"
              % (e.name, e.filename, bb[2], bb[3],
                 frag.count("<n"), frag.count("<b"), frag.count("<curve")))
        for s in sorted(e.canonical):
            print("               matches", s)
