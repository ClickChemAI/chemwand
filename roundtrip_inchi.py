"""roundtrip_inchi.py — read the file back and check the molecules survived.

    from roundtrip_inchi import check
    check("/tmp/scheme.cdxml", {"1": "c1ccccc1", "2": "CCO"})
    -> {"checked": 2, "matched": 2, "failed": []}

WHAT THIS VERIFIES, and it is narrower than it sounds: that the connectivity
written is the connectivity given. It says nothing about whether the
input was right, and nothing about whether the figure is readable -- a scheme
with every structure correct and every label on top of the wrong arrow passes
this cleanly.

WHY RDKit'S PARSER RATHER THAN A HOME-GROWN ONE. Writing a second CDXML reader
would mean the check shares its assumptions with the writer, which is how a round
trip comes back clean because both ends are wrong the same way. RDKit's reader
was written by people who do not know what this engine's writer does.

MATCHED BY MULTISET, NOT BY POSITION. The parser returns molecules in document
order and drops the ones it cannot build, so position i in the output is not
cell i in the input -- comparing pairwise would report the wrong cell as failed
as soon as one earlier fragment was skipped. Comparing InChI multisets answers
the question actually being asked: did each molecule sent end up in the file.

KNOWN PARSER LIMITS, none of them the engine's fault and all of them counted as
`unparsed` rather than as mismatches:

  * a contracted superatom whose nested fragment RDKit cannot zip raises
    "molzip: bond info not properly setup" and that fragment is dropped from
    the returned list, silently. A three-structure scheme with such a
    superatom comes back with two.
  * drawing objects -- arrows, captions, brackets -- are not molecules and are
    correctly absent.
"""

from __future__ import annotations

import contextlib
import io

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def _inchi(mol):
    try:
        s = Chem.MolToInchi(mol)
        return s or None
    except Exception:
        return None


def _key(mol):
    """InChI, or a canonical-SMILES key for the structures InChI cannot hold.

    A Markush structure carries a dummy atom and InChI answers "" for those --
    the SAME empty string an unparseable molecule gives, so on InChI alone the
    two would be indistinguishable and every generic structure filed under
    `input_unparsed` and never compared. Scope tables are generic by definition,
    so the check would be off for exactly the figure type it matters most in.

    Nothing about the drawing is the problem: RDKit's CDXML reader brings
    these back complete -- dummy atoms, map numbers and wedges -- so a canonical
    SMILES compares them as exactly as an InChI would. Only the InChI FORMAT
    cannot represent one.

    MAP NUMBERS ARE STRIPPED, on both sides. Renaming a placeholder through
    rlabels ("R1" -> "Ar") loses the number on the way into the file, so keeping
    it would fail every renamed structure for a difference in a label. What stays
    compared is which ATOM carries a placeholder, which is the connectivity
    question this file exists to ask; R1 and R2 swapping names does not show up.
    """
    if mol is None:
        return None
    if any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        m = Chem.Mol(mol)
        for a in m.GetAtoms():
            a.SetAtomMapNum(0)
        try:
            return "SMILES:" + Chem.MolToSmiles(m)
        except Exception:
            return None
    return _inchi(mol)


def _key_of_smiles(smi):
    m = Chem.MolFromSmiles(smi)
    return _key(m) if m is not None else None


def read_mols(path):
    """Molecules RDKit can build from the file. Never raises."""
    err = io.StringIO()
    try:
        with contextlib.redirect_stderr(err):
            return list(Chem.MolsFromCDXMLFile(path) or [])
    except Exception:
        return []


def check(path, built_from: dict, skip: dict | None = None) -> dict:
    """`built_from` is {cell id: input SMILES}. Returns the report fragment.

    "THE FILE WAS NOT THERE" IS NOT "THE MOLECULES ARE WRONG". Without the guard
    below, an unreadable path makes `read_mols` return [] and every structure
    comes out in `failed` -- indistinguishable from a real mismatch. A failure
    goes onto `warnings`, where it would accuse a correct figure.
    """
    import os
    if path and not os.path.exists(path):
        return {"checked": 0, "matched": 0, "failed": [],
                "unreadable": "no file at %s -- nothing was checked" % path}
    want = {}
    unparsed_input = []
    skip = dict(skip or {})
    for cid, smi in (built_from or {}).items():
        if cid in skip:
            continue
        key = _key_of_smiles(smi)
        if key is None:
            unparsed_input.append(cid)
            continue
        want.setdefault(key, []).append(cid)

    got: dict = {}
    for m in read_mols(path):
        k = _key(m)
        if k:
            got[k] = got.get(k, 0) + 1

    failed = []
    matched = 0
    for key, cids in want.items():
        have = got.get(key, 0)
        matched += min(have, len(cids))
        for cid in cids[have:]:
            failed.append(cid)

    rep = {"checked": sum(len(v) for v in want.values()),
           "matched": matched,
           "failed": sorted(failed)}
    if unparsed_input:
        # the INPUT SMILES did not parse -- nothing to compare against, and not
        # a round-trip failure
        rep["input_unparsed"] = sorted(unparsed_input)
    if skip:
        # NOT a pass and not a failure. Named, because "checked: 6" over eight
        # cells is the shape this whole file exists to avoid.
        rep["not_checked"] = {c: skip[c] for c in sorted(skip)}
    return rep
