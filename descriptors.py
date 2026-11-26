"""descriptors.py — the numbers a paper reports about a molecule.

Pure RDKit; no ChemDraw, no network, no API. One call is microseconds, which is
why this is worth having at all: every one of these is already in the SI of every
med-chem paper, and a chemist who has the structure should not have to leave for
them.

TWO NAMES ARE DELIBERATELY NOT "cLogP" AND NOT "MW".

`logp_crippen` is RDKit's Wildman-Crippen estimate. It is NOT ChemDraw's cLogP,
not BioByte's ClogP, and not an experimental logD; the three disagree, routinely
by 0.3-0.5 and sometimes by more than a unit, because they are different fitted
models and not different implementations of one. Labelling it "cLogP" would put a
number in someone's SI that a referee recomputing it in ChemDraw would not get.
The field name carries the method so the reply cannot lose it.

`mw` is the average molecular weight, the one that goes on a bottle and into a
stoichiometry calculation. `exact_mass` is monoisotopic, the one HRMS is compared
against ("HRMS calcd for ..."). They differ by enough to matter -- 180.16 against
180.0423 for aspirin -- and a paper needs both, for different sentences.

tPSA IS ERTL'S, AND IT STILL DOES NOT ALWAYS AGREE WITH PUBCHEM. Measured on
caffeine: RDKit 61.82, PubChem 58.44. The 3.38 is not rounding and not a bug --
the two amide nitrogens of the pyrimidinedione ring are perceived as AROMATIC by
RDKit, which scores an aromatic N with three connections at 4.93, while Cactvs
(what PubChem runs) scores them as tertiary amide N at 3.24. 2 x (4.93 - 3.24) =
3.38 exactly. Every atom contribution was read off to confirm that is the whole
difference. So the number is right for its method and a chemist checking it
against PubChem will see a mismatch; say which method, always.
"""
from __future__ import annotations

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")


def describe(smiles: str) -> dict:
    """Descriptors for one SMILES. Never raises; a bad SMILES comes back as
    {"smiles": ..., "error": ...} so one bad entry in a list does not lose the
    rest, which is the whole point of being handed a list."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"smiles": smiles, "error": "SMILES did not parse"}
    return {
        "smiles": Chem.MolToSmiles(mol),          # canonical, so two spellings
        "formula": rdMolDescriptors.CalcMolFormula(mol),   # of one molecule agree
        "mw": round(Descriptors.MolWt(mol), 2),
        "exact_mass": round(Descriptors.ExactMolWt(mol), 4),
        "tpsa": round(Descriptors.TPSA(mol), 2),
        "logp_crippen": round(Crippen.MolLogP(mol), 2),
        "hbd": Descriptors.NumHDonors(mol),
        "hba": Descriptors.NumHAcceptors(mol),
        "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "rings": rdMolDescriptors.CalcNumRings(mol),
    }


def fill_mw(d: dict, key: str | None = None) -> str | None:
    """Fill `<key>_mw` from `<key>_smiles`, in place. Returns an error or None.

    SHARED because more than one caller needs it (reaction_amounts among them),
    and each would otherwise carry its own copy.
    It lives here because this module already owns SMILES -> molecular weight;
    stoichiometry.py deliberately has no RDKit in it, so it could not go there.

    `key` is None for a plain {"smiles", "mw"} pair and "product"/"substrate" for
    the prefixed ones a route step carries."""
    smi_k = ("%s_smiles" % key) if key else "smiles"
    mw_k = ("%s_mw" % key) if key else "mw"
    if d.get(mw_k) in (None, "") and d.get(smi_k):
        got = describe(str(d[smi_k]))
        if "error" in got:
            return "%s: %s" % (d[smi_k], got["error"])
        d[mw_k] = got["mw"]
    return None


def describe_all(smiles_list) -> list[dict]:
    return [describe(s) for s in (smiles_list or [])]
