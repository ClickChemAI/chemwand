"""
stereo_tools.py — read / audit stereochemistry for reaction schemes.

DESIGN RULE (the whole point): NEVER invent stereochemistry. These functions read
what ChemDraw already perceived from the wedges you drew (via an SDF export or
Copy-as-SMILES), report it, and audit it. The drawing is the source of truth; the
tool only reads + verifies. RDKit reads SDF/MOL (not CDXML), so ingest takes an
SDF exported from ChemDraw (File > Save As > SDF) or a list of isomeric SMILES.

Runs on the Mac (needs rdkit). Imported by chemwand_mcp.py and wrapped as
tools; safe to import alongside the existing modules (no name clashes).
"""

from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import rdFMCS


# ---------------------------------------------------------------- stereo census
def _centers(mol):
    """Return (assigned, unassigned):
        assigned   = [(atomIdx, 'R'|'S'), ...]   defined tetrahedral centers
        unassigned = [atomIdx, ...]              potential centers with NO config
                                                 (wavy bond / left undefined)
    Uses the modern CIP labeller, falling back to the legacy one on old rdkit."""
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    try:
        cc = Chem.FindMolChiralCenters(mol, includeUnassigned=True,
                                       useLegacyImplementation=False)
    except TypeError:  # very old rdkit without the kwarg
        cc = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    assigned = [(i, lab) for i, lab in cc if lab in ("R", "S")]
    unassigned = [i for i, lab in cc if lab not in ("R", "S")]
    return assigned, unassigned


# --------------------------------------------------------------------- ingest
def ingest(sdf_path: str | None = None, smiles: list | None = None):
    """Read structures and return records with canonical ISOMERIC SMILES + census.

    sdf_path : an SDF exported from ChemDraw (perceives stereo from 2D wedges).
    smiles   : a list of isomeric SMILES (e.g. ChemDraw Copy-as-SMILES).
    Returns list of dicts: {id, ok, smiles, assigned, unassigned}."""
    mols = []
    if sdf_path:
        supp = Chem.SDMolSupplier(sdf_path)           # 2D wedges -> chirality on read
        for i, m in enumerate(supp):
            name = m.GetProp("_Name") if (m is not None and m.HasProp("_Name")
                                          and m.GetProp("_Name").strip()) else str(i + 1)
            mols.append((name, m))
    if smiles:
        for i, smi in enumerate(smiles):
            mols.append((str(i + 1), Chem.MolFromSmiles(smi)))

    records = []
    for name, m in mols:
        if m is None:
            records.append(dict(id=name, ok=False, smiles=None,
                                assigned=[], unassigned=[]))
            continue
        assigned, unassigned = _centers(m)
        records.append(dict(id=name, ok=True,
                            smiles=Chem.MolToSmiles(m, isomericSmiles=True),
                            assigned=assigned, unassigned=unassigned))
    return records


def format_ingest(records) -> str:
    """Human-readable table for the ingest records."""
    out = ["id     ok    #ctr  undef(atomIdx)        isomeric SMILES",
           "-" * 92]
    for r in records:
        if not r["ok"]:
            out.append(f"{r['id']:6s} FAIL  -     -                     <SMILES/mol did not parse>")
            continue
        n = len(r["assigned"]) + len(r["unassigned"])
        und = ",".join(str(i) for i in r["unassigned"]) or "-"
        out.append(f"{r['id']:6s} ok    {n:<5d} {und:21s} {r['smiles']}")
    return "\n".join(out)


# ----------------------------------------------------------------------- lint
def lint(entries):
    """Per-structure stereo audit.  entries = [(id, smiles), ...].

    Flags (a) parse failures and (b) UNDEFINED stereocenters — centers that still
    need a source before the scheme is complete (a deliberately-racemic center,
    like a wavy bond, shows here too so it's never silently dropped).

    Returns (report_str, n_problems).

    NOTE: 'is this center the SAME one, conserved down the route?' needs
    atom-mapping / MCS across consecutive structures — that's a separate, guarded
    layer (propagate), not this function. lint only reports what each structure
    independently says, which is 100% reliable."""
    out = ["id     status  #ctr  undefined(atomIdx)   CIP labels (atomIdx:R/S)",
           "-" * 88]
    problems = 0
    for cid, smi in entries:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            out.append(f"{cid:6s} FAIL    -     -                    -")
            problems += 1
            continue
        assigned, unassigned = _centers(m)
        if unassigned:
            problems += 1
        labs = " ".join(f"{i}:{l}" for i, l in assigned) or "-"
        und = ",".join(str(i) for i in unassigned) or "-"
        status = "ok    " if not unassigned else "UNDEF "
        out.append(f"{cid:6s} {status}  {len(assigned) + len(unassigned):<5d} "
                   f"{und:20s} {labs}")
    return "\n".join(out), problems


# --------------------------------------------------------------------- gallery
def gallery_cdxml(entries, header_prefix, per_row: int = 4):
    """Render every structure WITH WEDGES into one CDXML, in a grid, captioned by
    id — for a visual diff against your drawn scheme. entries = [(id, smiles), ...].
    Reuses the same pipeline (smiles_to_fragment wedge=True) so what you see here is
    exactly what draw_linear_scheme would emit.

    Returns (document, skipped); skipped is [(id, smiles, reason), ...] for
    entries that could not be built. One production caller, so the tuple
    costs nothing and needs no module-level list to reset."""
    import cdxml_core as C
    from smiles_to_frag import smiles_to_fragment

    gx, gy, base = 175.0, 150.0, 1000
    parts, col, row, skipped = [], 0, 0, []
    for cid, smi in entries:
        try:
            frag, _ = smiles_to_fragment(smi, base_id=1, wedge=True)
        except Exception as e:
            # NAMED, NOT SKIPPED. A gallery exists to be read by eye, so a
            # structure that vanishes from it is one the reader scans past and
            # counts as checked -- the one failure this tool must not have.
            skipped.append((cid, smi, "%s: %s" % (type(e).__name__, e)))
            continue
        px, py = col * gx, row * gy
        parts.append(C.place(frag, base, px, py))
        parts.append(C.number(base + 480, px + 40, py + 95, str(cid)))
        base += 600
        col += 1
        if col >= per_row:
            col, row = 0, row + 1
    return C.document(header_prefix, "\n".join(parts)), skipped


# ------------------------------------------------------------------ propagate
def _perm_sign(perm) -> int:
    """Sign (+1/-1) of a permutation given as a list of distinct ints (inversion parity)."""
    inv = 0
    for i in range(len(perm)):
        for j in range(i + 1, len(perm)):
            if perm[i] > perm[j]:
                inv += 1
    return -1 if (inv & 1) else 1


def _mcs_amap(anchor, target):
    """Map anchor atomIdx -> target atomIdx over their maximum common substructure
    (elements + bond orders match; rings matched as whole rings). {} if none."""
    res = rdFMCS.FindMCS([anchor, target],
                         atomCompare=rdFMCS.AtomCompare.CompareElements,
                         bondCompare=rdFMCS.BondCompare.CompareOrder,
                         ringMatchesRingOnly=True, completeRingsOnly=True,
                         timeout=10)
    if not res.smartsString:
        return {}
    patt = Chem.MolFromSmarts(res.smartsString)
    am = anchor.GetSubstructMatch(patt)
    tm = target.GetSubstructMatch(patt)
    if not am or not tm or len(am) != len(tm):
        return {}
    return {am[i]: tm[i] for i in range(len(am))}


def _reapply_chirality(anchor, target, amap):
    """Copy each anchor stereocenter's configuration onto its amap'd target atom,
    correcting for the two atoms' differing neighbor orders via permutation parity.

    A center is left UNTOUCHED (its target idx returned) when a substituent AT the
    center changed — i.e. an anchor neighbor isn't in the shared core (e.g. a
    benzylic OH->Cl). That's a center needing its own source, never a silent guess.
    Returns the list of such unresolved target idxs."""
    unresolved = []
    for a_idx, t_idx in amap.items():
        aa = anchor.GetAtomWithIdx(a_idx)
        tag = aa.GetChiralTag()
        if tag not in (Chem.ChiralType.CHI_TETRAHEDRAL_CW,
                       Chem.ChiralType.CHI_TETRAHEDRAL_CCW):
            continue
        a_nbr = [n.GetIdx() for n in aa.GetNeighbors()]
        if not all(n in amap for n in a_nbr):        # a substituent at the center changed
            unresolved.append(t_idx)
            continue
        a_nbr_as_t = [amap[n] for n in a_nbr]        # anchor's neighbor order, as target idxs
        ta = target.GetAtomWithIdx(t_idx)
        t_nbr = [n.GetIdx() for n in ta.GetNeighbors()]
        if sorted(a_nbr_as_t) != sorted(t_nbr):      # neighbor SETS differ -> can't transfer
            unresolved.append(t_idx)
            continue
        # both centers carry one implicit H in the same slot convention, so aligning the
        # heavy neighbors + a parity correction reproduces the 3D configuration.
        ta.SetChiralTag(tag)
        perm = [t_nbr.index(x) for x in a_nbr_as_t]
        if _perm_sign(perm) == -1:
            ta.InvertChirality()
    return unresolved


def propagate(anchor_smiles, target_smiles):
    """Transfer perceived chirality from anchors onto flat targets.

    anchor_smiles : list of ISOMERIC SMILES (defined-stereo reference structures).
    target_smiles : list of (id, flat_smiles).
    Each target uses the anchor sharing the largest scaffold. Conserved centers are
    filled; centers whose substituent changed, or that are new, are left UNDEFINED
    and reported. Returns list of dicts {id, smiles, unresolved}."""
    anchors = []
    for s in anchor_smiles:
        m = Chem.MolFromSmiles(s)
        if m is not None:
            Chem.AssignStereochemistry(m, cleanIt=True, force=True)
            anchors.append(m)
    out = []
    for cid, tsmi in target_smiles:
        t = Chem.MolFromSmiles(tsmi)
        if t is None:
            out.append(dict(id=cid, smiles=None, unresolved=[]))
            continue
        best_amap, best_anchor = {}, None
        for a in anchors:
            amap = _mcs_amap(a, t)
            if len(amap) > len(best_amap):
                best_amap, best_anchor = amap, a
        unres = _reapply_chirality(best_anchor, t, best_amap) if best_anchor else []
        Chem.AssignStereochemistry(t, cleanIt=True, force=True)
        out.append(dict(id=cid, smiles=Chem.MolToSmiles(t, isomericSmiles=True),
                        unresolved=unres))
    return out


def format_propagate(results):
    """Report propagated SMILES, then AUTO-RUN the lint census on the output (guardrail)."""
    lines = ["id     filled  isomeric SMILES", "-" * 92]
    lint_pairs = []
    for r in results:
        if r["smiles"] is None:
            lines.append(f"{r['id']:6s} -       <target did not parse>")
            continue
        lint_pairs.append((r["id"], r["smiles"]))
        lines.append(f"{r['id']:6s} yes     {r['smiles']}")
    report, n = lint(lint_pairs)     # GUARDRAIL: the propagated output is audited
    return ("\n".join(lines)
            + "\n\n--- auto stereo_lint on the propagated output ---\n" + report
            + f"\n\n{n} propagated structure(s) still have an undefined center "
              "(a new center — give it its own source anchor, then re-run).")


if __name__ == "__main__":
    # ---- read/audit demo ----
    print("### ingest + lint ###")
    demo = ["OC[C@@H]([C@H](C1=CC=C(C=C1)Br)O)N", "C[C@H](N)C(=O)O", "CC(N)C(=O)O"]
    print(format_ingest(ingest(smiles=demo)))
    print()
    rep, n = lint([(str(i + 1), s) for i, s in enumerate(demo)])
    print(rep)
    print(f"\nproblems: {n}  (expected 1 — the 3rd has an undefined center)")

    # ---- propagate self-test: anchors -> targets ----
    print("\n\n### propagate ###")
    anchors = [
        # 2 (aminodiol) -> S,S ; 8 (azetidine) -> R,R,S  (perceived earlier via stereo_ingest)
        "OC[C@@H]([C@H](C1=CC=C(C=C1)[N+]([O-])=O)O)N",
        "O=S(C1=CC=CC=C1[N+]([O-])=O)(N(C[C@H]2[C@H]([C@H](N2CC=C)COC(C3=CC=CC=C3)"
        "(C4=CC=CC=C4)C5=CC=CC=C5)C6=CC=C(C=C6)Br)CC=C)=O",
    ]
    targets = [
        ("3",     "OCC(N)C(O)c1ccc(Br)cc1"),                                  # aminodiol, NO2->Br
        ("4",     "OCC(N(CC=C)CC#N)C(O)c1ccc(Br)cc1"),                        # N-alkylated
        ("6",     "N#CC1C(c2ccc(Br)cc2)C(COC(c3ccccc3)(c4ccccc4)c5ccccc5)N1CC=C"),  # azetidine
        ("5core", "ClC(c1ccc(Br)cc1)C(N)CO"),                                 # benzylic OH->Cl (subst. changed)
    ]
    print(format_propagate(propagate(anchors, targets)))
    print("\nEXPECT: 3 & 4 come back S,S; 6 comes back R,R,S; 5core fills the amino")
    print("center (S) but FLAGS the benzylic as undefined (its substituent changed).")

