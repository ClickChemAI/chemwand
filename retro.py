"""retro.py — the one thing about a retrosynthetic proposal that code can settle.

A retrosynthesis is a guess. Claude can make a good one, and it will make a bad
one in exactly the same confident register, which is a property of the model and
not something a prompt fixes. So the useful division of labour is not "let the
model be careful"; it is to find the part of the claim that is mechanically
falsifiable and falsify it.

"Target <- precursor A + precursor B" contains one such part: EVERY HEAVY ATOM OF
THE TARGET HAS TO COME FROM SOMEWHERE. Atoms cannot appear. If the proposal is
short an atom, the route is wrong no matter how plausible the reaction names are,
and that is checkable without knowing any chemistry at all.

WHAT THIS DOES NOT CHECK, AND IT IS THE IMPORTANT PART. Atom accounting says
nothing about connectivity. Measured: 4-methylbiphenyl from 2-bromotoluene plus
phenylboronic acid accounts for all 13 heavy atoms and passes cleanly -- and it is
wrong, because the methyl ends up ortho, not para. Regiochemistry, stereochemistry,
chemoselectivity, whether the reaction is known, whether the functional groups
survive: none of it is here.

So `not_checked` is part of the return value, not a footnote in a docstring. A
caller that reports "checked, looks fine" without passing that on has overstated
what happened, and the whole reason this file exists is to not do that.

THE PLACEMENTS MUST BE DISJOINT. Taking each precursor's first substructure
match independently lets two precursors both land on the same ring: the correct
disconnection above would score 7 of 13 and look identical to a proposal missing
an entire precursor. Two precursors cannot supply
the same atom.
"""
from __future__ import annotations

import itertools

from rdkit import Chem, RDLogger
from rdkit.Chem import rdFMCS

RDLogger.DisableLog("rdApp.*")

MAX_PLACEMENTS = 60          # per precursor, before falling back to greedy
MAX_COMBINATIONS = 200_000   # total, same

NOT_CHECKED = [
    "Connectivity / regioselectivity -- matching atom counts do not mean the "
    "pieces are joined at the right positions (e.g. 4-methylbiphenyl from "
    "2-bromotoluene + phenylboronic acid accounts for every atom, but that "
    "coupling gives the ortho isomer)",
    "Stereochemistry",
    "Chemoselectivity, functional-group compatibility, and whether protecting groups are needed",
    "Whether this reaction has literature precedent",
    "Reagents, conditions and yield",
]


def _placements(target, prec_smiles: str):
    """Every way this precursor's largest common substructure can sit on the target."""
    m = Chem.MolFromSmiles(prec_smiles)
    if m is None:
        return None
    return _mcs_hits(target, m)


def _mcs_hits(target, m):
    # NO RING RULES. Forming a ring is turning chain atoms into ring atoms, so
    # "ring matches ring only" failed every ring-forming step outright: CuAAC,
    # RCM, lactonisation, a morpholinone closure.
    res = rdFMCS.FindMCS([target, m], timeout=5,
                         atomCompare=rdFMCS.AtomCompare.CompareElements,
                         bondCompare=rdFMCS.BondCompare.CompareAny,
                         ringMatchesRingOnly=False, completeRingsOnly=False)
    if not res.smartsString:
        return []
    # ANY BOND MEANS ANY BOND HERE TOO. The MCS SMARTS records the bond orders it
    # happened to see, so matching it back re-imposed them: in a Diels-Alder the
    # diene's C=C bonds become single and its middle bond double, and both
    # precursors were squeezed onto the product's one C=C.
    q = Chem.RWMol(Chem.MolFromSmarts(res.smartsString))
    for b in list(q.GetBonds()):
        q.ReplaceBond(b.GetIdx(), Chem.BondFromSmarts("~"), preserveProps=False)
    hits = target.GetSubstructMatches(q, uniquify=False,
                                      maxMatches=MAX_PLACEMENTS)
    return [frozenset(h) for h in hits]


def _fill(target, prec_smiles: str, missing: set) -> frozenset:
    """The most of `missing` this precursor can supply, matched against those
    atoms alone.

    A precursor is placed by its LARGEST common piece or not at all, and that
    piece can be the wrong one: Ph3P=CH2's phenyl lands on cyclohexanone's ring,
    collides, and the ylide is dropped -- with it the CH2 it really brings.
    Matching what is still unaccounted for asks the right question for it.

    ONLY THE PRECURSOR'S NON-RING ATOMS. Otherwise any carbon can supply any
    missing carbon: cyclohexanone + PPh3, no ylide at all, passed with a phenyl
    carbon standing in for the CH2. What a left-out precursor really brings is
    the part outside its rings (the ylide's CH2, pinacol's carbons)."""
    m = Chem.RWMol(Chem.MolFromSmiles(prec_smiles))
    ri = m.GetRingInfo()
    for i in sorted((a.GetIdx() for a in m.GetAtoms() if ri.NumAtomRings(a.GetIdx())),
                    reverse=True):
        m.RemoveAtom(i)
    if m.GetNumAtoms() == 0:
        return frozenset()
    m = m.GetMol()
    m.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(m)
    keep = sorted(missing)
    sub = Chem.RWMol(target)
    for i in sorted(set(range(target.GetNumAtoms())) - missing, reverse=True):
        sub.RemoveAtom(i)
    sub = sub.GetMol()
    sub.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(sub)
    hits = _mcs_hits(sub, m)
    best = max(hits, key=len, default=frozenset())
    return frozenset(keep[i] for i in best)


def _best_disjoint(opts):
    """Largest total coverage with no two precursors supplying the same atom.

    Exhaustive while the product of the option counts is small, greedy after that.
    Greedy is flagged in the output rather than passed off as exact: it can
    under-count, which would accuse a correct proposal."""
    # EVERY PRECURSOR MAY ALSO BE UNPLACED. With placement compulsory, a pair
    # whose only placements overlap has NO valid assignment at all, so coverage
    # collapses to zero and the reply accuses a merely incomplete proposal of
    # accounting for none of the target. Measured:
    # 4-methylbiphenyl from 4-bromotoluene plus toluene reported 0 of 13 where the
    # honest answer is 7 of 13 -- one of them fits, the other has nowhere to go.
    slots = [list(o) + [frozenset()] for o in opts]
    total = 1
    for sl in slots:
        total *= len(sl)
        if total > MAX_COMBINATIONS:
            break
    if total <= MAX_COMBINATIONS:
        best, best_pick = frozenset(), [frozenset()] * len(opts)
        for pick in itertools.product(*slots):
            used = frozenset()
            for s in pick:
                if s & used:
                    used = None
                    break
                used |= s
            if used is not None and len(used) > len(best):
                best, best_pick = used, list(pick)
        return best, True, best_pick

    used, pick = frozenset(), [frozenset()] * len(opts)
    for i in sorted(range(len(opts)), key=lambda k: -max((len(s) for s in opts[k]), default=0)):
        for s in sorted(opts[i], key=len, reverse=True):
            if not (s & used):
                used |= s
                pick[i] = s
                break
    return used, False, pick


def check_disconnection(target: str, precursors: list) -> dict:
    """Does every heavy atom of `target` come from one of `precursors`?"""
    t = Chem.MolFromSmiles(target or "")
    if t is None:
        return {"ok": False, "error": "The target SMILES could not be parsed: %r" % target,
                "not_checked": list(NOT_CHECKED)}
    precursors = [str(p) for p in (precursors or [])]
    if not precursors:
        return {"ok": False, "error": "No precursors given.", "not_checked": list(NOT_CHECKED)}

    opts, bad = [], []
    for p in precursors:
        pl = _placements(t, p)
        if pl is None:
            bad.append(p)
            opts.append([])
        else:
            opts.append(pl)
    if bad:
        return {"ok": False,
                "error": "Precursor SMILES could not be parsed: %s" % ", ".join(bad),
                "not_checked": list(NOT_CHECKED)}

    covered, exact, pick = _best_disjoint(opts)
    n = t.GetNumHeavyAtoms()
    # A precursor the search had to leave out may still bring part of what is
    # missing (see _fill). Only those: one already placed has spent its atoms.
    for p, s in zip(precursors, pick):
        left = set(range(n)) - set(covered)
        if not left:
            break
        if not s:
            covered = covered | _fill(t, p, left)
    missing = sorted(set(range(n)) - set(covered))

    per = []
    for p, pl in zip(precursors, opts):
        m = Chem.MolFromSmiles(p)
        best = max((len(s) for s in pl), default=0)
        per.append({"smiles": Chem.MolToSmiles(m),
                    "heavy_atoms": m.GetNumHeavyAtoms(),
                    "supplies_up_to": best,
                    "leaving_atoms": m.GetNumHeavyAtoms() - best})

    out = {
        "ok": not missing,
        "target_heavy_atoms": n,
        "atoms_accounted_for": len(covered),
        "atoms_from_nowhere": len(missing),
        "unaccounted_atom_indices": missing,
        "precursors": per,
        "search_exact": exact,
        "not_checked": list(NOT_CHECKED),
    }
    # ONE PRECURSOR MEANS A FUNCTIONAL-GROUP INTERCONVERSION, AND THEN THE MISSING
    # ATOMS ARE THE REAGENT'S. Measured on a live proposal: Vilsmeier formylation
    # of 4-methylbiphenyl came back "2 atoms from nowhere" -- the CHO carbon and
    # oxygen, which of course arrive from DMF/POCl3 and not from the substrate.
    # The chemistry was right and a fragment count alone called it impossible.
    #
    # Over-accusing is the worse failure. A chemist told that a sound proposal is
    # absurd stops reading the output, and then the check catches nothing at all.
    # So the same count is reported either way and only the reading changes: with
    # two or more fragments, atoms that no fragment supplies are a real problem;
    # with one, they are what the reagent brings.
    out["shape"] = "fgi" if len(precursors) == 1 else "coupling"
    if missing:
        syms = [t.GetAtomWithIdx(i).GetSymbol() for i in missing]
        if out["shape"] == "fgi":
            # WITH ONE PRECURSOR, CODE CANNOT TELL THE TWO READINGS APART, and
            # the verdict says so instead of choosing. Vilsmeier formylation is
            # short a C and an O because DMF brings them; a Suzuki written with
            # only one partner is short a whole ring because a precursor was
            # forgotten. Both are "one precursor, atoms missing". Counting the
            # missing atoms would separate them, but picking the number where a
            # reagent stops and an omission starts is a guess dressed as a rule.
            out["verdict"] = ("%d heavy atom(s) of the target (%s) do not come from this precursor. "
                              "With a single precursor there are two possible readings and the check "
                              "cannot tell them apart: in a functional-group interconversion the "
                              "reagent supplies them (e.g. the CHO of a formylation); in a coupling, a "
                              "precursor is missing. The reaction conditions decide which."
                              % (len(missing), ", ".join(syms)))
        else:
            out["verdict"] = ("This disconnection does not work: %d heavy atom(s) of the target (%s) "
                              "have no source among the given precursors. Atoms cannot appear from "
                              "nothing."
                              % (len(missing), ", ".join(syms)))
    else:
        out["verdict"] = ("Atom accounting passed: all %d heavy atoms have a source. This does not "
                          "mean the reaction is feasible -- see not_checked." % n)
    if not exact:
        out["verdict"] += " (Too many ways to place the precursors, so a greedy search was used; the coverage may be under-counted.)"
    return out
