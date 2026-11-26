"""
smiles_to_frag.py — SMILES -> 2D coordinates -> cdxml_core Atom/Bond lists.

Needs RDKit (`pip install rdkit`).

This is the offline replacement for the web library's `CDW.smilesToCDXML`
(the in-page converter). Because layout happens in rdkit, the Si-SMILES hang
that plagues the web parser does not occur here.

Pipeline per molecule:
  MolFromSmiles -> Kekulize (explicit single/double, NEVER aromatic 1.5)
  -> Compute2DCoords -> scale bonds to 14.4 units -> flip y to ChemDraw's
  y-grows-down convention -> emit Atom/Bond lists -> cdxml_core.fragment().
"""

from __future__ import annotations

import re
from statistics import median

import cdxml_core as C
# The open-source release carries only the hand-written table below.
_IMPORTED_ABBREVIATIONS: list = []

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem import rdDepictor
    try:
        from rdkit.Chem import rdCoordGen           # Schrodinger CoordGen (better fused/spiro layout)
    except ImportError:
        rdCoordGen = None
    try:
        from rdkit.Chem import rdFMCS               # largest common core, for auto_core_patt
    except ImportError:
        rdFMCS = None
except ImportError as e:  # pragma: no cover - environment-dependent
    raise SystemExit(
        "rdkit not found. On the Mac: `pip install rdkit` (or conda). "
        "This module only runs where rdkit is installed."
    ) from e


# Hand-drawn CDXML, consulted before rdkit. A broken library must not stop a figure
# from being drawn -- but it must not disappear either: a bad SMILES spelling in a
# new entry would otherwise take the whole library out, every molecule would quietly
# go back to being computed, and the symptom would be a ferrocene that looks wrong
# for a reason nowhere near the ferrocene. So the reason is kept and reported.
try:
    import depiction_library as _DL
    _DL_ERROR = ""
except Exception as _e:                  # noqa: BLE001
    _DL, _DL_ERROR = None, "%s: %s" % (type(_e).__name__, _e)


_DIR_TO_DISPLAY = {
    Chem.BondDir.BEGINWEDGE: "WedgeBegin",
    Chem.BondDir.BEGINDASH: "WedgedHashBegin",
}


# --- abbreviations: collapse pendant groups to ChemDraw superatoms ------------
# label -> SMARTS for the pendant group. The attachment atom is inferred as the
# matched atom bonded to something OUTSIDE the match. List longest/most-specific
# first so e.g. "OMe" is claimed before a bare "Me" inside it.
#
# THIS TABLE IS ABOUT CONVENTION, NOT WIDTH. Worth stating because the opposite
# is the natural assumption and it is measurably false. Contracting a group only
# narrows a structure when the group is DRAWN wide -- a ring, a chain. When the
# group draws compact and its label is long, contracting makes the cell WIDER,
# because a drawn group can lie vertically and a text label is always horizontal.
# Measured, horizontal extent, same scaffold:
#
#     OBn  +31.0 pt   Bn   +26.6   NHBoc +19.1   Ph  +14.2   Piv +11.9
#     2-Py  +4.8      SiH2Ph +4.7  iPr    +0.8
#     tBu   -1.6      OTBS   -2.9  CO2H   -3.4   SiMe3 -15.0
#
# **OTBS costs 2.9 pt.** So "abbreviate to fit the column" is not what this
# table does, for any row. It exists so a figure reads the way a chemist writes:
# nobody draws CO2H out in a scope table, and published tables print SiMe3 as
# three characters. Width is an occasional side effect.
#
# DO NOT "fix" the negative rows by deleting them or reordering the list. If a
# figure is too wide, the answer is in the rows at the TOP of that table -- rings
# and chains -- or somewhere else entirely.
_HAND_WRITTEN = [
    # ===== groups that CONTAIN a smaller labeled group -> must come FIRST =====
    # nitro-bearing sulfonyls (before NO2). TWO ISOMERS, TWO NAMES, AND THEY ARE
    # CHEMDRAW'S NAMES. ChemDraw's own nickname library ships Ns.cdxml as the
    # ORTHO isomer only and the para one as Nis.cdxml, so a p-nosyl labelled
    # "Ns" carries a label whose ChemDraw expansion is the wrong ring. The
    # structure in the file would still be right -- a Superatom embeds its own
    # atoms -- but the printed name would not.

    # ===== silyl ethers: the TMS pattern is a subset of the rest -> TMS LAST ==
    ("OTBDPS", "[OX2][Si](c1ccccc1)(c1ccccc1)C(C)(C)C"),
    ("OTBS",   "[OX2][Si](C)(C)C(C)(C)C"),
    ("OTIPS",  "[OX2][Si](C(C)C)(C(C)C)C(C)C"),
    ("OTMS",   "[OX2][Si](C)(C)C"),

    # ===== trityl / PMB ethers (PMB has a methoxy -> before OMe) ==============
    ("OTr",    "[OX2]C(c1ccccc1)(c1ccccc1)c1ccccc1"),
    ("OPMB",   "[OX2][CH2]c1ccc(OC)cc1"),

    # ===== carbamates (attach N for secondary NH-, carbonyl C for tertiary) ===
    # Cbz BEFORE OBn: its O-CH2-Ph also matches OBn, so it must claim it first
    ("NHCbz",  "[NX3;H1]C(=O)OCc1ccccc1"),
    ("NHFmoc", "[NX3;H1]C(=O)OCC1c2ccccc2-c2ccccc21"),
    ("NHAlloc","[NX3;H1]C(=O)OCC=C"),
    ("NHTroc", "[NX3;H1]C(=O)OCC(Cl)(Cl)Cl"),

    # ===== the same four on a TERTIARY nitrogen, and the guard is the point ===
    # ChemDraw ships Cbz / Fmoc / Alloc / Troc unguarded, and imported as-is
    # they cost a documented behaviour: an ordinary BENZYL ESTER is a
    # benzyloxycarbonyl group by connectivity, so unguarded Cbz (10 atoms)
    # outranks OBn (8) and `CCC(=O)OCc1ccccc1` comes back labelled Cbz -- a
    # protecting-group name on a plain ester. That is the same defect that got
    # Eoc and COOCH3 excluded from the import; these rows are the guarded form
    # of it. $([CX3][#7]) says out loud what a carbamate is: a carbonyl bonded
    # to nitrogen.
    ("Cbz",    "[C;$([CX3][#7])](=O)OCc1ccccc1"),
    ("Fmoc",   "[C;$([CX3][#7])](=O)OCC1c2ccccc2-c2ccccc21"),
    ("Alloc",  "[C;$([CX3][#7])](=O)OCC=C"),
    ("Troc",   "[C;$([CX3][#7])](=O)OCC(Cl)(Cl)Cl"),

    ("OBn",    "[OX2][CH2]c1ccccc1"),

    # ===== acetal / MOM-type ethers: end in a methoxy/silyl -> before OMe =====
    ("OMEM",   "[OX2][CH2]OCCOC"),
    ("OMOM",   "[OX2][CH2]O[CH3]"),
    ("OSEM",   "[OX2][CH2]OCC[Si](C)(C)C"),
    ("OTHP",   "[OX2]C1CCCCO1"),

    # ===== ester O-groups (before generic methoxy) ===========================
    ("OPiv",   "[OX2]C(=O)C(C)(C)C"),                     # pivaloate
    # Piv NOT on oxygen -- N-Piv above all, which is what surveyed published
    # figures actually have (about 21 of them, acyl on nitrogen). Must follow
    # OPiv, which contains it.
    ("OBz",    "[OX2]C(=O)c1ccccc1"),                     # any O-benzoate (incl. esters)
    ("OAc",    "[OX2]C(=O)[CH3]"),                        # any O-acetate
    ("OMe",    "[OX2][CH3]"),                             # also matches an ester -OMe

    # ===== amide / amine N-groups (attach via N) =============================
    ("NHBoc",  "[NX3;H1]C(=O)OC(C)(C)C"),
    ("NHBz",   "[NX3;H1]C(=O)c1ccccc1"),
    ("NHAc",   "[NX3;H1]C(=O)[CH3]"),
    ("NMe2",   "[NX3]([CH3])[CH3]"),

    # ===== Boc vs CO2tBu: ONE STRUCTURE, TWO NAMES, TOLD APART BY THE PARENT ==
    # These two rows are the same seven atoms. ChemDraw has one nickname for them
    # (Boc / t-BOC) and no CO2tBu at all, so this pair CANNOT be generated from
    # ChemDraw's library -- it is hand-written and has to stay that way.
    #
    # The discriminator is the carbonyl's third neighbour, and the two guards are
    # mutually exclusive by construction: a Boc carbonyl carries =O, -O- and N,
    # so it is bonded to no carbon; a tert-butyl ester's carries =O, -O- and C.
    # A tert-butyl CARBONATE (RO-CO-OtBu) has neither and matches neither row,
    # which is correct -- OtBu below then claims the tBu-O half of it.
    ("Boc",    "[C;$([CX3][#7])](=O)OC(C)(C)C"),
    ("CO2tBu", "[C;$([CX3][#6])](=O)OC(C)(C)C"),

    # ===== O-linked tBu / iPr, AFTER the two rows that contain them ===========
    # Built the way the other 21 O/NH composites here are built: ChemDraw's own
    # nickname is the group WITHOUT the linking heteroatom (its library ships tBu
    # and iPr, and no OtBu or OiPr), so these labels are written here, but every
    # token in them is one ChemDraw parses -- the same shape as OTBS, which works.
    #
    # POSITION IS LOAD-BEARING: Boc and CO2tBu both contain this row's whole
    # match ([OX2] + the tBu), so either one placed after OtBu would be shredded
    # into OtBu + a bare carbonyl.
    ("OtBu",   "[OX2]C(C)(C)C"),
    ("OiPr",   "[OX2][CH1](C)C"),

    # ===== sulfonyls (attach via S) ==========================================

    # ===== phosphine oxide (attach via P) ====================================
    # Written attachment-atom-LAST, unlike every other row here: the P carries the
    # bond, and a chemist writes this substituent Ph2P=O, not P(=O)Ph2. ChemDraw
    # renders the digit as a subscript from the plain text, the same way it draws
    # NO2 as NO~2~, so the label needs no markup.
    ("Ph2P=O", "P(=O)(c1ccccc1)c1ccccc1"),

    # ===== small distinctive groups (attach via their own atom) ==============
    ("NO2",    "[NX3+](=O)[O-]"),                          # after Ns
    ("CF3",    "[CX4]([F])([F])[F]"),
    ("CN",     "[CX2]#[NX1]"),
    # The ACID, not an ester: [OX2H1] is an -OH, so OAc / OBz / OPiv above (whose
    # oxygens carry two heavy neighbours) cannot be claimed by this row.
    # COSTS 3.4 pt of width and is kept anyway -- see the header. A scope table
    # writes CO2H; drawing it out is the thing that would look wrong.
    ("CO2H",   "C(=O)[OX2H1]"),

    # ===== silicon and benzyl: each contains something claimed above ==========
    # SiMe3 must follow OTMS *and* OSEM, both of which contain [Si](C)(C)C, so it
    # cannot live in the silyl-ether block up top where it would outrank OSEM.
    #
    # THE WORST ROW IN THE TABLE BY WIDTH: -15.0 pt, i.e. contracting it makes the
    # cell wider by more than a bond length. A silicon with three methyls draws
    # compact and can lie vertically; "SiMe3" is six characters and always
    # horizontal. It stays because published tables print those six characters
    # -- and because the alternative reading, that this row is a bug someone should
    # optimise away, is exactly what the header exists to prevent.
    ("SiMe3",  "[Si]([CH3])([CH3])[CH3]"),
    ("SiH2Ph", "[SiH2]c1ccccc1"),                          # before Ph
    # Bn NOT on oxygen. Surveyed figures have 36/36 on a sulfoximine SULFUR
    # (S-CH2Ph), which no O-linked row can reach. It sits inside OBn, NHCbz/Cbz
    # and OPMB, so it has to come after all three -- and before Ph, which sits
    # inside IT.
    #
    # ON A HETEROATOM ONLY, and that guard is load-bearing rather than tidy. A
    # bare [CH2]c1ccccc1 outranks Ph, so butylbenzene would come out "Pr-Bn"
    # instead of "Bu-Ph". "Bn" is what a chemist writes
    # for a benzyl ON a heteroatom (O-Bn, N-Bn, S-Bn); a phenyl at the end of an
    # alkyl chain is a Ph, and the ring is the part worth contracting.
    ("Bn",     "[CH2;$([CH2][!#6;!#1])]c1ccccc1"),
    # 2-pyridyl, substituted at the position ortho to N and nowhere else: the
    # attachment atom is INFERRED from which matched atom reaches outside, so a
    # bare pyridine SMARTS would label a 3- or 4-pyridyl "2-Py". The explicit
    # [cH] on the other four ring carbons is what pins the regiochemistry.
    ("2-Py",   "[cH0]1[cH][cH][cH][cH]n1"),

    # ===== bare alkyls and carbocycles: SECOND-TO-LAST, ALL OF THEM ==========
    # Every row in this block sits inside something above it, which is why the
    # whole block is down here rather than filed by chemistry:
    #
    #   tBu     inside OTBS, OTBDPS, OPiv, Piv, NHBoc, Boc, CO2tBu, OtBu
    #   iPr     inside OTIPS, OiPr, iBu
    #   Mes     inside Ph (the ring)
    #   p-Tolyl inside Ts (a tosyl IS a p-tolyl on a sulfonyl) and Ph
    #   1-/2-Naph  inside Ph
    #
    # NAMES ARE CHEMDRAW'S, not invented here: Ad (not "1-adamantyl"), p-Tolyl (not
    # "Tol"), 1-Naph / 2-Naph (not "Np" -- ChemDraw's Np is 4-nitrophenyl).
    # p-Tolyl costs 16.1 pt of label width against "Tol" at size 10 (30.00 vs
    # 13.89, Arial, measured) and is kept anyway, because a label ChemDraw does
    # not know is a latent typography bug and a wide label is only wide.
    #
    # EVERY RING ROW LOCKS TOTAL-H AND DEGREE ON ITS NON-ATTACHMENT ATOMS.
    # Without it Cy claims a decalin ring and a 4-methylcyclohexyl, Ad claims a
    # 3-substituted adamantyl, and p-Tolyl claims a 2,4-xylyl -- the attachment
    # atom is INFERRED from whichever matched atom reaches outside, so an
    # under-specified ring pattern does not fail to match, it matches the wrong
    # ring and names it confidently. (For neutral carbon the H count already
    # forces the degree; `D` is written out because it is the property being
    # asserted, and the next element added here may not be carbon.)
    # alpha before beta is arbitrary (the two are mutually exclusive), but both
    # must precede Ph, which sits inside either one.
    # iBu BEFORE iPr: an isobutyl's CH(CH3)2 is an iPr, so the shorter row would
    # claim it and leave a stray CH2 drawn between two labels.

    # ===== LAST, and the position is load-bearing ============================
    # A phenyl sits inside a dozen rows above (OBn, OBz, NHCbz, OTr, OPMB, Ts,
    # Ph2P=O, NHBz, Fmoc, Mes, p-Tolyl, 1-/2-Naph...), so Ph has to be the last
    # thing that gets to claim a ring. Its guard in _NEEDS_SUBSTANTIVE_PARENT
    # also READS `claimed`, which is only complete once every other row has had
    # its turn.
    ("Ph",     "c1ccccc1"),
]


# --- the full table: hand-written rows + the imported abbreviation table -------
# ORDER IS THE ALGORITHM. _find_groups walks this list ONCE and marks every
# matched atom unavailable, so a row placed too early eats the group a later row
# was for. Three tiers, applied in this order, and the reason for the tiers is
# that **largest-first is not the whole rule**.
#
#   TIER 1  functional groups -- anything containing a heteroatom.
#           Sorted largest-first, which for these IS the containment order: a
#           substructure cannot have more atoms than what contains it, so sorting
#           by heavy-atom count descending is a valid topological order and needs
#           no hand-maintenance across 250 rows.
#
#   TIER 2  bare hydrocarbon skeletons -- groups of only C and H.
#           AFTER every functional group REGARDLESS OF SIZE, and this is the part
#           that size cannot express. A skeleton pattern matches a molecule's
#           PRINCIPAL ring or chain as readily as a pendant one.
#
#   TIER 3  Ph, alone.
#
# **PH-LAST IS POLICY, NOT SIZE. DO NOT FOLD IT INTO THE SORT.** This is written
# out because the size sort looks like it subsumes the tiers and it does not:
# Ph is 6 atoms and OMe is 2, so a pure largest-first table puts Ph FIRST, and
# anisole comes out as the two labels "Ph" and "OMe" with nothing drawn between
# them -- measured, not hypothetical. Ph has to lose to a row a third its size.
# The same holds for every skeleton row against every functional row.
#
# Ph's guard in _NEEDS_SUBSTANTIVE_PARENT also READS `claimed`, so it needs every
# other row to have had its turn before it runs at all.


def _is_skeleton(smarts: str) -> bool:
    """True if the pattern's atoms are all carbon -- i.e. tier 2.

    Read off the SMARTS rather than declared per row: the imported rows are
    generated and nobody is going to hand-maintain a parallel list of which of
    them are hydrocarbons.
    """
    patt = Chem.MolFromSmarts(smarts)
    if patt is None:
        return False
    for a in patt.GetAtoms():
        # generated rows say "#6"; hand-written ones say C/c/CH2/cH0/...
        d = a.GetSmarts()
        if not ("#6" in d or "C" in d.upper().replace("CH", "C")):
            return False
        if any(sym in d for sym in ("#7", "#8", "#14", "#15", "#16", "#9",
                                    "#17", "#35", "#53", "N", "O", "S", "P",
                                    "F", "Cl", "Br", "I", "Si")):
            return False
    return True


def _build_table(hand, imported):
    """hand: [(label, smarts)]; imported: [(label, smarts, heavy)] -> ordered."""
    rows = []
    for label, sma in hand:
        patt = Chem.MolFromSmarts(sma)
        if patt is None:
            raise ValueError("hand-written row %r has unparsable SMARTS" % label)
        rows.append((label, sma, patt.GetNumAtoms()))
    seen = {l for l, _, _ in rows}
    for label, sma, heavy in imported:
        if label in seen:            # a hand-written row wins its own name
            continue
        rows.append((label, sma, heavy))

    tier3 = [r for r in rows if r[0] == "Ph"]
    rest = [r for r in rows if r[0] != "Ph"]
    tier2 = [r for r in rest if _is_skeleton(r[1])]
    tier1 = [r for r in rest if not _is_skeleton(r[1])]
    key = lambda r: (-r[2], r[0])
    return [(l, s) for l, s, _ in
            sorted(tier1, key=key) + sorted(tier2, key=key) + tier3]


_ABBREV_SMARTS = _build_table(_HAND_WRITTEN, _IMPORTED_ABBREVIATIONS)

# Compiled once. 250 rows x every structure in a scope table was 8000+ SMARTS
# compilations per figure; the pattern for a row never changes.
_PATT_CACHE: dict[str, "Chem.Mol"] = {}


def _patt(smarts):
    if smarts not in _PATT_CACHE:
        _PATT_CACHE[smarts] = Chem.MolFromSmarts(smarts)
    return _PATT_CACHE[smarts]


# --- which rows need more than a pendant bond to earn their contraction -------
# `c1ccccc1` matches a molecule's PRINCIPAL ring as readily as a pendant one, and
# contracting that ring is never right: it turns anisole into two labels with
# nothing drawn between them, and a scope table's benzene core -- the thing the
# figure exists to show -- into the word "Ph".
#
# The test is TOPOLOGY, not size: the ring must hang off an atom
# that is outside the ring, that is DRAWN rather than swallowed by another
# contraction, and that carries something else of its own.
#
#   CCCCc1ccccc1   parent CH2, other heavy neighbours -> Ph          (contract)
#   Cc1ccccc1      parent CH3, nothing else on it -> toluene stays whole
#   COc1ccccc1     parent O, but the O is inside the OMe superatom -> ring drawn
#   CC(=O)c1ccccc1 parent C=O, carries a methyl and an oxygen -> Ph   (contract)
#
# A size rule ("leave >= 2 heavy atoms") was rejected for standing in for this
# question: it misfires on small pendants and large cores alike.
#   PhCH2CH2Ph     both CH2s match Bn; without this rule the whole molecule
#                  contracts to "Bn-Bn" and nothing is drawn between them. The
#                  parent-already-claimed half of the test is what stops the
#                  second one.
#   Cc1ccccn1      2-methylpyridine: the ring IS the molecule, so "Me-2-Py" is
#                  not an abbreviation, it is a refusal to draw the compound.
#   Cc1ccc(C)cc1   p-xylene: parent is a CH3 with nothing else -> stays whole
#   CC(C)(C)C      neopentane is not "Me-tBu"
#   CC1CCCCC1      methylcyclohexane is not "Me-Cy"
#
# EVERY BARE ALKYL AND CARBOCYCLE ROW IS IN HERE, and the reason is
# the one above rather than "rings are risky": each of them can be the whole
# molecule with a single methyl hanging off it, and for each of those molecules
# the two-label reading is a refusal to draw the compound. Ad is in for the same
# reason in a less obvious shape -- adamantan-1-ol matches the Ad row exactly and
# correctly, and contracting it leaves the figure with "Ad-OH" and no structure.
# Bz, from the imported table, is in here too, and it is NOT a ring row -- it is
# here because acetophenone is exactly "a benzoyl with one methyl on it".
# Without the guard PhC(=O)CH3 contracts to "Bz" with a lone methyl drawn beside
# it; with it, Bz is refused, Ph then takes the ring off a carbonyl that carries a
# methyl and an oxygen, and the cell reads the way a chemist draws acetophenone.
_NEEDS_SUBSTANTIVE_PARENT = frozenset({"Ph", "Bn", "2-Py", "Bz"})

# EVERY TIER-2 ROW JOINS IT, computed rather than listed. A bare hydrocarbon can
# be the entire molecule with one small substituent on it, and for each of those
# molecules the two-label reading is a refusal to draw the compound: neopentane
# is not "Me-tBu", methylcyclohexane is not "Me-Cy", p-xylene is not
# "Me-p-Tolyl". The 214 imported rows are generated, so a hand-maintained
# parallel list of which of them are hydrocarbons would go stale on the first
# regeneration -- _is_skeleton answers it from the pattern instead.
_NEEDS_SUBSTANTIVE_PARENT = frozenset(
    _NEEDS_SUBSTANTIVE_PARENT
    | {l for l, sma in _ABBREV_SMARTS if _is_skeleton(sma)})


# --- groups this table CANNOT collapse, but that are worth naming -------------
# Recognition only. Nothing here is ever contracted; these patterns exist so that
# "abbreviate was on and this came out drawn in full" can name the GROUP instead
# of the molecule.
#
# WHY THE TWO LISTS ARE SEPARATE. A group can be recognisable and still be the
# wrong thing to collapse. `c1ccccc1` matches a molecule's PRINCIPAL ring as
# readily as a pendant one, so contracting it turns anisole into two labels with
# nothing drawn between them, and turns a scope table's benzene core -- the thing
# the figure exists to show -- into the word "Ph". Recognising it costs nothing
# and tells the caller exactly what the width went on.
#
# WHAT WOULD STILL BE WRONG TO PUT HERE: Me and Et. A methyl is one line segment
# and an ethyl is two; they cost almost no width, most schemes draw them as lines,
# and listing them would put a note on essentially every figure, which is how
# a report stops being read. This is a rule about the REPORT, not
# about the table -- it is the reason to leave a group out of both lists rather
# than to move it between them.
# EMPTY, AND KEPT EMPTY RATHER THAN DELETED. Groups such as Mes, Tol, Np, Cy,
# tBu and iPr belong in _ABBREV_SMARTS, because "recognise it but never contract
# it" is not the intent for any of them; this list is where a group waits while
# its SMARTS is still too loose to trust with a contraction.
#
# THE MECHANISM STILL EARNS ITS KEEP WITH NO ROWS IN IT. _expanded_groups reads
# both lists, and a row of _ABBREV_SMARTS can be left expanded on its own -- it
# may lose the greedy race, or not be a pendant group in this molecule, or be
# refused by _NEEDS_SUBSTANTIVE_PARENT. That reporting path is what turns "the
# figure came out wide" into "Ph x3". This list is the slot for the NEXT group
# that is worth naming before it is worth collapsing.
_REPORT_ONLY_SMARTS = [
]


# --- depiction templates -----------------------------------------------------
# A substrate-scope table only reads if the scaffold is drawn IDENTICALLY in every
# cell -- the eye is supposed to jump straight to the substituent. CoordGen picks
# an orientation per molecule with no knowledge of the series, so the same core
# comes out mirrored, rotated or flipped from one product to the next.
# GenerateDepictionMatching2DStructure fixes the shared substructure to a
# reference depiction and lays out only the rest.
#
# Any product of the series works as the template: as a substructure query, a
# phenyl matches a tolyl, and a methyl matches the first carbon of an ethyl.

_TEMPLATE_CACHE: dict[str, "Chem.Mol"] = {}

# Set on a reference whose coordinates came out of the depiction library, so
# the members templated on it are laid out the same way it was.
CANNED_REF_PROP = "_from_depiction_library"

# Set by _free_coords on a molecule it sent back to RDKit's own layout because
# CoordGen drew a stated E/Z the wrong way round; read only on a template
# reference. Like a canned chair, that shape is one CoordGen would not choose --
# run the members' constrained layout with CoordGen on and it squares the ring
# again, which is what it did to trans-cyclooctene -- so it is switched off.
EZ_REDRAWN_REF_PROP = "_ez_redrawn"

# --- what this batch could not do -------------------------------------------
# Three separate reports, all reset by clear_template_cache() at the start of a
# batch and read by the tool layer afterwards. A figure with a flagged molecule is
# still worth opening, so nothing here raises -- it is instrumentation, and the
# rule is the same as the Allocator's: if nobody reads the instrument, the defect
# ships.

# THE DENOMINATOR. The three lists below say WHICH molecules lost the template and
# why; none of them says how many were ASKED. "4 structures did not match" reads
# as a disaster in a series of 4 and as a footnote in a series of 40, and
# without these two counters the return value could not tell those apart -- see
# the note on `template_match` in chemwand_mcp._figure_report.
#
# Counted, not derived: `placed` is the flag the layout actually ended on, after
# the match AND both quality measures, so attempted - placed is every way a
# template can be lost including ones nothing lists yet.
TEMPLATE_ATTEMPTED: list[str] = []
TEMPLATE_PLACED: list[str] = []

# Molecules whose core did NOT match the template, so they fell back to plain
# CoordGen and may be oriented differently.
TEMPLATE_MISSES: list[str] = []

# Molecules that DID match but came out of the constrained layout with bonds so
# uneven that the template was thrown away. (smiles, templated ratio, free ratio).
TEMPLATE_DISTORTED: list[tuple[str, float, float]] = []

# Molecules whose core came out UNDEFORMED but not in the same place as the rest of
# the series -- rotated, mirrored or shifted as a rigid body. (smiles, bond lengths
# of drift). The cell looks fine on its own and crooked next to its neighbours.
TEMPLATE_DRIFTED: list[tuple[str, float]] = []

# Molecules that MATCHED, but whose copied layout drew a stated E/Z the wrong way
# round -- a trans ring double bond squared into a cis one -- so the template was
# dropped for them and _free_coords redrew them.
TEMPLATE_EZ_REDRAWN: list[str] = []

# Molecules drawn with abbreviate=True where NOT ONE group in _ABBREV_SMARTS
# matched, so the whole structure was drawn out. Without it, silence here would
# look identical to success.
ABBREV_MISSES: list[str] = []

# Per-GROUP, not per-structure: (smiles, label, count) for every recognisable
# group that came out drawn in full. ABBREV_MISSES only fires when a structure
# matches ZERO groups, so a molecule that collapsed its Boc and drew three phenyl
# rings out in full reports nothing there at all -- and three expanded rings is what
# takes a scheme from 552 pt to 708 pt. The count matters: "Ph x3" is a width
# budget, "a structure had a Ph" is not.
ABBREV_EXPANDED: list[tuple[str, str, int]] = []

# Placeholders actually drawn, (smiles, [labels]). Reported on EVERY figure that
# has one, not only on failure, and that is deliberate: a Markush drawing is not a
# molecule. It cannot be searched, its formula is meaningless, and a reader who
# copies it into a query gets nothing. Unlabelled, a dummy atom comes out as an
# unlabelled vertex -- which ChemDraw reads as a methyl -- so the figure is not
# merely un-queryable, it is WRONG, and nothing says so.
PLACEHOLDERS: list[tuple[str, list]] = []

# rlabels keys that no atom in the SMILES claimed: (smiles, key). Almost always a
# typo ({"1": "Ar"} against a SMILES that says [*:2]), and without this the label
# silently falls back to the default and the figure looks plausible.
RLABEL_UNUSED: list[tuple[str, str]] = []

# Depiction templates that are themselves Markush. Only interesting next to
# TEMPLATE_MISSES, which it explains (see template_mol).
TEMPLATE_MARKUSH: list[str] = []

# Placeholders redrawn as floating substituents, (smiles, [labels]). Reported
# because a floating bond is a DRAWING convention: the bond in the file no longer
# runs to the atom the SMILES said it did, so the picture and the connection
# table have deliberately been allowed to disagree.
FLOATING: list[tuple[str, list]] = []

# Asked to float and could not, (smiles, label, why). Falling back to an ordinary
# bond is the right drawing, and an indistinguishable one -- the same silence the
# abbreviation table's misses would otherwise ship in.
FLOAT_SKIPPED: list[tuple[str, str, str]] = []

# Anchor queries that hit nothing, (smiles, smarts). A decoration whose anchor
# matches nothing draws nothing, and a figure missing a highlight looks exactly
# like a figure that was never asked for one.
ANCHOR_MISSES: list[tuple[str, str]] = []

# style "fill" asked for on a match that encloses no ring, (smiles, smarts).
# A native fill is a coloured area bounded by the bonds it names, so a chain or a
# lone atom has nothing to enclose; those fall back to "color". Reported rather
# than done quietly, because the caller asked for one thing and got another and
# the difference is visible on paper.
FILL_NOT_RING: list[tuple[str, str]] = []

# An anchor that matched MORE THAN ONCE, (smiles, smarts, style, matches, atoms).
# A substructure anchor is re-runnable by design -- one spec, every fluorine --
# so N matches is often exactly what was asked for and this is NOT an error. It is
# reported because the other reading is just as common and looks identical on
# paper: "the amide" written as [#6](=[#8])[#7] matches every amide in the
# molecule, and style "dot" draws a disc per ATOM of every match, so one line of
# spec puts NINE discs on a tripeptide, and without this nothing records that.
#
# Both numbers are kept because there are two multipliers and they fail
# differently: matches>1 means the pattern is not as specific as it reads
# ([$(...)] or "maps" narrows it), atoms>matches means each match names a group
# rather than an atom (which "halo" wants and "dot" usually does not).
#
# "color" never lands here: it unions its matches into ONE region a few lines
# above, so by this point it has exactly one.
ANCHOR_MULTI: list[tuple[str, str, str, int, int]] = []

# A squiggle whose ink reaches outside its structure's own bounding box,
# (smiles, smarts, how far). A squiggle DELIBERATELY overhangs the bond it crosses
# -- that is what makes it read as a cut rather than as a decoration -- and the
# box a packer reserves for a cell is the structure's, so the overhang is spent
# out of the gap to the next cell. Reported rather than reserved: widening the box
# would move the structure inside its cell and change every figure that has one,
# and reserving separately would need the extent threaded out to five packers.
# Told the number, the caller can widen a gap; told nothing, a squiggle touching
# its neighbour looks like a packing bug.
SQUIGGLE_OVERFLOW: list[tuple[str, str, float]] = []

# Enhanced-stereo groups that could not be written onto a node,
# (smiles, atom index, "And:1", why). An ABS/AND/OR group is the ONLY thing in
# this pipeline that carries "drawn as one diastereomer, actually racemic" -- lose
# it and the picture still looks right while the embedded data quietly asserts a
# single enantiomer, which is the exact defect enhanced-stereo support exists to
# remove. Dropping one silently would build that defect back in.
#
# Per ATOM and with the reason, not a count: a list which looks complete and is
# not misleads worse than no list --
# "2 groups dropped" sends you looking at the SMILES, "atom 7, And:1, collapsed
# into superatom Bn" sends you at the abbreviation that ate it.
STEREO_GROUP_DROPPED: list[tuple[str, int, str, str]] = []

# A bond style that overwrote something the bond was already saying,
# (smiles, anchor, style, "Display WedgeBegin -> Wavy"). CDXML gives a bond ONE
# Display and ONE Order, so restyling is destructive exactly when the bond was
# already carrying meaning -- and both losses are chemistry, not decoration:
#
#   Display  a wedge or hash IS the stereocentre's configuration. Overwrite it
#            and the drawing stops asserting one, silently.
#   Order    "dative" replaces "2". The connection table now says single.
#
# Reported rather than refused, per the rule all these reports follow: the
# opposite of silent is LOUD, not forbidden. Restyling a wedge to wavy is what a
# scope table does when it cuts the shared fragment away at a stereocentre, so a
# blanket ban would close a legitimate drawing, the way a blanket "reject
# duplicate ids" rule would.
BOND_STYLE_OVERWROTE: list[tuple[str, str, str, str]] = []


# --- how uneven a templated layout is allowed to be --------------------------
# Measured as max bond / min bond, which is scale-free, so it survives the
# rescaling to BOND_LENGTH further down.
#
# _TEMPLATE_FLAT_ENOUGH -- below this, accept without a second thought.
#   Window: every planar layout in the survey (benzene through steroid, spiro-
#   oxindole, 12-crown-4, strychnine) came out at exactly 1.000, and the mildest
#   REAL artifact measured in a real aryl scope table was 1.23 (an 11.67-unit bond
#   among 14.40s). 1.05 sits inside (1.000, 1.23]. It is also the point where the
#   difference stops being visible: 5% of a 14.40 bond is 0.72 units, about the
#   0.60 LineWidth of the line that draws the bond.
#
# _TEMPLATE_DISPERSION_FACTOR -- above this MULTIPLE of the molecule's own free
#   layout, the unevenness is the template's fault, not the molecule's.
#   An absolute limit cannot do this job, and that is the whole finding: a cube
#   or a cyclophane CANNOT be drawn flat with equal bonds, and their honest free
#   layouts measure 1.58 and 2.39 -- worse than the distorted biphenyl (1.80) that
#   motivates this guard. So the question has to be asked relative to the same
#   molecule laid out with no template at all.
#   Window: 1.005 is the worst ratio seen when templating was doing its job
#   (benzyl macrolactam; adamantane, cubane, cyclophane and the aryl series all
#   sat at 1.000-1.003), and 1.23 is the mildest real artifact, as above. 1.10 is
#   inside (1.005, 1.23].
_TEMPLATE_FLAT_ENOUGH = 1.05
_TEMPLATE_DISPERSION_FACTOR = 1.10


# --- how far the core may sit from where the rest of the series put it --------
# The dispersion guard above asks whether the core got BENT. This one asks whether
# it got MOVED, which is the other half of the only promise a depiction template
# makes: the shared scaffold is drawn identically in every cell. Without this,
# nothing checks that promise, and the failure is invisible one cell at a time -- an
# unbent core rotated 40 degrees looks like a perfectly good drawing until it is
# next to its neighbours.
#
# Measured in BOND LENGTHS (the molecule's own median bond), so it is scale-free.
#
# _CORE_DRIFT_TOL window:
#   Ceiling (must not fire) -- 0.001. Series that template correctly agree with
#     each other to essentially nothing: the 6-member aryl series 0.000, the
#     adamantane series 0.001, the cubane series 0.001.
#   Floor (must fire) -- 0.201, o-xylene as a template for indane, where the two
#     frozen methyl carbons have to become a fused ring and CoordGen slides the
#     core instead of stretching a bond. The louder cases are the benzyl
#     macrolactam at 1.765, the cyclophane at 1.734 and biphenyl-as-fluorene at
#     1.898.
#   0.05 sits inside (0.001, 0.201], with 50x margin below and 4x above, and it has
#   the same physical derivation as the dispersion gate: 0.05 of a 14.40 bond is
#   0.72 units, about the 0.60 LineWidth of the line that draws the atom's bonds.
#   A core atom that moves less than the width of its own line cannot be seen.
_CORE_DRIFT_TOL = 0.05

# Reading a core off a symmetric molecule can be done several ways (adamantane's
# amide leaves 4, cubane's 6). Every one of them is the same drawing, so the guard
# takes the reading that agrees BEST with the series and asks its question of that.
# Capped because a highly symmetric core could otherwise enumerate for a long time.
_CORE_MATCH_CAP = 64

# The yardstick is the FIRST molecule of the batch that templated successfully --
# NOT the reference depiction. This is the whole finding behind this guard: the
# constrained layout may place the core in a different frame than the reference mol
# was drawn in, and do it CONSISTENTLY for every member. Measured against the
# reference, the adamantane series drifts 3.455 bond lengths and the cubane series
# 6.864, while agreeing with each other to 0.001 -- those tables are perfect, and a
# guard anchored on the reference would have thrown away every one of their
# templates. What the eye compares is cell against cell, so that is what is
# measured. Stored on the reference Mol itself, so it lives exactly as long as the
# cached reference does: for a SMILES template (what the tools pass)
# clear_template_cache() drops the Mol and the yardstick with it. A caller that
# hands in its OWN Mol keeps the yardstick across batches, which is the same answer
# -- the constrained layout is deterministic, so the next batch's first member
# reproduces those coordinates and measures 0 drift against them.
_YARDSTICK_PROP = "_cdxmlCoreYardstick"


def _bond_dispersion(mol) -> float:
    """max bond / min bond in the current conformer. 1.0 == every bond equal."""
    conf = mol.GetConformer()
    lengths = []
    for b in mol.GetBonds():
        p = conf.GetAtomPosition(b.GetBeginAtomIdx())
        q = conf.GetAtomPosition(b.GetEndAtomIdx())
        lengths.append(((p.x - q.x) ** 2 + (p.y - q.y) ** 2) ** 0.5)
    if not lengths or min(lengths) <= 0:
        return float("inf")            # a zero-length bond is a broken layout
    return max(lengths) / min(lengths)


def _free_coords(mol, coordgen: bool = True):
    """Lay `mol` out with no template -- the same call the fallback path makes."""
    if coordgen and rdCoordGen is not None:
        rdCoordGen.AddCoords(mol)
        # CoordGen returns NaN, silently, for a molecule made only of hydrogen
        # (H+, H-, H2, D2); every other species measured comes back finite. One
        # NaN width turns every later x in the row into NaN. RDKit's own
        # layout handles these, so use it whenever CoordGen did not.
        conf = mol.GetConformer()
        if not all(v == v and abs(v) != float("inf")
                   for i in range(mol.GetNumAtoms())
                   for v in (conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y)):
            AllChem.Compute2DCoords(mol)
        # CoordGen also lays a ring out as a regular polygon and does not always
        # honour a ring double bond's E/Z: trans-cyclooctene came out cis in all
        # three spellings tried, and so did one E bond of humulene (11-ring) and
        # of (E,Z)-cycloocta-1,5-diene. RDKit's own layout drew all 32 measured
        # cases right, so it redraws any molecule whose stated E/Z the drawing
        # contradicts.
        elif _ez_drawn_wrong(mol):
            AllChem.Compute2DCoords(mol)
            mol.SetBoolProp(EZ_REDRAWN_REF_PROP, True)
    else:
        AllChem.Compute2DCoords(mol)


def _ez_drawn_wrong(mol) -> bool:
    """True if any double bond with stated E/Z is drawn the other way round.

    Read back the way every reader of the file will: from the 2D coordinates,
    through a mol block. Atom and bond order survive that, so each bond is
    compared with itself."""
    want = {b.GetIdx(): b.GetStereo() for b in mol.GetBonds()
            if b.GetBondType() == Chem.BondType.DOUBLE
            and b.GetStereo() != Chem.BondStereo.STEREONONE}
    if not want:
        return False
    try:
        back = Chem.MolFromMolBlock(Chem.MolToMolBlock(mol), removeHs=False)
    except Exception:
        return False
    if back is None or back.GetNumBonds() != mol.GetNumBonds():
        return False
    Chem.AssignStereochemistry(back, cleanIt=True, force=True)
    return any(back.GetBondWithIdx(i).GetStereo() != s for i, s in want.items())


def _template_worth_keeping(mol, smiles: str, coordgen: bool) -> bool:
    """Measure the templated layout and say whether to keep it.

    The constrained layout has to attach whatever the template does not cover to a
    core it is not allowed to move, and when it cannot, it stretches a bond instead
    of failing. A 36-unit bond among 14.40s is what that looks like on paper.

    False means: throw the template away and lay this molecule out on its own.
    """
    ratio = _bond_dispersion(mol)
    if ratio <= _TEMPLATE_FLAT_ENOUGH:
        return True
    probe = Chem.Mol(mol)                       # same molecule, no template
    _free_coords(probe, coordgen)
    free = _bond_dispersion(probe)
    if ratio > _TEMPLATE_DISPERSION_FACTOR * free:
        TEMPLATE_DISTORTED.append((smiles, ratio, free))
        return False
    return True                                 # the molecule is just like that


def _median_bond(mol) -> float:
    conf = mol.GetConformer()
    lengths = sorted(
        ((conf.GetAtomPosition(b.GetBeginAtomIdx()).x - conf.GetAtomPosition(b.GetEndAtomIdx()).x) ** 2
         + (conf.GetAtomPosition(b.GetBeginAtomIdx()).y - conf.GetAtomPosition(b.GetEndAtomIdx()).y) ** 2) ** 0.5
        for b in mol.GetBonds())
    return lengths[len(lengths) // 2] if lengths else 1.0


def _core_reading(mol, ref, yard):
    """The core's coordinates, in reference-atom order, as (points, drift).

    A symmetric core can be read off the molecule several ways and they are all the
    same drawing, so every reading is tried and the one that agrees best with the
    yardstick wins. drift is in bond lengths; it is 0.0 when there is no yardstick
    yet, because the first member cannot disagree with itself.
    """
    matches = mol.GetSubstructMatches(ref, uniquify=False, maxMatches=_CORE_MATCH_CAP)
    if not matches:
        return None, 0.0
    conf = mol.GetConformer()
    med = _median_bond(mol)
    best = None
    for m in matches:
        pts = [(conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y) for i in m]
        if yard is None:
            return pts, 0.0
        worst = max(((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2) ** 0.5
                    for p, q in zip(pts, yard)) / med
        if best is None or worst < best[1]:
            best = (pts, worst)
    return best


def _core_placed_alike(mol, ref, smiles: str) -> bool:
    """Is this cell's core in the same place as the rest of the series'?

    False means: the core came through unbent but rotated, mirrored or shifted, so
    the template did not deliver the one thing it is for. Throw it away, draw the
    molecule on its own terms, and say so.
    """
    yard = None
    if ref.HasProp(_YARDSTICK_PROP):
        flat = [float(t) for t in ref.GetProp(_YARDSTICK_PROP).split()]
        yard = list(zip(flat[0::2], flat[1::2]))
    reading = _core_reading(mol, ref, yard)
    if reading[0] is None:
        return True                     # no match to measure against; no verdict
    pts, drift = reading
    if yard is None:
        ref.SetProp(_YARDSTICK_PROP,
                    " ".join("%.6f" % v for p in pts for v in p))
        return True
    if drift > _CORE_DRIFT_TOL:
        TEMPLATE_DRIFTED.append((smiles, drift))
        return False
    return True


# The smallest core that can fix an ORIENTATION. Two atoms are one bond: a
# direction, but nothing that stops the rest flipping across it. Three atoms not
# in a line are a frame. Derived, not tasted -- and in practice the cases this
# rules out come back with 0 or 1 atoms anyway (a series with nothing in common).
_AUTO_CORE_MIN_ATOMS = 3

# An upper bound on FindMCS, not an expected cost. Measured: 24 large, mutually
# different molecules take 0.004 s, and DRAWING that table takes 0.10 s -- so the
# step is ~4% of the figure at its worst. 2 s is 500x the measured worst case and
# 20x the whole figure; a run that hits it returns its best-so-far, which is still
# a usable core, so the ceiling costs correctness nothing.
#
# AN INT, because FindMCS's `timeout` is an unsigned int and 2.0 raises
# ArgumentError -- which the broad `except` below would then swallow, so the
# whole feature would return None in silence.
_AUTO_CORE_TIMEOUT = 2


def auto_core_patt(smiles_list, template=None):
    """The SMARTS this series shares, for `template_patt`. None if there isn't one.

    WHAT IT IS FOR. `template="auto"` takes the FIRST PRODUCT as the reference,
    and a reference carries its own substituent -- so only members whose
    substituent starts with the same kind of atom can match it. Measured on a
    17-member para-substituted benzamide series, the textbook scope table:
    7 of 17 placed. The misses were every member whose para substituent is not
    carbon (F, Cl, Br, OMe, NO2), which is most of a real table.

    The fix is not a different reference -- the reference has to be a drawable
    molecule, and every member is equally lopsided. It is to say WHICH PART of it
    has to match, which is exactly `template_patt`, and the largest common
    substructure is that part by definition. Same series: 16 of 17.

    KEKULIZED FIRST, and this is the whole trick. The pattern is matched against
    a reference that has had clearAromaticFlags run on it, so an MCS computed on
    aromatic molecules comes back with ':,-' bonds that cannot hit a written-out
    double bond -- measured, 0 of 17. Running FindMCS on kekulized copies with
    CompareAny yields '~'-style bonds that survive it.

    NEVER WORSE, by construction: the caller uses this only when no explicit
    `template_patt` was given, and the result is checked against the reference
    before it is used. A series with nothing in common returns None and the
    figure behaves exactly as it would without it.
    """
    if rdFMCS is None or len(smiles_list) < 2:
        return None
    mols = []
    for s in smiles_list:
        m = Chem.MolFromSmiles(s)
        if m is None:
            return None                   # one unparseable member: no opinion
        try:
            Chem.Kekulize(m, clearAromaticFlags=True)
        except Exception:
            return None
        mols.append(m)
    try:
        res = rdFMCS.FindMCS(mols, timeout=_AUTO_CORE_TIMEOUT,
                             bondCompare=rdFMCS.BondCompare.CompareAny)
    except Exception:
        return None
    if not res.smartsString or res.numAtoms < _AUTO_CORE_MIN_ATOMS:
        return None
    patt = Chem.MolFromSmarts(res.smartsString)
    if patt is None:
        return None
    # ...and it has to hit the reference it will be used with. The MCS covers
    # every member INCLUDING the first, so for template="auto" this holds by
    # construction -- but a section that INHERITS the table's reference is drawn
    # on a molecule its own members' MCS never saw. Checking here keeps the
    # promise above (never worse) instead of handing smiles_to_atoms_bonds a
    # pattern it will refuse.
    if template is not None:
        ref = template_mol(template)
        if ref is None or not ref.HasSubstructMatch(patt):
            return None
    return res.smartsString


def template_mol(spec):
    """Build (and cache) a reference molecule that carries 2D coordinates.

    `spec` may be a SMILES string or an rdkit Mol. A Mol without a conformer
    gets one; a Mol with one is used as-is.
    """
    if spec is None:
        return None
    if not isinstance(spec, str):
        if spec.GetNumConformers() == 0:
            _free_coords(spec)
        return spec
    hit = _TEMPLATE_CACHE.get(spec)
    if hit is not None:
        return hit
    ref = Chem.MolFromSmiles(spec)
    if ref is None:
        raise ValueError(f"template SMILES did not parse: {spec!r}")
    if any(a.GetAtomicNum() == 0 for a in ref.GetAtoms()):
        # A Markush parent CANNOT template explicit products. rdkit matches the
        # reference as a MOLECULE, not a query, so a * matches nothing but another
        # * -- measured: c1ccccc1[*:1] fails against toluene, anisole, and even
        # plain benzene, with and without ConstrainedDepictionParams.allowRGroups.
        # Markush against Markush is fine, which is the case that matters for a
        # sectioned Markush table. Recorded so TEMPLATE_MISSES can say WHY.
        TEMPLATE_MARKUSH.append(spec)
    Chem.Kekulize(ref, clearAromaticFlags=True)
    # THE REFERENCE ITSELF GETS THE CANNED DRAWING WHERE THERE IS ONE. Otherwise
    # a series containing a sugar comes out flat even though every member of it
    # could have had the chair: `template="auto"` takes the first cell's SMILES,
    # lays it out with CoordGen -- which flattens a pyranose -- and every other
    # member is then pinned to that flat ring. Measured on a published
    # carbohydrate figure: 7 of the 11 structures the canned core could reach
    # were lost exactly here.
    #
    # Doing it on the reference rather than letting the canned core outrank the
    # caller's template is what keeps the series CONSISTENT: one chair, and every
    # member aligned to it, instead of some cells chairs and some cells flat.
    if not _canned_into(ref):
        _free_coords(ref)       # the NaN guard lives there
    _TEMPLATE_CACHE[spec] = ref
    return ref


def _canned_into(ref) -> bool:
    """Lay `ref` out from a canned depiction if the library has one. In place."""
    if _DL is None:
        return False
    try:
        smi = Chem.MolToSmiles(ref)
        hit = _DL.core_lookup(smi)
        if hit is None:
            return False
        _entry, canned_ref, patt = hit
        prefer = rdDepictor.GetPreferCoordGen()
        rdDepictor.SetPreferCoordGen(False)     # see core_lookup: on, it flattens
        try:
            AllChem.GenerateDepictionMatching2DStructure(
                ref, canned_ref, refPatt=patt, acceptFailure=False)
        finally:
            rdDepictor.SetPreferCoordGen(prefer)
        # Marked so the members templated on it get the same treatment. Without
        # this the reference is a chair and every cell pinned to it is flat
        # again, because the member's own call still had CoordGen on.
        ref.SetBoolProp(CANNED_REF_PROP, True)
        return ref.GetNumConformers() > 0
    except Exception:
        # A reference that will not take the canned core is a reference that gets
        # laid out the ordinary way. Never a failure of the figure.
        return False


def clear_template_cache():
    """Start of a batch: drop cached reference depictions AND every report."""
    _TEMPLATE_CACHE.clear()          # takes the core yardstick with it
    TEMPLATE_ATTEMPTED.clear()
    TEMPLATE_PLACED.clear()
    TEMPLATE_MISSES.clear()
    TEMPLATE_DISTORTED.clear()
    TEMPLATE_DRIFTED.clear()
    TEMPLATE_EZ_REDRAWN.clear()
    ABBREV_MISSES.clear()
    ABBREV_EXPANDED.clear()
    PLACEHOLDERS.clear()
    RLABEL_UNUSED.clear()
    TEMPLATE_MARKUSH.clear()
    FLOATING.clear()
    FLOAT_SKIPPED.clear()
    ANCHOR_MISSES.clear()
    FILL_NOT_RING.clear()
    ANCHOR_MULTI.clear()
    SQUIGGLE_OVERFLOW.clear()
    STEREO_GROUP_DROPPED.clear()
    BOND_STYLE_OVERWROTE.clear()
    # The library keeps its own two reports, and they have to be reset here too:
    # otherwise, once one figure used a canned depiction, every later figure in
    # the same process would claim to have used one too -- naming a molecule
    # from somebody else's call. That matters in the MCP server, which never
    # restarts.
    if _DL is not None:
        _DL.reset_reports()


def _subject(smiles: str, where: str | None) -> str:
    """What a miss report calls this structure.

    The tuple stays three fields wide on purpose -- callers and tests unpack it as
    (subject, smarts, why) -- so the cell label rides in the subject rather than
    widening the shape. "16 products, one anchor missed" is not a report; "cell
    '3a'" is.
    """
    return "%s (%s)" % (where, smiles) if where else smiles


_EXPLICIT_H = re.compile(r"\[\d*H(?![a-z])")


def parse_subject(smiles: str):
    """Chem.MolFromSmiles for the SMILES a figure is ABOUT -- in ONE place.

    Four functions parse the subject: the two anchor resolvers, the bond
    resolver, and the atom/bond builder. They hand each other atom INDICES
    through `index_map`. If any two of them disagree about how many atoms the
    molecule has, every anchor lands on the wrong atom and nothing says so. So
    the decision lives here and all four call it.

    THE DECISION. RDKit drops a written-out hydrogen on the floor: `[H]Br`
    comes back as one Br carrying an implicit H, and there is then no H ATOM to
    point an arrow at. Harmless in every other figure -- and fatal in a
    mechanism, where a third of the curved arrows ever drawn end on a proton
    (31.1% of mech-USPTO-31K's 224,767 arrows, measured).

    Keeping hydrogens ALWAYS would grow an atom on every SMILES anyone has
    already written. So: keep them WHEN THE CALLER WROTE ONE OUT. A standalone
    bracket hydrogen (`[H]` `[2H]` `[H+]` `[H-]`) is not something you type by
    accident. An H inside another atom's bracket (`[nH]` `[NH2]` `[C@@H]`) is
    the ordinary hydrogen COUNT and stays implicit exactly as before -- the
    negative lookahead is what separates the two, and it is also what keeps
    `[Hg]` `[He]` `[Hf]` `[Ho]` `[Hs]` from reading as hydrogen.

    Nothing in this code base writes a standalone [H] (checked: the only source
    hits are the Python name `[HALO]`), so this changes no existing figure.
    """
    if _EXPLICIT_H.search(smiles or ""):
        ps = Chem.SmilesParserParams()
        ps.removeHs = False
        mol = Chem.MolFromSmiles(smiles, ps)
    else:
        mol = Chem.MolFromSmiles(smiles)
    if mol is not None and mol.GetNumHeavyAtoms() > MAX_ATOMS:
        raise ValueError(
            "This structure has %d atoms; the limit is %d. Laying out a structure "
            "that size takes gigabytes of memory on the shared server. For a "
            "polymer, draw the repeat unit; for a long peptide or chain, draw it "
            "in segments." % (mol.GetNumHeavyAtoms(), MAX_ATOMS))
    return mol


# Layout memory and time grow with the SQUARE of the atom count (RDKit's 2D
# coordinate generation): on the Mac Mini 1000 atoms take 0.4 s / 206 MB, 6000
# take 13 s / 3.1 GB, reached within 2 s. The largest real molecule in the
# 2026-10-29 survey, maitotoxin, has 236 heavy atoms; 1000 is about 4x that.
MAX_ATOMS = 1000


def match_node_ids(smiles: str, smarts: str, index_map: dict, base_id: int = 1,
                   where: str | None = None):
    """CDXML node ids for every atom matching `smarts`, one tuple per match.

    THE BRIDGE. Everything decorative -- a highlighted core, a coloured region, a
    dot on an atom -- needs to name atoms, and the two obvious ways of naming them
    are both wrong:

      coordinates die the moment a substituent changes and the layout moves;
      atom INDICES die the moment the SMILES is rewritten, which is a rewrite that
      changes nothing chemically and everything positionally.

    A substructure query survives both, because it names the atoms the way a
    chemist would: "the coumarin core", not "atoms 4, 5, 6, 9, 10, 11". That is
    what makes a decoration RE-RUNNABLE -- one spec over sixteen scope cells, each
    landing on its own copy of the core. Hand-drawing loses this: Cmd+G remembers
    a set of coordinates and forgets that they were ever the coumarin.

    The machinery is already in the repo twice over (_ABBREV_SMARTS matches 36
    groups, draw_scope_table matches a depiction template); this only points it at
    a different job -- addressing instead of typesetting.

    `index_map` is the {mol atom index -> position in the atom list} filled in by
    smiles_to_atoms_bonds for the SAME molecule. Node id is base_id + position,
    which is the rule fragment() emits by, so the ids are valid for the fragment
    string as returned -- BEFORE C.place() renumbers it. Decorate first, place
    after.

    Matched on the molecule as written, before Kekulize, for the same reason
    _find_groups is: an aromatic SMARTS has to see an aromatic ring.
    """
    patt = Chem.MolFromSmarts(smarts)
    if patt is None:
        raise ValueError(f"anchor SMARTS did not parse: {smarts!r}")
    mol = parse_subject(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")
    out = []
    for match in mol.GetSubstructMatches(patt):
        # dict.fromkeys: an atom collapsed into a superatom shares that node with
        # its neighbours, so a match reaching into one would otherwise name the
        # same id several times.
        ids = list(dict.fromkeys(base_id + index_map[a]
                                 for a in match if a in index_map))
        if ids:
            out.append(tuple(ids))
    if not out:
        ANCHOR_MISSES.append((_subject(smiles, where), smarts,
                              _why_no_match(mol, patt)))
    return out


def enclosing_ring_node_ids(smiles: str, node_ids, index_map: dict,
                            base_id: int = 1):
    """Node ids of the SMALLEST ring that contains every one of `node_ids`, or [].

    For the delocalized depiction. IUPAC GR-5.4 puts a delocalized charge symbol
    "centered within the delocalized system", which for an arenium is the five
    carbons and NOT the six-ring -- the two centres are 0.3 bond apart. The
    default here is the ring's centre. Both are one lookup away once you can name
    the ring, and neither can be derived from the anchor alone, so this is that
    lookup.

    Empty when the anchored atoms are not all in one ring (an allyl cation, a
    non-cyclic delocalized system), which is also the case where "the ring's
    centre" has no meaning and the caller must fall back.
    """
    mol = parse_subject(smiles)
    if mol is None:
        return []
    want = {int(i) for i in node_ids}
    best = None
    for ring in mol.GetRingInfo().AtomRings():
        ids = {base_id + index_map[a] for a in ring if a in index_map}
        if want <= ids and (best is None or len(ids) < len(best)):
            best = ids
    return sorted(best) if best else []


def match_node_ids_by_map(smiles: str, smarts: str, index_map: dict,
                          base_id: int = 1, where: str | None = None):
    """Like match_node_ids, but keyed by the anchor's OWN :n map numbers.

    match_node_ids returns ids in SMARTS-match order, which is right for "colour
    these atoms" and wrong the moment the caller has to say WHICH end is which --
    an electron-pushing arrow has a tail and a head, and match order is an
    accident of how the query happened to be written. match_bond_ids already
    reads :1 / :2 for exactly this reason (line 908); this is the atom-list
    sibling of that, so an annotate anchor can name four atoms and say which two
    are the tail bond and which two the head.

    Returns one dict per match: {map number -> node id}. Unmapped query atoms are
    absent, so a caller that wants them still uses match_node_ids. Empty list
    when the SMARTS carries no maps at all -- the caller then falls back.
    """
    patt = Chem.MolFromSmarts(smarts)
    if patt is None:
        raise ValueError(f"anchor SMARTS did not parse: {smarts!r}")
    qmaps = {a.GetAtomMapNum(): a.GetIdx() for a in patt.GetAtoms()
             if a.GetAtomMapNum()}
    if not qmaps:
        return []
    mol = parse_subject(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")
    out = []
    for match in mol.GetSubstructMatches(patt):
        got = {}
        for n, qi in qmaps.items():
            a = match[qi]
            if a in index_map:
                got[n] = base_id + index_map[a]
        if got:
            out.append(got)
    if not out:
        ANCHOR_MISSES.append((_subject(smiles, where), smarts,
                              _why_no_match(mol, patt)))
    return out


def node_ids_by_atom_map(smiles: str, index_map: dict, base_id: int = 1,
                         where: str | None = None) -> dict:
    """{atom map number -> CDXML node id}, read off the SUBJECT's own :n.

    THE OTHER WAY TO NAME AN ATOM, and for a mechanism the right one. Its
    sibling above names atoms through a QUERY: write a SMARTS, let RDKit match
    it, read the query's own :n off the match. That is correct for "colour every
    ketone in this scope table" and it is the wrong tool for "this arrow goes
    from THAT bromide to THAT carbon", because a query cannot say *which* copy:

      * `GetSubstructMatches` uniquifies by ATOM SET, so a SMARTS that spells out
        all four bromines of CBr4 has one match, not four -- `nth` cannot reach
        the near one because there is only ever one.
      * loosen the SMARTS until the four sets differ and you get four matches in
        RDKit's canonical order, which is not the drawn order, so `nth` becomes
        a number found by rendering and looking.

    Both turn up, in that order, on ONE arrow of the Appel mechanism.

    A map number written in the subject is none of that: it is unique by
    construction, it survives any layout, and it needs no matching at all --
    `index_map` already carries SMILES atom order -> node offset, so this is a
    lookup. `[Br-:1].[C:2]([O:3][P+:4])C` plus `tail=[1], head=[2]` says what it
    means and cannot select the wrong bromine.

    Unmapped atoms are simply absent. Raises on a duplicate map number, which
    RDKit accepts and which would make the arrow's end ambiguous -- exactly the
    failure this exists to remove.
    """
    mol = parse_subject(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")
    out = {}
    for a in mol.GetAtoms():
        n = a.GetAtomMapNum()
        if not n:
            continue
        if n in out:
            raise ValueError(
                "atom map :%d appears twice in %s, so an arrow end naming it "
                "would be ambiguous" % (n, _subject(smiles, where)))
        i = a.GetIdx()
        if i in index_map:
            out[n] = base_id + index_map[i]
    return out


def atom_map_classes(smiles: str, index_map: dict, base_id: int = 1) -> dict:
    """{map number -> node ids of every atom SYMMETRY-EQUIVALENT to it, itself first}.

    WHAT ADDRESSING CANNOT DECIDE, and this is the other half of the answer.
    `node_ids_by_atom_map` names one atom exactly, which kills the ambiguity; it
    cannot know that the packer put THAT one of CBr4's four bromines on the far
    side of the carbon, so the arrow flew 4.71 bond lengths across the molecule.
    The four bromines are the same atom chemically, so any of them is a correct
    answer and the readable one is the right answer -- but "readable" is a fact
    about the drawing, so the choice belongs downstream, next to the coordinates.
    This supplies the candidates; decorate picks among them.

    `includeAtomMaps=False` is load-bearing: with maps included (the default)
    RDKit ranks `[Br:2]` apart from its three unmapped twins and every class comes
    back a singleton, which silently loses the equivalence.
    """
    mol = parse_subject(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False,
                                         includeAtomMaps=False))
    by_rank: dict = {}
    for a in mol.GetAtoms():
        by_rank.setdefault(ranks[a.GetIdx()], []).append(a.GetIdx())
    out = {}
    for a in mol.GetAtoms():
        n = a.GetAtomMapNum()
        if not n:
            continue
        i = a.GetIdx()
        sibs = [i] + [j for j in by_rank[ranks[i]] if j != i]
        out[n] = [base_id + index_map[j] for j in sibs if j in index_map]
    return out


def match_bond_ids(smiles: str, smarts: str, index_map: dict, base_id: int = 1,
                   where: str | None = None, style: str = "squiggle"):
    """CDXML node ids of the two atoms of the BOND an anchor names, per match.

    The sibling of match_node_ids, and it has to be a sibling rather than a filter
    over it: that one dedupes with dict.fromkeys and drops atoms missing from
    index_map, which is right for "colour these atoms" and destroys exactly what is
    needed here -- the position-by-position correspondence between the query's
    atoms and the molecule's.

    WHICH TWO ATOMS:

      1. the ones the SMARTS maps :1 and :2. THE DOCUMENTED FORM, because it is
         order-independent and says out loud which bond is meant.
      2. failing that, matched atoms 0 and 1 in SMARTS order.

    Why 2 is a fallback and not the rule: [CX3](=O)[CX3] matches in the order
    C, O, C, so the positional reading lands on the C=O and not on the C-C anyone
    writing that query meant. The maps remove the trap; the fallback is kept only
    so the one-bond case ("[#6][#6]") needs no ceremony.

    Returns [(id_a, id_b), ...]. An anchor that matches but names something that is
    not a drawable bond is REPORTED (ANCHOR_MISSES), not raised: one bad cell out
    of sixteen should not cost the whole figure, and the report is what the other
    anchor failures already use. A malformed SMARTS still raises -- that is the
    caller's spec being wrong rather than this molecule being unusual.
    """
    patt = Chem.MolFromSmarts(smarts)
    if patt is None:
        raise ValueError(f"anchor SMARTS did not parse: {smarts!r}")
    mol = parse_subject(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")

    maps = {a.GetAtomMapNum(): a.GetIdx() for a in patt.GetAtoms()
            if a.GetAtomMapNum()}
    if 1 in maps and 2 in maps:
        qi, qj = maps[1], maps[2]
    elif patt.GetNumAtoms() >= 2:
        qi, qj = 0, 1
    else:
        raise ValueError(
            "style %r marks a BOND, so its anchor needs at least two "
            "atoms; %r matches one. Map the two you mean, e.g. "
            "'[C:1]-[C:2]'." % (style, smarts))

    out, why = [], ""
    for match in mol.GetSubstructMatches(patt):
        ai, aj = match[qi], match[qj]
        if ai not in index_map or aj not in index_map:
            why = ("the anchor names an atom that the drawing does not draw as "
                   "its own node")
            continue
        ida, idb = base_id + index_map[ai], base_id + index_map[aj]
        if ida == idb:
            # Both ends collapsed into ONE node. Distinct from "not bonded": the
            # bond is real and is simply not drawn, because abbreviate=True turned
            # the group it lives in into a single contracted label. Naming that is
            # the difference between "fix your query" and "this bond is inside an
            # abbreviation, expand it or anchor outside".
            why = ("both ends of that bond collapsed into one abbreviated group, "
                   "so there is no drawn bond between them -- anchor a bond "
                   "outside the abbreviation, or drop abbreviate")
            continue
        if mol.GetBondBetweenAtoms(ai, aj) is None:
            why = ("the two atoms the anchor names are not bonded to each other "
                   "(with no ':1'/':2' maps the first two atoms IN SMARTS ORDER "
                   "are taken, and for '[CX3](=O)[CX3]' those are the C and the "
                   "O) -- write the maps, e.g. '[CX3:1](=O)[CX3:2]'")
            continue
        out.append((ida, idb))
    if not out:
        ANCHOR_MISSES.append((_subject(smiles, where), smarts,
                              why or _why_no_match(mol, patt)))
    return out


def _why_no_match(mol, patt) -> str:
    """A guess at why an anchor hit nothing, or "" if there is no obvious one.

    Worth the six lines because ONE cause accounts for most of it and it is
    invisible: rdkit's aromaticity model is more generous than a chemist's pen.
    Coumarin comes back as O=c1ccc2ccccc2o1 -- the lactone ring aromatic, carbonyl
    included -- so the natural query O=C-O matches nothing at all, and neither
    does C=C. Told "no match", anyone would go looking at their ring perception;
    told this, they write [#8]=[#6][#8] and move on.
    """
    try:
        kek = Chem.Mol(mol)
        Chem.Kekulize(kek, clearAromaticFlags=True)
        if kek.GetSubstructMatches(patt):
            return ("the query is aliphatic but rdkit reads this ring as "
                    "aromatic; write it element-wise, e.g. [#6]=[#6] or "
                    "[#8]=[#6][#8]")
    except Exception:                      # noqa: BLE001 - a hint, never a failure
        pass
    return ""


def _placeholder_labels(mol, rlabels, smiles):
    """atom index -> label, for every dummy atom in `mol`.

    WHERE a placeholder sits is carried by the SMILES ([*:1] is ordinary SMILES,
    and rdkit gives it real 2D coordinates, aromaticity and stereo for free); WHAT
    it says is carried by `rlabels`, keyed by the same map number. The two are
    joined here and nowhere else.

    Defaults are chosen so the common figure needs no rlabels at all: a bare * is
    R, and [*:n] is R^n^ -- written with the superscript MARKER rather than as
    "Rn", because ChemDraw's auto-formula mode would set a trailing digit as a
    SUBSCRIPT (see cdxml_core._label_runs).
    """
    rl = {str(k): v for k, v in (rlabels or {}).items()}
    used, out = set(), {}
    for a in mol.GetAtoms():
        if a.GetAtomicNum() != 0:
            continue
        n = a.GetAtomMapNum()
        key = str(n)
        used.add(key)
        out[a.GetIdx()] = rl.get(key, "R" if n == 0 else "R^%d^" % n)
    for key in rl:
        if key not in used:
            RLABEL_UNUSED.append((smiles, key))
    if out:
        PLACEHOLDERS.append((smiles, sorted(set(out.values()))))
    return out


def _float_placeholders(mol, atoms, bonds, slot, plabels, floating, smiles,
                        attach_set=False):
    """Redraw the named placeholders as FLOATING substituents.

    A floating bond says the position is not being specified: it crosses a ring
    edge and stops inside instead of ending on one atom. The SMILES still has to
    attach the placeholder SOMEWHERE -- that attachment stops meaning "the
    substituent is here" and starts meaning "it comes in from this side", which
    is the same choice a hand-drawer makes when they pick which edge to cross.

    Three things are decided here, all of them from floating bonds drawn by hand:

    WHICH DIRECTION. Towards the midpoint of a ring EDGE, not out through a
    vertex. The hand-drawn bonds all cross the middle of an edge, and a bond aimed
    at a vertex would run along the very bond it is replacing. The edge chosen is
    the one whose direction best matches where rdkit already put the placeholder, so
    moving the attachment in the SMILES still steers it.

    HOW DEEP. C.FLOAT_DEPTH of the ring's inradius (see the constant).

    HOW LONG. One ordinary bond, C.BOND_LENGTH. All three hand-drawn ones
    measured 14.40 exactly -- a normal bond, moved.

    Mutates `atoms` and `bonds` in place; returns the labels it floated.
    """
    ring_info = mol.GetRingInfo()
    conf_rings = [set(r) for r in ring_info.AtomRings()]
    floated = []

    for idx, label in plabels.items():
        n = mol.GetAtomWithIdx(idx).GetAtomMapNum()
        if not (floating is True or n in floating):
            continue
        nbrs = [b.GetOtherAtomIdx(idx) for b in mol.GetAtomWithIdx(idx).GetBonds()]
        if len(nbrs) != 1:
            FLOAT_SKIPPED.append((smiles, label, "not a pendant placeholder"))
            continue
        rings = [r for r in conf_rings if nbrs[0] in r]
        if not rings:
            FLOAT_SKIPPED.append((smiles, label, "its attachment atom is in no ring"))
            continue
        ring = min(rings, key=len)

        ph = atoms[slot[idx]]
        pts = [(atoms[slot[a]].x, atoms[slot[a]].y) for a in ring if a in slot]
        if len(pts) < 3:
            FLOAT_SKIPPED.append((smiles, label, "ring atoms were collapsed away"))
            continue
        cx = sum(p[0] for p in pts) / len(pts)
        cy = sum(p[1] for p in pts) / len(pts)

        # edge midpoints of the ring, in ring order
        mids = [((pts[i][0] + pts[(i + 1) % len(pts)][0]) / 2.0,
                 (pts[i][1] + pts[(i + 1) % len(pts)][1]) / 2.0)
                for i in range(len(pts))]
        inradius = min(_hypot(m[0] - cx, m[1] - cy) for m in mids)

        want = _unit(ph.x - cx, ph.y - cy)
        best = max(mids, key=lambda m: (_unit(m[0] - cx, m[1] - cy)[0] * want[0]
                                        + _unit(m[0] - cx, m[1] - cy)[1] * want[1]))
        ux, uy = _unit(best[0] - cx, best[1] - cy)

        stop = (cx + ux * C.FLOAT_DEPTH * inradius,
                cy + uy * C.FLOAT_DEPTH * inradius)
        ph.x, ph.y = stop[0] + ux * C.BOND_LENGTH, stop[1] + uy * C.BOND_LENGTH

        # A PLAIN NODE, which ChemDraw reads as a carbon -- so the ring quietly
        # gains a methyl at a point that draws no ink. That is a real defect and
        # it is chosen anyway, because the alternatives are worse ON PAPER: the
        # only node types that are not atoms and not that methyl are
        # ExternalConnectionPoint and MultiAttachment, and ChemDraw marks the
        # latter with a visible asterisk (compared side by side). It is also
        # exactly what floating bonds drawn by hand in ChemDraw do. The figure is
        # right, the connection table is not, and _chem_report says so out loud
        # rather than letting it pass -- which is the same bargain the floating
        # bond itself is: a drawing that deliberately contradicts the structure.
        atoms.append(C.MultiAttachment(stop[0], stop[1],
                                       tuple(slot[a] for a in ring if a in slot))
                     if attach_set else C.Atom(stop[0], stop[1]))
        inner = len(atoms) - 1
        # The bond no longer belongs to the ring atom: it runs from the label to
        # a point with no atom on it. Rewritten rather than added, or the ring
        # would keep a second, real substituent underneath the floating one.
        # DROPPED and re-APPENDED, not edited where it sat. Two reasons, and the
        # second one is the one that shows on paper:
        #
        #   the bond no longer belongs to the ring atom, so leaving the original
        #   in place would keep a real substituent underneath the floating one;
        #
        #   and DRAW ORDER decides which of two crossing bonds gets the break.
        #   Rewritten in place, the floating bond comes out before the ring bond it
        #   crosses -- so the ring edge draws on top and cuts the floating bond in
        #   half. Hand-drawn ones are the other way round (ring fragment Z=165,
        #   floating R fragment Z=247 in a ChemDraw-drawn reference): the
        #   floating bond runs over, and the RING edge carries the break. Last in
        #   the list is last drawn, which is the only way this writer can say
        #   "in front" -- it emits no Z.
        me = slot[idx]
        bonds[:] = [b for b in bonds if me not in (b.i, b.j)]
        bonds.append(C.Bond(me, inner))
        floated.append(label)

    if floated:
        FLOATING.append((smiles, floated))
    return floated


def _hypot(dx, dy):
    return (dx * dx + dy * dy) ** 0.5


def _unit(dx, dy):
    d = _hypot(dx, dy)
    return (dx / d, dy / d) if d else (1.0, 0.0)


def _pendant(mol, ms):
    """(attach, parent) if `ms` is a single-attachment pendant group, else None."""
    attach = [ai for ai in ms
              if any(nb.GetIdx() not in ms
                     for nb in mol.GetAtomWithIdx(ai).GetNeighbors())]
    if len(attach) != 1:
        return None
    parents = [nb.GetIdx() for nb in mol.GetAtomWithIdx(attach[0]).GetNeighbors()
               if nb.GetIdx() not in ms]
    if len(parents) != 1:
        return None
    return attach[0], parents[0]


def _eligible(mol, label, ms, claimed):
    """(attach, parent) if `ms` may be contracted under `label`, else None.

    THE one definition of "a group that can be contracted", shared by the
    collapser and the reporter. If they used different rules the report would name
    groups that could never have been contracted in the first place -- a
    para-phenylene is not a Ph that got away, it is two attachment points and
    therefore not that group.

    Two tests: pendant shape, then _NEEDS_SUBSTANTIVE_PARENT's topology rule for
    the rows that ask for it.
    """
    pend = _pendant(mol, ms)
    if pend is None:
        return None
    if label in _NEEDS_SUBSTANTIVE_PARENT:
        attach, parent = pend
        # DRAWN, not swallowed: a parent already inside another contraction is a
        # label, and hanging a second label off the first leaves a figure with no
        # structure in it (anisole -> Ph-OMe).
        if parent in claimed:
            return None
        # ...and carrying something of its own, so the ring is a substituent
        # rather than the whole molecule (toluene).
        others = [nb.GetIdx() for nb in mol.GetAtomWithIdx(parent).GetNeighbors()
                  if nb.GetIdx() != attach]
        if not others:
            return None
    return pend


def _find_groups(mol):
    """Non-overlapping pendant-group matches on an aromaticity-perceived mol.
    Returns a list of dicts: {label, atoms(set), attach(idx), parent(idx)}.
    Only single-attachment pendant groups are collapsed."""
    claimed, groups = set(), []
    for label, sma in _ABBREV_SMARTS:
        patt = _patt(sma)
        if patt is None:
            continue
        for match in mol.GetSubstructMatches(patt):
            ms = set(match)
            if ms & claimed:
                continue
            pend = _eligible(mol, label, ms, claimed)
            if pend is None:
                continue
            claimed |= ms
            groups.append({"label": label, "atoms": ms,
                           "attach": pend[0], "parent": pend[1]})

    # NOTHING LEFT TO DRAW. Two rows can tile a molecule between them and the
    # figure is then LABELS WITH NO BOND, because there is no atom left for a
    # bond to join: measured on 2-nitrobenzoic acid (OBz + NO2), on
    # 3,3'-bipyridine (3Py + 3Py) and on N-(quinolin-8-yl)benzamide (Qu + NHBz),
    # each of which comes out as two floating words. RDKit cannot
    # read such a file back at all ("molzip: bond info not properly setup"), so
    # the round trip reports a mismatch -- correctly, and for a reason that reads
    # exactly like the parser limits it is not.
    #
    # _eligible ALREADY refuses one row that swallows everything: a whole molecule
    # has no bond leaving it, so _pendant finds no attachment and returns None.
    # What it cannot see is TWO rows doing it together, because it is asked once
    # per row and each of them does leave something behind.
    #
    # Give the BIGGEST one back, not the last one: it is the row whose absence
    # leaves the most structure drawn, and on all three cases it is also the one a
    # chemist would not have contracted -- the ring-bearing OBz / Qu / 3Py rather
    # than the NO2 / NHBz. Dropping one is always enough, because that row held at
    # least one atom.
    #
    # THE CALLER IS TOLD FOR FREE. The dropped row's atoms go back to unclaimed,
    # so _expanded_groups picks it up and it comes out in the existing
    # "these groups were drawn in FULL anyway" note.
    if groups and sum(len(g["atoms"]) for g in groups) >= mol.GetNumAtoms():
        groups.remove(max(groups, key=lambda g: len(g["atoms"])))
    return groups


def _expanded_groups(mol, claimed):
    """Recognisable groups that came out drawn IN FULL -> [(label, count)].

    Runs after the collapsing pass and skips anything already inside a contracted
    group, so the methoxy inside an OPMB is not reported as an expanded OMe.

    Both tables are searched. A row of _ABBREV_SMARTS can be left expanded too --
    its match may lose the greedy race or not be a pendant group in this
    particular molecule -- and _REPORT_ONLY_SMARTS is everything the table cannot
    contract at all. The caller does not care which list a group came from; it
    cares that the figure got wider than it asked for.
    """
    counts: dict[str, int] = {}
    seen: set = set()
    for label, sma in list(_ABBREV_SMARTS) + list(_REPORT_ONLY_SMARTS):
        patt = _patt(sma)
        if patt is None:
            continue
        for match in mol.GetSubstructMatches(patt):
            ms = set(match)
            # `seen` is what makes the ORDER of the two lists meaningful here:
            # a mesityl is reported once as Mes, not also as the Tol and the Ph
            # that sit inside the same ring.
            if ms & claimed or ms & seen:
                continue
            # The SAME eligibility the collapser uses: a ring the pendant rule
            # deliberately keeps drawn (a principal ring) is not "left expanded",
            # it is drawn on purpose, and reporting it would be noise.
            if _eligible(mol, label, ms, claimed) is None:
                continue
            seen |= ms
            counts[label] = counts.get(label, 0) + 1
    return sorted(counts.items())


# --- enhanced stereochemistry: ABS / AND / OR groups -------------------------
# RDKit's names on the left, CDXML's EnhancedStereoType on the right. Absolute
# takes no group number; the CDX spec says And and Or require one.
_RDKIT_STEREO = {"STEREO_ABSOLUTE": "Absolute",
                 "STEREO_AND": "And",
                 "STEREO_OR": "Or"}


def _stereo_of(mol) -> dict:
    """{mol atom index -> (EnhancedStereoType, group number or None)}.

    Read STRAIGHT off what the caller wrote, and only off that. `@`/`@@` alone
    means absolute configuration and is left alone; a group appears here only
    because the SMILES carried a CXSMILES extension saying so:

        C[C@H](O)[C@@H](C)Cl |&1:1,3|     -> And, group 1, atoms 1 and 3
        ...                  |o1:1,3|     -> Or, group 1
        ...                  |a:5|        -> Absolute

    **Nothing is inferred.** A molecule labelled "dr 3:1" in the paper is not
    guessed into an And group, because the engine cannot tell dr from er by
    looking at a structure, and a wrong guess here does not produce a wrong
    PICTURE -- it produces a file that asserts the wrong chemistry while looking
    correct, which is the failure mode this whole feature exists to remove.
    Guessing would build it back in.
    """
    out = {}
    for g in mol.GetStereoGroups():
        kind = _RDKIT_STEREO.get(str(g.GetGroupType()))
        if kind is None:                      # a future RDKit group type
            continue
        num = g.GetReadId() if kind in ("And", "Or") else None
        for at in g.GetAtoms():
            out[at.GetIdx()] = (kind, num)
    return out


def _stereo_dropped(smiles, idx, spec, why):
    """Record a group that had nowhere to land. See STEREO_GROUP_DROPPED."""
    kind, num = spec
    STEREO_GROUP_DROPPED.append(
        (smiles, idx, kind if num is None else "%s:%d" % (kind, num), why))


def _ate_atom(groups, idx) -> str:
    """Which abbreviation swallowed this atom -- named, so the report is actionable."""
    for g in groups:
        if idx in g["atoms"]:
            return "collapsed into superatom %r" % g["label"]
    return "not written as a plain atom"


def _electron_state(a):
    """(Radical, NumHydrogens) for one RDKit atom, or (None, None) for an ordinary one.

    THE FAILURE THIS CLOSES. `[C]1SC(C)=C(C)N1C` parses, and the canonical SMILES
    comes back with `[C]` intact -- but the CDXML said nothing about the free
    valence, so reading it back gave `CC1=C(C)N(C)CS1`: ChemDraw filled the carbene
    in with two hydrogens and nothing anywhere said so. Measured with RDKit's own
    CDXML parser, which is what `roundtrip_inchi` uses.

    RDKit counts UNPAIRED ELECTRONS, so the mapping is arithmetic, not taste:
    one is a doublet, two is a triplet. Singlet is NOT reachable from SMILES --
    `[C]` and a singlet carbene are the same string -- so this never emits it;
    a caller who means a singlet has to say so some other way, and currently
    cannot.

    NumHydrogens rides along on the SAME atoms and only those. It is what actually
    stops the fill (either attribute alone restores the round trip); it is emitted
    with the radical rather than on every atom because emitting it everywhere would
    move every existing figure's bytes.
    """
    n = a.GetNumRadicalElectrons()
    if not n:
        # A CHARGED ATOM NEEDS THE SAME PIN, and without it every ion comes back
        # neutral. ChemDraw recomputes the hydrogen count from
        # its own valence model when NumHydrogens is absent, and for an ion it
        # reaches a different answer and DROPS the charge to make it fit: an
        # alkoxide came back as HO, a phosphonium as PH. Measured, not guessed --
        # of 6724 charged heteroatoms in a set of ChemDraw files, 6550 carry
        # NumHydrogens and the 102 that carry only Charge + Element were this
        # converter's own output. This DOES move bytes, for any figure with an
        # ion in it. That is the point; without it the bytes draw the wrong
        # species.
        if a.GetFormalCharge():
            return None, a.GetTotalNumHs()
        return None, None
    if n == 1:
        return "Doublet", a.GetTotalNumHs()
    if n == 2:
        return "Triplet", a.GetTotalNumHs()
    # Three or more unpaired electrons has no CDXML spelling: the DTD stops at
    # Triplet. Refused rather than rounded down -- a quartet drawn as a triplet is
    # a different species, and silently drawing one is the failure these pins
    # exist to prevent.
    raise ValueError(
        "%s has %d unpaired electrons; CDXML's Radical attribute enumerates "
        "Singlet / Doublet / Triplet and has no spelling past three."
        % (a.GetSymbol(), n))


def smiles_to_atoms_bonds(smiles: str, wedge: bool = True, coordgen: bool = True,
                          abbreviate: bool = False, template=None,
                          template_patt=None,
                          strict_template: bool = False,
                          template_coordgen: bool | None = None,
                          rlabels: dict | None = None,
                          floating=None, float_attach_set: bool = False,
                          index_map: dict | None = None):
    """Return (atoms, bonds) lists ready for cdxml_core.fragment().

    coordgen=True uses RDKit's CoordGen (Schrodinger layout) for clean fused/spiro
    geometry - fixes the spiro-oxindole self-overlap at the source, no ChemDraw
    cleanup needed. Falls back to Compute2DCoords if CoordGen isn't available.

    template= (SMILES or Mol) fixes the shared substructure to a reference
    depiction so a whole series is drawn in one orientation -- what a scope table
    needs. Molecules whose core does not match fall back to CoordGen and are
    recorded in TEMPLATE_MISSES; pass strict_template=True to raise instead.

    IT IS A COORDINATE DONOR, NOT A QUERY. `template_mol` computes 2D coordinates
    for it and `GenerateDepictionMatching2DStructure` COPIES those coordinates onto
    each member, so it has to be a real, drawable structure (SMILES) -- a SMARTS
    has no coordinates to give. Measured: an aromatic scaffold written as SMARTS
    makes every aromatic member fall back, because the query's aromatic atoms no
    longer match the Kekulized molecule; 7 of 16 tested (template, member)
    pairs regress that way. The tool layer refuses SMARTS-only syntax here
    (preflight.check_template) and points at the argument below.

    template_patt is an optional SMARTS narrowing which part of the reference has
    to match, for when the reference itself is bigger than the shared core.
    A templated layout is MEASURED before it is accepted, against both halves of
    what a template promises. If its bonds came out materially more uneven than the
    same molecule's own untemplated layout the template is discarded and the
    molecule recorded in TEMPLATE_DISTORTED; if the core came through unbent but
    rotated or shifted away from where the rest of the batch drew it, the template
    is discarded and the molecule recorded in TEMPLATE_DRIFTED.

    abbreviate=True collapses pendant groups in _ABBREV_SMARTS to ChemDraw
    superatoms (contracted labels that Expand/Contract in ChemDraw). Matching runs
    BEFORE Kekulize so aromatic SMARTS (e.g. Ts) still match; the whole molecule is
    laid out first, so each collapsed group reuses its in-place coordinates as the
    embedded expansion (expanded view looks right, no Clean Up needed). A molecule
    where nothing matched is drawn in full and recorded in ABBREV_MISSES.

    `smiles` ACCEPTS CXSMILES, and that is worth saying out loud because nothing in
    the parameter name suggests it. RDKit's MolFromSmiles reads the extension block
    itself and this function passes the string through untouched, so enhanced
    stereochemistry needs no new argument and no new syntax:

        "C[C@H](O)[C@@H](C)Cl |&1:1,3|"   AND group 1 -- drawn as one
                                           diastereomer, data says racemic with
                                           the relative configuration defined
        "...                  |o1:1,3|"   OR group 1
        "...                  |a:5|"      absolute

    That is the difference between a figure whose wedges are right and whose
    embedded data says "single enantiomer" anyway, and one where both agree. An
    AND group whose atom cannot be written as a plain node -- collapsed by
    abbreviate=, or drawn as a Markush placeholder -- is reported in
    STEREO_GROUP_DROPPED rather than lost.

    Nothing is inferred: see _stereo_of.
    """
    mol = parse_subject(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")

    # Resolved BEFORE Kekulize, alongside the abbreviation match, for the same
    # reason: it reads the molecule as written. Nothing here depends on
    # aromaticity, but the report it fills does depend on being called exactly
    # once per molecule per batch.
    plabels = _placeholder_labels(mol, rlabels, smiles)

    # match groups while aromaticity is still perceived (Kekulize clears it below)
    groups = _find_groups(mol) if abbreviate else []
    grouped = set().union(*(g["atoms"] for g in groups)) if groups else set()
    if abbreviate:
        # PER GROUP. The strikeout report below fires only when a structure
        # matched nothing at all, which would let a molecule that contracted one
        # group and drew three rings out in full say nothing.
        for label, n in _expanded_groups(mol, grouped):
            ABBREV_EXPANDED.append((smiles, label, n))
    if abbreviate and not groups:
        # Asked for contracted labels and got none. Drawing the whole structure is
        # the right thing to DO -- it is just not something to do quietly, because
        # the figure then looks the same as one where abbreviate was never set, and
        # a missing table row and a molecule with no abbreviable group are
        # indistinguishable on paper.
        ABBREV_MISSES.append(smiles)

    # Kekule form => explicit 1/2/3 bond orders (avoids Order=1.5 -> ring circle).
    Chem.Kekulize(mol, clearAromaticFlags=True)

    if wedge:
        # assign stereo from the SMILES and let rdkit choose wedge directions
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

    placed = False
    ref = template_mol(template)
    if ref is not None:
        # A SMARTS THAT DOES NOT PARSE MUST NOT BECOME `None`. MolFromSmarts
        # returns None on bad input, and `refPatt=None` means "use the whole
        # reference" -- so a typo would quietly widen the match instead of
        # narrowing it, and the figure would come back looking fine. Exposing
        # `template_patt` at the tool layer is what makes this reachable
        # from a recipe, so it is refused here rather than in preflight: this is
        # the module that has rdkit, and preflight deliberately imports no
        # engine.
        patt = None
        if isinstance(template_patt, Chem.Mol):
            # Already a query. The canned-core path (depiction_library.core_lookup)
            # hands one over rather than a SMARTS string, because writing this
            # pattern as SMARTS would mean typing a sugar's stereochemistry out by
            # hand -- error-prone enough to ban outright. It is built from the
            # entry's own molecule instead.
            patt = template_patt
        elif template_patt:
            patt = Chem.MolFromSmarts(template_patt)
            if patt is None:
                raise ValueError(
                    "template_patt did not parse as SMARTS: %r. It is a MATCH "
                    "PATTERN, not a structure -- `template` is the one that has "
                    "to be a drawable molecule." % (template_patt,))
            # ...AND A SMARTS THAT PARSES BUT CAN NEVER MATCH IS THE SAME BUG ONE
            # STEP LATER. The guard above catches a typo; this one catches the
            # pattern everybody actually writes.
            #
            # `ref` was Kekulized with clearAromaticFlags (template_mol, and the
            # member likewise above), so the reference the pattern is tested
            # against has NO aromatic atoms and NO aromatic bonds. Measured on a
            # 17-member para-substituted benzamide series:
            #   c1ccccc1                        0 of 17   lowercase c wants an
            #                                             aromatic flag that is gone
            #   [#6]1[#6][#6][#6][#6][#6]1      0 of 17   an OMITTED SMARTS bond
            #                                             means "single or aromatic"
            #                                             and the ring is now single
            #                                             AND double
            #   [#6]1~[#6]~[#6]~[#6]~[#6]~[#6]1 15 of 17
            #   [#8]=[#6][#7]                   16 of 17  (no ring, so no trap)
            #
            # Unchecked, every one of those zeroes would arrive as 17 ordinary
            # template misses -- indistinguishable from "these molecules genuinely
            # share no core". And the caller is STEERED here twice: preflight
            # refuses a SMARTS in `template=` and points at this argument, and the
            # fallback warning recommends it by name. So the natural thing to
            # write is the thing that cannot work, and it fails silently.
            #
            # Raised, not warned: a pattern matching nothing means the whole
            # series is about to fall back, which is not a figure anyone wants
            # written to disk while a message scrolls past.
            if not ref.HasSubstructMatch(patt):
                raise ValueError(
                    "template_patt %r matches nothing in the template. It is "
                    "tested against the KEKULIZED reference -- aromatic rings "
                    "have been written out as alternating single/double bonds -- "
                    "so a pattern that says 'aromatic' (c1ccccc1) or that leaves "
                    "its ring bonds unwritten ([#6]1[#6][#6][#6][#6][#6]1) can "
                    "never match. Write the ring bonds as ~ (any bond): "
                    "[#6]1~[#6]~[#6]~[#6]~[#6]~[#6]1. A pattern that avoids the "
                    "ring needs no special care: [#8]=[#6][#7] is fine."
                    % (template_patt,))
        # CoordGen has to drive the CONSTRAINED layout too, and it does not by
        # default. GenerateDepictionMatching2DStructure honours a global flag: with
        # it off, the core is pinned and the REST is built by rdkit's legacy 2D
        # embedder, which pays for the pinned core by stretching a bond -- 1.80x on
        # a biphenyl, 2.50x on a naphthalene (measured on a real aryl scope table).
        # With it on, the same two molecules come out at 1.004 with the core still
        # within 0.006 bond lengths of the reference. Saved and restored because it
        # is process-global state that belongs to rdkit, not to this module.
        #
        # EXCEPT WHEN THE TEMPLATE'S WHOLE POINT IS A GEOMETRY COORDGEN WOULD NOT
        # CHOOSE. A canned sugar's reference is a CHAIR, and flattening a ring to
        # its own idea of one is exactly what CoordGen is for: on dapagliflozin,
        # with it on the ring comes out with 0.3 degrees of interior-angle spread
        # against the reference's 105.4, i.e. flat. With it off the spread is
        # 105.4 to the decimal and the bonds land in 1.49-1.55, so the reason it
        # was turned on in the first place does not bite here. `template_coordgen`
        # is separate from `coordgen` because the FREE layout, which is where this
        # falls back to if the gates below throw the template away, still wants it.
        prefer = rdDepictor.GetPreferCoordGen()
        _cg = coordgen if template_coordgen is None else template_coordgen
        if ref is not None and (ref.HasProp(CANNED_REF_PROP)
                                or ref.HasProp(EZ_REDRAWN_REF_PROP)):
            _cg = False                  # a canned drawing, or a ring CoordGen
                                         # would square again
        rdDepictor.SetPreferCoordGen(bool(_cg))
        try:
            AllChem.GenerateDepictionMatching2DStructure(
                mol, ref, refPatt=patt, acceptFailure=False)
            placed = True
        except Exception as exc:
            if strict_template:
                raise ValueError(
                    f"template did not match {smiles!r}: {exc}") from exc
            TEMPLATE_MISSES.append(smiles)
        finally:
            rdDepictor.SetPreferCoordGen(prefer)
        # ...and then MEASURE it, twice, because a depiction template promises two
        # things and can break either one on its own: the core comes through UNBENT
        # (bonds stay even) and IN THE SAME PLACE as the rest of the series. Bend
        # first -- it is the cheaper question, and a bent core is usually a moved
        # one too, so reporting it twice would only be noise.
        if placed and not _template_worth_keeping(mol, smiles, coordgen):
            placed = False
        if placed and not _core_placed_alike(mol, ref, smiles):
            placed = False
        # ...and the chemistry: a stated E/Z the copied layout drew the other way
        # round. On a CoordGen reference a trans ring double bond comes out
        # squared into a cis one; _free_coords below redraws it right.
        if placed and _ez_drawn_wrong(mol):
            placed = False
            TEMPLATE_EZ_REDRAWN.append(smiles)
        # AFTER the measures, because "placed" before them is not the answer the
        # caller gets: a molecule can match the template and still be laid out
        # freely because the constrained layout was worse. This is the one spot
        # where the final verdict for this molecule is known.
        TEMPLATE_ATTEMPTED.append(smiles)
        if placed:
            TEMPLATE_PLACED.append(smiles)

    if not placed:
        _free_coords(mol, coordgen)      # clean layout for complex ring systems
    conf = mol.GetConformer()

    if wedge:
        try:
            Chem.WedgeMolBonds(mol, conf)
        except Exception:
            pass  # wedging is best-effort; connectivity still correct without it

    # scale so the median bond length == ChemDraw's 14.4 units
    lengths = []
    for b in mol.GetBonds():
        p = conf.GetAtomPosition(b.GetBeginAtomIdx())
        q = conf.GetAtomPosition(b.GetEndAtomIdx())
        lengths.append(((p.x - q.x) ** 2 + (p.y - q.y) ** 2) ** 0.5)
    scale = C.BOND_LENGTH / median(lengths) if lengths else 1.0

    def xy(idx):
        p = conf.GetAtomPosition(idx)
        return p.x * scale, -p.y * scale     # flip y: rdkit up -> ChemDraw down

    # 1) ungrouped atoms -> plain Atoms; remember each atom's slot in the list
    # Enhanced-stereo groups are read here, from the molecule as the caller wrote
    # it. Everything that follows can only LOSE one -- an atom that collapses into
    # a superatom or turns into a Markush placeholder stops being a plain <n> --
    # so each of those two exits reports rather than shrugs.
    stereo = _stereo_of(mol)
    atoms, slot = [], {}
    for a in mol.GetAtoms():
        i = a.GetIdx()
        if i in grouped:
            if i in stereo:
                _stereo_dropped(smiles, i, stereo[i], _ate_atom(groups, i))
            continue
        x, y = xy(i)
        slot[i] = len(atoms)
        if i in plabels:
            if i in stereo:
                _stereo_dropped(smiles, i, stereo[i],
                                "drawn as Markush placeholder %r" % plabels[i])
            # GetSymbol() on a dummy returns "*", which _ELEMENTS does not know, so
            # without this it would fall through to the .get(..., 6) default and
            # emit Element="6" -- a carbon, drawn label-free, i.e. a methyl.
            # Verified on the real writer.
            atoms.append(C.placeholder(x, y, plabels[i]))
            continue
        st, gnum = stereo.get(i, (None, None))
        _rad, _nh = _electron_state(a)
        # GetIsotope() is 0 for an unlabelled atom, which cdxml_core reads as
        # "emit nothing". Dropping it draws `[2H]c1ccccc1`, and reads it back,
        # as plain benzene.
        atoms.append(C.Atom(x, y, a.GetSymbol(), charge=a.GetFormalCharge(),
                            stereo_type=st, stereo_group=gnum,
                            radical=_rad, n_hydrogens=_nh,
                            isotope=a.GetIsotope() or None))

    # 2) each matched group -> one Superatom (reusing its in-place coordinates)
    gslot = {}
    for gi, g in enumerate(groups):
        sub_idx = sorted(g["atoms"])
        pos = {ai: j for j, ai in enumerate(sub_idx)}
        sub = []
        for ai in sub_idx:
            x, y = xy(ai)
            at = mol.GetAtomWithIdx(ai)
            _r2, _n2 = _electron_state(at)
            sub.append(C.Atom(x, y, at.GetSymbol(), charge=at.GetFormalCharge(),
                              radical=_r2, n_hydrogens=_n2))
        sub_bonds = [C.Bond(pos[b.GetBeginAtomIdx()], pos[b.GetEndAtomIdx()],
                            int(round(b.GetBondTypeAsDouble())))
                     for b in mol.GetBonds()
                     if b.GetBeginAtomIdx() in g["atoms"]
                     and b.GetEndAtomIdx() in g["atoms"]]
        ax, ay = xy(g["attach"])
        px, py = xy(g["parent"])
        gslot[gi] = len(atoms)
        atoms.append(C.Superatom(ax, ay, g["label"], sub, sub_bonds,
                                 pos[g["attach"]], px, py))

    # 3) bonds: skip group-internal ones; reroute an attachment bond to its superatom
    bonds = []
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        ig, jg = i in grouped, j in grouped
        if ig and jg:
            continue                              # internal to a group -> in sub_bonds
        if not ig and not jg:
            disp = _DIR_TO_DISPLAY.get(b.GetBondDir()) if wedge else None
            bonds.append(C.Bond(slot[i], slot[j],
                                int(round(b.GetBondTypeAsDouble())), disp))
            continue
        core = i if not ig else j
        grp = j if jg else i
        gi = next(k for k, g in enumerate(groups) if grp in g["atoms"])
        # THE WEDGE IS THE CONFIGURATION, and contracting the group must not drop
        # it. C.Bond's fourth argument is `display`; if only the both-ends-plain
        # branch above passed one, every bond that reached a superatom would go
        # out flat. Measured on three cells drawn that way: kainic acid
        # read back /t9? where it was 9-, compound 11 /t19? where it was 19+,
        # periconianone 7 /t16? where it was 16+ -- connectivity identical every
        # time, one centre ERASED rather than inverted, which is why the drawing
        # looked merely unspecified instead of looking wrong.
        #
        # Carried only when the narrow end is the CORE atom. RDKit puts the
        # narrow end at the stereocentre, and every group in the table above is
        # pendant -- none holds a stereocentre -- so that is the case that
        # arises. The other direction would need "WedgedHashEnd", a spelling
        # this file has never written and no reader here is known to accept, so
        # it stays flat rather than being guessed at.
        disp = _DIR_TO_DISPLAY.get(b.GetBondDir()) if wedge else None
        if disp and b.GetBeginAtomIdx() != core:
            disp = None
        bonds.append(C.Bond(slot[core], gslot[gi], 1, disp))  # parent -> superatom

    # LAST: it moves a placeholder and rewrites the bond that reaches it, so it
    # has to see the finished lists. Per placeholder, not per molecule -- one
    # structure routinely mixes a floating R with a fixed one.
    if floating:
        _float_placeholders(mol, atoms, bonds, slot, plabels, floating, smiles,
                            attach_set=float_attach_set)

    if index_map is not None:
        # Handed OUT rather than recomputed by whoever needs it: this mapping is
        # the one place that knows a grouped atom lives on its superatom's node
        # and an ungrouped one on its own, and a second implementation of that
        # rule would drift from this one the first time either changed.
        index_map.update(slot)
        for gi, g in enumerate(groups):
            for ai in g["atoms"]:
                index_map[ai] = gslot[gi]

    return atoms, bonds


def smiles_to_fragment(smiles: str, base_id: int = 1, wedge: bool = True,
                       coordgen: bool = True, abbreviate: bool = False,
                       template=None, template_patt=None,
                       strict_template: bool = False, canned: bool = True,
                       template_coordgen: bool | None = None,
                       rlabels: dict | None = None, floating=None,
                       float_attach_set: bool = False,
                       index_map: dict | None = None,
                       highlight=None, decorations: list | None = None,
                       where: str | None = None):
    """Convenience: SMILES straight to a (fragment_str, bbox) pair.

    A CANNED depiction wins. depiction_library holds hand-drawn CDXML for molecules
    that SMILES cannot express -- ferrocene's iron is not bonded to ten carbons --
    and for those, computing a layout is the wrong move no matter how good the
    layout is. Set canned=False to force the converter (e.g. to see what it does).

    DECORATION HAPPENS HERE, and here is the only place it can. Every one of the
    six figure tools funnels through this function, so a highlight applied here is
    a highlight every tool accepts; a highlight applied in the tool layer is one
    that five of the six will silently drop, which is exactly what doing it here
    prevents. It also has to be here for a second reason: an anchor names
    atoms by the ids fragment() has just handed out, BEFORE C.place() renumbers
    them, and this is the last moment those ids are known to be the ones the
    anchor resolved against.

    `highlight` is the spec list (see decorate.STYLES).

    `decorations` is an out-parameter, the same idiom `index_map` uses. Styles that
    only rewrite the fragment ("color", "fill") need nothing from the caller and
    work everywhere. Styles that produce SIBLING objects -- halo, dot, locants --
    have nowhere to live inside a (fragment, bbox) pair, so the caller must hand in
    a list to receive them. A caller that does not is one that cannot place them,
    and gets a raise rather than a silent drop.

    `where` labels this structure in the miss report ("cell '3a'"). Without it a
    scope table with sixteen products says "1 anchor matched nothing" and leaves
    you to find which one.
    """
    if canned:
        entry = _DL.lookup(smiles) if _DL is not None else None
        if entry is not None:
            frag, bbox = entry.fragment()
            if highlight:
                # A library entry never reaches smiles_to_atoms_bonds, so no
                # index_map is ever built and no anchor can resolve. Say THAT --
                # the aromaticity guess _why_no_match makes would send you looking
                # at your SMARTS, which is not where the problem is.
                for h in highlight:
                    if h.get("on"):
                        ANCHOR_MISSES.append((
                            _subject(smiles, where), h["on"],
                            "this molecule is drawn from the hand-drawn depiction "
                            "library, which carries no atom mapping, so no anchor "
                            "can be resolved against it"))
            # base_id is where the caller wants this fragment's ids to start; the
            # entry hands them out from 1, so shift by the difference.
            return (C.renumber(frag, base_id - 1) if base_id != 1 else frag), bbox

        # NOT the whole molecule, but it CONTAINS one. The canned drawing then
        # becomes a depiction template: its coordinates pin the core and rdkit
        # builds the rest around them. Everything downstream -- wedges, labels,
        # anchors, the drift and evenness gates -- is the ordinary converted
        # path, so a highlight still resolves and a distorted result is still
        # thrown away.
        #
        # This is what makes the sugars worth having. Whole-molecule matching
        # reached NOTHING in a real carbohydrate figure (a published one: 18
        # structures, 16 with a pyranose ring, 0 hits) because a free sugar is
        # not what those figures draw -- a glycoside or a C-glycoside is.
        #
        # The caller's own template wins. Someone who passed `template=` is
        # lining a series up on a shared core and that intent outranks this.
        if template is None and _DL is not None:
            core = _DL.core_lookup(smiles)
            if core is not None:
                _entry, template, template_patt = core
                template_coordgen = False       # see the note at the call site
    # An anchor needs the mapping whether or not the caller wanted one back.
    imap = {} if index_map is None else index_map
    atoms, bonds = smiles_to_atoms_bonds(smiles, wedge=wedge, coordgen=coordgen,
                                         abbreviate=abbreviate, template=template,
                                         template_patt=template_patt,
                                         template_coordgen=template_coordgen,
                                         strict_template=strict_template,
                                         rlabels=rlabels, floating=floating,
                                         float_attach_set=float_attach_set,
                                         index_map=imap)
    frag, bbox = C.fragment(atoms, bonds, base_id=base_id)
    if highlight:
        import decorate as _DEC        # deferred: decorate imports this module
        deco, frag = _DEC.apply(frag, smiles, imap, highlight,
                                base_id=base_id, where=where)
        if deco:
            if decorations is None:
                # NOT "this cannot be done" -- it can, everywhere.
                # cdxml_core.place_with() carries a decoration by the same shift
                # as its structure, and every layout engine calls it. So reaching
                # here means the CALLER did not hand in a sink to receive the
                # decoration, which is a wiring omission in one tool, not a limit
                # of the engine.
                #
                # The distinction is the whole message. "The capability does not
                # exist" tells a caller to give up and use 'fill'; "this path is
                # not wired" tells them to use another tool or to report it.
                # Saying the first would be a false statement dressed as a
                # limitation -- the most expensive kind of error text, because
                # nobody re-tests a "no".
                raise ValueError(
                    "%s: styles halo / dot / locants draw objects that sit BESIDE the "
                    "structure rather than inside it, and this tool cannot place them. Use "
                    "style 'fill' (a native ChemDraw highlight) or 'color', which work "
                    "anywhere."
                    % (where or smiles))
            decorations.append(deco)
    return frag, bbox


if __name__ == "__main__":
    # quick self-check when run directly
    for smi in ["C1=CC=CC=C1", "O=S=NC1=CC=CC=C1"]:
        frag, bb = smiles_to_fragment(smi)
        n_nodes = frag.count("<n ")
        n_bonds = frag.count("<b ")
        print(f"{smi:28s} -> {n_nodes} atoms, {n_bonds} bonds, "
              f"bbox w={bb[2]-bb[0]:.1f} h={bb[3]-bb[1]:.1f}")
