"""pubchem.py — a named compound's structure, read from PubChem, never from memory.

WHY THIS EXISTS. The InChI gate reads a figure back and asks whether the drawing
is the SMILES it was given. Nothing asks whether the SMILES it was given is the
compound the name means -- and for a named compound (a reagent, a drug, a natural
product) that SMILES is otherwise written from a model's memory. Looking named
compounds up in PubChem by hand catches exactly this: SuFEx-IT, written from
memory, had one methyl too few. This is that step as one call.

WHAT LEAVES THE MACHINE: the name, and nothing else. It goes to PUG-REST as POST
form data -- in the body, not the URL, so "Pd/C" is not read as a path -- and no
structure is ever sent. There is deliberately no structure -> name direction; it
would put an unpublished structure in NCBI's logs. A compound number ("3a") sent
as a name leaves the machine too, and comes back not found; keeping those out is
the caller's rule, and the tool description says so.

WHAT PUBCHEM DOES WITH A NAME, measured on 48 names. Each warning
below exists because of a row in this table:

  (S)-BINAP        -> "Binap, (A+-)-" (PubChem's own    Axial chirality has no
                      spelling of (+-)), no stereo.       SMILES notation: asked
  (R)-BINOL        -> "Binol", no stereo at all.          for one enantiomer,
                                                          handed the racemate,
                                                          and no error.
  tartaric acid    -> L-tartaric acid          A name that does not say which
  lutidine         -> 2,6-lutidine             isomer gets one chosen for it, and
  TMS              -> 2,3',4,5'-tetramethoxystilbene      an abbreviation can land
                                                          on an unrelated compound.
                                                          Only the title shows it,
                                                          so a title that is not
                                                          the name asked for is
                                                          said out loud.
  menthol          -> 3 stereocentres undefined ("(-)-menthol": none)
  DIAD             -> two records: the (E) one, and one with N=N unassigned
  n-butyllithium   -> "[Li].CCC[CH2]" -- a lithium atom beside a butyl radical
  Pd(PPh3)4        -> a Pd atom and four unbonded PPh3
  Grubbs II        -> Ru bonded to the NHC, PCy3 floating free
  Pd/C             -> a bare Pd atom
  DMP, TEA, PCC,   -> not found
  CAN, DME, 3a

THE SMILES IS PASSED THROUGH VERBATIM. It is the string the PubChem page shows,
so a chemist can check it against the link, and rewriting a SMILES is how
chirality gets flipped without any gate noticing. PubChem renamed the property in
2025 (IsomericSMILES -> SMILES), and asking for the old name is not an error: the
field is silently absent. So a reply without "SMILES" is refused here, and
ConnectivitySMILES -- which carries no stereochemistry -- is never used in its
place.

NOT DONE HERE, on purpose: no cache and no retry. Unassigned double bonds are
not counted, because stereo_lint does not count them either; DIAD's surfaced as
two records instead.

THE THROTTLE IS MACHINE-WIDE. One batch runs a name at a time, 0.2-0.55 s each,
which alone stays under the 5 a second PubChem asks for -- and this file used to
rely on that. The server then got a pool of worker processes, each looking names
up without knowing about the others: 8 at once is about 26 a second, and a block
would cut the whole institute's address off. So every request first takes the
next free slot from one stamp file shared by all processes on the machine.

    lookup("DIAD") -> {"query", "found", "records": [...], "warnings": [...]}
"""
from __future__ import annotations

try:                       # Mac / Linux
    import fcntl
    def _lock(f):
        fcntl.flock(f, fcntl.LOCK_EX)
except ImportError:        # Windows has no fcntl; lock byte 0 with msvcrt
    import msvcrt
    def _lock(f):
        f.seek(0)
        msvcrt.locking(f.fileno(), msvcrt.LK_LOCK, 1)
import json
import re
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request

from rdkit import Chem, RDLogger

import stereo_tools as ST

RDLogger.DisableLog("rdApp.*")

URL = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/property/"
       "SMILES,InChIKey,MolecularFormula,Title/JSON")
PAGE = "https://pubchem.ncbi.nlm.nih.gov/compound/%d"

# Chosen against the measured 0.2-0.55 s per request.
TIMEOUT_S = 10.0

# Against PubChem's 5 a second, with room for it lowering its limit under load.
RATE_PER_S = 4.0
# Against the measured 40 names in 12.3 s; a large scope table fits.
MAX_NAMES = 50
# Against the server's 45 s per call: a call it kills loses every answer it
# had, so past this the names not yet sent come back marked busy instead.
BUDGET_S = 30.0
# Every worker process on the machine reads and writes this one file.
_STAMP = __import__("pathlib").Path(tempfile.gettempdir()) / "chemwand_pubchem_slot"


class Busy(Exception):
    """The next free slot is past the caller's deadline; nothing was sent."""


def _take_slot(deadline: float | None) -> None:
    """Wait for this machine's next request slot, RATE_PER_S apart across every
    process. The slot is reserved under the lock and waited for outside it, so
    one waiting worker does not hold the others."""
    with open(_STAMP, "a+") as f:
        _lock(f)
        f.seek(0)
        try:
            last = float(f.read() or 0)
        except ValueError:
            last = 0.0
        at = max(time.time(), last + 1.0 / RATE_PER_S)
        if deadline is not None and at > deadline:
            raise Busy()
        f.seek(0)
        f.truncate()
        f.write(repr(at))
    time.sleep(max(0.0, at - time.time()))


def request_for(name: str) -> urllib.request.Request:
    """The one request a lookup makes. The name rides in the POST body."""
    return urllib.request.Request(
        URL, data=urllib.parse.urlencode({"name": name}).encode("utf-8"))


def fetch(name: str, deadline: float | None = None) -> tuple[int, str]:
    """(HTTP status, body). An HTTP error status is an answer and comes back as
    one; a network failure raises OSError, because it is not an answer. Raises
    Busy, before sending anything, if the machine's next slot is past `deadline`."""
    _take_slot(deadline)
    try:
        with urllib.request.urlopen(request_for(name), timeout=TIMEOUT_S) as r:
            return r.status, r.read().decode("utf-8", errors="replace")
    except urllib.error.HTTPError as e:
        return e.code, e.read().decode("utf-8", errors="replace")


# ------------------------------------------------ does the name ask for stereo?
# A racemate is asked for by saying so, and a structure without stereo is the
# right answer to it. (+-), (RS), rac- and DL- match none of the patterns below,
# so nothing has to strip them first -- stripping them changes no answer.
# "rac-(1R,2S)" still asks for a RELATIVE configuration and is caught by its
# descriptors.
_DASHES = str.maketrans({"−": "-", "–": "-", "—": "-"})
_CIP = r"(?:\d+[a-z]?\s*)?(?:[RSEZ]\*?|Ra|Sa|Rp|Sp|P|M)"
_ASKS = (
    re.compile(r"\(\s*%s(?:\s*,\s*%s)*\s*\)" % (_CIP, _CIP)),  # (S) (1R,2S) (2E,4Z)
    re.compile(r"\(\s*[+-]\s*\)"),                               # (+) (-)
    re.compile(r"(?<![A-Za-z])[DL]-(?=[A-Za-z(\[])"),             # L-proline, D-glucose
    re.compile(r"(?<![A-Za-z])(?:cis|trans|meso|endo|exo)-", re.I),
)


def asks_for_stereo(name: str) -> bool:
    """Does the name specify a stereoisomer? Oxidation states -- palladium(0),
    copper(I) -- are not descriptors and do not match."""
    s = (name or "").translate(_DASHES)
    return any(p.search(s) for p in _ASKS)


# ------------------------------------------------------------------- metals
# Everything that is not a non-metal or a metalloid (B, Si, Ge, As, Sb, Te stay
# out: TBSCl is not an organometallic).
_NOT_METAL = frozenset({0, 1, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18,
                        32, 33, 34, 35, 36, 51, 52, 53, 54, 85, 86})
# ...except a lone, charged alkali or alkaline-earth ion. Na+ in sodium ascorbate,
# K+ in K2CO3, Li+ in LDA: for a salt, the ionic drawing is the right one.
_COUNTERION = frozenset({3, 4, 11, 12, 19, 20, 37, 38, 55, 56, 87, 88})


def metals_to_check(mol) -> list[str]:
    """Element symbols of the metals whose depiction PubChem may have mangled.
    A Pd2+ beside two acetates is NOT exempt: Pd(OAc)2 is written as text."""
    out = set()
    for a in mol.GetAtoms():
        z = a.GetAtomicNum()
        if z in _NOT_METAL:
            continue
        if z in _COUNTERION and a.GetDegree() == 0 and a.GetFormalCharge() != 0:
            continue
        out.add(a.GetSymbol())
    return sorted(out)


# ------------------------------------------------------------------ lookup
W_UNDEFINED = (
    "%d stereocenter(s) have no configuration in this record. PubChem does not "
    "know it either: if one stereoisomer is meant, look it up by a name that says "
    "which (\"(-)-menthol\", not \"menthol\"). Never add @/@@ by hand.")
W_ASKED = (
    "The name asks for a specific stereoisomer, but this record's SMILES carries "
    "no stereochemistry at all. Axial chirality (BINAP, BINOL) has no SMILES "
    "notation; otherwise PubChem matched a record without stereo -- read the "
    "title. A drawing from this SMILES will not show the configuration the name "
    "asks for.")
W_METAL = (
    "Contains %s. PubChem's SMILES often breaks metal-ligand bonds "
    "(n-butyllithium comes back as a lone Li atom beside a butyl radical, "
    "Pd(PPh3)4 as a Pd atom and four free PPh3). Look at the structure before "
    "drawing it; in a scheme such reagents are usually written as text.")
W_UNPARSED = ("RDKit cannot parse PubChem's SMILES, so the engine cannot draw "
              "this record.")
W_TITLE = ("PubChem took \"%s\" to mean \"%s\". Confirm that is the compound "
           "meant.")


def _folded(s: str) -> str:
    """The rule: the title and the name are compared literally, case
    and whitespace aside. Nothing cleverer -- "DMAP" against
    "4-(Dimethylamino)pyridine" is reported too, and that is the price of TMS
    against a tetramethoxystilbene never going unsaid."""
    return re.sub(r"\s+", "", s).casefold()


def _record(row: dict, query: str, asks: bool) -> dict:
    smi = row["SMILES"]
    cid = int(row["CID"])
    rec = {"cid": cid, "title": row.get("Title"), "smiles": smi,
           "inchikey": row.get("InChIKey"),
           "formula": row.get("MolecularFormula"), "url": PAGE % cid,
           "undefined_stereocenters": None, "warnings": []}
    title = rec["title"] or ""
    if _folded(title) != _folded(query):
        rec["warnings"].append(W_TITLE % (query, title or "a record with no title"))
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        rec["warnings"].append(W_UNPARSED)
        return rec
    # the same census stereo_lint reports, so the two tools cannot disagree
    unassigned = ST.ingest(smiles=[smi])[0]["unassigned"]
    rec["undefined_stereocenters"] = len(unassigned)
    if unassigned:
        rec["warnings"].append(W_UNDEFINED % len(unassigned))
    if asks and not any(c in smi for c in "@/\\"):
        rec["warnings"].append(W_ASKED)
    metals = metals_to_check(mol)
    if metals:
        rec["warnings"].append(W_METAL % ", ".join(metals))
    return rec


def lookup(name, deadline: float | None = None) -> dict:
    """One name or CAS number -> what PubChem holds under it. Never raises: a name
    that fails comes back as {"found": False, "error": ...} and carries no SMILES,
    so there is nothing to draw by mistake."""
    q = str(name if name is not None else "").strip()
    out = {"query": q, "found": False}
    if not q:
        out["error"] = "Empty name; nothing was sent."
        return out
    try:
        status, body = fetch(q, deadline)
    except Busy:
        out["busy"] = True
        out["error"] = ("Not looked up: other lookups on this server are using "
                        "PubChem's request allowance. Nothing was sent -- look "
                        "this name up again.")
        return out
    except OSError as e:
        out["unreachable"] = True
        out["error"] = ("PubChem could not be reached (%s: %s). No structure is "
                        "returned, and none should be written from memory in "
                        "its place." % (type(e).__name__, e))
        return out
    try:
        reply = json.loads(body)
    except ValueError:
        out["error"] = ("PubChem's reply (HTTP %d) was not JSON -- an error page "
                        "or a proxy, not an answer." % status)
        return out
    if status != 200:
        fault = reply.get("Fault", {}) if isinstance(reply, dict) else {}
        if fault.get("Code") == "PUGREST.NotFound":
            out["error"] = ("PubChem has no compound under this name. "
                            "Abbreviations often fail (DMP, TEA and PCC were not "
                            "found); try the full name or the CAS number.")
        else:
            out["error"] = "PubChem answered HTTP %d: %s" % (
                status, fault.get("Message") or fault.get("Code") or "no message")
        return out
    rows = (reply.get("PropertyTable") or {}).get("Properties") or []
    if not rows:
        out["error"] = "PubChem answered with no records."
        return out
    missing = [str(r.get("CID")) for r in rows if "SMILES" not in r]
    if missing:
        out["error"] = (
            "PubChem's reply has no SMILES field (CID %s). PubChem renamed this "
            "property once already (IsomericSMILES -> SMILES, 2025), and a stale "
            "name comes back as silence, not as an error. Please report this to the "
            "server administrator. ConnectivitySMILES is not used instead -- it "
            "carries no stereochemistry." % ", ".join(missing))
        return out

    asks = asks_for_stereo(q)
    out["found"] = True
    out["records"] = [_record(r, q, asks) for r in rows]
    out["warnings"] = []
    if len(rows) > 1:
        skeletons = {str(r.get("InChIKey", "")).split("-")[0] for r in rows}
        if len(skeletons) == 1:
            out["warnings"].append(
                "%d PubChem records carry this name and share one skeleton (the "
                "first InChIKey block): they differ in stereochemistry, isotopes "
                "or charge. Do not pick one silently -- say which is meant, or "
                "ask." % len(rows))
        else:
            out["warnings"].append(
                "%d PubChem records carry this name and they are DIFFERENT "
                "compounds (%d skeletons). Ask which one is meant."
                % (len(rows), len(skeletons)))
    return out


def lookup_all(names) -> list[dict]:
    """A list of names, one result each, in order. One that fails does not lose
    the others. More than MAX_NAMES is refused; names still unsent after
    BUDGET_S come back marked busy rather than holding the call past its limit."""
    names = list(names or [])
    if len(names) > MAX_NAMES:
        raise ValueError("%d names in one lookup; the limit is %d -- split them "
                         "into several lookups" % (len(names), MAX_NAMES))
    deadline = time.time() + BUDGET_S
    out = []
    for n in names:
        r = lookup(n, deadline)
        out.append(r)
        if r.get("busy"):
            deadline = 0.0          # the rest would only wait to be refused too
    return out
