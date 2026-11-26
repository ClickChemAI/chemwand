"""rxn_actions.py — a controlled vocabulary for experimental steps, and prose from it.

The vocabulary is IBM RXN's, from paragraph2actions:

    Copyright 2020 RXN for Chemistry team. MIT licence.
    https://github.com/rxn4chemistry/paragraph2actions
    (full licence text in THIRD_PARTY_NOTICES.txt)

Their action names and fields are distilled from about two million mined
procedures, which is why they are worth borrowing rather than inventing: the set
is what synthetic write-ups actually contain. Checked before use -- both
paragraph2actions and smiles2actions carry the verbatim MIT grant, "sell copies"
included, so a commercial product may use them with the notice retained.

WHAT IS NOT BORROWED, AND WHY. Their ActionSequenceValidator turned out to be a
TRAINING-DATA FILTER, not a procedure checker -- its own docstring says it
"defines which action sequences can be used for training". Three of its four
checks drop text-mining debris (InvalidAction, FollowOtherProcedure,
OtherLanguage) or reject sequences under five actions as uninformative, and one
counts "SLN" mentions against MakeSolution actions to police their own
representation. Nothing in it looks at temperature contradictions or addition
order, which a procedure checker would have to do. Only its multi-step detection
carries chemical meaning, and that one is reimplemented here.

Their four artefact types are therefore absent from ACTIONS: InvalidAction,
FollowOtherProcedure, OtherLanguage and NoAction exist because a pipeline was
reading paragraphs, and nothing here reads paragraphs.

THE MODEL NAMES CHEMICALS; IT DOES NOT COUNT THEM. paragraph2actions' Chemical
carries a `quantity`, and letting a model fill that would put its arithmetic back
into the procedure through a side door. Here a chemical is a NAME, and the amount
is substituted at render time from the stoichiometry calculator -- the same rule
that governs the charge sentence.
"""
from __future__ import annotations

# field -> required?  Fields not listed are rejected, so a model inventing
# `Add(pressure=...)` is told rather than silently ignored.
ACTIONS: dict = {
    "Add":            {"material": True, "dropwise": False, "temperature": False,
                       "atmosphere": False, "duration": False},
    "MakeSolution":   {"materials": True},
    "SetTemperature": {"temperature": True},
    "Stir":           {"duration": False, "temperature": False, "atmosphere": False},
    "Reflux":         {"duration": False, "dean_stark": False, "atmosphere": False},
    "Wait":           {"duration": True, "temperature": False},
    "Microwave":      {"duration": False, "temperature": False},
    "Sonicate":       {"duration": False, "temperature": False},
    # NOT IBM's, added here. Their vocabulary has Microwave and Sonicate and no
    # light source at all, so in stress tests every photochemical step (3 of 8)
    # smuggled "390 nm" into Stir's temperature field.
    "Irradiate":      {"wavelength": False, "light_source": False, "duration": False,
                       "temperature": False, "atmosphere": False},
    "Degas":          {"gas": False, "duration": False},
    "Quench":         {"material": True, "dropwise": False, "temperature": False},
    "PH":             {"material": False, "ph": False, "dropwise": False,
                       "temperature": False},
    "Extract":        {"solvent": True, "repetitions": False},
    "Wash":           {"material": True, "repetitions": False},
    "Partition":      {"material_1": True, "material_2": True},
    "PhaseSeparation": {},
    "CollectLayer":   {"layer": True},
    "Filter":         {"phase_to_keep": False},
    "DrySolution":    {"material": False},
    "DrySolid":       {"duration": False, "temperature": False, "atmosphere": False},
    "Concentrate":    {},
    "Purify":         {},
    "Recrystallize":  {"solvent": True},
    "Triturate":      {"solvent": True},
    "Yield":          {"material": True},
}

# Between two Yield actions, only these may appear -- anything else means the
# sequence has run two reactions together. This is the one check in IBM's
# validator with chemical meaning, kept.
_BETWEEN_YIELDS = {"Yield", "Purify", "Concentrate", "Recrystallize"}

_CHEMICAL_FIELDS = ("material", "solvent", "material_1", "material_2")


def validate(actions: list) -> dict:
    """Structural check. Returns {"ok", "problems": [...]}.

    NOT A CHEMISTRY CHECK. It catches a misspelled action, an invented field, a
    missing required one, and a sequence that is really two reactions. It says
    nothing about whether adding that reagent at that temperature is sensible --
    IBM's validator did not do that either, and rules for it would be chemistry
    nobody has specified."""
    problems = []
    if not actions:
        return {"ok": False, "problems": ["No actions: the action sequence is empty."]}

    for i, a in enumerate(actions, 1):
        if not isinstance(a, dict):
            problems.append("Action %d is not an object: %r" % (i, a))
            continue
        name = a.get("action")
        # THE KEY IS "action". Nearly every run in stress tests first wrote
        # {"type": "Add", ...} and was told "None is not in the vocabulary",
        # which names neither the mistake nor the fix.
        if name is None and "type" in a:
            problems.append(
                "Action %d uses the key 'type'; the key here is 'action', e.g. "
                "{\"action\": \"Add\", \"material\": {\"name\": \"K2CO3\"}}" % i)
            continue
        if name not in ACTIONS:
            problems.append(
                "Action %d: %r is not in the vocabulary. Available: %s"
                % (i, name, ", ".join(sorted(ACTIONS))))
            continue
        spec = ACTIONS[name]
        for k in a:
            if k != "action" and k not in spec:
                problems.append(
                    "Action %d (%s) has no field %r. It accepts: %s"
                    % (i, name, k, ", ".join(sorted(spec)) or "(no fields)"))
        for k, required in spec.items():
            if required and a.get(k) in (None, "", [], {}):
                problems.append("Action %d (%s) is missing %r" % (i, name, k))

    yields = [i for i, a in enumerate(actions)
              if isinstance(a, dict) and a.get("action") == "Yield"]
    if len(yields) > 1:
        between = [a.get("action") for a in actions[yields[0]:yields[-1] + 1]
                   if isinstance(a, dict)]
        if any(n not in _BETWEEN_YIELDS for n in between):
            problems.append(
                "There are %d Yield actions with %s between them -- this sequence "
                "describes more than one reaction step; write each step separately."
                % (len(yields),
                   ", ".join(sorted(set(n for n in between if n not in _BETWEEN_YIELDS)))))
    return {"ok": not problems, "problems": problems}


def chemicals_used(actions: list) -> list:
    """Every chemical NAME the sequence mentions, in order of first appearance.

    The caller needs this to check the sequence against the reagent table: a
    workup that mentions a reagent nobody weighed out is worth flagging, and one
    that never mentions a reagent that was weighed out probably forgot it."""
    seen, out = set(), []
    for a in actions or []:
        if not isinstance(a, dict):
            continue
        vals = [a.get(f) for f in _CHEMICAL_FIELDS]
        vals += list(a.get("materials") or [])
        for v in vals:
            nm = v.get("name") if isinstance(v, dict) else v
            if isinstance(nm, str) and nm.strip() and nm not in seen:
                seen.add(nm)
                out.append(nm)
    return out


def charged_names(actions: list) -> list:
    """The chemicals the sequence itself PUTS IN -- Add, MakeSolution, Quench, PH.

    procedure.render leaves these out of its opening charge sentence, so each
    reagent is written once, where it goes in. Workup chemicals (Extract, Wash,
    ...) are not here: they wash things out, they are not part of the charge."""
    out = []
    for a in actions or []:
        if not isinstance(a, dict):
            continue
        if a.get("action") in ("Add", "Quench", "PH"):
            vals = [a.get("material")]
        elif a.get("action") == "MakeSolution":
            vals = list(a.get("materials") or [])
        else:
            continue
        for v in vals:
            nm = v.get("name") if isinstance(v, dict) else v
            if isinstance(nm, str) and nm.strip() and nm.strip() not in out:
                out.append(nm.strip())
    return out


# ============================== prose ========================================
#
# ACS supporting-information voice: passive, past tense, one clause per action.
#
# `Purify` HAS NO FIELDS IN THEIR VOCABULARY, which is a real gap rather than an
# omission in this module -- "purified" with no method is all their mined
# representation carried. It renders as a bare clause and the caller is told, so
# the chemist can say how.

def _amt(name, amounts):
    """The computed quantity for a named chemical, or nothing.

    A chemical the reagent table never mentioned gets NO invented amount -- it is
    reported instead, because a workup that washes with a reagent nobody weighed
    out is either a missing table row or a hallucinated reagent."""
    got = (amounts or {}).get(name)
    return " (%s)" % got if got else ""


def _bits(*pairs):
    return "".join(p for p in pairs if p)


def render(actions: list, amounts: dict = None, labels: dict = None) -> dict:
    """Actions -> sentences, with quantities substituted from the calculator.

    labels -- {any name a chemical goes by: the name to print}. Optional; without
              it a chemical prints as the action named it."""
    amounts = amounts or {}
    labels = labels or {}
    out, unknown = [], []

    def chem(v, with_amount=True):
        nm = v.get("name") if isinstance(v, dict) else v
        nm = str(nm or "").strip()
        if nm and nm not in amounts and nm not in unknown:
            unknown.append(nm)
        # A WORKUP SOLVENT CARRIES NO AMOUNT, EVEN WHEN IT SHARES A NAME WITH THE
        # REACTION SOLVENT. Extract(CH2Cl2) picked up the reaction's computed
        # 50.0 mL -- "extracted twice with CH2Cl2 (50.0 mL)", where 30 mL each
        # was used (seen in a stress test). No workup volume is ever computed
        # here; the calculator's number belongs to the charge alone.
        return labels.get(nm, nm) + (_amt(nm, amounts) if with_amount else "")

    for a in actions or []:
        n = a.get("action")
        t = a.get("temperature")
        d = a.get("duration")
        atm = a.get("atmosphere")
        at = " at %s" % t if t else ""
        during = " for %s" % d if d else ""
        under = " under %s" % atm if atm else ""

        if n == "Add":
            out.append("%s was added%s%s%s%s." % (
                chem(a["material"]), " dropwise" if a.get("dropwise") else "",
                at, " over %s" % d if d else "", under))
        elif n == "MakeSolution":
            mats = [chem(m) for m in a.get("materials") or []]
            out.append("A solution of %s was prepared." % (
                " in ".join(mats) if len(mats) == 2 else ", ".join(mats)))
        elif n == "SetTemperature":
            out.append("The mixture was brought to %s." % a["temperature"])
        elif n == "Stir":
            out.append("The mixture was stirred%s%s%s." % (at, during, under))
        elif n == "Reflux":
            out.append("The mixture was heated at reflux%s%s%s." % (
                during, " with a Dean-Stark trap" if a.get("dean_stark") else "",
                under))
        elif n == "Wait":
            out.append("The mixture was left%s for %s." % (at, a["duration"]))
        elif n == "Microwave":
            out.append("The mixture was irradiated in a microwave reactor%s%s."
                       % (at, during))
        elif n == "Sonicate":
            out.append("The mixture was sonicated%s%s." % (at, during))
        elif n == "Irradiate":
            wl = str(a.get("wavelength") or "").strip()
            src = str(a.get("light_source") or "").strip()
            light = (("%s %s" % (wl, src)) if wl and src
                     else src or ("%s light" % wl if wl else ""))
            out.append("The mixture was irradiated%s%s%s%s." % (
                " with %s" % light if light else "", at, during, under))
        elif n == "Degas":
            out.append("The mixture was degassed%s%s." % (
                " with %s" % a["gas"] if a.get("gas") else "", during))
        elif n == "Quench":
            out.append("The reaction was quenched%s by%s addition of %s." % (
                at, " dropwise" if a.get("dropwise") else "", chem(a["material"])))
        elif n == "PH":
            with_what = " with %s" % chem(a["material"]) if a.get("material") else ""
            out.append("The pH was adjusted%s%s%s%s." % (
                " to %s" % a["ph"] if a.get("ph") else "", at, with_what,
                " dropwise" if a.get("dropwise") else ""))
        elif n == "Extract":
            out.append("The mixture was extracted%s with %s." % (
                _reps(a.get("repetitions")), chem(a["solvent"], False)))
        elif n == "Wash":
            out.append("The organic layer was washed%s with %s." % (
                _reps(a.get("repetitions")), chem(a["material"], False)))
        elif n == "Partition":
            out.append("The residue was partitioned between %s and %s."
                       % (chem(a["material_1"], False), chem(a["material_2"], False)))
        elif n == "PhaseSeparation":
            out.append("The layers were separated.")
        elif n == "CollectLayer":
            out.append("The %s layer was collected." % a["layer"])
        elif n == "Filter":
            keep = a.get("phase_to_keep")
            out.append("The mixture was filtered%s."
                       % (", keeping the %s" % keep if keep else ""))
        elif n == "DrySolution":
            out.append("The organic layer was dried%s."
                       % (" over %s" % chem(a["material"], False) if a.get("material") else ""))
        elif n == "DrySolid":
            out.append("The solid was dried%s%s%s." % (at, during, under))
        elif n == "Concentrate":
            out.append("The solvent was removed under reduced pressure.")
        elif n == "Purify":
            out.append("The residue was purified.")
        elif n == "Recrystallize":
            out.append("The solid was recrystallised from %s."
                       % chem(a["solvent"], False))
        elif n == "Triturate":
            out.append("The solid was triturated with %s." % chem(a["solvent"], False))
        elif n == "Yield":
            # Attaches to the sentence before it, as in a real write-up.
            tail = "to afford %s" % chem(a["material"])
            if out:
                out[-1] = out[-1].rstrip(".") + " " + tail + "."
            else:
                out.append(tail.capitalize() + ".")

    res = {"sentences": out, "text": " ".join(out)}
    if unknown:
        res["chemicals_not_in_the_table"] = unknown
        # PHRASED SO IT IS NOT A FALSE ALARM. Calling every unmatched name
        # suspicious would flag the extraction EtOAc and the wash water --
        # which belong in nobody's stoichiometry table. A
        # workup solvent being absent is normal; a REAGENT being absent is not,
        # and only the reader can tell which this is.
        res["note"] = ("These names are not in the reagent table, so they have no amounts: %s. "
                       "Workup solvents are never in the amounts table, so that is normal for "
                       "them; but if any of these should have been weighed out, a row is "
                       "missing from the table."
                       % ", ".join(unknown))
    if any(a.get("action") == "Purify" for a in actions or []):
        res["purify_is_vague"] = ("Purify has no fields in this vocabulary -- it can only say "
                                  "\"purified\", not how. Add the method yourself (e.g. flash "
                                  "chromatography on silica gel, eluent).")
    return res


# AN ADVERB, NOT "(2 x)". An SI writes "washed with water (2 x 10 mL)" -- the
# repetition and the volume together. There is no volume here (the model supplies
# the count and nobody weighed the wash water), so "(2 x)" comes out as a dangling
# multiplier: "washed with water (2 x)". In English the count goes before the
# solvent as an adverb, and reads like a procedure.
_ADVERBS = {2: " twice", 3: " three times", 4: " four times", 5: " five times"}


def _reps(n):
    try:
        n = int(n)
    except (TypeError, ValueError):
        return ""
    if not n or n < 2:
        return ""
    return _ADVERBS.get(n, " %d times" % n)
