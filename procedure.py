"""procedure.py — an ACS-style experimental paragraph, with the numbers put in by code.

The design's first principle, and the reason this file exists at all:

    Numbers are computed by code, prose is written by Claude, and the two never mix.

A language model's arithmetic is wrong in exactly the same voice as its arithmetic
is right, and "weigh 1.38 g" in someone's supporting information is not a place to
be approximately right. So THE MODEL IS NEVER SHOWN A MASS. It supplies the shape
of the paragraph -- what goes in first, at what temperature, how it was worked up
-- and every quantity is substituted here from stoichiometry.compute(). A wrong
amount is then not unlikely; it is structurally impossible.

WHAT IS STILL THE MODEL'S, AND SAID SO. Addition order, temperature, time, workup
and purification are judgement, and by design they are generated rather than
typed. The workup sentence in particular carries numbers the model
DID write -- "washed with 3 x 10 mL water" -- so the return value separates the
computed quantities from the free text rather than handing back one paragraph that
looks uniformly authoritative.

EQUIVALENTS, NOT MOL%, UNLESS ASKED. Catalysts are conventionally reported in
mol%, and the practical boundary sits somewhere near 0.2 equiv -- which is exactly
the kind of number this project does not invent. `catalyst: true` on a reagent
turns on mol%; nothing is inferred from its size.
"""
from __future__ import annotations

import re

import stoichiometry as STOICH

# Actions that state the reaction conditions themselves; with one of these in
# the sequence the generated "The mixture was stirred at ..." would say it twice.
_CONDITION_ACTIONS = {"Stir", "Reflux", "Wait", "Microwave", "Sonicate",
                      "Irradiate", "SetTemperature"}
# Actions that are still part of charging the flask; the reaction sentence goes
# after the run of these at the start of the sequence.
_CHARGE_ACTIONS = {"Add", "MakeSolution", "Degas"}


def _sentence(s: str) -> str:
    """Capitalise the first letter -- but not "pH" or "n-BuLi"."""
    if len(s) > 1 and s[0].islower() and s[1].islower():
        return s[0].upper() + s[1:]
    return s


def _purified_by(p) -> str:
    """The "how" of a purification, from whatever the caller wrote:
    "flash chromatography (...)", "purified by flash ...", or a full sentence."""
    p = str(p or "").strip().rstrip(".。")
    return re.sub(r"(?i)^(the (residue|crude( product)?) was )?purified (by|via|using) ",
                  "", p)

# WHO WROTE THE PROSE IS A PARAMETER, BECAUSE IT DIFFERS BY CALLER AND GETTING IT
# WRONG IS THE EXACT ERROR THIS MODULE EXISTS TO PREVENT. Through the MCP tool
# Claude supplies the conditions; on the /stoich page the chemist types them into a
# form. One note saying "written by the model" in both cases would tell a chemist
# that their own typed temperature had been invented by a model.
#
# TWO SENTENCES, NOT ONE WITH A WORD SWAPPED. Substituting the subject gives
# "the amounts were computed by the calculator; the person filling in the form
# took no part" -- true of the arithmetic and absurd as a
# whole, since the equivalents and molecular weights are exactly what the chemist
# did supply. The clause that matters differs by caller, so it is written twice.
_NOTES = {
    # "WRITTEN BY THE MODEL" OVERSTATED IT. Through the tool Claude passes these
    # on, often straight from what the chemist said, and nothing here can tell
    # dictated from proposed -- so the note says both and asks for a check. And
    # the product's mass and yield sit in parentheses too but are measured, not
    # computed.
    "model": ("The amounts in parentheses were computed by the calculator -- except "
              "the product's mass and yield, which are the measured values given. "
              "The addition order, temperatures, times, and the workup and "
              "purification wording were filled in by the model: they may be what "
              "you told it or its own proposal, so check them -- including any "
              "volumes in the workup sentences (e.g. \"washed with 3 × 10 mL\"), "
              "which were not calculated."),
    "user": ("The amounts in parentheses were computed by the calculator from the "
             "equivalents and molecular weights you entered (the product's mass and "
             "yield are your own values); no model was involved. "
             "The addition order, temperatures, times, workup and purification were "
             "entered by you -- including any volumes in the workup sentences (e.g. "
             "\"washed with 3 × 10 mL\"), which were not calculated."),
}

# ACS supporting-information house style, as used in J. Org. Chem. / Org. Lett.:
# the charge sentence names each species with (mass, mmol, equiv) in that order,
# solvents get a volume only, and solutions state their concentration first.
_JOIN_LAST = " and "


def _equiv_txt(eq: float) -> str:
    """Two decimals where the second one is real, one where it is not.

    "%.1f" PRINTS 0.25 equiv AS "0.2" AND 2.25 AS "2.2". 0.25 is exact in
    binary, so the halfway case rounds to even and a dictated amount comes back
    different from what was weighed. A
    trailing zero is dropped so 2.0 and 1.5 still read the way an SI writes
    them."""
    s = "%.2f" % eq
    return s[:-1] if s.endswith("0") else s


def _fmt_amount(r: dict, limiting: bool) -> str:
    """The parenthetical after a reagent name. Every number here is computed."""
    bits = []
    if r["form"] == "solution" and r.get("molarity"):
        bits.append("%g M" % r["molarity"])
    if r.get("display"):
        bits.append(r["display"])
    if r.get("mmol") is not None:
        bits.append("%s mmol" % STOICH.sig(r["mmol"]))
    if r.get("equiv") is not None:
        if r.get("catalyst"):
            bits.append("%.3g mol%%" % (r["equiv"] * 100))
        elif limiting:
            bits.append("1.0 equiv")
        else:
            # "2.0 equiv", not "2" -- ACS style keeps the decimal, and a bare
            # integer beside "1.0 equiv" on the limiting reagent reads as sloppy.
            eq = r["equiv"]
            bits.append(("%s equiv" % _equiv_txt(eq)) if eq >= 0.1
                        else ("%.3g equiv" % eq))
    if r.get("purity"):
        bits.append("%.0f wt%%" % (r["purity"] * 100))
    return ", ".join(bits)


def _name(r: dict, names: dict) -> str:
    """A reagent's written name, with its compound number if it has one."""
    nm = names.get(r["name"], r["name"])
    num = (r.get("compound_number") or "").strip() if isinstance(
        r.get("compound_number"), str) else ""
    return "%s (**%s**" % (nm, num) if num else nm


def _series(parts: list) -> str:
    if not parts:
        return ""
    if len(parts) == 1:
        return parts[0]
    return ", ".join(parts[:-1]) + _JOIN_LAST + parts[-1]


def render(amounts: dict, *, limiting: str = None, order: list = None,
           names: dict = None, temperature: str = None, time: str = None,
           addition_temperature: str = None,
           atmosphere: str = None, workup: str = None,
           purification: str = None, product: dict = None,
           catalysts: list = None, prose_author: str = "model",
           actions: list = None) -> dict:
    """Write the paragraph.

    amounts     -- the dict stoichiometry.compute() returned. Not recomputed here.
    limiting    -- name of the limiting reagent. REQUIRED, never guessed,
                   because it is not always the substrate.
    order       -- reagent names in addition order. The model's call.
    names       -- {reagent name: written name}, e.g. {"K2CO3": "potassium carbonate"}
    catalysts   -- names (or {"name": ...}) to report in mol% instead of equiv; one
                   that is not a row comes back in `catalysts_not_in_the_table`.
    actions     -- a row an Add / MakeSolution / Quench / PH action names is left
                   out of the opening charge and written where the action adds it.
    product     -- {"name", "mass_mg"?, "yield_pct"?} for the closing sentence.
    """
    if not amounts or not amounts.get("reagents"):
        return {"ok": False, "error": "No reagent amounts to work from: `reagents` is empty or nothing could be computed."}
    if amounts.get("ask"):
        return {"ok": False, "error": "The reagent data is incomplete, so no procedure was written. Answer the questions in `ask` first.",
                "ask": list(amounts["ask"])}
    if not limiting:
        return {"ok": False,
                "error": "`limiting` is required: name the limiting reagent. It is usually the "
                         "substrate, but an expensive coupling partner or a scarce ligand can be "
                         "limiting, so the first reagent is never assumed."}

    prose_author = str(prose_author or "model")
    if prose_author not in _NOTES:
        prose_author = "model"
    names = dict(names or {})
    # A CHEMICAL MAY BE NAMED BY ITS KEY OR BY ITS WRITTEN NAME -- in `catalysts`
    # as in an action. One map resolves the written name back to the key.
    key_of = {v: k for k, v in names.items()}
    # NAMES OR {"name": ...}, BOTH. Catalysts passed the way chemicals are passed
    # everywhere else, [{"name": "Pd"}], would otherwise crash the whole call with
    # "unhashable type: 'dict'".
    cats = []
    for c in catalysts or []:
        nm = c.get("name") if isinstance(c, dict) else c
        if isinstance(nm, str) and nm.strip():
            cats.append(nm.strip())
    cat_keys = {key_of.get(c, c) for c in cats}
    rows = {r["name"]: dict(r, catalyst=(r["name"] in cat_keys))
            for r in amounts["reagents"]}
    # A CATALYST THAT IS NOT A ROW MUST NOT BE DROPPED WITHOUT A WORD. "[Ni]" and
    # "BP8" have no molecular weight, so they never become rows, and the paragraph
    # would leave them out while the reply said nothing. So they are reported.
    cats_missing = [c for c in cats if key_of.get(c, c) not in rows]
    if limiting not in rows:
        return {"ok": False,
                "error": "The limiting reagent %r is not in the reagent table. The table has: %s"
                         % (limiting, ", ".join(rows))}

    seq = [n for n in (order or list(rows)) if n in rows]
    for n in rows:                          # anything the order forgot still goes in
        if n not in seq:
            seq.append(n)

    # A ROW AN ACTION PUTS IN IS NOT ALSO IN THE OPENING CHARGE. A charge
    # sentence listing every row puts a catalyst added after degassing, or a base
    # added dropwise at 0 °C, in twice -- once in "was added ..." and again in
    # the action -- while the reply still says ok: true. Seen with Boc2O, n-BuLi,
    # a quench HCl and a Pd catalyst. Each row appears once, where it actually
    # goes in.
    later = set()
    if actions:
        import rxn_actions as RXA
        for nm in RXA.charged_names(actions):
            k = nm if nm in rows else key_of.get(nm)
            if k in rows:
                later.add(k)

    solvents = [rows[n] for n in seq
                if rows[n]["form"] == "solvent" and n not in later]
    others = [rows[n] for n in seq
              if rows[n]["form"] != "solvent" and n not in later]
    # THE LIMITING REAGENT LEADS THE CHARGE UNLESS AN ACTION ADDS IT -- an ester
    # dropped into freshly made LDA, say. Then the charge opens with whatever is
    # in the flask first; if nothing is, there is no charge sentence at all.
    if limiting not in later:
        first = rows[limiting]
    else:
        first = others[0] if others else None
    rest = [r for r in others if first is None or r["name"] != first["name"]]

    computed = []
    def _cite(r):
        txt = _fmt_amount(r, r["name"] == limiting)
        # THE WRITTEN NAME, WHERE THERE IS ONE. Otherwise the amounts table shows
        # "O=Cc1ccc(Br)cc1  385 mg" because the /retro handoff fills the reagent
        # key with the SMILES -- readable in the prose, which uses `names`, and
        # unreadable in the table beside it.
        computed.append({"reagent": names.get(r["name"], r["name"]),
                         "key": r["name"], "amount": txt})
        nm = _name(r, names)
        close = ", %s)" % txt if nm.endswith("**") else " (%s)" % txt
        return nm + close

    # 1. the charge sentence
    # A SOLVENT NAME MAY ALREADY END IN A PARENTHESIS. "toluene/EtOH (4:1)" plus
    # " (10.4 mL)" gives "toluene/EtOH (4:1) (10.4 mL)"; an SI writes the ratio and
    # the volume in one bracket.
    def _solv(r):
        nm = names.get(r["name"], r["name"]).strip()
        if nm.endswith(")"):
            return "%s, %s)" % (nm[:-1].rstrip(), r["display"])
        return "%s (%s)" % (nm, r["display"])

    solv_txt = _series([_solv(r) for r in solvents])
    for r in solvents:
        computed.append({"reagent": names.get(r["name"], r["name"]),
                         "key": r["name"], "amount": r["display"]})
    # WITH NOTHING TO ADD, "To a solution of X in Y" IS NOT A SENTENCE. It has no
    # verb -- a step whose only non-solvent reagent IS the limiting one would
    # come out as "To a solution of 4-bromobenzaldehyde (385 mg, 2.08 mmol,
    # 1.0 equiv) in toluene (10.4 mL)." and stop there. An SI writes that case
    # the other way round: the substrate was dissolved.
    dissolve_only = not rest
    if first is None:
        lead = None
    elif dissolve_only:
        # "was dissolved" needs something to dissolve INTO. With no solvent row
        # the sentence would come out as "4-bromobenzaldehyde (385 mg...) was
        # dissolved." -- in what? A neat charge is "charged".
        if solv_txt:
            lead = "%s was dissolved in %s" % (_cite(first), solv_txt)
        else:
            lead = "%s was charged" % _cite(first)
    else:
        lead = ("To a solution of %s" % _cite(first) if solvents
                else "To %s" % _cite(first))
        if solv_txt:
            lead += " in %s" % solv_txt
    # THE ADDITION TEMPERATURE IS NOT THE REACTION TEMPERATURE. Putting the one
    # `temperature` in both sentences produces "To a solution ... at 85 °C was
    # added ... The mixture was stirred at 85 °C", which says the same thing twice
    # and is usually wrong besides: most charges go in cold or at room temperature
    # and are heated afterwards. Only an explicitly given addition temperature
    # appears here.
    if first is None:
        # EVERYTHING WEIGHED GOES IN BY AN ACTION. A solvent that no action names
        # still has a computed volume, and dropping it would lose the one number
        # that says how dilute the reaction was.
        sentences = ["The reaction was run in %s." % solv_txt] if solv_txt else []
    else:
        if addition_temperature:
            lead += " at %s" % addition_temperature
        sentences = [
            (lead + " was added " + _series([_cite(r) for r in rest]) + ".")
            if rest else (lead + ".")]

    # THE ROWS AN ACTION ADDS STILL NEED THEIR AMOUNTS -- in the map the actions
    # read below, and in the table shown beside the paragraph. Leaving them out
    # of the charge sentence must not leave them out of either.
    for n in seq:
        if n in later:
            r = rows[n]
            txt = (r["display"] if r["form"] == "solvent"
                   else _fmt_amount(r, n == limiting))
            computed.append({"reagent": names.get(n, n), "key": n,
                             "amount": txt})

    # 2. THE REACTION SENTENCE, UNLESS THE ACTIONS COVER IT. An action sequence
    # carries its own Stir, so generating one here as well produces "The mixture
    # was stirred. The mixture was stirred at 85 °C for 6 h under N2." -- and the
    # first of those is a bare sentence saying nothing, because with actions in
    # play nobody passes `temperature` or `time` separately.
    # BUT ONLY IF THEY DO COVER IT. A sequence that starts at the workup has no
    # Stir in it, and the 90 °C / 12 h / N2 the caller gave vanished: the
    # paragraph went straight from the charge to the extraction.
    covered = bool(actions) and any(
        isinstance(a, dict) and a.get("action") in _CONDITION_ACTIONS for a in actions)
    stir = None
    if not actions or (not covered and (temperature or time or atmosphere)):
        stir = "The mixture was stirred"
        if temperature:
            stir += " at %s" % temperature
        if time:
            stir += " for %s" % time
        if atmosphere:
            stir += " under %s" % atmosphere
        if not actions:
            sentences.append(stir + ".")

    # 3-4. WHAT HAPPENS AFTER THE REAGENTS ARE IN, one of two ways.
    #
    # `actions` is the structured path: IBM RXN's controlled vocabulary, checked
    # for shape and rendered to prose here, with quantities substituted from the
    # calculator exactly as the charge sentence does. It exists because a model
    # writing free prose cannot be checked at all, while a model choosing from 24
    # action types can be told it misspelled one or invented a field.
    #
    # `workup`/`purification` stay for the /stoich page, where a chemist types
    # their own sentences and there is nothing to check -- they wrote it.
    #
    # THE ACTIONS START AFTER THE CHARGE. The charge sentence is generated from
    # the amounts, so a sequence that begins by adding the reagents again would
    # say it twice. Add actions are still allowed -- a reagent added later, or a
    # quench, is exactly that -- but the opening charge is not the sequence's job.
    action_report = None
    if actions:
        import rxn_actions as RXA

        checked = RXA.validate(actions)
        if not checked["ok"]:
            return {"ok": False,
                    "error": "The action sequence is invalid, so no procedure was written; see `problems`.",
                    "problems": checked["problems"]}
        # BOTH NAMES MAP TO THE AMOUNT. `computed` reports the WRITTEN name so
        # the table beside the paragraph is readable, but an action may name a
        # chemical either way -- and keying only on one of them produces
        # "K2CO3 was added" with no quantity, plus a complaint that K2CO3 was not
        # in the reagent table when it was, under "potassium carbonate".
        amt_map = {}
        for c in computed:
            amt_map[c["reagent"]] = c["amount"]
            if c.get("key"):
                amt_map[c["key"]] = c["amount"]
        # THE PRODUCT'S OWN NUMBERS HAVE TO REACH THE Yield ACTION. Its material
        # is the product, which is not a row in the reagent table, so without
        # this the paragraph ends "...to afford the biaryl." -- the mass and the
        # yield silently dropped, which is the one number a reader looks for.
        yielded = next((a.get("material") for a in actions
                        if isinstance(a, dict) and a.get("action") == "Yield"), None)
        if yielded is not None and product:
            key = (yielded.get("name") if isinstance(yielded, dict)
                   else str(yielded)).strip()
            pbits = []
            if product.get("mass_mg") is not None:
                mg = float(product["mass_mg"])
                pbits.append(STOICH._fmt_mass(mg))
            if product.get("yield_pct") is not None:
                pbits.append("%.0f%%" % float(product["yield_pct"]))
            if key and pbits:
                amt_map[key] = ", ".join(pbits)
                computed.append({"reagent": key, "amount": ", ".join(pbits)})
        # THE WRITTEN NAME IN THE PROSE, WHICHEVER NAME THE ACTION USED. An action
        # naming the key would print the key -- "alcohol (3.00 g ...) was added"
        # -- while the charge sentence beside it uses `names`. Both print what the
        # chemist would write.
        labels = {}
        for c in computed:
            if c.get("key"):
                labels[c["key"]] = c["reagent"]
                labels[c["reagent"]] = c["reagent"]
        action_report = RXA.render(actions, amt_map, labels=labels)
        free_parts = [x.rstrip(".") for x in action_report["sentences"]]
        # THE REACTION GOES IN AFTER THE LAST THING ADDED, NOT AFTER THE CHARGE.
        # Put straight after the opening sentence it read "stirred at 90 °C for
        # 12 h. The mixture was degassed. Pd(PPh3)4 was added." -- heated before
        # the catalyst was in. Each of these actions renders one sentence.
        if stir:
            k = 0
            for a in actions:
                if not (isinstance(a, dict) and a.get("action") in _CHARGE_ACTIONS):
                    break
                k += 1
            free_parts.insert(k, stir)
        # A Purify action has no fields, so the dictated column conditions in
        # `purification` were dropped and the paragraph said only "purified".
        how = _purified_by(purification)
        if how:
            i = next((k for k, s in enumerate(free_parts)
                      if s.startswith("The residue was purified")), None)
            if i is not None:
                free_parts[i] = free_parts[i].replace(
                    "The residue was purified", "The residue was purified by " + how, 1)
            else:
                free_parts.append("The residue was purified by " + how)
    else:
        # The caller's own sentences, made into sentences: "the mixture was
        # diluted..." began lowercase after a full stop, and a bare
        # "flash chromatography (...)" had no verb at all.
        free_parts = [_sentence(str(workup).strip().rstrip(".。"))
                      ] if workup and str(workup).strip() else []
        how = _purified_by(purification)
        if how:
            p = str(purification).strip().rstrip(".。")
            free_parts.append(_sentence(p) if re.search(r"\b(was|were)\b", p)
                              else "The residue was purified by " + how)

    # 5. THE YIELD CLAUSE ATTACHES, IT DOES NOT STAND ALONE. Appending it as its
    # own sentence produces "...chromatography (hexanes/EtOAc 9:1). to afford the
    # product (340 mg, 83%)." -- a fragment starting with a lowercase "to" after a
    # full stop. In a real SI the yield hangs off the purification clause, so it
    # goes onto the last sentence there is.
    # A Yield ACTION HAS ALREADY ATTACHED THE PRODUCT. Appending `product` on top
    # would name it twice in one sentence.
    tail = ""
    if actions and any(a.get("action") == "Yield" for a in actions
                       if isinstance(a, dict)):
        product = None
    if product:
        bits = []
        if product.get("mass_mg") is not None:
            mg = float(product["mass_mg"])
            bits.append(STOICH._fmt_mass(mg))
        if product.get("yield_pct") is not None:
            bits.append("%.0f%%" % float(product["yield_pct"]))
        tail = "to afford %s" % (product.get("name") or "the product")
        if bits:
            tail += " (%s)" % ", ".join(bits)
            computed.append({"reagent": product.get("name") or "product",
                             "amount": ", ".join(bits)})

    if free_parts and tail:
        free_parts[-1] = free_parts[-1] + " " + tail
    elif tail:
        sentences[-1] = sentences[-1].rstrip(".") + " " + tail
    sentences.extend(f + "." for f in free_parts)

    out_extra = {}
    if cats_missing:
        out_extra["catalysts_not_in_the_table"] = cats_missing
        out_extra["catalysts_note"] = (
            "These catalysts are not in the reagent table and are missing from the "
            "paragraph: %s. Usually this means no molecular weight was available -- "
            "give a structure (smiles) or an mw so they can be added to the table "
            "and their amounts computed." % ", ".join(cats_missing))
    if action_report:
        for k in ("chemicals_not_in_the_table", "purify_is_vague"):
            if action_report.get(k):
                out_extra[k] = action_report[k]
        if action_report.get("note"):
            out_extra["actions_note"] = action_report["note"]
    return {
        "ok": True,
        "text": " ".join(sentences),
        **out_extra,
        "scale_mmol": amounts.get("scale_mmol"),
        "limiting": limiting,
        "computed_amounts": computed,
        "model_wrote": [k for k, v in (("order", order), ("temperature", temperature),
                                       ("time", time), ("atmosphere", atmosphere),
                                       ("workup", workup),
                                       ("purification", purification)) if v],
        "prose_author": prose_author,
        "note": _NOTES.get(prose_author, _NOTES["model"]),
        "from_actions": bool(actions),
    }
