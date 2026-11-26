"""stoichiometry.py — how much of each thing to weigh out, and nothing else.

Pure arithmetic. No RDKit, no network, no model. Molecular weights arrive as
numbers; resolving them from a SMILES is the caller's job, which keeps this file
testable against hand-checked values and keeps chemistry perception out of the one
place that must never be approximately right.

WHY THERE IS NO REAGENT TABLE. One was considered -- a hundred common reagents
with densities and molarities filled in once -- and rejected for a better reason
than any in its favour: concentrations are not fixed. A table that says n-BuLi
is 2.5 M, next to a bottle that says 2.2 M, does not save the user a lookup; it
computes a wrong number with the same confidence as a right one, and nothing in the
output says which happened. The only correct source for what is in someone's bottle
is that someone.

SO THE ERRORS ARE THE PRODUCT. `ask` comes back with one sentence per missing
number, phrased as a question. A caller that cannot compute an amount is not
supposed to fall back on a plausible default -- it is supposed to ask. That is the
whole reason this returns {"ok": False, "ask": [...]} instead of raising or
guessing.

WHAT THE GUESS WOULD COST, measured, because "be careful with units" is not an
argument: n-BuLi at 1.2 equiv on 5 mmol is 2.40 mL as its usual 2.5 M solution and
0.57 mL treated as a neat liquid at d 0.68. A factor of 4.2. Someone follows the
second number and the reaction is dead, and the procedure reads perfectly.
"""
from __future__ import annotations

import math

# Every form, and exactly what it needs. A form that is not in here is refused by
# name rather than silently treated as a solid.
REQUIRED = {
    "solid":    ("mw",),                 # -> mass
    "neat":     ("mw", "density"),       # -> mass and volume
    "solution": ("molarity",),           # -> volume; mw is not needed and not used
    "solvent":  ("concentration",),      # -> volume, from the limiting reagent's mmol
}

_LABEL = {"mw": "molecular weight (mw, or the reagent's smiles)", "density": "density (g/mL)", "molarity": "molarity (mol/L)",
          "concentration": "reaction concentration (mol/L, relative to the limiting reagent)"}

# BOTH NAMES, BECAUSE `ask` HAS TWO READERS. A chemist reads these strings
# directly, where a bare "neat" means nothing; Claude meets them through the MCP
# tool and needs the key it must actually put in the call. Printing both "neat
# liquid" and "neat" costs one word and serves both, instead of picking an
# audience to fail.
_FORM = {"solid": "solid", "neat": "neat liquid",
         "solution": "solution", "solvent": "solvent"}


SIG = 3            # significant figures on every amount printed (ACS SI style)
_SIG_DATA = 4      # one more in the returned data, so nothing is lost there


def _round_sig(x: float, n: int = _SIG_DATA) -> float:
    """x to n significant figures, as a number."""
    return float("%.*g" % (n, x)) if x else 0.0


def sig(x: float, n: int = SIG) -> str:
    """x to n significant figures, as written in an SI: trailing zeros kept,
    never in exponent form. 2.4 -> "2.40", 0.163 -> "0.163", 333.3 -> "333".

    FIXED DECIMAL PLACES WERE THE BUG. "%.1f mg" printed 0.163 mg of a catalyst
    as 0.2 mg (+23%), "%.0f µL" printed 0.34 µL of Et3N as 0 µL, and "%.2f M"
    printed 0.005 M as 0.01 M. Precision belongs to the number, not its band."""
    x = _round_sig(x, n)         # round first: 9.996 is "10.0", not "10.00"
    if x == 0:
        return "0"
    places = n - 1 - math.floor(math.log10(abs(x)))
    return "%.*f" % (places, x) if places > 0 else "%d" % round(x)


def _fmt_mass(mg: float) -> str:
    """A procedure says 1.38 g, not 1382 mg. Below a gram it says mg."""
    mg = _round_sig(mg, SIG)
    if mg >= 1000:
        return "%s g" % sig(mg / 1000)
    return "%s mg" % sig(mg)


def _fmt_vol(ml: float) -> str:
    ml = _round_sig(ml, SIG)
    if ml >= 1:
        return "%s mL" % sig(ml)
    return "%s µL" % sig(ml * 1000)


def scale_to_mmol(scale: dict) -> tuple[float | None, str | None]:
    """Accept either mmol directly or a mass plus the limiting reagent's MW.

    Both, because the scale is asked for in conversation and a chemist says
    whichever they have -- "5 mmol" or "start from 1 g of the bromide"."""
    if not isinstance(scale, dict):
        return None, "scale must be an object, e.g. {\"mmol\": 5} or {\"mass_mg\": 1000, \"mw\": 171.04}"
    if scale.get("mmol") is not None:
        v = float(scale["mmol"])
        return (v, None) if v > 0 else (None, "scale: mmol must be greater than 0")
    mg, mw = scale.get("mass_mg"), scale.get("mw")
    if mg is not None and mw:
        if float(mg) <= 0 or float(mw) <= 0:
            return None, "scale: mass_mg and mw must both be greater than 0"
        return float(mg) / float(mw), None
    if mg is not None and not mw:
        return None, "a scale given as mass_mg also needs the limiting reagent's molecular weight (mw, or its smiles)"
    return None, "no scale given: pass mmol, or mass_mg plus the limiting reagent's mw (or smiles)"


def compute(scale: dict, reagents: list) -> dict:
    """Amounts for one step.

    scale     -- {"mmol": 5.0} or {"mass_mg": 1000, "mw": 171.04}
    reagents  -- [{"name", "equiv", "form", ...form-specific data...}, ...]
                 a "solvent" entry takes `concentration` instead of `equiv`.

    Returns {"ok", "scale_mmol", "reagents": [...], "ask": [...]}. `ok` is False
    when anything is missing, and `ask` then holds one question per gap. Entries
    that COULD be computed are still returned, so a single unknown density does
    not cost the other six rows."""
    mmol0, err = scale_to_mmol(scale)
    if err:
        return {"ok": False, "scale_mmol": None, "reagents": [], "ask": [err]}

    rows, ask = [], []
    for i, r in enumerate(reagents or []):
        name = str(r.get("name") or "reagent %d" % (i + 1))
        form = r.get("form")
        if form not in REQUIRED:
            ask.append("%s: form must be one of %s (got %r). Solid, neat liquid, solution and "
                       "solvent are computed differently, so there is no default."
                       % (name, "/".join(REQUIRED), form))
            continue

        missing = [k for k in REQUIRED[form] if r.get(k) in (None, "")]
        if form != "solvent" and r.get("equiv") in (None, ""):
            missing.append("equiv")
        if missing:
            ask.append("%s (%s): also needs %s."
                       % (name, _FORM.get(form, form),
                          ", ".join(_LABEL.get(k, k) for k in missing)))
            continue

        row: dict = {"name": name, "form": form}

        if form == "solvent":
            conc = float(r["concentration"])
            if conc <= 0:
                ask.append("%s: concentration must be greater than 0" % name); continue
            ml = mmol0 / 1000.0 / conc * 1000.0     # mmol -> L -> mL
            row.update(concentration=conc, volume_ml=_round_sig(ml),
                       display=_fmt_vol(ml),
                       note="at %g M (relative to the limiting reagent)" % conc)
            rows.append(row); continue

        eq = float(r["equiv"])
        if eq <= 0:
            ask.append("%s: equiv must be greater than 0" % name); continue
        mmol = mmol0 * eq
        row.update(equiv=eq, mmol=_round_sig(mmol))

        if form == "solution":
            molar = float(r["molarity"])
            if molar <= 0:
                ask.append("%s: molarity must be greater than 0" % name); continue
            ml = mmol / molar
            row.update(molarity=molar, volume_ml=_round_sig(ml),
                       display=_fmt_vol(ml),
                       note="%g M solution" % molar)
            rows.append(row); continue

        mw = float(r["mw"])
        if mw <= 0:
            ask.append("%s: molecular weight must be greater than 0" % name); continue
        mg = mmol * mw
        row["mw"] = mw

        purity = r.get("purity")
        if purity not in (None, ""):
            purity = float(purity)
            if not 0 < purity <= 1:
                ask.append("%s: purity must be a fraction between 0 and 1 (write 60%% as 0.6)" % name)
                continue
            weigh = mg / purity
            row.update(mass_mg=_round_sig(mg), weigh_mg=_round_sig(weigh),
                       purity=purity, display=_fmt_mass(weigh),
                       note="purity %.0f%%: weigh out %s (%s of active reagent)"
                            % (purity * 100, _fmt_mass(weigh), _fmt_mass(mg)))
        else:
            row.update(mass_mg=_round_sig(mg), display=_fmt_mass(mg))

        if form == "neat":
            dens = float(r["density"])
            if dens <= 0:
                ask.append("%s: density must be greater than 0" % name); continue
            ml = (row.get("weigh_mg", mg)) / 1000.0 / dens
            row.update(density=dens, volume_ml=_round_sig(ml),
                       display=_fmt_vol(ml),
                       note=(row.get("note", "") + "; " if row.get("note") else "")
                            + "d %g g/mL: %s (%s)" % (dens, _fmt_vol(ml), _fmt_mass(mg)))
        rows.append(row)

    # scale_mw IS RETURNED EVEN THOUGH NOTHING HERE NEEDS IT. Measured on a
    # live run: with only scale_mmol in the reply, Claude wanted to write
    # "4-bromotoluene (MW 171.xx)" in the reagent table, had nothing to quote, and
    # recalled 171.03 against the real 171.04. Harmless that time -- the mmol it
    # used came from here -- but it is the one number in the write-up that did not
    # come from a tool, and the whole design is that none of them do. A field
    # costs nothing; a remembered molecular weight is a number nobody checked.
    return {"ok": not ask, "scale_mmol": _round_sig(mmol0),
            "scale_mw": (round(float(scale["mw"]), 2)
                         if isinstance(scale, dict) and scale.get("mw") else None),
            "reagents": rows, "ask": ask}
