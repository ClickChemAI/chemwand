"""preflight.py — refuse a bad spec BEFORE any geometry runs.

WHY THIS IS A LAYER AND NOT SIX MORE `if`s. Every engine already checks the keys
it reads, and those checks are good; what none of them can do is check the spec as
a WHOLE before the first coordinate is computed. So each engine grew its own
partial answer at its own depth, and the gaps between them all have the same
shape:

    the figure is drawn, it is missing what the caller asked for, and
    nothing anywhere says so.

Measured without this module -- each of these drew a figure and returned
`ok: true`:

  * `draw_linear_scheme(connectors=[{"from": "a", "to": "zzz"}])` -- layout_linear
    does `struct_pos.get(...)` then `if fpos is None: continue`, so the arrow is
    simply not drawn (layout_linear.py, the connector loop)
  * `draw_molecule(smiles="c1ccccc1[*:1]", rlabels={"1": "Ar", "9": "X"})` -- "9"
    names a placeholder that is not in the SMILES
  * `draw_scope_table(header={"items": [{"arrow": {"reagent_smiles": ...}}]})` --
    the reagent is accepted and never drawn; the same key on a `draw_linear_scheme`
    arrow DOES draw one, so the two tools disagree in silence
  * `arc_style={"hed": "Full"}` -- `layout_cycle._arc_kw` filters to `_ARC_KEYS`
  * `page="lettr"` on all six -- `cdxml_core.finalize` has `fit` as its `else`

and two more that crashed instead of complaining, which is the same defect wearing
a different hat: `KeyError: 'ZZ'` for a radial slot that is not one of the eight,
`KeyError: 'nope'` for `edges[].ring` naming a ring that is not there.

THE MCP SCHEMA CANNOT CLOSE THIS. `chemwand_mcp` publishes every object as
`additionalProperties: false`, but FastMCP validates a call against a pydantic
model built from the SIGNATURE, and pydantic's default for an extra field is to
IGNORE it -- the note above `_obj` in that file says so and names this module as
the answer.

WHAT THIS MODULE IS NOT. It is not a second copy of the engines' key whitelists:
those live in the layer that reads the keys and stay there. It checks the three
things no single-level whitelist can see -- an id used twice, a name that refers to
nothing, and a value outside a closed set -- and it does that on the caller's own
spec, before a fragment is built.
"""

from __future__ import annotations

import re


class SpecError(Exception):
    """The spec is wrong, and the message can name the field and the way out.

    THREE FIELDS RATHER THAN A STRING, because two readers want different things
    out of one error. A person wants the sentence; a test (and, later, a client
    that renders the failure) wants to key off `field` without parsing English
    that will be reworded. `valid_options` is separate for the same reason and
    carries the whole of the answer to "well, what SHOULD I have written" --
    which for a dangling reference is the ids the figure actually has, and for a
    closed set is the set.

    NOT a ValueError. The engines raise ValueError for everything they check, and
    keeping the two apart is what lets a test say "this was refused BEFORE the
    geometry" rather than "something, somewhere, was unhappy".
    """

    def __init__(self, field: str, reason: str, valid_options=None):
        self.field = str(field)
        self.reason = str(reason)
        self.valid_options = (None if valid_options is None
                              else [str(o) for o in valid_options])
        msg = "%s: %s" % (self.field, self.reason)
        if self.valid_options:
            msg += " Valid: %s." % ", ".join(self.valid_options)
        elif self.valid_options == []:
            # An EMPTY list is a real answer and a different one from "no list":
            # "nothing is valid here" is what a dangling reference into an empty
            # figure means, and printing "Valid: ." would be worse than saying it.
            msg += " (this figure defines none)"
        super().__init__(msg)


# ---------------------------------------------------------------------------
# the three checks
# ---------------------------------------------------------------------------

def check_unique_ids(entries, field: str) -> None:
    """One id, one molecule. Raise when two entries share an id and DIFFER.

    `entries` is [(path, id, smiles)] in the caller's own order; `path` is what
    goes in the message, so it must be the recipe path (`items[2]`), not a count.

    SHARING AN ID IS LEGITIMATE and stays legitimate: a struck-out failed target
    and the actual product are routinely both labelled 8a, and `_register` builds
    that molecule once and puts it in both places on purpose. There is a byte-gate
    regression case that does exactly this across two containers, and it exists to
    prove an address cannot be keyed on the structures key.

    WHAT IS NOT LEGITIMATE is two DIFFERENT molecules under one id. `structures`
    is a dict keyed by id, so the second one is never built -- and without this
    check the figure comes back with every position filled, a note in the notes
    field and `ok: false`. A report on a figure that is already wrong is the weaker
    answer: refusing is the stronger one, because the caller cannot miss it.

    An id that is blank or absent means "no id given" and never collides
    (`"id": ""` on every header item is the natural thing to write, and `dict.get`
    falls back only on a MISSING key, so without this rule those would all pile
    onto one entry).
    """
    seen: dict = {}
    for path, sid, smiles in entries:
        sid = str(sid or "").strip()
        if not sid:
            continue
        if sid not in seen:
            seen[sid] = (path, smiles)
            continue
        first_path, first_smiles = seen[sid]
        if first_smiles == smiles:
            continue                    # the same molecule twice: fine, on purpose
        raise SpecError(
            field,
            "id %r is used twice for DIFFERENT structures -- %s has %r and %s has "
            "%r. A figure holds one structure per id, so the second one would "
            "never be drawn and every position would still have something in it. "
            "Give the second its own id, or drop its 'smiles' if it was meant to "
            "repeat the first."
            % (sid, first_path, first_smiles, path, smiles))


def check_references(used, defined, field: str, what: str,
                     note: str = "") -> None:
    """Every name in `used` must be in `defined`. Raise on the first that is not.

    `used` is [(path, name)] and `defined` is the set of names that exist. `what`
    names the kind of thing being referred to ("cell", "edge", "ring"), because
    "'zzz' is not a cell in this figure" reads as a sentence and "'zzz' is not in
    the set" does not.

    THE OPTIONS ARE THE FIGURE'S OWN NAMES, always, even when there are none. A
    dangling reference is nearly always a typo, and the fix is one of the names
    already on the page -- so the message carries them rather than making the
    caller go and look.

    `note` is for the case where the option list alone would mislead: a name can be
    absent from it because it was never written OR because the caller is naming
    something that is not addressable at all. Same parameter, same purpose, as
    check_choice's.
    """
    have = {str(d) for d in defined}
    for path, name in used:
        if str(name) in have:
            continue
        raise SpecError(
            path,
            "%r is not %s %s in this figure.%s"
            % (str(name), _article(what), what, (" " + note) if note else ""),
            sorted(have))


def check_choice(value, options, field: str, note: str = "") -> None:
    """`value` must be one of `options`. Raise naming the set.

    THE SET COMES FROM THE CODE THAT BRANCHES ON IT, never from a docstring. Every
    one of these is silently accepted without it: `page_mode` falls through to `fit`
    in `cdxml_core.finalize`'s `else`, and the rest are read with `.get` or
    compared with `==` and simply never match.

    None passes. Every one of these arguments has a default, and "not given" is a
    different statement from "given something outside the set".
    """
    if value is None:
        return
    if str(value) in {str(o) for o in options}:
        return
    raise SpecError(field, "%r is not one of the values this reads.%s"
                    % (value, (" " + note) if note else ""), options)


def _article(word: str) -> str:
    return "an" if word[:1].lower() in "aeiou" else "a"


# ---------------------------------------------------------------------------
# what every engine shares
# ---------------------------------------------------------------------------

# cdxml_core.finalize's own branches, in its own order. `fit` is the `else`, which
# is exactly why an unchecked unknown value is indistinguishable from asking for it.
PAGE_MODES = ("letter", "single", "double", "auto", "fit")


def _cells(rows):
    """(path, cell) for the two engines whose cells live in `rows`.

    scope_table and snaking, so `rows[0][1]` in a message is built once and reads
    the same in both. linear has its own walker because its `items` is EITHER a
    flat list of cells or a list of rows, and cycle's cells are `nodes`, which is
    neither shape.

    A row that is a dict rather than a list is a scope table's `{"section": ...}`
    marker; it is yielded as `rows[i]` so a stray key on one still gets a path.
    """
    out = []
    for i, row in enumerate(rows or []):
        if isinstance(row, list):
            for j, cell in enumerate(row):
                if isinstance(cell, dict):
                    out.append(("rows[%d][%d]" % (i, j), cell))
        elif isinstance(row, dict):
            out.append(("rows[%d]" % i, row))
    return out


def _id_entries(pairs):
    """[(path, id, smiles)] for the cells that name a structure.

    Cells with no `smiles` are skipped: an id on a written-out species names no
    molecule, so it cannot collide with one in the structures dict.
    """
    return [(path, c.get("id"), c.get("smiles"))
            for path, c in pairs if c.get("smiles")]


# The eight slots layout_radial.POS names -- the 3x3 grid around the hub, without
# the hub itself, which is `center_smiles` and not a product. Read off that dict:
# unchecked, a slot outside it comes back as `KeyError: 'ZZ'` from the geometry.
RADIAL_SLOTS = ("TL", "T", "TR", "L", "R", "BL", "B", "BR")

# layout_cycle._ARC_KEYS, and the reason it is named here rather than imported is
# in `preflight_cycle`.
ARC_KEYS = ("head", "tail", "head_type", "head_size", "head_width",
            "head_center", "line")

# ONE HIGHLIGHT SPEC, and the keys `decorate.apply()` actually reads.
#
# THE DEFINITION LIVES HERE. The set names a dict its own module reads with
# .get() and checks NOT AT ALL, so a schema that publishes it as closed is not
# enough: at runtime a typo goes straight through -- `{"style": "color",
# "colour": "red"}` draws an uncoloured figure and says nothing, on all six tools.
# chemwand_mcp's schema aliases this name, so there is one list rather than
# two that can drift.
HIGHLIGHT_KEYS = frozenset({
    "on", "style", "color",                          # every style
    "labels",                                        # "locants"
    "text", "offset", "arrow", "bulge", "size",      # "annotate"
    "tail", "head",                                  # "annotate" ends by :n map
    "tail_dir", "head_dir",                          # "annotate" Bezier tangents
    "maps",                                          # address by the subject's :n
    "charge", "center", "radius",                    # "delocalized"
    "nth",                                           # "annotate" one match only
    "amp", "wavelength", "overhang",                 # "squiggle"
    "n",                                             # "lonepair" / "radical"
})


def check_closed(obj, keys, field: str, note: str = "") -> None:
    """Every key of `obj` must be one `keys` names. Raise on the first that is not.

    THE ONE PLACE a closed nested object is enforced. Every set this module holds
    describes a dict its own module reads with `.get()` -- which is precisely why
    none of them could be enforced by the reader: `.get()` on a key nobody wrote
    is indistinguishable from `.get()` on a key nobody spelled right. So the check
    has to be here, ahead of the reader.

    `note` says what the silence COSTS, because that is the sentence a caller
    needs: "unknown key" tells them they typed something wrong, and "so the
    decoration comes out with its defaults" tells them why their figure looks
    the way it does.
    """
    if not isinstance(obj, dict):
        return
    for k in obj:
        check_choice(k, sorted(keys), "%s key" % field, note)


def check_highlights(specs, field: str) -> None:
    """Every key of every highlight in `specs` must be one decorate.apply() reads.

    A LIST OF SPECS, because that is the shape every tool takes it in -- the
    figure-level `highlight` argument and each cell's own are the same object, and
    checking them in one place is what makes "all six tools" true rather than
    "five of six and whichever was remembered".
    """
    for i, hl in enumerate(specs or []):
        check_closed(hl, HIGHLIGHT_KEYS, "%s[%d]" % (field, i),
                     "An unknown key here is read by nobody, so the "
                     "decoration comes out with its defaults.")


# ---- four more closed sets ---------------------------------------------------
#
# Same story as HIGHLIGHT_KEYS above, four more times: chemwand_mcp's schema
# needs a frozenset for `_obj` to point at, and each names a dict whose own module
# reads it with `.get()` and checks it not at all. Measured: twelve places publish
# one of these four, and unchecked, a typo'd key comes back BYTE-IDENTICAL to the
# correct spec in all twelve.
#
# The definitions live here and that file aliases them, so there is one list.

# arrows.apply_style(). Every one of these is read with `style.get(...)`, so an
# unknown key leaves the arrow at its default and says nothing.
ARROW_STYLE_KEYS = frozenset({
    "kind", "heads", "headType", "shaftSpacing", "headSize", "width", "barb",
    "curve", "noGo", "line", "color", "dipole", "hashSpacing", "lineWidth",
    "boldWidth", "equilibriumRatio",
})
# arrows._resolve_arcs(), one entry of an arrow's `arcs`. NOT layout_cycle's
# `_ARC_KEYS` -- the names are alike and the sets are not: that one is the cycle's
# arc STYLE (head / tail / head_type / ...) and lives in ARC_KEYS above.
ARC_ENTRY_KEYS = frozenset({"at", "off", "head", "tail", "line", "dir", "bulge"})
# legend.measure() / compose_legend(), and one CELL of a legend row.
LEGEND_KEYS = frozenset({"rows", "header", "align", "cols", "fill"})
LEGEND_CELL_KEYS = frozenset({"swatch", "text"})

# THE SHAFT STYLES AN ARROW CAN ACTUALLY HAVE.
#
# `arrows.apply_style` maps this onto `LineType`, and the map ends in
# `.get(line, "Dashed")` -- so without this set, EVERY unrecognised value comes
# out dashed. Measured: `line="dotted"` wrote LineType="Dashed" with no error and
# no warning. That, not a dead key, is the fault in this family; a dead-key
# comparison reads it as dead only when its own baseline (`"plain"`) is not a
# legal value either, so it compares Dashed against Dashed.
#
# `wavy` IS NOT HERE, and it is the one entry that needs a measurement rather than
# a rule: `arrows.py` records that ChemDraw STORES LineType="Wavy" on an arrow --
# it survives a round trip -- and renders it identically to plain, confirmed
# visually. So it is a published value the format cannot express. Left out pending
# a hand-drawn ChemDraw sample. Note that a byte-comparison criterion cannot see
# this case at all: the bytes move and the picture does not.
ARROW_LINE_TYPES = ("solid", "dashed", "bold")

_STYLE_NOTE = ("An unknown key here is read by nobody, so the arrow comes out "
               "with its defaults and nothing says why.")
_ARC_NOTE = ("An unknown key here is read by nobody, so the arc comes out with "
             "its defaults and nothing says why.")
_LEGEND_NOTE = ("An unknown key here is read by nobody, so the legend comes out "
                "with its defaults and nothing says why.")


def check_line_type(value, field: str) -> None:
    """An arrow's shaft style. One place, because two code paths write it.

    `arrows.apply_style` (for a `style`) and `cdxml_core.arc_arrow` (for an
    `arcs[]` entry) are different functions that end up at the same `LineType`
    attribute, so the value they accept has to be decided once.
    """
    if value is None:
        return
    if str(value).lower() == "wavy":
        raise SpecError(
            field,
            "'wavy' is not available on an arrow. ChemDraw STORES LineType=\"Wavy\" "
            "-- it survives a round trip -- and draws the shaft identically to a "
            "plain one, so it would change nothing you can see. For a heavier shaft "
            "use 'bold', or set 'lineWidth' / 'boldWidth' on the same style for "
            "thicker still.",
            list(ARROW_LINE_TYPES))
    check_choice(value, ARROW_LINE_TYPES, field,
                 "Any other value would be drawn as a dashed line without warning.")


def check_arrow_style(style, field: str) -> None:
    check_closed(style, ARROW_STYLE_KEYS, field, _STYLE_NOTE)
    if isinstance(style, dict):
        check_line_type(style.get("line"), "%s.line" % field)


def check_arcs(arcs, field: str) -> None:
    for i, a in enumerate(arcs or []):
        check_closed(a, ARC_ENTRY_KEYS, "%s[%d]" % (field, i), _ARC_NOTE)
        if isinstance(a, dict):
            check_line_type(a.get("line"), "%s[%d].line" % (field, i))


def check_legend(legend, field: str) -> None:
    """A legend, and every CELL of every row of it.

    The rows are a ragged nest -- a row may be a bare string (a NOTE row) or a
    list whose entries are strings or `{"swatch", "text"}` dicts -- so the walk
    checks what is a dict and leaves the rest to the schema. A string cell has no
    keys to spell wrong.
    """
    check_closed(legend, LEGEND_KEYS, field, _LEGEND_NOTE)
    if not isinstance(legend, dict):
        return
    for r, row in enumerate(legend.get("rows") or []):
        if not isinstance(row, list):
            continue
        for c, cell in enumerate(row):
            check_closed(cell, LEGEND_CELL_KEYS,
                         "%s.rows[%d][%d]" % (field, r, c), _LEGEND_NOTE)


def refuse_unread_arrow(arrow, field: str) -> None:
    """`style` / `arcs` on a HEADER or PANEL arrow: refused whole, not spell-checked.

    THE OBJECT IS NOT READ AT ALL HERE, which is a different fault from a typo and
    takes a different answer. `layout_grid._emit_header` draws this arrow with a
    bare `C.arrow(...)` and reads only `above` / `below`; `panels.py` borrows that
    same function (it builds `{"header": {"items": panel["items"]}}`). Measured:
    changing a REAL key -- `headSize`, `bulge` -- moves not one byte.

    So spell-checking the keys would be false comfort: spelled right, it still
    does nothing. It is the same treatment `reagent_smiles` gets -- same fault,
    same cause (the scope table reuses
    `_S_LINEAR_ARROW` whole while its renderer implements a subset of it).
    """
    if not isinstance(arrow, dict):
        return
    for key, what in (("style", "an arrow style"), ("arcs", "attached arcs")):
        if arrow.get(key) is None:
            continue
        raise SpecError(
            "%s.%s" % (field, key),
            "%s is not built for this arrow -- it is accepted and never drawn. "
            "The arrow here is a plain shaft plus the two caption lines 'above' "
            "and 'below', and nothing else of it is read. An arrow that must "
            "carry a style or an attached arc belongs in `draw_linear_scheme`, "
            "whose renderer does read both." % what,
            ["above", "below"])


_MAPNO = re.compile(r"\[\*:(\d+)\]")


# SMARTS-ONLY TOKENS, for `template`.
#
# `template` is a COORDINATE DONOR: `smiles_to_frag.template_mol` computes 2D
# coordinates for it and `GenerateDepictionMatching2DStructure` copies them onto
# every member. A SMARTS has no coordinates to give, so writing one there is not a
# stricter query -- it is a silent, near-total fallback. Measured:
# `[#6]S(=N[#6]=O)[#6]` produced 35/36 fallbacks with nothing said, and the
# engine still returns `ok: true` for it (measured on a four-cell table:
# `warnings: []`, one sentence in `notes`). `template_patt` is the SMARTS slot
# and takes exactly this.
#
# EACH PATTERN IS HERE BECAUSE IT WAS TESTED, not remembered. The criterion was
# "does `MolFromSmiles` refuse it, or does it mean something else there":
#
#   [c,n]  [!C]  [C;X4]  [C&X4]  [CX4]  [CR1]  [D2]  [$(CC)]
#       MolFromSmiles refuses all of these outright.
#   [#6]
#       PARSES as SMILES (it means carbon), which is exactly why that case is silent
#       rather than an error. Nobody writes `[#6]` in a SMILES on purpose, so it
#       is listed -- this one is a judgement about intent, not a parse failure.
#   ~
#       any-bond. RDKit's SMILES parser was lenient about it in this build, so it
#       is listed on the same intent grounds.
#
# THE WILDCARD IS SPLIT, and the split is the whole subtlety of this check.
#
# In SMARTS a bare `*` is "any atom" -- the most permissive thing you can write,
# and as a template it is meaningless: a reference whose atoms are all wildcards
# has nothing to donate. So a `*` OUTSIDE brackets is refused.
#
# A BRACKETED one is not. `[*:1]` is this project's Markush attachment point --
# draw_molecule documents it, `c1ccccc1[*:1]` is a live template value in
# test_placeholders, and `template_mol` has a whole branch
# (TEMPLATE_MARKUSH) for a Markush parent templating a Markush series. `[*]`
# without a map number is bracketed too and passes for the same reason; it is the
# same dummy atom, just unnamed.
#
# Told apart STRUCTURALLY rather than with a lookbehind: strip every bracketed
# atom expression, then ask whether a `*` is left. `c1ccccc1[*:1]` strips to
# `c1ccccc1` and passes; `c1ccccc1*` keeps its star and is refused.
_BRACKETED = re.compile(r"\[[^]]*\]")

_SMARTS_ONLY = (
    (re.compile(r"\[[^]]*,"), "an atom list like [c,n]"),
    (re.compile(r"\[[^]]*[;&]"), "a logical ; or & inside brackets"),
    (re.compile(r"\[[^]]*!"), "a negation like [!C]"),
    (re.compile(r"\$\("), "a recursive SMARTS $(...)"),
    (re.compile(r"\[#\d"), "an atomic number like [#6]"),
    (re.compile(r"\[[^]]*[XDR]\d"), "a query like [CX4] / [D2] / [CR1]"),
    (re.compile(r"~"), "the any-bond operator ~"),
)


def check_template(value, field: str) -> None:
    """`template` must be a drawable structure, not a match pattern.

    "auto" / "off" are the two sentinels the tool reads before anything is built,
    so they pass untouched.
    """
    if value is None:
        return
    s = str(value)
    if s.lower() in ("auto", "off"):
        return
    if "*" in _BRACKETED.sub("", s):
        raise SpecError(
            field,
            "%r contains a bare `*`. In SMARTS that is \"any atom\", and a "
            "reference made of wildcards has no coordinates to donate -- "
            "`template` has to be a real drawable structure (SMILES). A "
            "BRACKETED placeholder is different and is accepted: `[*:1]` is this "
            "project's Markush attachment point, and a Markush parent may "
            "template a Markush series. If you meant a match pattern, it goes in "
            "`template_patt`." % s,
            ["auto", "off", "<SMILES>"])
    for rx, what in _SMARTS_ONLY:
        if rx.search(s):
            raise SpecError(
                field,
                "%r contains %s, which is SMARTS. `template` is a COORDINATE "
                "DONOR -- its 2D layout is computed once and copied onto every "
                "member -- so it has to be a real drawable structure, written as "
                "SMILES (any member of the series works). A pattern goes in "
                "`template_patt`, which narrows WHICH PART of the template is "
                "matched and is parsed as SMARTS. Left as it is, this does not "
                "match more strictly: it makes nearly every cell fall back to its "
                "own layout, which the figure does not show."
                % (s, what),
                ["auto", "off", "<SMILES>"])


def _rlabel_refs(smiles, rlabels, field: str) -> None:
    """`rlabels` keys must name a placeholder the SMILES actually has.

    `rlabels: {"1": "Ar"}` renames `[*:1]`, so a key with no matching map number
    renames nothing. Unchecked, it is accepted in silence, the worst outcome for
    this key in particular: the figure comes back with the DEFAULT label ("R1")
    on every placeholder and looks like a rendering bug rather than a typo.

    A bare `*` carries no map number and is always drawn "R", so it is not
    addressable here and contributes nothing to the set -- which is why the message
    can say "this SMILES has none" and mean it.
    """
    if not rlabels:
        return
    have = set(_MAPNO.findall(str(smiles or "")))
    check_references([("%s[%r]" % (field, k), k) for k in rlabels],
                     have, field, "placeholder")


# ---------------------------------------------------------------------------
# one per engine
# ---------------------------------------------------------------------------
#
# EACH TAKES THE TOOL'S OWN ARGUMENTS as a dict, and is called from the tool
# before anything is built. Not the engine recipe: by the time a recipe exists the
# SMILES have been parsed and the structures dict has already collapsed two cells
# onto one id, which is the very thing check_unique_ids is here to catch.

def preflight_molecule(spec: dict) -> None:
    check_choice(spec.get("page"), PAGE_MODES, "page")
    check_highlights(spec.get("highlight"), "highlight")
    check_legend(spec.get("legend"), "legend")
    _rlabel_refs(spec.get("smiles"), spec.get("rlabels"), "rlabels")
    # NO ID CHECK, and that is a statement rather than an omission: this tool draws
    # exactly one structure and has no id space to collide in. It is the one
    # documented exception to going through `_register` for the same reason
    # (pinned by test_no_tool_builds_a_structure_behind_the_registrar).


_RESERVE_KEYS = frozenset({"after", "h"})


def check_reserve(reserve, rows_of_ids, field: str = "reserve") -> None:
    """`reserve` asks the packer to leave `h` points under the row a cell sits in.

    REFUSED RATHER THAN IGNORED, on preflight's terms, and the failure it prevents
    is a quiet one: a caller who asks for room and silently gets none drops a
    sibling figure straight onto the next row, and the only symptom is an overlap
    that reads as the packer's fault.

    `rows_of_ids` is one set of cell ids per row, in row order -- the last row is
    refused because nothing sits below it to be pushed down, so the request would
    be granted and change nothing.
    """
    have = sorted({i for r in rows_of_ids for i in r})
    for i, rv in enumerate(reserve or []):
        path = "%s[%d]" % (field, i)
        if not isinstance(rv, dict):
            raise SpecError(path, "each entry is a dict, e.g. "
                                  '{"after": "8-2", "h": 208}.', [])
        unknown = set(rv) - _RESERVE_KEYS
        if unknown:
            raise SpecError(path, "%s is not a reserve key."
                            % ", ".join(sorted(unknown)), sorted(_RESERVE_KEYS))
        try:
            h = float(rv.get("h"))
        except (TypeError, ValueError):
            h = None
        if h is None or h <= 0:
            raise SpecError("%s.h" % path,
                            "the room to leave, in points, as a positive number; "
                            "got %r." % (rv.get("h"),), [])
        check_references([("%s.after" % path, rv.get("after"))], have,
                         field, "cell",
                         "'after' names the cell whose ROW gets the room below it.")
        ri = next(r for r, ids in enumerate(rows_of_ids)
                  if str(rv.get("after")) in ids)
        if ri == len(rows_of_ids) - 1:
            raise SpecError(
                "%s.after" % path,
                "%r is in the last row, and nothing sits below it to be pushed "
                "down -- the room would be granted and change nothing. Put the "
                "sibling below the figure instead." % str(rv.get("after")), have)


def check_span_width(span, field: str = "span") -> None:
    """`span` tells an engine how wide it may be, instead of asking the page.

    NAMED IN FULL to keep it apart from `check_cell_spans`, which validates the
    OTHER span -- a distance between two named cells. The two take different
    shapes (a number vs a list of dicts), refuse different things, and live on
    different tools; `check_span` / `check_spans` would have differed by one
    letter at every call site.

    REFUSED RATHER THAN IGNORED, on preflight's terms and for check_reserve's
    reason.
    The number is read three times -- the column count, the cell spread, and the
    "too wide" note -- and a nonsense one does not crash: `int(-200 // cell_w)`
    is negative, `max(1, ...)` rescues it, and the figure comes out in ONE column
    looking like the packer's own judgement. Zero does the same thing, and a
    string is the likelier typo of the two.
    """
    if span is None:
        return
    try:
        w = float(span)
    except (TypeError, ValueError):
        w = None
    if w is None or w <= 0:
        raise SpecError(field,
                        "the width to lay the figure out within, in points, as a "
                        "positive number; got %r." % (span,), [])


_CELL_SPAN_KEYS = frozenset({"from", "to", "w"})


def check_cell_spans(spans, rows, field: str = "span",
                     arrows_between: bool = False) -> None:
    """`span` asks the packer to put two cells exactly `w` points apart.

    `rows` is one list of cell dicts per row, in row order -- the SHAPE of the row
    is needed here, not just its ids, because two of the three refusals below are
    about what sits between the endpoints.

    REFUSED RATHER THAN IGNORED, on preflight's terms and for check_reserve's
    reason, and the three refusals are the three ways the request cannot be met at
    all:

      DIFFERENT ROWS   a span is a horizontal distance and rows are packed
                       independently, so there is no arrow the stretch could
                       spend. The parent moves its own part instead.
      NOTHING BETWEEN  only an ARROW stretches. Two adjacent cells, or a pair with
                       only a "+" between them, are as far apart as their own
                       widths make them, and the request would change nothing.
      NO SUCH CELL     via check_references, which lists the ids that do exist.

    `arrows_between` IS FOR SNAKING AND IS NOT A CONVENIENCE. A snaking row is
    written as structures and nothing else -- its cells take no `arrow` key at all
    -- and the tool puts one between every adjacent pair on the way to the linear
    packer. So the "nothing between" test, run on what the caller wrote, would
    refuse every span on that engine while the engine could have granted it.

    What is NOT refused here is a `w` SMALLER than the cells already are. That one
    needs the widths, which this layer does not have and will not guess -- the
    engine reports it instead (layout_linear.SPAN_SHORT).
    """
    have = sorted({str(c.get("id")).strip() for row in rows for c in row
                   if isinstance(c, dict) and c.get("id") is not None})
    for i, sv in enumerate(spans or []):
        path = "%s[%d]" % (field, i)
        if not isinstance(sv, dict):
            raise SpecError(path, "each entry is a dict, e.g. "
                                  '{"from": "III", "to": "II", "w": 308.7}.', [])
        unknown = set(sv) - _CELL_SPAN_KEYS
        if unknown:
            raise SpecError(path, "%s is not a span key." % ", ".join(sorted(unknown)),
                            sorted(_CELL_SPAN_KEYS))
        try:
            w = float(sv.get("w"))
        except (TypeError, ValueError):
            w = None
        if w is None or w <= 0:
            raise SpecError("%s.w" % path,
                            "how far apart to put them, in points, as a positive "
                            "number; got %r." % (sv.get("w"),), [])
        check_references([("%s.from" % path, sv.get("from")),
                          ("%s.to" % path, sv.get("to"))], have, field, "cell",
                         "'from' and 'to' name two cells in ONE row.")
        a, b = str(sv.get("from")).strip(), str(sv.get("to")).strip()
        if a == b:
            raise SpecError(path, "'from' and 'to' are both %r, and a cell is "
                                  "zero from itself." % a, have)
        where = {}
        for ri, row in enumerate(rows):
            for ci, c in enumerate(row):
                if isinstance(c, dict) and str(c.get("id")).strip() in (a, b):
                    where[str(c["id"]).strip()] = (ri, ci)
        if where[a][0] != where[b][0]:
            raise SpecError(
                path,
                "%r is in row %d and %r is in row %d. A span is a horizontal "
                "distance and rows are packed independently, so there is no arrow "
                "between them to lengthen -- move the sibling part instead."
                % (a, where[a][0], b, where[b][0]), have)
        row = rows[where[a][0]]
        lo, hi = sorted((where[a][1], where[b][1]))
        if arrows_between:
            continue
        if not any(isinstance(c, dict) and "arrow" in c for c in row[lo + 1:hi]):
            raise SpecError(
                path,
                "nothing between %r and %r can stretch: only an ARROW takes up the "
                "slack, and the cells between them are %s. Their distance is fixed "
                "by their own widths, so the span would be granted and change "
                "nothing."
                % (a, b, ", ".join(sorted({k for c in row[lo + 1:hi]
                                           if isinstance(c, dict) for k in c
                                           if k in ("s", "smiles", "plus", "text")}))
                   or "nothing at all"), have)


def check_one_under_per_row(rows, field: str = "items") -> None:
    """`under` x-aligns a structure with one in an EARLIER row, by shifting the
    WHOLE row. A row therefore has exactly one of them.

    REFUSED RATHER THAN IGNORED, and this is the quietest of the four such cases:
    layout_linear walks the row and does `shifts[ri] += ...` per `under`, each
    time recomputing the offset against the shift the previous one just applied.
    So the LAST one lands exactly and every earlier one is undone -- a row that
    asked to sit under two things sits under one of them, the figure is drawn,
    and nothing anywhere says which one won.

    (The DERIVED `under` -- read off the connectors -- already guards this:
    "one row is pulled at most once". Only the hand-written key needs this check.)
    """
    for ri, row in enumerate(rows):
        if not isinstance(row, list):
            continue
        got = [(ci, c) for ci, c in enumerate(row)
               if isinstance(c, dict) and c.get("under")]
        if len(got) > 1:
            raise SpecError(
                "%s[%d]" % (field, ri),
                "%d cells in this row carry 'under' (%s). `under` shifts the whole "
                "row, so a row can only sit under ONE structure -- the others would "
                "be silently overruled. Keep the one that matters (usually the row's "
                "FIRST structure, which is what makes a connector into it vertical) "
                "and drop the rest."
                % (len(got), ", ".join("%s[%d][%d] -> %r"
                                       % (field, ri, ci, c["under"]) for ci, c in got)),
                [])


def preflight_linear(spec: dict) -> None:
    check_choice(spec.get("page"), PAGE_MODES, "page")
    check_template(spec.get("template"), "template")
    check_choice(spec.get("align"), ("left", "center"), "align",
                 "Only 'left' and 'center' are supported; any other value would be drawn "
                 "as 'center'.")
    # `items` IS EITHER SHAPE. The tool takes a single row (a flat list of cells)
    # or a list of rows, and decides with `isinstance(items[0], list)` -- the same
    # test, spelled the same way, because reading it only one way here tells a
    # multi-row figure that it "defines no cells".
    items = spec.get("items") or []
    if items and isinstance(items[0], list):
        cells = [("items[%d][%d]" % (r, c), it)
                 for r, row in enumerate(items) if isinstance(row, list)
                 for c, it in enumerate(row) if isinstance(it, dict)]
    else:
        cells = [("items[%d]" % i, it) for i, it in enumerate(items)
                 if isinstance(it, dict)]
    entries = _id_entries(cells)
    # A REAGENT ABOVE AN ARROW IS A STRUCTURE LIKE ANY OTHER and goes into the same
    # structures dict under `reagent_id`, so it collides with a cell of that id.
    # Left out, the one id space would be checked in two halves.
    for path, it in cells:
        arrow = it.get("arrow")
        if isinstance(arrow, dict) and arrow.get("reagent_smiles"):
            entries.append(("%s.arrow" % path, arrow.get("reagent_id"),
                            arrow.get("reagent_smiles")))
    check_unique_ids(entries, "items")
    for path, it in cells:
        _rlabel_refs(it.get("smiles"), it.get("rlabels"), "%s.rlabels" % path)
        check_highlights(it.get("highlight"), "%s.highlight" % path)
        arrow = it.get("arrow")
        if isinstance(arrow, dict):
            check_highlights(arrow.get("reagent_highlight"),
                             "%s.arrow.reagent_highlight" % path)
            check_arrow_style(arrow.get("style"), "%s.arrow.style" % path)
            check_arcs(arrow.get("arcs"), "%s.arrow.arcs" % path)
    # `connectors` name two PLACED structures by id. layout_linear looks each up in
    # `struct_pos` and does `if fpos is None or tpos is None: continue` -- so a typo
    # is an arrow that is quietly not drawn.
    #
    # THE CELL IDS, NOT `entries`. `entries` also holds the reagents above arrows,
    # and `struct_pos` is filled only for `it["type"] == "s"` -- a row cell. So a
    # connector naming a reagent is dangling too, and taking the wider set would
    # have let exactly that one through.
    have = {str(c.get("id") or "").strip() for _p, c in cells
            if c.get("smiles") and str(c.get("id") or "").strip()}
    # THE RIGHT REASON, NOT THE FIRST ONE THAT FITS. A `text` cell has an id and is
    # not in `have`, so check_references alone would reject it with "a connector
    # endpoint must be a cell that gave itself an 'id'" -- which is false of the
    # cell in front of the caller, who then goes looking for the id they can see
    # they wrote. What is actually true is narrower and is a limit of the engine:
    # `struct_pos` is filled only for drawn structures, so a written-out species
    # has no geometry for an arrow to end on.
    _texts = {str(c.get("id") or "").strip(): p for p, c in cells
              if c.get("text") is not None and str(c.get("id") or "").strip()}
    used = []
    for ci, con in enumerate(spec.get("connectors") or []):
        if not isinstance(con, dict):
            continue
        check_arrow_style(con.get("style"), "connectors[%d].style" % ci)
        for end in ("from", "to"):
            if con.get(end) is not None:
                used.append(("connectors[%d].%s" % (ci, end), con[end]))
    for _path, _ref in used:
        _r = str(_ref).strip()
        if _r in _texts and _r not in have:
            raise SpecError(
                _path,
                "%r is a `text` cell, and a connector cannot end on one: the "
                "packer records a position only for DRAWN structures, so there is "
                "nothing for the arrow to reach. Draw that species with `smiles` "
                "and the connector works, or write the step as a row arrow instead."
                % _r, sorted(have))
    check_references(
        used, have, "connectors", "cell",
        "A connector endpoint must be a cell that gave itself an 'id'. The keys "
        "this layer assigns to an unnamed cell (s1, s2, ...) are internal and are "
        "NOT usable as a target -- name the cell you want to point at.")
    # `items` is either shape here too -- same test, same spelling, same reason as
    # the one at the top of this function.
    _rows = items if (items and isinstance(items[0], list)) else [items]
    check_reserve(spec.get("reserve"),
                  [{str(c.get("id")).strip() for c in row
                    if isinstance(c, dict) and c.get("smiles") and c.get("id")}
                   for row in _rows if isinstance(row, list)])
    # THE ROWS THEMSELVES, not a set of ids: two of check_cell_spans' three
    # refusals are about what sits BETWEEN the endpoints, which a set has thrown
    # away.
    check_cell_spans(spec.get("span"),
                     [row for row in _rows if isinstance(row, list)])
    check_one_under_per_row(_rows)


    # A FRAMELESS LONE LIGAND BESIDE THE ROW BELONGS IN THE ROW. Same rule as
    # draw_scope_table's (see preflight_scope_table), measured here too: the same
    # Suzuki row with one phosphine came out 419 pt wide in a slot and 373 pt with
    # the ligand as a trailing item -- 11%, the same order as measured on scope
    # tables.
    #
    # ONLY WHEN IT HAS NO FRAME AND NO FILL, and that is the whole difference from
    # the scope-table rule. A BOXED ligand beside the arrow is a published idiom
    # and the box is the point of it; a frameless structure parked to one side is
    # drawn exactly as an inline item would be, so the 46 pt buys nothing. On a
    # scope table the cost is a column of the GRID, so there the frame does not
    # earn it back and the rule fires either way.
    if items:
        for _i, _sl in enumerate(spec.get("slots") or []):
            if not isinstance(_sl, dict):
                continue
            if _sl.get("frame") or _sl.get("fill"):
                continue
            if str(_sl.get("side") or "right") not in ("right", "left"):
                continue
            _st = [x for x in (_sl.get("structures") or []) if isinstance(x, dict)]
            if len(_st) != 1:
                continue
            _id = _st[0].get("id") or _st[0].get("label") or "L1"
            raise SpecError(
                "slots[%d]" % _i,
                "one structure, no frame, parked to the side -- that is a ligand, and a "
                "frameless slot draws it exactly as a trailing `items` cell would while "
                "making the figure about 10%% wider. Put it at the end of `items`:\n    "
                "items=[..., {\"smiles\": ..., \"id\": %r}]\nKeep the slot if the box IS "
                "the point -- give it a `frame` (that is the published ligand box and it "
                "is not refused), or put two structures in it." % _id)


def preflight_scope_table(spec: dict) -> None:
    check_choice(spec.get("page"), PAGE_MODES, "page")
    check_span_width(spec.get("span"))
    check_choice(spec.get("struct_valign"), ("top", "center", "bottom"),
                 "struct_valign")
    # BOTH LEVELS: a section may re-declare its own `template`, and
    # `_section_template` hands whatever it finds straight to `template_mol`, so
    # checking only the table's would leave the same silence one level down.
    check_template(spec.get("template"), "template")
    for i, r in enumerate(spec.get("rows") or []):
        if isinstance(r, dict) and r.get("section") is not None:
            check_template(r.get("template"), "rows[%d].template" % i)
    cells = _cells(spec.get("rows"))
    header = spec.get("header")
    hitems = (header.get("items") or []) if isinstance(header, dict) else []
    hcells = [("header.items[%d]" % i, it) for i, it in enumerate(hitems)
              if isinstance(it, dict)]
    check_unique_ids(_id_entries(cells) + _id_entries(hcells), "rows / header.items")
    for path, c in cells + hcells:
        _rlabel_refs(c.get("smiles"), c.get("rlabels"), "%s.rlabels" % path)
        check_highlights(c.get("highlight"), "%s.highlight" % path)
        check_legend(c.get("legend"), "%s.legend" % path)
    for ci, con in enumerate(spec.get("connectors") or []):
        if isinstance(con, dict):
            check_arrow_style(con.get("style"), "connectors[%d].style" % ci)
    # DECLARED AND NOT BUILT. `reagent_smiles` on a header arrow is accepted by the
    # renderer and the reagent is never drawn -- while the SAME key on a
    # draw_linear_scheme arrow does draw one, so the two tools disagree and neither
    # says so. Refused rather than implemented: this closes the silence, it does not
    # add the feature.
    # ...and `style` / `arcs` on the same arrow, for the same reason and by the
    # same rule. PANELS TOO: `panels.py` renders its `items` through
    # `layout_grid._emit_header`, so a panel arrow is a header arrow and inherits
    # every one of its limits. (Nothing else about a panel item is preflighted --
    # its ids are outside the unique-id check above. That is a separate gap;
    # widening the check is a separate change.)
    for _pi, _pan in enumerate(spec.get("panels") or []):
        for _ii, _it in enumerate((_pan or {}).get("items") or []
                                  if isinstance(_pan, dict) else []):
            if isinstance(_it, dict):
                refuse_unread_arrow(_it.get("arrow"),
                                    "panels[%d].items[%d].arrow" % (_pi, _ii))
    for i, it in enumerate(hitems):
        arrow = it.get("arrow") if isinstance(it, dict) else None
        refuse_unread_arrow(arrow, "header.items[%d].arrow" % i)
        if isinstance(arrow, dict) and arrow.get("reagent_smiles"):
            raise SpecError(
                "header.items[%d].arrow.reagent_smiles" % i,
                "a reagent drawn ABOVE the header arrow is not built for a scope "
                "table -- it is accepted and never drawn. Write the reagent as its "
                "own `header.items` cell before the arrow, or name it in the "
                "arrow's `above` caption; a structure that must sit over the arrow "
                "belongs in `draw_linear_scheme`, which does draw one.")

    # A LONE LIGAND BESIDE THE REACTION ROW OF A PRODUCT GRID BELONGS IN THE ROW.
    # `slots` and `panels` place their content OUTSIDE the reaction row. On a
    # figure that HAS a product grid that costs the figure a column of its own,
    # while a header item shares the row's width and alignment. The docstring says
    # so under `boxes` ("THE LIGAND BESIDE THE ARROW is hdr:<id>", measured 827 pt
    # inline against 987), and documentation alone does not prevent it: measured on
    # published figures, a ligand placed outside the row cost 10-17% of the
    # figure's width. A note in the return value would not catch it: those figures
    # came back clean. So it is refused.
    #
    # FOUR CONDITIONS, and each one is a shape that must keep working:
    #   a header to move it into      -- an optimisation figure with no header has
    #                                    nowhere to put it.
    #   a product grid                -- WITHOUT one, a right-hand panel on the
    #                                    scheme is the engineered shape: its width
    #                                    comes out of that block's band, so it does
    #                                    not cost a column (test_panel_band,
    #                                    test_panels_no_grid -- about half of the
    #                                    surveyed optimisation figures).
    #   anchored to the SCHEME, at its side -- `of: "table"` or `side: "below"` is
    #                                    a block of its own, not part of the row.
    #   one structure and no reaction -- a panel holding an arrow or a plus is the
    #                                    lettered panel and is none of this.
    if hitems and cells:
        for _kind, _list in (("slots", spec.get("slots") or []),
                             ("panels", spec.get("panels") or [])):
            for _i, _entry in enumerate(_list):
                if not isinstance(_entry, dict):
                    continue
                if str(_entry.get("of") or "") != "scheme":
                    continue
                if str(_entry.get("side") or "right") not in ("right", "left"):
                    continue
                _items = _entry.get("structures" if _kind == "slots" else "items")
                _items = [x for x in (_items or []) if isinstance(x, dict)]
                if any("arrow" in x or x.get("plus") for x in _items):
                    continue
                _structs = [x for x in _items if x.get("smiles")]
                if len(_structs) != 1:
                    continue
                _id = _structs[0].get("id") or _structs[0].get("label") or "L1"
                raise SpecError(
                    "%s[%d]" % (_kind, _i),
                    "this holds ONE structure and no reaction and hangs off the side of the "
                    "scheme, so it is a ligand or catalyst standing beside the reaction -- "
                    "and on a figure WITH a product grid a %s puts it in a column of its "
                    "own, which is what makes the figure too wide. Put it in `header.items` "
                    "instead, where it shares the reaction row's width and alignment:\n    "
                    "header={\"items\": [..., {\"smiles\": ..., \"id\": %r}]}\nand if it "
                    "needs a frame, box THAT: boxes=[{\"around\": \"hdr:%s\", \"style\": "
                    "\"Plain\"}].\nA %s stays right for what is NOT part of the reaction "
                    "row: a lettered panel holding its own arrow, `of: \"table\"` or `side: "
                    "\"below\"`, or an optimisation figure with no product grid -- on that "
                    "last one the panel's width comes out of the scheme's own band and costs "
                    "nothing."
                    % (_kind[:-1], _id, _id, _kind[:-1]))


def preflight_cycle(spec: dict) -> None:
    check_choice(spec.get("page"), PAGE_MODES, "page")
    # THE `steps` FORM DOES NOT CHECK THIS. layout_cycle._direction is called
    # only from the `rings` path; the `steps` path copies cfg["direction"] into the
    # ring row, where `== "cw"` is simply False and the ring comes out turning the
    # other way.
    check_choice(spec.get("direction"), ("cw", "ccw"), "direction")
    check_choice(spec.get("label_pos"), ("below", "outward", "auto"), "label_pos")
    nodes = [n for n in (spec.get("nodes") or []) if isinstance(n, dict)]
    entries, rl = [], []
    for i, nd in enumerate(nodes):
        path = "nodes[%d]" % i
        if nd.get("smiles"):
            entries.append((path, nd.get("id"), nd["smiles"]))
            rl.append((path, nd))
        for j, it in enumerate(nd.get("content") or []):
            if isinstance(it, dict) and it.get("smiles"):
                entries.append(("%s.content[%d]" % (path, j), it.get("id"),
                                it["smiles"]))
                rl.append(("%s.content[%d]" % (path, j), it))
    for i, st in enumerate(spec.get("steps") or []):
        if not isinstance(st, dict):
            continue
        for side in ("in", "out"):
            for j, e in enumerate(st.get(side) or []):
                if isinstance(e, dict) and e.get("smiles"):
                    entries.append(("steps[%d].%s[%d]" % (i, side, j),
                                    e.get("id"), e["smiles"]))
                    rl.append(("steps[%d].%s[%d]" % (i, side, j), e))
    check_unique_ids(entries, "nodes / steps")
    for path, d in rl:
        _rlabel_refs(d.get("smiles"), d.get("rlabels"), "%s.rlabels" % path)
        check_highlights(d.get("highlight"), "%s.highlight" % path)
    # `edges[].ring` names a row of `rings`. Unchecked, it reaches the geometry as
    # `KeyError: 'nope'` -- true, and not a sentence that names the field or the
    # rings the figure has.
    rings = [r for r in (spec.get("rings") or []) if isinstance(r, dict)]
    if rings:
        have = {str(r.get("as")) for r in rings if r.get("as")}
        used = [("edges[%d].ring" % i, e["ring"])
                for i, e in enumerate(spec.get("edges") or [])
                if isinstance(e, dict) and e.get("ring")]
        check_references(used, have, "edges", "ring")
    # `arc_style` / a step's `style`: layout_cycle._arc_kw takes the keys it knows
    # and DROPS the rest, so `{"hed": "Full"}` draws a default arrow. The MCP schema
    # already publishes this object as closed (`_obj(CY._ARC_KEYS)`); pydantic
    # ignores extras, so the schema alone cannot refuse it.
    #
    # ARC_KEYS IS RESTATED HERE, NOT IMPORTED, and that is on purpose: this module
    # must not import an engine, or the six tools would each drag a layout module
    # in before they had validated anything. test_reject_bad_input asserts
    # the two lists are equal, so a key added to one and not the other fails there.
    # THREE PLACES, not two. The schema publishes `_obj(CY._ARC_KEYS)` at
    # `arc_style`, `steps[].style` AND `edges[].style`; checking only the first two
    # leaves `edges[0].style = {"hed": "Full"}` drawing a default arrow in silence
    # (measured).
    for field, style in ([("arc_style", spec.get("arc_style"))]
                         + [("steps[%d].style" % i, (st or {}).get("style"))
                            for i, st in enumerate(spec.get("steps") or [])
                            if isinstance(st, dict)]
                         + [("edges[%d].style" % i, (e or {}).get("style"))
                            for i, e in enumerate(spec.get("edges") or [])
                            if isinstance(e, dict)]):
        check_closed(style, ARC_KEYS, field,
                     "An unknown key here is dropped, so the arrow comes out "
                     "with the engine default and nothing says why.")


def preflight_radial(spec: dict) -> None:
    check_choice(spec.get("page"), PAGE_MODES, "page")
    products = spec.get("products") or {}
    if isinstance(products, dict):
        # Unchecked, a slot outside the eight is `KeyError: 'ZZ'` from inside the
        # geometry -- the position is the whole address of a radial product, so a
        # typo here is the most likely mistake this tool has.
        check_references([("products[%r]" % k, k) for k in products],
                         RADIAL_SLOTS, "products", "slot")
        cells = [("products[%r]" % k, v) for k, v in products.items()
                 if isinstance(v, dict)]
        check_unique_ids(_id_entries(cells), "products")
        for path, c in cells:
            _rlabel_refs(c.get("smiles"), c.get("rlabels"), "%s.rlabels" % path)
            check_highlights(c.get("highlight"), "%s.highlight" % path)
    # The FIGURE-LEVEL highlight decorates the CENTRE; a product carries its own,
    # checked in the loop above. No centre `rlabels` check, because the tool has no
    # such argument -- draw_radial_figure's placeholders are named per product.
    check_highlights(spec.get("highlight"), "highlight")


def preflight_snaking(spec: dict) -> None:
    check_choice(spec.get("page"), PAGE_MODES, "page")
    cells = _cells(spec.get("rows"))
    check_unique_ids(_id_entries(cells), "rows")
    for path, c in cells:
        check_choice(c.get("align"), ("L", "R"), "%s.align" % path,
                     "A cell is pinned to an edge with 'L' or 'R'; leave align out for the "
                     "default placement.")
        _rlabel_refs(c.get("smiles"), c.get("rlabels"), "%s.rlabels" % path)
        check_highlights(c.get("highlight"), "%s.highlight" % path)
    check_reserve(spec.get("reserve"),
                  [{str(x.get("id")).strip() for x in row
                    if isinstance(x, dict) and x.get("id")}
                   for row in (spec.get("rows") or []) if isinstance(row, list)])
    check_cell_spans(spec.get("span"),
                     [row for row in (spec.get("rows") or [])
                      if isinstance(row, list)], arrows_between=True)
