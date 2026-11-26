"""address.py — figure addressing: where a thing IS.

BOTH LEVELS. Level one answers "where is this object" off a box; level two
refines that to an atom by FORWARDING to the existing resolver. `layout_cycle`
fills the registry as it places (registration is the only way in).

WHY A NEW MODULE RATHER THAN A WIDER Anchors. The existing path is
SMARTS -> RDKit atom index -> CDXML node id, and every step of it presumes an
ATOM. A text node has none: `Anchors.record` returns on its first line when there
is no SMILES, and `node_points` reads `<n>` and cannot see a `<t>`. So the thing
that was missing was never a bigger atom resolver, it was a layer ABOVE atoms that
can say where a BOX is. That is this module. `anchors.py` is unchanged -- not one
line -- and serves as the second level unchanged; see `Registry.atom_point`,
which is only a forward.

WHAT THIS MODULE DOES NOT DO, all three deliberate:
  - it does not place anything. Coordinates come from the caller.
  - it does not draw. No arrow, no curve, no reservation.
  - it does not touch the other five engines, and there is no abstraction here
    reaching for them.

THE ONE IMPORT THAT IS DEFERRED, and why it matters: `anchors` is imported inside
`_anchors()`, on the first atom address and never before. Atom refinement is
lazy: a figure that writes no `atom` must not pay `node_points()`, and that path
pulls in RDKit.

THREE DESIGN DECISIONS, each a deliberate choice rather than an accident of the
implementation:

  1. `ink` on the record. The parametric form, {"of": "e3", "at": 0.5}, is "0..1
     measured along the DRAWN INK", and an arc's ink is not recoverable from its
     bounding rect -- the chord is not the arc. Hence a stored polyline rather
     than a parametric function.
  2. The trim. Left to the implementation, on the principle that reading it off
     an existing emitter beats measuring; it is `arc_margin`. See ATTACH_TRIM
     below.
  3. The reference point a paired end presents to the other end is its own box
     CENTRE. This is a KNOWN LIMITATION -- a wide box attaches from the middle of
     its edge where the eye wants the near end. Deliberate; see `Registry._ref`.

The third one is the instructive one: it is as much a choice as the other two, and
an unlisted decision buried in an implementation becomes a bug report later.

FRAME: LAYOUT, always, and negative coordinates are normal (layout_cycle puts the
ring centre at the origin). The LAYOUT / PAGE / CONTENT discipline is unchanged and
this module never converts between them -- see `_reconcile`'s docstring in
chemwand_mcp.py, which names the three frames.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from allocator import Rect

# Every kind here came out of a real published figure, so this set is closed on
# purpose: an unknown kind is a typo, and a typo that is accepted becomes an
# object with no address.
KINDS = frozenset({
    "structure",     # a drawn species
    "text",          # a written species: L~n~Pd, [Fe^III^], Cp*RhX~2~
    "marker",        # SET box, Bond Formation, e- circle, kcal/mol cloud
    "edge",          # a step arrow; the host a feeder hangs off
    "edgegroup",     # a pair of edges named together
    "ring",          # the circle itself, for a ring label or centre note
    "offring",       # a species off the circle
    "subgraph",      # a nested figure: equilibrium pair, coordination slot
})

# The only kind an atom address may refine (the atom level). Written as a set
# rather than `== "structure"` so that the error message below can list it, and so
# making `subgraph` refinable later changes one line.
ATOM_KINDS = frozenset({"structure"})

# Kinds that cannot be an endpoint. A group of edges has no single
# attachment point, so it may HOST a caption and may not terminate an arrow.
HOST_ONLY_KINDS = frozenset({"edgegroup"})

_FORMS = ("object", "atom", "param", "group")
_DICT_KEYS = frozenset({"of", "atom", "at"})


# ATTACH_TRIM: how far outside the box the attachment point sits, in page units.
#
# 12.0, and it is NOT a new number -- it is layout_cycle's `arc_margin`, which is
# already defined as exactly this quantity ("how far an arc backs off from a node
# before starting", CYCLE_DEFAULTS, and the feeder code measures it "from the edge
# of the box along" the tangent). Reading it off the emitter that already draws
# this gap is the third option in the house rule for pinning a constant -- derive
# in closed form, else read the emitter, and only measure as a fallback -- and it
# is the one that applies here, because the gap a reader sees on an existing cycle
# IS arc_margin.
#
# NOT imported from layout_cycle. This module has no engine dependency and is not
# going to grow one for a float; a test asserts the two agree, so a drift is
# caught rather than assumed away.
#
# It is ONE constant for both ends of an edge, and that is load-bearing rather than
# tidy: a per-role trim is exactly what the attachment rule forbids, and direction
# invariance is what would break.
ATTACH_TRIM = 12.0


class AddressError(ValueError):
    """A bad address, or one that names nothing.

    ValueError, so the tool layer's existing `except ValueError` paths and the
    house style of `raise ValueError("path: what went wrong")` both keep working;
    a distinct type so a test can tell an addressing failure from a chemistry one.

    Addressing failure ALWAYS raises. Comparing the existing engines is the whole
    argument -- on linear a text cell fails honestly ("cell 'T' is not in this
    figure, the ids are ...") and a reagent fails silently (accepted, then no arrow
    drawn). Never the second one.
    """


class MissingHandleError(AddressError):
    """An atom address on a structure whose placer supplied no handle.

    Its own type because the fix is not the caller's: the ADDRESS is fine, the
    placement did not pass `smiles` / `index_map` through. Distinguishing them
    matters -- "your SMARTS is wrong" and "the engine did not hand over what it
    knew" send you to different files.
    """


@dataclass(frozen=True)
class Record:
    """One house number. The object record.

    `handle` is a HANDLE, not resolved atom points, and that is deliberate: an
    object address costs nothing (the packer computed the rect anyway) while atom
    refinement must stay lazy, because `node_points()` on every figure is what the
    current engine pays for annotations nobody asked for. Nothing in this module
    looks inside it.
    """
    aid: str
    kind: str
    rect: Rect                      # LAYOUT frame; may be negative
    handle: object | None = None    # opaque; only the atom level reads it
    ink: tuple = ()                 # ((x, y), ...) LAYOUT; only `at` reads it
    where: str = ""                 # provenance, for the duplicate-id message


@dataclass(frozen=True)
class Resolved:
    """What an address came out as. One shape for all four forms, tagged by `form`.

    `records` is a tuple because the group form names several; the other three name
    exactly one, and `record` is the convenience for them.
    """
    form: str
    records: tuple
    atom: str | None = None         # the SMARTS, handed on to Anchors
    at: float | None = None

    @property
    def record(self) -> Record:
        if len(self.records) != 1:
            raise AddressError(
                "this address names %d objects, so it has no single record"
                % len(self.records))
        return self.records[0]

    @property
    def handle(self):
        """The structure handle, for atom refinement. None when it has none."""
        return self.record.handle


def _num(v, what: str, where: str) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        raise AddressError("%s: %s is a number, got %r" % (where, what, v))


def parse(addr, where: str = "?") -> tuple:
    """Shape-check one address. -> (form, ids, atom, at). No registry involved.

    The four address forms and nothing else. Separate from lookup so that "this is
    not an address" and "this address names nothing here" are two different
    messages -- the caller fixes them differently.

    `where` is the recipe path, e.g. "edges[2].to". It goes into every message,
    in the style the KEYS ARE CHECKED errors already use.
    """
    if isinstance(addr, str):
        aid = addr.strip()
        if not aid:
            raise AddressError(
                "%s: an address cannot be the empty string. Give an id, or one of "
                '{"of": id, "atom": SMARTS} / {"of": edge, "at": 0..1} / '
                '{"of": [id, id]}.' % where)
        return "object", (aid,), None, None

    if not isinstance(addr, dict):
        raise AddressError(
            "%s: an address is an id string or a dict with 'of', got %s (%r)"
            % (where, type(addr).__name__, addr))

    unknown = set(addr) - _DICT_KEYS
    if unknown:
        raise AddressError(
            "%s: %s is not an address key. Known: %s."
            % (where, ", ".join(sorted(map(repr, unknown))),
               ", ".join(sorted(_DICT_KEYS))))
    if "of" not in addr:
        raise AddressError(
            "%s: an address dict needs 'of' (the object it is about), got keys %s."
            % (where, ", ".join(sorted(map(repr, addr))) or "none"))
    if "atom" in addr and "at" in addr:
        raise AddressError(
            "%s: 'atom' and 'at' are mutually exclusive -- 'atom' names an atom "
            "INSIDE an object, 'at' names a point ALONG an edge's ink." % where)

    of = addr["of"]
    if isinstance(of, (list, tuple)):
        # The group form. A refinement of a group is not a thing -- a
        # group has no single ink and no single structure.
        if "atom" in addr or "at" in addr:
            raise AddressError(
                "%s: 'of' is a list (an edge group), so it takes neither 'atom' "
                "nor 'at' -- a group has no single atom and no single ink. Name "
                "one member instead." % where)
        ids = tuple(str(v).strip() for v in of)
        if not ids:
            raise AddressError("%s: 'of' is an empty list; name at least one "
                               "object." % where)
        if any(not i for i in ids):
            raise AddressError("%s: 'of' holds an empty id: %r" % (where, list(of)))
        if len(set(ids)) != len(ids):
            raise AddressError(
                "%s: 'of' repeats an id (%s). A group is a set of distinct "
                "objects." % (where, ", ".join(sorted(
                    {i for i in ids if ids.count(i) > 1}))))
        return "group", ids, None, None

    if not isinstance(of, str) or not of.strip():
        raise AddressError(
            "%s: 'of' is an id string or a list of them, got %r" % (where, of))
    aid = of.strip()

    if "atom" in addr:
        sma = addr["atom"]
        if not isinstance(sma, str) or not sma.strip():
            raise AddressError(
                "%s: 'atom' is a SMARTS string, got %r" % (where, sma))
        return "atom", (aid,), sma.strip(), None

    if "at" in addr:
        t = _num(addr["at"], "'at'", where)
        if not (0.0 <= t <= 1.0):
            raise AddressError(
                "%s: 'at' is 0..1 along the edge's drawn ink, got %g." % (where, t))
        return "param", (aid,), None, t

    return "object", (aid,), None, None


def _polyline_point(ink, t: float, aid: str, where: str) -> tuple:
    """The point at fraction `t` of a polyline's ARC LENGTH.

    Arc length rather than vertex count, because that is what "measured along the
    drawn ink" means and because curve_polyline samples a cubic at a fixed count
    per segment, not at a fixed spacing -- so the vertices are not evenly spread
    and counting them would put `at: 0.5` somewhere that is not the middle.
    """
    pts = [(float(x), float(y)) for x, y in ink]
    if len(pts) < 2:
        raise AddressError(
            "%s: the ink of %r is %d point(s); a point along it needs at least 2."
            % (where, aid, len(pts)))
    seg = [math.dist(pts[i], pts[i + 1]) for i in range(len(pts) - 1)]
    total = math.fsum(seg)
    if total <= 0.0:
        raise AddressError(
            "%s: the ink of %r has zero length, so 'at' cannot name a point on it."
            % (where, aid))
    want = t * total
    run = 0.0
    for i, d in enumerate(seg):
        if run + d >= want or i == len(seg) - 1:
            f = 0.0 if d <= 0.0 else (want - run) / d
            f = min(1.0, max(0.0, f))
            (x0, y0), (x1, y1) = pts[i], pts[i + 1]
            return (x0 + (x1 - x0) * f, y0 + (y1 - y0) * f)
        run += d
    return pts[-1]                                      # unreachable


def _on_boundary(rect: Rect, toward: tuple) -> tuple:
    """Where the ray from `rect`'s centre to `toward` leaves the rect. -> (x, y, ux, uy)

    The attachment rule, and the whole of it: WHICH EDGE is decided by the incoming
    direction, not by whether this end is the source or the target. There is no
    role parameter in this function and there is not going to be one -- that
    absence is what makes direction invariance true by construction rather
    than by a matching pair of adjustments.
    """
    cx, cy = rect.cx, rect.cy
    dx, dy = float(toward[0]) - cx, float(toward[1]) - cy
    d = math.hypot(dx, dy)
    if d <= 0.0:
        raise AddressError(
            "the other end of this connection is at the centre of the box it is "
            "connecting to, so there is no direction to attach from. Two objects "
            "cannot be addressed against each other at the same point.")
    ux, uy = dx / d, dy / d
    hw, hh = rect.w / 2.0, rect.h / 2.0
    # Scale the direction until it lands on the nearer of the two edge pairs. A
    # zero component means that pair is never reached, hence the infinities.
    sx = hw / abs(dx) if dx else math.inf
    sy = hh / abs(dy) if dy else math.inf
    s = min(sx, sy)
    if s is math.inf:                                   # zero-area box, no edges
        return cx, cy, ux, uy
    return cx + dx * s, cy + dy * s, ux, uy


@dataclass
class Registry:
    """Every house number in one figure. One flat global namespace.

    FLAT, including the insides of a nested subgraph, and that is a simplification
    rather than a shortcut: with ids globally unique, "how do I address something
    inside a subgraph" stops being a question -- no path syntax, no scope
    resolution, no shadowing. The cost is that the caller must not repeat an id,
    and a repeat is caught here naming both origins.

    An id is an INSTANCE, not a chemical identity. The same compound drawn
    twice is two ids. That is the opposite of `_register`'s rule, where one id
    means one molecule built once -- construction may be shared, an address may
    not -- so this registry deliberately does not consult, extend or resemble
    `structures`.

    REGISTRATION IS THE ONLY WAY IN, and that claim is worth stating
    precisely rather than broadly, because one half of it is enforceable in Python
    and the other half is not:

      CLOSED. `_by_id` is `init=False`, so there is no constructor argument that
              takes a filled table -- `Registry({...})` would be a way to install
              records that never saw `register`. It is also
              private, and `ids()` hands back a sorted tuple rather than the dict,
              so there is no writable handle to it either.
      NOT CLOSED, and not worth forcing: `Record(...)` is a public dataclass, so
              anyone can build one, and `reg._by_id[k] = rec` reaches the store
              through the underscore. Python has no way to prevent that which is
              not either trivially defeated or a nuisance to the engine that has
              to use this module.

    The second one does not weaken the rule, because the rule is not about hostile
    or accidental hand-construction. It is about ENGINE-INTERNAL PLACEMENT drifting
    out of step with registration -- the layout_linear hole, where two
    `place_with` call sites existed and only one recorded. That failure needs no
    smuggled Record: it happens when a placement path simply never calls the
    registrar. What defends against it is the engine calling one wrapper that
    places AND registers in the same operation; this class is what that wrapper
    writes into.

    So: the door a normal caller could walk through by accident is shut, the door
    that needs an underscore is not, and neither is the door the rule is about.

    ONE REGISTRY PER FIGURE. NOT PER SESSION, NOT PER ENGINE, NOT REUSED.

    The flat namespace is flat WITHIN one figure. A Registry handed to two
    renders raises on the second, because the second figure re-registers `node#0`
    and every other name the first one took -- and the message will say "duplicate
    id 'node#0'" and name two provenances that look identical, which is a confusing
    way to find out. Make a new one, or let `render_cycle` make it (it does, when
    the caller passes none).

    This is written down because for a multi-ring figure, reusing one registry
    across the rings of one page looks like the obvious thing to do -- and it is
    exactly the thing that breaks. If several rings ever need to share a
    namespace, they have to share ONE render, not one Registry passed twice;
    otherwise the second ring's `node#0` is the first ring's. A test pins the
    raise so nobody discovers it in a figure.
    """
    _by_id: dict = field(default_factory=dict, init=False)
    # The atom-level resolver, None until the first atom address asks for it.
    # init=False for the same reason _by_id is: there is no constructor argument
    # that installs a pre-loaded one.
    _anch: object = field(default=None, init=False)

    # ---- filling it in ---------------------------------------------------

    def register(self, aid, kind: str, rect: Rect, *,
                 handle: object | None = None, ink=(),
                 where: str = "") -> Record:
        """File one object. Raises on a repeat or an unknown kind."""
        aid = str(aid).strip()
        if not aid:
            raise AddressError(
                "%s: an object needs a non-empty id; ids are how everything else "
                "in the figure refers to it." % (where or "register"))
        if kind not in KINDS:
            raise AddressError(
                "%s: %r is not an object kind. Known: %s."
                % (where or ("object %r" % aid), kind, ", ".join(sorted(KINDS))))
        if not isinstance(rect, Rect):
            raise AddressError(
                "%s: rect is an allocator Rect in the LAYOUT frame, got %s."
                % (where or ("object %r" % aid), type(rect).__name__))
        old = self._by_id.get(aid)
        if old is not None:
            # BOTH origins named, which is the whole point of keeping `where`: one
            # of the two is the one the author meant to rename, and a message that
            # names only the loser cannot say which.
            raise AddressError(
                "duplicate id %r. Ids are unique across the WHOLE figure, including "
                "inside a nested subgraph. It is already the %s registered at %s, and is "
                "being registered again as a %s at %s. The same compound drawn twice "
                "needs two ids -- an id is an instance, not a chemical identity."
                % (aid, old.kind, old.where or "(no origin recorded)",
                   kind, where or "(no origin recorded)"))
        rec = Record(aid=aid, kind=kind, rect=rect, handle=handle,
                     ink=tuple((float(x), float(y)) for x, y in (ink or ())),
                     where=where)
        self._by_id[aid] = rec
        return rec

    def grow(self, aid, rect: Rect, where: str = "?") -> Record:
        """Widen an existing record's rect to include `rect`. Cannot create one.

        WHY THIS EXISTS AND WHY IT IS NOT A SECOND DOOR. A node's attachment
        box is the STRUCTURE union its LABEL, and the label's position is
        only known after it has been searched for -- so the record is filed when the
        structure is placed (placement and registration stay ONE operation) and
        widened when the label lands. The only-way-in rule is about a PLACEMENT
        that never registers; this creates nothing, so it cannot be that.

        Refuses an unknown id for the same reason: a `grow` that quietly registered
        would be exactly the door that rule closes.

        Only ever widens. A rect already inside the record is a no-op rather than a
        shrink, so the order two calls arrive in cannot change the answer.
        """
        rec = self._find(str(aid).strip(), where)
        self._by_id[rec.aid] = Record(
            aid=rec.aid, kind=rec.kind, rect=rec.rect.union(rect),
            handle=rec.handle, ink=rec.ink, where=rec.where)
        return self._by_id[rec.aid]

    # ---- reading it ------------------------------------------------------

    def __contains__(self, aid) -> bool:
        return str(aid).strip() in self._by_id

    def __len__(self) -> int:
        return len(self._by_id)

    def ids(self) -> tuple:
        return tuple(sorted(self._by_id))

    def __getitem__(self, aid) -> Record:
        return self._by_id[str(aid).strip()]

    def _find(self, aid: str, where: str) -> Record:
        rec = self._by_id.get(aid)
        if rec is None:
            raise AddressError(
                "%s: %r is not an object in this figure. The ids are %s."
                % (where, aid, ", ".join(map(repr, self.ids())) or "(none)"))
        return rec

    # ---- shape, existence, attachment -------------------------------------

    def resolve(self, addr, where: str = "?") -> Resolved:
        """Address -> Resolved. Shape then existence, then per-form legality.

        Does NOT resolve atoms and does NOT compute a point: the first is the atom
        level (`atom_point`), the second is `attach`.
        """
        form, ids, atom, at = parse(addr, where)
        recs = tuple(self._find(i, where) for i in ids)

        if form == "atom":
            # A text node has no atoms; that is a definition, not a
            # missing feature, so it is refused rather than approximated by the
            # box centre.
            if recs[0].kind not in ATOM_KINDS:
                raise AddressError(
                    "%s: 'atom' refines an object that HAS atoms, and %r is a %s. "
                    "Only %s takes an 'atom'. A written species has no atoms to "
                    "name -- address the object itself and the attachment point "
                    "comes off its box."
                    % (where, recs[0].aid, recs[0].kind,
                       "/".join(sorted(ATOM_KINDS))))
            h = recs[0].handle
            if h is None:
                raise MissingHandleError(
                    "%s: %r is a structure but was registered without a structure "
                    "handle, so there is nothing to resolve %r against. Whatever "
                    "placed it did not pass one."
                    % (where, recs[0].aid, atom))
            missing = [k for k in ("placed", "smiles", "index_map")
                       if not h.get(k)]
            if missing:
                raise MissingHandleError(
                    "%s: %r has a handle with no %s, so %r cannot be resolved. "
                    "The engine placed the structure but did not pass what the "
                    "atom lookup needs -- on the cycle engine that comes from "
                    "recipe['anchor_src'], which the MCP layer fills in."
                    % (where, recs[0].aid, " / ".join(missing), atom))

        if form == "param" and not recs[0].ink:
            raise AddressError(
                "%s: 'at' measures along an object's DRAWN INK and %r (a %s) was "
                "registered without any. Only an object placed with its ink -- an "
                "edge -- can take an 'at'."
                % (where, recs[0].aid, recs[0].kind))

        return Resolved(form=form, records=recs, atom=atom, at=at)

    # ---- the second level, forwarded --------------------------------------

    def _anchors(self):
        """The lazily-built Anchors, made on the FIRST atom address and no sooner.

        THE LAZINESS IS THIS METHOD. `Anchors.record` calls `C.node_points()`, which
        regexes a fragment string; doing that for every structure of every figure
        is exactly the cost layout_linear avoids by accident (its tool layer only
        fills `anchor_src` when annotations are asked for). Here the object
        addresses are built eagerly because they are free, and this -- the
        expensive half -- waits to be asked. A figure that never writes an `atom`
        never imports anchors and never parses a fragment.

        DEFERRED IMPORT for the same reason layout_linear defers its own: anchors
        pulls in decorate and smiles_to_frag, and smiles_to_frag pulls in RDKit.
        """
        if self._anch is None:
            import anchors as _AN         # deferred: RDKit lives down this path
            self._anch = _AN.Anchors()
        return self._anch

    def atom_point(self, addr, where: str = "?"):
        """(x, y) of the atom an `{"of": .., "atom": ..}` address names.

        THE FORWARD, and it is only a forward: the SMARTS and the handle go to the
        existing `Anchors` / `match_node_ids`, which are not modified by one line --
        that level does not need changing. What was missing was never the atom
        resolver; it was a layer above it that could say where a BOX is.

        The mean of the first match, via `Anchors.centre`: a whole-ring SMARTS
        arrives at the ring centre, which is what "attacks the ring" means, and a
        single-atom one lands on that atom.

        A MISS RAISES rather than returning None. Anchors reports it through
        S.ANCHOR_MISSES, which is right for a highlight -- a missing decoration is
        cosmetic -- and wrong here: an address that resolves to nothing would draw
        a figure with a piece missing, and the always-raise rule (see AddressError)
        is exactly about not doing that.
        """
        got = self.resolve(addr, where)
        if got.form != "atom":
            raise AddressError(
                "%s: atom_point needs an address with an 'atom', got the %s form."
                % (where, got.form))
        rec = got.record
        an = self._anchors()
        if rec.aid not in an:
            h = rec.handle
            # Keyed by the ADDRESS id, not the structures key: two instances of one
            # compound are two addresses and must not share a cache slot --
            # they sit at different places on the ring.
            an.record(rec.aid, h["placed"], h["base"],
                      smiles=h["smiles"], index_map=h["index_map"],
                      base_id=h.get("base_id", 1))
        pt = an.centre(rec.aid, got.atom, where=where)
        if pt is None:
            raise AddressError(
                "%s: %r matched nothing in %r. An address that resolves to nothing would "
                "leave a hole in the drawing, so it raises rather than being skipped."
                % (where, got.atom, rec.handle.get("s") or rec.aid))
        return pt

    def host_box(self, addr, where: str = "?") -> Rect:
        """The box an address offers as a HOST -- for a caption, a name, a feeder.

        Defined for every form including the group (its actual use): a group's
        host box is the union of its members. Note this is a box, not a point --
        hosting and terminating are different questions and a group answers only
        the first.
        """
        got = self.resolve(addr, where)
        box = got.records[0].rect
        for r in got.records[1:]:
            box = box.union(r.rect)
        return box

    def attach(self, addr, toward, where: str = "?") -> tuple:
        """(x, y) in LAYOUT: where a connection from `toward` meets this object.

        `toward` is the OTHER end's reference point, and passing a point rather
        than a role is the mechanism for direction invariance: flipping an edge
        swaps which call is made first and changes neither call's arguments, so
        both coordinates come out bit-identical.

        The point sits ATTACH_TRIM outside the box, on the edge the incoming
        direction actually hits.
        """
        got = self.resolve(addr, where)

        if got.form == "group":
            # A group has no single attachment point, stated as the reason
            # rather than as a refusal.
            raise AddressError(
                "%s: an edge group cannot be an endpoint -- a set of edges has no "
                "single attachment point. A group hosts a caption, a name or a "
                "feeder; name one of its members (%s) to terminate a connection."
                % (where, ", ".join(map(repr, [r.aid for r in got.records]))))

        rec = got.record
        if rec.kind in HOST_ONLY_KINDS:
            raise AddressError(
                "%s: %r is a %s, which hosts things and cannot be an endpoint." % (where, rec.aid, rec.kind))

        if got.form == "atom":
            # SECOND LEVEL. The point IS the atom -- no trim, for the same
            # reason a parametric point is not trimmed: an electron-pushing arrow's
            # tail sits ON the lone pair, which is what the verified behaviour
            # already does (measured at 0.0000 pt from the target atom).
            return self.atom_point(addr, where)

        if got.form == "param":
            # A point ON the ink, not beside a box: an edge is a line, so there is
            # no edge of a box to hit and nothing to trim. "The parametric point
            # on the edge" is the whole definition.
            return _polyline_point(rec.ink, got.at, rec.aid, where)

        x, y, ux, uy = _on_boundary(rec.rect, toward)
        return (x + ux * ATTACH_TRIM, y + uy * ATTACH_TRIM)

    def endpoints(self, a, b, *, where: str = "?") -> tuple:
        """Both ends of one connection, as ((ax, ay), (bx, by)).

        The pairing helper, and it exists so that direction invariance has
        something to test that can actually go red: each end is attached toward
        the OTHER's reference point, so the two calls are symmetric, and an
        implementation that trimmed by role would show up here as a flip changing
        the coordinates.

        The reference point is the object's box centre -- or, for a parametric
        address, the point on the ink, since an `at` point is already a point.
        """
        pa, pb = self._ref(a, where + ".from"), self._ref(b, where + ".to")
        return (self.attach(a, pb, where + ".from"),
                self.attach(b, pa, where + ".to"))

    def _ref(self, addr, where: str) -> tuple:
        """The point an address presents to the OTHER end as "where I am".

        The box CENTRE -- a decision with a known limitation recorded alongside
        it, not an unconsidered default. A wide box (published cycles have
        several: a long ketone label, a side-by-side `2a + 2b`) attaches from the
        middle of the edge the ray crosses, where the
        eye wants the end nearer the other object. Adjacent species on the ring
        are barely affected; an off-ring node or an interior chord will visibly
        skew.

        Left alone on purpose until arrows are actually drawn, at which
        point "aim at the other object's nearest boundary point" is the obvious
        thing to evaluate. Until then, skew on a wide box is KNOWN and INTENDED
        rather than a defect.

        A parametric address is already a point, so it presents itself.
        """
        got = self.resolve(addr, where)
        if got.form == "param":
            return _polyline_point(got.record.ink, got.at, got.record.aid, where)
        if got.form == "atom":
            # An atom address is already a point, so it presents itself -- exactly
            # as a parametric one does. This is also what keeps direction
            # invariance true for the refined form: neither end's reference
            # depends on the other, so flipping the edge changes nothing.
            return self.atom_point(addr, where)
        box = self.host_box(addr, where)
        return (box.cx, box.cy)
