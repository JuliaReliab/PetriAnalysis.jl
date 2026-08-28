# PetriAnalysis 0.2.1

- `reachability_graph` is substantially faster and allocates less: about 6x the
  speed and a third of the memory on a closed ring net of 38,760 states, mostly
  from the concrete types added in PetriStructure 1.3.0. It now also
  - calls `isenabled`/`fire` instead of `enablefunc`/`firingfunc`, which built a
    closure per transition per marking;
  - visits the enabled transitions instead of collecting them with `filter`,
    which allocated an array and a closure at every state;
  - skips building a GenVec for a net with no general transitions, rather than
    allocating an empty vector per marking.
  The marking graph, generator and MRSPN blocks are unchanged.
- Declare `[compat]` for `PetriStructure`, `SparseArrays` and `LinearAlgebra`.
  Only `julia` was bounded, so nothing recorded which versions of the
  dependencies this package is known to work with.

# PetriAnalysis 0.2.0

- **MRSPN support.** `reachability_graph` is now GEN-aware: it fires general
  (`GenTrans`) transitions as regeneration jumps and records each marking's
  GenVec (the `GenStatus` — `GENABLE`/`GDISABLE`/`GPREEMPT` — of every general
  transition, following gospn's preemption rule). Earlier versions silently
  ignored GEN transitions, which dropped reachable states.
- `mrspn(mg | pn)` — Markov regenerative decomposition (gospn `mark` style):
  groups markings into regeneration classes by `(StateType, GenVec)` and emits
  sparse block matrices — subordinated-CTMC generator blocks (`exp_block`),
  immediate branch-probability blocks (`imm_block`), and general-transition jump
  blocks (`gen_block`) — plus `active_gens` (the aging GEN + its distribution per
  group), `group_markings`, and `initial_vectors`. The MRGP numeric solve is left
  to the caller.
- `generator` now refuses a net containing GEN transitions, directing the user to
  `mrspn`; `has_gen(mg)` reports whether a graph is an MRSPN.

# PetriAnalysis 0.1.0

- initial release
- `reachability_graph(pn)` — depth-first enumeration of reachable markings,
  classifying each as tangible / vanishing / absorbing (GSPN semantics: an
  immediate transition suppresses timed transitions); place capacities bound the
  graph
- `generator(mg | pn; tangible=true)` — CTMC infinitesimal generator `Q` as a
  sparse matrix, with vanishing (immediate) markings eliminated by the standard
  GSPN reduction; timeless traps are detected
- `exp_rate_matrix`, `imm_prob_matrix` — underlying sparse rate / branching matrices
- `markgraph_todot(mg)` — Graphviz DOT export of the marking graph
