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
