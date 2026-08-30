# PetriAnalysis 0.2.5

- **`gospn mark -t` is crosschecked too.** That flag runs a different search on gospn's
  side -- it vanishes the immediate markings during the reachability walk rather than
  leaving them for the elimination -- so it reaches the generator by a third route.
  `spnp_example2` is now compared from both files, and the assembly was generalised for
  it: which blocks a file contains depends on the search (`-t` has no `I0I0I` and gains
  a direct `G0A0E`), so blocks are read as present-or-zero rather than by a fixed list.

  All three agree to 1e-12: gospn's plain output, gospn's `-t` output, and this
  package's own elimination. Perturbing a rate fails both gospn comparisons.

# PetriAnalysis 0.2.4

- **The MRSPN crosscheck now covers the distributions, not just the block structure.**
  gospn 0.22.0 records which general transition each `P<k>` block is and what governs
  each regeneration group (`gentrans`, `groupgen`), which is what the previous version of
  this test asked for: it compared every block and still passed when `det(5)` was changed
  to `det(99)`. Both are now compared, and both perturbations fail.

- A fourth fixture, `two_gen`, has **two** general transitions with different
  distributions — the case where the `P<k>` numbering has something to get wrong
  (`P0` is `Trebuild` in gospn's `raid6.spn` and `Trecon` in `raid10.spn`). Swapping the
  two distributions fails the test.

  Requires fixtures from gospn 0.22.0 or later; `test/data/README.md` says so.

# PetriAnalysis 0.2.3

- **The MRSPN block matrices are now checked against gospn's** as well as the generator.
  `test_gospn_crosscheck.jl` reads a `.npz` from `gospn mark` and compares every
  exponential, immediate and general block, for every pair of regeneration groups, on
  three nets: the textbook fail/repair, one whose general transition stays enabled across
  an EXP firing (so its group spans two markings), and one with a vanishing state between
  the general and exponential transitions. Initial vectors are compared too. They agree
  to 1e-12, and **no semantic difference was found**.

  Neither the group order nor the marking order agrees between the two implementations,
  so groups are matched by their *set* of markings and rows by the marking vector. The
  test also asserts that every block gospn wrote was one the comparison looked at:
  without that, a name wrong on both sides at once passes, since an absent block reads
  as zero.

  **Known limitation, and it is gospn's file rather than this test**: a general block is
  a 0/1 jump matrix, and the result file records neither which general transition each
  `P<k>` is nor the distribution governing each group. Changing `det(5)` to `det(99)`
  leaves every matrix in the file identical. The fixtures therefore have exactly one
  general transition, and what is verified is the *structure* of the regenerative
  process.

# PetriAnalysis 0.2.2

- **The generator is now checked against gospn's.** The two implement the same
  semantics twice, in two languages, and nothing compared them.
  `test/test_gospn_crosscheck.jl` reads a `.npz` written by `gospn mark` and checks
  that the generator it implies is the one this package builds, for an SPN
  (`spnp_example1`) and a GSPN whose vanishing states have to be eliminated
  (`spnp_example2`). They agree to 1e-12.

  The two enumerate the state space in different orders, so nothing can be compared
  position by position: every comparison is keyed on the marking vector, using the
  `place` and `mark<G>` elements gospn writes from its 0.20.0. The fixtures and the
  command that produced them are in `test/data/`.

  NPZ is a test-only dependency; the package itself gains none.

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
