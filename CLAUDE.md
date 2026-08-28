# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Run the full test suite (from the PetriTools workspace root)
julia --project=PetriAnalysis.jl -e 'using Pkg; Pkg.test()'

# Run a single test file interactively
julia --project=PetriAnalysis.jl PetriAnalysis.jl/test/test_generator.jl

# Instantiate dependencies (first time / after Project.toml changes)
julia --project=PetriAnalysis.jl -e 'using Pkg; Pkg.instantiate()'

# Run the examples
julia --project=PetriAnalysis.jl PetriAnalysis.jl/examples/mm1_queue.jl
julia --project=PetriAnalysis.jl PetriAnalysis.jl/examples/gspn_example.jl
```

This package depends on the sibling `PetriStructure.jl` **by relative path**,
declared in `Project.toml` via a `[sources]` entry (`path = "../PetriStructure.jl"`,
Julia ≥ 1.11). On older Julia, `Pkg.develop(path="../PetriStructure.jl")`.

## Architecture

`PetriAnalysis.jl` adds a state-space layer on top of `PetriStructure.jl`, which
is structural-only. It ports the marking-graph construction of the Go tool
`gospn` (`../gospn/pkg/petrinet/dfs.go`, `genvec.go`, `markinggraph.go`) and
covers SPN, GSPN **and MRSPN**:

- **SPN/GSPN** (immediate + exponential): the marking graph reduces to a CTMC
  after eliminating immediate ("vanishing") markings; `generator` returns `Q`.
- **MRSPN** (with general `GenTrans` transitions): the process is a Markov
  regenerative process. `mrspn` groups markings by GenVec into regeneration
  classes and returns the block matrices (it does not solve the MRGP).

Four source files:

- **`reachability.jl`** — `reachability_graph(pn)` runs a depth-first search from
  `initial(pn)`, interning markings in a `Dict{Vector{Int},Int}` (Julia hashes
  vectors by value, so no separate key has to be built). Each marking is
  classified into the `@enum StateType` (`VANISHING` / `TANGIBLE` / `ABSORBING`),
  its GenVec (`@enum GenStatus` per general transition, via `gen_status`) is
  recorded, and every firing is a `MarkEdge` (`:imm` / `:exp` / `:gen`). The search
  **reuses PetriStructure's token game** (`isenabled`, `fire`) and fires
  immediate, exponential **and** general transitions. Result is a `MarkingGraph`.

  The loop is written to allocate as little as it can, because it runs per
  transition per marking: `isenabled`/`fire` rather than
  `enablefunc`/`firingfunc`, which build a closure per call; the enabled
  transitions are **visited, not collected** (`filter` allocated an array and a
  closure at every state); and a GenVec is built only for a net that has general
  transitions, so an SPN or GSPN does not allocate an empty vector per marking.
  Keep it that way — see the typing design point in `PetriStructure.jl/CLAUDE.md`,
  which is where most of the cost used to be.

- **`generator.jl`** — SPN/GSPN path. `exp_rate_matrix` (timed rates) and
  `imm_prob_matrix` (weight-normalised immediate branches) are the building
  blocks. `generator(mg)` returns `(Q, states)`: a plain generator for a pure SPN;
  with vanishing markings it performs the GSPN reduction
  `R_eff = R_TT + R_TV·(I − P_VV)⁻¹·P_VT`, returning `Q` over tangible states only.
  It **errors on a net with GEN transitions** (`has_gen(mg)`), directing to `mrspn`.

- **`mrspn.jl`** — MRSPN path (port of gospn's `genvec.go` + `markinggraph.go`
  `TransMatrix`). `mrspn(mg)` partitions markings into regeneration groups
  (`MRGroup`) keyed by `(StateType, GenVec)` and assembles three sparse block
  families keyed by group pair: `expblocks` (subordinated-CTMC generators — self
  blocks carry `−(total exp exit rate)` on the diagonal), `immblocks` (branch
  probabilities), `genblocks` (0/1 regeneration jumps, keyed also by gen `trid`).
  `active_gens(an, g)` returns the aging GEN transitions (with `.dist`).

- **`dot.jl`** — `markgraph_todot(mg)` emits Graphviz DOT (style by `StateType`,
  `:gen` edges labelled `gen`), matching `PetriStructure.jl/src/dot.jl`.

## Key design points

- **GSPN semantics: immediate suppresses timed.** In a marking where any
  immediate transition is enabled, only immediate transitions fire and the
  marking is `VANISHING` — exactly gospn's `visitImmMark`-before-`visitGenMark`
  ordering. Exponential/general edges are not generated from vanishing markings.
- **GEN status / preemption.** `gen_status(tr, m)` mirrors gospn's GEN rule: a
  general transition is `GDISABLE` if an input arc is under-supplied, else
  `GENABLE` if its guard holds, else `GDISABLE` under `:prd` / `GPREEMPT` under
  `:prs`/`:pri`. PetriStructure has no inhibitor arcs, so only arcs and guards
  count. The GenVec (these statuses over all GEN transitions) is the regeneration
  class key; a `:gen` firing is a regeneration jump.
- **Constant rates/weights.** Unlike gospn (marking-dependent `ratefunc`
  closures), `PetriStructure` transition rates (`ExpTrans.rate`) and weights
  (`ImmTrans.weight`) are constant `Float64` fields, so edge values are read
  directly. There is also no transition `priority`, so gospn's priority cutoff is
  dropped.
- **Capacity bounds make the graph finite.** A firing whose result violates a
  place domain (`minmark .<= m .<= maxmark`) is treated as disabled. Unbounded
  behaviour must be modelled by a large place `max`; otherwise `maxstates`
  guards against runaway enumeration.
- **Generator is over tangible states.** `generator(...; tangible=false)` errors
  when vanishing markings exist — a CTMC generator is only defined once they are
  eliminated. A timeless trap (immediate cycle with no timed exit) makes
  `(I − P_VV)` singular and raises an informative error.
- **Returns `(Q, states)`.** `states` maps each row of `Q` back to `mg.states`,
  since DFS discovery order (not marking order) determines indexing.

## Dependencies

- `PetriStructure` (≥ 1.3) — net types and the `isenabled`/`fire` token game
- `SparseArrays`, `LinearAlgebra` (stdlib) — matrices and the elimination solve
