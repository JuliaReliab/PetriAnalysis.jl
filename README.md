# PetriAnalysis.jl

A Julia package that builds the **reachability (marking) graph** and the
**CTMC infinitesimal generator** of a Petri net defined with
[`PetriStructure.jl`](../PetriStructure.jl).

`PetriStructure.jl` is intentionally structural-only (incidence matrix,
P/T-invariants) and does not enumerate states. `PetriAnalysis.jl` adds that
state-space layer. The algorithm follows the Go tool `gospn`'s marking-graph construction:

- **SPN/GSPN** (immediate + exponential transitions) — the marking graph reduces
  to a continuous-time Markov chain once immediate ("vanishing") markings are
  eliminated; `generator` returns the CTMC generator `Q`.
- **MRSPN** (with general `GenTrans` transitions) — the process is a Markov
  regenerative process, not a CTMC. `mrspn` groups the markings into regeneration
  classes (by GenVec) and returns the block matrices that feed an MRGP solve.

## What it does

- **Reachability graph** — depth-first enumeration of all reachable markings,
  each classified as:
  - **tangible** — only timed (exponential) transitions enabled,
  - **vanishing** — at least one immediate transition enabled (left in zero time;
    GSPN semantics suppress timed transitions here),
  - **absorbing** — no transition enabled.
- **CTMC generator `Q`** (SPN/GSPN) — a `SparseArrays.SparseMatrixCSC{Float64}`
  whose rows sum to zero. Vanishing markings are removed by the standard GSPN
  reduction, folding immediate branch probabilities (transition weights) into the
  effective timed rates.
- **MRSPN block matrices** (nets with general transitions) — markings are grouped
  into regeneration classes by GenVec, and the dynamics are returned as sparse
  blocks: subordinated-CTMC generators (`exp_block`), immediate branch
  probabilities (`imm_block`), and general-transition jumps (`gen_block`), with
  each group's aging GEN distribution available via `active_gens`.
- **DOT export** of the marking graph for Graphviz.

It is intentionally *not* a solver: it returns `Q` / the MRGP blocks, leaving the
steady-state / transient solve to the caller (`LinearAlgebra`, etc.). MATLAB
`.mat` export and JSPNL `.spn` input are out of scope.

## Installation

This package depends on `PetriStructure.jl` by relative path (declared in
`Project.toml` via a `[sources]` entry). From the `PetriTools` workspace:

```julia
using Pkg
Pkg.activate("PetriAnalysis.jl")
Pkg.instantiate()
```

On Julia < 1.11 (no `[sources]` support), develop the dependency explicitly:

```julia
Pkg.develop(path = "PetriStructure.jl")
```

## Quick start

```julia
using PetriStructure
using PetriAnalysis

# M/M/1/K queue: a buffer of capacity K, arrival and service transitions.
K = 4
pn = petri()
buf = place(pn, "buf", 0, K)
arr = exptrans(pn, "arrival", 0.6)
srv = exptrans(pn, "service", 1.0)
arc(pn, arr, buf)   # arrival produces a token (capacity K bounds the graph)
arc(pn, buf, srv)   # service consumes a token

mg = reachability_graph(pn)
# MarkingGraph(5 states: 5 tangible, 0 vanishing, 0 absorbing; 8 edges)

Q, states = generator(mg)   # SparseMatrixCSC; states maps rows -> mg.states
```

For a GSPN with immediate transitions, the generator is built over the tangible
markings only:

```julia
Q, states = generator(mg; tangible = true)   # vanishing markings eliminated
print(markgraph_todot(mg))                    # Graphviz DOT of the full graph
```

For an **MRSPN** — a net with general (`gentrans`) transitions — `generator` does
not apply; use `mrspn`:

```julia
# failure (EXP) / deterministic repair (GEN): a Markov regenerative process
pn = petri()
up   = place(pn, "up", 1, 1)
down = place(pn, "down", 0, 1)
Tf = exptrans(pn, "Tfail", 0.1)
Tr = gentrans(pn, "Trepair", detdist(5.0))   # deterministic repair time
arc(pn, up, Tf);   arc(pn, Tf, down)
arc(pn, down, Tr); arc(pn, Tr, up)

an = mrspn(pn)
for g in an.groups
    @show g.label, group_markings(an, g.id), active_gens(an, g.id)
end
an.expblocks   # subordinated-CTMC generator blocks, keyed (src_group, dst_group)
an.genblocks   # regeneration jumps, keyed (src_group, dst_group, gen_trid)
```

## API

| Function | Description |
|---|---|
| `reachability_graph(pn; maxstates=10_000_000)` | Build the full GEN-aware `MarkingGraph` |
| `nstates(mg)` | Number of reachable markings |
| `tangible_states(mg)` / `vanishing_states(mg)` | State indices by class |
| `has_gen(mg)` | Whether the net is an MRSPN (contains GEN firings) |
| `generator(mg; tangible=true)` | CTMC generator `Q` and its `states` map (SPN/GSPN) |
| `generator(pn; tangible=true)` | Convenience: graph + generator in one call |
| `exp_rate_matrix(mg)` / `imm_prob_matrix(mg)` | Sparse rate / immediate-branch matrices |
| `mrspn(mg)` / `mrspn(pn)` | MRSPN regeneration-block decomposition (`MRSPNGraph`) |
| `ngroups(an)` / `group_markings(an, g)` | Regeneration groups and their markings |
| `active_gens(an, g)` | Aging GEN transitions (with `.dist`) of group `g` |
| `exp_block` / `imm_block` / `gen_block` | Block matrices of the MRSPN decomposition |
| `initial_vectors(an)` | Per-group initial probability vectors |
| `markgraph_todot(mg)` | Graphviz DOT string of the marking graph |

`generator` returns `(Q, states)` where `mg.states[states[k]]` is the marking of
row `k` of `Q`. For `mrspn`, a group's `marks` order is the local row/column
index of its blocks.

### Notes

- **Bounded nets only.** A firing that would push a place beyond its `max` is
  treated as disabled, which keeps the graph finite. Model an unbounded place by
  giving it a sufficiently large `max`.
- **Timeless traps.** An immediate cycle with no timed exit makes the reduction
  singular; `generator` raises an informative error.

## Examples

See [`examples/`](examples/):

| File | Topic |
|---|---|
| `mm1_queue.jl` | Pure SPN; generator vs analytic M/M/1/K distribution |
| `gspn_example.jl` | GSPN; vanishing-state elimination and DOT export |
| `mrspn_example.jl` | MRSPN; regeneration groups and block matrices (fail/repair) |

## Testing

```bash
julia --project=PetriAnalysis.jl -e 'using Pkg; Pkg.test()'
```

## Dependencies

- [`PetriStructure.jl`](../PetriStructure.jl) — net definition and the
  `enablefunc` / `firingfunc` token game reused to drive the search
- `SparseArrays`, `LinearAlgebra` (standard library)

## License

MIT — see [LICENSE](LICENSE).

## Author

Hiroyuki Okamura <okamu@hiroshima-u.ac.jp>
