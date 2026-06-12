"""
    StateType

Classification of a reachable marking.

- `VANISHING`  — one or more immediate (`ImmTrans`) transitions are enabled; the
  marking is left in zero time. Following GSPN semantics, timed (exponential and
  general) transitions are *not* fired from a vanishing marking.
- `TANGIBLE`   — no immediate transition is enabled, but at least one timed
  transition (exponential, or general with `GENABLE` status) is. The sojourn is
  governed by those timed transitions.
- `ABSORBING`  — no transition is enabled.
"""
@enum StateType VANISHING TANGIBLE ABSORBING

"""
    GenStatus

Per-marking status of a general (`GenTrans`) transition, mirroring `gospn`:

- `GENABLE`  — enabled (input arcs satisfied and guard holds); it is aging.
- `GDISABLE` — disabled (an input arc is under-supplied, or the guard fails under
  the `:prd` policy, which resets the age).
- `GPREEMPT` — preempted (the guard fails under `:prs`/`:pri`, which retains the
  age).

The vector of `GenStatus` over all general transitions — the *GenVec* — labels
the regeneration class a marking belongs to (see [`mrspn`](@ref)).
"""
@enum GenStatus GDISABLE GENABLE GPREEMPT

"""
    MarkEdge

A directed edge of the marking graph, produced by firing transition `trid`
(a `tr.id` in the source `PN`) from state `src` to state `dst`. `kind` is one of
`:imm`, `:exp`, `:gen`; `value` is the immediate weight, the exponential rate, or
`1.0` for a general-transition jump.
"""
struct MarkEdge
    src::Int
    dst::Int
    trid::Int
    kind::Symbol
    value::Float64
end

"""
    MarkingGraph

The reachability graph of a Petri net. Markings are interned: each distinct
marking vector appears once in `states`, and `index` maps a marking back to its
1-based position. `types[i]` classifies `states[i]`; `genvecs[i]` is its GenVec
(the `GenStatus` of every general transition, in `pn.gentrans` order); `edges`
holds every firing.
"""
struct MarkingGraph
    pn::PN
    states::Vector{Vector{Int}}
    index::Dict{Vector{Int},Int}
    types::Vector{StateType}
    genvecs::Vector{Vector{GenStatus}}
    edges::Vector{MarkEdge}
    initial_state::Int
end

"""
    gen_status(tr, m) -> GenStatus

Status of general transition `tr` in marking `m`, following `gospn`'s GEN
enabling rule (`firing.go`): disabled if any input arc is under-supplied;
otherwise enabled if the guard holds; otherwise disabled under `:prd` or
preempted under `:prs`/`:pri`. (`PetriStructure` has no inhibitor arcs, so only
arc supply and guards are considered.)
"""
function gen_status(tr, m)
    for a in tr.inarcs
        if m[a.src.id] < a.mul
            return GDISABLE
        end
    end
    guard_ok = all(g -> evaluate(g, m), tr.guard)   # empty guard list -> true
    guard_ok && return GENABLE
    return tr.policy == :prd ? GDISABLE : GPREEMPT
end

"""
    reachability_graph(pn; maxstates = 10_000_000) -> MarkingGraph

Enumerate the reachable markings of `pn` by depth-first search from the initial
marking, building the full marking graph. Immediate, exponential, **and general**
transitions are all expanded, so the graph is correct for SPN, GSPN and MRSPN.

A candidate marking produced by a firing is admitted only if it lies within every
place's domain (`minmark(pn) .<= m .<= maxmark(pn)`); a firing that would exceed a
place capacity is treated as disabled. This keeps the graph finite for bounded
nets — model an unbounded place by giving it a sufficiently large `max`.

GSPN semantics are enforced: in a marking where any immediate transition is
enabled, only immediate transitions fire (the marking is `VANISHING`). Otherwise
every enabled exponential transition and every general transition with `GENABLE`
status fires; a general-transition firing is a regeneration jump (a `:gen` edge).

Throws if more than `maxstates` distinct markings are discovered.
"""
function reachability_graph(pn::PN; maxstates::Int = 10_000_000)
    m0 = initial(pn)
    mmax = maxmark(pn)
    mmin = minmark(pn)

    states = Vector{Vector{Int}}()
    index = Dict{Vector{Int},Int}()
    types = Vector{StateType}()
    genvecs = Vector{Vector{GenStatus}}()
    edges = Vector{MarkEdge}()

    # Intern a marking, returning its index and whether it is newly seen.
    function intern!(m)
        i = get(index, m, 0)
        i != 0 && return i, false
        push!(states, m)
        push!(types, ABSORBING)            # provisional; set when the state is expanded
        push!(genvecs, GenStatus[])        # filled when the state is expanded
        i = length(states)
        index[m] = i
        return i, true
    end

    enabled(tr, m) = enablefunc(pn, tr)(m)
    indomain(m) = all(mmin[k] <= m[k] <= mmax[k] for k in eachindex(m))

    i0, _ = intern!(copy(m0))
    stack = [i0]
    while !isempty(stack)
        i = pop!(stack)
        m = states[i]
        genvecs[i] = [gen_status(tr, m) for tr in pn.gentrans]

        en_imm = filter(tr -> enabled(tr, m), pn.immtrans)
        if !isempty(en_imm)
            types[i] = VANISHING
            for tr in en_imm
                m2 = firingfunc(pn, tr)(m)
                indomain(m2) || continue
                j, isnew = intern!(m2)
                isnew && push!(stack, j)
                push!(edges, MarkEdge(i, j, tr.id, :imm, tr.weight))
            end
        else
            en_exp = filter(tr -> enabled(tr, m), pn.exptrans)
            en_gen = filter(tr -> gen_status(tr, m) == GENABLE, pn.gentrans)
            if !isempty(en_exp) || !isempty(en_gen)
                types[i] = TANGIBLE
                for tr in en_exp
                    m2 = firingfunc(pn, tr)(m)
                    indomain(m2) || continue
                    j, isnew = intern!(m2)
                    isnew && push!(stack, j)
                    push!(edges, MarkEdge(i, j, tr.id, :exp, tr.rate))
                end
                for tr in en_gen
                    m2 = firingfunc(pn, tr)(m)
                    indomain(m2) || continue
                    j, isnew = intern!(m2)
                    isnew && push!(stack, j)
                    push!(edges, MarkEdge(i, j, tr.id, :gen, 1.0))
                end
            else
                types[i] = ABSORBING
            end
        end

        length(states) > maxstates && error(
            "reachability_graph: exceeded maxstates = $maxstates distinct markings; " *
            "the net may be unbounded or place capacities (max) too large.")
    end

    MarkingGraph(pn, states, index, types, genvecs, edges, i0)
end

"""
    nstates(mg) -> Int

Number of reachable markings.
"""
nstates(mg::MarkingGraph) = length(mg.states)

"""
    tangible_states(mg) -> Vector{Int}

Indices of markings that are `TANGIBLE` or `ABSORBING` (i.e. not vanishing) —
the states over which a CTMC generator is defined.
"""
tangible_states(mg::MarkingGraph) = findall(t -> t != VANISHING, mg.types)

"""
    vanishing_states(mg) -> Vector{Int}

Indices of `VANISHING` markings.
"""
vanishing_states(mg::MarkingGraph) = findall(t -> t == VANISHING, mg.types)

"""
    has_gen(mg) -> Bool

Whether the graph contains any general-transition (`:gen`) firing, i.e. the net
is an MRSPN that must be analysed with [`mrspn`](@ref) rather than [`generator`](@ref).
"""
has_gen(mg::MarkingGraph) = any(e -> e.kind === :gen, mg.edges)

function Base.show(io::IO, mg::MarkingGraph)
    nv = count(==(VANISHING), mg.types)
    nt = count(==(TANGIBLE), mg.types)
    na = count(==(ABSORBING), mg.types)
    ng = count(e -> e.kind === :gen, mg.edges)
    print(io, "MarkingGraph($(nstates(mg)) states: $nt tangible, $nv vanishing, ",
          "$na absorbing; $(length(mg.edges)) edges, $ng gen)")
end
