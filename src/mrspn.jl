"""
    MRGroup

A regeneration class of the MRSPN marking graph: the set of markings that share a
`(StateType, GenVec)` pair. `marks` holds the member state indices (into the
`MarkingGraph`) in local order — a marking's position in `marks` is its **local
index**, the row/column index used by the block matrices. `label` is a short tag
(`"G0"` tangible, `"I0"` vanishing, `"A0"` absorbing).
"""
struct MRGroup
    id::Int
    gtype::StateType
    genvec::Vector{GenStatus}
    marks::Vector{Int}
    label::String
end

"""
    MRSPNGraph

The MRSPN (Markov regenerative) decomposition of a marking graph, mirroring the
block matrices that `gospn`'s `mark` command produces. The markings are
partitioned into regeneration groups ([`MRGroup`]); the dynamics between/within
groups are given as sparse blocks:

- `expblocks[(i, j)]` — exponential block from group `i` to group `j`. The
  self-block `(i, i)` of a tangible group is the **generator block of the
  subordinated CTMC**: off-diagonal entries are the timed rates between
  same-group markings, and each diagonal entry is the negated *total* exponential
  exit rate of that marking (so its row sum equals the rate leaving the group).
- `immblocks[(i, j)]` — immediate block from a vanishing group `i`, row-normalised
  to **branch probabilities** (each vanishing marking's immediate out-weights sum
  to 1 across all destination groups).
- `genblocks[(i, j, trid)]` — the **regeneration jump** of general transition
  `trid`: a 0/1 matrix whose every active row has a single 1, mapping a marking to
  its successor when that general transition fires.

The blocks, together with each tangible group's governing distribution(s) (see
[`active_gens`](@ref)), are the inputs to a Markov regenerative process solve,
which is left to the caller.
"""
struct MRSPNGraph
    mg::MarkingGraph
    groups::Vector{MRGroup}
    groupof::Vector{Int}    # state index -> group id
    localof::Vector{Int}    # state index -> local index within its group
    expblocks::Dict{Tuple{Int,Int},SparseMatrixCSC{Float64,Int}}
    immblocks::Dict{Tuple{Int,Int},SparseMatrixCSC{Float64,Int}}
    genblocks::Dict{Tuple{Int,Int,Int},SparseMatrixCSC{Float64,Int}}
end

const _COO = Vector{Tuple{Int,Int,Float64}}

"""
    mrspn(mg) -> MRSPNGraph
    mrspn(pn; maxstates = 10_000_000) -> MRSPNGraph

Build the MRSPN block-matrix decomposition from a marking graph (or directly from
a net). Markings are grouped by `(StateType, GenVec)` into regeneration classes,
and the exponential / immediate / general blocks are assembled (see
[`MRSPNGraph`]). Works for any net; for a pure SPN/GSPN (no general transitions)
prefer [`generator`](@ref).
"""
function mrspn(mg::MarkingGraph)
    n = nstates(mg)

    # --- 1. group states by (type, genvec) -----------------------------------
    keyof(i) = (mg.types[i], mg.genvecs[i])
    groupmap = Dict{Tuple{StateType,Vector{GenStatus}},Vector{Int}}()
    for i in 1:n
        push!(get!(groupmap, keyof(i), Int[]), i)
    end
    keys_sorted = sort!(collect(keys(groupmap)); by = k -> (Int.(k[2]), Int(k[1])))

    groups = Vector{MRGroup}()
    groupof = zeros(Int, n)
    localof = zeros(Int, n)
    counter = Dict(VANISHING => 0, TANGIBLE => 0, ABSORBING => 0)
    for (gid, k) in enumerate(keys_sorted)
        gtype, gv = k
        marks = sort(groupmap[k])
        for (li, s) in enumerate(marks)
            groupof[s] = gid
            localof[s] = li
        end
        prefix = gtype == TANGIBLE ? "G" : gtype == VANISHING ? "I" : "A"
        push!(groups, MRGroup(gid, gtype, gv, marks, string(prefix, counter[gtype])))
        counter[gtype] += 1
    end
    gsize(gid) = length(groups[gid].marks)

    # --- 2. per-marking out-sums ---------------------------------------------
    expsum = zeros(Float64, n)
    immsum = zeros(Float64, n)
    for e in mg.edges
        e.kind === :exp && (expsum[e.src] += e.value)
        e.kind === :imm && (immsum[e.src] += e.value)
    end

    # --- 3. COO accumulation per block ---------------------------------------
    expCOO = Dict{Tuple{Int,Int},_COO}()
    immCOO = Dict{Tuple{Int,Int},_COO}()
    genCOO = Dict{Tuple{Int,Int,Int},_COO}()
    for e in mg.edges
        gs = groupof[e.src]; gd = groupof[e.dst]
        li = localof[e.src]; lj = localof[e.dst]
        if e.kind === :exp
            push!(get!(expCOO, (gs, gd), _COO()), (li, lj, e.value))
        elseif e.kind === :imm
            push!(get!(immCOO, (gs, gd), _COO()), (li, lj, e.value / immsum[e.src]))
        else # :gen
            push!(get!(genCOO, (gs, gd, e.trid), _COO()), (li, lj, 1.0))
        end
    end

    # self exponential block of every tangible group must exist and carry the
    # negated total exit rate on its diagonal (the subordinated-CTMC generator).
    for g in groups
        g.gtype == TANGIBLE || continue
        v = get!(expCOO, (g.id, g.id), _COO())
        for (li, s) in enumerate(g.marks)
            push!(v, (li, li, -expsum[s]))
        end
    end

    build(coo, m, k) = sparse(getindex.(coo, 1), getindex.(coo, 2), getindex.(coo, 3), m, k, +)
    expblocks = Dict(k => build(v, gsize(k[1]), gsize(k[2])) for (k, v) in expCOO)
    immblocks = Dict(k => build(v, gsize(k[1]), gsize(k[2])) for (k, v) in immCOO)
    genblocks = Dict(k => build(v, gsize(k[1]), gsize(k[2])) for (k, v) in genCOO)

    MRSPNGraph(mg, groups, groupof, localof, expblocks, immblocks, genblocks)
end

mrspn(pn::PN; kwargs...) = mrspn(reachability_graph(pn; kwargs...))

"""
    ngroups(an) -> Int

Number of regeneration groups.
"""
ngroups(an::MRSPNGraph) = length(an.groups)

"""
    group_markings(an, gid) -> Vector{Vector{Int}}

The marking vectors belonging to group `gid`, in local-index order.
"""
group_markings(an::MRSPNGraph, gid::Int) = [an.mg.states[s] for s in an.groups[gid].marks]

"""
    active_gens(an, gid) -> Vector

The general transitions that are `GENABLE` (aging) in group `gid` — the
transitions whose firing-time distributions govern the group's sojourn. Each
returned object is a `GenTrans` carrying `.label`, `.dist` and `.policy`.
"""
function active_gens(an::MRSPNGraph, gid::Int)
    gv = an.groups[gid].genvec
    [an.mg.pn.gentrans[k] for k in eachindex(gv) if gv[k] == GENABLE]
end

"""
    exp_block(an, i, j) / imm_block(an, i, j) / gen_block(an, i, j, trid)

Return the requested block, or an all-zero sparse matrix of the correct size when
that block is absent.
"""
exp_block(an::MRSPNGraph, i::Int, j::Int) =
    get(an.expblocks, (i, j), spzeros(length(an.groups[i].marks), length(an.groups[j].marks)))
imm_block(an::MRSPNGraph, i::Int, j::Int) =
    get(an.immblocks, (i, j), spzeros(length(an.groups[i].marks), length(an.groups[j].marks)))
gen_block(an::MRSPNGraph, i::Int, j::Int, trid::Int) =
    get(an.genblocks, (i, j, trid), spzeros(length(an.groups[i].marks), length(an.groups[j].marks)))

"""
    initial_vectors(an) -> Dict{Int,Vector{Float64}}

Per-group initial probability vectors: the group containing the net's initial
marking has a 1 at that marking's local index; all other entries are 0.
"""
function initial_vectors(an::MRSPNGraph)
    s0 = an.mg.initial_state
    g0 = an.groupof[s0]
    Dict(g.id => (g.id == g0 ? begin
                      v = zeros(Float64, length(g.marks)); v[an.localof[s0]] = 1.0; v
                  end : zeros(Float64, length(g.marks))) for g in an.groups)
end

function Base.show(io::IO, g::MRGroup)
    print(io, "MRGroup(", g.label, ", ", g.gtype, ", ", length(g.marks), " marks)")
end

function Base.show(io::IO, an::MRSPNGraph)
    nt = count(g -> g.gtype == TANGIBLE, an.groups)
    nv = count(g -> g.gtype == VANISHING, an.groups)
    na = count(g -> g.gtype == ABSORBING, an.groups)
    print(io, "MRSPNGraph(", ngroups(an), " groups: $nt tangible, $nv vanishing, $na absorbing; ",
          length(an.expblocks), " exp / ", length(an.immblocks), " imm / ",
          length(an.genblocks), " gen blocks)")
end
