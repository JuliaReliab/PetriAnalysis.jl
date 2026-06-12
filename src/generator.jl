"""
    exp_rate_matrix(mg) -> SparseMatrixCSC{Float64}

Off-diagonal rate matrix `R` over all `nstates(mg)` markings, where `R[i, j]` is
the sum of the rates of the exponential (`:exp`) edges from state `i` to state
`j`. Diagonal entries are zero.
"""
function exp_rate_matrix(mg::MarkingGraph)
    n = nstates(mg)
    I = Int[]; J = Int[]; V = Float64[]
    for e in mg.edges
        e.kind === :exp || continue
        push!(I, e.src); push!(J, e.dst); push!(V, e.value)
    end
    sparse(I, J, V, n, n, +)
end

"""
    imm_prob_matrix(mg) -> SparseMatrixCSC{Float64}

Row-stochastic branching matrix `P` over all `nstates(mg)` markings. For each
vanishing state, the immediate (`:imm`) out-edges are normalised by total weight
so that `P[i, j]` is the probability of moving from `i` to `j`. Rows for
non-vanishing states are empty.
"""
function imm_prob_matrix(mg::MarkingGraph)
    n = nstates(mg)
    rowsum = zeros(Float64, n)
    for e in mg.edges
        e.kind === :imm || continue
        rowsum[e.src] += e.value
    end
    I = Int[]; J = Int[]; V = Float64[]
    for e in mg.edges
        e.kind === :imm || continue
        push!(I, e.src); push!(J, e.dst); push!(V, e.value / rowsum[e.src])
    end
    sparse(I, J, V, n, n, +)
end

"""
    generator(mg; tangible = true) -> (Q, states)
    generator(pn; tangible = true) -> (Q, states)

Build the CTMC infinitesimal generator `Q` (a `SparseMatrixCSC{Float64}` whose
rows sum to zero) from a marking graph (or directly from a net).

`states` is the vector of original marking-graph indices corresponding to the
rows/columns of `Q`, so `mg.states[states[k]]` is the marking of row `k`.

- If the graph has no vanishing markings (a pure SPN), `Q` covers every state.
- Otherwise, with `tangible = true` (default), vanishing (immediate) markings are
  eliminated by the standard GSPN reduction and `Q` covers only the tangible /
  absorbing states. With `tangible = false` an error is raised, since a CTMC
  generator is only defined once immediate markings are removed (the full graph
  remains available on `mg`).

The reduction forms, from the exponential rate matrix `R` and the immediate
branching matrix `P`, the block decomposition over tangible (`T`) and vanishing
(`V`) states and computes the effective tangible-to-tangible rates

    R_eff = R_TT + R_TV * (I - P_VV)^{-1} * P_VT

An immediate cycle with no exponential exit ("timeless trap") makes `(I - P_VV)`
singular and raises an error.
"""
function generator(mg::MarkingGraph; tangible::Bool = true)
    has_gen(mg) && error(
        "generator: the net contains general (GEN) transitions, so it is an MRSPN " *
        "whose process is not a CTMC. Use `mrspn(mg)` to obtain the regenerative " *
        "block matrices instead.")
    R = exp_rate_matrix(mg)
    vstates = vanishing_states(mg)

    if isempty(vstates)
        states = collect(1:nstates(mg))
        return _generator_from_rates(R), states
    end

    tangible || error(
        "generator: the marking graph contains $(length(vstates)) vanishing " *
        "(immediate) markings; a CTMC generator requires tangible=true to " *
        "eliminate them. The full marking graph is available on the MarkingGraph.")

    T = tangible_states(mg)
    V = vstates
    P = imm_prob_matrix(mg)

    R_TT = R[T, T]
    R_TV = R[T, V]
    P_VV = P[V, V]
    P_VT = P[V, T]

    # Solve (I - P_VV) X = P_VT for X (probability of eventually exiting V into
    # each tangible state). Dense solve: the vanishing subgraph is typically small.
    A = Matrix(I - P_VV)
    local X
    try
        X = A \ Matrix(P_VT)
    catch err
        error("generator: failed to eliminate vanishing markings — the net " *
              "likely contains a timeless trap (an immediate cycle with no " *
              "timed exit). ($(err))")
    end

    R_eff = R_TT + R_TV * sparse(X)
    return _generator_from_rates(R_eff), T
end

generator(pn::PN; kwargs...) = generator(reachability_graph(pn); kwargs...)

# Turn an off-diagonal rate matrix into a generator by setting the diagonal to
# the negated row sums.
function _generator_from_rates(R::AbstractMatrix)
    n = size(R, 1)
    Q = SparseMatrixCSC{Float64,Int}(R)
    Q[diagind(Q)] .= 0.0          # ensure no stray diagonal contribution
    d = vec(sum(Q, dims = 2))
    Q + sparse(1:n, 1:n, -d, n, n)
end
