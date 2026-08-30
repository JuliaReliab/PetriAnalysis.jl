# gospn and PetriAnalysis.jl implement the same semantics twice, in two languages, with
# nothing comparing them. These tests read a result file gospn wrote and check that the
# generator it implies is the one PetriAnalysis builds.
#
# The two enumerate the state space in different orders -- gospn in its search order,
# this package in DFS discovery order -- so nothing can be compared position by
# position. Every comparison here is keyed on the marking vector, using the `place` and
# `mark<G>` elements gospn has written since 0.20.0.
#
# The fixtures in test/data were produced with:
#   gospn mark -i spnp_example1.spn -o spnp_example1.npz

using NPZ

datadir = joinpath(@__DIR__, "data")

# gospn stores a sparse matrix as four CSC members with a 0-origin.
function gospn_sparse(z, name)
    rows, cols = z[name*".shape"]
    SparseMatrixCSC(Int(rows), Int(cols),
        Int.(z[name*".indptr"]) .+ 1,
        Int.(z[name*".indices"]) .+ 1,
        z[name*".data"])
end

gospn_places(z) = split(String(z["place"]), "\n")

# rowperm maps a row of gospn's matrix to a row of ours, by marking. `theirs` is the
# mark<G> matrix, `ours` the marking of each of our rows; `places` is gospn's column
# order, which need not be ours.
function rowperm(theirs, places, ours, ourplaces)
    col = Dict(p => i for (i, p) in enumerate(ourplaces))
    index = Dict{Vector{Int},Int}()
    for k in axes(theirs, 1)
        m = zeros(Int, length(ourplaces))
        for (i, p) in enumerate(places)
            m[col[p]] = Int(theirs[k, i])
        end
        index[m] = k
    end
    @assert length(index) == size(theirs, 1) "gospn markings are not distinct"
    [index[m] for m in ours]
end

@testset "gospn crosscheck: SPN generator (spnp_example1)" begin
    # The same net as test/data/spnp_example1.spn, by hand.
    pn = @petrinet begin
        p0[3]; p1[0]; p2[0]; p3[0]; p4[0]
        exp(1.0): t0
        exp(3.0): t1
        exp(7.0): t2
        exp(9.0): t3
        exp(5.0): t4
        p0 => t0; t0 => p1; t0 => p2
        p1 => t1; t1 => p3
        p2 => t2; t2 => p4
        p3 => t3; t3 => p1
        p3 => t4; p4 => t4; t4 => p0
    end

    mg = reachability_graph(pn)
    Q, st = generator(mg; tangible = true)

    z = npzread(joinpath(datadir, "spnp_example1.npz"))
    G = gospn_sparse(z, "G0G0E")
    @test size(G) == size(Q)

    # A pure-EXP net has one group, and its matrix is the generator, diagonal included.
    perm = rowperm(z["markG0"], gospn_places(z),
        [mg.states[i] for i in st], [p.label for p in pn.places])
    @test Matrix(Q) ≈ Matrix(G[perm, perm]) rtol = 1e-12

    # The initial state is the same one, not merely a state of the same degree.
    @test findfirst(!iszero, z["initG0"][perm]) == findfirst(==(mg.initial_state), st)
end

@testset "gospn crosscheck: GSPN generator (spnp_example2)" begin
    # test/data/spnp_example2.spn. All five IMM transitions share one priority there,
    # which is what makes this net expressible here: PetriStructure has no priorities.
    pn = @petrinet begin
        p0[1]; p1[0]; p2[0]; p3[0]; p4[0]; p5[0]; p6[0]; p7[0]; p8[0]
        exp(1.0): A
        exp(0.3): B1
        exp(0.2): C
        exp(7.0): D
        imm(0.4): t2
        imm(0.6): t3
        imm(0.05): t6
        imm(0.95): t7
        imm(1.0): t8
        p0 => A; A => p1; A => p3
        p1 => B1; B1 => p2
        p3 => t2; t2 => p4
        p3 => t3; t3 => p5
        p4 => C; C => p6
        p5 => D; D => p7
        p7 => t6; t6 => p6
        p7 => t7; t7 => p5
        p2 => t8; p6 => t8; t8 => p8
    end

    mg = reachability_graph(pn)
    Q, st = generator(mg; tangible = true)

    z = npzread(joinpath(datadir, "spnp_example2.npz"))

    # gospn keeps the groups apart: tangible (G0), vanishing (I0) and absorbing (A0).
    # Eliminating the vanishing states is the same computation generator() performs,
    # written over gospn's blocks -- which is the point: the two arrive here by
    # different routes.
    nG, nI, nA = size(z["markG0"], 1), size(z["markI0"], 1), size(z["markA0"], 1)
    QTT = blockdiag(gospn_sparse(z, "G0G0E"), spzeros(nA, nA))
    RTV = vcat(gospn_sparse(z, "G0I0E"), spzeros(nA, nI))
    PVV = gospn_sparse(z, "I0I0I")
    PVT = hcat(gospn_sparse(z, "I0G0I"), gospn_sparse(z, "I0A0I"))
    G = QTT + RTV * ((I - PVV) \ Matrix(PVT))

    @test size(G) == size(Q)
    @test all(abs.(vec(sum(G, dims = 2))) .< 1e-12)     # still a generator

    perm = rowperm(vcat(z["markG0"], z["markA0"]), gospn_places(z),
        [mg.states[i] for i in st], [p.label for p in pn.places])
    @test Matrix(Q) ≈ G[perm, perm] rtol = 1e-12
end

# -- MRSPN --------------------------------------------------------------------

# gospn labels a group G0/I0/A0 and names a block <src><dest><kind>, where kind is "E",
# "I" or "P<k>".
#
# Two things the block matrices alone cannot say, and which gospn's result file does not
# record either:
#
#   - which general transition P<k> is. The k is a counter over the general transitions
#     in the order their blocks were found, so the fixtures here have exactly one, where
#     P0 is unambiguous.
#   - the distribution governing each group's sojourn. A gen block is a 0/1 jump matrix,
#     so changing Trepair from det(5) to det(99) leaves every matrix in the file
#     identical -- this test was checked, and it does not fail on that.
#
# So this compares the *structure* of the regenerative process, which is what the file
# holds. An MRP solved from the file needs the distributions from somewhere else.
gospn_groups(z) = sort([k[5:end] for k in keys(z) if startswith(k, "mark")])

# Match gospn's groups to ours by the set of markings they hold: neither the group order
# nor the marking order agrees between the two, and nothing else identifies a group.
# Returns gospn label -> (our group id, our-row -> gospn-row permutation).
function match_groups(z, an, pn)
    ourplaces = [p.label for p in pn.places]
    places = gospn_places(z)
    out = Dict{String,Tuple{Int,Vector{Int}}}()
    for label in gospn_groups(z)
        theirs = z["mark"*label]
        gid = findfirst(g -> Set(group_markings(an, g)) ==
                             Set(gospn_markings(theirs, places, ourplaces)),
                        1:ngroups(an))
        @assert gid !== nothing "no group of ours holds the markings of gospn's $label"
        out[label] = (gid, rowperm(theirs, places, group_markings(an, gid), ourplaces))
    end
    out
end

# The rows of a mark<G> matrix as marking vectors in our place order.
function gospn_markings(theirs, places, ourplaces)
    col = Dict(p => i for (i, p) in enumerate(ourplaces))
    [[Int(theirs[k, findfirst(==(p), places)]) for p in ourplaces] for k in axes(theirs, 1)]
end

# The block gospn wrote for this pair, or a zero matrix when it wrote none -- an absent
# block and an all-zero block mean the same thing, and the two sides disagree about
# which of them to store.
function gospn_block(z, an, m, src, dst, kind)
    (i, pi), (j, pj) = m[src], m[dst]
    name = src * dst * kind
    if !haskey(z, name * ".shape")
        return zeros(length(an.groups[i].marks), length(an.groups[j].marks))
    end
    Matrix(gospn_sparse(z, name))[pi, pj]
end

function check_mrspn(name, npz, pn)
    @testset "gospn crosscheck: MRSPN blocks ($name)" begin
        z = npzread(joinpath(datadir, npz))
        an = mrspn(pn)
        m = match_groups(z, an, pn)
        @test length(m) == ngroups(an)
        @test length(pn.gentrans) == 1          # P0 is only unambiguous with one
        trid = pn.gentrans[1].id

        # Every block gospn wrote must be one this loop actually looks at. Without
        # this the test passes when a name is wrong on both sides at once: an absent
        # block reads as zero, and zero equals zero.
        written = Set(k[1:end-6] for k in keys(z) if endswith(k, ".shape"))
        visited = Set{String}()

        for src in keys(m), dst in keys(m)
            i, j = m[src][1], m[dst][1]
            for kind in ("E", "I", "P0")
                push!(visited, src * dst * kind)
            end
            @test gospn_block(z, an, m, src, dst, "E") ≈ Matrix(exp_block(an, i, j)) rtol = 1e-12
            @test gospn_block(z, an, m, src, dst, "I") ≈ Matrix(imm_block(an, i, j)) rtol = 1e-12
            @test gospn_block(z, an, m, src, dst, "P0") ≈ Matrix(gen_block(an, i, j, trid)) rtol = 1e-12
        end
        @test isempty(setdiff(written, visited))

        # And the initial marking is the same one, not merely a marking of some group.
        ours = initial_vectors(an)
        for (label, (gid, perm)) in m
            @test z["init"*label][perm] ≈ ours[gid]
        end
    end
end

check_mrspn("fail/repair", "mrspn_fail_repair.npz", fail_repair())
check_mrspn("subordinated group", "mrspn_subordinated.npz", subordinated_pair())
check_mrspn("gen -> imm -> exp", "mrspn_gen_imm_exp.npz", gen_imm_exp())
