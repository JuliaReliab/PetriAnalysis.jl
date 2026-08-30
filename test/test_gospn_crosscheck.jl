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

# The rows of a mark<G> matrix as marking vectors in our place order.
function gospn_markings(theirs, places, ourplaces)
    col = Dict(p => i for (i, p) in enumerate(ourplaces))
    [[Int(theirs[k, findfirst(==(p), places)]) for p in ourplaces] for k in axes(theirs, 1)]
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

# The net of test/data/spnp_example2.spn, by hand. All five of its IMM transitions
# share one priority, which is what makes it expressible here: PetriStructure has none.
function spnp_example2()
    @petrinet begin
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
end

# A block gospn wrote, or a zero matrix of the right size when it wrote none. Which
# blocks exist depends on the net *and* on the search: `mark -t` vanishes immediate
# markings as it goes, so its file has no I0I0I and gains a direct G0A0E.
blockOrZero(z, name, r, c) =
    haskey(z, name * ".shape") ? Matrix(gospn_sparse(z, name)) : zeros(r, c)

# The generator over the tangible markings, assembled from gospn's blocks by eliminating
# the vanishing ones: Q_TT + R_TV (I - P_VV)^-1 P_VT. This is the same computation
# generator() performs, written over the file -- which is the point, the two reach it by
# different routes. Rows are G0 then A0.
function gospn_tangible_generator(z)
    nG, nA = size(z["markG0"], 1), size(z["markA0"], 1)
    nI = haskey(z, "markI0") ? size(z["markI0"], 1) : 0
    QTT = [blockOrZero(z, "G0G0E", nG, nG) blockOrZero(z, "G0A0E", nG, nA)
           blockOrZero(z, "A0G0E", nA, nG) blockOrZero(z, "A0A0E", nA, nA)]
    nI == 0 && return QTT
    RTV = vcat(blockOrZero(z, "G0I0E", nG, nI), blockOrZero(z, "A0I0E", nA, nI))
    PVV = blockOrZero(z, "I0I0I", nI, nI)
    PVT = hcat(blockOrZero(z, "I0G0I", nI, nG), blockOrZero(z, "I0A0I", nI, nA))
    QTT + RTV * ((I - PVV) \ PVT)
end

# The tangible markings of a gospn file, in the row order of the assembled generator.
gospn_tangible_marks(z, places, ourplaces) =
    vcat(gospn_markings(z["markG0"], places, ourplaces),
         gospn_markings(z["markA0"], places, ourplaces))

function check_gspn_generator(name, npz, pn)
    @testset "gospn crosscheck: GSPN generator ($name)" begin
        mg = reachability_graph(pn)
        Q, st = generator(mg; tangible = true)
        z = npzread(joinpath(datadir, npz))

        G = gospn_tangible_generator(z)
        @test size(G) == size(Q)
        @test all(abs.(vec(sum(G, dims = 2))) .< 1e-12)     # still a generator

        ourplaces = [p.label for p in pn.places]
        theirs = gospn_tangible_marks(z, gospn_places(z), ourplaces)
        ours = [mg.states[i] for i in st]
        index = Dict(m => k for (k, m) in enumerate(theirs))
        @assert length(index) == length(theirs) "gospn markings are not distinct"
        perm = [index[m] for m in ours]

        @test Matrix(Q) ≈ G[perm, perm] rtol = 1e-12
    end
end

check_gspn_generator("spnp_example2", "spnp_example2.npz", spnp_example2())

# The same net through `gospn mark -t`, which runs a different search: it vanishes the
# immediate markings during the reachability walk rather than leaving them for the
# elimination. Fewer states reach the file -- 10 rather than 11, and no I0I0I block --
# and the generator over the tangible markings has to be the same one.
check_gspn_generator("spnp_example2 via mark -t", "spnp_example2_tangible.npz", spnp_example2())

# -- MRSPN --------------------------------------------------------------------

# gospn labels a group G0/I0/A0 and names a block <src><dest><kind>, where kind is "E",
# "I" or "P<k>".
#
# The block matrices alone cannot say which general transition a P<k> block is, nor what
# distribution governs a group: a gen block is a 0/1 jump matrix, so changing det(5) to
# det(99) leaves every matrix identical. gospn 0.22.0 added two text elements for exactly
# that -- `gentrans` and `groupgen` -- and they are checked here too. Before them this
# test could compare every block and still miss a wrong distribution; it was the reason
# they exist.

# gentrans is "P<k>\t<transition>\t<distribution>" per line.
gospn_gentrans(z) = Dict(p[1] => (p[2], p[3]) for p in
    (split(l, "\t") for l in split(String(z["gentrans"]), "\n")))

# groupgen is "<group>\t<transition>\t<status>\t<distribution>" per line; status is E
# for aging and P for preempted.
function gospn_groupgen(z)
    out = Dict{String,Vector{Tuple{String,String,String}}}()
    haskey(z, "groupgen") || return out
    for l in split(String(z["groupgen"]), "\n")
        g, tr, st, dist = split(l, "\t")
        push!(get!(out, String(g), Tuple{String,String,String}[]), (String(tr), String(st), String(dist)))
    end
    out
end

# The same rendering gospn uses, so the two can be compared as strings.
gospn_dist(d::DetDist) = "det(" * fmtnum(d.value) * ")"
gospn_dist(d::UnifDist) = "unif(" * fmtnum(d.a) * "," * fmtnum(d.b) * ")"
gospn_dist(d::ExpDist) = "expdist(" * fmtnum(d.rate) * ")"
fmtnum(x::Float64) = isinteger(x) ? string(Int(x)) : string(x)
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

        # Which transition each P<k> is, and with which distribution. Without the
        # gentrans element a net with two general transitions could not be checked at
        # all -- P0 is Trebuild in raid6.spn and Trecon in raid10.spn.
        gt = gospn_gentrans(z)
        byname = Dict(t.label => t for t in pn.gentrans)
        @test Set(keys(gt)) ⊆ Set("P" * string(k) for k in 0:length(pn.gentrans)-1)
        for (_, (label, dist)) in gt
            @test haskey(byname, label)
            @test dist == gospn_dist(byname[label].dist)
        end
        kinds = vcat("E", "I", collect(keys(gt)))

        # Every block gospn wrote must be one this loop actually looks at. Without
        # this the test passes when a name is wrong on both sides at once: an absent
        # block reads as zero, and zero equals zero.
        written = Set(k[1:end-6] for k in keys(z) if endswith(k, ".shape"))
        visited = Set{String}()

        for src in keys(m), dst in keys(m)
            i, j = m[src][1], m[dst][1]
            for kind in kinds
                push!(visited, src * dst * kind)
            end
            @test gospn_block(z, an, m, src, dst, "E") ≈ Matrix(exp_block(an, i, j)) rtol = 1e-12
            @test gospn_block(z, an, m, src, dst, "I") ≈ Matrix(imm_block(an, i, j)) rtol = 1e-12
            for (pk, (label, _)) in gt
                trid = byname[label].id
                @test gospn_block(z, an, m, src, dst, pk) ≈ Matrix(gen_block(an, i, j, trid)) rtol = 1e-12
            end
        end
        @test isempty(setdiff(written, visited))

        # Which general transitions govern each group. active_gens returns the aging
        # ones; gospn also reports the preempted ones, which is more than we model.
        gg = gospn_groupgen(z)
        for (label, (gid, _)) in m
            theirs = sort([tr for (tr, st, _) in get(gg, label, []) if st == "E"])
            @test theirs == sort([t.label for t in active_gens(an, gid)])
            for (tr, _, dist) in get(gg, label, [])
                @test dist == gospn_dist(byname[tr].dist)
            end
        end

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

# Two independent general transitions with different distributions: the P<k> numbering
# has something to get wrong, and the groups cover several GenVec combinations.
function two_gen()
    pn = petri()
    a1 = place(pn, "a1", 1, 1); b1 = place(pn, "b1", 0, 1)
    a2 = place(pn, "a2", 1, 1); b2 = place(pn, "b2", 0, 1)
    T1 = gentrans(pn, "T1", detdist(1.0))
    T2 = gentrans(pn, "T2", unifdist(2.0, 4.0))
    arc(pn, a1, T1); arc(pn, T1, b1)
    arc(pn, a2, T2); arc(pn, T2, b2)
    pn
end

check_mrspn("two general transitions", "mrspn_two_gen.npz", two_gen())
