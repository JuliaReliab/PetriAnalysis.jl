# Net builders -----------------------------------------------------------------

# Failure (EXP) / deterministic repair (GEN, DET): the textbook MRSPN.
function fail_repair(; lam = 0.1, tau = 5.0)
    pn = petri()
    up = place(pn, "up", 1, 1)
    down = place(pn, "down", 0, 1)
    Tf = exptrans(pn, "Tfail", lam)
    Tr = gentrans(pn, "Trepair", detdist(tau))
    arc(pn, up, Tf); arc(pn, Tf, down)
    arc(pn, down, Tr); arc(pn, Tr, up)
    pn
end

# A general transition stays enabled across an EXP firing, so its regeneration
# group spans two markings (subordinated CTMC with a real off-diagonal rate).
function subordinated_pair(; r = 2.0, tau = 1.0)
    pn = petri()
    g = place(pn, "g", 1, 1); a = place(pn, "a", 1, 1); b = place(pn, "b", 0, 1)
    Te = exptrans(pn, "Te", r)
    Tg = gentrans(pn, "Tg", detdist(tau))
    arc(pn, a, Te); arc(pn, Te, b)      # EXP moves a -> b, leaving g untouched
    arc(pn, g, Tg)                       # GEN consumes g (ends the group)
    pn
end

# GEN -> vanishing (IMM) -> EXP -> back: exercises immediate blocks.
function gen_imm_exp(; r = 1.0, w = 1.0, tau = 1.0)
    pn = petri()
    a = place(pn, "a", 1, 1); b = place(pn, "b", 0, 1); c = place(pn, "c", 0, 1)
    Tg = gentrans(pn, "Tg", detdist(tau))
    Ti = immtrans(pn, "Ti", w)
    Te = exptrans(pn, "Te", r)
    arc(pn, a, Tg); arc(pn, Tg, b)
    arc(pn, b, Ti); arc(pn, Ti, c)
    arc(pn, c, Te); arc(pn, Te, a)
    pn
end

# Invariant helpers (independent of group ordering) ----------------------------

group_exp_rowsum(an, g) =
    sum((vec(sum(exp_block(an, g, j), dims = 2)) for j in 1:ngroups(an)); init = zeros(length(an.groups[g].marks)))
group_imm_rowsum(an, g) =
    sum((vec(sum(imm_block(an, g, j), dims = 2)) for j in 1:ngroups(an)); init = zeros(length(an.groups[g].marks)))
group_gen_rowsum(an, g, trid) =
    sum((vec(sum(gen_block(an, g, j, trid), dims = 2)) for j in 1:ngroups(an)); init = zeros(length(an.groups[g].marks)))

# -----------------------------------------------------------------------------

@testset "MRSPN: reachability is GEN-aware" begin
    mg = reachability_graph(fail_repair())
    @test nstates(mg) == 2
    @test has_gen(mg)
    @test all(t -> t == TANGIBLE, mg.types)
    # the down marking has Trepair aging; the up marking does not
    up_i = findfirst(s -> s == [1, 0], mg.states)
    dn_i = findfirst(s -> s == [0, 1], mg.states)
    @test mg.genvecs[up_i] == [GDISABLE]
    @test mg.genvecs[dn_i] == [GENABLE]
    @test count(e -> e.kind === :gen, mg.edges) == 1
    @test count(e -> e.kind === :exp, mg.edges) == 1
end

@testset "MRSPN: fail/repair block matrices" begin
    lam, tau = 0.1, 5.0
    an = mrspn(fail_repair(; lam = lam, tau = tau))
    @test ngroups(an) == 2

    # the group whose GenVec has Trepair enabled
    gid = findfirst(g -> GENABLE in g.genvec, an.groups)
    ags = active_gens(an, gid)
    @test length(ags) == 1
    @test ags[1].label == "Trepair"
    @test ags[1].dist == DetDist(tau)

    # subordinated generator: each tangible group's exp rows sum to 0
    for g in 1:ngroups(an)
        an.groups[g].gtype == TANGIBLE || continue
        @test all(abs.(group_exp_rowsum(an, g)) .< 1e-12)
    end
    # the "up" group exits at rate lam into the "down" group
    up_gid = findfirst(g -> g.genvec == [GDISABLE], an.groups)
    @test exp_block(an, up_gid, up_gid)[1, 1] ≈ -lam
    # GEN jump is row-stochastic out of the aging group
    @test all(group_gen_rowsum(an, gid, ags[1].id) .≈ 1.0)
end

@testset "MRSPN: subordinated group spans two markings" begin
    r = 2.0
    an = mrspn(subordinated_pair(; r = r))
    g2 = findfirst(g -> length(g.marks) == 2, an.groups)
    @test g2 !== nothing
    @test an.groups[g2].gtype == TANGIBLE
    @test length(active_gens(an, g2)) == 1

    # within-group EXP: [-r r; 0 0] (m0 ->(r) m1; m1 has no EXP)
    B = exp_block(an, g2, g2)
    @test Matrix(B) ≈ [-r r; 0.0 0.0]
    @test all(abs.(group_exp_rowsum(an, g2)) .< 1e-12)

    # both markings enable the GEN, so every gen row is stochastic
    tg = active_gens(an, g2)[1]
    @test all(group_gen_rowsum(an, g2, tg.id) .≈ 1.0)
end

@testset "MRSPN: immediate blocks are branch probabilities" begin
    an = mrspn(gen_imm_exp())
    vgid = findfirst(g -> g.gtype == VANISHING, an.groups)
    @test vgid !== nothing                       # the IMM marking forms a vanishing group
    @test all(group_imm_rowsum(an, vgid) .≈ 1.0) # immediate out-weights normalise to 1
    @test any(g -> g.gtype == TANGIBLE, an.groups)
end

@testset "MRSPN: generator refuses a net with GEN" begin
    mg = reachability_graph(fail_repair())
    @test_throws ErrorException generator(mg)
end

@testset "MRSPN: pure SPN/GSPN still works via generator" begin
    # no GEN -> has_gen false, generator unaffected
    pn = petri()
    buf = place(pn, "buf", 0, 3)
    arr = exptrans(pn, "arrival", 0.5)
    srv = exptrans(pn, "service", 1.0)
    arc(pn, arr, buf); arc(pn, buf, srv)
    mg = reachability_graph(pn)
    @test !has_gen(mg)
    Q, _ = generator(mg)
    @test all(abs.(vec(sum(Q, dims = 2))) .< 1e-12)
end
