# Helpers ---------------------------------------------------------------------

# M/M/1/K queue: single place of capacity K, an always-enabled arrival
# transition (capacity bound stops it at K) and a service transition.
function mm1k(K; lam = 0.6, mu = 1.0)
    pn = petri()
    buf = place(pn, "buf", 0, K)
    arr = exptrans(pn, "arrival", lam)
    srv = exptrans(pn, "service", mu)
    arc(pn, arr, buf)
    arc(pn, buf, srv)
    pn
end

# GSPN: p0 -exp T1-> p1 (vanishing); p1 -imm I1(w1)-> p2, -imm I2(w2)-> p3;
# p2 -exp T2-> p0; p3 -exp T3-> p0.
function gspn_branch(; a = 2.0, b = 1.0, c = 3.0, w1 = 1.0, w2 = 3.0)
    pn = petri()
    place(pn, "p0", 1, 1); place(pn, "p1", 0, 1)
    place(pn, "p2", 0, 1); place(pn, "p3", 0, 1)
    T1 = exptrans(pn, "T1", a); T2 = exptrans(pn, "T2", b); T3 = exptrans(pn, "T3", c)
    I1 = immtrans(pn, "I1", w1); I2 = immtrans(pn, "I2", w2)
    p0, p1, p2, p3 = pn.places
    arc(pn, p0, T1); arc(pn, T1, p1)
    arc(pn, p1, I1); arc(pn, I1, p2)
    arc(pn, p1, I2); arc(pn, I2, p3)
    arc(pn, p2, T2); arc(pn, T2, p0)
    arc(pn, p3, T3); arc(pn, T3, p0)
    pn
end

# -----------------------------------------------------------------------------

@testset "reachability: M/M/1/K enumeration" begin
    K = 4
    mg = reachability_graph(mm1k(K))
    @test nstates(mg) == K + 1
    @test sort(mg.states) == [[i] for i in 0:K]
    @test all(t -> t == TANGIBLE, mg.types)        # every state has a timed move
    @test isempty(vanishing_states(mg))
    # birth/death edges only, with the right rates
    @test all(e -> e.kind === :exp, mg.edges)
end

@testset "reachability: capacity bound keeps graph finite" begin
    # arrival with no service: without the max bound this would be unbounded.
    pn = petri()
    buf = place(pn, "buf", 0, 3)
    arr = exptrans(pn, "arrival", 1.0)
    arc(pn, arr, buf)
    mg = reachability_graph(pn)
    @test nstates(mg) == 4                          # markings 0..3 only
    @test maximum(only.(mg.states)) == 3            # never exceeds capacity
end

@testset "reachability: GSPN classification & IMM suppresses EXP" begin
    mg = reachability_graph(gspn_branch())
    @test nstates(mg) == 4
    @test count(==(VANISHING), mg.types) == 1
    @test count(==(TANGIBLE), mg.types) == 3
    # the vanishing marking [0,1,0,0] emits only immediate edges
    vi = only(vanishing_states(mg))
    @test mg.states[vi] == [0, 1, 0, 0]
    out = filter(e -> e.src == vi, mg.edges)
    @test !isempty(out)
    @test all(e -> e.kind === :imm, out)
end

@testset "reachability: absorbing state" begin
    # a single timed move into a dead marking
    pn = petri()
    a = place(pn, "a", 1, 1)
    b = place(pn, "b", 0, 1)
    t = exptrans(pn, "t", 1.0)
    arc(pn, a, t); arc(pn, t, b)
    mg = reachability_graph(pn)
    @test nstates(mg) == 2
    @test count(==(ABSORBING), mg.types) == 1
    @test mg.states[only(findall(==(ABSORBING), mg.types))] == [0, 1]
end
