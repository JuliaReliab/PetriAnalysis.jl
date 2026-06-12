# Reuse the net builders defined in test_reachability.jl (included earlier).

# Solve the stationary distribution of a CTMC generator Q (rows sum to 0).
function stationary(Q)
    n = size(Q, 1)
    A = Matrix(transpose(Q))
    A[end, :] .= 1.0
    b = zeros(n); b[end] = 1.0
    A \ b
end

@testset "generator: M/M/1/K matches analytic steady state" begin
    K = 5; lam = 0.7; mu = 1.3
    Q, st = generator(mm1k(K; lam = lam, mu = mu))
    @test size(Q) == (K + 1, K + 1)
    @test all(abs.(vec(sum(Q, dims = 2))) .< 1e-12)     # generator rows sum to 0
    @test st == collect(1:K+1)

    # rows/cols are ordered by discovery; reorder to marking order 0..K
    mm1 = reachability_graph(mm1k(K; lam = lam, mu = mu))
    order = sortperm(only.(mm1.states))                  # marking value -> position
    p = stationary(Q)[order]

    rho = lam / mu
    analytic = [rho^i for i in 0:K]; analytic ./= sum(analytic)
    @test maximum(abs.(p .- analytic)) < 1e-10
end

@testset "generator: GSPN vanishing elimination (probabilistic branch)" begin
    a, b, c, w1, w2 = 2.0, 1.0, 3.0, 1.0, 3.0
    mg = reachability_graph(gspn_branch(; a = a, b = b, c = c, w1 = w1, w2 = w2))
    Q, st = generator(mg; tangible = true)

    @test size(Q) == (3, 3)                              # vanishing state removed
    @test st == tangible_states(mg)
    @test all(abs.(vec(sum(Q, dims = 2))) .< 1e-12)

    # locate rows by their markings
    row(mk) = findfirst(i -> mg.states[i] == mk, st)
    r0, r2, r3 = row([1,0,0,0]), row([0,0,1,0]), row([0,0,0,1])
    @test Q[r0, r2] ≈ a * w1 / (w1 + w2)                 # 0.5
    @test Q[r0, r3] ≈ a * w2 / (w1 + w2)                 # 1.5
    @test Q[r2, r0] ≈ b
    @test Q[r3, r0] ≈ c
end

@testset "generator: tangible=false errors when vanishing states exist" begin
    mg = reachability_graph(gspn_branch())
    @test_throws ErrorException generator(mg; tangible = false)
end

@testset "generator: timeless trap is detected" begin
    # Two immediate transitions forming a cycle with no timed exit:
    # p0 -exp T-> p1; p1 -imm I1-> p2; p2 -imm I2-> p1  (trap among {p1,p2}).
    pn = petri()
    place(pn, "p0", 1, 1); place(pn, "p1", 0, 1); place(pn, "p2", 0, 1)
    T = exptrans(pn, "T", 1.0)
    I1 = immtrans(pn, "I1", 1.0); I2 = immtrans(pn, "I2", 1.0)
    p0, p1, p2 = pn.places
    arc(pn, p0, T); arc(pn, T, p1)
    arc(pn, p1, I1); arc(pn, I1, p2)
    arc(pn, p2, I2); arc(pn, I2, p1)
    mg = reachability_graph(pn)
    @test !isempty(vanishing_states(mg))
    @test_throws ErrorException generator(mg; tangible = true)
end

@testset "generator: building blocks" begin
    mg = reachability_graph(mm1k(3))
    R = exp_rate_matrix(mg)
    @test size(R) == (4, 4)
    @test all(diag(R) .== 0.0)                           # off-diagonal only

    P = imm_prob_matrix(mg)                              # no IMM here -> empty
    @test nnz(P) == 0

    Pg = imm_prob_matrix(reachability_graph(gspn_branch()))
    vi = only(vanishing_states(reachability_graph(gspn_branch())))
    @test sum(Pg[vi, :]) ≈ 1.0                           # branch row is stochastic
end
