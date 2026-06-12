"""
M/M/1/K queue (pure SPN)

A single buffer place of capacity K, an always-enabled arrival transition
(the place's `max` bounds it at K) and a service transition. Builds the
reachability graph and the CTMC generator, then compares the stationary
distribution against the analytic M/M/1/K result.
"""

using PetriStructure
using PetriAnalysis
using LinearAlgebra

K = 6
lam = 0.6   # arrival rate
mu = 1.0    # service rate

pn = petri()
buf = place(pn, "buf", 0, K)
arr = exptrans(pn, "arrival", lam)
srv = exptrans(pn, "service", mu)
arc(pn, arr, buf)   # arrival produces a token
arc(pn, buf, srv)   # service consumes a token

mg = reachability_graph(pn)
println(mg)

Q, states = generator(mg)
println("\nGenerator Q (", size(Q, 1), "×", size(Q, 2), "):")
display(Matrix(Q)); println()

# Stationary distribution: solve πQ = 0 with Σπ = 1.
n = size(Q, 1)
A = Matrix(transpose(Q)); A[end, :] .= 1.0
b = zeros(n); b[end] = 1.0
prob = A \ b

order = sortperm([only(mg.states[s]) for s in states])   # reorder to 0..K
println("\nStationary distribution (marking 0..K):")
println(prob[order])

rho = lam / mu
analytic = [rho^i for i in 0:K]; analytic ./= sum(analytic)
println("Analytic M/M/1/K:")
println(analytic)
println("\nMax abs error: ", maximum(abs.(prob[order] .- analytic)))
