"""
GSPN with immediate transitions (vanishing-state elimination)

A timed transition moves a token into a marking where two *immediate*
transitions compete. That marking is vanishing — it is left in zero time — so the
CTMC generator is defined only over the tangible markings. `generator(mg)`
eliminates the vanishing marking, folding the immediate branch probabilities
(weight w1 vs w2) into the effective timed rates.
"""

using PetriStructure
using PetriAnalysis

a, b, c = 2.0, 1.0, 3.0    # timed rates
w1, w2 = 1.0, 3.0          # immediate weights (branch 1:3)

pn = petri()
place(pn, "p0", 1, 1); place(pn, "p1", 0, 1)
place(pn, "p2", 0, 1); place(pn, "p3", 0, 1)
T1 = exptrans(pn, "T1", a); T2 = exptrans(pn, "T2", b); T3 = exptrans(pn, "T3", c)
I1 = immtrans(pn, "I1", w1); I2 = immtrans(pn, "I2", w2)
p0, p1, p2, p3 = pn.places
arc(pn, p0, T1); arc(pn, T1, p1)     # p0 -> p1 (vanishing)
arc(pn, p1, I1); arc(pn, I1, p2)     # immediate split ...
arc(pn, p1, I2); arc(pn, I2, p3)
arc(pn, p2, T2); arc(pn, T2, p0)
arc(pn, p3, T3); arc(pn, T3, p0)

mg = reachability_graph(pn)
println(mg)
println("\nMarkings:")
for i in 1:nstates(mg)
    println("  s$i = ", mg.states[i], "  ", mg.types[i])
end

Q, states = generator(mg; tangible = true)
println("\nTangible CTMC generator (vanishing marking eliminated):")
println("rows correspond to markings: ", [mg.states[s] for s in states])
display(Matrix(Q)); println()

println("\nEffective rate p0 -> p2 = a·w1/(w1+w2) = ", a * w1 / (w1 + w2))
println("Effective rate p0 -> p3 = a·w2/(w1+w2) = ", a * w2 / (w1 + w2))

println("\nDOT of the full marking graph:\n")
println(markgraph_todot(mg))
