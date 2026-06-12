"""
MRSPN: failure / deterministic repair (block matrices)

A component fails at exponential rate `lam` and is repaired in a *deterministic*
time `tau` (a general transition). The process is a Markov regenerative process,
not a CTMC, so `generator` does not apply. `mrspn` groups the markings into
regeneration classes (by GenVec) and emits the block matrices — the
subordinated-CTMC generator blocks, immediate branch-probability blocks, and
general-transition jump blocks — that feed an MRGP solve.
"""

using PetriStructure
using PetriAnalysis

lam, tau = 0.1, 5.0

pn = petri()
up = place(pn, "up", 1, 1)
down = place(pn, "down", 0, 1)
Tf = exptrans(pn, "Tfail", lam)
Tr = gentrans(pn, "Trepair", detdist(tau))   # deterministic repair time
arc(pn, up, Tf); arc(pn, Tf, down)
arc(pn, down, Tr); arc(pn, Tr, up)

mg = reachability_graph(pn)
println(mg)

an = mrspn(mg)
println(an)
println()

for g in an.groups
    ags = active_gens(an, g.id)
    gens = isempty(ags) ? "—" : join(["$(t.label)~$(t.dist)" for t in ags], ", ")
    println("group $(g.label) [$(g.gtype)]: markings ", group_markings(an, g.id),
            "  aging GEN: ", gens)
end

println("\nEXP blocks (subordinated-CTMC generators; rows sum to the group exit rate):")
for k in sort(collect(keys(an.expblocks)))
    println("  groups $k => ", Matrix(an.expblocks[k]))
end

println("\nGEN jump blocks (regeneration; row-stochastic 0/1), keyed (src, dst, trid):")
for k in sort(collect(keys(an.genblocks)))
    trlabel = pn.trans[k[3]].label
    println("  $k  [$trlabel] => ", Matrix(an.genblocks[k]))
end

println("\nInitial per-group vectors: ", initial_vectors(an))
