module PetriAnalysis

using PetriStructure
using PetriStructure: PN
using SparseArrays
using LinearAlgebra

export MarkingGraph, MarkEdge, StateType, VANISHING, TANGIBLE, ABSORBING,
       GenStatus, GDISABLE, GENABLE, GPREEMPT,
       reachability_graph, nstates, tangible_states, vanishing_states,
       has_gen, gen_status,
       generator, exp_rate_matrix, imm_prob_matrix,
       MRSPNGraph, MRGroup, mrspn, ngroups, group_markings, active_gens,
       exp_block, imm_block, gen_block, initial_vectors,
       markgraph_todot

include("reachability.jl")
include("generator.jl")
include("mrspn.jl")
include("dot.jl")

end
