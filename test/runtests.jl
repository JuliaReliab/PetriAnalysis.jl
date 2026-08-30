using PetriStructure
using PetriAnalysis
using SparseArrays
using LinearAlgebra
using Test

@testset verbose=true "PetriAnalysis.jl" begin
    include("test_reachability.jl")
    include("test_generator.jl")
    include("test_mrspn.jl")
    include("test_gospn_crosscheck.jl")
end
