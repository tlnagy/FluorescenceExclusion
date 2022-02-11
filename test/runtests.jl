using FileIO
using FluorescenceExclusion
using Test
using Documenter

@testset "Doctest" begin
    doctest(FluorescenceExclusion; manual=false)
end

include("flatfield.jl")
include("segment.jl")