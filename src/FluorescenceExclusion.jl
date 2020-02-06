module FluorescenceExclusion

using Images
using ImageFiltering
using ScatteredInterpolation
using ImageBinarization
using Statistics
using Distributions
using ProgressMeter

include("flatfield.jl")
include("denoise.jl")
include("segment.jl")

export identify, denoise, compute_flatfield

end # module
