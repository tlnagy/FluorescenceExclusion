module FluorescenceExclusion

using DataStructures
using Distributions
using ImageBinarization
using ImageFiltering
using ImageSegmentation
using Images
using ProgressMeter
using ScatteredInterpolation
using Statistics

include("flatfield.jl")
include("denoise.jl")
include("segment.jl")
include("correct.jl")
include("tracking.jl")

export identify, segment!, correct!, link, denoise, compute_flatfield

end # module
