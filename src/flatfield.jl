const BUFFER = 5

"""
    compute_flatfield(slice, mask; len) -> flatfield

Given a 2D fluorescence image and a boolean mask where all foreground objects
are labeled true, this function returns a interpolated flatfield image to
account for inhomogeneities in illumination. It does this by fitting a
Multiquadratic Radial Basis function to a set of background grid points and
interpolating the rest. This works well for gradual changes in topology so any
sharp transitions in the background illumination will increase the error of the
interpolant.
"""
function compute_flatfield(slice::AbstractArray{T, 2}, mask::AbstractArray{Bool, 2}; len=10) where {T}

    points = generate_sample_grid(mask; len=len)

    values = vec(Float64.(slice[points]))
    spoints = hcat([[pos[1], pos[2]] for pos in points]...)

    itp = ScatteredInterpolation.interpolate(Multiquadratic(), spoints, values);

    # evaluate the interpolant at every point
    eval_points = hcat([[i, j] for i in 1:1024, j in 1:1024]...)

    flatfield = reshape(ScatteredInterpolation.evaluate(itp, eval_points), 1024, :)

    # make extra smooth
    imfilter!(flatfield, flatfield, Kernel.gaussian(30))

    flatfield
end

"""
    generate_sample_grid(mask)

Create a sinusoidally spaced grid that avoiding areas labeled true in `mask`.
The higher density of sampling points near the edges helps with the increased
steepness in signal loss and interpolation error that occurs at the image
boundaries.
"""
function generate_sample_grid(mask::AbstractArray{Bool, 2}; len = 10)
    width, height = size(mask)
    sinusoidal_spacing = sin.(range(-π/2, stop=π/2, length=len))
    xrange = round.(Int, sinusoidal_spacing .* (width/2 - BUFFER) .+ width/2)
    yrange = round.(Int, sinusoidal_spacing .* (height/2 - BUFFER) .+ height/2)
    gridpoints = Set(CartesianIndex(x, y) for x in xrange, y in yrange)

    # only select positions on the grid that are also far enough way from objects
    safe_areas = distance_transform(feature_transform(mask)) .> 10
    collect(intersect(gridpoints, Set(findall(safe_areas))))
end

"""
    display_grid(arr)

Takes the output of `generate_sample_grid` and creates a simple visualization of
the sample points
"""
function display_grid(arr::AbstractArray{CartesianIndex{2}})
    eval_points = falses((maximum(arr) + CartesianIndex(BUFFER, BUFFER)).I)
    eval_points[arr] .= true
    Gray.(distance_transform(feature_transform(eval_points)) .< 5)
end


function get_medians(img::AbstractArray{T, 2}, centers::AbstractArray{Bool, 2}) where {T}
    pillar_medians = Dict{Int, Float64}()
    pillar_labels = label_components(copy(centers))

    for i in unique(pillar_labels)
        (i == 0) && continue # ignore the background
        val = Float64(median(view(img, pillar_labels .== i)))
        pillar_medians[i] = val
    end
    pillar_medians
end


function get_medians(img::AbstractArray{T, 3}, centers::AbstractArray{Bool, 3}) where {T}
    pillar_medians = DefaultDict{Int, Vector{Float64}}(Vector{Float64})

    @showprogress for t in 1:size(img, 3)
        medians = get_medians(view(img, :, :, 1), view(centers, :, :, 1))
        for (k, v) in medians
            push!(pillar_medians[k], v)
        end
    end
    pillar_medians
end