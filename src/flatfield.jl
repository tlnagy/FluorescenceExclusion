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

    xrange = round.(Int, range(5, stop=size(slice, 1)-5, length=len))
    yrange = round.(Int, range(5, stop=size(slice, 2)-5, length=len))
    gridpoints = Set(CartesianIndex(x, y) for x in xrange, y in yrange)

    # only select positions on the grid that are also far enough way from objects
    safe_areas = distance_transform(feature_transform(mask)) .> 5
    points = collect(intersect(gridpoints, Set(findall(safe_areas))))

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