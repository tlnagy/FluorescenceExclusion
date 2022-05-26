using StatsBase

function correct!(data::AbstractArray{<: Colorant, 3}; mask = nothing)
    labels = Array{Int}(undef, size(dn_img))
    (mask !== nothing) && (@assert size(data) == size(mask))

    @showprogress for (t, slice) in enumerate(eachslice(data, dims=3))
        if mask !== nothing
            labels[:, :, t] = correct!(slice; mask=view(mask, :, :, t))
        else
            labels[:, :, t] = correct!(slice)
        end
    end
    labels
end

"""
    correct!(data; mask)

Given a FxM image `data`, corrects for inhomogeneities both in illumination and
in warping. 

It does this by sampling the foreground signal of the excluded dye
and then interpolating into the "holes" in the signal like the pillars and
cells. Passing an additional optional parameter `mask` can augment the
built-in detection of cells to avoid using those areas to estimate the
foreground. Additionally, it estimates the true background by interpolating
between pillars.

The returned values should be centered around 0 for the centers of the pillars
and around 1 for the foregound areas.

```jldoctest; setup = :(using TiffImages, Pkg, FluorescenceExclusion; path = joinpath(Pkg.pkgdir(FluorescenceExclusion), "test", "testdata", "220125_lane2_fxmraw.tif"))
julia> img = TiffImages.load(path);

julia> fimg = float.(img); #convert to Gray{Float32}

julia> correct!(fimg); # correct FxM channel
```
"""
function correct!(data::A, dn_img) where {A <: AbstractMatrix{<: AbstractGray{<: AbstractFloat}}}
    pillar_mask, pillar_centers = identify(Pillars(), dn_img)
    cell_mask = identify(Cells(), dn_img, pillar_mask)

    bkg = darkfield(data, pillar_centers)

    data .-= bkg

    flatfield = compute_flatfield(data, cell_mask .| pillar_mask, len=30);

    labeled = label_components(cell_mask)

    data .= @. clamp(data / flatfield, Gray(-0.1), Gray(1.1))

    labeled
end

function fix_labels(labels)
    prev_max = 0
    for t in 1:size(labels, 3)
        labeled = view(labels, :, :, t)
        labeled[labeled .> 0] .+= prev_max
        prev_max = maximum(labeled)
    end
end