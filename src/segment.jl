
abstract type AbstractSegmentable end

struct Pillars <: AbstractSegmentable end

"""
    identify(Pillars(), img; dist)

Use a simple algorithm to identify pillars and their centers. The optional
keyword `dist` allows varying the number of pixels from the edge of the
pillar that a pixel must be to be considered a center pixel. This helps correct
for incomplete sealing at the edges of pillars.
"""
function identify(::Pillars, img::AbstractArray{T, 2}; dist=30) where {T}
    M₁ = opening(binarize(Bool, img, ImageBinarization.Balanced()))

    h, w = size(img)

    # prevent edge pixels from being included in pillars
    M₁[1:h, [1, w]] .= true
    M₁[[1, h], 1:w] .= true

    dtrans = distance_transform(feature_transform(M₁))

    M₁ .= .~ M₁

    M₁, dtrans .> dist
end

function identify(::Pillars, img::AbstractArray{T, 3}; dist=30, verbose=true) where {T}
    pillar_masks = Array{Bool}(undef, size(img))
    pillar_centers = Array{Bool}(undef, size(img))

    wait_time = verbose ? 1 : Inf

    @showprogress wait_time for t in 1:size(img, 3)
        I = @view img[:, :, t]
        out = identify(Pillars(), I)

        pillar_masks[:, :, t] .= out[1]
        pillar_centers[:, :, t] .= out[2]
    end
    pillar_masks, pillar_centers
end

struct Cells <: AbstractSegmentable end

"""
    identify(Cells(), img, pillar_masks)

Identify cells in the FxM channel using a custom edge-based detection system that
thresholds using the distribution of spatial gradient sizes across cells.
"""
function identify(::Cells, 
                 img::AbstractArray{T, 2},
                 pillar_masks::AbstractArray{Bool, 2}) where {T}

    grad_y, grad_x, mag, orient = imedge(img, KernelFactors.scharr);

    # remove the pillars and their vicinities from the calculations since they
    # can skew the results
    foreground = distance_transform(feature_transform((pillar_masks))) .> 30 
    mag .*= foreground

    flattened = filter(x->x > 0.0, reshape(mag, :))
    lo, hi = quantile(flattened, (0.01, 0.99))

    # the magnitudes are almost log normal distributed, so we fit the
    # distribution then use 1 σ above the mean as the cutoff. A quick comparison
    # suggests that this approach catches the thin edges of cells much better
    # than common algorithms like Otsu or Yen.
    normfit = fit_mle(LogNormal, view(flattened, lo .< flattened .< hi))
    
    thres = Images.opening(.~imfill(mag .< exp(normfit.μ + 1*normfit.σ), (0, 500)))
    imfill(thres, (0, 50))
end

function identify(::Cells, 
                 img::AbstractArray{T, 3},
                 pillar_masks::AbstractArray{Bool, 3}) where {T}

    @assert size(img) == size(pillar_masks)
    masks = Array{Bool}(undef, size(img))
    @showprogress for i in 1:size(img, 3)
        masks[:, :, i] .= identify(Cells(), view(img, :, :, i), view(pillar_masks, :, :, i))
    end
    masks
end

function remove_small!(mask::AbstractArray{Bool, 2}, labels::AbstractArray{Int, 2}; lim=100) where {T}
    for label in findall(component_lengths(labels) .< lim) .- 1
        mask[labels .== label] .= false
    end
end

function remove_small!(mask::AbstractArray{Bool, 2}; lim=100) where {T}
    # TODO: we shouldn't copy here but it needs to wait till
    # https://github.com/JuliaImages/ImageMorphology.jl/issues/21 is resolved
    labels = label_components(copy(mask))
    remove_small!(mask, labels; lim=lim)
end

function segment!(img::AbstractArray{T, 2}, label::AbstractArray{Int, 2}) where {T}
    seeds = label
    bkg = label .== 0
    # seeds will be computed by eroding 8 pixels in from the background to get
    # an idea of where the cells are 
    seeds[distance_transform(feature_transform(bkg)) .<= 8] .= 0

    segments = ImageSegmentation.watershed(img, seeds; mask=.~ bkg, compactness=0.01)
    label .= labels_map(segments)
end