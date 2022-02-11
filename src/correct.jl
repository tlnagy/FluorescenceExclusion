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

function correct!(data::A; mask=zeros(Bool, size(data)...)) where {A <: AbstractMatrix{<: AbstractGray{<: AbstractFloat}}}
    pillar_mask, pillar_centers = identify(Pillars(), data)
    cell_mask = identify(Cells(), data, pillar_mask) .| mask

    bkg = darkfield(data, pillar_centers)

    flatfield = compute_flatfield(data, cell_mask .| pillar_mask, len=20) .- bkg;

    labeled = label_components(cell_mask)

    data .= @. clamp((data - bkg) / flatfield, Gray(-0.05), Gray(1.05))

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