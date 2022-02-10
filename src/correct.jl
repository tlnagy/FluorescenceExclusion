using StatsBase

function correct!(data::AbstractArray{<: Colorant, 3}, dn_img::AbstractArray{<: Colorant, 3}; mask = nothing)
    labels = Array{Int}(undef, size(dn_img))
    (mask !== nothing) && (@assert size(data) == size(mask))

    @showprogress for t in 1:size(dn_img, 3)
        slice = view(data, :, :, t)
        dn_slice = view(dn_img, :, :, t)

        if mask !== nothing
            labels[:, :, t] = correct!(slice, dn_slice; mask=mask[:, :, t])
        else
            labels[:, :, t] = correct!(slice, dn_slice)
        end
    end
    labels
end

function correct!(data::AbstractArray{<: Colorant, 2}, dn_img::AbstractArray{<: Colorant, 2}; mask=zeros(Bool, size(data)...)) 
    pillar_mask, pillar_centers = identify(Pillars(), dn_img)
    cell_mask = identify(Cells(), dn_img, pillar_mask) .| mask

    pillar_medians = get_medians(data, pillar_centers)
    bkg = mean(values(pillar_medians))

    flatfield = compute_flatfield(data, cell_mask .| pillar_mask, len=20) .- bkg;

    labeled = label_components(cell_mask)


    # remove_small!(cell_mask, labeled)
    # segment!(data, labeled)

    data .= (data ./ flatfield)
    data ./= percentile(vec(data), 99.9)

    pillar_medians = get_medians(data, pillar_centers)
    bkg = mean(values(pillar_medians))

    # subtract the background and clamp values so that the median background
    # value is now 0
    data .= clamp01.(data .- bkg)

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