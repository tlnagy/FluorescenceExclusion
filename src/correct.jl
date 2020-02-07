function correct!(img::AbstractArray{T, 3}) where {T}
    flatfield = similar(view(img, :, :, 1))
    labels = Array{Int}(undef, size(img))
    n_pillars = 0
    prev_max::Int = 0

    @showprogress for t in 1:size(img, 3)
        slice = view(img, :, :, t)
        
        pillar_mask, pillar_centers = identify(Pillars(), slice)
        cell_mask = identify(Cells(), slice, pillar_mask)

        pillar_medians = get_medians(slice, pillar_centers)
        
        if t == 1
            flatfield .= compute_flatfield(slice, cell_mask .| pillar_mask, len=10);
            n_pillars = length(keys(pillar_medians))
            @info "Found $n_pillars pillars"
        elseif length(keys(pillar_medians)) != n_pillars
            @error "Number of pillars changed at frame $t, this is likely a major problem"
        end

        labeled = label_components(cell_mask)

        remove_small!(cell_mask, labeled)
        
        segment!(slice, labeled)

        slice .= imadjustintensity(slice ./ (flatfield .* mean(values(pillar_medians))))

        labeled[labeled .> 0] .+= prev_max
        prev_max = maximum(labeled)
        labels[:, :, t] .= labeled
    end

    labels
end