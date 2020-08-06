using DataFrames

"""
    build_tp_df(img, components; dist)

Build a `DataFrames.DataFrame` that is compatible with trackpys `link_df`
function. Needs to be converted to a `Pandas.DataFrame` before passing to
trackpy. `dist` is a 2-tuple of integers indicating the minimum and maximum
distance away in pixels from each cell to include in its local background
calculation.
"""
function build_tp_df(img::AxisArray{T1, 4},
                     components::AxisArray{T2, 3}; dist=(2, 10)) where {T1, T2 <: Integer}

    particles = DataFrames.DataFrame[]

    for (idx, timepoint) in enumerate(timeaxis(img))
        component_slice = view(components, :, :, idx)

        # get all ids present in the current slice
        slice_ids = sort(unique(component_slice))
        # get the number of pixels for each id in the slice
        lengths = component_lengths(component_slice)[slice_ids .+ 1]

        correct_size = trues(length(lengths))
        correct_size[1] = false

        ids = slice_ids[correct_size]
        # get the centroids for all ids and then select only the ids that appear
        # in the current slice and *also* have the correct size
        centroids = component_centroids(component_slice)[slice_ids[correct_size] .+ 1]

        n = length(centroids)
        ys = map(f->f[1], centroids) # y corresponds to rows
        xs = map(f->f[2], centroids) # x corresponds to columns
        frames = fill(idx-1, n)

        # select the current time slice and enforce storage order to match
        # components
        slice = view(img, Axis{:y}(:), Axis{:x}(:), Axis{:channel}(:), Axis{:time}(timepoint))
        localities = identify(Locality(), component_slice, ids, dist=dist)
        # dictionary of ids to areas
        data = OrderedDict(:x=>xs,
                           :y=>ys,
                           :frame=>frames,
                           :id=>ids,
                           :area=>lengths[2:end],
                           :footprint=>map(x->CartesianIndex.(x), (component_subscripts(component_slice)[ids .+ 1])),
                           :locality=>collect(values(localities)))

        cax = AxisArrays.axes(img, Axis{:channel})
        for c in cax
            channelslice = view(slice, Axis{:channel}(c))
            tfs = Float64[]
            medbkgs = Measurement{Float64}[]

            for (id, indices) in localities
                # total fluorescence is the sum of all signal in the actual
                # footprint of the object. We have to do a copy operation here
                # due to https://github.com/JuliaArrays/AxisArrays.jl/issues/179
                push!(tfs, sum(channelslice[component_slice .== id]))

                # median background is the median of background signal in the
                # locality of object
                bkg = Float64.(channelslice[indices])
                push!(medbkgs, median(bkg) ± (1.253 * std(bkg) / sqrt(length(bkg))))
            end
            data[Symbol("tf_", c)] = tfs
            data[Symbol("medbkg_", c)] = medbkgs
        end
        push!(particles, DataFrames.DataFrame(data))
    end
    vcat(particles...)
end

function build_tp_df(img::AxisArray{T1, 3},
                     thresholds::AxisArray{T2, 3}; dist=(2, 10)) where {T1, T2 <: Bool}
    build_tp_df(AxisArray(reshape(img, size(img)..., 1),
                          AxisArrays.axes(img)..., Axis{:channel}([:FxM])),
                thresholds;
                dist=dist
               )
end

function build_tp_df(img::AxisArray{T1, 4},
                     thresholds::AxisArray{T2, 3}; dist=(2, 10)) where {T1, T2 <: Bool}

    # we have to pass the underlying array due to
    # https://github.com/JuliaImages/ImageMorphology.jl/issues/21
    components = Images.label_components(thresholds.data, [axisdim(thresholds, Axis{:y}), axisdim(thresholds, Axis{:x})])
    build_tp_df(img, AxisArray(components, AxisArrays.axes(thresholds)); dist=dist)
end

function build_tp_df(img::AxisArray{T1, 3},
                     thresholds::AxisArray{T2, 3}; dist=(2, 10)) where {T1, T2 <: Integer}

    build_tp_df(AxisArray(reshape(img, size(img)..., 1),
                          AxisArrays.axes(img)..., Axis{:channel}([:FxM])),
                thresholds;
                dist=dist
               )
end

"""
Given an `img` with at least `y`, `x`, and `t` axes and a 3 dimensional boolean
array, `thresholds`, in yxt order.
"""
function build_tp_df(img::AxisArray{T1, 4},
                     thresholds::BitArray{3}; dist=(2, 10)) where {T1}
    _axes = Tuple(AxisArrays.axes(img, Axis{ax}) for ax in (:y, :x, :time))
    build_tp_df(img, AxisArray(Bool.(thresholds), _axes...), dist=dist)
end