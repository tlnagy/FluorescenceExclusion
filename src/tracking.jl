using SegmentationTools
using Suppressor
using DataFrames
using Unitful: μm, ustrip
using AxisArrays
using Query
using PyCall
import Pandas

function link(img, labels; dist=(3, 12), chamber_height=12.96μm)
    particle_df = SegmentationTools.build_tp_df(img, labels, dist=dist);
    return particle_df

    tp = pyimport("trackpy")

    # convert particles to a Pandas DataFrame and call trackpy

    @info "Linking..."
    t = @suppress_out tp.link_df(Pandas.DataFrame(particle_df), 25, memory=2);
    linked = DataFrame(Pandas.DataFrame(t));

    xstepsize = step(AxisArrays.axes(img, Axis{:x}).val)
    ystepsize = step(AxisArrays.axes(img, Axis{:y}).val)

    # compute the maximum intensity value for each cell's background
    Imax = linked[:, :medbkg_slice]

    linked[!, :rel_volumes] = (linked[!, :medbkg_slice] .* linked[!, :area]) .- linked[!, :tf_slice];

    # α is the maximum signal per voxel
    α = (Imax .- 0.0) ./ (chamber_height .* xstepsize .* ystepsize)

    # equation 2 from Cadart 2017
    linked[!, :abs_volumes] .= linked[!, :rel_volumes] ./ α

    linked[!, :area] = linked[!, :area] .* xstepsize .* ystepsize

    linked
end