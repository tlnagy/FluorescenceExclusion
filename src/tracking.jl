using SegmentationTools
using Suppressor
using DataFrames
using Unitful: μm, ustrip
using AxisArrays
using Query
using PyCall
import Pandas

function link(img, labels; dist=(5, 50), chamber_height=12.96μm)
    particle_df = SegmentationTools.build_tp_df(img, labels, dist=dist);

    tp = pyimport("trackpy")

    # convert particles to a Pandas DataFrame and call trackpy

    @info "Linking..."
    particle_df[!, :area] .= ustrip.(particle_df[!, :area])
    t = @suppress_out tp.link_df(Pandas.DataFrame(particle_df), 25, memory=2);
    linked = DataFrame(Pandas.DataFrame(t));
    linked[!, :area] .= linked[!, :area] .* μm^2
    volumes = linked[!, :bkg_slice] .- linked[!, :tf_slice];

    xstepsize = step(AxisArrays.axes(img, Axis{:x}).val)
    ystepsize = step(AxisArrays.axes(img, Axis{:x}).val)

    # compute the maximum intensity value for each cell's background
    Imax = linked[:, :bkg_slice] ./ (linked[:, :area] ./ (xstepsize * ystepsize));

    # equation 2 from Cadart 2017
    linked[!, :abs_volumes] .= volumes .* (Imax .- 0.0) .* chamber_height .* xstepsize .* ystepsize

    linked
end