using Suppressor
using DataFrames 
using Unitful: μm, ustrip
using AxisArrays
using PyCall
using Pandas

function link(img, labels; dist=(3, 12), chamber_height=12.96μm)
    particle_df = build_tp_df(img, labels, dist=dist);

    to_track = particle_df[!, [:frame, :id, :x, :y]]

    select!(particle_df, Not([:x, :y]))

    tp = pyimport("trackpy")

    # convert particles to a Pandas DataFrame and call trackpy

    @info "Linking..."
    t = @suppress_out tp.link_df(Pandas.DataFrame(to_track), 25, memory=2);
    linked = join(DataFrames.DataFrame(Pandas.DataFrame(t)), particle_df, on=[:frame, :id])

    xstepsize = step(AxisArrays.axes(img, Axis{:x}).val)
    ystepsize = step(AxisArrays.axes(img, Axis{:y}).val)

    # compute the maximum intensity value for each cell's background
    Imax = linked[:, :medbkg_FxM]

    rel_vols = (linked[!, :medbkg_FxM] .* linked[!, :area]) .- linked[!, :tf_FxM];
    insertcols!(linked, 6, :rel_volume => rel_vols)

    # α is the maximum signal per voxel
    α = (Imax .- 0.0) ./ (chamber_height .* xstepsize .* ystepsize)

    # equation 2 from Cadart 2017
    abs_vols = linked[!, :rel_volume] ./ α
    insertcols!(linked, 7, :abs_volume => abs_vols)

    linked[!, :area] = linked[!, :area] .* xstepsize .* ystepsize

    linked
end