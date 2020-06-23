using AxisArrays
using Distributions
using FixedPointNumbers
using Unitful: μm, ms, @u_str

@testset "Extract information from images" begin
    img = load(joinpath("testdata", "191108_fxm_high_spatiotemporal_res_gfp_corrected_t1.tif"));

    labels_raw = FileIO.load(joinpath("testdata", "191108_fxm_high_spatiotemporal_res_labels_t1.tif"))
    # convert labels to an integer matrix
    labels = Int.(reinterpret.(UInt8, convert.(N0f8, labels_raw)))

    axs = (Axis{:y}(0.0μm:0.6518μm:666.7914000000001μm), Axis{:x}(0.0μm:0.6518μm:666.7914000000001μm), Axis{:time}(0.0ms:2000.0ms:0.0ms))

    linked = FluorescenceExclusion.link(AxisArray(reshape(img, size(img)..., 1), axs...), AxisArray(reshape(labels, size(labels)..., 1), axs...))

    # select non-dead cells
    nondead = linked[linked[!, :abs_volumes] .> 250u"μm^3", :abs_volumes]

    mean_vol = mean(nondead)

    # check if the mean volume of cells is close to 510 microns
    @test isapprox(mean_vol, 510μm^3, rtol=0.01)
end