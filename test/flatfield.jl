using Distributions
using ImageContrastAdjustment
using ImageMorphology
using FluorescenceExclusion: identify, Pillars, get_medians
using ColorTypes

@testset "Flatfield Correction" begin

    # lets create an off-center non-flat light source
    norm = pdf.(Normal(512, 512), 1:1024)
    img = adjust_histogram(norm * norm', LinearStretching(nothing => (0,1)))

    ranges = (x,y) -> (x-50:x+50, y-50:y+50)

    mask = falses(size(img))

    mask[ranges(512, 70)...] .= true
    mask[ranges(70, 512)...] .= true
    mask[ranges(512, 1024-70)...] .= true
    mask[ranges(1024-70, 512)...] .= true

    img_w_pillars = copy(img)
    img_w_pillars[mask] .= 0.0

    flatfield = FluorescenceExclusion.compute_flatfield(img_w_pillars, mask, len=20);

    # the very outer pixels contribute quite a bit of error
    @test all(isapprox.(flatfield, img, atol=0.03))

    # however, the error in the center should be very low
    @test all(isapprox.(flatfield[256:768, 256:768], img[256:768, 256:768], atol=0.006))

end

@testset "Temporal drift correction" begin
    img = zeros(10,10,3)
    img[2,2,1] = 2.0
    img[4,4,1] = 2.5
    img[6:9, 6, 1:3] .= [3.0, 4.0, 8.0, 9.0]
    centers = img .> 0.0

    vals = get_medians(img, centers; verbose=false)

    @test vals[1] == [2.0, 6.0, 6.0]
    @test vals[2] == [2.5]
    @test vals[3] == [6.0]

    centers = repeat(img[:, :, 1] .> 0, 1, 1, 3)

    vals = get_medians(img, centers; verbose=false)

    @test vals[1] == [2.0, 0.0, 0.0]
    @test vals[2] == [2.5, 0.0, 0.0]
    @test vals[3] == [6.0, 6.0, 6.0]
end

@testset "Darkfield correction" begin
    pillar_centers = Gray.(falses(512, 256))
    pillar_centers[128 .+ (-5:5), 128 .+ (-5:5)] .= true
    pillar_centers[384 .+ (-5:5), 128 .+ (-5:5)] .= true

    img = zeros(Gray{Float32}, size(pillar_centers)...)
    img[1:256, 1:256] .= Gray(1.0)

    bkg = FluorescenceExclusion.darkfield(img, pillar_centers)

    # make sure the pillars are correct
    @test bkg[128, 128] == Gray(1.0)
    @test bkg[384, 128] == Gray(0.0)

    # test that halfway between the two "pillars" it's exactly half of the signal
    @test bkg[256, 128] ≈ Gray(0.5)
end