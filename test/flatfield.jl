using Distributions
using ImageContrastAdjustment
using ImageMorphology

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

    vals = FluorescenceExclusion.get_medians(img, centers; verbose=false)

    @test vals[1] == [2.0, 6.0, 6.0]
    @test vals[2] == [2.5]
    @test vals[3] == [6.0]

    centers = repeat(img[:, :, 1] .> 0, 1, 1, 3)

    vals = FluorescenceExclusion.get_medians(img, centers; verbose=false)

    @test vals[1] == [2.0, 0.0, 0.0]
    @test vals[2] == [2.5, 0.0, 0.0]
    @test vals[3] == [6.0, 6.0, 6.0]
end