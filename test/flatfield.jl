using Distributions

@testset "Flatfield Correction" begin

    # lets create an off-center non-flat light source
    norm = pdf.(Normal(512, 512), 1:1024)
    img = imadjustintensity(norm * norm')

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