
@testset "Segmenting pillars" begin
    img = ones(30,40, 2)

    img[5:10, 5:10, 1] .= 0.0

    # edge pillar, which should get clipped
    img[25:30, 35:40, 1] .= 0.0
    pillars, centers = identify(FluorescenceExclusion.Pillars(), img[:, :, 1], dist=1)

    # pillars should be the same as image except for the 11 pixels
    # on the bottom and right of the edge pillar
    @test sum(pillars .!== (img[:, :, 1] .== 0.0)) == 11

    # we should erode the non-edge pillar by two to 4x4, while the edge pillar
    # also loses an extra pixel due to the edge
    @test sum(centers) == 4*4 + 3*3

    pillars, centers = identify(FluorescenceExclusion.Pillars(), img, dist=1, verbose=false)

    @test sum(centers[:, :, 2]) == 0
end