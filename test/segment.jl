using ImageDraw
using FluorescenceExclusion: isencapsulated

@testset "Segmenting pillars" begin
    img = ones(30,40, 2)

    img[5:10, 5:10, 1] .= 0.0

    # edge pillar, which should get clipped
    img[25:30, 35:40, 1] .= 0.0
    pillars, centers = identify(FluorescenceExclusion.Pillars(), img[:, :, 1], dist=1, min_area = 0)

    # pillars should be the same as image except for the 11 pixels
    # on the bottom and right of the edge pillar
    @test sum(pillars .!== (img[:, :, 1] .== 0.0)) == 11

    # we should erode the non-edge pillar by two to 4x4, while the edge pillar
    # also loses an extra pixel due to the edge
    @test sum(centers) == 4*4 + 3*3

    pillars, centers = identify(FluorescenceExclusion.Pillars(), img, dist=1, min_area=0, verbose=false)

    @test sum(centers[:, :, 2]) == 0
end

@testset "Segmenting cells" begin
    img = load(joinpath("testdata", "191108_fxm_high_spatiotemporal_res_gfp_corrected_t1.tif"));

    pillars, centers = identify(FluorescenceExclusion.Pillars(), img)

    cell_mask = identify(FluorescenceExclusion.Cells(), img, pillars)

    labels = label_components(cell_mask)

    FluorescenceExclusion.remove_small!(cell_mask, labels)

    localities = distance_transform(feature_transform(cell_mask)) .< 15;

    locality_labels = label_components(localities)

    indices = component_indices(locality_labels)
    deleteat!(indices, 1)

    pillar_idxs = Set(findall(vec(pillars)))

    for idxs in indices
        # if there are more than one cells in this locality than remove it by setting it equal to the
        # background
        unq_idxs = unique(labels[idxs])
        if length(unq_idxs) != 2 # 2 because the background is also a label
            cell_mask[idxs] .= 0
        # make sure that the expanded area doesn't touch a pillar
        elseif length(intersect(pillar_idxs, Set(idxs))) > 0
            println("Touching a pillar!")
            cell_mask[idxs] .= 0
        end
    end

end

@testset "Encapsulation" begin
    img = Gray.(falses(100, 100))

    # first draw two circles so that we have a hole (a cell) that is fully
    # encapsulated by the true signal (aka the locality)
    draw!(img, Ellipse(CirclePointRadius(50, 50, 40)))
    draw!(img, Ellipse(CirclePointRadius(50, 50, 20)), Gray.(false))

    imgr = reinterpret(Bool, img)

    # should be encapsulated
    @test isencapsulated(findall(imgr))

    # Let's add a thin connector to connect the "hole" with the outside so that
    # it is no longer encapsulated
    imgr[50, 1:50] .= false

    @test !isencapsulated(findall(imgr)) # make sure we aren't encapsulated

    img .= false

    # Now, let's test a locality with a thickness of a single pixel, the
    # calculation of the convex hull can have a couple pixel error so this makes
    # sure it can handle thin regions
    draw!(img, Ellipse(CirclePointRadius(50, 50, 40)))
    draw!(img, Ellipse(CirclePointRadius(50, 50, 39)), Gray.(false))

    imgr = reinterpret(Bool, img)

    @test isencapsulated(findall(imgr))

    # Let's add a thin connector to connect the "hole" with the outside so that
    # it is no longer encapsulated
    imgr[50, 1:50] .= false

    @test !isencapsulated(findall(imgr)) # make sure we aren't encapsulated

    # if there are <3 points in the locality than the convex hull calculation
    # fails, lets make sure we handle this properly
    img .= false
    @test !isencapsulated(findall(imgr))


    draw!(img, Ellipse(CirclePointRadius(90, 10, 2)))
    draw!(img, Ellipse(CirclePointRadius(50, 80, 1)))

    # handle cases where the area of the convex hull is negative (i.e.
    # counterclockwise) 
    @test !isencapsulated(findall(imgr))
end