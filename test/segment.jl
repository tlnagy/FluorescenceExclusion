
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