# FluorescenceExclusion.jl

A metapackage for analyzing fluorescence exclusion microscopy (FxM) images in
Julia. Very much a WIP.

## TL;DR

```@setup 1
using Pkg
using TiffImages
using FluorescenceExclusion

path = joinpath(Pkg.pkgdir(FluorescenceExclusion), "test", "testdata", "220125_lane2_fxmraw.tif")
```

This package provides tools to segment and correct fluorescence exclusion
microscopy images:

```@example 1
img = TiffImages.load(path) # load a FxM image

fimg = float.(img) # convert to floating point

correct!(fimg, fimg)

hcat(img, fimg)
```

## Public

```@docs
correct!
identify
build_tp_df
```

## Internal

```@docs
FluorescenceExclusion.compute_flatfield
FluorescenceExclusion.darkfield
FluorescenceExclusion.generate_sample_grid
FluorescenceExclusion.get_locality
FluorescenceExclusion.get_medians
FluorescenceExclusion.isencapsulated
```