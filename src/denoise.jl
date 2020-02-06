
"""
    denoise(filepath; cores)

Denoise a 3D TIFF image using the ndsafir tool from PRIISM. Assumes dimension
order is XYT.
"""
function denoise(filepath; cores=Sys.CPU_THREADS)
    filename = basename(filepath)
    external_progs = ["tiff2mrc", "ndsafir_priism", "mrc2tiff"]
    if any(isnothing(Sys.which.(external_progs)))
        @error "Make sure that $external_progs are all in the PATH"
    end
    if !endswith(filepath, r".tif|.tiff")
        @error "Given file must be a tif file"
    end
    # switch file ending to MRC
    mrc_filepath = reduce(replace, [".tiff"=>".mrc", ".tif"=>".mrc"], init=filepath)
    outfile = replace(mrc_filepath, ".mrc" => ".out")

    @info "Converting $filename for processing "
    run(`tiff2mrc -in=$filepath -palette=gray -multi=t $mrc_filepath`)

    # denoise image using NDSAFIR
    dn_mrc = replace(mrc_filepath, ".mrc" => "_dn.mrc")
    @info "Denoising $filename)"
    cmd = `ndsafir_priism "$mrc_filepath" "$dn_mrc" -iter=4 -np=$cores -3d=t -noise=poisson -adapt=10.0 -island=4.0 -p=1 -sampling=-1`
    run(pipeline(cmd, stdout=outfile, stderr=outfile, append=true))

    dn_tiff = replace(dn_mrc, "_dn.mrc" => "_dn.tiff")
    @info "Converting $filename back to tif"
    run(`mrc2tiff "$dn_mrc" -single -out="$dn_tiff"`)
end