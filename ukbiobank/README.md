## Scripts for UK Biobank experiments

- Launch julia with `--project=<path-to-this directory>` or activate the environment with `using Pkg; Pkg.activate(normpath("/path/to/this/directory"))`.
- Then run `using Pkg; Pkg.instantiate()`. 

Note: The package was named `vGWAS` at the time of experiments. Using the commands above sets the name of the package as `vGWAS`. The package names should be changed to `TrajGWAS`, and the function `vgwas` should be changed to `trajgwas` for the use with the current version of `TrajGWAS.jl`.

