# TrajGWAS reproducibility 

Scripts to reproduce the experiments performed for the TrajGWAS manuscript:

_Ko S, German CA, Jensen A, Shen J, Wang A, Mehrotra DV, Sun YV, Sinsheimer JS, Zhou H, Zhou JJ. GWAS of longitudinal trajectories at biobank scale. The American Journal of Human Genetics. In press. [doi:10.1016/j.ajhg.2022.01.018](https://doi.org/10.1016/j.ajhg.2022.01.018)._

- `simulation` : Scripts for reproducing simulation studies. Experiments are conducted with an older version of `TrajGWAS`, named `vGWAS`. This version of package is kept on the branch [`gwas_dev`](https://github.com/OpenMendel/TrajGWAS.jl/tree/gwas_dev) of the package repository. You can load this version following the instruction inside the [directory](https://github.com/kose-y/TrajGWAS-reproducibility/tree/main/simulation).
- `ukbiobank` : Scripts for UK Biobank data analyses. It can run with `TrajGWAS` v0.1.2. Its optimization interface has changed to MathOptInterface since v0.4, and `KNITRO.KnitroSolver` no longer works. Use
```julia
solver = KNITRO.Optimizer()
solver_config = Dict("outlev" => 3)
```
as keyword arguments for the `trajgwas()` function instead to run with newer version of TrajGWAS and KNITRO.

