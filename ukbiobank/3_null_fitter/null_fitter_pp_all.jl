using DataFrames, CSV, Dates, CategoricalArrays
using Statistics
using TrajGWAS
using Ipopt, WiSER
using LinearAlgebra
using KNITRO
using BGEN
# fit the null model
BLAS.set_num_threads(1)
solver = KNITRO.KnitroSolver(outlev=3) # outlev 0-6
# can use other solvers e.g., 
# Ipopt.IpoptSolver(print_level=0, mehrotra_algorithm = "yes", warm_start_init_point="yes", max_iter=100)
# (for Ipopt version 0.8)
# KNITRO (commercial software) is more efficient/stable/robust than other solvers

bp_data = CSV.read("bp_all.csv", DataFrame)

@time nm = trajgwas(@formula(std_pbp ~ 1 + SEX + std_age + std_age_sq + std_age&SEX +
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10),
    @formula(std_pbp ~ 1 + std_age),
    @formula(std_pbp ~ 1 + SEX + std_age + std_age_sq + std_age&SEX),
    :IID,
    bp_data,
    nothing;
    nullfile="pbp.run5_final.all.null.txt",
    solver=solver,
    runs=10
)

println(nm)
using Serialization
open("fittednullmodel.pbp.run5_final.all.jls", "w") do io
    Serialization.serialize(io, nm)
end
