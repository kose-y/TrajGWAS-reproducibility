using DataFrames, CSV, Dates, CategoricalArrays
using Statistics
using vGWAS
using Ipopt, WiSER
using LinearAlgebra
using KNITRO
using BGEN
# fit the null model
BLAS.set_num_threads(1)
solver = KNITRO.KnitroSolver(outlev=3) # outlev 0-6

bp_data = CSV.read("bp_all.csv", DataFrame)

@time nm = vgwas(@formula(std_pbp ~ 1 + SEX + std_age + std_age_sq + std_age&SEX +
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
