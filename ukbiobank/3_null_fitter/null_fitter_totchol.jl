using DataFrames, CSV
using Statistics
using Ipopt, WiSER
using LinearAlgebra
using KNITRO
using BGEN
# fit the null model
BLAS.set_num_threads(1)
solver = KNITRO.KnitroSolver(outlev=0) # outlev 0-6



totchol_data_d = CSV.read("totchol_diabetics.csv", DataFrame)
totchol_data_nd = CSV.read("totchol_nondiabetics.csv", DataFrame)
totchol_data = vcat(totchol_data_d, totchol_data_nd)

@time nm = WSVarLmmModel(@formula(std_totchol ~ 1 + SEX + std_age + std_age_sq + std_bmi + SEX & std_age + std_bmi & std_age +
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10),
    @formula(std_totchol ~ 1 + std_age),
    @formula(std_totchol ~ 1 + SEX + std_age + std_age_sq +
        std_bmi + SEX & std_age + std_bmi & std_age),
    :IID,
    totchol_data)

fit!(nm,
solver,
runs=5)

println(nm)
using Serialization
open("fittednullmodel.totchol.test.allsubjs.jls", "w") do io
    Serialization.serialize(io, nm)
end