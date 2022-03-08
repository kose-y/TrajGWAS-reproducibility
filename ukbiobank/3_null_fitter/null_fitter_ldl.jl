using DataFrames, CSV
using Statistics
using Ipopt, WiSER
using LinearAlgebra
using KNITRO
using BGEN
# fit the null model
BLAS.set_num_threads(1)
solver = KNITRO.KnitroSolver(outlev=0) # outlev 0-6
# can use other solvers e.g., 
# Ipopt.IpoptSolver(print_level=0, mehrotra_algorithm = "yes", warm_start_init_point="yes", max_iter=100)
# (for Ipopt version 0.8)
# KNITRO (commercial software) is more efficient/stable/robust than other solvers



ldl_data_d = CSV.read("ldl_diabetics.csv", DataFrame)
ldl_data_nd = CSV.read("ldl_nondiabetics.csv", DataFrame)
ldl_data = vcat(ldl_data_d, ldl_data_nd)
genetic_iids_subsample = unique(ldl_data.FID)
@time nm = WSVarLmmModel(@formula(std_ldl ~ 1 + SEX + std_age + std_age_sq + std_bmi + SEX & std_age + std_bmi & std_age +
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10),
    @formula(std_ldl ~ 1 + std_age),
    @formula(std_ldl ~ 1 + SEX + std_age + std_age_sq +
        std_bmi + SEX & std_age + std_bmi & std_age),
    :IID,
    ldl_data)

fit!(nm,
solver,
runs=5)


println(nm)
using Serialization
open("fittednullmodel.ldl.test.allsubjs.jls", "w") do io
    Serialization.serialize(io, nm)
end