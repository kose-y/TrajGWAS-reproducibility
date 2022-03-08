using DataFrames, CSV
using Statistics, WiSER
using LinearAlgebra
using KNITRO
import StatsBase: countmap
# fit the null model
solver = KNITRO.KnitroSolver(outlev=0) # outlev 0-6
# can use other solvers e.g., 
# Ipopt.IpoptSolver(print_level=0, mehrotra_algorithm = "yes", warm_start_init_point="yes", max_iter=100)
# (for Ipopt version 0.8)
# KNITRO (commercial software) is more efficient/stable/robust than other solvers

hba1c_data = CSV.read("a1c_all.csv", DataFrame)
diabetics = CSV.read("../diabetics.txt", DataFrame)[!, :FID]
using DataStructures
diabetics_map = DefaultDict{Int, Bool}(false)
for v in diabetics
    diabetics_map[v] = true
end
hba1c_data.diabetes_status = map(x -> Int64.(diabetics_map[x.IID]), eachrow(hba1c_data))

# self insulin and beta and tau 
@time nm = WSVarLmmModel(@formula(std_a1c ~ 1 + SEX + std_age + std_age_sq + std_bmi + std_bmi & std_age +
        self_insulin + diabetes_status + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10),
    @formula(std_a1c ~ 1 + std_age),
    @formula(std_a1c ~ 1 + SEX + std_age + std_age_sq + diabetes_status + self_insulin +
        std_bmi + std_bmi & std_age),
    :FID,
    hba1c_data)

@time fit!(nm, solver, init = init_ls!(nm, gniters = 0), runs = 10) 
println(nm)
using Serialization
open("fittednullmodel.hba1c.indicatorboth.allsubjs.jls", "w") do io
    Serialization.serialize(io, nm)
end