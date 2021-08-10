using DataFrames, CSV
using Statistics, WiSER
using LinearAlgebra
using KNITRO
import StatsBase: countmap
# fit the null model
solver = KNITRO.KnitroSolver(outlev=0) # outlev 0-6

fg_data = CSV.read("glucose_fasting.csv", DataFrame)
diabetics = CSV.read("../diabetics.txt", DataFrame)[!, :FID]
using DataStructures
diabetics_map = DefaultDict{Int, Bool}(false)
for v in diabetics
    diabetics_map[v] = true
end
fg_data.diabetes_status = map(x -> Int64.(diabetics_map[x.FID]), eachrow(fg_data))
#self insulin just in beta
@time nm = WSVarLmmModel(@formula(std_glucose_fast ~ 1 + SEX + std_age + std_age_sq + std_bmi + SEX & std_age + std_bmi & std_age + self_insulin + diabetes_status +
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10),
    @formula(std_glucose_fast ~ 1 + std_age),
    @formula(std_glucose_fast ~ 1 + SEX + std_age + std_age_sq +
        std_bmi + SEX & std_age + std_bmi & std_age + diabetes_status),
    :FID,
    fg_data)
solver = KNITRO.KnitroSolver(outlev=0,ftol=1e-5)
@time fit!(nm,
solver,
runs=10, init = init_ls!(nm, gniters=2))

println(nm)
using Serialization
open("fittednullmodel.glucose_fast.indicatorbeta_withdiabetesstatus.allsubjs.jls", "w") do io
    Serialization.serialize(io, nm)
end

nfg_data = CSV.read("glucose_nonfasting.csv", DataFrame)
diabetics = CSV.read("../diabetics.txt", DataFrame)[!, :FID]
using DataStructures
diabetics_map = DefaultDict{Int, Bool}(false)
for v in diabetics
    diabetics_map[v] = true
end
nfg_data.diabetes_status = map(x -> Int64.(diabetics_map[x.FID]), eachrow(nfg_data))
#self insulin just in beta
@time nm = WSVarLmmModel(@formula(std_glucose_nonfast ~ 1 + SEX + std_age + std_age_sq + std_bmi + SEX & std_age + std_bmi & std_age + self_insulin + diabetes_status +
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10),
    @formula(std_glucose_nonfast ~ 1 + std_age),
    @formula(std_glucose_nonfast ~ 1 + SEX + std_age + std_age_sq +
        std_bmi + SEX & std_age + std_bmi & std_age + diabetes_status),
    :FID,
    nfg_data)
solver = KNITRO.KnitroSolver(outlev=0,ftol=1e-5)
@time fit!(nm,
solver,
runs=10)
println(nm)
using Serialization
open("fittednullmodel.glucose_nonfast.indicatorbeta_withdiabetesstatus.allsubjs.jls", "w") do io
    Serialization.serialize(io, nm)
end