using DataFrames, CSV
using Statistics, WiSER
using LinearAlgebra
using KNITRO
import StatsBase: countmap
# fit the null model
solver = KNITRO.KnitroSolver(outlev=0) # outlev 0-6

tg_data = CSV.read("triglycerides_all.csv", DataFrame)
#self insulin  in beta and tau
@time nm = WSVarLmmModel(@formula(std_triglycerides ~ 1 + SEX + std_age + std_age_sq + std_bmi + SEX & std_age + std_bmi & std_age +
        + self_cholesteroldrugs + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10),
    @formula(std_triglycerides ~ 1 + std_age),
    @formula(std_triglycerides ~ 1 + SEX + std_age + std_age_sq +
        std_bmi + SEX & std_age + std_bmi & std_age + self_cholesteroldrugs),
    :FID,
    tg_data)

@time fit!(nm,
solver,
runs=10)

println(nm)
using Serialization
open("fittednullmodel.triglycerides.indicatorboth.allsubjs.jls", "w") do io
    Serialization.serialize(io, nm)
end