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

bmi_data = CSV.read("BMI_PC.csv", DataFrame)

@time nm = vgwas(@formula(std_BMI ~ 1 + SEX  + std_age + std_age_sq + std_age&SEX +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10),
    @formula(std_BMI ~ 1 + std_age),
    @formula(std_BMI ~ 1 + SEX + std_age + std_age_sq + std_age&SEX),
    :FID,
    bmi_data,
    nothing;
    nullfile="bmi.all.null.txt",
    solver=solver,
    runs=10, 
    init_mom=true
)

println(nm)
using Serialization
open("fittednullmodel.bmi.all.jls", "w") do io
    Serialization.serialize(io, nm)
end
