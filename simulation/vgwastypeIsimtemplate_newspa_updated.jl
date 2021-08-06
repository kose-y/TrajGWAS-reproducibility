#Code to do type I error simulation for vGWAS.jl

# Covariates input by runtypeIsims.jl
# index = index for that run 
# maf= minor allele freq p 
# ni = # obs per person, a range
# m = # number of individuals
println("index = $index, maf = $maf, ni = $ni, m = $m")
using vGWAS, WiSER, Distributions, Random, LinearAlgebra, StatsBase
using CodecZlib, Interpolations, Roots
import Optim: minimizer, optimize, LBFGS, NelderMead
BLAS.set_num_threads(1)

nistring = string(minimum(ni)) * "to" * string(maximum(ni))

# set parameters

pA = 1 - maf
paa = maf^2
pAa = 2 * maf *  pA
pAA = pA^2
alleles = [0.0, 1.0, 2.0]
pAlleles_p = [pAA, pAa, paa]

βtrue = [10.0; 5.0; 0.5; -0.3] 
τtrue = [0.25; 0.3; -0.15; 0.1] 
Σγ    = Matrix(Diagonal([2.0]))
δγω   = [0.0]
σω    = [0.1]

p, l, q = length(βtrue), length(τtrue), size(Σγ, 2)

Σγω   = [Σγ δγω; δγω' σω]
Lγω   = cholesky(Symmetric(Σγω), check = false).L
Lγ    = Lγω[1:q, 1:q]
lγω   = Lγω[q + 1, 1:q]
lω    = Lγω[q + 1, q + 1]

γω = Vector{Float64}(undef, q + 1)
z  = similar(γω) # hold vector of iid std normal

obsvec = Vector{WSVarLmmObs{Float64}}(undef, m)

# open p-value file 
gzio = GzipCompressorStream(open("pvals_type1_newspa_maf$(maf)_$(index)_ni$(nistring)_m$(m)_.csv.gz", "w"))
println(gzio, "betapval,taupval,jointpval,betapval_spa,taupval_spa,jointpval_spa")
seed = 818 + Int(maf * 1_000_000) + index + maximum(ni) * 100 + m
Random.seed!(seed)
Xbank = Dict{Int, Matrix{Float64}}()
Zbank = Dict{Int, Matrix{Float64}}()
Wbank = Dict{Int, Matrix{Float64}}()

for ni_i in ni
    Xbank[ni_i] = Matrix{Float64}(undef, ni_i, p)
    Zbank[ni_i] = Matrix{Float64}(undef, ni_i, q)
    Wbank[ni_i] = Matrix{Float64}(undef, ni_i, l)
    (Xbank[ni_i])[:, 1] .= 1
    (Zbank[ni_i])[:, 1] .= 1
    (Wbank[ni_i])[:, 1] .= 1
end

for i in 1:m
    ni_i = rand(ni)
    sex = rand(0:1)
    age_std = randn() #time-invariant
    bmi_std = randn(ni_i) #time-varying
    bp_std = randn(ni_i) #time-varying

    X = Xbank[ni_i]
    Z = Zbank[ni_i]
    W = Wbank[ni_i]

    X[:, 2] .= sex
    X[:, 3] .= age_std
    X[:, 4] = bmi_std

    W[:, 2] .= sex
    W[:, 3] .= age_std
    W[:, 4] = bp_std


    # generate random effects: γω = Lγω * z
    mul!(γω, Lγω, randn!(z))
    # generate y
    μy = X * βtrue + Z * γω[1:q]
    @views ysd = exp.(0.5 .* (W * τtrue .+ dot(γω[1:q], lγω) .+ γω[end]))
    y = ysd .* randn(ni_i) .+ μy
    # form a VarLmmObs instance
    obsvec[i] = WSVarLmmObs(y, X, Z, W)
end

# fit null model and create testing object
nm = WSVarLmmModel(obsvec)
WiSER.fit!(nm)
show(nm)
st = WSVarScoreTestInvariant(nm, 1, 1)

Ks = vGWAS.ecgf(st)

fakesnp = zeros(m)
tmp = Array{Float64}(undef, m)
tmp2 = Array{Float64}(undef, 3)

for j in 1:10
    global fakesnp
    wsample!(alleles, pAlleles_p, fakesnp)
    ps, zs = vGWAS.test!(st, fakesnp, fakesnp)
    bpval, tpval, jpval = ps
    bzval, tzval = zs

    #SPA for tau 
    bpspa, tpspa, jpspa = vGWAS.spa(fakesnp, st, ps, Ks; tmp_g=tmp, tmp_g2=tmp2)
    #fakesnp = normalize(fakesnp)
    #bpspa = vGWAS.spa(fakesnp, beta_pre_vec, bpval, 
    #        K0_beta, K1_beta, K2_beta; tmp_g = tmp, tmp_g2 = tmp2, r_var=r_var_beta) 
    #tpspa = vGWAS.spa(fakesnp, tau_pre_vec, tpval, 
    #        K0_tau, K1_tau, K2_tau; tmp_g = tmp, tmp_g2 = tmp2, r_var=r_var_tau) 
    #jpspa = vGWAS.spa(fakesnp, joint_pre_vec, jpval, 
    #        K0_joint, K1_joint, K2_joint; tmp_g = tmp, tmp_g2 = tmp2, r_var=r_var_joint) 
    
    println(gzio, "$bpval,$tpval,$jpval,$bpspa,$tpspa,$jpspa")
end

close(gzio)














