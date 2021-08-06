#Code to do type I error simulation for vGWAS.jl
# after reworking SPA 

using vGWAS, WiSER, Distributions, Random, LinearAlgebra, StatsBase
using CodecZlib, Interpolations, Roots
import Optim: minimizer, optimize, LBFGS, NelderMead
BLAS.set_num_threads(1)

# Covariates input by runtypeIsims.jl
# maf = minor allele freq p 
# ni = # obs per person, a range
# m = # number of individuals
# effectsize = # effect size to evaluate
println("effectsize = $effectsize, maf = $maf, ni = $ni, m = $m")


nistring = string(minimum(ni)) * "to" * string(maximum(ni))

# set parameters

nreps = 50
pA = 1 - maf
paa = maf^2
pAa = 2 * maf *  pA
pAA = pA^2
alleles = [0.0, 1.0, 2.0]
pAlleles_p = [pAA, pAa, paa]

βtrue = [10.0; 5.0; 0.5; -0.3; fill(effectsize, 20)] 
τtrue = [0.25; 0.3; -0.15; 0.1; fill(effectsize, 20)] 
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
gzio = GzipCompressorStream(open("pvals_power_maf$(maf)_effectsize$(effectsize)_ni$(nistring)_m$(m)__.csv.gz", "w"))
#println(gzio, "betapval,taupval,jointpval,beta_z,tau_z,taupval_spa")
println(gzio, "betapval,taupval,jointpval,betapval_spa,taupval_spa,jointpval_spa")
seed = 190 + Int(maf * 1_000_000) + round(Int, effectsize * 10000) + maximum(ni) * 100 + m
Random.seed!(seed)

tmp = Array{Float64}(undef, m)
tmp2 = Array{Float64}(undef, 3)
pvs = Array{Float64}(undef, 20)

for j in 1:nreps 
    snps = wsample(alleles, pAlleles_p, (m, 20))
    snp_mat_norm = normalize(snps);
    for i in 1:m
        ni_i = rand(ni)
        sex = rand(0:1)
        age_std = randn() #time-invariant
        bmi_std = randn(ni_i) #time-varying
        bp_std = randn(ni_i) #time-varying

        X = Matrix{Float64}(undef, ni_i, p - 20)
        Z = Matrix{Float64}(undef, ni_i, q)
        W = Matrix{Float64}(undef, ni_i, l - 20)
        Xtrue = Matrix{Float64}(undef, ni_i, p)
        Wtrue = Matrix{Float64}(undef, ni_i, l)
        X[:, 1] .= 1
        Z[:, 1] .= 1
        W[:, 1] .= 1
        Xtrue[:, 1] .= 1
        Wtrue[:, 1] .= 1

        X[:, 2] .= sex
        X[:, 3] .= age_std
        X[:, 4] = bmi_std

        W[:, 2] .= sex
        W[:, 3] .= age_std
        W[:, 4] = bp_std

        Xtrue[:, 2] .= sex
        Xtrue[:, 3] .= age_std
        Xtrue[:, 4] = bmi_std

        Wtrue[:, 2] .= sex
        Wtrue[:, 3] .= age_std
        Wtrue[:, 4] = bp_std
        # add snps to W and X 
        for ij in 1:20
            Xtrue[:, 4 + ij] .= snps[i, ij]
            Wtrue[:, 4 + ij] .= snps[i, ij]
        end

        # generate random effects: γω = Lγω * z
        mul!(γω, Lγω, randn!(z))
        # generate y
        μy = Xtrue * βtrue + Z * γω[1:q]
        @views ysd = exp.(0.5 .* (Wtrue * τtrue .+ dot(γω[1:q], lγω) .+ γω[end]))
        y = ysd .* randn(ni_i) .+ μy
        # form a VarLmmObs instance
        obsvec[i] = WSVarLmmObs(y, X, Z, W)
    end
    # fit null model and create testing object
    nm = WSVarLmmModel(obsvec)
    WiSER.fit!(nm)
    show(nm)
    st = vGWAS.WSVarScoreTestInvariant(nm, 1, 1)
    
    Ks = vGWAS.ecgf(st)

    for snpind in 1:20
        @views ps, zs = vGWAS.test!(st, snps[:, snpind], snps[:, snpind])
        bpval, tpval, jpval = ps
        bzval, tzval = zs

	snp = view(snps, :, snpind)

	bpspa, tpspa, jpspa = vGWAS.spa(snp, st, ps, Ks; tmp_g=tmp, tmp_g2=tmp2)

    	#bpspa = vGWAS.spa(snp, beta_pre_vec, tpval, 
        #    K0_beta, K1_beta, K2_beta; tmp_g = tmp, tmp_g2 = tmp2, r_var=r_var_beta) 
    	#tpspa = vGWAS.spa(snp, tau_pre_vec, tpval, 
        #    K0_tau, K1_tau, K2_tau; tmp_g = tmp, tmp_g2 = tmp2, r_var=r_var_tau) 
    	#jpspa = vGWAS.spa(snp, joint_pre_vec, tpval, 
        #    K0_joint, K1_joint, K2_joint; tmp_g = tmp, tmp_g2 = tmp2, r_var=r_var_joint) 

        #print results 
    	println(gzio, "$bpval,$tpval,$jpval,$bpspa,$tpspa,$jpspa")
    end
end

close(gzio)














