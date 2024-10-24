
##16/05/2024
##Sylvain Moinard, LECA

## Inference for simulated communities with different total initial amount of molecules
##____________________________________________________________________________________________________

commu1 = Float64.(Matrix(CSV.read("data/export_to_julia/sim100_m0tot_1.csv", DataFrame)))
commu2 = Float64.(Matrix(CSV.read("data/export_to_julia/sim100_m0tot_2.csv", DataFrame)))
commu3 = Float64.(Matrix(CSV.read("data/export_to_julia/sim100_m0tot_3.csv", DataFrame)))

K = 1.204e13

#pour chaque commu
#pour chaque sim
#une inference

ps_infer1 = zeros(size(commu1))
ps_infer2 = zeros(size(commu2))
ps_infer3 = zeros(size(commu3))

for sim in axes(commu1, 1)
  opt = optim_ps(reshape(commu1[sim,:], 1, size(commu1, 2)),
  Lambda_infer,
  ninfer = 1,
  randomTheta0 = true,
  previousTheta0 = false, 
  m0tot = 2.5e4,
  ncycles = 40,
  K = K,
  nsim = 190)
  ps_infer1[sim,:] .= opt.minimizer[1,:]
end

ps_infer1 = ps_infer1  ./ sum(ps_infer1, dims = 2)

mean(commu1, dims = 1)
mean(ps_infer1, dims = 1)

for sim in axes(commu2, 1)
  opt = optim_ps(reshape(commu2[sim,:], 1, size(commu2, 2)),
  Lambda_infer,
  ninfer = 1,
  randomTheta0 = true,
  previousTheta0 = false, 
  m0tot = 2.5e5,
  ncycles = 40,
  K = K,
  nsim = 190)
  ps_infer2[sim,:] .= opt.minimizer[1,:]
end

ps_infer2 = ps_infer2  ./ sum(ps_infer2, dims = 2)
mean(commu2, dims = 1)
mean(ps_infer2, dims = 1)

for sim in axes(commu3, 1)
  opt = optim_ps(reshape(commu3[sim,:], 1, size(commu3, 2)),
  Lambda_infer,
  ninfer = 1,
  randomTheta0 = true,
  previousTheta0 = false, 
  m0tot = 2.5e6,
  ncycles = 40,
  K = K,
  nsim = 190)
  ps_infer3[sim,:] .= opt.minimizer[1,:]
end

ps_infer3 = ps_infer3  ./ sum(ps_infer3, dims = 2)
mean(commu3, dims = 1)
mean(ps_infer3, dims = 1)

CSV.write("data/export_to_r/res_commu1.csv", Tables.table(ps_infer1 ./ sum(ps_infer1, dims = 2)))
CSV.write("data/export_to_r/res_commu2.csv", Tables.table(ps_infer2 ./ sum(ps_infer2, dims = 2)))
CSV.write("data/export_to_r/res_commu3.csv", Tables.table(ps_infer3 ./ sum(ps_infer3, dims = 2)))


##Below: Tests
#____________________________________________________________________________________________________

##Simulate another MU community
##____________________________________________________________________________________________________


simUm = simu_pcr(qty_init_U ./ 10. ;
  Lambda = Lambda_infer,
  K = K,
  n_reads = 1e5,
  ncycles = 40,
  nreplicate = 10)

simU = simu_pcr(qty_init_U ;
  Lambda = Lambda_infer,
  K = K,
  n_reads = 1e5,
  ncycles = 40,
  nreplicate = 10)

  simUp = simu_pcr(qty_init_U .* 10 ;
  Lambda = Lambda_infer,
  K = K,
  n_reads = 1e5,
  ncycles = 40,
  nreplicate = 10)

mean(simUm ./ sum(simUm, dims = 2), dims = 1)
std(simUm ./ sum(simUm, dims = 2), dims = 1)

mean(simU ./ sum(simU, dims = 2), dims = 1)
std(simU ./ sum(simU, dims = 2), dims = 1)

mean(simUp ./ sum(simUp, dims = 2), dims = 1)
std(simUp ./ sum(simUp, dims = 2), dims = 1)


x = mean(simU ./ sum(simU, dims = 2), dims = 1)
x[1]/x[13]

x = mean(simUm ./ sum(simUm, dims = 2), dims = 1)
x[1]/x[13]

x = mean(simUp ./ sum(simUp, dims = 2), dims = 1)
x[1]/x[13]

## Tests
##____________________________________________________________________________________________________

Random.seed!(1234)
Q = rand(190, 13*42)

J_efficiencies(optU_eff.minimizer, Q, reads_data = reads_U, n_reads = sum(reads_U, dims = 2))
J_efficiencies(fill(0.825, 12) , Q, reads_data = reads_U, n_reads = sum(reads_U, dims = 2))

a = simu_pcr_normalQ_eff(Lambda_infer, Q, K = K, n_reads = sum(reads_U, dims = 2))
b = simu_pcr_normalQ_eff(fill(0.825, 13), Q, K = K, n_reads = sum(reads_U, dims = 2))

mean(reads_U ./ sum(reads_U, dims = 2), dims = 1)
mean(a ./ sum(a, dims = 2), dims = 1)
mean(b ./ sum(b, dims = 2), dims = 1)

#J_ps([1., 1., 1.], rand(190, 1000), reads_data = reads_dataU, Lambda = Lambda, ncycles = 40)

println("___")
println(repl)
println(typeof(pdata))
println(typeof(psimu))
println(((psimu .- pdata).^2. ./ pdata))



##Serial inference for different values of lambda_max

Lambda_max = vcat((0.7:0.05:0.95), 0.999)

res_lambda = zeros(length(Lambda_max), 13)
scores = zeros(length(Lambda_max))
i = 1

nsim2 = 190
Q = rand(nsim2, 13*42)

for lambda_max in Lambda_max
    println(lambda_max)
    optU_eff2 = optim_efficiencies(reads_U, qty_init_U,
    ncycles = 40,
    dispersion = 1.,
    nsim = nsim2,
    lambda_max = lambda_max,
    quantiles = Q)
    res_lambda[i,:] .= vcat(lambda_max, vec(optU_eff2.minimizer))
    scores[i] = optU_eff2.minimum[1]
    i = i + 1
end

scores
res_lambda

# CSV.write("data/efficiencies_U.csv", Tables.table(res_lambda))


#mean(reads_U, dims = 1) ./ sum(mean(reads_U, dims = 1)) ./ (1 .+ optU_eff.minimizer).^26

sum(reads_U, dims = 1)



res_lambda[:,[9, 7, 8]]
Lambda

res_lambda ./ res_lambda[:,1]


optU_eff2 = optim_efficiencies(reads_U, qty_init_U,
  ncycles = 40,
  ninfer = 1,
  nsim = nsim,
  Kmult = 1e12,
  Kmin = 0.01,
  Kmax = 10)

  res = vcat(optU_eff.minimizer, optU_eff2.minimizer)

  median(res, dims = 1)


  ##Varying m0tot and K
##____________________________________________________________________________________________________

scoreT = zeros(10*10)
ps0T = zeros(10*10, 13)

scoreG = zeros(10*10)
ps0G = zeros(10*10, 13)

i = 0

for m0tot in (0.5:0.5:5)
  i = i+1
  optT_full = optim_ps(reads_T, Lambda_infer,
  ninfer = 10,
  randomTheta0 = true,
  previousTheta0 = false,
  m0tot = m0tot*1e5,
  ncycles = 40,
  K = K,
  nsim = 190)
  scoreT[i:(i+9)] = optT_full.minimum
  ps0T[i:(i+9),:] = optT_full.minimizer

  optG_full = optim_ps(reads_G, Lambda_infer,
  m0tot = m0tot*1e5,
  ncycles = 40,
  K = K,
  nsim = 190)
  scoreG[i] = optG_full.minimum[1]
  ps0G[i,:] = vec(optG_full.minimizer)
end

scoreT
ps0T = ps0T ./ sum(ps0T, dims = 2)
std(ps0T, dims = 1)
std(ps0T, dims = 1) ./mean(ps0T, dims = 1)

scoreG
ps0G = ps0G ./ sum(ps0G, dims = 2)
std(ps0G, dims = 1)
std(ps0G, dims = 1) ./mean(ps0G, dims = 1)
