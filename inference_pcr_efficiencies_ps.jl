
##Inference of PCR efficiencies in metabarcoding data
##Sylvain Moinard, LECA
## 2023

#From the final proportions of each species and the PCR amplifiation rates inferred or measured by Taqman qPCR,
#we infer the initial proportions of each species with the Fixed Landscape Inference MethOd (flimo)

#See https://git.metabarcoding.org/lecasofts/flimo

## Use functions in file metabar_bias_functions.jl

##____________________________________________________________________________________________________
# Setup: load data

# U community : reads

reads_U = Float64.(Matrix(CSV.read("data/reads_U.csv", DataFrame)))

# Initial number of molecules, assayed by ddPCR

qty_init_U = fill(19084., size(reads_U, 2))


##____________________________________________________________________________________________________
# Efficiencies inference

nreplicate = size(reads_U, 1)
n_reads = sum(reads_U, dims = 2)

## Optimize
#nsim proportional to nreplicates
nsim = 190

Random.seed!(21011996)
optU_eff = optim_efficiencies(reads_U, qty_init_U,
  ncycles = 40,
  ninfer = 101,
  nsim = nsim,
  Kmult = 1e13,
  Kmin_log = -3, #log10
  Kmax_log = 3) #log10

optU_eff.minimizer
optU_eff.minimum
optU_eff.time
optU_eff.initial


CSV.write("data/res_infer_Lambda_K.csv", Tables.table(optU_eff.minimizer))
#CSV.write("data/res_infer_J.csv", Tables.table(optU_eff.minimum))

#Stats
#K
mean(optU_eff.initial[:,13])
median(optU_eff.minimizer[:,13])
std(optU_eff.minimizer[:,13])
std(10 .^optU_eff.minimizer[:,13].*1e13)
10^median(optU_eff.minimizer[:,13])*1e13

#Lambda and K
mean(optU_eff.minimizer, dims = 1)
median(optU_eff.minimizer, dims = 1)
std(optU_eff.minimizer, dims = 1)

index = findall(x->x==median(optU_eff.minimizer, dims = 1)[1,13], optU_eff.minimizer[:,13])[1]

Lambda_infer = vcat(1., vec(optU_eff.minimizer[index,1:12]))
CSV.write("data/efficiencies_U.csv", Tables.table(Lambda_infer))
K = 10^optU_eff.minimizer[index,13]*1e13

simU = simu_pcr_normalQ_eff(vcat(1., vec(optU_eff.minimizer)),
    rand(190, 13*62),
    m0 = qty_init_U,
    K = K,
    n_reads = vec(repeat(sum(reads_U, dims = 2), 10)),
    ncycles = 40)

mean(simU, dims = 1) ./ sum(mean(simU, dims = 1))
mean(reads_U, dims = 1) ./ sum(mean(reads_U, dims = 1))

Lambda #Cbe, Cbp, Fex, file inference_metabarcoding.jl
optU_eff.minimizer[[8, 6, 7]] #Cbe, Cbp, Fe

optU_eff.minimizer[[8, 6, 7]] ./ Lambda


optU_eff.minimizer[[8]] / 0.972 #CbeA
optU_eff.minimizer[[8]] / 0.947 #CbeB


## Proportions inference for communities T and G with inferred Lambdas
##____________________________________________________________________________________________________

psU = mean(reads_U ./ sum(reads_U, dims = 2), dims = 1)

neffU = (mean(log.(psU))+log(K)-mean(log.(qty_init_U)))/mean(log.(1 .+ Lambda_infer))

reads_T = Float64.(Matrix(CSV.read("data/reads_T.csv", DataFrame)))
reads_G = Float64.(Matrix(CSV.read("data/reads_G.csv", DataFrame)))

psT = mean(reads_T ./ sum(reads_T, dims = 2), dims = 1)
psG = mean(reads_G ./ sum(reads_G, dims = 2), dims = 1)

m0totU = K .* sum(vec(psU) ./ (1 .+ Lambda_infer).^neffU)
m0totT = K .* sum(vec(psT) ./ (1 .+ Lambda_infer).^neffU)
m0totG = K .* sum(vec(psG) ./ (1 .+ Lambda_infer).^neffU)

optU_full = optim_Qmetabar(reads_U, Lambda_infer,
  ninfer = 10,
  randomTheta0 = true,
  previousTheta0 = false,
  m0tot = m0totU,
  ncycles = 40,
  K = K,
  nsim = 190)

optU_full.minimizer ./ sum(optU_full.minimizer, dims = 2)
mean(optU_full.minimizer ./ sum(optU_full.minimizer, dims = 2), dims = 1)
std(optU_full.minimizer ./ sum(optU_full.minimizer, dims = 2), dims = 1)

mean(reads_U ./ sum(reads_U, dims = 2), dims = 1)


optT_full = optim_Qmetabar(reads_T, Lambda_infer,
  ninfer = 10,
  randomTheta0 = true,
  previousTheta0 = false, 
  m0tot = m0totT,
  ncycles = 40,
  K = K,
  nsim = 190)

optT_full.minimizer ./ sum(optT_full.minimizer, dims = 2)
mean(optT_full.minimizer ./ sum(optT_full.minimizer, dims = 2), dims = 1)
std(optT_full.minimizer ./ sum(optT_full.minimizer, dims = 2), dims = 1)

std(optT_full.minimizer ./ sum(optT_full.minimizer, dims = 2), dims = 1) ./ mean(optT_full.minimizer ./ sum(optT_full.minimizer, dims = 2), dims = 1)

mean(reads_T ./ sum(reads_T, dims = 2), dims = 1)


optG_full = optim_Qmetabar(reads_G, Lambda_infer,
  ninfer = 10,
  randomTheta0 = true,
  previousTheta0 = false, 
  m0tot = m0totG,
  ncycles = 40,
  K = K,
  nsim = 190)

optG_full.minimizer ./ sum(optG_full.minimizer, dims = 2)
mean(optG_full.minimizer ./ sum(optG_full.minimizer, dims = 2), dims = 1)
std(optG_full.minimizer ./ sum(optG_full.minimizer, dims = 2), dims = 1)

std(optG_full.minimizer ./ sum(optG_full.minimizer, dims = 2), dims = 1) ./ mean(optG_full.minimizer ./ sum(optG_full.minimizer, dims = 2), dims = 1)

vec(optG_full.minimizer ./ sum(optG_full.minimizer))
mean(reads_G ./ sum(reads_G, dims = 2), dims = 1)

CSV.write("data/prop_inferU.csv", Tables.table(mean(optU_full.minimizer ./ sum(optU_full.minimizer, dims = 2), dims = 1)))
CSV.write("data/prop_inferT.csv", Tables.table(mean(optT_full.minimizer ./ sum(optT_full.minimizer, dims = 2), dims = 1)))
CSV.write("data/prop_inferG.csv", Tables.table(mean(optG_full.minimizer ./ sum(optG_full.minimizer, dims = 2), dims = 1)))


## Inference for simulated communities with different total initial amount of molecules
##____________________________________________________________________________________________________

commu1 = Float64.(Matrix(CSV.read("data/sim100_m0tot_1.csv", DataFrame)))
commu2 = Float64.(Matrix(CSV.read("data/sim100_m0tot_2.csv", DataFrame)))
commu3 = Float64.(Matrix(CSV.read("data/sim100_m0tot_3.csv", DataFrame)))

#pour chaque commu
#pour chaque sim
#une inference

ps_infer1 = zeros(size(commu1))
ps_infer2 = zeros(size(commu2))
ps_infer3 = zeros(size(commu3))

for sim in 1:size(commu1, 1)
  opt = optim_Qmetabar(reshape(commu1[sim,:], 1, size(commu1, 2)),
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

for sim in 1:size(commu2, 1)
  opt = optim_Qmetabar(reshape(commu2[sim,:], 1, size(commu2, 2)),
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

for sim in 1:size(commu3, 1)
  opt = optim_Qmetabar(reshape(commu3[sim,:], 1, size(commu3, 2)),
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

CSV.write("data/res_commu1.csv", Tables.table(ps_infer1 ./ sum(ps_infer1, dims = 2)))
CSV.write("data/res_commu2.csv", Tables.table(ps_infer2 ./ sum(ps_infer2, dims = 2)))
CSV.write("data/res_commu3.csv", Tables.table(ps_infer3 ./ sum(ps_infer3, dims = 2)))


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

#J_Qmetabar([1., 1., 1.], rand(190, 1000), reads_data = reads_dataU, Lambda = Lambda, ncycles = 40)

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
  optT_full = optim_Qmetabar(reads_T, Lambda_infer,
  ninfer = 10,
  randomTheta0 = true,
  previousTheta0 = false,
  m0tot = m0tot*1e5,
  ncycles = 40,
  K = K,
  nsim = 190)
  scoreT[i:(i+9)] = optT_full.minimum
  ps0T[i:(i+9),:] = optT_full.minimizer

  optG_full = optim_Qmetabar(reads_G, Lambda_infer,
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
