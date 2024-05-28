
##Inference of initial proportions in metabarcoding data from Taqman data
##Sylvain Moinard, LECA
## 2023

#run metabar_bias_functions.jl before this file

##____________________________________________________________________________________________________
## Pipeline

#From the final proportions of each species and the PCR amplifiation rates inferred or measured by Taqman qPCR,
#we infer the initial proportions of each species with the Fixed Landscape Inference MethOd (flimo)

#See https://git.metabarcoding.org/lecasofts/flimo


##____________________________________________________________________________________________________
##Inference from Taqman data

## Import data

#Species :
# Carpinus_betulus
# Capsella_bursa-pastoris
# Fraxinus_excelsior

nspecies = 3

reads_dataU = Float64.(Matrix(CSV.read("data/dfSper01U_taq.csv", DataFrame)[:,2:(nspecies+1)]))
reads_dataT = Float64.(Matrix(CSV.read("data/dfSper01T_taq.csv", DataFrame)[:,2:(nspecies+1)]))
reads_dataG = Float64.(Matrix(CSV.read("data/dfSper01G_taq.csv", DataFrame)[:,2:(nspecies+1)]))

mean(reads_dataU ./ sum(reads_dataU, dims = 2), dims = 1)
mean(reads_dataT ./ sum(reads_dataT, dims = 2), dims = 1)
mean(reads_dataG ./ sum(reads_dataG, dims = 2), dims = 1)

## Define Lambda
Lambda = [mean([0.972, 0.947]), 0.989, 0.940]

#Lambda_hybrid = mean(Lambda)

#Define K
K = 7.57e12

## Optimize
#nsim proportional to nreplicates
nsim = 1900

#_____
optU = optim_ps(reads_dataU, Lambda,
  m0tot = 1e5,
  ncycles = 40,
  K = K,
  dispersion = 1.,
  nsim = nsim)

optU.minimizer ./ sum(optU.minimizer)
reshape(fill(1, nspecies) ./ nspecies, (1,nspecies))
mean(reads_dataU ./ sum(reads_dataU, dims = 2), dims = 1)

inferU = vec(optU.minimizer ./ sum(optU.minimizer))
vraiU = fill(1, nspecies) ./ nspecies
finU = vec(mean(reads_dataU ./ sum(reads_dataU, dims = 2), dims = 1))

#Aboslute RMSE:
sqrt(sum((inferU .- vraiU) .^ 2)/nspecies)
sqrt(sum((finU .- vraiU) .^ 2)/nspecies)

#Relative RMSE:
sqrt(sum(((inferU .- vraiU) ./ vraiU) .^ 2)/nspecies)
sqrt(sum(((finU .- vraiU)./ vraiU) .^ 2)/nspecies)

##Exponential model
prop_exp(reads_dataU, Lambda)

#_____
optT = optim_ps(reads_dataT, Lambda,
  m0tot = 1e5,
  ncycles = 40,
  K = K,
  dispersion = 1.,
  nsim = nsim)

optT.minimizer ./ sum(optT.minimizer)
reshape([41, 30, 67] ./ sum([41, 30, 67]), (1,nspecies))
mean(reads_dataT ./ sum(reads_dataT, dims = 2), dims = 1)

inferT = vec(optT.minimizer ./ sum(optT.minimizer))
vraiT = [20, 14, 33] ./ sum([20, 14, 33])
finT = vec(mean(reads_dataT ./ sum(reads_dataT, dims = 2), dims = 1))

#Absolute RMSE:
sqrt(sum((inferT .- vraiT) .^ 2)/nspecies)
sqrt(sum((finT .- vraiT) .^ 2)/nspecies)

#Relative RMSE:
sqrt(sum(((inferT .- vraiT) ./ vraiT) .^ 2)/nspecies)
sqrt(sum(((finT .- vraiT)./ vraiT) .^ 2)/nspecies)

##Exponential model
prop_exp(reads_dataT, Lambda)

#_____
optG = optim_ps(reads_dataG, Lambda,
  m0tot = 2.5e5,
  ncycles = 40,
  K = K,
  dispersion = 1.,
  nsim = nsim)

optG.minimizer ./ sum(optG.minimizer)
reshape([1, 4, 16] ./ 21, (1,nspecies))
mean(reads_dataG ./ sum(reads_dataG, dims = 2), dims = 1)

inferG = vec(optG.minimizer ./ sum(optG.minimizer))
vraiG = [1, 4, 16] ./ 21
finG = vec(mean(reads_dataG ./ sum(reads_dataG, dims = 2), dims = 1))

#Absolute RMSE:
sqrt(sum((inferG .- vraiG) .^ 2)/nspecies)
sqrt(sum((finG .- vraiG) .^ 2)/nspecies)

#Relative RMSE:
sqrt(sum(((inferG .- vraiG) ./ vraiG) .^ 2)/nspecies)
sqrt(sum(((finG .- vraiG)./ vraiG) .^ 2)/nspecies)

##Exponential model
prop_exp(reads_dataT, Lambda)


##____________________________________________________________________________________________________
#Using either CbeA or CbeB

LambdaA = [0.972, 0.989, 0.940]
LambdaB = [0.947, 0.989, 0.940]

#Change LambdaA to LambdaB to compare

#_____
optU = optim_ps(reads_dataU, LambdaA,
  m0tot = 1e5,
  ncycles = 40,
  K = K,
  dispersion = 1.,
  nsim = nsim)

optU.minimizer ./ sum(optU.minimizer)
reshape(fill(1, nspecies) ./ nspecies, (1,nspecies))
mean(reads_dataU ./ sum(reads_dataU, dims = 2), dims = 1)

inferU = vec(optU.minimizer ./ sum(optU.minimizer))
vraiU = fill(1, nspecies) ./ nspecies
finU = vec(mean(reads_dataU ./ sum(reads_dataU, dims = 2), dims = 1))

#Aboslute RMSE:
sqrt(sum((inferU .- vraiU) .^ 2)/nspecies)
sqrt(sum((finU .- vraiU) .^ 2)/nspecies)

#Relative RMSE:
sqrt(sum(((inferU .- vraiU) ./ vraiU) .^ 2)/nspecies)
sqrt(sum(((finU .- vraiU)./ vraiU) .^ 2)/nspecies)

#_____
optT = optim_ps(reads_dataT, LambdaA,
  m0tot = 1e5,
  ncycles = 40,
  K = K,
  dispersion = 1.,
  nsim = nsim)

optT.minimizer ./ sum(optT.minimizer)
reshape([41, 30, 67] ./ sum([41, 30, 67]), (1,nspecies))
mean(reads_dataT ./ sum(reads_dataT, dims = 2), dims = 1)

inferT = vec(optT.minimizer ./ sum(optT.minimizer))
vraiT = [41, 30, 67] ./ sum([41, 30, 67])
finT = vec(mean(reads_dataT ./ sum(reads_dataT, dims = 2), dims = 1))

#Absolute RMSE:
sqrt(sum((inferT .- vraiT) .^ 2)/nspecies)
sqrt(sum((finT .- vraiT) .^ 2)/nspecies)

#Relative RMSE:
sqrt(sum(((inferT .- vraiT) ./ vraiT) .^ 2)/nspecies)
sqrt(sum(((finT .- vraiT)./ vraiT) .^ 2)/nspecies)

#_____
optG = optim_ps(reads_dataG, LambdaA,
  m0tot = 1e5,
  ncycles = 40,
  K = K,
  dispersion = 1.,
  nsim = nsim)

optG.minimizer ./ sum(optG.minimizer)
reshape([12, 24, 190] ./ sum([12, 24, 190]), (1,nspecies))
mean(reads_dataG ./ sum(reads_dataG, dims = 2), dims = 1)

inferG = vec(optG.minimizer ./ sum(optG.minimizer))
vraiG = [12, 24, 190] ./ sum([12, 24, 190])
finG = vec(mean(reads_dataG ./ sum(reads_dataG, dims = 2), dims = 1))

#Absolute RMSE:
sqrt(sum((inferG .- vraiG) .^ 2)/nspecies)
sqrt(sum((finG .- vraiG) .^ 2)/nspecies)

#Relative RMSE:
sqrt(sum(((inferG .- vraiG) ./ vraiG) .^ 2)/nspecies)
sqrt(sum(((finG .- vraiG)./ vraiG) .^ 2)/nspecies)


##____________________________________________________________________________________________________
#Check type stability
@code_warntype optim_ps(reads_dataU, Lambda,
    m0tot = 1e5,
    ncycles = 40,
    K = 1e11,
    dispersion = 1.,
    nsim = 1000)

##____________________________________________________________________________________________________
# Mechanistic model

optU_mec = optim_ps(reads_dataU, Lambda,
  m0tot = 1e5,
  ncycles = 40,
  K = K,
  dispersion = 1.,
  nsim = nsim,
  model = "",
  c = 0.9,
  e = 1.1)

optU_mec.minimizer ./ sum(optU_mec.minimizer)
optU.minimizer ./ sum(optU.minimizer)
reshape(fill(1, nspecies) ./ nspecies, (1,nspecies))
mean(reads_dataU ./ sum(reads_dataU, dims = 2), dims = 1)

optT_mec = optim_ps(reads_dataT, Lambda,
  m0tot = 1e5,
  ncycles = 40,
  K = K,
  dispersion = 1.,
  nsim = nsim,
  model = "",
  c = 0.9,
  e = 1.1)

optT_mec.minimizer ./ sum(optT_mec.minimizer)
optT.minimizer ./ sum(optT.minimizer)
reshape([41, 30, 67] ./ sum([41, 30, 67]), (1,nspecies))
mean(reads_dataT ./ sum(reads_dataT, dims = 2), dims = 1)


optG_mec = optim_ps(reads_dataG, Lambda,
  m0tot = 1e5,
  ncycles = 40,
  K = K,
  dispersion = 1.,
  nsim = nsim,
  model = "",
  c = 0.9,
  e = 1.1)

optG_mec.minimizer ./ sum(optG_mec.minimizer)
optG.minimizer ./ sum(optG.minimizer)
reshape([12, 24, 190] ./ sum([12, 24, 190]), (1,nspecies))
mean(reads_dataG ./ sum(reads_dataG, dims = 2), dims = 1)

##____________________________________________________________________________________________________
#Adding an "hybrid" species accounting for every other species
Lambda_hybrid = mean(Lambda)
Lambda_complete = vcat(Lambda, Lambda_hybrid)

reads_dataU_complete = Float64.(Matrix(CSV.read("data/dfSper01U_taq_complete.csv", DataFrame)[:,2:5]))
reads_dataT_complete = Float64.(Matrix(CSV.read("data/dfSper01T_taq_complete.csv", DataFrame)[:,2:5]))
reads_dataG_complete = Float64.(Matrix(CSV.read("data/dfSper01G_taq_complete.csv", DataFrame)[:,2:5]))


#_____
optU_complete = optim_ps(reads_dataU_complete, Lambda_complete,
  m0tot = 1e5,
  ncycles = 40,
  K = K,
  dispersion = 1.,
  nsim = nsim)

optU_complete.minimizer[:,1:3] ./ sum(optU_complete.minimizer[:,1:3])
reshape(fill(1, 3) ./ 3, (1,3))
mean(reads_dataU ./ sum(reads_dataU, dims = 2), dims = 1)

optU_complete.minimizer ./ sum(optU_complete.minimizer)
mean(reads_dataU_complete ./ sum(reads_dataU_complete, dims = 2), dims = 1)



inferU = vec(optU_complete.minimizer[:,1:3] ./ sum(optU_complete.minimizer[:,1:3]))
vraiU = fill(1, 3) ./ 3
finU = vec(mean(reads_dataU ./ sum(reads_dataU, dims = 2), dims = 1))

#Aboslute RMSE:
sqrt(sum((inferU .- vraiU) .^ 2)/3)
sqrt(sum((finU .- vraiU) .^ 2)/3)

#Relative RMSE:
sqrt(sum(((inferU .- vraiU) ./ vraiU) .^ 2)/3)
sqrt(sum(((finU .- vraiU)./ vraiU) .^ 2)/3)

#_____
optT_complete = optim_ps(reads_dataT_complete, Lambda_complete,
  m0tot = 2e5,
  ncycles = 20,
  K = K/7,
  dispersion = 1.,
  nsim = nsim)

optT_complete.minimizer[:,1:3] ./ sum(optT_complete.minimizer[:,1:3])
reshape([41, 30, 67] ./ sum([41, 30, 67]), (1,3))
mean(reads_dataT ./ sum(reads_dataT, dims = 2), dims = 1)

mean(reads_dataT_complete ./ sum(reads_dataT_complete, dims = 2), dims = 1)

inferT = vec(optT_complete.minimizer[:,1:3] ./ sum(optT_complete.minimizer[:,1:3]))
vraiT = [41, 30, 67] ./ sum([41, 30, 67])
finT = vec(mean(reads_dataT ./ sum(reads_dataT, dims = 2), dims = 1))

#Absolute RMSE:
sqrt(sum((inferT .- vraiT) .^ 2)/3)
sqrt(sum((finT .- vraiT) .^ 2)/3)

#Relative RMSE:
sqrt(sum(((inferT .- vraiT) ./ vraiT) .^ 2)/3)
sqrt(sum(((finT .- vraiT)./ vraiT) .^ 2)/3)

#_____
optG_complete = optim_ps(reads_dataG_complete, Lambda_complete,
  m0tot = 1e5,
  ncycles = 40,
  K = K,
  dispersion = 1.,
  nsim = nsim)

optG_complete.minimizer[:,1:3] ./ sum(optG_complete.minimizer[:,1:3])
reshape([12, 24, 190] ./ sum([12, 24, 190]), (1,3))
mean(reads_dataG ./ sum(reads_dataG, dims = 2), dims = 1)

mean(reads_dataG_complete ./ sum(reads_dataG_complete, dims = 2), dims = 1)

inferG = vec(optG_complete.minimizer[:,1:3] ./ sum(optG_complete.minimizer[:,1:3]))
vraiG = [12, 24, 190] ./ sum([12, 24, 190])
finG = vec(mean(reads_dataG ./ sum(reads_dataG, dims = 2), dims = 1))

#Absolute RMSE:
sqrt(sum((inferG .- vraiG) .^ 2)/3)
sqrt(sum((finG .- vraiG) .^ 2)/3)

#Relative RMSE:
sqrt(sum(((inferG .- vraiG) ./ vraiG) .^ 2)/3)
sqrt(sum(((finG .- vraiG)./ vraiG) .^ 2)/3)

x = finG ./ (1 .+ Lambda).^3
x./ sum(x)



##____________________________________________________________________________________________________
#Tests

optU.minimum

J_ps([1., 1., 1.], rand(190, 1000), reads_data = reads_dataU, Lambda = Lambda, ncycles = 40)
J_ps(optU.minimizer, rand(190, 1000), reads_data = reads_dataU, Lambda = Lambda, ncycles = 40)


sim = simu_pcr_normalQ(1e5 .* [1., 1., 1.], rand(1, 1000) ; Lambda = Lambda, n_reads = 5e4)
sim ./ sum(sim)

simu_pcr(1e5 .* [1., 1., 1.], Lambda = Lambda, n_reads = 5e4)

puc = mean(reads_dataU_complete ./ sum(reads_dataU_complete, dims = 2), dims = 1)
pu = mean(reads_dataU ./ sum(reads_dataU, dims = 2), dims = 1)

m0tot = 2.47e5

sim = simu_pcr_normalQ(m0tot .* puc, rand(1, 1000), K=K, Lambda = Lambda_hybrid)
sim = simu_pcr_normalQ(m0tot .* pu, rand(1, 1000), K=K, Lambda = Lambda)

sim = simu_pcr_normalQ(m0tot .* pu, rand(1, 1000), K=K, Lambda = Lambda, model = "", e = 1.1, c = 0.9) #mechanistic

sim ./ sum(sim)