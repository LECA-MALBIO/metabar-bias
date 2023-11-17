

##Inference of initial proportions in metabarcoding data
##Sylvain Moinard, LECA
## 2023

##____________________________________________________________________________________________________
## Pipeline

#From the final proportions of each species and the PCR amplifiation rates measured by Taqman qPCR,
#we infer the initial proportions of each species with the Fixed Landscape Inference MethOd (flimo)

#See https://git.metabarcoding.org/lecasofts/flimo


#TO DO 14/11/2023 : implement quantile multinomial (with cumsum)

##____________________________________________________________________________________________________
## Setup

using CSV
using DataFrames
using Distributions
using ForwardDiff
using Optim
using Random

#using Pkg
#Pkg.add(url ="https://git.metabarcoding.org/lecasofts/flimo/jflimo.git")
using Jflimo


##____________________________________________________________________________________________________
## Functions


#Function used in the mechanistic saturation model
function gamma_factor(x ; c = 1.)
  2/(1+sqrt(1-4*c^2*x*(1-x)))
end

#The simulators are adapted to the flimo method, i.e.
#each random draw (rand(...)) is replaced by a call to the appropriate quantile function (quantile(...))

#The normal approximation of the model makes it possible to obtain a continuous distribution, and therefore to use
#an automatic differentiation algorithm to minimise the objective function


## PCR Simulation

#The implemented model is the logistic model***
#
# Lambda is an array containing the amplification rates of each species (0 <= Lambda <= 1)
# K is the charge capacity
# n_reads is the number of reads after sequencing
# ncycles is the number of PCR cycles
# nreplicate is the number of replicates to simulate
# dispersion is the initial inter-replicate dispersion

#Theta: array containing the initial quantities
#quantiles: 2D array. Nrow = Number of simulations; Ncol = nspecies*(ncycles+2)

#Exact simulation:
function simu_pcr(Theta ;
  Lambda = fill(one(eltype(Theta)), length(Theta)),
  K = 1e11,
  n_reads = 1e5,
  ncycles = 40,
  nreplicate = 1,
  dispersion = 1.,
  model = "logistic",
  c = 1.,
  e = 1.)

#Average initial number of molecules of each species:
  m0 = Theta
  nspecies = length(m0) #Number of species
  kinetics = Array{Int64,2}(undef, nreplicate, nspecies) #initialize kinetics
  lam = Vector{Float64}(undef, nreplicate) #effective amplification rates
  crea = Vector{Int64}(undef, nreplicate) #new molecules create at each cycle for one species
  crea_tot = Vector{Int64}(undef, nreplicate) #total number of molecules created at each cycle
  kin_tot = zeros(eltype(kinetics), nreplicate) #total number of molecules created

#Initialisation
  @inbounds for species in 1:nspecies
    if dispersion <= 1.
      kinetics[:,species] .= rand.(Poisson(m0[species]))
    else
      if m0[species] == zero(eltype(m0))
        kinetics[:,species] .= zero(eltype(m0))
      else
        kinetics[:,species] .= rand.(NegativeBinomial(m0[species]/(dispersion^2-1),1/dispersion^2))
      end
    end
  end

#Amplification cycle by cycle
  @inbounds for cyc in 1:ncycles #for each cycle
    crea_tot[:] .= zero(eltype(m0))
    @inbounds for species in 1:nspecies #for each species
#effective rate
      if model == "logistic"
#Logistic model
        lam[:] .= Lambda[species].*(one(eltype(lam)) .- kin_tot ./ K)
      else
#mechanistic model
        lam[:] .= Lambda[species].* (one(eltype(lam)) .- kin_tot ./ K) .* c .* gamma_factor.(kin_tot ./ K, c = c) .* (e .- kin_tot ./ K)
        lam[lam .>= one(eltype(lam))] .= one(eltype(lam)) #happens if e.c > 1
      end
      lam[lam .<= zero(eltype(lam))] .= zero(eltype(lam))
#Amplification
      crea[:] .= rand.(Binomial.(view(kinetics,:,species), lam))
      broadcast!(+, crea_tot, crea_tot, crea) #update molecules created at this cycle
      kinetics[:,species] .= view(kinetics, :, species).+crea
    end
    broadcast!(+, kin_tot, kin_tot, crea_tot) #update molecules created
  end

#Sequencing
  if length(n_reads) == 1
    n_reads = fill(n_reads, nreplicate)
  end
  kinetics_tot = vec(sum(kinetics, dims = 2))
  for i in 1:nreplicate
    kinetics[i,:] .= vec(rand(Multinomial(Int64(n_reads[i]), vec(kinetics[i,:] ./ kinetics_tot[i]))))
  end

  return kinetics
end

function simu_pcrQ(Theta,
                   quantiles ;
                   Lambda = fill(one(eltype(Theta)), length(Theta)),
                   K = 1e11,
                   n_reads = 1e5,
                   ncycles = 40,
                   nreplicate = size(quantiles, 1),
                   dispersion = 1.,
                   model = "logistic",
                   c = 1.,
                   e = 1.)

                   if length(n_reads) == 1
                    n_reads = fill(1e5, nreplicate)
                  end
              
                  if length(n_reads) < nreplicate
                    n_reads = vec(repeat(n_reads, Int64(nreplicate/length(n_reads)))) #multiple needed
                  end

  #Average initial number of molecules of each species:
  m0 = Theta
  nspecies = length(m0) #Number of species
  kinetics = Array{Int64,2}(undef, nreplicate, nspecies) #initialize kinetics
  lam = Vector{Float64}(undef, nreplicate) #effective amplification rates
  crea = Vector{Int64}(undef, nreplicate) #new molecules create at each cycle for one species
  crea_tot = Vector{Int64}(undef, nreplicate) #total number of molecules created at each cycle
  kin_tot = zeros(eltype(kinetics), nreplicate) #total number of molecules created

  #Initialisation
  @inbounds for species in 1:nspecies
    if dispersion <= 1.
      kinetics[:,species] .= quantile.(Poisson(m0[species]), view(quantiles,1:nreplicate,species))
    else
      if m0[species] == zero(eltype(m0))
        kinetics[:,species] .= zero(eltype(m0))
      else
        kinetics[:,species] .= quantile.(NegativeBinomial(m0[species]/(dispersion^2-1),1/dispersion^2),
        view(quantiles,1:nrepli,species))
      end
    end
  end

  #Amplification cycle by cycle
  @inbounds for cyc in 1:ncycles #for each cycle
    crea_tot[:] .= zero(eltype(m0))
    @inbounds for species in 1:nspecies #for each species
      #effective rate
      if model == "logistic"
        #Logistic model
        lam[:] .= Lambda[species].*(one(eltype(lam)) .- kin_tot ./ K)
      else
        #mechanistic model
        lam[:] .= Lambda[species].* (one(eltype(lam)) .- kin_tot ./ K) .* c .* gamma_factor.(kin_tot ./ K, c = c) .* (e .- kin_tot ./ K)
        lam[lam .>= one(eltype(lam))] .= one(eltype(lam)) #happens if e.c > 1
      end
      lam[lam .<= zero(eltype(lam))] .= zero(eltype(lam))
      #Amplification
      crea[:] .= quantile.(Binomial.(view(kinetics,:,species), lam),
                          view(quantiles,1:nreplicate,nspecies*cyc+species))
      broadcast!(+, crea_tot, crea_tot, crea) #update molecules created at this cycle
      kinetics[:,species] .= view(kinetics, :, species).+crea
    end
    broadcast!(+, kin_tot, kin_tot, crea_tot) #update molecules created
  end

  #Sequencing -> no quantile function available for multinomial
  for species in 1:nspecies
    kinetics[:,species] .= quantile.(Binomial.(Int64.(n_reads),
                                      view(kinetics,:,species) ./ sum(kinetics, dims = 2)),
    view(quantiles,1:nreplicate,(nspecies*(ncycles+1)+species)))
  end
  return kinetics
end



#Normally approximated simulation:
function simu_pcr_normalQ(Theta,
    quantiles ;
    Lambda = fill(one(eltype(Theta)), length(Theta)),
    K = 1e11,
    nreplicate = size(quantiles, 1),
    n_reads = 1e5,
    ncycles = 40,
    dispersion = 1.,
    model = "logistic",
    c = 1.,
    e = 1.)

    if length(n_reads) == 1
      n_reads = fill(1e5, nreplicate)
    end

    if length(n_reads) < nreplicate
      n_reads = vec(repeat(n_reads, Int64(nreplicate/length(n_reads)))) #multiple needed
    end

    #Average initial number of molecules of each species:
    m0 = Theta
    nspecies = length(m0) #Number of species
    kinetics = Array{eltype(Theta),2}(undef, nreplicate, nspecies) #initialize kinetics
    lam = Vector{eltype(Theta)}(undef, nreplicate) #effective amplification rates
    crea = Vector{eltype(Theta)}(undef, nreplicate) #new molecules create at each cycle for one species
    crea_tot = Vector{eltype(Theta)}(undef, nreplicate) #total number of molecules created at each cycle
    kin_tot = zeros(eltype(Theta), nreplicate) #total number of molecules created

    #Initialisation
    @inbounds for species in 1:nspecies
        kinetics[:,species] .= quantile.(Normal(m0[species],
            dispersion*sqrt(m0[species])), view(quantiles,1:nreplicate,species))
    end
    kinetics[kinetics .< zero(eltype(m0))] .= zero(eltype(m0))
    #pc = true
    #Amplification cycle by cycle
    @inbounds for cyc in 1:ncycles #for each cycle
      #if (mean(kin_tot) >= mean(K)/10) & pc
      #  pc = false
      #end
      crea_tot[:] .= zero(eltype(m0))
      @inbounds for species in 1:nspecies #for each species
        #effective rate
        if model == "logistic"
          #Logistic model
          lam[:] .= Lambda[species].*(one(eltype(lam)) .- kin_tot ./ K)
        else
          #mechanistic model
          lam[:] .= Lambda[species].* (one(eltype(lam)) .- kin_tot ./ K) .* c .* gamma_factor.(kin_tot ./ K, c = c) .* (e .- kin_tot ./ K)
          lam[lam .>= one(eltype(lam))] .= one(eltype(lam)) #happens if e.c > 1
        end
        lam[lam .<= zero(eltype(lam))] .= zero(eltype(lam))
        #Amplification
        crea[:] .= quantile.(Normal.(view(kinetics,:,species) .* lam,
                            sqrt.(view(kinetics,:,species) .* lam .* (one(eltype(lam)) .- lam))),
                            view(quantiles,1:nreplicate,nspecies*cyc+species))
        crea[crea .< zero(eltype(crea))] .= zero(eltype(crea))
        broadcast!(+, crea_tot, crea_tot, crea) #update molecules created at this cycle
        kinetics[:,species] .= view(kinetics, :, species).+crea
      end
      broadcast!(+, kin_tot, kin_tot, crea_tot) #update molecules created
    end

  #Sequencing -> simulate marginal distributions
  kinetics_tot = sum(kinetics, dims = 2)

  for species in 1:nspecies
    #reuse crea array for pre-allocation
    crea[:] .= vec(view(kinetics,:, species) ./ kinetics_tot)
    kinetics[:,species] .= quantile.(Normal.(n_reads .* crea,
                                             sqrt.(n_reads .* crea .* (one(eltype(crea)) .- crea))),
                                    view(quantiles,1:nreplicate,(nspecies*(ncycles+1)+species)))
  end

  kinetics[kinetics .< zero(eltype(m0))] .= zero(eltype(m0))
  
  return kinetics .* n_reads ./ sum(kinetics, dims = 2) #***uncorrect conditionning
end



##Objective function
#To minimize to find optimal starting conditions

#reads_data: 2D array: Nrow : nreplicates Ncol : reads for each species


function J_Qmetabar(Theta,
    quantiles ;
    reads_data = rand(1,1)::Array{Float64, 2},
    Lambda = ones(eltype(Theta), length(Theta)),
    K = 1e11,
    n_reads = 1e5,
    ncycles = 40,
    dispersion = 1.,
    m0tot = 1e5,
    model = "logistic",
    c = 1.,
    e = 1.)

  nsim = size(quantiles, 1) #number of simulations to perform by flimo

  nreplicate = size(reads_data, 1)

  if length(n_reads) == 1
    n_reads = fill(1e5, nreplicate)
  end

  #nspecies = size(reads_data, 2)

  m0 = Theta .* m0tot #initial quantities at right scale

  n_reads_repli = Vector{eltype(m0)}(undef, nsim) #Number of reads for each simulation
  pdata = vec(mean(reads_data ./ sum(reads_data, dims = 2), dims = 1)) #average species proportions in data
  
  @inbounds for repl in 1:Int64(floor(nsim/nreplicate))
      #more simulations than replicates needed in this implementation
      n_reads_repli[(1+(repl-1)*nreplicate):(repl*nreplicate)] .= vec(n_reads)
    end
    if Int64(floor(nsim/nreplicate))*nreplicate < nsim
      n_reads_repli[(1+Int64(floor(nsim/nreplicate))*nreplicate):nsim] .= view(n_reads_repli,1:nsim-Int64(floor(nsim/nreplicate))*nreplicate)
    end
    Ksim = K

    kinetics = simu_pcr_normalQ(m0,
                                   quantiles ;
                                   Lambda = Lambda,
                                   K = Ksim,
                                   n_reads = n_reads_repli,
                                   ncycles = ncycles,
                                   nreplicate = nsim,
                                   dispersion = dispersion,
                                   model = model,
                                   c = c,
                                   e = e)
    
    @inbounds for sim in 1:nsim #Number of molecules for each simulation
      cin_tot = sum(view(kinetics, sim, :))
      kinetics[sim,:] .= view(kinetics, sim,:)./cin_tot #Relative abundance
    end
    
    psimu = vec(mean(kinetics, dims = 1))
    psimu[:] .= psimu ./ sum(psimu)
    
    score = sum(((psimu .- pdata).^2. ./ pdata))

  return score
end

##____________________________________________________________________________________________________
##Optim

#quantiles : Nrow = number replicates
#Ncol = nspecies * (ncycles+2)
#first column : initial randomness, then each cycle then sequencing

function optim_Qmetabar(reads_data::Array{Float64, 2}, Lambda::Array{Float64, 1};
  m0tot = 1e5::Float64,
  ncycles = 40::Int64,
  K = 1e11::Float64,
  nsim = 1000::Int64,
  Theta0 = nothing,
  ninfer = 1,
  randomTheta0 = false,
  previousTheta0 = true,
  maxit = 200::Int64,
  gtol = 1e-6::Float64,
  xtol=1e-32::Float64,
  dispersion = 1.::Float64,
  model = "logistic",
  c = 1.,
  e = 1.)

  #reads_data: 2D array. Rows = replicates, Col = species

  nreplicate = size(reads_data, 1)
  nspecies = size(reads_data, 2)

  n_reads_repli = sum(reads_data, dims = 2)

  pdata = vec(mean(reads_data ./ sum(reads_data, dims = 2), dims = 1)) #average species proportions in data

  if (isnothing(Theta0))
    Theta0 = pdata
  end
  obj = (Theta, quantiles) -> J_Qmetabar(Theta, quantiles,
                                        reads_data = reads_data,
                                        Lambda = Lambda,
                                        K = K,
                                        n_reads = n_reads_repli,
                                        ncycles = ncycles,
                                        dispersion = dispersion,
                                        m0tot = m0tot,
                                        model = model,
                                        c = c,
                                        e = e)

  opt = Jflimo.flimoptim(nspecies*(ncycles+2),
                     obj = obj,
                     nsim = nsim,
                     ninfer = ninfer,
                     lower = pdata ./ 5.,
                     upper = pdata .* 5.,
                     AD = true,
                     Theta0 = Theta0,
                     randomTheta0 = randomTheta0,
                     previousTheta0 = previousTheta0,
                     gtol = gtol,
                     xtol = xtol,
                     maxit = maxit,
                     show_trace = true,
                     show_every = 10)
    return opt
end


function prop_exp(reads_data, Lambda, low_c = 15, high_c = 30)
  #proportions with the exponential model for ncycles between low_c and high_c
  prop = zeros(length(low_c:high_c), size(reads_data, 2))
  for s in 1:axes(reads_data, 2)
    prop[:,s] = mean(reads_data[:,s]) ./ (1+Lambda[s]) .^(low_c:high_c) .* K ./ mean(sum(reads_data, dims = 2))
  end
  prop ./ sum(prop, dims = 2)
end


##____________________________________________________________________________________________________
##Inference

## Import data

#Species :
# Carpinus_betulus
# Capsella_bursa-pastoris
# Fraxinus_excelsior

#Change 4 by nspecies+1

reads_dataU = Float64.(Matrix(CSV.read("data/dfSper01U_taq.csv", DataFrame)[:,2:4]))
reads_dataT = Float64.(Matrix(CSV.read("data/dfSper01T_taq.csv", DataFrame)[:,2:4]))
reads_dataG = Float64.(Matrix(CSV.read("data/dfSper01G_taq.csv", DataFrame)[:,2:4]))

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
optU = optim_Qmetabar(reads_dataU, Lambda,
  m0tot = 1e5,
  ncycles = 40,
  K = K,
  dispersion = 1.,
  nsim = nsim)

optU.minimizer ./ sum(optU.minimizer)
reshape(fill(1, 3) ./ 3, (1,3))
mean(reads_dataU ./ sum(reads_dataU, dims = 2), dims = 1)

inferU = vec(optU.minimizer ./ sum(optU.minimizer))
vraiU = fill(1, 3) ./ 3
finU = vec(mean(reads_dataU ./ sum(reads_dataU, dims = 2), dims = 1))

#Aboslute RMSE:
sqrt(sum((inferU .- vraiU) .^ 2)/3)
sqrt(sum((finU .- vraiU) .^ 2)/3)

#Relative RMSE:
sqrt(sum(((inferU .- vraiU) ./ vraiU) .^ 2)/3)
sqrt(sum(((finU .- vraiU)./ vraiU) .^ 2)/3)

##Exponential model
prop_exp(reads_dataU, Lambda)

#_____
optT = optim_Qmetabar(reads_dataT, Lambda,
  m0tot = 1e5,
  ncycles = 40,
  K = K,
  dispersion = 1.,
  nsim = nsim)

optT.minimizer ./ sum(optT.minimizer)
reshape([41, 30, 67] ./ sum([41, 30, 67]), (1,3))
mean(reads_dataT ./ sum(reads_dataT, dims = 2), dims = 1)

inferT = vec(optT.minimizer ./ sum(optT.minimizer))
vraiT = [20, 14, 33] ./ sum([20, 14, 33])
finT = vec(mean(reads_dataT ./ sum(reads_dataT, dims = 2), dims = 1))

#Absolute RMSE:
sqrt(sum((inferT .- vraiT) .^ 2)/3)
sqrt(sum((finT .- vraiT) .^ 2)/3)

#Relative RMSE:
sqrt(sum(((inferT .- vraiT) ./ vraiT) .^ 2)/3)
sqrt(sum(((finT .- vraiT)./ vraiT) .^ 2)/3)

##Exponential model
prop_exp(reads_dataT, Lambda)

#_____
optG = optim_Qmetabar(reads_dataG, Lambda,
  m0tot = 2.5e5,
  ncycles = 40,
  K = K,
  dispersion = 1.,
  nsim = nsim)

optG.minimizer ./ sum(optG.minimizer)
reshape([1, 4, 16] ./ 21, (1,3))
mean(reads_dataG ./ sum(reads_dataG, dims = 2), dims = 1)

inferG = vec(optG.minimizer ./ sum(optG.minimizer))
vraiG = [1, 4, 16] ./ 21
finG = vec(mean(reads_dataG ./ sum(reads_dataG, dims = 2), dims = 1))

#Absolute RMSE:
sqrt(sum((inferG .- vraiG) .^ 2)/3)
sqrt(sum((finG .- vraiG) .^ 2)/3)

#Relative RMSE:
sqrt(sum(((inferG .- vraiG) ./ vraiG) .^ 2)/3)
sqrt(sum(((finG .- vraiG)./ vraiG) .^ 2)/3)

##Exponential model
prop_exp(reads_dataT, Lambda)


##____________________________________________________________________________________________________
#Using either CbeA or CbeB

LambdaA = [0.972, 0.989, 0.940]
LambdaB = [0.947, 0.989, 0.940]

#Change LambdaA to LambdaB to compare

#_____
optU = optim_Qmetabar(reads_dataU, LambdaA,
  m0tot = 1e5,
  ncycles = 40,
  K = K,
  dispersion = 1.,
  nsim = nsim)

optU.minimizer ./ sum(optU.minimizer)
reshape(fill(1, 3) ./ 3, (1,3))
mean(reads_dataU ./ sum(reads_dataU, dims = 2), dims = 1)

inferU = vec(optU.minimizer ./ sum(optU.minimizer))
vraiU = fill(1, 3) ./ 3
finU = vec(mean(reads_dataU ./ sum(reads_dataU, dims = 2), dims = 1))

#Aboslute RMSE:
sqrt(sum((inferU .- vraiU) .^ 2)/3)
sqrt(sum((finU .- vraiU) .^ 2)/3)

#Relative RMSE:
sqrt(sum(((inferU .- vraiU) ./ vraiU) .^ 2)/3)
sqrt(sum(((finU .- vraiU)./ vraiU) .^ 2)/3)

#_____
optT = optim_Qmetabar(reads_dataT, LambdaA,
  m0tot = 1e5,
  ncycles = 40,
  K = K,
  dispersion = 1.,
  nsim = nsim)

optT.minimizer ./ sum(optT.minimizer)
reshape([41, 30, 67] ./ sum([41, 30, 67]), (1,3))
mean(reads_dataT ./ sum(reads_dataT, dims = 2), dims = 1)

inferT = vec(optT.minimizer ./ sum(optT.minimizer))
vraiT = [41, 30, 67] ./ sum([41, 30, 67])
finT = vec(mean(reads_dataT ./ sum(reads_dataT, dims = 2), dims = 1))

#Absolute RMSE:
sqrt(sum((inferT .- vraiT) .^ 2)/3)
sqrt(sum((finT .- vraiT) .^ 2)/3)

#Relative RMSE:
sqrt(sum(((inferT .- vraiT) ./ vraiT) .^ 2)/3)
sqrt(sum(((finT .- vraiT)./ vraiT) .^ 2)/3)

#_____
optG = optim_Qmetabar(reads_dataG, LambdaA,
  m0tot = 1e5,
  ncycles = 40,
  K = K,
  dispersion = 1.,
  nsim = nsim)

optG.minimizer ./ sum(optG.minimizer)
reshape([12, 24, 190] ./ sum([12, 24, 190]), (1,3))
mean(reads_dataG ./ sum(reads_dataG, dims = 2), dims = 1)

inferG = vec(optG.minimizer ./ sum(optG.minimizer))
vraiG = [12, 24, 190] ./ sum([12, 24, 190])
finG = vec(mean(reads_dataG ./ sum(reads_dataG, dims = 2), dims = 1))

#Absolute RMSE:
sqrt(sum((inferG .- vraiG) .^ 2)/3)
sqrt(sum((finG .- vraiG) .^ 2)/3)

#Relative RMSE:
sqrt(sum(((inferG .- vraiG) ./ vraiG) .^ 2)/3)
sqrt(sum(((finG .- vraiG)./ vraiG) .^ 2)/3)


##____________________________________________________________________________________________________
#Check type stability
@code_warntype optim_Qmetabar(reads_dataU, Lambda,
    m0tot = 1e5,
    ncycles = 40,
    K = 1e11,
    dispersion = 1.,
    nsim = 1000)

##____________________________________________________________________________________________________
# Mechanistic model

optU_mec = optim_Qmetabar(reads_dataU, Lambda,
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
reshape(fill(1, 3) ./ 3, (1,3))
mean(reads_dataU ./ sum(reads_dataU, dims = 2), dims = 1)

optT_mec = optim_Qmetabar(reads_dataT, Lambda,
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
reshape([41, 30, 67] ./ sum([41, 30, 67]), (1,3))
mean(reads_dataT ./ sum(reads_dataT, dims = 2), dims = 1)


optG_mec = optim_Qmetabar(reads_dataG, Lambda,
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
reshape([12, 24, 190] ./ sum([12, 24, 190]), (1,3))
mean(reads_dataG ./ sum(reads_dataG, dims = 2), dims = 1)

##____________________________________________________________________________________________________
#Adding an "hybrid" species accounting for every other species
Lambda_hybrid = mean(Lambda)
Lambda_complete = vcat(Lambda, Lambda_hybrid)

reads_dataU_complete = Float64.(Matrix(CSV.read("data/dfSper01U_taq_complete.csv", DataFrame)[:,2:5]))
reads_dataT_complete = Float64.(Matrix(CSV.read("data/dfSper01T_taq_complete.csv", DataFrame)[:,2:5]))
reads_dataG_complete = Float64.(Matrix(CSV.read("data/dfSper01G_taq_complete.csv", DataFrame)[:,2:5]))


#_____
optU_complete = optim_Qmetabar(reads_dataU_complete, Lambda_complete,
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
optT_complete = optim_Qmetabar(reads_dataT_complete, Lambda_complete,
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
optG_complete = optim_Qmetabar(reads_dataG_complete, Lambda_complete,
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

J_Qmetabar([1., 1., 1.], rand(190, 1000), reads_data = reads_dataU, Lambda = Lambda, ncycles = 40)
J_Qmetabar(optU.minimizer, rand(190, 1000), reads_data = reads_dataU, Lambda = Lambda, ncycles = 40)


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