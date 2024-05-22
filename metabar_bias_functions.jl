
##28/03/24
#Functions for the paper
#Sylvain Moinard

#From the final proportions of each species and the PCR amplifiation rates inferred or measured by Taqman qPCR,
#we infer the initial proportions of each species with the Fixed Landscape Inference MethOd (flimo)

#See https://git.metabarcoding.org/lecasofts/flimo

#Run this file before running inferences
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


#Function used in the mechanistic saturation model, not presented in paper
function gamma_factor(x ; c = 1.)
  2/(1+sqrt(1-4*c^2*x*(1-x)))
end

#The simulators are adapted to the flimo method, i.e.
#each random draw (rand(...)) is replaced by a call to the appropriate quantile function (quantile(...))

#The normal approximation of the model makes it possible to obtain a continuous distribution, and therefore to use
#an automatic differentiation algorithm to minimise the objective function


## PCR Simulation

#The implemented model is the logistic model
#
# Lambda is an array containing the amplification rates of each species (0 <= Lambda <= 1)
# K is the charge capacity
# n_reads is the number of reads after sequencing
# ncycles is the number of PCR cycles
# nreplicate is the number of replicates to simulate
# dispersion is the initial inter-replicate dispersion

#Theta: array containing the initial quantities
#quantiles: 2D array. Nrow = Number of simulations; Ncol = nspecies*(ncycles+2)

#Exact simulation of a PCR:
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
#mechanistic model, not presented in paper
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

#Exact simulations of a PCR with fixed quantiles (randomness management for flimo)
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



#Normally approximated simulation with fixed quantiles:
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
  
  return kinetics .* n_reads ./ sum(kinetics, dims = 2) #uncorrect conditionning*** -> multinomial
end



##Objective function
#To minimize to find optimal starting conditions

#reads_data: 2D array: Nrow : nreplicates Ncol : reads for each species


function J_ps(Theta,
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
##Optimisation

#quantiles : Nrow = number replicates
#Ncol = nspecies * (ncycles+2)
#first column : initial randomness, then each cycle then sequencing

function optim_ps(reads_data::Array{Float64, 2}, Lambda::Array{Float64, 1};
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
    c = 1., #mechanistic model, not shown in paper
    e = 1., #mechanistic model, not shown in paper
    show_trace = true,
    show_every = 10) 
  
    #reads_data: 2D array. Rows = replicates, Col = species
  
    #nreplicate = size(reads_data, 1)
    nspecies = size(reads_data, 2)
  
    n_reads_repli = sum(reads_data, dims = 2)
  
    pdata = vec(mean(reads_data ./ sum(reads_data, dims = 2), dims = 1)) #average species proportions in data
  
    if (isnothing(Theta0))
      Theta0 = pdata
    end
    obj = (Theta, quantiles) -> J_ps(Theta, quantiles,
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
                       show_trace = show_trace,
                       show_every = show_every)
      return opt
  end
  
  
  function prop_exp(reads_data, Lambda, low_c = 15, high_c = 30)
    #proportions with the exponential model for ncycles between low_c and high_c
    prop = zeros(length(low_c:high_c), size(reads_data, 2))
    for s in axes(reads_data, 2)
      prop[:,s] = mean(reads_data[:,s]) ./ (1+Lambda[s]) .^(low_c:high_c) .* K ./ mean(sum(reads_data, dims = 2))
    end
    prop ./ sum(prop, dims = 2)
  end

  ##____________________________________________________________________________________________________
  ## Inference of efficiencies

  # Inference Function

#Normally approximated simulation of a PCR with fixed quantiles:
function simu_pcr_normalQ_eff(Theta,
    quantiles ;
    m0 = fill(2e4, length(Theta)),
    K = 1e11,
    n_reads = 1e5,
    ncycles = 40,
    nreplicate = size(quantiles, 1),
    dispersion = 1.,
    model = "logistic",
    c = 1.,
    e = 1.)

    Lambda = Theta

    if length(n_reads) == 1
      n_reads = fill(1e5, nreplicate)
    end

    if length(n_reads) < nreplicate
      n_reads = vec(repeat(n_reads, Int64(nreplicate/length(n_reads)))) #multiple needed
    end

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
      #  println(c)
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

  #Sequencing
  kinetics_tot = sum(kinetics, dims = 2)
  for species in 1:nspecies
    #reuse crea array for pre-allocation
    crea[:] .= vec(view(kinetics,:, species) ./ kinetics_tot)
    kinetics[:,species] .= quantile.(Normal.(n_reads .* crea,
                                             sqrt.(n_reads .* crea .* (one(eltype(crea)) .- crea))),
                                    view(quantiles,1:nreplicate,(nspecies*(ncycles+1)+species)))
  end
  kinetics[kinetics .< zero(eltype(m0))] .= zero(eltype(m0))
  
  return kinetics .* n_reads ./ sum(kinetics, dims = 2)
end

#Cost function to minimize to find optimal lambda_s and K
function J_efficiencies(Theta,
    quantiles ;
    reads_data = rand(1,1)::Array{Float64, 2},
    m0 = qty_init_U,
    n_reads = 1e5,
    ncycles = 40,
    dispersion = 1.,
    model = "logistic",
    c = 1.,
    e = 1.,
    lambda_max = 1.,
    Kmult = 1e14,
    Kfixed = false)

    #One should choose nsim = k * nreplicate, k integer
    if Kfixed
      K = Kmult
      if !isnothing(lambda_max)
        Theta = vcat(lambda_max, Theta) #most abundant species in first position !!!
      end
    else
      K = 10^Theta[length(Theta)]*Kmult #log10(K) stored in Theta
      if !isnothing(lambda_max)
        Theta = vcat(lambda_max, vec(Theta[1:(length(Theta)-1)])) #most abundant species in first position !!!
      end
    end


  

  nsim = size(quantiles, 1) #number of simulations to perform by flimo

  nreplicate = size(reads_data, 1)

  #Ksim = Vector{eltype(Theta)}(undef, nsim) #estimated value of K for each simulation
  pdata = vec(mean(reads_data ./ sum(reads_data, dims = 2), dims = 1)) #average species proportions in data
  #pdata = vec(mean(reads_data, dims = 1))

  nsimr = Int64(floor(nsim/nreplicate)) #number of simulations by replicate

  n_reads_repli = Vector{eltype(m0)}(undef, nsim) #Number of reads for each simulation
  
  @inbounds for repl in 1:nreplicate
      #more simulations than replicates needed in this implementation
      n_reads_repli[(1+(repl-1)*nsimr):(repl*nsimr)] .= n_reads[repl]
    end
    if nsimr*nreplicate < nsim #unused in the score computation
      n_reads_repli[(1+nsimr*nreplicate):nsim] .= view(n_reads,1:nsim-nsimr*nreplicate)
    end

    kinetics = simu_pcr_normalQ_eff(Theta,
                                   quantiles ;
                                   m0 = m0,
                                   K = K,
                                   n_reads = n_reads_repli,
                                   ncycles = ncycles,
                                   nreplicate = nsim,
                                   dispersion = dispersion,
                                   model = model,
                                   c = c,
                                   e = e)
    
    kinetics = kinetics ./ sum(kinetics, dims = 2)
    #@inbounds for sim in 1:nsim #Number of molecules for each simulation
    #  cin_tot = sum(view(kinetics, sim, :))
    #  kinetics[sim,:] .= view(kinetics, sim,:)./cin_tot #Relative abundance
    #end

    @inbounds for sim in 1:nsim #Number of molecules for each simulation
        cin_tot = sum(view(kinetics, sim, :))
        kinetics[sim,:] .= view(kinetics, sim,:)./cin_tot #Relative abundance
      end
      
      psimu = vec(mean(kinetics, dims = 1))
      psimu[:] .= psimu ./ sum(psimu)

      
      #pdata = reads_data ./ sum(reads_data, dims = 2)#without agregating ps_data
      score = sum(((psimu .- pdata).^2. ./ pdata))
      #score = sum((transpose(pdata) .- psimu).^2 ./ transpose(pdata))#without agregating ps_data

    #score = zero(eltype(Theta))

    #@inbounds for repl in 1:nreplicate
    #    pdata = view(reads_data, repl, :) ./ sum(view(reads_data, repl, :))
    #    psimu = vec(mean(view(kinetics, (1+(repl-1)*nsimr):(repl*nsimr), :), dims = 1))
    #   score = score + sum((psimu .- pdata).^2. ./ pdata)
    #  end
    
  return score
end



#Inference function to get the PCR efficiencies Lambda_s and K
function optim_efficiencies(reads_data::Array{Float64, 2}, m0::Array{Float64, 1};
    ncycles = 40::Int64,
    nsim = 1000::Int64,
    ninfer = 1,
    Theta0 = nothing,
    maxit = 200::Int64,
    gtol = 1e-6::Float64,
    xtol=1e-32::Float64,
    dispersion = 1.::Float64,
    model = "logistic",
    c = 1.,
    e = 1.,
    lambda_max = 1.,
    quantiles = nothing,
    Kfixed = false,
    K = 1e14,
    Kmin_log = -2,
    Kmax_log = 2,
    show_trace = true,
    show_every = 10)
  
    #reads_data: 2D array. Rows = replicates, Col = species
    nspecies = size(reads_data, 2)
  
    #Estimation of K for the data replicates
    n_reads_repli = sum(reads_data, dims = 2)
  
    #pdata = vec(mean(reads_data ./ sum(reads_data, dims = 2), dims = 1)) #average species proportions in data
    if Kfixed
      if (isnothing(Theta0))
        Theta0 = fill(lambda_max*0.999, nspecies-1)
      end
        theta_lower = fill(0.5, nspecies-1)
        theta_upper = fill(lambda_max*0.9999, nspecies-1)
    else  
      if (isnothing(Theta0))
        Theta0 = vcat(fill(lambda_max*0.999, nspecies-1), 0.) #log10(K)
      end
        theta_lower = vcat(fill(0.5, nspecies-1), Kmin_log)
        theta_upper = vcat(fill(lambda_max*0.9999, nspecies-1), Kmax_log)
    end

    obj = (Theta, quantiles) -> J_efficiencies(Theta, quantiles,
      reads_data = reads_data,
      m0 = m0,
      n_reads = n_reads_repli,
      ncycles = ncycles,
      dispersion = dispersion,
      model = model,
      c = c,
      e = e,
      lambda_max = lambda_max,
      Kfixed = Kfixed,
      Kmult = K)

    opt = Jflimo.flimoptim(nspecies*(ncycles+2),
                       obj = obj,
                       nsim = nsim,
                       ninfer = ninfer,
                       lower = theta_lower,
                       upper = theta_upper,
                       AD = true,
                       randomTheta0 = true,
                       previousTheta0 = false,
                       gtol = gtol,
                       xtol = xtol,
                       maxit = maxit,
                       show_trace = show_trace,
                       show_every = show_every,
                       quantiles = quantiles)
      return opt
  end

##____________________________________________________________________________________________________
## RMSE computation

function AbsErr(ps_observed, ps_expected)
  sqrt(sum((ps_observed .- ps_expected).^2) / length(ps_expected))
end

function RelErr(ps_observed, ps_expected)
  sqrt(sum(((ps_observed .- ps_expected) ./ ps_expected).^2) / length(ps_expected))
end

  #TO DO 14/11/2023 : implement quantile multinomial (with cumsum)