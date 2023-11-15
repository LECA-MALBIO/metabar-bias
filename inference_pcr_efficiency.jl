
##Inference of PCR efficiencies in metabarcoding data
##Sylvain Moinard, LECA
## 2023

## Uses function in file inference_metabarcoding.jl

##____________________________________________________________________________________________________
# Setup

# U community : reads

reads_U = Float64.(Matrix(CSV.read("data/reads_U.csv", DataFrame)))

# Initial number of molecules, assayed by ddPCR

qty_init_U = fill(19084., size(reads_U, 2))

# Inference Function

#Normally approximated simulation:
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


function J_efficiencies(Theta,
    quantiles ;
    reads_data = rand(1,1)::Array{Float64, 2},
    m0 = qty_init_U,
    K = 1e11,
    n_reads = 1e5,
    ncycles = 40,
    dispersion = 1.,
    model = "logistic",
    c = 1.,
    e = 1.,
    lambda_max = 1.)

    #One should choose nsim = k * nreplicate, k integer

  Theta = vcat(lambda_max, vec(Theta)) #most abundant species in first position !!!
  println(Theta)

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

      println(pdata)
      println(psimu)
      
      score = sum(((psimu .- pdata).^2. ./ pdata))

    #score = zero(eltype(Theta))

    #@inbounds for repl in 1:nreplicate
    #    pdata = view(reads_data, repl, :) ./ sum(view(reads_data, repl, :))
    #    psimu = vec(mean(view(kinetics, (1+(repl-1)*nsimr):(repl*nsimr), :), dims = 1))
    #   score = score + sum((psimu .- pdata).^2. ./ pdata)
    #  end
    
  return score
end


function optim_efficiencies(reads_data::Array{Float64, 2}, m0::Array{Float64, 1};
    ncycles = 40::Int64,
    K = 1e11::Float64,
    nsim = 1000::Int64,
    Theta0 = nothing,
    maxit = 200::Int64,
    gtol = 1e-6::Float64,
    xtol=1e-32::Float64,
    dispersion = 1.::Float64,
    model = "logistic",
    c = 1.,
    e = 1.,
    lambda_max = 1.,
    quantiles = nothing)
  
    #reads_data: 2D array. Rows = replicates, Col = species
    nspecies = size(reads_data, 2)
  
    #Estimation of K for the data replicates
    n_reads_repli = sum(reads_data, dims = 2)
  
    #pdata = vec(mean(reads_data ./ sum(reads_data, dims = 2), dims = 1)) #average species proportions in data
    if (isnothing(Theta0))
      Theta0 = fill(lambda_max*0.999, nspecies-1)
    end
    obj = (Theta, quantiles) -> J_efficiencies(Theta, quantiles,
                                          reads_data = reads_data,
                                          m0 = m0,
                                          K = K,
                                          n_reads = n_reads_repli,
                                          ncycles = ncycles,
                                          dispersion = dispersion,
                                          model = model,
                                          c = c,
                                          e = e,
                                          lambda_max = lambda_max)
    
    opt = Jflimo.flimoptim(nspecies*(ncycles+2),
                       obj = obj,
                       nsim = nsim,
                       ninfer = 1,
                       lower = fill(0.5, nspecies-1),
                       upper = fill(lambda_max*0.9999, nspecies-1),
                       AD = true,
                       Theta0 = Theta0,
                       randomTheta0 = false,
                       gtol = gtol,
                       xtol = xtol,
                       maxit = maxit,
                       show_trace = true,
                       show_every = 10,
                       quantiles = quantiles)
      return opt
  end

##____________________________________________________________________________________________________
# Efficiencies inference

#Define K
K = 7.57e12

nreplicate = size(reads_U, 1)

n_reads = sum(reads_U, dims = 2)

## Optimize
#nsim proportional to nreplicates
nsim = 1900

optU_eff = optim_efficiencies(reads_U, qty_init_U,
  ncycles = 40,
  K = K,
  dispersion = 1.,
  nsim = nsim)

optU_eff.minimizer
optU_eff.minimum
optU_eff.time

res_lambda = vcat(1, vec(optU_eff.minimizer))
# CSV.write("data/efficiencies_U.csv", Tables.table(res_lambda))

simU = simu_pcr_normalQ_eff(vcat(0.9, vec(optU_eff.minimizer)),
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

Lambda_infer = vcat(1., vec(optU_eff.minimizer))

reads_T = Float64.(Matrix(CSV.read("data/reads_T.csv", DataFrame)))
reads_G = Float64.(Matrix(CSV.read("data/reads_G.csv", DataFrame)))

optU_full = optim_Qmetabar(reads_U, Lambda_infer,
  m0tot = 1e5,
  ncycles = 40,
  K = K,
  dispersion = 1.,
  nsim = 190)

optU_full.minimizer ./ sum(optU_full.minimizer)
mean(reads_U ./ sum(reads_U, dims = 2), dims = 1)

optT_full = optim_Qmetabar(reads_T, Lambda_infer,
  m0tot = 1e5,
  ncycles = 40,
  K = K,
  dispersion = 1.,
  nsim = 190)

optT_full.minimizer ./ sum(optT_full.minimizer)
mean(reads_T ./ sum(reads_T, dims = 2), dims = 1)


optG_full = optim_Qmetabar(reads_G, Lambda_infer,
  m0tot = 1e5,
  ncycles = 40,
  K = K,
  dispersion = 1.,
  nsim = 190)

vec(optG_full.minimizer ./ sum(optG_full.minimizer))
mean(reads_G ./ sum(reads_G, dims = 2), dims = 1)

CSV.write("data/prop_inferU.csv", Tables.table(optU_full.minimizer ./ sum(optU_full.minimizer)))
CSV.write("data/prop_inferT.csv", Tables.table(optT_full.minimizer ./ sum(optT_full.minimizer)))
CSV.write("data/prop_inferG.csv", Tables.table(optG_full.minimizer ./ sum(optG_full.minimizer)))

## Tests
##____________________________________________________________________________________________________

Random.seed!(1234)
Q = rand(190, 13*42)

J_efficiencies(optU_eff.minimizer, Q, reads_data = reads_U, K = K, n_reads = sum(reads_U, dims = 2))
J_efficiencies(fill(0.825, 12) , Q, reads_data = reads_U, K = K, n_reads = sum(reads_U, dims = 2))

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
    K = K,
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
