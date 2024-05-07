
##16/04/2024
##Sylvain Moinard, LECA
##Inference from two reference communities

##____________________________________________________________________________________________________
##Import reads
reads_U = Float64.(Matrix(CSV.read("data/reads_U.csv", DataFrame)))
reads_T = Float64.(Matrix(CSV.read("data/reads_T.csv", DataFrame)))
reads_G = Float64.(Matrix(CSV.read("data/reads_G.csv", DataFrame)))

# Initial number of molecules, assayed by ddPCR

qty_init_U = vec(Float64.(Matrix(CSV.read("data/qty_initU.csv", DataFrame))))
qty_init_T = vec(Float64.(Matrix(CSV.read("data/qty_initT.csv", DataFrame))))
qty_init_G = vec(Float64.(Matrix(CSV.read("data/qty_initG.csv", DataFrame))))

##____________________________________________________________________________________________________
##Functions

#Inference function to get the PCR efficiencies Lambda_s and K
function optim_efficiencies_2commu(
    reads_data1::Array{Float64, 2},
    reads_data2::Array{Float64, 2},
    m01::Array{Float64, 1},
    m02::Array{Float64, 1};
    ncycles1 = 40::Int64,
    ncycles2 = 40::Int64, #warning: assumes here that ncycles1 = ncycles2
    nsim = 1000::Int64,
    ninfer = 1,
    Theta0 = nothing,
    maxit = 200::Int64,
    gtol = 1e-6::Float64,
    xtol=1e-32::Float64,
    dispersion = 1.::Float64,
    model = "logistic",
    c = 1., #mechanistic model, not shown in paper
    e = 1., #mechanistic model, not shown in paper
    lambda_ref = 0.9,
    quantiles = nothing,
    Kmin_log = -2,
    Kmax_log = 2,
    Kmult = 1e12)
  
    #reads_data: 2D array. Rows = replicates, Col = species
    nspecies = size(reads_data1, 2) #same for data2
  
    #Estimation of K for the data replicates
    n_reads_repli1 = sum(reads_data1, dims = 2)
    n_reads_repli2 = sum(reads_data2, dims = 2)
  
    #pdata = vec(mean(reads_data ./ sum(reads_data, dims = 2), dims = 1)) #average species proportions in data
    if (isnothing(Theta0))
      Theta0 = vcat(fill(lambda_ref, nspecies), 0.) #last: log10(K)
    end
    obj = (Theta, quantiles) -> J_efficiencies(Theta, quantiles[:,1:nspecies*(ncycles1+2)],
                                          reads_data = reads_data1,
                                          m0 = m01,
                                          n_reads = n_reads_repli1,
                                          ncycles = ncycles1,
                                          dispersion = dispersion,
                                          model = model,
                                          c = c,
                                          e = e,
                                          lambda_max = nothing,
                                          Kmult = Kmult)+
                                J_efficiencies(Theta, quantiles[:,(1+nspecies*(ncycles1+2):(2*nspecies*(ncycles1+2)))],
                                          reads_data = reads_data2,
                                          m0 = m02,
                                          n_reads = n_reads_repli2,
                                          ncycles = ncycles2,
                                          dispersion = dispersion,
                                          model = model,
                                          c = c,
                                          e = e,
                                          lambda_max = nothing,
                                          Kmult = Kmult)
    
    opt = Jflimo.flimoptim(2*nspecies*(ncycles1+2), #number of random draws per simulation (*2 commu)
                       obj = obj,
                       nsim = nsim,
                       ninfer = ninfer,
                       lower = vcat(fill(0.5, nspecies), Kmin_log),
                       upper = vcat(fill(0.9999, nspecies), Kmax_log),
                       AD = true,
                       randomTheta0 = true,
                       previousTheta0 = false,
                       gtol = gtol,
                       xtol = xtol,
                       maxit = maxit,
                       show_trace = true,
                       show_every = 10,
                       quantiles = quantiles)
      return opt
  end



##____________________________________________________________________________________________________
##M_U and M_T

## Optimize
#nsim proportional to nreplicates (19 each)
nsim = 190

optU_eff2UT = optim_efficiencies_2commu(reads_U, reads_T,
    qty_init_U, qty_init_T,
    ncycles1 = 40,
    ncycles2 = 40,
    ninfer = 11,
    nsim = nsim,
    Kmult = 1e13,
    Kmin_log = -13., #log10
    Kmax_log = 3.) #log10


optU_eff2UT.minimizer
optU_eff2UT.minimum



##Inference for M_G

##____________________________________________________________________________________________________
##M_U and M_G

##Inference for M_T