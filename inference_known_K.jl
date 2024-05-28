
##Inference of PCR efficiencies in metabarcoding data, known K
##Sylvain Moinard, LECA
## 28/05/2024

#From the final proportions of each species and the PCR amplifiation rates inferred or measured by Taqman qPCR,
#we infer the initial proportions of each species with the Fixed Landscape Inference MethOd (flimo)

#See https://git.metabarcoding.org/lecasofts/flimo

## Use functions in file metabar_bias_functions.jl

##____________________________________________________________________________________________________
# Setup: load data
reads_U = Float64.(Matrix(CSV.read("data/export_to_julia/reads_U.csv", DataFrame)))
reads_T = Float64.(Matrix(CSV.read("data/export_to_julia/reads_T.csv", DataFrame)))
reads_G = Float64.(Matrix(CSV.read("data/export_to_julia/reads_G.csv", DataFrame)))

# Initial number of molecules, assayed by ddPCR

qty_init_U = vec(Float64.(Matrix(CSV.read("data/export_to_julia/qty_initU.csv", DataFrame))))
qty_init_T = vec(Float64.(Matrix(CSV.read("data/export_to_julia/qty_initT.csv", DataFrame))))
qty_init_G = vec(Float64.(Matrix(CSV.read("data/export_to_julia/qty_initG.csv", DataFrame))))


##____________________________________________________________________________________________________
##Testing values of K 


K = 1.204e13

function full_inference_K(K, id ; nsim = 190)
    
    Random.seed!(21011996+id)
    optU_eff = optim_efficiencies(reads_U, qty_init_U,
      ncycles = 40,
      ninfer = 10,
      nsim = nsim,
      Kfixed = true,
      K = K,
      show_trace = false) 
    
    Lambda_infer = vcat(1., vec(mean(optU_eff.minimizer, dims  = 1)))
    CSV.write("data/export_to_r/efficiencies_U_K"*string(id)*".csv", Tables.table(optU_eff.minimizer))
    
    psU = mean(reads_U ./ sum(reads_U, dims = 2), dims = 1)
    
    neffU = (mean(log.(psU))+log(K)-mean(log.(qty_init_U)))/mean(log.(1 .+ Lambda_infer))
    
    psT = mean(reads_T ./ sum(reads_T, dims = 2), dims = 1)
    psG = mean(reads_G ./ sum(reads_G, dims = 2), dims = 1)
    
    m0totT = K .* sum(vec(psT) ./ (1 .+ Lambda_infer).^neffU)
    m0totG = K .* sum(vec(psG) ./ (1 .+ Lambda_infer).^neffU)
    
    optT_full = optim_ps(reads_T, Lambda_infer,
      ninfer = 10,
      randomTheta0 = true,
      previousTheta0 = false, 
      m0tot = m0totT,
      ncycles = 40,
      K = K,
      nsim = nsim,
      show_trace = false)

      CSV.write("data/export_to_r/prop_inferT_K"*string(id)*".csv", Tables.table(optT_full.minimizer ./ sum(optT_full.minimizer, dims = 2)))

    
    optG_full = optim_ps(reads_G, Lambda_infer,
      ninfer = 10,
      randomTheta0 = true,
      previousTheta0 = false, 
      m0tot = m0totG,
      ncycles = 40,
      K = K,
      nsim = nsim,
      show_trace = false)

    CSV.write("data/export_to_r/prop_inferG_K"*string(id)*".csv", Tables.table(optG_full.minimizer ./ sum(optG_full.minimizer, dims = 2)))

    
    #M_T
    at = AbsErr(vec(mean(optT_full.minimizer ./ sum(optT_full.minimizer, dims = 2), dims = 1)), qty_init_T ./ sum(qty_init_T))
    rt = RelErr(vec(mean(optT_full.minimizer ./ sum(optT_full.minimizer, dims = 2), dims = 1)), qty_init_T ./ sum(qty_init_T))
    
    #M_G
    ag = AbsErr(vec(mean(optG_full.minimizer ./ sum(optG_full.minimizer, dims = 2), dims = 1)), qty_init_G ./ sum(qty_init_G))
    rg = RelErr(vec(mean(optG_full.minimizer ./ sum(optG_full.minimizer, dims = 2), dims = 1)), qty_init_G ./ sum(qty_init_G))
    
    #CSV.write("data/at_rt_ag_rg_K"*string(id)*".csv", Tables.table(vcat(at, rt, ag, rg)))
    return vcat(at, rt, ag, rg)
end

rmseK1 = full_inference_K(K, 13)
rmseK2 = full_inference_K(10*K, 14)
rmseK2 = full_inference_K(100*K, 15)
rmseK3 = full_inference_K(K/10, 12)
rmseK4 = full_inference_K(K/100, 11)

rmseK1
rmseK2
rmseK3
rmseK4
rmseK5
