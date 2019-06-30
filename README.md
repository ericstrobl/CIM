# Causal Inference over Mixtures (CIM)

This repository contains code for an algorithm called Causal Inference over Mixtures (CIM) which relaxes the single DAG assumption by modeling causal processes using a mixture of DAGs so that the graph can change over time. CIM uses longitudinal data to improve the accuracy of causal discovery on both real and synthetic clinical datasets even when cycles, non-stationarity, non-linearity, latent variables and selection bias exist simultaneously.

# Installation

The package depends on the MASS, pcalg and igraph packages, so please install those first. Then:

> library(devtools)

> install_github("ericstrobl/CIM"); install_github("ericstrobl/RCIT")

> library(CIM); library(RCIT); library(pcalg)

# Generate Synthetic Data

> waves = list(w1=1:8,w2=9:16,w3=17:24) # create 3 waves containing 8 variables each

> mDAGs = generate_mix_DAGs(nIndep=sample(5:15,1),p=24,en=2,waves=waves) # generate a mixture of 5-15 DAGs with the waves; also include latent and selection variables

> synth_list= sample_mix_DAGs(mDAGs,1000);  suffStat=list(); suffStat$data = synth_list$data; # generate 1000 samples from the mixture of DAGs with latent and selection variables

> plot(as(synth_list$mDAGs$graph,"graphNEL")) # plot the ground truth father graph

# Run CIM on Synthetic Data

> out = CIM(suffStat, RCoT_wrap, alpha=0.01, p=ncol(suffStat$data), waves=mDAGs$waves) # run CIM

> colnames(out$pofaag) <- mDAGs$actual_indices; rownames(out$pofaag) <- mDAGs$actual_indices # modify indices to account for possible latent and selection variables

> print(out$pofaag) # print recovered partially oriented father AAG


# How to Interpret the Output

Let S denote the selection variables.

`out$pofaag[i,j] = 0` means that CIM could render i and j conditionally independent

`out$pofaag[i,j] = 1` means CIM could *not* render i and j conditionally independent, and CIM does *not* know if j is an ancestor or not an ancestor of i or S

`out$pofaag[i,j] = 2` means j is *not* an ancestor of i in the father AAG

`out$pofaag[i,j] = 3` means j is an ancestor of i or S in the father AAG

