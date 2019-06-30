# Causal Inference over Mixtures (CIM)

Many real causal processes contain cycles and evolve. However, most causal discovery algorithms assume that the underlying causal process follows a single directed acyclic graph (DAG) that does not change over time. This repository contains code for a new algorithm called Causal Inference over Mixtures (CIM) which relaxes the single DAG assumption by modeling causal processes using a mixture of DAGs so that the graph can change over time. CIM uses longitudinal data to improves the accuracy of causal discovery on both real and synthetic clinical datasets even when cycles, non-stationarity, non-linearity, latent variables and selection bias exist simultaneously.

# Installation

The package depends on the MASS and pcalg packages on CRAN, so please install those first. Then:

> library(devtools)

> install_github("ericstrobl/CIM"); install_github("ericstrobl/RCIT")

> library(CIM); library(RCIT); library(pcalg)

# Generate Synthetic Data

> waves = list(w1=1:8,w2=9:16,w3=17:24) # create 3 waves containing 8 variables each

> mDAGs = generate_mix_DAGs(nIndep=sample(5:15,1),p=24,en=2,waves=waves) # generate a mixture of DAGs across the 3 waves with latent and selection variables

> waves2= list(w1=1:length(mDAGs$waves$w1),
               w2=(length(mDAGs$waves$w1)+1):(length(mDAGs$waves$w1)+length(mDAGs$waves$w2)),
               w3=(length(mDAGs$waves$w1)+length(mDAGs$waves$w2)+1):(length(mDAGs$waves$w1)+length(mDAGs$waves$w2)+length(mDAGs$waves$w3))
               ) # reorganize waves according to generated latent and selection variables

> plot(as(synth_list$mDAGs$graph,"graphNEL")) # plot the ground truth father graph

> synth_list= sample_mix_DAGs(mDAGs,1000);  suffStat=list(); suffStat$data = synth_list$data; # generate 1000 samples from the mixture of DAGs with latent and selection variables

# Run CIM on Synthetic Data

> out = CIM(suffStat, RCoT_wrap, alpha=0.01, p=ncol(suffStat$data), waves=waves) # run CIM

> colnames(out$maag) <- mDAGs$actual_indices; rownames(out$maag) <- mDAGs$actual_indices # modify indices to account for possible latent and selection variables

> print(out$maag) # print recovered partially oriented father graph

