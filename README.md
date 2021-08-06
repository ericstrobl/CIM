# Causal Inference over Mixtures (CIM)

This repository contains code for an algorithm called Causal Inference over Mixtures (CIM) which relaxes the single DAG assumption by modeling causal processes using a mixture of DAGs, so that the graph and causal relations can change over time and sub-populations. CIM uses longitudinal data to improve the accuracy of causal discovery on both real and synthetic clinical datasets. Each time step in the longitudinal dataset may correspond to a mixture of multiple DAGs. CIM accurately recovers causal relations even when cycles, non-stationarity, non-linearity, latent variables and selection bias exist simultaneously.

The Experiments folder contains code needed to replicate the synthetic data results.

# Installation

> library(devtools)

> install_github("ericstrobl/CIM")

> library(CIM)

# Sample from a Mixture of DAGs

> waves = list(w1=1:8,w2=9:16,w3=17:24) # create 3 waves, aka time steps, containing 8 variables each

> mixDAG = generate_mix_DAGs2(24,en=2,waves) # generate a mixture DAGs, also include latent and selection variables

> resort_p = sample(c(mixDAG$waves_L$w1,mixDAG$waves_L$w2,mixDAG$waves_L$w3),24-length(mixDAG$L),replace=FALSE) # remove latent variables and randomize variable order

> waves = list(w1 = match(mixDAG$waves_L$w1,resort_p), w2 = match(mixDAG$waves_L$w2,resort_p), w3 = match(mixDAG$waves_L$w3,resort_p)) # wave prior knowledge

> synth_data = sample_mix_DAGs2(mixDAG,samps) # sample from the mixture of DAGs

# Run CIM on the Data

> suffStat = list(); suffStat$data = synth_data[,resort_p];

> out = CIM(suffStat, GCM_wrap, alpha=0.01, p=ncol(suffStat$data), waves=waves) # run CIM

> print(out$f_star) # print F*


# How to Interpret the Output

Let S denote the selection variables.

`out$f_star[i,j] = 0` means that CIM could find a set rendering i and j conditionally independent

`out$f_star[i,j] = 1` means CIM could *not* find a set rendering i and j conditionally independent, and CIM does *not* know if j is an ancestor or not an ancestor of i or S in F*

`out$f_star[i,j] = 2` means j is *not* an ancestor of i in F*

`out$f_star[i,j] = 3` means j is an ancestor of i or S in F*

