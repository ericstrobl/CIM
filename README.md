# Causal Inference over Mixtures (CIM)

This repository contains code for an algorithm called Causal Inference over Mixtures (CIM) which relaxes the single DAG assumption by modeling causal processes using a mixture of DAGs so that the graph can change over time. CIM uses longitudinal data to improve the accuracy of causal discovery on both real and synthetic clinical datasets even when cycles, non-stationarity, non-linearity, latent variables and selection bias exist simultaneously.

# Installation

The package depends on the MASS, pcalg and igraph packages, so please install those first. Then:

> library(devtools)

> install_github("ericstrobl/CIM")

> library(CIM)

# Generate Synthetic Data

> waves = list(w1=1:8,w2=9:16,w3=17:24) # create 3 waves containing 8 variables each

> mixDAG = generate_mix_DAGs2(24,en=2,waves) # generate a mixture DAGs, also include latent and selection variables

> synth_data = sample_mix_DAGs2(mixDAG,samps) # sample from the mixture of DAGs

# Run CIM on Synthetic Data

> out = cim_out = CIM(suffStat, GCM_wrap, alpha=0.01, p=ncol(suffStat$data), waves=waves) # run CIM

> print(out$f_star) # print F*


# How to Interpret the Output

Let S denote the selection variables.

`out$pofaag[i,j] = 0` means that CIM could render i and j conditionally independent

`out$pofaag[i,j] = 1` means CIM could *not* render i and j conditionally independent, and CIM does *not* know if j is an ancestor or not an ancestor of i or S

`out$pofaag[i,j] = 2` means j is *not* an ancestor of i in F*

`out$pofaag[i,j] = 3` means j is an ancestor of i or S in F*

