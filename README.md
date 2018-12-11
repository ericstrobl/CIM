# Causal Inference over Mixtures (CIM)

CIM is a causal discovery algorithm that can infer causal structure even when cycles, latent variables, selection bias, non-linearity and non-stationarity exist simultaneously. 

CIM assumes that the joint distribution can be modeled as a mixture of DAGs. Why? Because it turns out that most other causal discovery algorithms orient too many arrowheads in practice, particularly if they allow latent variables and selection bias. See it for yourself:



A possible explanation for this phenomenon is that real joint distributions actually arise from multiple causal processes (i.e., a mixture of DAGs) rather than just a single causal process (i.e., a single directed graph).


CIM was designed to accurately discover causal relationships from clinical data, where causal processes often contain cycles and evolve over time. I suspect that CIM works well with other data types as well.

# Causal Inference over Mixtures (CIM)


# Run the Algorithm on Clinical Datasets
