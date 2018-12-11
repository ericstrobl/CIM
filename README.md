# Causal Inference over Mixtures (CIM)

CIM is a causal discovery algorithm that can infer causal structure even when cycles, latent variables, selection bias, non-linearity and non-stationarity exist simultaneously. CIM assumes that the joint distribution can be modeled as a mixture of DAGs. The algorithm was designed to accurately discover causal relationships from clinical data, where causal processes often contain cycles and evolve over time. I suspect that the method works well with other data types as well.
