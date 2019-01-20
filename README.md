# Causal Inference over Mixtures (CIM)

Real causal processes often do not follow a single directed graph because causal processes usually change over time. CIM is an algorithm for inferring causation from non-experimental data when the causal process changes. 

CIM requires longitudinal data or some other prior knowledge to rule out causal relationships. CIM outperforms PC, FCI, RFCI and CCI by a large margin on average, even if we give the other algorithms the same time information or prior knowledge. The algorithm was inspired by the observation that X and Y are often independent in real non-experimental data when X directly causes Y.

CIM can also handle cycles, latent variables, selection bias, and non-linearity simultaneously; every effort was made to make the algorithm actually work on real data. This repository allows you install and then run CIM with a few lines of codes.

# Installation


# Run the Algorithm on Real Longitudinal Data
- Framingham Heart Study:

- Mayo Clinic Primary Biliary Cirrhosis Study:

- 
