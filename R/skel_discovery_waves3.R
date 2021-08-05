#' Skeleton discovery with longitudinal data, wrapper function

skel_discovery_waves3 <- function (suffStat, indepTest, alpha, p, waves, max.cs = Inf) 
{
  time_start <- proc.time()
  
  skel <- skeleton_waves(suffStat, indepTest, alpha, waves, p = p, m.max = max.cs)
  G_sk <- as(skel@graph, "matrix")
  sepset_sk <- skel@sepset
  time_skel = proc.time() - time_start

  res_sk <- list(G_sk = G_sk, sepset_sk = sepset_sk,
                time_skel = time_skel, skel = skel)
  return(res_sk)
}