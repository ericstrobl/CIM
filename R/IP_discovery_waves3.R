IP_discovery_waves3 <- function (suffStat, skel, indepTest, alpha, p, waves, max.cs = Inf) 
{
  time_start <- proc.time()
  # skel <- skeleton_waves1(suffStat, indepTest, alpha, waves, p = p, m.max = max.cs,
  #                        method = "stable")
  G_sk <- as(skel@graph, "matrix")
  sepset_sk <- skel@sepset
  # time_skel = proc.time() - time_start
  pdsepRes <- pdsep_waves(skel@graph, suffStat, indepTest, p = p, 
                          sepset_sk, alpha, skel@pMax, waves, m.max = 2)
  G <- pdsepRes$G
  sepset <- pdsepRes$sepset
  time_pdsep = proc.time() - time_start
  # resIP <- list(G = G, G_sk = G_sk, sepset = sepset, sepset_sk = sepset_sk, 
  #               time_skel = time_skel, time_pdsep = time_pdsep, skel = skel, 
  #               pdsepRes = pdsepRes)
  resIP <- list(G = G, sepset = sepset, G_sk = G_sk, sepset_sk = sepset_sk,
                time_pdsep=time_pdsep)
  return(resIP)
}