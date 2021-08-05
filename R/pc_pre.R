pc_pre <- function (suffStat, indepTest, alpha, labels, p, skeleton_pre, 
                    waves, prior_know, verbose = FALSE) 
{
  if (verbose) 
    cat("Compute Skeleton\n================\n")
  
  if (is.null(skeleton_pre)){
    pdsepRes = IP_discovery(suffStat, indepTest = indepTest, 
                            alpha = alpha, p = p)
    G_sk <- pdsepRes$G_sk
    sepset_sk <- pdsepRes$sepset_sk
  } else {
    G_sk <- skeleton_pre$G_sk
    sepset_sk <- skeleton_pre$sepset_sk
  }
  
  res=udag2pdag_new(G_sk, sepset_sk, waves, prior_know, verbose)
  
  G_sk = matrix(0,p,p)
  idx = which(res>0,arr.ind=TRUE);
  for (i in seq_len(nrow(idx))){
   if (res[idx[i,1],idx[i,2]]==1 & res[idx[i,2],idx[i,1]]==0){
     G_sk[idx[i,1],idx[i,2]]=2;
     G_sk[idx[i,2],idx[i,1]]=3;
   } else if (res[idx[i,1],idx[i,2]]==1 & res[idx[i,2],idx[i,1]]==1){
     G_sk[idx[i,1],idx[i,2]]=1;
     G_sk[idx[i,2],idx[i,1]]=1;
   }
  }
  
  return(list(maag = G_sk))
}