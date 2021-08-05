sample_mix_DAGs2<-function(mixDAG, nsamps){
  
  nsamps = nsamps*6;
  
  p=nrow(mixDAG$graph)
  data = matrix(0,nsamps,p)
  
  n = 0;
  for (s in 1:nsamps){
    
    if (length(mixDAG$T_prob)>1){
      Ts = c()
      for (t in 1:length(mixDAG$T_prob)){
         Ts = c(Ts, sample(0:1,1,prob=c(1-mixDAG$T_prob[t],mixDAG$T_prob[t])))
      }
      
      Tp=which(Ts==0)
      for (tp in Tp){
        A = mixDAG$graph;
        A[mixDAG$DAGs[[tp]]] = 0;
      }
      
    } else{
      A = mixDAG$graph;
    }
    Cov = ginv(diag(nrow(A)) - t(A))
    data[s,] = mvrnorm(1, rep(0,p), Cov %*% t(Cov));
    
  }
  nsamps = nsamps/6
  
  sample_cut = c();
  for (s in seq_len(length(mixDAG$S))){
    s_cutoff = quantile(data[,mixDAG$S[s]], mixDAG$S_prob[s]);
    sample_cut = c(sample_cut, which(data[,mixDAG$S[s]]<=s_cutoff))
  }
  
  if (length(sample_cut)>0){
    data = data[-sample_cut,];
  }
  
  # if (length(c(mixDAG$L))>0){
  #   data=data[,-c(mixDAG$L)];
  # }
  
  samps = sample(nrow(data),nsamps, replace=FALSE)
  data = data[samps,,drop=FALSE]
  
  # colnames(data) = setdiff(1:ncol(mixDAG$graph), mixDAG$L);
  colnames(data) = 1:ncol(mixDAG$graph)
  
  return(data)
  
  
}