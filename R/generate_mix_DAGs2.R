generate_mix_DAGs2 <- function(p, en, waves){
  
  N = p*p - p;
  samplesB = rbinom(N/2,1, en/(p-1) );
  graph = matrix(0,p,p)
  graph[upper.tri(graph, diag=FALSE)] <- samplesB;
  
  for (w1 in 1:(length(waves)-1)){
    graph[waves[[w1]],waves[[w1+1]]]=graph[waves[[w1]],waves[[w1+1]]] + diag(length(waves[[1]]))
  }
  graph = graph>0
  
  nS = sample(0:2,1);
  nL = sample(0:2,1);
  # nS = 0
  # nL = 0
  
  pL = which(rowSums(graph)>=2); #variables with >=2 children
  L = sample(pL, min(length(pL),nL));
  actual_indices=setdiff(1:p,L)

  waves_L = waves;
  for (ll in seq_len(length(L))){
    wave_idx = which(sapply(waves_L, FUN=function(X) L[ll] %in% X))
    waves_L[[wave_idx]]=setdiff(waves_L[[wave_idx]],L[ll])
  }
  
  pS = setdiff(which(colSums(graph)>=2),L); #variables with >=2 parents that are not in L
  S = sample(pS, min(length(pS),nS));
  S_prob = runif(length(S),0.1,0.5);
  
  nT = sample(5:15,1)
  T_prob = runif(nT,0.3,0.7);
  # T_prob=1
  
  weights = matrix((0.75*runif(p^2)+0.25)*sample(c(-1,1),p^2,replace=TRUE),p,p)
  graph = graph*weights;
  
  is = which(graph!=0)
  DAGs = approx_block_random(is,g=length(T_prob))
  # DAGs = c()
  
  return(list(graph=graph,T_prob=T_prob,S_prob=S_prob,waves_L=waves_L,waves=waves,S=S,L=L,DAGs=DAGs))
  
}

