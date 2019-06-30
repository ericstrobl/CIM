generate_mix_DAGs <- function(nIndep, p, en, waves){
  
  enn=en/(2*nIndep);
  N = p*p - p;
  graphs=vector("list",nIndep)
  graphs_sum = matrix(0,p,p)
  for (In in 1:(2*nIndep)){
    samplesB = rbinom(N,1, enn/(p-1) );
    
    graph = matrix(0,p,p);
    pt = sample(1:p,p,replace=FALSE)
    graphs[[In]]=graph;
     graphs[[In]][upper.tri(graph, diag=FALSE)] <- samplesB[1:(N/2)];
     # graphs[[In]][lower.tri(graph, diag=FALSE)] <- samplesB[(N/2+1):N];
    graphs[[In]] = graphs[[In]][pt,pt]
    
    for (w2 in length(waves):2){
      for (w1 in (w2-1):1){
        graphs[[In]][waves[[w2]],waves[[w1]]]=0
      }
    }
    
    for (w1 in 1:(length(waves)-1)){
      # dd = rbinom(length(waves[[1]]),1,0.5) ###
      if (In %% 2 == 0 & w1 == 2){
        graphs[[In]][waves[[w1]],waves[[w1+1]]]=graphs[[In]][waves[[w1]],waves[[w1+1]]]+diag(length(waves[[1]]))
      } else if (In %% 2 != 0 & w1 == 1){
        graphs[[In]][waves[[w1]],waves[[w1+1]]]=graphs[[In]][waves[[w1]],waves[[w1+1]]]+diag(length(waves[[1]]))
      }
    }
    graphs[[In]][which(graphs[[In]]>1)]=1;
    
    graphs_sum = graphs_sum + graphs[[In]]
    
  }
  graphs_sum = graphs_sum>0;
  
  mDAGs=list();
  mDAGs$father_graph = graphs_sum
  mDAGs$graphs=graphs;
  mDAGs$weights = matrix((0.75*runif(p^2)+0.25)*sample(c(-1,1),p^2,replace=TRUE),p,p)
  
  nS = sample(0:2,1);
  nL = sample(0:2,1);

  pL = which(rowSums(graphs_sum)>=2); #variables with >=2 children
  mDAGs$L = sample(pL, min(length(pL),nL));
  mDAGs$actual_indices=setdiff(1:p,mDAGs$L)
  
  for (L in seq_len(length(mDAGs$L))){
    wave_idx = which(sapply(waves, FUN=function(X) mDAGs$L[L] %in% X))
    waves[[wave_idx]]=setdiff(waves[[wave_idx]],mDAGs$L[L])
  }
  mDAGs$waves=waves;
  
  pS = setdiff(which(colSums(graphs_sum)>=2),mDAGs$L); #variables with >=2 parents that are not in L
  mDAGs$S = sample(pS, min(length(pS),nS));
  
  graph = matrix(0,p+length(mDAGs$S), p+length(mDAGs$S));
  graph[1:p,1:p] = graphs_sum;
  
  for (s in seq_len(length(mDAGs$S))){
    graph[mDAGs$S[s],p+s]=1;
  }
  
  mDAGs$graph=graph;
  
  S_prob = runif(length(mDAGs$S),0.1,0.5);
  
  if (length(mDAGs$S)>0){
    mDAGs$elim_prop = runif(1,0.1,0.5);
  } else{
    mDAGs$elim_prop = 0;
  }
  
  mDAGs$S_prob = mDAGs$elim_prop*(S_prob / sum(S_prob));
  
  mDAGs$indep_prob = runif(nIndep)
  
  # weights = runif(ncol(mDAGs$indep_prob))
  # mDAGs$weights = weights/sum(weights);
  
  return(mDAGs)
  
}

