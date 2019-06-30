sample_mix_DAGs <-function(mDAGs,nsamps_o){

  p=nrow(mDAGs$graphs[[1]])
  nsamps = ceiling(nsamps_o * 4)
  data = matrix(0,nsamps,p)

  cycle_father = isCyclic(mDAGs$father_graph>0)

  cnt=0;
  father_graph_approx = matrix(0,p,p)
  for (s in 1:nsamps){

      Cs = runif(length(mDAGs$graphs)/2)

      Cs = Cs < mDAGs$indep_prob

      graph_sample = matrix(0,p,p);
      for (c in 1:length(Cs)){
          if (Cs[c]){
            graph_sample = graph_sample + mDAGs$graphs[[2*c-1]]
          } else{
            graph_sample = graph_sample + mDAGs$graphs[[2*c]]
          }
      }
      if (cycle_father){
        if (isCyclic(graph_sample>0)){
          cnt=cnt+1;
          # print(cnt)
          graph_sample = make_acyclic(graph_sample>0)
        }
      }
      father_graph_approx = (father_graph_approx + graph_sample)>0;

      data[s,] = sample_DAG_simple(1,graph_sample)
  }
  mDAGs$father_graph_approx = father_graph_approx

  sample_cut = c();
  for (s in seq_len(length(mDAGs$S))){
    s_cutoff = quantile(data[,mDAGs$S[s]], mDAGs$S_prob[s]);
    sample_cut = c(sample_cut, which(data[,mDAGs$S[s]]<=s_cutoff))
    #print(length(which(data_d[,mDAGs$S[s]]<=s_cutoff)))
  }

  if (length(sample_cut)>0){
    data = data[-sample_cut,];
  }

  if (length(c(mDAGs$L))>0){
    data=data[,-c(mDAGs$L)];
  }

  samps = sample(nrow(data),nsamps_o, replace=FALSE)

  colnames(data) = setdiff(1:ncol(mDAGs$graphs[[1]]), mDAGs$L);

  return(list(data=data[samps,],mDAGs=mDAGs))


}
