make_acyclic <- function(DCG){
  
  cycles = get_all_cycles(DCG);
  for (c in seq_len(length(cycles))){
    # del_nodes = sample(seq_len(length(cycles[[c]])), 2, replace=FALSE);
    del_nodes = 1:2;
    cycle = c(cycles[[c]], cycles[[c]][1])
    del_edges = cbind( cycle[del_nodes], cycle[del_nodes+1] );
    DCG[del_edges[1,1],del_edges[1,2]]=0
  }
  
  return(DCG)
}
