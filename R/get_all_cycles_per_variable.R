get_all_cycles_per_variable <- function(adj, v, list){

  pa = which(adj[,v]);
  for (p in pa){
    cycles = all_simple_paths(graph_from_adjacency_matrix(adj),v,p, mode="out");
    
    for (c in seq_len(length(cycles))){
      cycles[[c]]=as.vector(cycles[[c]]);
      #cycles[[c]]=c(cycles[[c]], v);
    }
    
    list = c(list, cycles)
  }
  
  return(list)

}