get_all_cycles <- function(adj){
  
  p=nrow(adj);
  list=list();
  for (v in 1:p){
    list = get_all_cycles_per_variable(adj, v, list)
  }
  
  list_n = list;
  for (l in seq_len(length(list))){
    list_n[[l]]=sort(list[[l]]);
  }
  
  list=list[which(!duplicated(list_n))]
  
  return(list)
  
  
}