approx_block_random <- function(is,g){
  
  assign = vector(mode = "list", length = g)
  
  rem = is
  for (i in seq_len(length(is))){
     ls = lengths(assign)
     lmin = which(ls == min(ls))[1]
    
     if (length(rem)>1){
       pick = sample(rem,1,replace=FALSE);
       assign[[lmin]] = c(assign[[lmin]], pick)
     } else{
       assign[[lmin]] = c(assign[[lmin]], rem)
     }
     
     rem = rem[-which(rem==pick)]
     
  }
  
  order = sample(1:g,g,replace=FALSE)
  assign  = assign[order]
  
  return(assign)
  
}