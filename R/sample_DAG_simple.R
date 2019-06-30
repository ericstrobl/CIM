sample_DAG_simple <- function(nsamps, weights){
  
  A=weights;
  
  A=t(A);
  
  p=nrow(weights);
  
  Cov = ginv(diag(p) - A)
  
  data = mvrnorm(nsamps, rep(0,p), Cov %*% t(Cov));
  
  return(data)
}