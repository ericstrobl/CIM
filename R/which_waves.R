which_waves <- function(x,y,waves){
  
  xi = which(sapply(waves, FUN=function(X) x %in% X))
  
  if (!is.null(y)){
    yi = which(sapply(waves, FUN=function(X) y %in% X))
    
    wmax = max(xi,yi)
    wmin = min(xi,yi)
    
    is = c()
    for (w in wmin:wmax){
      is = c(is,waves[[w]])
    }
  } else{
    is = waves[[xi]]
  }
  is = setdiff(is,c(x,y))
  
  
  # ##############
  # for (w in 1:length(waves)){
  #   is = c(is,waves[[w]])
  # }
  # is = setdiff(is,c(x,y))
  # ###########
  
  return(is)
}