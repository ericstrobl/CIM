matrix2 <- function(mat)
  
  if (is.matrix(mat)) {
    return(mat);
  } else {
    mat = matrix(mat);
    return(mat); }