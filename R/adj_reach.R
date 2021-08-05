adj_reach <- function (x, amat, verbose = FALSE) 
{
  stopifnot((ncol(amat) == nrow(amat)), x <= ncol(amat), all(amat %in% 
                                                               c(0, 1, 2)), all((amat != 0) == (t(amat != 0))))
  A. <- (amat != 0)
  PSEP <- which(A.[x, ])
  sort(setdiff(PSEP, x))
}