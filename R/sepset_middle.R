sepset_middle <- function (pag, suffStat, indepTest, sepset, alpha, verbose = FALSE) 
{
  p = nrow(pag)
  allPdsep.tmp <- vector("list", p)
  sepset2 <- lapply(1:p, function(.) lapply(1:p, function(.) vector("list", 
                                                                       p)))
  allPdsep <- lapply(1:p, qreach, amat = pag)
  ind_a = which(pag == 2, arr.ind = TRUE)
  for (a in ind_a[, 1]) {
    ind_c = setdiff(ind_a[, 1], a)
    for (c in ind_c) {
      if (pag[a, c] == 0) {
        ind_b = which(pag[a, ] == 2 & pag[c, ] == 1)
        ind_b = setdiff(ind_b,sepset[[a]][[c]]) ###
        for (b in ind_b) {
          local_a = setdiff(allPdsep[[a]], c(a, b, c))
          if (length(sepset2[[a]][[b]][[c]]) == 0 & length(sepset[[a]][[c]])>0) {
            
            m=length(sepset[[a]][[c]])-1
            
            if (length(local_a) >= m) {
              if (length(local_a) > 1 | m == 0) {
                tmp.combn <- combn(local_a, m)
              }
              else if (length(local_a) == 1) {
                tmp.combn = matrix(local_a, 1, 1)
              }
              for (k in seq_len(ncol(tmp.combn))) {
                T <- tmp.combn[, k]
                T1 = union(T,b)
                pval <- indepTest(a, c, T1, suffStat)
                if (pval > alpha) {
                  sepset2[[a]][[b]][[c]] = T1
                  sepset2[[c]][[b]][[a]] = T1
                  break
                }
              }
            }

          }
        }
      }
    }
  }
  return(sepset2)
}