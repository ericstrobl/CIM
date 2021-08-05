rfci_pre <- function (suffStat, indepTest, alpha, p, skel.method = c("stable", 
                                                                 "original", "stable.fast"), fixedGaps = NULL, fixedEdges = NULL, 
          NAdelete = TRUE, m.max = Inf, rules = rep(TRUE, 10), conservative = FALSE, 
          maj.rule = FALSE, numCores = 1, verbose = FALSE, skeleton_pre=NULL) 
{
  if (verbose) 
    cat("Compute Skeleton\n================\n")
  
  if (is.null(skeleton_pre)){
    pdsepRes = IP_discovery(suffStat, indepTest = indepTest, 
                            alpha = alpha, p = p)
    G_sk <- pdsepRes$G_sk
    sepset_sk <- pdsepRes$sepset_sk
  } else {
    G_sk <- skeleton_pre$G_sk
    sepset_sk <- skeleton_pre$sepset_sk
  }
  sk.A <- G_sk
  sepset <- sepset_sk
  u.t <- find.unsh.triple(sk.A, check = FALSE)
  r.v. <- rfci.vStruc(suffStat, indepTest, alpha, sepset, sk.A, 
                      unshTripl = u.t$unshTripl, unshVect = u.t$unshVect, conservative = (conservative || 
                                                                                            maj.rule), version.unf = c(1, 1), maj.rule = maj.rule, 
                      verbose = verbose)
  A <- r.v.$amat
  sepset <- r.v.$sepset
  if (verbose) 
    cat("\nDirect egdes:\n=============\nUsing rules:", which(rules), 
        "\n")
  res <- udag2apag(A, suffStat, indepTest, alpha, sepset, rules = rules, 
                   unfVect = r.v.$unfTripl, verbose = verbose)
  Amat <- res$graph
  colnames(Amat) <- rownames(Amat) <- 1:p
  return(list(maag=Amat))
}