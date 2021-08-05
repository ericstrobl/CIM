fci_pre <- function (suffStat, indepTest, alpha, p, skeleton_pre=NULL, skel.method = c("stable",
                                                                                            "original", "stable.fast"), type = c("normal", "anytime",
                                                                                                                                 "adaptive"), fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE,
                  m.max = Inf, pdsep.max = Inf, rules = rep(TRUE, 10), doPdsep = TRUE,
                  biCC = FALSE, conservative = FALSE, maj.rule = FALSE, verbose = FALSE)
{

  if (verbose)
    cat("Compute Skeleton\n================\n")
  
  if (is.null(skeleton_pre)){
    pdsepRes = IP_discovery(suffStat, indepTest = indepTest, 
                            alpha = alpha, p = p)
    G <- pdsepRes$G
    sepset <- pdsepRes$sepset
  } else {
    G <- skeleton_pre$G
    sepset <- skeleton_pre$sepset
  }

  if (verbose)
    cat("\nDirect egdes:\n=============\nUsing rules:", which(rules),
        "\nCompute collider:\n")
  res <- udag2pag(pag = G, sepset, rules = rules, verbose = verbose)
  colnames(res) <- rownames(res) <- 1:p
  return(list(maag = res))
}
