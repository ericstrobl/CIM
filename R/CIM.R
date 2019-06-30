#' Causal Inference over Mixtures (CIM) algorithm
#' @param suffStat Sufficient statistics for conditional independence test
#' @param indepTest Conditional independence test
#' @param alpha alpha level
#' @param p Total number of variables
#' @param waves List containing wave information
#' @param verbose
#' @param prior_know List containing any other prior knowledge (default: NULL)
#' @param resIP Precomputed skeleton if available (default: NULL)
#' @export

CIM <- function (suffStat, indepTest, alpha, p, waves,
                 verbose = FALSE, prior_know = NULL, resIP = NULL)
{
  if (verbose)
    cat("Compute Skeleton\n================\n")

  if (is.null(resIP)){
    resIP = skel_discovery_waves3(suffStat, indepTest, alpha=alpha,
                                              p=p, waves)
    G <- resIP$G_sk
    sepset <- resIP$sepset_sk
  } else {
    G <- resIP$G_sk
    sepset <- resIP$sepset_sk
  }

  tripleList <- NULL
  if (verbose)
    cat("\nDirect egdes:\n=============\n")

   G1=G

  for (w1 in seq_len(length(waves)-1)){
    for (w2 in (w1+1):length(waves)){
      G1[waves[[w1]],waves[[w2]]]=G[waves[[w1]],waves[[w2]]]*2;
    }
  }

  G2 = G;
  for (p in seq_len(length(prior_know))){
    G2[prior_know[[p]][[1]],] = G[prior_know[[p]][[1]],]*2;
    G2[prior_know[[p]][[1]],prior_know[[p]][[2]]]=G[prior_know[[p]][[1]],prior_know[[p]][[2]]];
  }
  G1 = pmax(G1,G2)

  sepset2 <- sepset_middle(G1, suffStat, indepTest,
                                 sepset, alpha, verbose = verbose)

  res <- udag2pag4_simple(pag = G1, sepset, sepset2, c(), rules = rules,
                          unfVect = NULL, verbose = verbose, rules_used = c())

  return(list(maag = res$pag, pre_OR_res = G1, rules_used = res$rules_used, sepset=sepset))
}
