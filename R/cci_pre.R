cci_pre <- function (suffStat, indepTest, alpha, p, skeleton_pre = NULL, 
          rules = rep(TRUE, 7), verbose = FALSE) 
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
  tripleList <- NULL
  if (verbose) 
    cat("\nDirect egdes:\n=============\n")
  G <- CCI:::v_struc(pag = G, sepset, unfVect = tripleList, verbose)
  list_pre_v <- CCI:::after_v_struc(pag = G, sepset, suffStat, indepTest, 
                              alpha, verbose = verbose)

  sup_sepset <- extra_dsep_waves(list_pre_v$G, suffStat, indepTest, 
                           sepset, alpha, waves, verbose = verbose)
  list_e <- CCI:::step_5(list_pre_v$G, sepset, sup_sepset, suffStat, 
                   indepTest, alpha, verbose = verbose, rules_used = list_pre_v$rules_used)
  list_f <- CCI:::step_6(list_e$pag, sepset, sup_sepset, suffStat, 
                   indepTest, alpha, verbose = verbose, list_e$rules_used)
  res <- CCI:::udag2pag4(pag = list_f$pag, sepset, rules = rules, 
                   unfVect = tripleList, verbose = verbose, rules_used = list_f$rules_used)
  return(list(maag = res$pag, pre_OR_res = list_f$pag, rules_used = res$rules_used, 
              sepset = sepset))
}