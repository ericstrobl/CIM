pdsep_waves <- function (skel, suffStat, indepTest, p, sepset, alpha, pMax, waves,
          m.max = Inf, pdsep.max = Inf, NAdelete = TRUE, unfVect = NULL, 
          biCC = FALSE, verbose = FALSE) 
{
  G <- (as(skel, "matrix") != 0)
  n.edgetests <- rep(0, 1000)
  ord <- 0L
  allPdsep.tmp <- vector("list", p)
  if (biCC) 
    conn.comp <- lapply(biConnComp(skel), as.numeric)
  if (any(G)) {
    amat <- G
    ind <- which(G, arr.ind = TRUE)
    storage.mode(amat) <- "integer"
    if (verbose) 
      cat("\nCompute collider:\n")
    for (i in seq_len(nrow(ind))) {
      x <- ind[i, 1]
      y <- ind[i, 2]
      allZ <- setdiff(which(amat[y, ] != 0), x)
      for (z in allZ) {
        if (amat[x, z] == 0 && !(y %in% sepset[[x]][[z]] || 
                                 y %in% sepset[[z]][[x]])) {
          if (length(unfVect) == 0) {
            amat[x, y] <- amat[z, y] <- 2
            if (verbose) 
              cat("\n", x, "*->", y, "<-*", z, "\n")
          }
          else {
            if (!any(unfVect == triple2numb(p, x, y, 
                                            z), na.rm = TRUE) && !any(unfVect == triple2numb(p, 
                                                                                             z, y, x), na.rm = TRUE)) {
              amat[x, y] <- amat[z, y] <- 2
              if (verbose) 
                cat("\n", x, "*->", y, "<-*", z, "\n")
            }
          }
        }
      }
    }
    allPdsep <- lapply(1:p, qreach, amat = amat)
    allPdsep.tmp <- vector("list", p)
    
    # G.w = matrix(TRUE, nrow = p, ncol = p);
    # for (w1 in seq_len(length(waves)-1)){
    #   for (w2 in (w1+1):length(waves)){
    #     G.w[waves[[w1]],waves[[w2]]]=FALSE;
    #   }
    # }
    # G.w = matrix(FALSE, nrow = p, ncol = p);
    # for (w1 in seq_len(length(waves))){
    #   G.w[waves[[w1]],waves[[w1]]]=TRUE;
    # }
    
    for (x in 1:p) {
      if (verbose) 
        cat("\nPossible D-Sep of", x, "is:", allPdsep[[x]], 
            "\n")
      if (any(an0 <- amat[x, ] != 0)) {
        # allPdsep[[x]] = intersect(allPdsep[[x]])
        tf1 <- setdiff(allPdsep[[x]], x)
        adj.x <- which(an0)
        for (y in adj.x) {
          if (verbose) 
            cat(sprintf("\ny = %3d\n.........\n", y))
          tf <- setdiff(tf1, y)
          diff.set <- setdiff(tf, adj.x)
          diff.set = intersect(diff.set, which_waves(x,y,waves))####
         
          allPdsep.tmp[[x]] <- c(tf, y)
          
          if (length(diff.set) > 0) {
            done <- FALSE
            ord <- 0L
            while (!done && ord < min(length(tf), m.max)) {
              ord <- ord + 1L
              if (verbose) 
                cat("ord = ", ord, "\n")
              if (ord == 1) {
                for (S in diff.set) {
                  pval <- indepTest(x, y, S, suffStat)
                  n.edgetests[ord + 1] <- n.edgetests[ord + 
                                                        1] + 1
                  if (is.na(pval)) 
                    pval <- as.numeric(NAdelete)
                  if (pval > pMax[x, y]) 
                    pMax[x, y] <- pval
                  if (pval >= alpha) {
                    amat[x, y] <- amat[y, x] <- 0
                    sepset[[x]][[y]] <- sepset[[y]][[x]] <- S
                    done <- TRUE
                    if (verbose) 
                      cat("x=", x, " y=", y, " S=", S, 
                          ": pval =", pval, "\n")
                    break
                  }
                }
              }
              else {
                tmp.combn <- combn(tf, ord)
                if (ord <= length(adj.x)) {
                  for (k in seq_len(ncol(tmp.combn))) {
                    S <- tmp.combn[, k]
                    if (!all(S %in% adj.x)) {
                      n.edgetests[ord + 1] <- n.edgetests[ord + 
                                                            1] + 1
                      pval <- indepTest(x, y, S, suffStat)
                      if (is.na(pval)) 
                        pval <- as.numeric(NAdelete)
                      if (pMax[x, y] < pval) 
                        pMax[x, y] <- pval
                      if (pval >= alpha) {
                        amat[x, y] <- amat[y, x] <- 0
                        sepset[[x]][[y]] <- sepset[[y]][[x]] <- S
                        done <- TRUE
                        if (verbose) 
                          cat("x=", x, " y=", y, " S=", 
                              S, ": pval =", pval, "\n")
                        break
                      }
                    }
                  }
                }
                else {
                  for (k in seq_len(ncol(tmp.combn))) {
                    S <- tmp.combn[, k]
                    n.edgetests[ord + 1] <- n.edgetests[ord + 
                                                          1] + 1
                    pval <- indepTest(x, y, S, suffStat)
                    if (is.na(pval)) 
                      pval <- as.numeric(NAdelete)
                    if (pMax[x, y] < pval) 
                      pMax[x, y] <- pval
                    if (pval >= alpha) {
                      amat[x, y] <- amat[y, x] <- 0
                      sepset[[x]][[y]] <- sepset[[y]][[x]] <- S
                      done <- TRUE
                      if (verbose) 
                        cat("x=", x, " y=", y, " S=", 
                            S, ": pval =", pval, "\n")
                      break
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    G[amat == 0] <- FALSE
    G[amat == 1] <- TRUE
    G[amat == 2] <- TRUE
  }
  list(G = G, sepset = sepset, pMax = pMax, allPdsep = allPdsep.tmp, 
       max.ord = ord, n.edgetests = n.edgetests[1:(ord + 1)])
}