udag2pag4_simple <- function (pag, indepTest, sepset, sepset2, sepset_e, rules = rep(TRUE, 8), unfVect = NULL, 
                           verbose = FALSE, rules_used = c(), alpha = 0.01) 
{
  rules = rep(FALSE,8)
  rules[c(1,5)]=TRUE
  stopifnot(is.logical(rules), length(rules) == 8)
  if (any(pag != 0)) {
    p <- as.numeric(dim(pag)[1])
    old_pag1 <- matrix(0, p, p)
    while (any(old_pag1 != pag)) {
      old_pag1 <- pag
      if (rules[1]) {
        ind <- which((pag == 2 & t(pag) != 0), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          b <- ind[i, 2]
          indC <- which((pag[b, ] != 0 & pag[, b] == 
                           1) & (pag[a, ] == 0 & pag[, a] == 0))
          indC <- setdiff(indC, a)

          for (c in indC) {
            if (b %in% sepset[[a]][[c]]){
              # pval1 <- indepTest(a, b, setdiff(sepset[[a]][[c]],b), suffStat)
              # pval2 <- indepTest(c, b, setdiff(sepset[[a]][[c]],b), suffStat)
              # if ((pval1 < alpha) & (pval2 < alpha)) {
              #   pag[c, b] <- 3
              # }
              pag[c, b] <- 3
              rules_used = unique(c(rules_used, 1))
            } else if (length(sepset2[[a]][[b]][[c]])>0){
              # pval1 <- indepTest(a, b, setdiff(sepset2[[a]][[b]][[c]],b), suffStat)
              # pval2 <- indepTest(c, b, setdiff(sepset2[[a]][[b]][[c]],b), suffStat)
              # if ((pval1 < alpha) & (pval2 < alpha)) {
              #   pag[c, b] <- 3
              # }
              pag[c, b] <- 3
              rules_used = unique(c(rules_used, 1))
            }  else if (length(sepset2[[c]][[b]][[a]])>0){
              # pval1 <- indepTest(a, b, setdiff(sepset2[[c]][[b]][[a]],b), suffStat)
              # pval2 <- indepTest(c, b, setdiff(sepset2[[c]][[b]][[a]],b), suffStat)
              # if ((pval1 < alpha) & (pval2 < alpha)) {
              #   pag[c, b] <- 3
              # }
              pag[c, b] <- 3
              rules_used = unique(c(rules_used, 1))
            }
            if (verbose) 
              cat("\nRule 1a", "\nOrient:", a, "*->", 
                  b, "o-*", c, "as:", b, "-*", c, 
                  "\n")
          }
        }
      }
      if (rules[2]) {
        ind <- which((pag != 0 & t(pag) == 1), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          b <- ind[i, 1]
          c <- ind[i, 2]
          indA <- which((pag[b, ] == 3 & pag[, b] != 
                           0) & (pag[c, ] == 0 & pag[, c] == 0))
          indA <- setdiff(indA, c)
          if (length(indA) > 0) {
            if (length(select_not_2triangle(pag, c, b, 
                                            indA, 3)) > 0) {
              if (length(unfVect) == 0) {
                pag[c, b] <- 3
                rules_used = unique(c(rules_used, 2))
                if (verbose) 
                  cat("\nRule 2", "\nOrient:", indA, 
                      "-o", b, "o-*", c, "as", b, "-*", 
                      c, "\n")
              }
            }
          }
        }
      }
      if (rules[3]) {
        ind <- which((pag == 2 & t(pag) != 0), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          b <- ind[i, 2]
          indC <- which(pag[b, ] == 3 & pag[, b] == 3 & 
                          pag[a, ] == 0 & pag[, a] == 0, arr.ind = TRUE)
          for (c in seq_len(length(indC))) {
            indDf = count_2triangle(pag, a, b, indC[c])
            if (length(indDf) == 1) {
              if (pag[a, indDf[1]] == 1) {
                pag[a, indDf[1]] = 2
                rules_used = unique(c(rules_used, 31))
                if (verbose) {
                  cat("\nRule 3a", "\nOrient:", a, "*-o", 
                      indDf[1], "as", a, "*->", indDf[1], 
                      "\n")
                }
              }
              if (pag[b, indDf[1]] == 1) {
                pag[b, indDf[1]] = 3
                rules_used = unique(c(rules_used, 32))
                if (verbose) {
                  cat("\nRule 3b", "\nOrient:", b, "*-o", 
                      indDf[1], "as", b, "*-", indDf[1], 
                      "\n")
                }
              }
              if (pag[indDf[1], b] == 1) {
                pag[indDf[1], b] = 3
                rules_used = unique(c(rules_used, 33))
                if (verbose) {
                  cat("\nRule 3c", "\nOrient:", indDf[1], 
                      "*-o", b, "as", indDf[1], "*-", b, 
                      "\n")
                }
              }
              p_und <- is_one_undirected_path(pag, indDf[1], 
                                              indC[c], b)
              if (length(p_und) == 1) {
                for (j in 2:length(p_und[[1]])) {
                  if (pag[p_und[[1]][j], p_und[[1]][j - 
                                                    1]] == 1) {
                    pag[p_und[[1]][j], p_und[[1]][j - 
                                                    1]] = 3
                    rules_used = unique(c(rules_used, 
                                          34))
                    if (verbose) 
                      cat("\nRule 3d", "\nOrient:", p_und[[1]][j], 
                          "*-o", p_und[[1]][j - 1], "as", 
                          p_und[[1]][j], "*-", p_und[[1]][j - 
                                                            1], "\n")
                  }
                  if (pag[p_und[[1]][j - 1], p_und[[1]][j]] == 
                      1) {
                    pag[p_und[[1]][j - 1], p_und[[1]][j]] = 3
                    rules_used = unique(c(rules_used, 
                                          34))
                    if (verbose) 
                      cat("\nRule 3d", "\nOrient:", p_und[[1]][j - 
                                                                 1], "*-o", p_und[[1]][j], "as", 
                          p_und[[1]][j - 1], "*-", p_und[[1]][j], 
                          "\n")
                  }
                }
              }
            }
          }
        }
      }
      if (rules[4]) {
        ind <- which((pag == 2 & t(pag) != 0), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          c <- ind[i, 2]
          pag1 = pag
          pag1[a, c] = 0
          pag1[c, a] = 0
          tail_paths <- is_tail_path_one_unknown(pag1, 
                                                 c, a)
          if (length(tail_paths) == 1) {
            for (j in seq_len(length(tail_paths[[1]]) - 
                              1)) {
              if (pag[tail_paths[[1]][j + 1], tail_paths[[1]][j]] == 
                  1) {
                pag[tail_paths[[1]][j + 1], tail_paths[[1]][j]] <- 2
                rules_used = unique(c(rules_used, 4))
                if (verbose) {
                  cat("\nRule 4", "\n")
                  cat("Orient:", tail_paths[[1]][j + 
                                                   1], "*-o", tail_paths[[1]][j], "as:", 
                      tail_paths[[1]][j + 1], "*->", tail_paths[[1]][j], 
                      "\n")
                }
              }
            }
          }
        }
      }
      if (rules[5]) {
        ind <- which((pag != 0 & t(pag) == 1), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          c <- ind[i, 2]
          if (is_dir_undirected(pag, a, c)) {
            pag[c, a] <- 3
            rules_used = unique(c(rules_used, 6))
            if (verbose) 
              cat("\nRule 5", "\nOrient:", a, "o-*", 
                  c, "as", a, "-*", c, "\n")
          }
        }
      }
      if (rules[6]) {
        ind <- which((pag != 0 & t(pag) == 1), arr.ind = TRUE)
        while (length(ind) > 0) {
          a <- ind[1, 1]
          c <- ind[1, 2]
          ind <- ind[-1, , drop = FALSE]
          indB <- which((pag[a, ] != 0) & (pag[, a] != 
                                             0) & (pag[c, ] == 0 & pag[, c] == 0))
          indB <- setdiff(indB, c)
          indB <- select_not_2triangle(pag, c, a, indB, 
                                       3)
          while ((length(indB) > 0) && (pag[c, a] == 
                                        1)) {
            b <- indB[1]
            indB <- indB[-1]
            upd <- minUncovPdPath3(p, pag, a, b, c, unfVect = unfVect, 
                                   verbose = verbose)
            if (length(upd) > 1) {
              pag[c, a] <- 3
              rules_used = unique(c(rules_used, 7))
              if (verbose) 
                cat("\nRule 6", "\nThere exists an uncovered potentially directed path between", 
                    a, "and", c, ". Orient:", a, " -*", 
                    c, "\n")
            }
          }
        }
      }
      if (rules[7]) {
        ind <- which((pag != 0 & t(pag) == 1), arr.ind = TRUE)
        while (length(ind) > 0) {
          a <- ind[1, 1]
          c <- ind[1, 2]
          ind <- ind[-1, , drop = FALSE]
          indB <- which((pag[c, ] == 3 & pag[, c] != 
                           0))
          if (length(indB) >= 2) {
            counterB <- 0
            while (counterB < length(indB) && (pag[c, 
                                                   a] == 1)) {
              counterB <- counterB + 1
              b <- indB[counterB]
              indD <- setdiff(indB, b)
              counterD <- 0
              while ((counterD < length(indD)) && (pag[c, 
                                                       a] == 1)) {
                counterD <- counterD + 1
                d <- indD[counterD]
                if ((pag[a, b] != 0) && (pag[b, a] != 
                                         0) && (pag[a, d] != 0) && (pag[d, a] != 
                                                                    0) && pag[d, b] == 0 && pag[b, d] == 
                    0 && a %in% sepset[[b]][[d]]) {
                  if (length(unfVect) == 0) {
                    pag[c, a] <- 3
                    rules_used = unique(c(rules_used, 
                                          81))
                    if (verbose) 
                      cat("\nRule 7 [easy]", "\nOrient:", 
                          a, "-*", c, "\n")
                  }
                  else if (!any(unfVect == triple2numb(p, 
                                                       b, a, d), na.rm = TRUE) && !any(unfVect == 
                                                                                       triple2numb(p, d, a, b), na.rm = TRUE)) {
                    pag[c, a] <- 3
                    rules_used = unique(c(rules_used, 
                                          81))
                    if (verbose) 
                      cat("\nRule 7 [easy]", "\nConservatively orient:", 
                          a, "-*", c, "\n")
                  }
                }
                else {
                  indX <- which(pag[a, ] != 0 & pag[, 
                                                    a] != 0, arr.ind = TRUE)
                  indX <- setdiff(indX, c)
                  if (length(indX >= 2)) {
                    counterX1 <- 0
                    while (counterX1 < length(indX) && 
                           pag[c, a] == 1) {
                      counterX1 <- counterX1 + 1
                      first.pos <- indX[counterX1]
                      indX2 <- setdiff(indX, first.pos)
                      counterX2 <- 0
                      while (counterX2 < length(indX2) && 
                             pag[c, a] == 1) {
                        counterX2 <- counterX2 + 1
                        sec.pos <- indX2[counterX2]
                        if (pag[first.pos, sec.pos] == 
                            0 & !(pag[first.pos, a] == 
                                  2 & pag[sec.pos, a] == 2)) {
                          if (!is_2triangle(pag, c, a, 
                                            first.pos) & !is_2triangle(pag, 
                                                                       c, a, sec.pos)) {
                            t1 <- minUncovPdPath3(p, 
                                                  pag, a, first.pos, b, unfVect = unfVect, 
                                                  verbose = verbose)
                            if (length(t1) > 1) {
                              t2 <- minUncovPdPath3(p, 
                                                    pag, a, sec.pos, d, unfVect = unfVect, 
                                                    verbose = verbose)
                              if (length(t2) > 1 && first.pos != 
                                  sec.pos) {
                                if (length(unfVect) == 
                                    0) {
                                  pag[c, a] <- 3
                                  rules_used = unique(c(rules_used, 
                                                        82))
                                  if (verbose) 
                                    cat("\nRule 7", "\nOrient:", 
                                        a, "-*", c, "\n")
                                }
                                else if (!any(unfVect == 
                                              triple2numb(p, first.pos, 
                                                          a, sec.pos), na.rm = TRUE) && 
                                         !any(unfVect == triple2numb(p, 
                                                                     sec.pos, a, first.pos), 
                                              na.rm = TRUE)) {
                                  pag[c, a] <- 3
                                  rules_used = unique(c(rules_used, 
                                                        82))
                                  if (verbose) 
                                    cat("\nRule 7", "\nConservatively orient:", 
                                        a, "-*", c, "\n")
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      if (rules[8]) {
        ind <- which((pag == 3 & t(pag) != 0), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          b <- ind[i, 2]
          indC <- which((pag[b, ] != 0 & pag[, b] == 
                           1) & (pag[a, ] == 0 & pag[, a] == 0))
          indC <- setdiff(indC, a)
          
          for (c in indC) {
            if (!(b %in% sepset[[a]][[c]]) | length(sepset_e[[a]][[b]][[c]])>0 ){
              pag[c, b] <- 2
              rules_used = unique(c(rules_used, 1))
              if (verbose) 
                cat("\nRule 8", "\nOrient:", a, "*-", 
                    b, "o-*", c, "as:", b, "<-*", c, 
                    "\n")
            }
          }
        }
        
      }
    }
  }
  return(list(pag = pag, rules_used = rules_used))
}