udag2pdag_new <- function (G_sk, sepset_sk, waves, prior_know, verbose = FALSE) 
{

  if (sum(G_sk) > 0) {
    g <- as(G_sk, "matrix")
    
    gp = matrix(0,ncol(g),ncol(g));
    for (w1 in seq_len(length(waves)-1)){
      for (w2 in (w1+1):length(waves)){
        gp[waves[[w1]],waves[[w2]]] = 1;
      }
    }
    for (p in seq_len(length(prior_know))){
      gp[prior_know[[p]][[1]],]=1;
      gp[prior_know[[p]][[1]],prior_know[[p]][[2]]]=1;
    }
    # gp = gp*-1 + 1; # can be ancestors
    # print(gp)
    
    p <- as.numeric(dim(g)[1])
    pdag <- g
    ind <- which(g == 1, arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      x <- ind[i, 1]
      y <- ind[i, 2]
      allZ <- setdiff(which(g[y, ] == 1), x)
      for (z in allZ) {
        if (g[x, z] == 0 && !(y %in% sepset_sk[[x]][[z]] || 
                              y %in% sepset_sk[[z]][[x]])) {
          
          if (verbose) {
            cat("\n", x, "->", y, "<-", z, "\n")
            cat("Sxz=", sepset_sk[[z]][[x]], "Szx=", 
                sepset_sk[[x]][[z]])
          }

          if (gp[y,x] != 1 & gp[y,z]!=1){
            pdag[x, y] <- pdag[z, y] <- 1
            pdag[y, x] <- pdag[y, z] <- 0
          }
        }
      }
    }
    
    res2 <- pdag2dag(as(pdag, "graphNEL"))
    if (res2$success) {
      old_pdag <- matrix(0, p, p)
      while (!all(old_pdag == pdag)) {
        old_pdag <- pdag
        ind <- which((pdag == 1 & t(pdag) == 0), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          b <- ind[i, 2]
          indC <- which((pdag[b, ] == 1 & pdag[, b] == 
                           1) & (pdag[a, ] == 0 & pdag[, a] == 0))
          if (length(indC) > 0) {
            for (c in 1:length(indC)){
              if(gp[c,b]!=1){
                pdag[b, c] <- 1
                pdag[c, b] <- 0
              }
            }
            if (verbose) 
              cat("\nRule 1:", a, "->", b, " and ", b, 
                  "-", indC, " where ", a, " and ", indC, 
                  " not connected: ", b, "->", indC, "\n")
          }
        }
        ind <- which((pdag == 1 & t(pdag) == 1), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          b <- ind[i, 2]
          indC <- which((pdag[a, ] == 1 & pdag[, a] == 
                           0) & (pdag[, b] == 1 & pdag[b, ] == 0))
          if (length(indC) > 0) {
            
            if(gp[b,a]!=1){
              pdag[a, b] <- 1
              pdag[b, a] <- 0
            }
            
            if (verbose) 
              cat("\nRule 2: Kette ", a, "->", indC, 
                  "->", b, ":", a, "->", b, "\n")
          }
        }
        ind <- which((pdag == 1 & t(pdag) == 1), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          b <- ind[i, 2]
          indC <- which((pdag[a, ] == 1 & pdag[, a] == 
                           1) & (pdag[, b] == 1 & pdag[b, ] == 0))
          if (length(indC) >= 2) {
            g2 <- pdag[indC, indC]
            if (length(g2) <= 1) {
              g2 <- 0
            }
            else {
              diag(g2) <- rep(1, length(indC))
            }
            if (any(g2 == 0)) {
              if(gp[b,a]!=1){
                pdag[a, b] <- 1
                pdag[b, a] <- 0
              }
              if (verbose) 
                cat("\nRule 3:", a, "->", b, "\n")
            }
          }
        }
      }
      graph <- pdag
      
    }
    else {
      graph <- res2$graph
      graph <- dag2cpdag(as(pdag, "graphNEL"))
    }
  }
  return(as(graph,"matrix"))
}