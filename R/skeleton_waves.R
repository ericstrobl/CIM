#' Skeleton discovery with longitudinal data


skeleton_waves <- function (suffStat, indepTest, alpha, waves, labels = NULL, p, m.max = Inf) 
{

  seq_p <- seq_len(p)
  
  G <- matrix(TRUE, nrow = p, ncol = p)
  diag(G) <- FALSE
  
  # if (length(waves)>=2){
  #  for (w1 in seq_len(length(waves)-2)){
  #    for (w2 in (w1+2):length(waves)){
  #      G[waves[[w1]],waves[[w2]]]=FALSE;
  #      G[waves[[w2]],waves[[w1]]]=FALSE;
  #    }
  #  }
  # }

  # G.w = matrix(FALSE, nrow = p, ncol = p);
  # for (w1 in seq_len(length(waves))){
  #     G.w[waves[[w1]],waves[[w1]]]=TRUE;
  # }
  pval <- NULL
  sepset <- lapply(seq_p, function(.) vector("list", p))
  pMax <- matrix(-Inf, nrow = p, ncol = p)
  diag(pMax) <- 1
  done <- FALSE
  ord <- 0L
  n.edgetests <- numeric(1)
  while (!done && any(G) && ord <= m.max) {
    n.edgetests[ord1 <- ord + 1L] <- 0
    done <- TRUE
    ind <- which(G, arr.ind = TRUE)
    ind <- ind[order(ind[, 1]), ]
    remEdges <- nrow(ind)

    G.l <- split(G, gl(p, p))
    
    for (i in 1:remEdges) {

      x <- ind[i, 1]
      y <- ind[i, 2]
      
      if (G[y, x]) {
        nbrsBool <- G.l[[x]] # get neighbors of x, but also in the wave of x or wave of y
        nbrsBool[y] <- FALSE # remove y from neighbors of x
        nbrs <- intersect(seq_p[nbrsBool],which_waves(x,y,waves)) ### convert neighbors to set of numbers
        length_nbrs <- length(nbrs) # find number of neighbors
        if (length_nbrs >= ord) {
          if (length_nbrs > ord) 
            done <- FALSE
          S <- seq_len(ord)
          repeat {
            n.edgetests[ord1] <- n.edgetests[ord1] + 
              1
            pval <- indepTest(x, y, nbrs[S], suffStat)
            
            # print(c(x,y,nbrs[S]))
            
            if (is.na(pval)) 
              pval <- as.numeric(NAdelete)
            if (pMax[x, y] < pval) 
              pMax[x, y] <- pval
            if (pval >= alpha) {
              G[x, y] <- G[y, x] <- FALSE
              sepset[[x]][[y]] <- sepset[[y]][[x]] <- nbrs[S]
              break
            }
            else {
              nextSet <- getNextSet(length_nbrs, ord, 
                                    S)
              if (nextSet$wasLast) 
                break
              S <- nextSet$nextSet
            }
          }
        }
      }
    }
    ord <- ord + 1L
  }
  for (i in 1:(p - 1)) {
    for (j in 2:p) pMax[i, j] <- pMax[j, i] <- max(pMax[i, 
                                                        j], pMax[j, i])
  }
  
  Gobject <- if (sum(G) == 0) {
    new("graphNEL", nodes = labels)
  }
  else {
    colnames(G) <- rownames(G) <- labels
    as(G, "graphNEL")
  }
  cl <- match.call()
  new("pcAlgo", graph = Gobject, call = cl, n = integer(0), 
      max.ord = as.integer(ord - 1), n.edgetests = n.edgetests, 
      sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1))
}