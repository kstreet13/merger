mergeManyPairwise <- function(clusteringMatrix, nCores = 3) {
  # Turn the matrix into a numeric matrix
  clusMat <- apply(clusteringMatrix, 2, function(x) {
    x[x != "-1"] <- as.numeric(factor(x[x != "-1"]))
    x[x == "-1"] <- -1
    x <- as.integer(x)
  })
  
  # Initialize the values
  clusters <- apply(clusMat, 2, unique)
  currentMat <- clusMat
  baseARI <- apply(clusMat, 2, function(x) {
    apply(clusMat, 2, function(y) {
      mclust::adjustedRandIndex(x, y)
    })
  })
  bestARI <- baseARI
  working <- TRUE
  merges <- NULL
  
  # Try to see if any merge would increse 
  while (working) {
    # Test all pairwise clusters to merge
    # For every cluster label list
    mergeResults <- mclapply(1:ncol(currentMat), function(whClus) {
      clus <- currentMat[, whClus]
      clusternames <- clusters[[whClus]]
      clusPairs <- combn(clusternames[clusternames != -1], 2)
      
      # For every pair of labels in that list
      deltaARI <- apply(clusPairs, 2, function(pair) {
        sapply((1:ncol(clusMat))[-whClus], function(otherClus) {
          # This will modify clus?
          clus[clus %in% pair] <- max(clus) + 1
          # Shouldn't it be currentMat here
          mclust::adjustedRandIndex(clus, clusMat[, otherClus])
        })
      }) - bestARI[whClus, -whClus]
      
      return(colMeans(deltaARI))
    }, mc.cores = nCores)
    
    # Find best pair to merge
    maxs <- sapply(mergeResults, max)
    
    # Only merge if it improves ARI
    if (max(maxs) > 0) {
      whClus <- which.max(maxs)
      # update clusters
      clusternames <- clusters[[whClus]]
      clusPairs <- combn(clusternames[clusternames != -1], 2)
      pair <- clusPairs[, which.max(mergeResults[[whClus]])]
      indsPair <- which(currentMat[, whClus] %in% pair)
      currentMat[indsPair, whClus] <- min(pair)
      clusters[[whClus]] <- unique(currentMat[, whClus])
      
      # update bestARI
      newARIs <- sapply((1:ncol(currentMat))[-whClus], function(j) {
        mclust::adjustedRandIndex(currentMat[, whClus], currentMat[, j])
      })
      bestARI[whClus, -whClus] <- newARIs
      bestARI[-whClus, whClus] <- newARIs
      # tracking
      merges <- rbind(merges, c(whClus, pair))
      print(c(whClus, pair))
    } else {
      working <- FALSE
    }
    
    # If no more to merge in any of them, stop
    if (sum(sapply(clusters, length) == 1) == length(clusters)) stop()
  }
  
  colnames(merges) <- c("clustering", "cluster1", "cluster2")
  return(list("currentMat" = currentMat,
              "merges" = merges))
}
