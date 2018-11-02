################################################################################
### pairwise ARIs ###
################################################################################
require(clusterExperiment)
clusMat <- readRDS('data/clusMat.rds')
clusMat <- apply(clusMat,2,function(x){
    x[x != '-1'] <- as.numeric(factor(x[x != '-1']))
    x[x == '-1'] <- -1
    x <- as.integer(x)
}) # convert to integer
#coClus <- makeConsensus(clusMat, proportion = 0.8)$percentageShared
#denom <- (clusMat!=-1) %*% t(clusMat!=-1)
#num <- coClus * denom
#rm(coClus)
clusters <- apply(clusMat,2,unique)

#base.ARI <- mean(arimps)
#best.ARI <- base.ARI

require(parallel)

#mergeManyPairwise <- function(clusMat){
require(matrixStats)
n <- nrow(clusMat)

clusMat.m <- clusMat

base.ARI <- apply(clusMat,2,function(x){
    apply(clusMat,2,function(y){
        mclust::adjustedRandIndex(x,y)
    })
})
best.ARI <- base.ARI

working <- TRUE
merges <- NULL
while(working){
    # test all pairwise clusters to merge
    mergeResults <- mclapply(1:ncol(clusMat.m),function(wh.clus){
        clus <- clusMat.m[,wh.clus]
        clusternames <- clusters[[wh.clus]]
        clusPairs <- combn(clusternames[clusternames!=-1], 2)#[,1:20] #########
        
        deltaARI <- apply(clusPairs,2,function(pair){
            sapply((1:ncol(clusMat))[-wh.clus], function(j.clus){
                clus[clus %in% pair] <- max(clus)+1
                mclust::adjustedRandIndex(clus, clusMat[,j.clus])
            })
        }) - best.ARI[wh.clus, -wh.clus]
        
        return(colMeans(deltaARI))
    }, mc.cores = 3)
    
    # find best pair to merge, only merge if it improves ARI
    maxes <- sapply(mergeResults,max)
    if(max(maxes) > 0){
        wh.clus <- which.max(maxes)
        # update clusters
        clusternames <- clusters[[wh.clus]]
        clusPairs <- combn(clusternames[clusternames!=-1], 2)
        pair <- clusPairs[,which.max(mergeResults[[wh.clus]])]
        ind.ii <- which(clusMat.m[,wh.clus] == pair[1])
        ind.jj <- which(clusMat.m[,wh.clus] == pair[2])
        clusMat.m[c(ind.ii,ind.jj), wh.clus] <- min(pair)
        clusters[[wh.clus]] <- unique(clusMat.m[,wh.clus])
        # update best.ARI
        newARIs <- sapply((1:ncol(clusMat.m))[-wh.clus], function(j){
            mclust::adjustedRandIndex(clusMat.m[,wh.clus], clusMat.m[,j])
        })
        best.ARI[wh.clus,-wh.clus] <- newARIs
        best.ARI[-wh.clus,wh.clus] <- newARIs
        # tracking
        merges <- rbind(merges, c(wh.clus, pair))
        print(c(wh.clus, pair))
    }else{
        working <- FALSE
    }
}
colnames(merges) <- c('clustering', 'cluster1','cluster2')


save(clusMat.m, merges, file='mergeManyPairwise_out.RData')
#}

#
