################################################################################
### ARImp ###
################################################################################
require(clusterExperiment)
NCORES = 2
clusMat <- readRDS('data/clusMat.rds')
clusMat <- apply(clusMat,2,function(x){
    x[x != '-1'] <- as.numeric(factor(x[x != '-1']))
    x[x == '-1'] <- -1
    x <- as.integer(x)
}) # convert to integer
coClus <- makeConsensus(clusMat, proportion = 0.8)$percentageShared
denom <- (clusMat!=-1) %*% t(clusMat!=-1)
num <- coClus * denom
rm(coClus)
clusters <- apply(clusMat,2,unique)

clusMat.m <- clusMat

# n x n similarity matrix, M / partition Q
require(Rcpp)
sourceCpp('./src/ARImp.cpp')

arimps <- apply(clusMat,2,ARImp,num/denom)

base.ARI <- mean(arimps)
best.ARI <- base.ARI

require(parallel)
working <- TRUE
merges <- NULL
while(working){
    # test all pairwise clusters to merge
    mergeResults <- mclapply(1:ncol(clusMat.m),function(wh.clus){
        clus <- clusMat.m[,wh.clus]
        clusternames <- clusters[[wh.clus]]
        clusPairs <- combn(clusternames[clusternames!=-1], 2)[,1:20] #########
        test_pairs(clus, clusPairs, num, denom) - arimps[wh.clus]
    }, mc.cores = NCORES)
    
    # find best pair to merge, only merge if it improves ARI
    maxes <- sapply(mergeResults,max)
    if(max(maxes) > 0){
        # update arimps
        wh.clus <- which.max(maxes)
        arimps[wh.clus] <- arimps[wh.clus] + max(maxes)
        # update clusters
        clusternames <- clusters[[wh.clus]]
        clusPairs <- combn(clusternames[clusternames!=-1], 2)
        pair <- clusPairs[,which.max(mergeResults[[wh.clus]])]
        ind.ii <- which(clusMat.m[,wh.clus] == pair[1])
        ind.jj <- which(clusMat.m[,wh.clus] == pair[2])
        clusMat.m[c(ind.ii,ind.jj), wh.clus] <- min(pair)
        clusters <- apply(clusMat.m,2,unique)
        # update num
        num[ind.ii,ind.jj] <- num[ind.ii,ind.jj] + 1
        num[ind.jj,ind.ii] <- num[ind.jj,ind.ii] + 1
        
        merges <- rbind(merges, c(wh.clus, pair))
        print(pair)
    }else{
        working <- FALSE
    }
}
colnames(merges) <- c('clustering', 'cluster1','cluster2')
#
