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

clusMat.m <- clusMat

#base.ARI <- mean(arimps)
#best.ARI <- base.ARI

require(parallel)
working <- TRUE
merges <- NULL


.ARI <- function(mat,matchoose2=NULL){
    if(class(mat)=='list'){
        matchoose2 <- mat[[2]]
        mat <- mat[[1]]
    }
    ei <- sum(choose(rowSums(mat),2))*sum(choose(colSums(mat),2)) / choose(sum(mat),2)
    (sum(matchoose2) - ei) / (.5*(sum(choose(rowSums(mat),2)) + sum(choose(colSums(mat),2))) - ei)
}


mergeManyPairwise <- function(clusMat){
    require(matrixStats)
    n <- nrow(clusMat)
    
    confusionMats <- list()
    for(ci1 in 1:(ncol(clusMat)-1)){
        for(ci2 in (ci1+1):ncol(clusMat)){
            mat.ij <- table(clusMat[,ci1],clusMat[,ci2])
            confusionMats[[paste0(c(ci1,ci2),collapse = '')]] <- list(mat = mat.ij, mat2 = choose(mat.ij,2))
        }
    }
    
    clusMat.m <- clusMat

    base.ARI <- sapply(1:ncol(clusMat),function(ci1){
        sapply(1:ncol(clusMat),function(ci2){
            if(ci1==ci2){ return(1) }
            return(.ARI(confusionMats[[paste(c(min(c(ci1,ci2)),max(c(ci1,ci2))),collapse='')]]))
        })
    })
    best.ARI <- base.ARI
    
    working <- TRUE
    while(working){
        
        # test all pairs of clusters to merge
        mergeResults <- lapply(1:ncol(clusMat.m),function(wh.clus){
            aris <- best.ARI[,wh.clus][-wh.clus]
            clus <- clusMat.m[,wh.clus]
            clusternames <- clusters[[wh.clus]]
            clusPairs <- combn(clusternames[clusternames!=-1], 2)
            
            apply(clusPairs,2,function(clp){
                aris.m <- sapply((1:ncol(clusMat.m))[-wh.clus],function(ci2){
                    idx <- paste0(c(wh.clus,ci2),collapse = '')
                    matlist <- confusionMats[[idx]]
                    matlist$mat
                })
                
                # if(jj > ii & all(clusters1[c(ii,jj)] != '-1')){
                #     newrow <- colSums(mat.m[c(ii,jj),])
                #     mat.k <- mat.m[-c(ii,jj),]
                #     mat.k <- rbind(mat.k, newrow)
                #     mat2.k <- mat2.m[-c(ii,jj),]
                #     mat2.k <- rbind(mat2.k, choose(newrow,2))
                #     return(.ARI(mat.k, mat2.k))
                # }else{
                #     return(-1)
                # }
            })
            
            test_pairs(clus, clusPairs, num, denom) - arimps[wh.clus]
        })
            
        # find best pair to merge, only merge if it improves ARI
        merged.aris <- unlist(c(mergeResults1,mergeResults2))
        if(max(merged.aris) > best.ARI){
            best.ARI <- max(merged.aris)
            if(max(mergeResults1) > max(mergeResults2)){
                ii <- which.max(rowMaxs(mergeResults1))
                jj <- which.max(colMaxs(mergeResults1))
                tomerge <- c(clusters1[ii], clusters1[jj])
                clus1.m[clus1.m %in% tomerge] <- paste0(tomerge,collapse=',')
                newrow <- colSums(mat.m[c(ii,jj),])
                mat.m <- mat.m[-c(ii,jj),]
                mat.m <- rbind(mat.m, newrow)
                mat2.m <- mat2.m[-c(ii,jj),]
                mat2.m <- rbind(mat2.m, choose(newrow,2))
                rownames(mat.m)[nrow(mat.m)] <- paste0(tomerge,collapse=',')
            }else{
                ii <- which.max(rowMaxs(mergeResults2))
                jj <- which.max(colMaxs(mergeResults2))
                
                tomerge <- c(clusters2[ii], clusters2[jj])
                clus2.m[clus2.m %in% tomerge] <- paste0(tomerge,collapse=',')
                newcol <- rowSums(mat.m[,c(ii,jj)])
                mat.m <- mat.m[,-c(ii,jj)]
                mat.m <- cbind(mat.m, newcol)
                mat2.m <- mat2.m[,-c(ii,jj)]
                mat2.m <- cbind(mat2.m, choose(newcol,2))
                colnames(mat.m)[ncol(mat.m)] <- paste0(tomerge,collapse=',')
            }
            #print(tomerge)
        }else{
            working <- FALSE
        }
    }
    return(cbind(clus1.m,clus2.m))
}


colnames(merges) <- c('clustering', 'cluster1','cluster2')
#
