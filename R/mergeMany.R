################################################################################
### ARImp ###
################################################################################
require(clusterExperiment)
clusMat <- readRDS('data/clusMat.rds')
clusMat <- apply(clusMat,2,function(x){
    x[x != '-1'] <- as.numeric(factor(x[x != '-1']))
    x[x == '-1'] <- -1
    x <- as.integer(x)
}) # convert to integer
coClus <- makeConsensus(clusMat, proportion = 0.8)$percentageShared
denom <- (clusMat!=-1) %*% t(clusMat!=-1)
num <- coClus * denom
clusters <- apply(clusMat,2,unique)

clusMat.m <- clusMat
coClus.m <- coClus

# n x n similarity matrix, M / partition Q
require(Rcpp)
sourceCpp('./src/ARImp.cpp')

arimps <- apply(clusMat,2,ARImp,coClus)

base.ARI <- mean(arimps)
best.ARI <- base.ARI


working <- TRUE
while(working){
    
    # test all pairwise clusters to merge
    mergeResults <- lapply(1:ncol(clusMat.m),function(wh.clus){
        clus <- clusMat.m[,wh.clus]
        clusters <- unique(clus)
        clusPairs <- combn(clusters[clusters!=-1], 2) #########
        
        # test_pairs(clus, clusPairs, num, denom)
        
        mcmapply(function(pp){
            ind.ii <- which(clus == clusPairs[1,pp])
            ind.jj <- which(clus == clusPairs[2,pp])
            # clus and coClus.m will reset with each iteration
            clus[c(ind.ii,ind.jj)] <- max(clusters)+1
            #coClus.m[ind.ii, ind.jj] <- (num[ind.ii,ind.jj] + 1) / denom[ind.ii,ind.jj]
            #coClus.m[ind.jj, ind.ii] <- (num[ind.jj,ind.ii] + 1) / denom[ind.jj,ind.ii]
            #ARImp(clus, coClus.m)
            update_arimp(Q=clus, indii=ind.ii, indjj=ind.jj, M=coClus.m,
                update=(num[ind.ii,ind.jj] + 1) / denom[ind.ii,ind.jj]) - arimps[wh.clus]
        }, 1:ncol(clusPairs), mc.cores = 2) #1:ncol(clusPairs)
    })
    
    mergeResults1 <- sapply(1:nclus1,function(ii){
        sapply(1:nclus1,function(jj){
            if(jj > ii & all(clusters1[c(ii,jj)] != '-1')){
                newrow <- colSums(mat.m[c(ii,jj),])
                mat.k <- mat.m[-c(ii,jj),]
                mat.k <- rbind(mat.k, newrow)
                mat2.k <- mat2.m[-c(ii,jj),]
                mat2.k <- rbind(mat2.k, choose(newrow,2))
                return(ARI(mat.k, mat2.k))
            }else{
                return(-1)
            }
        })
    })
    # find best pair to merge, only merge if it improves ARI
    merged.aris <- unlist(c(mergeResults1,mergeResults2))
    best.ARI <- max(c(mergeResults1,mergeResults2))
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
    print(tomerge)
    
    # split as needed
    split <- TRUE
    while(split){
        clusters1 <- unique(clus1.m)
        clusters2 <- unique(clus2.m)
        nclus1 <- length(unique(clus1.m))
        nclus2 <- length(unique(clus2.m))
        
        map1 <- sapply(clusters1,function(c1){
            unique(clus1[clus1.m==c1])
        })
        map1 <- map1[sapply(map1,length) > 1]
        map2 <- sapply(clusters2,function(c2){
            unique(clus2[clus2.m==c2])
        })
        map2 <- map2[sapply(map2,length) > 1]
        
        # make maximal table of possible subsets (2^max(K)-2 by K)
        maxK <- max(unlist(c(sapply(c(map1,map2),length),1)))
        submat <- matrix(c(0,1), nrow=2)
        while(ncol(submat) < maxK){
            submat <- rbind(cbind(submat,0), cbind(submat,1))
        }
        # find best split (splitting merged clusters into groups of original clusters)
        splitResults1 <- lapply(map1,function(set){
            k <- length(set)
            submat.k <- submat[1:2^(k-1), 1:k]
            submat.k <- submat.k[-1, , drop=FALSE]
            aris.k <- apply(submat.k,1,function(subset){
                clus.split1 <- set[as.logical(subset)]
                clus.split2 <- set[!as.logical(subset)]
                
                foo <- clus1.m
                foo[clus1 %in% clus.split1] <- 'foo.clus.split1'
                foo[clus1 %in% clus.split2] <- 'foo.clus.split2'
                
                adjustedRandIndex(foo, clus2.m)
            })
            return(aris.k)
        })
        splitResults2 <- lapply(map2,function(set){
            k <- length(set)
            mat.k <- mat[1:2^(k-1), 1:k]
            mat.k <- mat.k[-1, , drop=FALSE]
            aris.k <- apply(mat.k,1,function(subset){
                clus.split1 <- set[as.logical(subset)]
                clus.split2 <- set[!as.logical(subset)]
                
                foo <- clus2.m
                foo[clus2 %in% clus.split1] <- 'foo.clus.split1'
                foo[clus2 %in% clus.split2] <- 'foo.clus.split2'
                
                adjustedRandIndex(clus1.m, foo)
            })
            return(aris.k)
        })
        # only split if it improves ARI
        split.aris <- unlist(c(-Inf,splitResults1,splitResults2))
        if(max(split.aris) > best.ARI){
            print('split')
            best.ARI <- max(split.aris)
            if(max(c(-Inf,unlist(splitResults1))) > max(c(-Inf,unlist(splitResults2)))){
                tosplit <- names(which.max(sapply(splitResults1,max)))
                set <- map1[[tosplit]]
                k <- length(set)
                mat.k <- mat[1:2^(k-1), 1:k]
                mat.k <- mat.k[-1, , drop=FALSE]
                split.idx <- as.logical(mat.k[which.max(splitResults1[[tosplit]]),])
                clus1.m[clus1 %in% set[split.idx]] <- paste(set[split.idx],collapse=',')
                clus1.m[clus1 %in% set[!split.idx]] <- paste(set[!split.idx],collapse=',')          
            }else{
                tosplit <- names(which.max(sapply(splitResults2,max)))
                set <- map2[[tosplit]]
                k <- length(set)
                mat.k <- mat[1:2^(k-1), 1:k]
                mat.k <- mat.k[-1, , drop=FALSE]
                split.idx <- as.logical(mat.k[which.max(splitResults2[[tosplit]]),])
                clus2.m[clus2 %in% set[split.idx]] <- paste(set[split.idx],collapse=',')
                clus2.m[clus2 %in% set[!split.idx]] <- paste(set[!split.idx],collapse=',')
            }
        }else{
            split <- FALSE
        }
    }
    
    if(all(clus1.m == clus1.m.0) & all(clus2.m == clus2.m.0)){
        working <- FALSE
    }
}
#
