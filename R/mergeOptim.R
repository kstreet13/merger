

clus1 <- clus1.m <- allen
clus2 <- clus2.m <- rsec
#clus1 <- clus1.m <- sample(allen)
#clus2 <- clus2.m <- sample(rsec)

clusters1.orig <- unique(clus1)
clusters2.orig <- unique(clus2)

ARI <- function(mat,matchoose2){
    ei <- sum(choose(rowSums(mat),2))*sum(choose(colSums(mat),2)) / choose(sum(mat),2)
    (sum(matchoose2) - ei) / (.5*(sum(choose(rowSums(mat),2)) + sum(choose(colSums(mat),2))) - ei)
}

n <- length(clus1)
mat <- mat.m <- table(clus1,clus2)
mat2 <- mat2.m <- choose(mat,2)

# check equal length
require(mclust)
base.ARI <- ARI(mat,mat2)
best.ARI <- base.ARI


# MERGE TWICE, SPLIT (if improves ARI)

working <- TRUE
while(working){
    clus1.m.0 <- clus1.m
    clus2.m.0 <- clus2.m
    # merge twice
    for(lcv in 1:2){
        clusters1 <- rownames(mat.m)
        clusters2 <- colnames(mat.m)
        nclus1 <- nrow(mat.m)
        nclus2 <- ncol(mat.m)
        # test all pairwise clusters to merge
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
        mergeResults2 <- sapply(1:nclus2,function(ii){
            sapply(1:nclus2,function(jj){
                if(jj > ii & all(clusters2[c(ii,jj)] != '-1')){
                    newcol <- rowSums(mat.m[,c(ii,jj)])
                    mat.k <- mat.m[,-c(ii,jj)]
                    mat.k <- cbind(mat.k, newcol)
                    mat2.k <- mat2.m[,-c(ii,jj)]
                    mat2.k <- cbind(mat2.k, choose(newcol,2))
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
    }
    
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
















################################################################################
### TWO-STEP VERSION ###
################################################################################




clus1 <- clus1.m <- allen
clus2 <- clus2.m <- rsec
#clus1 <- clus1.m <- sample(allen)
#clus2 <- clus2.m <- sample(rsec)

ARI <- function(mat,matchoose2){
    ei <- sum(choose(rowSums(mat),2))*sum(choose(colSums(mat),2)) / choose(sum(mat),2)
    (sum(matchoose2) - ei) / (.5*(sum(choose(rowSums(mat),2)) + sum(choose(colSums(mat),2))) - ei)
}

n <- length(clus1)
mat <- mat.m <- table(clus1,clus2)
mat2 <- mat2.m <- choose(mat,2)

# check equal length
require(mclust)
base.ARI <- ARI(mat,mat2)
best.ARI <- base.ARI

# THIS IS THE TWO-STEP VERSION
# (it doesn't quite work)
merge <- TRUE
while(merge){
    print(dim(mat.m))
    clusters1 <- rownames(mat.m)
    clusters2 <- colnames(mat.m)
    nclus1 <- nrow(mat.m)
    nclus2 <- ncol(mat.m)
    # test all pairwise clusters to merge
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
    mergeResults1 <- list(max = max(mergeResults1), 
        clusters = clusters1[c(which.max(rowMaxs(mergeResults1)),which.max(colMaxs(mergeResults1)))])
    mergeResults2 <- sapply(1:nclus2,function(ii){
        sapply(1:nclus2,function(jj){
            if(jj > ii & all(clusters2[c(ii,jj)] != '-1')){
                newcol <- rowSums(mat.m[,c(ii,jj)])
                mat.k <- mat.m[,-c(ii,jj)]
                mat.k <- cbind(mat.k, newcol)
                mat2.k <- mat2.m[,-c(ii,jj)]
                mat2.k <- cbind(mat2.k, choose(newcol,2))
                return(ARI(mat.k, mat2.k))
            }else{
                return(-1)
            }
        })
    })
    mergeResults2 <- list(max = max(mergeResults2), 
        clusters = clusters2[c(which.max(rowMaxs(mergeResults2)),which.max(colMaxs(mergeResults2)))])
    mergeResults12 <- list(max = -1, clusters = NA)
    for(i1 in 1:nclus1){
        for(j1 in 1:nclus1){
            if(j1 > i1 & all(clusters1[c(i1,j1)] != '-1')){
                # search over pairs of clusters in clusters2 such that one cluster has non-zero overlap with clusters1[i1] or clusters1[j1]
                clusters2.k <- colnames(mat.m)[which(colSums(mat.m[c(i1,j1),]) > 0)]
                clusters2.k <- clusters2.k[clusters2.k != '-1']
                
                newrow <- colSums(mat.m[c(ii,jj),])
                mat.k <- mat.m[-c(ii,jj),]
                mat.k <- rbind(mat.k, newrow)
                mat2.k <- mat2.m[-c(ii,jj),]
                mat2.k <- rbind(mat2.k, choose(newrow,2))
                
                for(i2 in 1:length(clusters2.k)){
                    for(j2 in 1:length(clusters2.k)){
                        if(j2 > i2 & all(clusters2.k[c(i2,j2)] != '-1')){
                            newcol <- rowSums(mat.k[,c(i2,j2)])
                            mat.kk <- mat.k[,-c(i2,j2)]
                            mat.kk <- cbind(mat.kk, newcol)
                            mat2.kk <- mat2.k[,-c(ii,jj)]
                            mat2.kk <- cbind(mat2.kk, choose(newcol,2))
                            ari <- ARI(mat.kk, mat2.kk)
                            if(ari > mergeResults12$max){
                                mergeResults12 <- list(max = ari, clusters = c(clusters1[c(i1,j1)],clusters2.k[c(i2,j2)]))
                            }
                        }
                    }
                }
                
            }
        }
    }
    # find best pair to merge, only merge if it improves ARI
    bestmerge <- which.max(sapply(list(mergeResults1,mergeResults2,mergeResults12),function(x){x$max}))
    if(max(sapply(list(mergeResults1,mergeResults2,mergeResults12),function(x){x$max})) > best.ARI){
        if(bestmerge == 1){
            best.ARI <- mergeResults1$max
            tomerge <- mergeResults1$clusters
            clus1.m[clus1.m %in% tomerge] <- paste0(tomerge,collapse=',')
            ii <- which(clusters1 == tomerge[1])
            jj <- which(clusters1 == tomerge[2])
            newrow <- colSums(mat.m[c(ii,jj),])
            mat.m <- mat.m[-c(ii,jj),]
            mat.m <- rbind(mat.m, newrow)
            mat2.m <- mat2.m[-c(ii,jj),]
            mat2.m <- rbind(mat2.m, choose(newrow,2))
            rownames(mat.m)[nrow(mat.m)] <- paste0(tomerge,collapse=',')
        }
        if(bestmerge == 2){
            best.ARI <- mergeResults2$max
            tomerge <- mergeResults2$clusters
            clus2.m[clus2.m %in% tomerge] <- paste0(tomerge,collapse=',')
            ii <- which(clusters2 == tomerge[1])
            jj <- which(clusters2 == tomerge[2])
            newcol <- rowSums(mat.m[,c(ii,jj)])
            mat.m <- mat.m[,-c(ii,jj)]
            mat.m <- cbind(mat.m, newcol)
            mat2.m <- mat2.m[,-c(ii,jj)]
            mat2.m <- cbind(mat2.m, choose(newcol,2))
            colnames(mat.m)[ncol(mat.m)] <- paste0(tomerge,collapse=',')
        }
        if(bestmerge == 3){
            best.ARI <- mergeResults2$max
            # merge rows
            tomerge <- mergeResults12$clusters[1:2]
            clus1.m[clus1.m %in% tomerge] <- paste0(tomerge,collapse=',')
            i1 <- which(clusters1 == tomerge[1])
            j1 <- which(clusters1 == tomerge[2])
            newrow <- colSums(mat.m[c(i1,j1),])
            mat.m <- mat.m[-c(i1,j1),]
            mat.m <- rbind(mat.m, newrow)
            mat2.m <- mat2.m[-c(i1,j1),]
            mat2.m <- rbind(mat2.m, choose(newrow,2))
            rownames(mat.m)[nrow(mat.m)] <- paste0(tomerge,collapse=',')
            # merge cols
            tomerge <- mergeResults12$clusters[3:4]
            clus2.m[clus2.m %in% tomerge] <- paste0(tomerge,collapse=',')
            i2 <- which(clusters2 == tomerge[1])
            j2 <- which(clusters2 == tomerge[2])
            newcol <- rowSums(mat.m[,c(i2,j2)])
            mat.m <- mat.m[,-c(i2,j2)]
            mat.m <- cbind(mat.m, newcol)
            mat2.m <- mat2.m[,-c(i2,j2)]
            mat2.m <- cbind(mat2.m, choose(newcol,2))
            colnames(mat.m)[ncol(mat.m)] <- paste0(tomerge,collapse=',')
        }
    }else{
        merge <- FALSE
    }
}










################################################################################
### SPLITTING ###
################################################################################


split <- TRUE
while(split){
    print(1)
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
    mat <- matrix(c(0,1), nrow=2)
    while(ncol(mat) < maxK){
        mat <- rbind(cbind(mat,0), cbind(mat,1))
    }
    # find best split (splitting merged clusters into groups of original clusters)
    splitResults1 <- lapply(map1,function(set){
        k <- length(set)
        mat.k <- mat[1:2^(k-1), 1:k]
        mat.k <- mat.k[-1, , drop=FALSE]
        aris.k <- apply(mat.k,1,function(subset){
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






################################################################################
### ARImp ###
################################################################################
require(clusterExperiment)
clusMat <- readRDS('~/clusMat.rds')
#clusMat <- clusterings
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
sourceCpp('./ARImp.cpp')

arimps <- apply(clusMat,2,ARImp,coClus)

base.ARI <- mean(arimps)
best.ARI <- base.ARI


working <- TRUE
while(working){

    # test all pairwise clusters to merge
    mergeResults <- lapply(1:ncol(clusMat.m),function(wh.clus){
        clus <- clusMat.m[,wh.clus]
        clusters <- unique(clus)
        clusPairs <- combn(clusters[clusters!=-1], 2)[,1:10] #########
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

