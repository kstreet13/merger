
.ARI <- function(mat,matchoose2){
    ei <- sum(choose(rowSums(mat),2))*sum(choose(colSums(mat),2)) / choose(sum(mat),2)
    (sum(matchoose2) - ei) / (.5*(sum(choose(rowSums(mat),2)) + sum(choose(colSums(mat),2))) - ei)
}


mergeTwo <- function(clus1, clus2){
    require(matrixStats)
    n <- length(clus1)
    mat <- mat.m <- table(clus1,clus2)
    mat2 <- mat2.m <- choose(mat,2)
    
    clus1.m <- clus1
    clus2.m <- clus2
    
    if(n != length(clus2)){
        stop('Input clusterings have different lengths, but should represent clusterings of the same data.')
    }
    base.ARI <- .ARI(mat,mat2)
    best.ARI <- base.ARI
    
    working <- TRUE
    while(working){
        clus1.m.0 <- clus1.m
        clus2.m.0 <- clus2.m
        
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
                    return(.ARI(mat.k, mat2.k))
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
                    return(.ARI(mat.k, mat2.k))
                }else{
                    return(-1)
                }
            })
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


