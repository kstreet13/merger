################################################################################
### pairwise ARIs ###
################################################################################
library(clusterExperiment)
library(parallel)
library(matrixStats)
library(here)
library(mclust)
source(here("R", "helper.R"))

clusMat <- readRDS(here("data", "clusMat.rds"))
InitialARI <- apply(clusMat, 2, function(x) {
  apply(clusMat, 2, function(y) {
    mclust::adjustedRandIndex(x, y)
  })
})
plotARIs(InitialARI)

mergers <- mergeManyPairwise(clusteringMatrix = clusMat, nCores = 3)

FinalARI <- apply(mergers$currentMat, 2, function(x) {
  apply(mergers$currentMat, 2, function(y) {
    mclust::adjustedRandIndex(x, y)
  })
})

plotARIs(FinalARI)