#!/bin/bash
#SBATCH --job-name=mergeMany
#SBATCH --cpus-per-task 4
#SBATCH -o mergeMany.out
#SBATCH -e mergeMany.err

R-3.5.0 CMD BATCH --vanilla R/mergeManyPairwise.R mergeMany.Rout
