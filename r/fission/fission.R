rm(list=ls())

library(mixMVPLN)
library(fission)

data("fission")

sel <- cbind(fission@colData$strain, fission@colData$minute, fission@colData$replicate)
sel_wt_r1 <- fission@colData$strain == "wt" & fission@colData$replicate == "r1"
sel_wt_r2 <- fission@colData$strain == "wt" & fission@colData$replicate == "r2"
sel_wt_r3 <- fission@colData$strain == "wt" & fission@colData$replicate == "r3"
sel_mut_r1 <- fission@colData$strain == "mut" & fission@colData$replicate == "r1"
sel_mut_r2 <- fission@colData$strain == "mut" & fission@colData$replicate == "r2"
sel_mut_r3 <- fission@colData$strain == "mut" & fission@colData$replicate == "r3"

counts <- fission@assays$data$counts

counts_wt_r1 <- t(matrix(as.numeric(counts[,sel_wt_r1]), ncol=7039, nrow=6)) 
counts_wt_r2 <- t(matrix(as.numeric(counts[,sel_wt_r2]), ncol=7039, nrow=6))
counts_wt_r3 <- t(matrix(as.numeric(counts[,sel_wt_r3]), ncol=7039, nrow=6)) 
counts_mut_r1 <- t(matrix(as.numeric(counts[,sel_mut_r1]), ncol=7039, nrow=6)) 
counts_mut_r2 <- t(matrix(as.numeric(counts[,sel_mut_r2]), ncol=7039, nrow=6))
counts_mut_r3 <- t(matrix(as.numeric(counts[,sel_mut_r3]), ncol=7039, nrow=6))

lista <- list(counts_wt_r1,
              counts_wt_r2,
              counts_wt_r3,
              counts_mut_r1,
              counts_mut_r2,
              counts_mut_r3)

mat_wt <- matrix(NA, ncol=6, nrow=7039)
for (i in 1:7039){
    for (j in 1:6){
        mat_wt[i,j] <- median(c(lista[[1]][i,j], lista[[2]][i,j], lista[[3]][i,j]))
    }
}
mat_mut <- matrix(NA, ncol=6, nrow=7039)
for (i in 1:7039){
    for (j in 1:6){
        mat_mut[i,j] <- median(c(lista[[4]][i,j], lista[[5]][i,j], lista[[6]][i,j]))
    }
}
lista <- list()
for (i in 1:7039){
    lista[[i]] <- rbind(mat_wt[i,], mat_mut[i,]) |> t()
}

lista

mvplnVGAclus(lista,
             normalize = "Yes",
             nInitIterations = 10,
             gmin = 1,
             gmax = 2)


source("customvga.R")
customvga(lista,
          normalize = "Yes",
          nInitIterations = 10,
          gmin = 1,
          gmax = 2)
