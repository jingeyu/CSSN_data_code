##############################################################
##################### Two-step algorithm #####################
##############################################################
rm(list = ls())
setwd("/Users/jinge.yu/Desktop/code_and_data/Real_Application")
load("ID35_processed.RData")
source("cs_network.R")
library(WGCNA)

G <- nrow(X)
n <- ncol(X)
nu <- rep(2*G, ncol(X))
c.thre <- 0.1
cell.den <- 70

#------Two step algorithm-------
Result <- CellGeneCov(X, cell.info, cell.den = cell.den,
                      nu = nu, is.scale = TRUE, rho = c.thre)
Sparse.Corr <- Result$`Sparse Correlation Matrix`
save(Sparse.Corr, file = "Sparse_est_70_0.1.RData")


#------WGCNA-------


# Notice that cell types in cell.info as factor
cell.type <- as.vector(cell.info[, 1])
# Cell Type number
K <- length(unique(cell.type))
# Transform factor cell types to numeric ones.
ct <- names(table(cell.type))
for(i in 1:length(ct)){
  cell.type <- gsub(ct[i], i, cell.type)
}
cell.type <- as.numeric(cell.type)
ind.cell.type <- list()

for(k in 1:K){
  tmp <- which(cell.type == k)
  ind.cell.type[[k]] <- tmp
}

softPower <- c(4,3,3,3,2,2,5,4,6,2,2,2,2)
TOM <- array(NA, dim = c(G,G,K))
for(k in 1:K){
  datExpr <- t(X[,ind.cell.type[[k]]])
  adjacency <- adjacency(datExpr, power = softPower[k])
  TOM[,,k] <- TOMsimilarity(adjacency)
}

wgcna.Corr <- array(NA, dim = c(G, G, n))
for(i in 1:n){
  wgcna.Corr[,,i] <- TOM[,,cell.type[i]]
}

wgcna.Corr[wgcna.Corr < 0.1] <- 0
save(wgcna.Corr, file = "WGCNA_result.RData")


