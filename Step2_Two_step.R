##############################################################
##################### Two-step algorithm #####################
##############################################################
rm(list = ls())
# setwd("/Users/jinge.yu/Desktop/code_and_data/Real_Application")
load("ID35_processed.RData")
source("cs_network.R")

G <- nrow(X)
nu <- rep(2*G, ncol(X))
c.thre <- 0.1
cell.den <- 70

Result <- CellGeneCov(X, cell.info, cell.den = cell.den,
                      nu = nu, is.scale = TRUE, rho = c.thre)
Sparse.Corr <- Result$`Sparse Correlation Matrix`
save(Sparse.Corr, file = "Sparse_est_70_0.1.RData")

