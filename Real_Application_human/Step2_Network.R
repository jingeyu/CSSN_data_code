##############################################################
#################### Two-step algorithm ######################
##############################################################
rm(list = ls())
setwd("/Users/jinge.yu/Desktop/code_and_data/Real_Application_human")

load("data.RData")
G <- nrow(X)
n <- ncol(X)
K <- 5
source("cs_network.R")

cell.den <- 70
c.thre <- 0.3
nu <- rep(G+G, n)

t1 <- proc.time()
Result <- CellGeneCov(X, cell.info, cell.den = cell.den,
                      nu = nu, is.scale = TRUE, rho = c.thre)
print(proc.time() - t1)

Sparse.Corr <- Result$`Sparse Correlation Matrix`
Sparse.Corr[Sparse.Corr != 0] <- 1
save(Sparse.Corr, file = "Result_70_0.3.RData")




