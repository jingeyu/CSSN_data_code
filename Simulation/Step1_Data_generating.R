###############################################################
###################### Data Generation ########################
###############################################################
rm(list = ls())
setwd("/Users/jinge.yu/Desktop/code_and_data/Simulation")
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
library(MASS)
library(Matrix)
library(CholWishart)

colors_func <- colorRampPalette(c('white', "black"))
colors <- colors_func(2)

#set seed to reproduce results
sch <- seq.Date(from = as.Date("20201205",format = "%Y%m%d"), by = "day", length.out = 10)
sch <- as.numeric(gsub("-", "", sch))

# gene number
G <- 100

# cell-type number
K <- 5

# the area of whole tissue cryosection
L <- 1000 # length, micrometre
H <- 750 # width, micrometre
area  <- L * H

for(seed in sch){
  set.seed(seed)
  
  #---- Step1 generate spatial coordinates of cells ----
  # Poisson process
  lambda <- c(2,3,4,5,6) / 100
  # 1. divide the area R into n.region(10000) sub-regions
  sub.num <- 100
  n.region <- sub.num * sub.num
  N <- matrix(NA, K, n.region)
  for(i in 1:n.region){
    N[, i] <- rpois(K, lambda)
  }
  
  # total cell numbers in one tissue cryosection 
  n <- sum(N)
  print(n)
  
  # Random and independently place N[,i] points in n.region sub-regions.
  N.i <- colSums(N)
  xlab <- seq(0, H, length.out = sub.num + 1)
  ylab <- seq(0, L, length.out =sub.num + 1)
  centroid.x <- NULL
  centroid.y <- NULL
  for(i in 1:sub.num){
    for(j in 1:sub.num){
      centroid.x <- c(centroid.x, runif(N.i[(i-1) * sub.num + j], xlab[i], xlab[i+1]))
      centroid.y <- c(centroid.y, runif(N.i[(i-1) * sub.num + j], ylab[j], ylab[j+1]))
    }
  }
  
  #---- set spatial pattern manually----
  cell.type <- NULL
  ct <- 1:K
  for(i in 1:n.region){
    cell.type <- c(cell.type, rep(ct, N[,i]))
  }
  
  cell.type[1:500] <- sample(1:K, 500, prob = c(0.2, 0.01, 0.05, 0.2, 0.64), replace = TRUE)
  cell.type[501:1000] <- sample(1:K, 500, prob = c(0.1, 0.05, 0.35, 0.3, 0.2), replace = TRUE)
  cell.type[801:1000] <- sample(4:5, 200, prob = c(0.5, 0.5), replace = TRUE)
  cell.type[1001:1500] <- sample(1:K, 500, prob = c(0.05, 0.65, 0.1, 0.1, 0.1), replace = TRUE)
  cell.type[1501:n] <- sample(1:K, n-1501+1, prob = c(0.3, 0.3, 0.2, 0.15, 0.05), replace = TRUE)
  cell.type[1601:1700] <- sample(c(2,3), prob = c(0.7,0.3), replace = TRUE)
  cell.type[1800:n] <- sample(1:3, n-1800+1, prob = c(0.8, 0.1, 0.2), replace = TRUE)
  length(cell.type)
  table(cell.type)
  
  t1 <- which(centroid.x < 710 & centroid.x > 650 & centroid.y < 850 & centroid.y > 750)
  t2 <- which(centroid.x < 500 & centroid.x > 350 & centroid.y < 550 & centroid.y > 230)
  t3 <- which(centroid.x < 600 & centroid.x > 400 & centroid.y < 820 & centroid.y > 600)
  t4 <- which(centroid.x < 200 & centroid.x > 0 & centroid.y < 250 & centroid.y > 0)
  t5 <- which(centroid.x < 300 & centroid.x > 0 & centroid.y < 750 & centroid.y > 250)
  t6 <- which(centroid.x < 710 & centroid.x > 650 & centroid.y < 850 & centroid.y > 750)
  t7 <- which(centroid.x < 710 & centroid.x > 650 & centroid.y < 850 & centroid.y > 750)
  
  cell.type[t1] <- sample(c(1,4, 5), length(t1), prob = c(0.8, 0.15, 0.05), replace = TRUE)
  cell.type[t2] <- sample(c(3,4,5), length(t2), prob = c(0.8, 0.15, 0.05), replace = TRUE)
  cell.type[t3] <- sample(c(2,4,5), length(t3), prob = c(0.7, 0.2, 0.1), replace = TRUE)
  cell.type[t4] <- sample(c(1,5,3,2), length(t4), prob = c(0.7, 0.1,0.1,0.1), replace = TRUE)
  cell.type[t5] <- sample(c(5,1,3,4,2), length(t5), prob = c(0.7, 0.15, 0.07, 0.05, 0.03), replace = TRUE)
  cell.type[t6] <- sample(c(3,4), length(t6), prob = c(0.8, 0.2), replace = TRUE)
  cell.type[t7] <- sample(c(3,4), length(t7), prob = c(0.8, 0.2), replace = TRUE)
  
  
  cell.info <- data.frame(cell.type, centroid.x, centroid.y)
  colnames(cell.info) <- c("CT", "X", "Y")
  
  #---- Step 2 Generate true Sigma.k (Sparse + Positive define) ----
  Indicator <- function(i, j, number){
    if(abs(i - j) == number){
      return(1)
    }else{
      return(0)
    }
  }
  
  isPosDef <- function(M) { 
    if (isSymmetric(M)) {  
      # first test symmetricity
      if (all(eigen(M)$values > 0)) {
        return(TRUE)
      }else {
        return(FALSE)
      } 
    }else{
      return(FALSE)
    }  # not symmetric
  }
  
  # K blocks
  Sigma.k <- array(0, dim = c(G, G, K))
  sub.piece <- G / K
  subSgm.1 <- matrix(0, sub.piece, sub.piece)
  subSgm.2 <- matrix(0, sub.piece, sub.piece)
  subSgm.3 <- matrix(0, sub.piece, sub.piece)
  subSgm.4 <- matrix(0, sub.piece, sub.piece)
  
  rho <- 0.7
  for(i in 1:sub.piece){
    for(j in 1:sub.piece){
      subSgm.1[i, j] <- rho^(abs(i - j))
      subSgm.2[i, j] <- 1 - abs(i - j) / 10
      subSgm.3[i, j] <- -0.3 * Indicator(i, j, 1) + 1.3 * Indicator(i, j, 0)
      subSgm.4[i, j] <- 1 - 2 * abs(i - j) / G
    }
  }
  
  diag(subSgm.1) <- diag(subSgm.1) + 0.5
  subSgm.2[subSgm.2 < 0] <- 0
  subSgm.4[subSgm.4 < 0] <- 0
  subSgm.5 <- matrix(NA, sub.piece, sub.piece)
  b <- matrix(NA, sub.piece, sub.piece)
  b[upper.tri(b, diag = TRUE)] <- runif(sub.piece * (sub.piece + 1) / 2, -0.2, 0.8) * rbinom(sub.piece * (sub.piece + 1) / 2, 1, 0.2)
  B <- t(b)
  B[upper.tri(t(B))] <- b[upper.tri(b)]
  epsilon <- max(c(0, -min(eigen(B)$values))) + 0.01
  subSgm.5 <- B + epsilon * diag(1, sub.piece)
  
  Sigma.k[,,1] <- as.matrix(bdiag(list(subSgm.1, subSgm.2, subSgm.3, subSgm.4,subSgm.5)))
  Sigma.k[,,2] <- as.matrix(bdiag(list(subSgm.1, subSgm.3, subSgm.2, subSgm.4,subSgm.5)))
  Sigma.k[,,3] <- as.matrix(bdiag(list(subSgm.1, subSgm.3, subSgm.2, subSgm.5,subSgm.4)))
  Sigma.k[,,4] <- as.matrix(bdiag(list(subSgm.3, subSgm.1, subSgm.2, subSgm.5,subSgm.4)))
  Sigma.k[,,5] <- as.matrix(bdiag(list(subSgm.3, subSgm.2, subSgm.5, subSgm.1,subSgm.4)))
  
  #---- Step 3 Generate true Sigma.ki (Sparse + Positive definite) ----
  # given the freedom of Inver-Wishart distribution:
  nu <- rep(G + 50, n)
  r = 80
  NeiFind <- function(i, r){
    ind.i <- as.numeric(cell.info[i,2:3])
    nei.inx <- which(abs(cell.info[,2] - ind.i[1]) < r & abs(cell.info[,3] - ind.i[2]) < r)
    # Remove the i-th cell
    cell.index <- nei.inx[nei.inx != i]
    cell.type <- cell.info[cell.index, 1]
    # Return a matrix with index of neighborhood cells and corresponding cell type.
    return(cbind(cell.index, cell.type))
  }
  
  #----- Expectation of Sigma.i(Inverse-Wishart) ----
  ExpSigma <- function(i, r, nu.i){
    nei.mat <- data.frame(NeiFind(i, r))
    colnames(nei.mat) <- c("cell.index", "cell.type")
    ni <- nrow(nei.mat)
    if(ni == 0){ 
      Lambda.i <- (nu.i - G - 1) * Sigma.k[,,cell.type[[i]]]
    }else{
      cell.label <- as.integer(names(table(nei.mat$cell.type)))
      nei.nk <- as.numeric(table(nei.mat$cell.type))
      weight <- nei.nk / ni
      tmp <- 0
      for(j in 1:length(cell.label)){
        tmp <- tmp + Sigma.k[,,cell.label[j]] * weight[j]
      }
      Lambda.i <- (nu.i - G - 1) * tmp 
    }
  }
  
  Sigma.i <- array(NA, dim = c(G, G, n))
  X <- matrix(NA, G, n)
  Lambda <- array(NA, dim = c(G, G, n))
  c.thre <- 0.5
  Corr.true <- array(NA, dim = c(G, G, n))
  for(i in 1:n){
    Lambda[,,i] <- ExpSigma(i, r, nu[i])
    Sigma.i[,,i] <- rInvWishart(1, nu[i], Lambda[,,i])[,,1]
    Sigma.i[,,i][abs(Sigma.i[,,i]) < c.thre] <- 0
    # ensure Sigma_i are positive definite
    diag(Sigma.i[,,i]) <- diag(Sigma.i[,,i]) + 5
    Corr.true[,,i] <- diag(diag(Sigma.i[,,i])^(-0.5)) %*% Sigma.i[,,i] %*% diag(diag(Sigma.i[,,i])^(-0.5))
    X[,i] <- mvrnorm(1, mu = rep(0, G), Sigma = Sigma.i[,,i])
  }
  
  
  #-----Save data--------
  save(list = c("X", "cell.info"), file = paste0("RData/", seed, "_Sim_network.RData"))
  save(list = c("Sigma.k", "Corr.true"), file = paste0("RData/", seed, "_Sigma.RData"))
  
  
}