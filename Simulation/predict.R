###############################################################
######################### Prediction ##########################
###############################################################
rm(list = ls())
library(CholWishart)
library(MASS)

set.seed(20201231)
# only 20201205
load(paste0("RData/", 20201205, "_Sigma.RData"))
load(paste0("RData/", 20201205, "_Sim_network.RData"))
load("RData/Result_20201205_0.1_70.RData")
Corr.true[abs(Corr.true) != 0] <- 1
Sparse.Corr[Sparse.Corr != 0] <- 1
G <- nrow(X)
n <- ncol(X)

# set number of missing genes
miss.num <- 50
# generate coordinates of missing genes
miss.x <- runif(miss.num, 0, 750)
miss.y <- runif(miss.num, 0, 1000)
# whether there coordinates in cell.info already
sum(miss.x %in% cell.info[,2])
sum(miss.y %in% cell.info[,3])
miss.indx <- cbind(miss.x, miss.y)
# Square neighborhood radius
r <- 80

NeiFind <- function(miss.indx){
  nei.indx <- which(abs(cell.info[,2] - miss.indx[1]) < r & abs(cell.info[,3] - miss.indx[2]) < r)
  cell.type.nei <- cell.type[nei.indx]
  return(cbind(nei.indx, cell.type.nei))
}

cell.type <- cell.info[,1]
ExpSigma <- function(miss.indx, r, nu.i){
  nei.mat <- data.frame(NeiFind(i))
  colnames(nei.mat) <- c("cell.index", "cell.type")
  ni <- nrow(nei.mat)
  if(ni == 0){ 
    Lambda.i <- (nu.i - G - 1) * Sigma.k[,,cell.type[i]]
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

nu <- rep(G + 50, miss.num)
Sigma.miss <- array(NA, dim = c(G, G, miss.num))
X.miss <- matrix(NA, G, miss.num)
Lambda.miss <- array(NA, dim = c(G, G, miss.num))
c.thre <- 0.5
Corr.miss <- array(NA, dim = c(G, G, miss.num))
for(i in 1:miss.num){
  Lambda.miss[,, i] <- ExpSigma(i, r, nu[i])
  Sigma.miss[,, i] <- rInvWishart(1, nu[i], Lambda.miss[,, i])[,,1]
  Sigma.miss[,, i][abs(Sigma.miss[,, i]) < c.thre] <- 0
  # ensure Sigma_i are positive definite
  diag(Sigma.miss[,,i]) <- diag(Sigma.miss[,,i]) + 5
  Corr.miss[,,i] <- diag(diag(Sigma.miss[,,i])^(-0.5)) %*% Sigma.miss[,,i] %*% diag(diag(Sigma.miss[,,i])^(-0.5))
  X.miss[,i] <- mvrnorm(1, mu = rep(0, G), Sigma = Sigma.miss[,,i])
}

Corr.miss[Corr.miss != 0] <- 1

#-------- Predictions of missing cells--------
est.miss <- array(NA, dim = c(G, G, miss.num))
for(i in 1:miss.num){
  miss.nei <- NeiFind(miss.indx[i,])
  tmp <- Sparse.Corr[,, miss.nei[,1]]
  tmp1 <- apply(tmp, 1:2, mean)
  tmp1[tmp1 < 0.5] <- 0
  tmp1[tmp1 >= 0.5] <- 1
  est.miss[,, i] <- tmp1
}

pre.error <- rep(0, miss.num)
for(i in 1:miss.num){
  pre.error[i] <- sum(abs(est.miss[,,i][upper.tri(est.miss[,,i])] - Corr.miss[,,i][upper.tri(Corr.miss[,,i])]))
}
sum(pre.error) / miss.num
# 347.84
