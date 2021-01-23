###############################################################
########################## ROC Data ###########################
###############################################################
rm(list = ls())
setwd("/Users/jinge.yu/Desktop/code_and_data/Simulation")

c.thre <- seq(0,0.3,0.01)
load(paste0("RData/", 20201205, "_Sigma.RData"))
load(paste0("RData/", 20201205, "_Sim_network.RData"))
Corr.true[Corr.true != 0] <- 1
n <- ncol(X)
# Gene number
G <- nrow(X)

true.edge <- list()
true.zero <- list()
edge.true <- rep(0,n)

# true situation
for(i in 1:n){
  edge.true[i] <- sum(Corr.true[,,i][upper.tri(Corr.true[,,i])])
  true.edge[[i]] <- which(Corr.true[,,i][upper.tri(Corr.true[,,i])] == 1)
  true.zero[[i]] <- which(Corr.true[,,i][upper.tri(Corr.true[,,i])] == 0)
}
non.edge <- G * (G - 1) / 2 - edge.true

#-----CTN---------
# Notice that cell types in cell.info are factor
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

#---- Step1 Centralization and Scaling ----
is.scale <- TRUE
for(k in 1:K){
  tmp <- which(cell.type == k)
  X[, tmp] <- X[, tmp] - rowMeans(X[,tmp])
  if(is.scale == TRUE){
    X[,tmp] <- X[, tmp] / apply(X[, tmp], 1, sd)
  }
  ind.cell.type[[k]] <- tmp
}

NA.num <- sum(is.na(X))
if(NA.num > 0){
  print("Waring! remove genes having zero values across one cell-type!")
  # delete the genes have all 0 expression across one cell-type for the sake of normalization
  ind.na <- which(is.na(X))
  coor.na <- NULL
  for(i in 1:length(ind.na)){
    if(ind.na[i] %% nrow(X) == 0){
      indrow <- nrow(X)
      indcol <- ind.na[i] / nrow(X)
    }else{
      indrow <- ind.na[i] %% nrow(X)
      indcol <- ceiling(ind.na[i] / nrow(X))
    }
    tmp <- c(indrow, indcol)
    coor.na <- rbind(coor.na, tmp)
  }
  del.gen <- unique(coor.na[,1])
  del.gen.num <- length(del.gen)
  print(paste("Delete", del.gen.num, "genes in total."))
  X <- X[-del.gen, ]
  G <- nrow(X)
  # stop("Please remove the genes having zero expression across spots or cells!")
}

Sgm.hat <- array(NA, dim = c(G,G,K))
for(k in 1:K){
  nk <- length(ind.cell.type[[k]])
  Sgm.hat[,,k] <- X[, ind.cell.type[[k]]] %*% t(X[, ind.cell.type[[k]]]) / nk
}

ct.Sigma <- array(NA, dim = c(G, G, n))
ct.Corr <- array(NA, dim = c(G, G, n))
for(i in 1:n){
  ct.Sigma[,,i] <- Sgm.hat[,,cell.type[i]]
  ct.Corr[,,i] <- cov2cor(ct.Sigma[,,i])
  # ct.Corr[,,i] <- diag(diag(ct.Sigma[,,i])^(-0.5)) %*% ct.Sigma[,,i] %*% diag(diag(ct.Sigma[,,i])^(-0.5))
}

#-------ERRORS----
total_FPR_ct <- c()
total_TPR_ct <- c()
FPR.ct <- list()
TPR.ct <- list()
est.edge.ct <- list()
for(s in 1:length(c.thre)){
  tmp <- ct.Corr
  tmp[abs(tmp) < c.thre[s]] <- 0
  tmp[tmp != 0] <- 1
  FPR.ct[[s]] <- rep(0, n)
  TPR.ct[[s]] <- rep(0, n)
  FP <- rep(0,n)
  TP <- rep(0,n)
  est.edge.ct[[s]] <- list()
  diff_ct <- NULL
  for(i in 1:n){
    est.edge.ct[[s]][[i]] <- which(tmp[,,i][upper.tri(tmp[,,i])] == 1)
    FP[i] <- length(intersect(est.edge.ct[[s]][[i]], true.zero[[i]]))
    FPR.ct[[s]][i] <- FP[i] / non.edge[i]
    TP[i] <- length(intersect(est.edge.ct[[s]][[i]], true.edge[[i]]))
    TPR.ct[[s]][i] <- TP[i] / edge.true[i]
  }
  total_FPR_ct[s] <- sum(FP) / sum(non.edge)
  total_TPR_ct[s] <- sum(TP) / sum(edge.true)
}

#--------Simualation----------
cell.den <- 70
source("cs_network.R")
nu <- rep(G+G, n)

FPR.sim <- list()
TPR.sim <- list()
est.edge.sim <- list()
total.FPR.sim <- c()
total.TPR.sim <- c()

for(s in 1:length(c.thre)){
  Result <- CellGeneCov(X, cell.info, cell.den = cell.den,
                        nu = nu, is.scale = TRUE, rho = c.thre[s])
  Sparse.Corr <- Result$`Sparse Correlation Matrix`
  tmp <- Sparse.Corr
  tmp[tmp != 0] <- 1
  
  #---- Estimation errors ----
  FPR.sim[[s]] <- rep(0, n)
  TPR.sim[[s]] <- rep(0, n)
  FP <- rep(0,n)
  TP <- rep(0,n)
  est.edge.sim[[s]] <- list()
  
  for(i in 1:n){
    est.edge.sim[[s]][[i]] <- which(tmp[,,i][upper.tri(tmp[,,i])] == 1)
    FP[i] <- length(intersect(est.edge.sim[[s]][[i]], true.zero[[i]]))
    FPR.sim[[s]][i] <- FP[i] / non.edge[i]
    TP[i] <- length(intersect(est.edge.sim[[s]][[i]], true.edge[[i]]))
    TPR.sim[[s]][i] <- TP[i] / edge.true[i]
  }
  
  total.FPR.sim[s] <- sum(FP) / sum(non.edge)
  total.TPR.sim[s] <- sum(TP) / sum(edge.true)

}


#---------CSN----------
load("RData/CSN_result_alpha1e-30--1e-4.RData")
load("RData/csct_result_alpha1e-30--1e-4.RData")
alpha1 <- c(1e-20, 1e-19, 1e-18, 1e-17, 1e-16, 1e-15, 
            1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 
            1e-8, 1e-7, 1e-6, 1e-5, 1e-4)
alpha2 <- c(1e-30, 1e-29, 1e-28, 1e-27, 1e-26, 1e-25, 1e-24, 1e-23, 1e-22,1e-21)
num.alpha1 <- length(alpha1)
num.alpha2 <- length(alpha2)
alpha <- c(alpha1, alpha2)
num.alpha <- num.alpha1 + num.alpha2

FPR.cs <- list()
TPR.cs <- list()
FPR.csct <- list()
TPR.csct <- list()
est.edge.cs <- list()
est.edge.csct <- list()
total_FPR_csn <- rep(0, num.alpha)
total_TPR_csn <- rep(0, num.alpha)
total_FPR_csnct <- rep(0, num.alpha)
total_TPR_csnct <- rep(0, num.alpha)

for(s in 1:num.alpha){
  FPR.cs[[s]] <- rep(0, n)
  TPR.cs[[s]] <- rep(0, n)
  FP <- rep(0,n)
  TP <- rep(0,n)
  est.edge.cs[[s]] <- list()
  for(i in 1:n){
    est.edge.cs[[s]][[i]] <- which(cs.Sgm[[s]][,,i][upper.tri(cs.Sgm[[s]][,,i])] == 1)
    FP[i] <- length(intersect(est.edge.cs[[s]][[i]], true.zero[[i]]))
    FPR.cs[[s]][i] <- FP[i] / non.edge[i]
    TP[i] <- length(intersect(est.edge.cs[[s]][[i]], true.edge[[i]]))
    TPR.cs[[s]][i] <- TP[i] / edge.true[i]
  }
  total_FPR_csn[s]<- sum(FP) / sum(non.edge)
  total_TPR_csn[s] <- sum(TP) / sum(edge.true)
}

total_FPR_csn <- c(0, total_FPR_csn, 1)
total_TPR_csn <- c(0, total_TPR_csn, 1)

#--------CSN-CT-------
for(s in 1:num.alpha){
  FPR.csct[[s]] <- rep(0, n)
  TPR.csct[[s]] <- rep(0, n)
  FP <- rep(0,n)
  TP <- rep(0,n)
  est.edge.csct[[s]] <- list()
  for(i in 1:n){
    est.edge.csct[[s]][[i]] <- which(cs.Sgm.celltype[[s]][,,i][upper.tri(cs.Sgm.celltype[[s]][,,i])] == 1)
    FP[i] <- length(intersect(est.edge.csct[[s]][[i]], true.zero[[i]]))
    FPR.csct[[s]][i] <- FP[i] / non.edge[i]
    TP[i] <- length(intersect(est.edge.csct[[s]][[i]], true.edge[[i]]))
    TPR.csct[[s]][i] <- TP[i] / edge.true[i]
  }
  total_FPR_csnct[s]<- sum(FP) / sum(non.edge)
  total_TPR_csnct[s] <- sum(TP) / sum(edge.true)
}

total_FPR_csnct <- c(0, total_FPR_csnct,1)
total_TPR_csnct <- c(0, total_TPR_csnct, 1)

#----save------

save(list = c("total.FPR.sim", "total.TPR.sim", "total_FPR_ct", "total_TPR_ct",
              "total_FPR_csn", "total_TPR_csn", "total_FPR_csnct", "total_TPR_csnct"),
     file = "RData/ROC_curve.RData")



