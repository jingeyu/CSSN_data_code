###############################################################
######################### Comparasion #########################
###############################################################
rm(list = ls())
# setwd("/Users/jinge.yu/Desktop/code_and_data/Simulation")

sch <- seq.Date(from = as.Date("20201205",format = "%Y%m%d"), by = "day", length.out = 10)
sch <- as.numeric(gsub("-", "", sch))
c.thre <- 0.1
G <- 100

#---Initialize ----
sum.diff.ct <- rep(0, 10)
sum.diff.sim <- rep(0, 10)
sum.diff.csn <- rep(0, 10)
sum.diff.csnct <- rep(0, 10)

time.ctn <- rep(0, 10)
time.sim <- rep(0, 10)
time.csn <- rep(0, 10)
time.csn.ct <- rep(0, 10)

#---- Cell-type specific results ------
for(s in 1:10){
  seed <- sch[s]
  # load data
  load(paste0("RData/", seed, "_Sigma.RData"))
  load(paste0("RData/", seed, "_Sim_network.RData"))
  Corr.true[Corr.true != 0] <- 1
  n <- dim(Corr.true)[3]
  
  t1 <- proc.time()
  # Cell number:
  n <- ncol(X)
  # Gene number
  G <- nrow(X)
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
  }
  
  ct.Corr[abs(ct.Corr) < c.thre] <- 0
  tmp <- ct.Corr
  tmp[tmp != 0] <- 1
  time.ctn[s] <- (proc.time() - t1)[3]
  print(time.ctn[s])
  save(ct.Corr, file = paste0("RData/CTS_", seed, ".RData"))
  #------Step 2-ERRORS----
  diff_ct <- NULL
  for(i in 1:n){
    diff_ct <- c(diff_ct, sum(abs(Corr.true[,,i][upper.tri(Corr.true[,,i])] - tmp[,,i][upper.tri(tmp[,,i])])))
  }
  sum.diff.ct[s] <- sum(diff_ct) / n
}


#-------Two-step algorithm Simulation result---------
cell.den <- 70
source("cs_network.R")
# sum.diff.my <- NULL

for(s in 1:10){
  seed <- sch[s]
  load(paste0("RData/", seed, "_Sim_network.RData"))
  load(paste0("RData/", seed, "_Sigma.RData"))
  
  Corr.true[Corr.true != 0] <- 1
  
  G <- nrow(X)
  n <- ncol(X)
  nu <- rep(G+G, n)
  t1 <- proc.time()
  Result <- CellGeneCov(X, cell.info, cell.den = cell.den,
                        nu = nu, is.scale = TRUE, rho = c.thre)
  time.sim[s] <- (proc.time() - t1)[3]
  print(time.sim[s])
  
  Sparse.Corr <- Result$`Sparse Correlation Matrix`
  tmp <- Sparse.Corr
  tmp[tmp != 0] <- 1
  save(Sparse.Corr, file = paste0("RData/Result_", seed, "_", c.thre, "_", cell.den,".RData"))
  
  #---- Estimation errors ----
  diff_my <- NULL
  for(i in 1:n){
    diff_my <- c(diff_my, sum(abs(Corr.true[,,i][upper.tri(Corr.true[,,i])] - tmp[,,i][upper.tri(tmp[,,i])])))
  }
  sum.diff.sim[s] <- sum(diff_my) / n
  
}


#---------Cell-specific gene network ------------
Corr.true.all <- list()
for(s in 1:10){
  load(paste0("RData/", sch[s], "_Sigma.RData"))
  Corr.true[Corr.true != 0] <- 1
  Corr.true.all[[s]] <- Corr.true
}
G <- ncol(Corr.true.all[[1]])
K <- 5
num.cell <- c()
for(i in 1:10){
  num.cell[i] <- dim(Corr.true.all[[i]])[3]
}

#---- load csn results ----
library(R.matlab)
csn <- readMat("RData/csn.mat")
#---- transform csn mat into ordinary matrix ----
cs.Sgm <- list()
for(num in 1:10){
  cs.Sgm[[num]] <- array(NA, dim = c(G, G, num.cell[num]))
  for(i in 1:num.cell[num]){
    cs.Sgm[[num]][,,i] <- as.numeric(as.matrix(csn[[1]][[50+num]][[1]][[i]][[1]]))
  }
}
save(cs.Sgm, file = "RData/CSN_result.RData")

cs.Sgm.ct <- list()
for(num in 1:10){
  cs.Sgm.ct[[num]] <- list()
  for(j in 1:K){
    n.num <- length(csn[[1]][[(num -1) * K + j]][[1]])
    cs.Sgm.ct[[num]][[j]] <- array(NA, dim = c(G, G, n.num))
    for(i in 1:n.num){
      cs.Sgm.ct[[num]][[j]][,,i] <- as.numeric(as.matrix(csn[[1]][[(num -1) *K + j]][[1]][[i]][[1]]))
    }
  }
}

#---- CSN-joint ----
for(num in 1:10){
  diff.csn1 <- NULL
  for(i in 1:num.cell[num]){
    diff.csn1 <- c(diff.csn1, sum(abs(Corr.true.all[[num]][,,i][upper.tri(Corr.true.all[[num]][,,i])] - cs.Sgm[[num]][,,i][upper.tri(cs.Sgm[[num]][,,i])])))
  }
  sum.diff.csn[num] <- sum(diff.csn1) / num.cell[num]
}


#------CSN-separate-------
cs.Sgm.celltype <- list()
for(num in 1:10){
  load(paste0("RData/", sch[num], "_Sim_network.RData"))
  cell.type.rep <- cell.info[,1]
  cs.Sgm.celltype[[num]] <- array(NA, dim = c(G, G, length(cell.type.rep)))
  for(k in 1:K){
    cs.Sgm.celltype[[num]][,,which(cell.type.rep == k)] <- cs.Sgm.ct[[num]][[k]]
  }
}

save(cs.Sgm.celltype, file = "RData/CSN_result_ct_integrate.RData")

for(num in 1:10){
  diff.csn2 <- NULL
  for(i in 1:dim(cs.Sgm[[num]])[3]){
    diff.csn2 <- c(diff.csn2, sum(abs(Corr.true.all[[num]][,,i][upper.tri(Corr.true.all[[num]][,,i])] - cs.Sgm.celltype[[num]][,,i][upper.tri(cs.Sgm.celltype[[num]][,,i])])))
  }
  sum.diff.csnct[num] <- sum(diff.csn2) / num.cell[num]
}

save(list = c('sum.diff.csnct',"sum.diff.csn", "sum.diff.sim", "sum.diff.ct"),
     file =  "RData/methods_error.RData")

