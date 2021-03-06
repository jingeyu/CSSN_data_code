###############################################################
################# Gene Network Function #######################
###############################################################
# packages needed 
# suppressMessages(library(dplyr))

CellGeneCov <- function(X, cell.info, 
                        nu, rho = 0.1,
                        cell.den = 70,
                        is.scale = TRUE){
  
  
  # 
  # Args:
  # X, he gene expression matrix (G, n). Each row of X is gene, and each column is cell.
  # cell.info, a (n, 3) matrix. The first column is cell type indicator, and the last two columns are X and Y centroid coordinate index of each cell respectively.
  # nu, n dimension vector, the degree of freedom of the prior Inverse-Wishart distribution of Sigma_i
  # cell.den, numeric parameter, stands for the average cell numbers in one neighborhood area. The default is 70.
  # rho, numeric thresholding parameter, the range is 0~1. The default is 0.1.
  # is.scale, bool, if TRUE, the variance of gene expression data in each cell type will be scaled to 1. The default is TRUE.
  
  #---- Step0 Preprocess of Data ------
  # Cell number:
  n <- ncol(X)
  # Gene number
  G <- nrow(X)
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
  
  cell.info <- as.data.frame(cell.info)
  colnames(cell.info) <- c("CT", "X", "Y")
  #---- Step2 Estimate Sigma_k ----
  Sgm.hat <- array(NA, dim = c(G,G,K))
  for(k in 1:K){
    nk <- length(ind.cell.type[[k]])
    Sgm.hat[,,k] <- X[, ind.cell.type[[k]]] %*% t(X[, ind.cell.type[[k]]]) / (nk - 1)
  }
  
  
  #---- Step3 Find the Neighbors of Cell i ----
  #--- Set neighborhood radius r ---
  # according to average density of cells in the whole tissue slice ---
  
  # Area of cell tissue slice
  W <- max(cell.info[,2]) - min(cell.info[,2])
  L <- max(cell.info[,3]) - min(cell.info[,3])
  
  # cell.info.upperleft <- filter(cell.info, cell.info[,2] < W / 10 & cell.info[,3] > 9 * L / 10)
  # cell.info.upperright <- filter(cell.info, cell.info[,2] > 9 * W / 10 & cell.info[,3] > 9 * L / 10)
  # cell.info.lowerleft <- filter(cell.info, cell.info[,2] < W / 10 & cell.info[,3] < L / 10)
  # cell.info.lowerright <- filter(cell.info, cell.info[,2] > 9 * W / 10 & cell.info[,3] < L / 10)
  
  # indx <- rbind(c(0,L), c(L, W), c(0,0), c(W,0))
  # # index of four vertices
  # ul <- cell.info.upperleft[which.min(colSums((t(cell.info.upperleft[, 2:3]) - indx[1,])^2)),2:3]
  # ur <- cell.info.upperright[which.min(colSums((t(cell.info.upperright[, 2:3]) - indx[2,])^2)),2:3]
  # ll <- cell.info.lowerleft[which.min(colSums((t(cell.info.lowerleft[, 2:3]) - indx[3,])^2)),2:3]
  # lr <- cell.info.lowerright[which.min(colSums((t(cell.info.lowerright[, 2:3]) - indx[4,])^2)),2:3]
  # indx.quad <- rbind(as.numeric(ul), as.numeric(ur), as.numeric(lr), as.numeric(ll))
  # # diagonal lines of the quadrangle
  # AC <- indx.quad[3,] - indx.quad[1,]
  # BD <- indx.quad[4,] - indx.quad[2,]
  # AreaQuad <- function(x, y){
  #   # x, y are diagonal lines of quadrangle
  #   return(abs(y[1] * x[2] - y[2] * x[1]) / 2)
  # }
  # # area of the quadrangle
  # S <- AreaQuad(AC, BD)
  
  S <- L * W
  # Square neighborhood radius
  r <- sqrt(S * cell.den / n) / 2

  NeiFind <- function(i){
    ind.i <- as.numeric(cell.info[i,2:3])
    nei.inx <- which(abs(cell.info[,2] - ind.i[1]) < r & abs(cell.info[,3] - ind.i[2]) < r)
    # Remove the i-th cell
    cell.index <- nei.inx[nei.inx != i]
    cell.type.nei <- cell.type[cell.index]
    # Return a matrix with index of neighborhood cells and corresponding cell type.
    return(cbind(cell.index, cell.type.nei))
  }
  
  
  #---- Step4 Estimate Posterior Expectation of Sigma_i ----
  SgmEst <- function(nu.i, i, K, Sgm.hat){
    nei.mat <- data.frame(NeiFind(i))
    colnames(nei.mat) <- c("cell.index", "cell.type")
    
    ni <- nrow(nei.mat)
    if(ni == 0){ 
      Lambda.i <- (nu.i - G - 1) * Sgm.hat[,,cell.type[[i]]]
    }else{
      cell.label <- as.integer(names(table(nei.mat$cell.type)))
      # cell.label <- sort(unique(nei.mat$cell.type))
      nei.nk <- as.numeric(table(nei.mat$cell.type))
      weight <- nei.nk / ni
      tmp <- 0
      for(j in 1:length(cell.label)){
        tmp <- tmp + Sgm.hat[,,cell.label[j]] * weight[j]
      }
      Lambda.i <- (nu.i - G - 1) * tmp 
    }
    # mu is known to be 0.
    # posterior mean
    S.i <- X[,i] %*% t(X[,i])
    Sgm.i.hat <- (S.i + Lambda.i) / (nu.i - G)
    
    return(Sgm.i.hat)
  }
  
  Sgm.est <- array(NA, dim = c(G, G, n))
  Corr.est <- array(NA, dim = c(G, G, n))
  for(i in 1:n){ 
    tmp <- SgmEst(nu[i], i, K, Sgm.hat)
    # transform covariance matrix into correlation matrix
    Corr.est[,,i] <- diag(diag(tmp)^(-0.5)) %*% tmp %*% diag(diag(tmp)^(-0.5))
    Corr.est[,,i] <- cov2cor(tmp)
    Sgm.est[,,i] <- tmp
    # write.table(i, "run.txt", append = TRUE, col.names = FALSE)
  }
  
  #---- Sparse estimation of correlation matrix ----
  Corr.est[abs(Corr.est) < rho] <- 0
  
  
  # The function Cell_Gene_Network returns the correlation index list,
  # The first one is posterior mean estimation of Sgm, and the second element is the sparse estimation of the posterior covariance matrix.
  return(list("Posterior Mean" = Sgm.est, 
              "Sparse Correlation Matrix" = Corr.est,
              "Neighborhood radius" = r))
}