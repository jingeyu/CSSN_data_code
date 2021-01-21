##############################################################
####################### Pre-processing #######################
##############################################################
rm(list = ls())
library(dplyr)

dat <- read.table("/Users/jinge.yu/Desktop/Gene_network/Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv", 
                  header = T, row.names = 1, sep = ",")

#-------data split--------
dat.i <- dat[dat$Animal_ID == 35, ]
dat.i <- dat.i[dat.i$Bregma == 0.26, ]
# delete cells belongs to 'Ambiguous'(Cell_class)
dat.i <- dat.i[-which(dat.i$Cell_class == "Ambiguous"), ]
# delete barcoded genes named as "Blank_"
Blank.ind <- paste0("Blank_",1:5)
dat.i <- dat.i[, -which(colnames(dat.i) %in% Blank.ind)]
cell.info.i <- data.frame(dat.i$Cell_class, 
                          as.numeric(dat.i$Centroid_X),
                          as.numeric(dat.i$Centroid_Y))
colnames(cell.info.i) <- c("CT", "X", "Y")

X <- t(dat.i[, -c(1:8)])
dim(X)
# sum(is.na(X))

#------ Filtering -----
# delete cell-types containing less than 10 cells!
ct.del <- names(which(table(cell.info.i$CT) < 10))
del.cell <- NULL
for(del in 1:length(ct.del)){
  del.cell <- c(del.cell, which(cell.info.i == ct.del[del])) 
}
cell.info <- cell.info.i[-del.cell, ]
rownames(cell.info) <- 1:nrow(cell.info)
X <- X[, -del.cell]

# Cell number:
n <- ncol(X)
# Gene number
G <- nrow(X)
# Notice that cell types in cell.info are factor
cell.type <- cell.info[,1]
K <- length(unique(cell.type))
K
ct <- names(table(cell.type))
ind.cell.type <- list()
X_tmp <- matrix(NA, G, K)
for(k in 1:K){
  tmp <- which(cell.type == ct[k])
  ind.cell.type[[k]] <- tmp
  X_tmp[,k] <- rowSums(X[, ind.cell.type[[k]]])
}


tmp <- rowSums(X_tmp == 0) == 0
X <- X[tmp, ]
dim(X)


# save gene-expression matrix X and cell information matrix cell.info
save(list = c("X", "cell.info"), file = "ID35_processed.RData")
