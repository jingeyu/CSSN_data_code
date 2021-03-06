##############################################################
####################### Preprocessing ########################
##############################################################
rm(list = ls())
setwd("/Users/jinge.yu/Desktop/code_and_data/Real_Application_human")
library(Seurat)

# total RNA counts
dat1 <- read.csv("pnas.1912459116.sd12.csv", header = TRUE)
rownames(dat1) <- dat1[,1]
dat1 <- dat1[,-1]
dim(dat1)

# Batch 1 coordinates
coor1 <- read.csv("Coor_B1.csv", header = TRUE)
rownames(coor1) <- coor1[,1]
coor1 <- coor1[,-1]

# Delete Blank genes
dat1 <- dat1[-c(953:3805),]
dat1_1 <- dat1[, 1:nrow(coor1)]

# ---- Choose Batch 1-----
#Normalization
tmp1 <- 1e6 / colSums(dat1_1)
RNA_norm <- sweep(dat1_1, 2, tmp1, '*')
RNA_sd <- apply(RNA_norm, 1, sd)
RNA_500 = RNA_norm[order(RNA_sd, decreasing = T),] %>% head(500)
RNA_c <- log2(1+RNA_500)



#############################
#Conduct Seurat analysis now
#############################
library(Seurat)
batch1 <- CreateSeuratObject(counts = dat1_1, 
                             min.cells = 3, min.features = 200,
                             project = "cell")
dim(batch1)

#normalization
batch1 <-NormalizeData(batch1, normalization.method = "LogNormalize", scale.factor = 10000)

#identify highly variable features
#the selected features are used in downstream PCA analysis
batch1 <- FindVariableFeatures(batch1, selection.method = "vst", nfeatures = 2000)

#scaling the data, each gene has mean zero and variance one
all.genes <- rownames(dat1_1)
batch1 <- ScaleData(batch1, features = all.genes)
#perform PCA (linear dimensional reduction)
batch1 <- RunPCA(batch1, features = VariableFeatures(object = batch1), verbose = F)
nPC <- 15
batch1 <- FindNeighbors(batch1, dims = 1:nPC)
batch1 <- FindClusters(batch1, resolution = 0.8)
#cluster id

table(Idents(batch1))

cell.type <- as.numeric(Idents(batch1))
cell.info <- as.data.frame(cbind(cell.type, coor1))

X <- as.matrix(RNA_c)

save(list = c("X", "cell.info"), file = "data.RData")

