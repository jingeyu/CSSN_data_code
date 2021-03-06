rm(list = ls())
setwd("/Users/jinge.yu/Desktop/code_and_data/Real_Application_mouse")
load("ID35_processed.RData")
load("Sparse_est_70_0.1.RData")
load("WGCNA_result.RData")
library(ggplot2)
library(pheatmap)
Sparse.Corr[Sparse.Corr != 0] <- 1
genes <- rownames(X)

##############################################################
########################### Figure 5 #########################
##############################################################

##############################################################
################### (A) Cell's Spatial Pattern ###############
##############################################################
cell.type <- as.vector(cell.info[, 1])
gg <- ggplot(cell.info, aes(x = X, y = Y, col = as.factor(cell.type)))

pl <- gg + geom_point(size = 4, alpha = 0.9) +
  # scale_color_manual(values = c(pal[1], pal[2], pal[3], pal[4], pal[5])) +
  theme_bw()+
  theme(legend.text=element_text(size=30,face = 'bold'),
        axis.title.x=element_text(size=50),
        axis.title.y=element_text(size=50),
        axis.text.x = element_text(size = 30,face = "bold"),
        axis.text.y = element_text(size = 30,face = "bold")
  ) + labs(x = "X", y = "Y") + 
  guides(color = guide_legend(title = "Cell class",
                              title.theme = element_text(size = 50,face = 'bold'),
                              override.aes = list(size = 20)
  ))

ggsave("pics/ID35.png", pl, width = 25, height = 20)


##############################################################
############ (B) Pak3-Grpr spatial pattern (Two-step )########
##############################################################
pal <- c(rgb(255,255,255, maxColorValue = 255),
         rgb(0,0,0, maxColorValue = 255))
filename <- paste0("pics/gene-gene-spatial_", paste0(genes[88],"_",genes[56]), ".png")
gg <- ggplot(cell.info, aes(x = X, y = Y, col = as.factor(Sparse.Corr[88, 56,])))
pl <- gg + geom_point(size = 4, alpha = 0.9) +
  scale_color_manual(values = c(pal[1], pal[2])) +
  theme_bw()+
  theme(legend.text=element_text(size=30, face = 'bold'),
        legend.key = element_rect(fill = "grey"),
        axis.title.x=element_text(size=50),
        axis.title.y=element_text(size=50),
        axis.text.x = element_text(size = 30,face = "bold"),
        axis.text.y = element_text(size = 30,face = "bold")
  ) + labs(x = "X", y = "Y") + 
  guides(color = guide_legend(title = 'Edge',
                              title.theme = element_text(size = 30,face = 'bold'),
                              override.aes = list(size = 20)
  ))
ggsave(filename, pl, width = 25, height = 20)

##############################################################
############ (C) Pak3-Grpr spatial pattern (Two-step )########
##############################################################
filename <- paste0("pics/gene-gene-spatial_WGCNA_0.1_",paste0(genes[88],"_",genes[56]), ".png")
wgcna.Corr[wgcna.Corr != 0] <- 1
gg <- ggplot(cell.info, aes(x = X, y = Y, col = as.factor(wgcna.Corr[88, 56,])))
pl <- gg + geom_point(size = 4, alpha = 0.9) +
  scale_color_manual(values = c(pal[1], pal[2])) +
  theme_bw()+
  theme(legend.text=element_text(size=30, face = 'bold'),
        legend.key = element_rect(fill = "grey"),
        axis.title.x=element_text(size=50),
        axis.title.y=element_text(size=50),
        axis.text.x = element_text(size = 30,face = "bold"),
        axis.text.y = element_text(size = 30,face = "bold")
  ) + labs(x = "X", y = "Y") + 
  guides(color = guide_legend(title = 'Edge',
                              title.theme = element_text(size = 30,face = 'bold'),
                              override.aes = list(size = 20)
  ))
ggsave(filename, pl, width = 25, height = 20)


##############################################################
######################## Figure 6 Heatmaps ###################
##############################################################
set.seed(1996)
cell.type <- cell.info[,1]
c1 <- sample(which(cell.type == "Inhibitory"), 1)
c2 <- sample(which(cell.type == "Excitatory"), 1)
cell <- c(c1, c2)
cell.info[cell,]
colors_func <- colorRampPalette(c('white', "black"))
colors <- colors_func(2)
genes <- rownames(X)

colnames(Sparse.Corr[,,cell[1]]) <- genes
rownames(Sparse.Corr[,,cell[1]]) <- genes
list = pheatmap(Sparse.Corr[,,cell[1]],color = colors,
                cluster_cols = T, cluster_rows = T,
)
newOrder <- Sparse.Corr[,,cell[1]][list$tree_row$order, list$tree_col$order]
colnames(newOrder) <- genes[list$tree_col$order]
rownames(newOrder) <- genes[list$tree_row$order]

png(filename = paste0("pics/Heatmap_clu_",cell[1], ".png"), width = 5400, height = 5000, res = 200)
pheatmap(newOrder, color = colors,
         legend_breaks = c(0,1),
         legend = FALSE,
         fontsize = 15.5,
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = F,
         cutree_rows = 3, cutree_cols = 3,
         labels_row = genes[list$tree_row$order]
)
dev.off()

newOrder2 <- Sparse.Corr[,,cell[2]][list$tree_row$order, list$tree_col$order]
rownames(newOrder2) <- genes[list$tree_row$order]
png(filename = paste0("pics/Heatmap_clu_v",cell[2], ".png"), width = 5400, height = 5000, res = 200)
pheatmap(newOrder2, color = colors,
         legend_breaks = c(0,1),
         legend = FALSE,
         fontsize = 15.5,
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = F,
         cutree_rows = 3, cutree_cols = 3,
         labels_row = genes[list$tree_row$order]
)
dev.off()


##############################################################
################### Figure 7 Violin plots ###################
##############################################################
n <- ncol(X)
G <- nrow(X)
degree.sim <- matrix(NA, G, n)
for(i in 1:n){
  degree.sim[,i] <- rowSums(Sparse.Corr[,,i]) - 1
}

cell.type <- as.vector(cell.info[, 1])
ct <- names(table(cell.type))
K <- length(table(cell.type))
for(i in 1:K){
  cell.type <- gsub(ct[i], i, cell.type)
}
cell.type <- as.numeric(cell.type)
ind.cell.type <- list()
nk <- rep(0, K)
for(k in 1:K){
  tmp <- which(cell.type == k)
  ind.cell.type[[k]] <- tmp
  nk[k] <- length(ind.cell.type[[k]])
}

indx1 <- order(apply(degree.sim, 1, sd), decreasing = TRUE)[1:2]
filename <- paste0("pics/degree_", genes[indx1], ".png")
Cell.type <- NULL
Degree.ct <- list()
for(i in 1:length(indx1)){
  Degree.ct[[i]] <- 0
  for(k in 1:K){
    Degree.ct[[i]] <- c(Degree.ct[[i]], degree.sim[indx1[i], ind.cell.type[[k]]])
  }
  Degree.ct[[i]] <- Degree.ct[[i]][-1]
}

for(k in 1:K){
  Cell.type <- c(Cell.type, rep(ct[k], nk[k]))
}

cols <- c(rgb(234,139, 129, maxColorValue = 255),
          rgb(217, 153, 66, maxColorValue = 255),
          rgb(191, 166, 66, maxColorValue = 255),
          rgb(157, 178, 67, maxColorValue = 255),
          rgb(104, 186, 66, maxColorValue = 255),
          rgb(99, 193, 134, maxColorValue = 255),
          rgb(100, 196, 180, maxColorValue = 255),
          rgb(97, 191, 218, maxColorValue = 255),
          rgb(91, 179, 246, maxColorValue = 255),
          rgb(152, 160, 248, maxColorValue = 255),
          rgb(205, 137, 247, maxColorValue = 255),
          rgb(233, 124, 219, maxColorValue = 255),
          rgb(235, 216, 179, maxColorValue = 255)
          )

for(i in 1:length(indx1)){
  p_data <- data.frame(Degree = Degree.ct[[i]],
                       Celltype = as.factor(Cell.type))
  p <- ggplot(data = p_data, aes(x=Celltype, y=Degree)) +
    geom_violin(aes(fill = Celltype)) + 
    ylab("Degree") + xlab("Cell Class") + 
    scale_color_manual(values = cols) + 
    scale_fill_discrete("")+
    theme(plot.title = element_text(hjust = 0.5),
          legend.position="none",
          axis.text.x = element_text(size = 13,face = 'bold'),
          axis.text.y = element_text(size = 18,face = 'bold'),
          axis.title.x=element_text(size=20,face = 'bold'),
          axis.title.y=element_text(size=20,face = 'bold')
          
    ) 
  ggsave(filename[i], p, width = 22, height = 6) 
}


#----- Gene list -------
# The most variable 15 and stable 15 genes (degree) in Inhibitory and Excitatory neurons
cell.type <- as.vector(cell.info[, 1])
Exict.cell <- which(cell.type == "Excitatory")
Inhi.cell <- which(cell.type == "Inhibitory")
inhi.indx <- order(apply(degree.sim[, Inhi.cell], 1, sd), decreasing = TRUE)
exict.indx <- order(apply(degree.sim[, Exict.cell], 1, sd), decreasing = TRUE)
tab.inhi <- cbind(genes[head(inhi.indx, 15)], genes[rev(tail(inhi.indx, 15))])
tab.exict <- cbind(genes[head(exict.indx,15)], genes[rev(tail(exict.indx, 15))])
colnames(tab.inhi) <- c("Most variable genes", "Most stable genes")
colnames(tab.exict) <- c("Most variable genes", "Most stable genes")
write.csv(tab.inhi, file = "Inhibitory degree genes.csv")
write.csv(tab.exict, file = "Excitatory degree genes.csv")






