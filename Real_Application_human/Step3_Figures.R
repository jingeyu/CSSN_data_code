##############################################################
######################### Figure 8 ###########################
##############################################################
rm(list = ls())
setwd("/Users/jinge.yu/Desktop/code_and_data/Real_Application_human")
library(dplyr)
library(ggplot2)
library(pheatmap)
load("Result_70_0.3.RData")

##############################################################
################### (A) Cells' spatial pattern ###############
##############################################################
cell.type <- cell.info[,1]
colnames(cell.info) <- c("Celltype", "X", "Y")
gg <- ggplot(cell.info, aes(x = X, y = Y, 
                            col = as.factor(cell.type), 
                            shape = as.factor(cell.type)))
pl <- gg + geom_point(size = 4, alpha = 0.9) +
  theme_bw()+
  theme(legend.text=element_text(size=30,face = 'bold'),
        axis.title.x=element_text(size=50),
        axis.title.y=element_text(size=50),
        axis.text.x = element_text(size = 30,face = "bold"),
        axis.text.y = element_text(size = 30,face = "bold")
  ) + labs(x = "X", y = "Y") + 
  guides(color = guide_legend(title = "Cell Class",
                              title.theme = element_text(size = 22,face = 'bold'),
                              override.aes = list(size = 10,face = 'bold')
  ),
  shape = guide_legend(title = "Cell Class",
                       title.theme = element_text(size = 22, face = 'bold'),
                       override.aes = list(size = 10,face = 'bold')))

ggsave("pics/spatial.png", pl, width = 15, height = 10)


##############################################################
####################### (B)-(F) Heatmaps #####################
##############################################################
set.seed(1996)
cell <- NULL
for(k in 1:K){
  cell <- c(cell, sample(which(cell.type == k), 1))
}
cell.info[cell,]
genes <- rownames(X)
colors_func <- colorRampPalette(c('white', "black"))
colors <- colors_func(2)

colnames(Sparse.Corr[,,cell[1]]) <- genes
rownames(Sparse.Corr[,,cell[1]]) <- genes
list = pheatmap(Sparse.Corr[,,cell[1]],color = colors,
                cluster_cols = T, cluster_rows = T,
)
newOrder <- Sparse.Corr[,,cell[1]][list$tree_row$order, list$tree_col$order]
colnames(newOrder) <- genes[list$tree_col$order]
rownames(newOrder) <- genes[list$tree_row$order]

png(filename = paste0("pics/Heatmap_clu_v_",cell[1], ".png"), width = 2300, height = 2000, res = 200)
pheatmap(newOrder, color = colors,
         # legend_breaks = c(0,1),
         legend = FALSE,
         fontsize = 15.5,
         cluster_cols = F, cluster_rows = F,
         show_rownames = F, show_colnames = F,
         cutree_rows = 3, cutree_cols = 3,
         labels_row = genes[list$tree_row$order]
)
dev.off()

for(i in 2:length(cell)){
  filenames <- paste0("pics/Heatmap_clu_v_",cell[i], ".png")
  newOrder2 <- Sparse.Corr[,,cell[i]][list$tree_row$order, list$tree_col$order]
  rownames(newOrder2) <- genes[list$tree_row$order]
  png(filename = filenames, width = 2300, height = 2000, res = 200)
  pheatmap(newOrder2, color = colors,
           legend = FALSE,
           fontsize = 15.5,
           cluster_cols = F, cluster_rows = F,
           show_rownames = F, show_colnames = F,
           cutree_rows = 3, cutree_cols = 3,
           labels_row = genes[list$tree_row$order]
  )
  dev.off()
}


##############################################################
#################### Figure 9 Violin plots ###################
##############################################################
#----degree------
cell.type <- as.vector(cell.info[,1])
genes <- rownames(X)
degree.sim <- matrix(NA, G, n)
for(i in 1:n){
  degree.sim[, i] <- rowSums(Sparse.Corr[,,i]) - 1
}
nk <- c()
ind.cell.type <- list()
for(k in 1:K){
  ind.cell.type[[k]] <- which(cell.type == k)
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
  Cell.type <- c(Cell.type, rep(k, nk[k]))
}

cols <- c(rgb(234,139, 129, maxColorValue = 255),
          rgb(173, 173, 66, maxColorValue = 255),
          rgb(99, 194, 143, maxColorValue = 255),
          rgb(168, 212, 246, maxColorValue = 255),
          rgb(239, 203, 247, maxColorValue = 255)
          
)

for(i in 1:length(indx1)){
  p_data <- data.frame(Degree = Degree.ct[[i]],
                       Celltype = as.factor(Cell.type))
  p <- ggplot(data = p_data, aes(x=Celltype, y=Degree)) +
    geom_violin(aes(fill = Celltype)) + 
    labs(x = "Cell Class", y = "Degree") + 
    scale_color_manual(values = cols) + 
    theme(plot.title = element_text(hjust = 0.5),
          # legend.text=element_text(size = 15,face = 'bold'),
          legend.position="none",
          # legend.title = element_text(size = 18, face = 'bold'),
          axis.text.x = element_text(size = 20,face = 'bold'),
          axis.text.y = element_text(size = 20,face = 'bold'),
          axis.title.x=element_text(size=25,face = 'bold'),
          axis.title.y=element_text(size=25,face = 'bold')
          
    ) 
  ggsave(filename[i], p, width = 16, height = 6) 
}







