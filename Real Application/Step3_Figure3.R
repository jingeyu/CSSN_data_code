##############################################################
########################### Figure 3 #########################
##############################################################
rm(list = ls())
# setwd("/Users/jinge.yu/Desktop/code_and_data/Real_Application")
load("ID35_processed.RData")
load("Sparse_est_70_0.1.RData")
library(ggplot2)
library(pheatmap)

##############################################################
################### (a) Cell's Spatial Pattern ###############
##############################################################
cell.type <- as.vector(cell.info[, 1])
gg <- ggplot(cell.info, aes(x = X, y = Y, col = as.factor(cell.type)))

pl <- gg + geom_point(size = 4, alpha = 0.9) +
  # scale_color_manual(values = c(pal[1], pal[2], pal[3], pal[4], pal[5])) +
  theme_bw()+
  theme(legend.text=element_text(size=30,face = 'bold'),
        # legend.key.size = unit(1,"inches"),
        # legend.key.width=unit(0.1,'cm'),
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
################ (b) Pak3-Grpr spatial pattern ###############
##############################################################
genes <- rownames(X)
Sparse.Corr[Sparse.Corr != 0] <- 1
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
######################## (c) Heatmaps ########################
##############################################################
set.seed(3)
c1 <- sample(which(cell.type == "Inhibitory"), 1)
c2 <- sample(which(cell.type == "Excitatory"), 1)
cell <- c(c1, c2)
cell.info[cell,]
filename <- paste0("pics/Heatmap_", cell, ".png")
colors_func <- colorRampPalette(c('white', "black"))
colors <- colors_func(2)
for(i in 1:length(cell)){
  colnames(Sparse.Corr[,,cell[i]]) <- genes
  rownames(Sparse.Corr[,,cell[i]]) <- genes
  list = pheatmap(Sparse.Corr[,,cell[i]],color = colors,
                  cluster_cols = T, cluster_rows = T,
  )
  newOrder=Sparse.Corr[,,cell[i]][list$tree_row$order, list$tree_col$order]
  colnames(newOrder) <- genes[list$tree_col$order]
  rownames(newOrder) <- genes[list$tree_row$order]
  
  png(filename = filename[i], width = 5400, height = 5000, res = 200)
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
}



