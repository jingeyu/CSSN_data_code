###############################################################
########################## Figure 2 ###########################
###############################################################
rm(list = ls())
setwd("/Users/jinge.yu/Desktop/code_and_data/Simulation")
library(ggplot2)
library(pheatmap)

###############################################################
################## (A) Cell's Spatial Pattern #################
###############################################################
load(paste0("RData/", 20201205, "_Sim_network.RData"))
pal <- c(rgb(221, 160, 221, maxColorValue = 255), 
         rgb(0, 206, 209, maxColorValue = 255),
         rgb(112, 173, 71, maxColorValue = 255), 
         rgb(47, 85, 151, maxColorValue = 255),
         rgb(237, 125, 49, maxColorValue = 240))
pal <- setNames(pal, c("1", "2", "3", "4","5"))
cell.type <- as.numeric(cell.info[,1])

gg <- ggplot(cell.info, aes(x = X, y = Y, col = as.factor(cell.type), shape = as.factor(cell.type)))
pl <- gg + geom_point(size = 2.5) +
  scale_color_manual(values = c(pal[1], pal[2], pal[3], pal[4], pal[5])) +
  theme_bw()+
  theme(legend.text=element_text(size=20, face = 'bold'),
        # legend.key.size = unit(1,"inches"),
        # legend.key.width=unit(0.1,'cm'),
        axis.title.x=element_text(size=25,face = "bold"),
        axis.title.y=element_text(size=25,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold"),
        axis.text.y = element_text(size = 20,face = "bold")
  ) + labs(x = "H", y = "L") + 
  guides(color = guide_legend(title = "Cell Type",
                              title.theme = element_text(size = 22,face = 'bold'),
                              override.aes = list(size = 10,face = 'bold')),
         shape = guide_legend(title = "Cell Type",
                              title.theme = element_text(size = 22, face = 'bold'),
                              override.aes = list(size = 5)))
ggsave(paste0("pics/sim_cell_view_manual_", 20201205, ".png"), pl, width = 9, height = 12)

###############################################################
#######################  (B) ROC Curve ########################
###############################################################
load("RData/ROC_curve.RData")

cols <- c(rgb(224, 190, 56, maxColorValue = 230),
          rgb(231, 125, 114, maxColorValue = 255),
          rgb(110, 156, 248, maxColorValue = 255),
          rgb(83, 183, 76, maxColorValue = 255),
          rgb(227, 112, 221, maxColorValue = 255))

png(filename = "pics/ROC.png", width = 1000, height = 1000, res = 200)
plot(total.FPR.sim, total.TPR.sim, type = "l", lwd = 1.8,
     col = cols[1], xlab = "FPR", ylab = "TPR",
     main = "",cex.lab = 1.2, cex.axis = 0.9)
lines(total_FPR_wgcna, total_TPR_wgcna, col = cols[2], lwd = 1.8)
lines(total_FPR_ct, total_TPR_ct, col = cols[3], lwd = 1.8)
lines(total_FPR_csn, total_TPR_csn, col = cols[4], lwd = 1.8)
lines(total_FPR_csnct, total_TPR_csnct, col = cols[5], lwd = 1.8)
legend('bottomright', c("Two-step", "WGCNA", "CTS", "CSN-joint", "CSN-separate"), 
       col = cols, lty = c(1,1,1,1,1), cex = 1, lwd = 1.8
)
dev.off()


###############################################################
#########################  Figure3   ##########################
###############################################################
###############################################################
######################## (D) Heatmaps #########################
###############################################################
load(paste0("RData/", 20201205, "_Sigma.RData"))
load(paste0("RData/", 20201205, "_Sim_network.RData"))
Corr.true[Corr.true != 0] <- 1
load("RData/Result_20201205_0.1_70.RData")
Sparse.Corr[Sparse.Corr != 0] <- 1
load("RData/CTS_20201205.RData")
ct.Corr[ct.Corr != 0] <- 1
load("RData/WGCNA_20201205.RData")
load("RData/CSN_result.RData")
cs.20201205 <- cs.Sgm[[1]]
load("RData/CSN_result_ct_integrate.RData")
csct.20201205 <- cs.Sgm.celltype[[1]]
n <- ncol(X)
indx <- c(545)
print(cell.info[indx,2:3])

colors_func <- colorRampPalette(c('white', "black"))
colors <- colors_func(2)
filenames1 <- paste0("pics/20201205_True_", indx, ".png")
filenames2 <- paste0("pics/20201205_Sim_", indx, ".png")
filenames3 <- paste0("pics/20201205_CTS_", indx, ".png")
filenames4 <- paste0("pics/20201205_CSN_", indx, ".png")
filenames5 <- paste0("pics/20201205_CSNCT_", indx, ".png")
filenames6 <- paste0("pics/20201205_WGCNA_", indx, ".png")

for(i in 1:length(indx)){
  p1 <- pheatmap(Corr.true[,,indx[i]],
                 color = colors,
                 legend_breaks = c(0,1),
                 fontsize = 20,
                 cluster_cols = F, cluster_rows = F,
                 show_rownames = F, show_colnames = F,
                 width = 4, height = 2.8,
                 filename = filenames1[i]
                 
  )
  
  p2 <- pheatmap(Sparse.Corr[,,indx[i]],
                 color = colors,
                 legend_breaks = c(0,1),
                 fontsize = 20,
                 cluster_cols = F, cluster_rows = F,
                 show_rownames = F, show_colnames = F,
                 width = 4, height = 2.8,
                 filename = filenames2[i]

  )

  p3 <- pheatmap(ct.Corr[,,indx[i]],
                 color = colors,
                 legend_breaks = c(0,1),
                 fontsize = 20,
                 cluster_cols = F, cluster_rows = F,
                 show_rownames = F, show_colnames = F,
                 width = 4, height = 2.8,
                 filename = filenames3[i]

  )

  p4 <- pheatmap(cs.20201205[,,indx[i]],
                 color = colors,
                 legend_breaks = c(0,1),
                 fontsize = 20,
                 cluster_cols = F, cluster_rows = F,
                 show_rownames = F, show_colnames = F,
                 width = 4, height = 2.8,
                 filename = filenames4[i]

  )

  p5 <- pheatmap(csct.20201205[,,indx[i]],
                 color = colors,
                 legend_breaks = c(0,1),
                 fontsize = 20,
                 cluster_cols = F, cluster_rows = F,
                 show_rownames = F, show_colnames = F,
                 width = 4, height = 2.8,
                 filename = filenames5[i]

  )

  p6 <- pheatmap(wgcna.Corr[,,indx[i]],
                 color = colors,
                 legend_breaks = c(0,1),
                 fontsize = 20,
                 cluster_cols = F, cluster_rows = F,
                 show_rownames = F, show_colnames = F,
                 width = 4, height = 2.8,
                 filename = filenames6[i]

  )
  
}


###############################################################
########################### Figure 4 ##########################
###############################################################

###############################################################
###################### (A)-(E) Violin plots ###################
###############################################################
G <- nrow(X)
n <- ncol(X)
K <- 5
degree.true <- matrix(NA, G, n)
degree.cts <- matrix(NA, G, n)
degree.wgcna <- matrix(NA, G, n)
degree.sim <- matrix(NA, G, n)
degree.csn <- matrix(NA, G, n)
degree.csnct <- matrix(NA, G, n)

for(i in 1:n){
  degree.true[,i] <- rowSums(Corr.true[,,i]) - 1
  degree.cts[,i] <- rowSums(ct.Corr[,,i]) - 1
  degree.wgcna[,i] <- rowSums(wgcna.Corr[,,i]) - 1
  degree.sim[,i] <- rowSums(Sparse.Corr[,,i]) - 1
  degree.csn[,i] <- rowSums(cs.20201205[,,i]) - 1
  degree.csnct[,i] <- rowSums(csct.20201205[,,i]) - 1
}

cell.type <- as.vector(cell.info[, 1])
ct <- names(table(cell.type))
for(i in 1:K){
  cell.type <- gsub(ct[i], i, cell.type)
}
cell.type <- as.numeric(cell.type)
ind.cell.type <- list()
for(k in 1:K){
  tmp <- which(cell.type == k)
  ind.cell.type[[k]] <- tmp
}

indx1 <- which.max(apply(degree.true, 1, sd))
filename1 <- paste0("pics/degree_ct_", 1:K, ".png")

cols <- c(rgb(85, 188, 194, maxColorValue = 255),
          rgb(224, 190, 56, maxColorValue = 230),
          rgb(231, 125, 114, maxColorValue = 255),
          rgb(110, 156, 248, maxColorValue = 255),
          rgb(83, 183, 76, maxColorValue = 255),
          rgb(227, 112, 221, maxColorValue = 255))

for(k in 1:K){
  nk <- length(ind.cell.type[[k]])
  p_data <- data.frame(Degree = c(degree.true[indx1, ind.cell.type[[k]]], degree.cts[indx1, ind.cell.type[[k]]],
                                  degree.wgcna[indx1, ind.cell.type[[k]]],degree.sim[indx1, ind.cell.type[[k]]],
                                  degree.csn[indx1, ind.cell.type[[k]]], degree.csnct[indx1, ind.cell.type[[k]]]),
                       Method = as.factor(c(rep("Truth", nk), rep("CTS", nk), rep("WGCNA", nk),
                                            rep("Two-step", nk), rep("CSN-joint", nk), 
                                            rep("CSN-separate", nk))))
  p_data$Method <- factor(p_data$Method, levels = c("Truth", "Two-step", "WGCNA","CTS", "CSN-joint", "CSN-separate"))
  p <- ggplot(data = p_data, aes(x=Method, y=Degree)) +
    geom_violin(aes(fill = Method)) + 
    ylab("Degree") + xlab("Method") + 
    scale_fill_manual(values = cols) + 
    theme(plot.title = element_text(hjust = 0.5),
          legend.text=element_text(size = 15,face = 'bold'),
          legend.title = element_text(size = 18, face = 'bold'),
          axis.text.x = element_text(size = 13,face = 'bold'),
          axis.text.y = element_text(size = 18,face = 'bold'),
          axis.title.x=element_text(size=20,face = 'bold'),
          axis.title.y=element_text(size=20,face = 'bold')
          
    )  
  ggsave(filename1[k], p, width = 10, height = 6) 
}

###############################################################
############### (F) Time computation Violin plots #############
###############################################################
load("RData/time_violin.RData")
# Comparasions
cols <- c(rgb(224, 190, 56, maxColorValue = 230),
          rgb(231, 125, 114, maxColorValue = 255),
          rgb(110, 156, 248, maxColorValue = 255),
          rgb(83, 183, 76, maxColorValue = 255),
          rgb(227, 112, 221, maxColorValue = 255))

p_data <- data.frame(Time = c(time.sim, time.wgcna, time.ctn, time.csn, time.csn.ct),
                     Method = as.factor(c(rep("Two-step", 10), rep("WGCNA", 10), rep("CTS", 10),
                                          rep("CSN-joint", 10), rep("CSN-separate", 10))))
p_data$Method <- factor(p_data$Method, levels = c("Two-step", "WGCNA", "CTS", "CSN-joint", "CSN-separate"))
p <- ggplot(data = p_data, aes(x=Method, y=Time)) +
  geom_violin(aes(fill = Method)) + 
  scale_fill_manual(values = cols) + 
  ylab("Time (Seconds)") + xlab("Method") + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size = 15,face = 'bold'),
        legend.title = element_text(size = 18, face = 'bold'),
        axis.text.x = element_text(size = 13,face = 'bold'),
        axis.text.y = element_text(size = 18,face = 'bold'),
        axis.title.x=element_text(size=20,face = 'bold'),
        axis.title.y=element_text(size=20,face = 'bold')
        
  )  

ggsave( "pics/time_violin.png", p, width = 9, height = 6) 


###############################################################
########################## Table 1 ############################
###############################################################

#------Mean error of and standard deviation of 10 replications -------
load("RData/methods_error.RData")
round(mean(sum.diff.sim),2)
round(sd(sum.diff.sim),2)
round(mean(sum.diff.wgcna),2)
round(sd(sum.diff.wgcna),2)
round(mean(sum.diff.ct),2)
round(sd(sum.diff.ct),2)
round(mean(sum.diff.csn),2)
round(sd(sum.diff.csn),2)
round(mean(sum.diff.csnct),2)
round(sd(sum.diff.csnct),2)




