## slingshot trajectory
library(slingshot)
library(Seurat)
sce<-as.SingleCellExperiment(neural)

sce <- slingshot(sce, reducedDim = 'PCA')  # no clusters

# Plot PC1 vs PC2 colored by Slingshot pseudotime.
colors <- rainbow(50, alpha = 1)
plot(reducedDims(sce)$PCA, col = colors[cut(sce$slingPseudotime_1,breaks=50)], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2)

library(ggplot2)
library(ggthemes)
library(ggbeeswarm)
library(dplyr)
library(cowplot)
library(RColorBrewer)
library(knitr)

library(ggbeeswarm)
library(ggthemes)
# Plot Slingshot pseudotime vs Phase. 
ggplot(as.data.frame(colData(sce)), aes(x = sce$slingPseudotime_1, y = Phase, 
                              colour = Phase)) +
    geom_quasirandom(groupOnX = FALSE) +
    scale_color_tableau() + theme_classic() +
    xlab("Slingshot pseudotime") + ylab("Timepoint") +
    ggtitle("Cells ordered by Slingshot pseudotime")

ggsave(paste0("pseudotime_slingshot.svg"))

library(gam)

# Only look at the 1,000 most variable genes when identifying temporally expressesd genes.
# Identify the variable genes by ranking all genes by their variance.
Y <- log2(counts(sce) + 1)
var1K <- names(sort(apply(Y, 1, var),decreasing = TRUE))[1:1000]
Y <- Y[var1K, ]  # only counts for variable genes

# Fit GAM for each gene using pseudotime as independent variable.
t <- sce$slingPseudotime_1
gam.pval <- apply(Y, 1, function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
})

# Identify genes with the most significant time-dependent model fit.
topgenes <- names(sort(gam.pval, decreasing = FALSE))[00:70]  

# Identify genes with the most significant time-dependent model fit.
topgenes <- names(sort(sce$slingPseudotime_1, decreasing = FALSE))[00:70]  

topgenes<-rownames(sorted)[1:70]
# Prepare and plot a heatmap of the top genes that vary their expression over pseudotime.

## best gam genes
library(scran)
library(scater)
on.first.path <- !is.na(neural.sce$slingPseudotime_1)
plotHeatmap(neural.sce[,on.first.path], order_columns_by=c("slingPseudotime_1","TSCAN.first","entropy"), 
    colour_columns_by=c("Phase","seurat_clusters"), features=topgenes,
    center=TRUE)

write.table(gam.pval,file="gam.tsv",row.names=T,sep="\t")


plotExpression(neural.sce, "SOX6", x = "slingPseudotime_1", 
               colour_by = "Phase", show_violin = FALSE,
               show_smooth = TRUE)+scale_color_tableau()


save(neural.sce,file="scegood.rda")