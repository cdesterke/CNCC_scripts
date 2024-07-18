


library(Seurat)
library(GSVA)

## download CNCC differentiation wikipathways murine data WP2074
library(escape)
set<-read.csv("wp2074.csv",h=T)
geneset<-split(set$hs,set$geneset)
library(janitor)
library(dplyr)
geneset%>%clean_names()->geneset
library(SingleCellExperiment)
load("scegood.rda")

matrix <- escape.matrix(neural.sce,method = "GSVA",gene.sets = geneset, make.positive = T,
min.size = 0)
matrix<-as.data.frame(matrix)
matrix$entropy<-neural.sce$entropy

matrix$slingshot<-neural.sce$slingPseudotime_1
cormat <- round(cor(matrix),2)
cormat
library(reshape2)
melted_cormat <- melt(cormat)
head(melted_cormat)
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}



# Reorder the correlation matrix
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
library(ggplot2)
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1,size = 12, hjust = 1),axis.text.y = element_text( vjust = 1,size = 12, hjust = 1))+
  coord_fixed()+xlab("")+ylab("")
# Print the heatmap
print(ggheatmap)

save(matrix,file="matrix.rda")
colnames(matrix)
neural.sce$anatomic_nervous_system<-matrix$anatomic_nervous_system
neural.sce$chondrocytes<-matrix$chondrocytes
neural.sce$epidermis<-matrix$epidermis
neural.sce$melanocytes<-matrix$melanocytes
neural.sce$mesoderm<-matrix$mesoderm
neural.sce$migratory_cncc<-matrix$migratory_cncc
neural.sce$neural_plate_border<-matrix$neural_plate_border
neural.sce$premigratory_cncc<-matrix$premigratory_cncc

colData(neural.sce)



colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)



save(neural.sce,file="goodsceescape.rda")
library(scater)
plotReducedDim(neural.sce, "PCA", colour_by="migratory_cncc") 
colnames(matrix)
matrix$slingshot<-neural.sce$slingPseudotime_1


plotExpression(neural.sce,features="CTNNB1",x="Phase",colour_by="Phase")

plotColData(neural.sce, y = "chondrocytes", x ="seurat_clusters" , colour_by = "Phase")

library(ggbeeswarm)
library(ggthemes
plotColData(neural.sce, "entropy", x = "slingPseudotime_1", 
               colour_by = "Phase", show_violin = FALSE,
               show_smooth = TRUE)+scale_color_tableau()




