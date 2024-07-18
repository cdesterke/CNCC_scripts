library(Seurat)
library(pals)

## load original Seurat object from GSE197513 dataset
data<-readRDS("GSE197513_CS20_craniofacial_snRNA.rds")

## dimension reduction in Seurat
DimPlot(data,reduction="tsne",group.by ="Idents",pt.size=1,cols=cols25(),label=F)
ElbowPlot(object = data,  ndims = 40)
data<-RunTSNE(data,dims=1:15,reduction="pca")

## metadata Seurat object
x<-data[[]]
head(x)

Idents(object = data) <- 'Idents'
Idents(object = data)
save(data,file="all.rda")

## find markers for each idents
markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
library(dplyr)
markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
library(dplyr)
markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top5

## doheatmap
library(ggplot2)
library(pals)
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
DoHeatmap(data, features = top5$gene,angle = 90,size = 1.5, group.colors=as.vector(cols25()))+
  scale_fill_gradientn(colours = rev(mapal))+ NoLegend()

write.table(markers,file="markersclusters.tsv",sep="\t",row.names=T)



## subset neural
neural<-subset(x = data, idents = "Neural Cells")
neural%>%NormalizeData()%>%FindVariableFeatures()%>%ScaleData()->neural

neural <- RunPCA(neural, features = VariableFeatures(object = neural))

DimPlot(neural,reduction="pca",group.by ="Phase",pt.size=1,cols=cols25(),label=F)
DimPlot(neural,reduction="pca",group.by ="Idents",pt.size=1,cols=cols25(),label=F)
ElbowPlot(object = data,  ndims = 40)

save(neural,file="neural.rda")

library(viridis)
library(RColorBrewer)
library(ggplot2)
FeaturePlot(neural,"HEY1",reduction = "pca",pt.size=2,min.cutoff = "q9") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "Spectral"))) + 
ggtitle("HEY1")


Idents(object = neural) <- 'seurat_clusters'
Idents(object = neural)


VlnPlot(neural, features = c("HEY1"), slot = "counts", log = TRUE,pt.size=1,split.by="Phase",
col=c("blue","red","green"))
ElbowPlot(object = neural,  ndims = 50)
neural <- FindNeighbors(object = neural, dims = 1:10)
neural <- FindClusters(object = neural)

DimPlot(neural,reduction="pca",pt.size=1,cols=cols25(),label=F)


markers <- FindAllMarkers(neural, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
library(dplyr)

library(dplyr)
markers %>%
    group_by(cluster) %>%
    top_n(n = 15, wt = avg_log2FC) -> top5

##custom doheatmap
library(ggplot2)
library(pals)
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
DoHeatmap(neural, features = top5$gene,angle = 90,size = 3, group.colors=as.vector(cols25()))+
  scale_fill_gradientn(colours = rev(mapal))+ NoLegend()

write.table(markers,file="markersneural4clusters.tsv",sep="\t",row.names=T)

Idents(object = neural) <- 'Phase'
Idents(object = neural)

markers <- FindAllMarkers(neural, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
library(dplyr)

library(dplyr)
markers %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = avg_log2FC) -> top5

##custom doheatmap
library(ggplot2)
library(pals)
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
DoHeatmap(neural, features = top5$gene,angle = 90,size = 3, group.colors=as.vector(cols25()))+
  scale_fill_gradientn(colours = rev(mapal))+ NoLegend()

write.table(markers,file="markersneuralCCphases.tsv",sep="\t",row.names=T)
