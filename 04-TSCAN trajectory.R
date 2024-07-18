## TSCAN cell trajectory
library(scater)
plotReducedDim(neural.sce, "PCA", colour_by="Phase") 
colLabels(neural.sce)<-neural.sce$Phase


by.cluster <- aggregateAcrossCells(neural.sce, ids=colLabels(neural.sce))
centroids <- reducedDim(by.cluster, "PCA")



library(TSCAN)
mst <- createClusterMST(centroids, clusters=NULL)
mst

line.data <- reportEdges(by.cluster, mst=mst, clusters=NULL, use.dimred="PCA")



colData(neural.sce)


## extract pca from neural seurat
y=Embeddings(neural, reduction = "pca")
y<-as.data.frame(y)
neural.sce$PCA1<-y$PC_1
neural.sce$PCA2<-y$PC_2
plotReducedDim(neural.sce,"PCA", colour_by="Phase") + 
    geom_line(data=line.data,mapping=aes(x=PC_2, y= PC_2, group=edge,size=0.05,alpha=1))

map.tscan <- mapCellsToEdges(neural.sce, mst=mst, use.dimred="PCA")
tscan.pseudo <- orderCells(map.tscan, mst)
head(tscan.pseudo)


common.pseudo <- averagePseudotime(tscan.pseudo) 
plotReducedDim(neural.sce,"PCA", colour_by=I(common.pseudo), 
        text_by="label", text_colour="red") +
    geom_line(data=line.data, mapping=aes(x=PC_1, y=PC_2, group=edge))





pseudo <- testPseudotime(neural.sce, pseudotime=tscan.pseudo[,1])[[1]]
pseudo$SYMBOL <- row.names(pseudo)
pseudo[order(pseudo$p.value),]

neural.sce$TSCAN.first <- pathStat(tscan.pseudo)[,1]

sorted <- pseudo[order(pseudo$p.value),]


best <- head(row.names(sorted), 12)
plotExpression(neural.sce, features=best, x="TSCAN.first", colour_by="label")


sorted%>%filter(FDR<=0.05)->sorted
up.left <- sorted[sorted$logFC < 0,]
head(up.left, 12)

best <- head(row.names(up.left), 20)
plotExpression(neural.sce, features=best,     x="TSCAN.first", colour_by="label")



up.right <- sorted[sorted$logFC > 0,]
head(up.right)

best <- head(row.names(up.right), 10)
plotExpression(neural.sce, features=best,  x="TSCAN.first", colour_by="label")

write.csv(pseudo,file="pseudotime_commun.csv")


sorted <- pseudo[order(pseudo$p.value),]
sorted<-as.data.frame(sorted)

write.table(sorted,file="pseudotime_sorted.tsv",row.names=T,sep="\t")

save(neural.sce,file="pseudotime.sce.rda")

# Identify genes with the most significant time-dependent model fit.
topgenes <- names(sort(gam.pval, decreasing = FALSE))[00:70]  






## best gam genes
library(scran)
library(scater)
on.first.path <- !is.na(sce$slingPseudotime_1)
plotHeatmap(neural.sce[,on.first.path], order_columns_by="slingPseudotime_1", 
    colour_columns_by="Phase", features=topgenes,
    center=TRUE)