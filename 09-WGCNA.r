library(SingleCellExperiment)
library(Seurat)
list.files()
library(TrajectoryUtils)
load("neural.rda")

meta<-neural[[]]

meta$slingshot<-neural.sce$slingPseudotime_1
meta$entropy<-neural.sce$entropy
meta$epidermis<-neural.sce$epidermis
meta$chondrocyte<-neural.sce$chondrocyte
meta$migratory_cncc<-neural.sce$migratory_cncc
meta$neural_plate_border<-neural.sce$neural_plate_border
meta$anatomic_nervous_system<-neural.sce$anatomic_nervous_system
meta$premigratory_cncc<-neural.sce$premigratory_cncc

meta$old.ident<-NULL
meta$orig.ident<-NULL
meta$RNA_snn_res.0.1<-NULL
meta$Idents<-NULL
meta$seurat_clusters<-NULL
str(meta)
library(fastDummies)
results<-fastDummies::dummy_cols(meta,remove_selected_columns=T)
head(results)
results<-as.data.frame(results)
row.names(results)<-row.names(meta)

save(results,file="wgcna_annot.rda")

data<-neural[["RNA"]]$data
dim(data)

load("topgenesEPI.rda")
epi<-topgenes
rm(topgenes)
load("topgenesTF.rda")
tf<-topgenes
rm(topgenes)

conca<-c(tf,epi)
length(conca)
conca<-unique(conca)
inter<-intersect(conca,row.names(data))

small<-data[row.names(data)%in%inter,]
trans<-t(small)

all(row.names(trans)==row.names(results))

save(trans,file="wgcna_data.rda")


library(WGCNA)
pheno<-results
matrix<-as.matrix(trans)

# Supprimez les gènes et les échantillons avec trop de valeurs manquantes
gsg <- goodSamplesGenes(matrix, verbose = 3)
gsg$allOK 
if (!gsg$allOK) {
  matrix <- matrix[gsg$goodSamples, gsg$goodGenes]
}
#matrix<-scale(matrix)
# Cluster des échantillons pour détecter les outliers
sampleTree <- hclust(dist(matrix), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")


save(matrix, pheno, file="SamplesAndTraits.RData")

# Choose a soft threshold power- USE A SUPERCOMPUTER IRL ------------------------------------

powers = c(c(1:10), seq(from =4, to=30, by=2)) #choosing a set of soft-thresholding powers
sft = pickSoftThreshold(matrix, powerVector=powers, verbose =5, networkType="signed") #call network topology analysis function

sizeGrWindow(9,5)
par(mfrow= c(1,2))
cex1=0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type= "n", main= paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
abline(h=0.90, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")


#build a adjacency "correlation" matrix
enableWGCNAThreads()
softPower =5
adjacency = adjacency(matrix, power = softPower, type = "signed") #specify network type
head(adjacency)

# Construct Networks- USE A SUPERCOMPUTER IRL -----------------------------
#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency, TOMType="signed") # specify network type
dissTOM = 1-TOM

# Generate Modules --------------------------------------------------------
# Generate Modules --------------------------------------------------------

library(flashClust)
# Generate a clustered gene tree
geneTree = flashClust(as.dist(dissTOM), method="average")
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)
#This sets the minimum number of genes to cluster into a module
minModuleSize = 15
dynamicMods = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=4, pamRespectsDendro= FALSE, minClusterSize = minModuleSize)
dynamicColors= labels2colors(dynamicMods)
MEList= moduleEigengenes(matrix, colors= dynamicColors,softPower = softPower)
MEs= MEList$eigengenes
MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")
save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_allSamples_signed_RLDfiltered.RData")
table(dynamicMods)

TOMplot(dissTOM,geneTree,dynamicColors)

## mds plot
cmd1=cmdscale(as.dist(dissTOM),2)

plot(cmd1,col=dynamicColors,main="MDS plot",xlab="Scaling Dimension 1",ylab="Scaling Dimension 2")


#plots tree showing how the eigengenes cluster together
#INCLUE THE NEXT LINE TO SAVE TO FILE
MEDissThres = 0.0
merge = mergeCloseModules(matrix, dynamicColors, cutHeight= MEDissThres, verbose =3)
mergedColors = merge$colors
mergedMEs = merge$newMEs


#plot dendrogram with module colors below it
#INCLUE THE NEXT LINE TO SAVE TO FILE
#pdf(file="cluster.pdf")
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
#INCLUE THE NEXT LINE TO SAVE TO FILE
#dev.off()



#Define number of genes and samples
nGenes = ncol(matrix)
nSamples = nrow(matrix)
#Recalculate MEs with color labels
MEs0 = moduleEigengenes(matrix, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, pheno, use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#Print correlation heatmap between modules and traits
textMatrix= paste(signif(moduleTraitCor, 2), "\n(",
                  signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)


sizeGrWindow(12, 6)
#par(mar= c(8, 8, 3, 3))
labeledHeatmap(Matrix= moduleTraitCor,
               xLabels= names(pheno),
               yLabels= names(MEs),
               ySymbols= names(MEs),
               colorLabels= FALSE,
               colors= blueWhiteRed(50),
               textMatrix= textMatrix,
               setStdMargins= FALSE,
               cex.text= 1,
               zlim= c(-1,1),
               main= paste("Module-trait relationships"))



## gene names from yellow modulle
yellow<-colnames(matrix)[moduleColors=="yellow"]
yellow<-data.frame(colnames(matrix)[moduleColors=="yellow"])
write.table(yellow,file="module_yellow.txt")

## gene names from blue modulle
blue<-colnames(matrix)[moduleColors=="blue"]
blue<-data.frame(colnames(matrix)[moduleColors=="blue"])
write.table(blue,file="module_blue.txt")

## gene names from turquoise modulle
turquoise<-colnames(matrix)[moduleColors=="turquoise"]
turquoise<-data.frame(colnames(matrix)[moduleColors=="turquoise"])
write.table(turquoise,file="module_turquoise.txt")

## gene names from brown modulle
brown<-colnames(matrix)[moduleColors=="brown"]
brown<-data.frame(colnames(matrix)[moduleColors=="brown"])
write.table(brown,file="module_brown.txt")

dftf<-as.data.frame(tf)
dftf$group<-"TF"
colnames(dftf)<-c("gene","group")

dfepi<-as.data.frame(epi)
dfepi$group<-"EpiFactor"
colnames(dfepi)<-c("gene","group")

all<-rbind(dftf,dfepi)
library(dplyr)
all%>%filter(gene%in%yellow)->cluster_yellow
all%>%filter(gene%in%blue)->cluster_blue
all%>%filter(gene%in%brown)->cluster_brown
all%>%filter(gene%in%turquoise)->cluster_turquoise

cluster_yellow$module<-"yellow"
cluster_blue$module<-"blue"
cluster_brown$module<-"brown"
cluster_turquoise$module<-"turquoise"

final<-rbind(cluster_yellow,cluster_blue,cluster_brown,cluster_turquoise)


library(ggplot2)
p=ggplot(final,aes(x=factor(module),fill=factor(group)))+
    geom_bar(position="fill")+
    geom_text(aes(label=..count..),stat='count',position=position_fill(vjust=0.5),size=6,colour="white")+
    scale_fill_brewer(palette="Dark2")+
    labs(fill = "")+  
    ggtitle("Module Stratification") +
    xlab("WGCNA modules") + ylab("relative proportions ")+theme_classic()

## custom size label == remove classic theme
p+theme(legend.position="bottom")
p+theme(plot.title = element_text(color="black", size=14, face="bold"),
    	axis.title.x = element_text(color="black", size=14, face="bold"),
    	axis.title.y = element_text(color="black", size=14, face="bold"))
p+coord_flip()  

write.table(final,file="moduleprogram.tsv",row.names=F,sep="\t")




# Calculer la connectivité intramodulaire
moduleColors <- dynamicColors
geneModuleMembership <- as.data.frame(cor(matrix, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
write.csv(geneModuleMembership,file="genemodulecor.csv",row.names=T)
write.csv(MMPvalue,file="genemodulepvalue.csv",row.names=T)

# Calculer la corrélation entre gènes et traits
geneTraitSignificance <- as.data.frame(cor(matrix, pheno, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
write.csv(geneTraitSignificance,file="genetraitcor.csv",row.names=T)
write.csv(GSPvalue,file="genetraitpvalue.csv",row.names=T)


# Identifier les gènes hubs
for (module in unique(moduleColors)) {
	df<-geneModuleMembership
	colnames(df)<-sub("^[^ME]*ME","",colnames(df))
	moduleGenes <- (moduleColors == module)
	cat("\nModule: ", module)
	cat("\nTop hub genes based on module membership and trait significance:\n")
 	hubGenes <- rownames(t(matrix))[moduleGenes]
	df<-df[row.names(df)%in%hubGenes,] 
	print(df%>%arrange(desc(!!sym(module))))
}



## connectivity for yellow modulle
library(WGCNA)
library(igraph)

# Exemple de matrice d'expression (genes x samples) et modules (dynamicColors)
# Assurez-vous que exprs_matrix et dynamicColors sont disponibles dans votre environnement
# exprs_matrix <- ...
# dynamicColors <- ...

# Choisir le module d'intérêt
module <- "yellow"
moduleGenes <- (moduleColors == module)
moduleGenes <- rownames(t(matrix))[moduleGenes]
moduleExprs <- matrix[,colnames(matrix)%in%moduleGenes]

# Calculer la similarité topologique (TOM)
TOM <- TOMsimilarityFromExpr(moduleExprs, power = 4)
names<-as.list(moduleGenes)
row.names(TOM) <- names
colnames(TOM)<- names

moduleTOM <- TOM[moduleGenes, moduleGenes]


# Convertir la matrice de similarité en un objet igraph
moduleGraph <- graph_from_adjacency_matrix(moduleTOM, mode = "undirected", weighted = TRUE, diag = FALSE)
moduleGraph <- simplify(moduleGraph, remove.multiple = TRUE, remove.loops = TRUE)

# Calculer la connectivité des gènes
geneConnectivity <- colSums(moduleTOM) - 1
topGenes <- sort(geneConnectivity, decreasing = TRUE)
print(head(topGenes))

# Ajouter des attributs de nœud
V(moduleGraph)$name <- colnames(moduleExprs)
V(moduleGraph)$size <- geneConnectivity * 5 / max(geneConnectivity)

# Visualiser le réseau avec des options d'igraph
plot(moduleGraph, vertex.label = V(moduleGraph)$name, vertex.label.cex = 1,vertex.label.dist=1, 
vertex.size = V(moduleGraph)$size, edge.width = E(moduleGraph)$weight * 2, 
layout = layout_nicely)
layout.fruchterman.reingold)

###blue
# Choisir le module d'intérêt
module <- "blue"
moduleGenes <- (moduleColors == module)
moduleGenes <- rownames(t(matrix))[moduleGenes]
moduleExprs <- matrix[,colnames(matrix)%in%moduleGenes]

# Calculer la similarité topologique (TOM)
TOM <- TOMsimilarityFromExpr(moduleExprs, power = 4)
names<-as.list(moduleGenes)
row.names(TOM) <- names
colnames(TOM)<- names

moduleTOM <- TOM[moduleGenes, moduleGenes]


# Convertir la matrice de similarité en un objet igraph
moduleGraph <- graph_from_adjacency_matrix(moduleTOM, mode = "undirected", weighted = TRUE, diag = FALSE)
moduleGraph <- simplify(moduleGraph, remove.multiple = TRUE, remove.loops = TRUE)

# Calculer la connectivité des gènes
geneConnectivity <- colSums(moduleTOM) - 1
topGenes <- sort(geneConnectivity, decreasing = TRUE)
print(head(topGenes))

# Ajouter des attributs de nœud
V(moduleGraph)$name <- colnames(moduleExprs)
V(moduleGraph)$size <- geneConnectivity * 5 / max(geneConnectivity)

# Visualiser le réseau avec des options d'igraph
plot(moduleGraph, vertex.label = V(moduleGraph)$name, vertex.label.cex = 1,vertex.label.dist=1, 
vertex.size = V(moduleGraph)$size, edge.width = E(moduleGraph)$weight * 2, 
layout = layout_nicely)



###turquoise
# Choisir le module d'intérêt
module <- "turquoise"
moduleGenes <- (moduleColors == module)
moduleGenes <- rownames(t(matrix))[moduleGenes]
moduleExprs <- matrix[,colnames(matrix)%in%moduleGenes]

# Calculer la similarité topologique (TOM)
TOM <- TOMsimilarityFromExpr(moduleExprs, power = 4)
names<-as.list(moduleGenes)
row.names(TOM) <- names
colnames(TOM)<- names

moduleTOM <- TOM[moduleGenes, moduleGenes]


# Convertir la matrice de similarité en un objet igraph
moduleGraph <- graph_from_adjacency_matrix(moduleTOM, mode = "undirected", weighted = TRUE, diag = FALSE)
moduleGraph <- simplify(moduleGraph, remove.multiple = TRUE, remove.loops = TRUE)

# Calculer la connectivité des gènes
geneConnectivity <- colSums(moduleTOM) - 1
topGenes <- sort(geneConnectivity, decreasing = TRUE)
print(head(topGenes))

# Ajouter des attributs de nœud
V(moduleGraph)$name <- colnames(moduleExprs)
V(moduleGraph)$size <- geneConnectivity * 5 / max(geneConnectivity)

# Visualiser le réseau avec des options d'igraph
plot(moduleGraph, vertex.label = V(moduleGraph)$name, vertex.label.cex = 1,vertex.label.dist=1, 
vertex.size = V(moduleGraph)$size, edge.width = E(moduleGraph)$weight * 2, 
layout = layout_nicely)


##top 40 turquoise module
topGenes <- names(sort(geneConnectivity, decreasing = TRUE)[1:40])

# Construire le réseau de co-expression pour les gènes les plus connectés
topModuleTOM <- moduleTOM[topGenes, topGenes]
topModuleGraph <- graph_from_adjacency_matrix(topModuleTOM, mode = "undirected", weighted = TRUE, diag = FALSE)
topModuleGraph <- simplify(topModuleGraph, remove.multiple = TRUE, remove.loops = TRUE)

# Ajouter des attributs de nœud
V(topModuleGraph)$name <- topGenes
V(topModuleGraph)$size <- geneConnectivity[topGenes] * 5 / max(geneConnectivity[topGenes])

# Visualiser le réseau avec des options d'igraph
# Visualiser le réseau avec des options d'igraph
plot(topModuleGraph, vertex.label = V(topModuleGraph)$name, vertex.label.cex = 1, vertex.size = V(topModuleGraph)$size, 
edge.width = E(topModuleGraph)$weight * 2, layout = layout_nicely)




###brown
# Choisir le module d'intérêt
module <- "brown"
moduleGenes <- (moduleColors == module)
moduleGenes <- rownames(t(matrix))[moduleGenes]
moduleExprs <- matrix[,colnames(matrix)%in%moduleGenes]

# Calculer la similarité topologique (TOM)
TOM <- TOMsimilarityFromExpr(moduleExprs, power = 4)
names<-as.list(moduleGenes)
row.names(TOM) <- names
colnames(TOM)<- names

moduleTOM <- TOM[moduleGenes, moduleGenes]


# Convertir la matrice de similarité en un objet igraph
moduleGraph <- graph_from_adjacency_matrix(moduleTOM, mode = "undirected", weighted = TRUE, diag = FALSE)
moduleGraph <- simplify(moduleGraph, remove.multiple = TRUE, remove.loops = TRUE)

# Calculer la connectivité des gènes
geneConnectivity <- colSums(moduleTOM) - 1
topGenes <- sort(geneConnectivity, decreasing = TRUE)
print(head(topGenes))

# Ajouter des attributs de nœud
V(moduleGraph)$name <- colnames(moduleExprs)
V(moduleGraph)$size <- geneConnectivity * 5 / max(geneConnectivity)

# Visualiser le réseau avec des options d'igraph
plot(moduleGraph, vertex.label = V(moduleGraph)$name, vertex.label.cex = 1,vertex.label.dist=1, 
vertex.size = V(moduleGraph)$size, edge.width = E(moduleGraph)$weight * 2, 
layout = layout_nicely)


topyellow<-as.data.frame(topyellow)
topyellow$group<-"yellow"
colnames(topyellow)<-c("connectivity","network")


topblue<-as.data.frame(topblue)
topblue$group<-"blue"
colnames(topblue)<-c("connectivity","network")

topbrown<-as.data.frame(topbrown)
topbrown$group<-"brown"
colnames(topbrown)<-c("connectivity","network")


topturquoise<-as.data.frame(topturquoise)
topturquoise$group<-"turquoise"
colnames(topturquoise)<-c("connectivity","network")

output<-rbind(topyellow,topblue,topbrown,topturquoise)
output$gene<-row.names(output)
all%>%left_join(output,by="gene")%>%arrange(desc(connectivity))->connectivity

library(ggpubr)
connectivity%>%filter(network=="yellow")->yellow
ggdensity(yellow, x = "connectivity",
   add = "mean", rug = TRUE,
   color = "group", palette = c("plum2", "turquoise"))



connectivity%>%filter(network=="turquoise")->turquoise
ggdensity(turquoise, x = "connectivity",
   add = "mean", rug = TRUE,
   color = "group", palette = c("plum2", "turquoise"))




connectivity%>%filter(network=="blue")->blue
ggdensity(blue, x = "connectivity",
   add = "mean", rug = TRUE,
   color = "group", palette = c("plum2", "turquoise"))



connectivity%>%filter(network=="brown")->brown
ggdensity(brown, x = "connectivity",
   add = "mean", rug = TRUE,
   color = "group", palette = c("plum2", "turquoise"))


write.table(connectivity,file="connectivity.tsv",row.names=F,sep="\t")


epi<-read.csv("epi.csv",h=T)
epi<-epi[complete.cases(epi$group),]

epi%>%inner_join(connectivity,by="gene")%>%filter(group.y=="EpiFactor")->epi

library(pals)
library(ggplot2)
p=ggplot(epi,aes(x=factor(network),fill=factor(group.x)))+
    geom_bar(position="fill")+
    geom_text(aes(label=..count..),stat='count',position=position_fill(vjust=0.5),size=4,colour="white")+
    scale_fill_manual(values=as.vector(cols25()))+
    labs(fill = "")+  
    ggtitle("Module Stratification") +
    xlab("WGCNA modules") + ylab("relative proportions ")+theme_classic()

## custom size label == remove classic theme

p+theme(plot.title = element_text(color="black", size=14, face="bold"),
    	axis.title.x = element_text(color="black", size=14, face="bold"),
    	axis.title.y = element_text(color="black", size=14, face="bold"))
p  

write.table(epi,file="epiconnectivity.tsv",row.names=F,sep="\t")


tf<-read.csv("tfpseudo.csv",h=T)
tf<-tf[complete.cases(tf$group),]
table(tf$group)
tf%>%inner_join(connectivity,by="gene")%>%filter(group.y=="TF")->tf


## tf yellow
tf%>%filter(network=="yellow")->yellow
yellow%>%mutate(threshold=ifelse(connectivity>=0.186,"up","bottom"))->yellow


library(ggplot2)
p=ggplot(yellow,aes(x=factor(threshold),fill=factor(group.x)))+
    geom_bar(position="fill")+
    geom_text(aes(label=..count..),stat='count',position=position_fill(vjust=0.5),size=4,colour="white")+
    scale_fill_manual(values=as.vector(cols25()))+
    labs(fill = "DBD")+  
    ggtitle("Yellow Module TF") +
    xlab("WGCNA connectivity") + ylab("relative proportions ")+theme_classic()
p+coord_flip()+theme(legend.position="bottom")




## tf blue
tf%>%filter(network=="blue")->blue
blue%>%mutate(threshold=ifelse(connectivity>=0.239,"up","bottom"))->blue


library(ggplot2)
p=ggplot(blue,aes(x=factor(threshold),fill=factor(group.x)))+
    geom_bar(position="fill")+
    geom_text(aes(label=..count..),stat='count',position=position_fill(vjust=0.5),size=4,colour="white")+
    scale_fill_manual(values=as.vector(cols25()))+
    labs(fill = "DBD")+  
    ggtitle("Blue Module TF") +
    xlab("WGCNA connectivity") + ylab("relative proportions ")+theme_classic()
p+coord_flip()+theme(legend.position="bottom")




## epi brown
epi%>%filter(network=="brown")->brown
brown%>%mutate(threshold=ifelse(connectivity>= 0.432,"up","bottom"))->brown


library(ggplot2)
p=ggplot(brown,aes(x=factor(threshold),fill=factor(group.x)))+
    geom_bar(position="fill")+
    geom_text(aes(label=..count..),stat='count',position=position_fill(vjust=0.5),size=4,colour="white")+
    scale_fill_manual(values=as.vector(cols25()))+
    labs(fill = "")+  
    ggtitle("Brown Module EpiFactors") +
    xlab("WGCNA connectivity") + ylab("relative proportions ")+theme_classic()
p




## epi turquoise
epi%>%filter(network=="turquoise")->turquoise
turquoise%>%mutate(threshold=ifelse(connectivity>= 0.185,"up","bottom"))->turquoise


library(ggplot2)
p=ggplot(turquoise,aes(x=factor(threshold),fill=factor(group.x)))+
    geom_bar(position="fill")+
    geom_text(aes(label=..count..),stat='count',position=position_fill(vjust=0.5),size=4,colour="white")+
    scale_fill_manual(values=as.vector(cols25()))+
    labs(fill = "")+  
    ggtitle("Turquoise Module EpiFactors") +
    xlab("WGCNA connectivity") + ylab("relative proportions ")+theme_classic()
p


## tf turquoise
tf%>%filter(network=="turquoise")->turquoise2
turquoise2%>%mutate(threshold=ifelse(connectivity>= 0.301,"up","bottom"))->turquoise2


library(ggplot2)
p=ggplot(turquoise2,aes(x=factor(threshold),fill=factor(group.x)))+
    geom_bar(position="fill")+
    geom_text(aes(label=..count..),stat='count',position=position_fill(vjust=0.5),size=4,colour="white")+
    scale_fill_manual(values=as.vector(cols25()))+
    labs(fill = "DBD")+  
    ggtitle("Turquoise Module TF") +
    xlab("WGCNA connectivity") + ylab("relative proportions ")+theme_classic()
p

write.table(tf,file="tfconnectivity.tsv",row.names=F,sep="\t")


