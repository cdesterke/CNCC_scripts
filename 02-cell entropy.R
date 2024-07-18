## single cell experiment
library(SingleCellExperiment)
neural.sce <- as.SingleCellExperiment(neural)

## cell entropy
library(TSCAN)
entropy<-perCellEntropy(neural.sce)


ent.data<-data.frame(cluster=neural$Phase,entropy=entropy)
head(ent.data)
library(dplyr)


ent.data$code<-as.numeric(ent.data$code)
ent.data$cluster<-as.factor(ent.data$Phase)
library(ggplot2)

ggplot(ent.data,aes(reorder(cluster,cluster),entropy))+
geom_violin(aes(fill = cluster), trim = FALSE) + 
geom_boxplot(width = 0.4,outlier.shape = NA)+
scale_fill_manual(values = c("blue", "red","green"))+
	#scale_fill_brewer(palette="Set1")+
  #geom_point(aes(fill=factor(cluster),size=0.5),shape = 21, alpha = .8, position = position_dodge2(width = .5))+
  theme_classic(base_size=18) +
  theme(legend.position = "none")+xlab("Clusters of CNCC")+ggtitle("Cell entropy")+ 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

neural$entropy<-ent.data$entropy
