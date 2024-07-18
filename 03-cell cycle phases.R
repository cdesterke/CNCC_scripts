## barplots metadata neural cells

x<-neural[[]]

library(ggplot2)
library(pals)
p=ggplot(x,aes(x=factor(seurat_clusters),fill=factor(Phase)))+
    geom_bar(position="fill")+
    #geom_text(aes(label=..count..),stat='count',position=position_fill(vjust=0.5),size=4,colour="white")+
#    scale_fill_brewer(palette="Dark2")+
	scale_fill_manual(values=as.vector(cols25()))+
    labs(fill = "orig.ident")+  
    ggtitle("Cell Stratification") +
    xlab("Groups of cells") + ylab("relative proportions ")+theme_classic()


p+coord_flip() 