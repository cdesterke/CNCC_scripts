list.files()

## download EpiGenes database
epi<-read.csv("EpiGenes_main.csv",h=T)

epi[epi=="#"]<-NA


library(dplyr)
epi%>%select(HGNC_symbol,Function)->epi

colnames(epi)<-c("gene","function")

save(epi,file="epi.rda")

## load pseudotime table results
pseudo<-read.table("pseudotime_sorted.tsv",h=T,sep="\t")

pseudo$gene<-row.names(pseudo)
pseudo%>%relocate(gene)->pseudo
pseudo%>%filter(FDR<=0.01)->sig

sig%>%inner_join(epi,by="gene")->df

library(writexl)
write_xlsx(df,"epifactorpseudo.xlsx")
save(df,file="epifactorpseudo.rda")

topgenes<-df$gene
save(topgenes,file="topgenes.rda")




