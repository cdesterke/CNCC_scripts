
## download human transciption factor database
tf<-read.csv("tf.csv",h=T)


tf%>%select(HGNC.symbol,DBD,Is.TF.)%>%filter(Is.TF.=="Yes")->tf

tf%>%dplyr::rename(gene="HGNC.symbol")->tf


## load pseudotime table results
pseudo<-read.table("pseudotime_sorted.tsv",h=T,sep="\t")

pseudo$gene<-row.names(pseudo)
pseudo%>%relocate(gene)->pseudo
pseudo%>%filter(FDR<=0.01)->sig

sig%>%inner_join(tf,by="gene")->df


library(writexl)
write_xlsx(df,"tfpseudo.xlsx")
save(df,file="tfpseudo.rda")

topgenes<-df$gene
save(topgenes,file="topgenes.rda")