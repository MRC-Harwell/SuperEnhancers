#!/bin/Rscript
args<-commandArgs(TRUE)   #1st arument = table , 2nd = output directory path

library(rGREAT)

df = read.table(args[1],header=F,sep="\t")
colnames(df) = c("chr","start","end","id")

#df$start = df$start + 100

job = submitGreatJob(df, species = "mm9",  version = "3.0.0", request_interval = 60)
#job = submitGreatJob(df, species = "mm10",  version = "3.0.0", request_interval = 60)

tb = getEnrichmentTables(job, ontology = c("GO Biological Process", "GO Molecular Function", "Mouse Phenotype", "Human Phenotype", "Disease Ontology", "PANTHER Pathway", "MSigDB Pathway"))
names(tb)
mp = tb[[3]]
head(mp[order(mp$Binom_Raw_PValue),])


par(mfrow = c(1, 3))
res = plotRegionGeneAssociationGraphs(job, type = NULL)


bed <- data.frame(chr=seqnames(res),
  start=start(res),
  end=end(res),
  gene = elementMetadata(res)$gene,
  dist = elementMetadata(res)$distTSS
)

setwd(args[2])
out = sub(".bed","",args[1])
file_name = paste(out,".rGreat.txt",sep="")
write.table(bed,file=file_name,sep="\t",quote=F,row.names=F)
