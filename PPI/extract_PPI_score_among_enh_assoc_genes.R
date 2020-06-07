#!/bin/Rscript
args<-commandArgs(TRUE)   #1st arument = list of enh associated genes in a tissue , 2nd = output directory path

library(STRINGdb)
library(igraph)

data = read.table(args[1], header = F)
colnames(data) = c("geneSymbol")

setwd(args[2])

cat("Initialising String object....\n")
string_db = STRINGdb$new( version="10", species=10090 , score_threshold = 900, input_directory="/NGS/users/Sid/StringDB_data")

cat("Mapping input genes to STRING IDs....\n")
d_map <- string_db$map( data, "geneSymbol", removeUnmappedRows = TRUE )
interactions = string_db$get_interactions(d_map$STRING_id)

scores = data.frame()
for( i in 1:nrow(d_map)){
  gene = d_map$geneSymbol[i]
  protein = d_map$STRING_id[i]
  sum = sum(interactions[interactions$from == protein | interactions$to == protein,16])
  scores[i,"Gene"] = gene
  scores[i,"Score"] = sum
}

out = sub(".all.genes.txt","",args[1])
file_name = paste(out,".ppiScore.txt",sep="")
write.table(scores, file = file_name, sep = "\t", col.names = T, row.names = F, quote = FALSE)
