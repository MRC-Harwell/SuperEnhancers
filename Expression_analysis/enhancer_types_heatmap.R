#!/bin/Rscript
args<-commandArgs(TRUE)   #1st arument = table , 2nd = output directory path
setwd(args[2])

library(pheatmap)

name = sub(".txt","",args[1])
file_name = paste(name,".png",sep="")
df = read.table(args[1], header=T, sep="\t")
rownames(df) = df$Gene
df = df[,-1]

df$sums = rowSums(df)

# write cell type specific table for each geneSymbol
table = df[,23,drop=FALSE]
table$Type = "Super"
write.table(table, file = "all.te.genes.cellTypeSpec.txt", sep="\t",row.names=T, col.names=F, quote=F)


############################################################
################ For Super enhancers #######################
############################################################

c1 = df[df$sums == 1,c(1:22)]
c2 = df[df$sums == 2,c(1:22)]
c3 = df[df$sums >= 3,c(1:22)]

s1 = c1[order(c1[,1],c1[,2],c1[,3],c1[,4],c1[,5],c1[,6],c1[7],c1[,8],c1[,9],c1[,10],c1[,11],c1[,12],c1[,13],c1[,14],c1[,15],c1[,16],c1[,17],c1[,18],c1[,19],c1[,20],c1[,21],c1[,22],decreasing=TRUE),]
s2 = c2[order(c2[,1],c2[,2],c2[,3],c2[,4],c2[,5],c2[,6],c2[7],c2[,8],c2[,9],c2[,10],c2[,11],c2[,12],c2[,13],c2[,14],c2[,15],c2[,16],c2[,17],c2[,18],c2[,19],c2[,20],c2[,21],c2[,22],decreasing=TRUE),]
s3 = c3[order(c3[,1],c3[,2],c3[,3],c3[,4],c3[,5],c3[,6],c3[7],c3[,8],c3[,9],c3[,10],c3[,11],c3[,12],c3[,13],c3[,14],c3[,15],c3[,16],c3[,17],c3[,18],c3[,19],c3[,20],c3[,21],c3[,22],decreasing=TRUE),]


g1 = nrow(c1)
g2 = g1 + nrow(c2)
gaps = c(g1,g2)

mat = rbind(s1,s2,s3)
mat=data.matrix(mat)

png(file_name,bg="transparent",units="in",width = 7, height= 10 ,res=600)
pheatmap(
  mat,
  color = colorRampPalette(c("black","forestgreen"))(50),
  #border_color = "white",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = F,
  show_colnames = T,
  main = NA,
  fontsize = 15,
  fontsize_row = 2, #10
  fontsize_col = 18, #13
  treeheight_row = 0,
  annotation_legend = FALSE,
  legend = FALSE,
  gaps_row = gaps
)
dev.off()

