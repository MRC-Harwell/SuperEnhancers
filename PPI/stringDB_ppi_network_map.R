#!/bin/Rscript
args<-commandArgs(TRUE)   #1st arument = table/tissue , 2nd = output directory path

library(STRINGdb)
library(igraph)

# gene files
se_file = paste(args[1],".se.genes.txt",sep='')
te_file = paste(args[1],".te.uniq.genes.txt",sep='')

se = read.table(se_file, header = F)
colnames(se) = c("geneSymbol")
se$type = "SE"
te = read.table(te_file, header = F)
colnames(te) = c("geneSymbol")
te$type = "TE"
data = rbind(se,te)

setwd(args[2])

cat("Initialising String object....\n")
string_db = STRINGdb$new( version="10", species=10090 , score_threshold = 900, input_directory="/NGS/users/Sid/StringDB_data")

cat("Mapping input genes to STRING IDs....\n")
d_map <- string_db$map( data, "geneSymbol", removeUnmappedRows = TRUE )
#interactions = string_db$get_interactions(d_map$STRING_id)

cat("\nFetching the sub network....\n")
network.t = string_db$get_subnetwork(d_map$STRING_id)
network = delete.vertices(network.t,which(degree(network.t)<1))

summary = string_db$get_summary(d_map$STRING_id)
summary = sub("\n",", ",summary)
summary = sub("\n",", ",summary)
summary = sub("\n","",summary)

### extracting gene names and adding as labels to the i graph
cat("Extracting genes names....\n")
names = V(network)$name
labels = data.frame()
for(i in 1:length(names)){
	l = d_map[d_map$STRING_id == names[i],,drop=T]
	labels[i,"geneSymbol"] = l$geneSymbol[1]
	labels[i,"STRING_id"] = l$STRING_id[1]
	#gene = l$geneSymbol[1]
	#type = data[data$geneSymbol == gene,]$type
	#labels[i,"type"] = type
	labels[i,"type"] = l$type[1]
}
### function to convert gene symbols UPPERCASE to Lowercase
r_ucfirst <- function (str) {
  paste(toupper(substring(str, 1, 1)), tolower(substring(str, 2)), sep = "")
}
labels$geneSymbol = r_ucfirst(labels$geneSymbol)
V(network)$label = labels$geneSymbol
V(network)$type = labels$type

### working on extracting phenotype, for adding attributes to the graph object
cat("Extracting the phenotypes....\n")
mp = "/NGS/users/Sid/ChromHmm/Posterior/State_6/POSTERIOR/Tau_enhancers/TauResults/Highly_specific/Pleitropy_mp/mgi_mp_matrix_update.txt"
df = read.table(mp,header=T,sep="\t")

phenotypes = data.frame()
for(i in 1:nrow(labels)){
	present = labels$geneSymbol[i] %in% df$Gene
	if(present == FALSE){
		phenotypes[i,"gene"] = labels$geneSymbol[i]
		phenotypes[i,"mp1"] = 0
		phenotypes[i,"mp2"] = 0
	}else{
		m = df[df$Gene == labels$geneSymbol[i],,drop=T]
		phenotypes[i,"gene"] = labels$geneSymbol[i]
		phenotypes[i,"mp1"] = m$MP.0005386[1]
		phenotypes[i,"mp2"] = m$MP.0003631[1]
	}
}
###### for single phenotype
#phenotypes$final = phenotypes$mp1
###### for multiple phenotypes
phenotypes$final = phenotypes$mp1 + phenotypes$mp2

V(network)$mp = phenotypes$final

out = sub(".genes.txt","",args[1])
file_name = paste(out,".png",sep="")
title = out
print(file_name)

#legend = "Adipose tissue phenotype" #5375
#legend = "Immune system and\nhematopoietic phenotype" #5387,5397
legend = "Behavioural/neurological and\nnervous system phenotype" #5386,3631
#legend = "Embryogenesis phenotype" #5380
#legend = "Cardiovascular system\nphenotype" #5385
#legend = "Renal/urinary system phenotype" #5367
#legend = "Limbs/digits/tail phenotype" #5371
#legend = "Liver/biliary system\nphenotype" #5370
#legend = "Respiratory system phenotype" #5388
#legend = "Digestive/alimentary\nphenotype" #5381
#legend = "Reproductive system\nphenotype" #5389

cat("Plotting network....\n")
png(file_name,bg="transparent",units="in",width = 7.25, height= 6.25 ,res=600)
lo = layout_with_fr(network, grid = "nogrid")
plot(network,
	vertex.label = NA,
	#vertex.size = 3.5,
	#layout=layout.kamada.kawai,
	vertex.size = 2.5,
	#vertex.size = 3.2,
	layout = lo,
	edge.width = 0.4,
	frame=TRUE,
	margin = c(0,0,0,0),
	main = title,
	#vertex.frame.color = adjustcolor("black",alpha.f=0.4),
	vertex.shape = c("circle","square")[1+(V(network)$type=="SE")],
	vertex.frame.color = adjustcolor(c("black","red")[1+(V(network)$type=="SE")],alpha.f=0.4),
	vertex.color=adjustcolor(c("hotpink","skyblue")[1+(V(network)$mp=="0")],alpha.f=0.8)
)
legend(x=0.76, y=-0.85, c("super-enh","typical-enh","known","novel"), pch=c(22,21,21,21), pt.bg=c("NA","NA","hotpink","skyblue"),col=c("red","black","black","black"),pt.cex=1, cex=.6, bty="n", ncol=2, title = legend)
text(0,1.15,summary, cex=0.7)
dev.off()


cat("Extracting all connected nodes....\n")
genesWithInteractions = V(network)$label
allConnecting = data.frame()
for(i in 1:length(genesWithInteractions)){
	l = phenotypes[phenotypes$gene == genesWithInteractions[i],,drop=T]
	allConnecting[i,"gene"] = l$gene[1]
	allConnecting[i,"mpFinal"] = l$final[1]
}

### extracting novel genes connecting to known phenotype genes
cat("Extracting novel nodes connected to phenotype nodes....\n")
edges_nk = as_ids(E(network)[ V(network)[mp=="0"] %->% V(network)[mp >= "1"] ])
k=1
novelConnecting = data.frame()
for(i in 1:length(edges_nk)){
	for(j in 1:2){
		n1 = labels[labels$STRING_id == unlist(strsplit(edges_nk[i],"[|]"))[j],,drop=T]$geneSymbol
		l = phenotypes[phenotypes$gene == n1,,drop=T]
		if(l$final == 0){
			novelConnecting[k,"gene"] = l$gene[1]
			novelConnecting[k,"mpFinal"] = l$final[1]
			k=k+1
		}
	}
}
novelConnecting = unique(novelConnecting)

file_name_2 = paste(out,"_allConnecting.txt",sep="")
write.table(allConnecting,file=file_name_2,sep="\t",quote=F,row.names=F)
file_name_3 = paste(out,"_novelConnecting.txt",sep="")
write.table(novelConnecting,file=file_name_3,sep="\t",quote=F,row.names=F)

### extracting no of genes with phenotype from original input data
cat("Calculating and writing stats....\n")
input_pheno = data.frame()
for(i in 1:nrow(data)){
	present = data$geneSymbol[i] %in% df$Gene
	if(present == FALSE){
		input_pheno[i,"mp1"] = 0
		input_pheno[i,"mp2"] = 0
	}else{
		m = df[df$Gene == as.character(data$geneSymbol[i]),,drop=T]
		input_pheno[i,"mp1"] = m$MP.0005386[1]
		input_pheno[i,"mp2"] = m$MP.0003631[1]
	}
}

input_genes = table(input_pheno)
network_genes = table(phenotypes$final)
stats_1 = rbind(input_genes,network_genes)
statsFile = paste(out,".stats.txt",sep="")
write.table(stats_1,file=statsFile,sep="\t",quote=F)

cat("writing mapping and Novel-known interactions....\n")
file_name_4 = paste(out,"_mapping.txt",sep="")
write.table(d_map,file=file_name_4,sep="\t",quote=F,row.names=F)
file_name_5 = paste(out,"_nkinteractions.txt",sep="")
write.table(edges_nk,file=file_name_5,sep="\t",quote=F,row.names=F)
