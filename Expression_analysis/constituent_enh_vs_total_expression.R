library(dplyr)
library(tidyr)
library(ggplot2)


# matrix of number of constituent enhancers for each gene
mat = read.table("/NGS/users/Sid/ENHANCERS/SuperEnh/AllGenes_distEnh/MATRIX_enhCount.txt", header=T, sep="\t", stringsAsFactor=F)
tissues = colnames(mat)[-1]

# file containing total-expression for each gene target
exp = read.table("/NGS/users/Sid/ENHANCERS/SuperEnh_27ac/Great_targets/Strict_targets/Expression_analysis/RESULTS.txt", header=F, sep="\t", stringsAsFactor=F)
colnames(exp) = c("Group","exp","Gene","tissue")


res = data.frame(Gene= character(0), exp= numeric(0), NumbEnh = numeric(0), Group = character(0))
for(i in 1:length(tissues)){
	tisName = tissues[i]
	# for exp
	tis.exp = exp[exp$tissue == tisName & exp$Group != "Absent" & exp$Group != "Weak",]
	data.exp = inner_join(tis.exp,mat)		
	#data <- gather(mat,tissue,enhCount, BAT:Wbrain, factor_key=TRUE)
	data.exp = subset(data.exp, select = c("Gene","exp",tisName,"Group"))
	colnames(data.exp) = c("Gene","exp","NumbEnh","Group")
	tisName = ifelse(tisName == "Es.E14","Es-E14",tisName)
	res = rbind(res,data.exp)
}
	
table(res$NumbEnh)
res = res[res$NumbEnh !=0,]
res$NumbEnh = ifelse(res$NumbEnh > 49, 50, res$NumbEnh)
r2.exp = round(cor(log10(res$exp+1),res$NumbEnh, method = "spearman"),2)
cor.test(log10(res$exp+1),res$NumbEnh, method = "spearman", exact=F)
text = paste("Spearman's rho =",r2.exp,sep="")
pval = "p < 2.2e-16" # outcome p-value is low, formatting it for the plot 

e.plot = ggplot(res, aes(x=factor(NumbEnh), y=log10(exp+1))) +
	 geom_boxplot(lwd = 0.35, width = 0.4, outlier.alpha=0.5, outlier.size = 0.12)+
	#stat_summary(fun.y=median, geom="line", aes(group=1), size = 0.1, color = "red", alpha = 1 )+
	geom_smooth(aes(group=1),size=0.4, method = "lm", fullrange=T, span = 10, color = "red") +
	theme_bw() +
	#theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
	theme(plot.title = element_text(size=12, hjust=0.5),
		legend.background = element_rect(fill = "transparent",colour = NA),
		legend.title=element_text(size=12),
		legend.text=element_text(size=11),
		legend.position="top",
		#legend.position= c(0.70,0.98),
		legend.key.size= unit(4,"mm"),
		axis.text.x = element_text(size=5),
		axis.text.y = element_text(size = 5),
		panel.border = element_rect(colour="BLACK",size=0.4),
		axis.title.x = element_text(size=9),
		axis.title.y = element_text(size=9,angle = 90),
		panel.background = element_rect(fill="transparent"),
		plot.background = element_rect(fill = "transparent",colour = NA)
	)+
	annotate("text",x=47,y=4.4,size=3,label= text, parse=F, color = "red") +
	annotate("text",x=47,y=4.2,size=3,label=pval, parse=T, color = "red") +
	scale_x_discrete(name = "Number of constituent enhancers", breaks = c(1:50), labels = c(1:49,"50+"))+
	scale_y_continuous(name = "Log10 expression")

file1 = "all.exp.png"
png(file1,bg="transparent",units="in",width = 6.25, height= 3.75 ,res=600)
plot(e.plot)
dev.off()
	
