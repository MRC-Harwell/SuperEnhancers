library(dplyr)
library(tidyr)
library(ggplot2)


# matrix of number of constituent enhancers for each gene
mat = read.table("/NGS/users/Sid/ENHANCERS/SuperEnh/AllGenes_distEnh/MATRIX_enhCount.txt", header=T, sep="\t", stringsAsFactor=F)
tissues = colnames(mat)[-1]

# file containing tau value for each gene target
tau = read.table("/NGS/users/Sid/ENHANCERS/SuperEnh_27ac/Great_targets/Strict_targets/Tau_analysis/RESULTS.txt", header=F, sep="\t", stringsAsFactor=F)
colnames(tau) = c("Group","tau","Gene","tissue")

# file containing total-expression for each gene target
exp = read.table("/NGS/users/Sid/ENHANCERS/SuperEnh_27ac/Great_targets/Strict_targets/Expression_analysis/RESULTS.txt", header=F, sep="\t", stringsAsFactor=F)
colnames(exp) = c("Group","exp","Gene","tissue")


# changing es-e14 to es.e14 from comparison
tau$tissue = ifelse(tau$tissue == "Es-E14","Es.E14",tau$tissue)
exp$tissue = ifelse(exp$tissue == "Es-E14","Es.E14",exp$tissue)

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
	text = paste("r =",r2.exp,sep="")
	pval = "p < 2.2e-16"

	e.plot = ggplot(res, aes(x=factor(NumbEnh), y=log10(exp+1))) +
	  	geom_boxplot(lwd = 0.35, width = 0.4, outlier.alpha=0.5, outlier.size = 0.12)+
		#stat_summary(fun.y=median, geom="line", aes(group=1), size = 0.1, color = "red", alpha = 1 )+
		#geom_smooth(aes(group=1),size=0.4, method = "lm",linetype="dashed", fullrange=TRUE, span = 10) +
		geom_smooth(aes(group=1),size=0.4, method = "lm", fullrange=T, span = 10, color = "red") +
		#geom_jitter(position=position_jitter(width=.1, height=0),size=0.1,alpha=0.3,fill="white", color = "blue", pch=21)+
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
		#ggtitle(text)+
		#annotate("text",x=47,y=4.4,size=3,label=paste("r^2 ==",,r2.exp,sep=""), parse=T, color = "red") +
		annotate("text",x=47,y=4.4,size=3,label= text, parse=F, color = "red") +
		annotate("text",x=47,y=4.2,size=3,label=pval, parse=T, color = "red") +
		scale_x_discrete(name = "Number of constituent enhancers", breaks = c(1:50), labels = c(1:49,"50+"))+
		scale_y_continuous(name = "Log10 expression")

	file1 = "all.exp.png"
	png(file1,bg="transparent",units="in",width = 6.25, height= 3.75 ,res=600)
	plot(e.plot)
	dev.off()
	
	
	# ### super enh
# 	res = res[res$NumbEnh !=0,]
# 	res$NumbEnh = ifelse(res$NumbEnh > 49, 50, res$NumbEnh)
# 	se = res[res$Group == "Super",]
# 	r2.exp = round(cor(log10(se$exp+1),se$NumbEnh),2)
# 	#text = paste("corr = ",r2.exp,sep="")
# 	text = paste("r^2 =",r2.exp,sep="")
#
# 	e.plot = ggplot(se, aes(x=factor(NumbEnh), y=log10(exp+1))) +
#   	geom_boxplot(lwd = 0.15, width = 0.4, outlier.alpha=0.5, outlier.size = 0.12)+
# 		stat_summary(fun.y=median, geom="line", aes(group=1), size = 0.1, color = "red", alpha = 1 )+
# 		#geom_smooth(aes(group=1),size=0.4, method = "lm",linetype="dashed", fullrange=TRUE, span = 10) +
# 		theme_bw() +
# 		#theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
# 		theme(plot.title = element_text(size=12, hjust=0.5),
# 			legend.background = element_rect(fill = "transparent",colour = NA),
# 			legend.title=element_text(size=12),
# 			legend.text=element_text(size=11),
# 			legend.position="top",
# 			#legend.position= c(0.70,0.98),
# 			legend.key.size= unit(4,"mm"),
# 			axis.text.x = element_text(size=6),
# 			axis.text.y = element_text(size = 6),
# 			panel.border = element_rect(colour="BLACK",size=0.4),
# 			axis.title.x = element_text(size=10),
# 			axis.title.y = element_text(size=10,angle = 90),
# 			panel.background = element_rect(fill="transparent"),
# 			plot.background = element_rect(fill = "transparent",colour = NA)
# 		)+
# 		ggtitle("Super-enhancer associated genes")+
# 		annotate("text",x=47,y=4.2,size=3,label=paste("r^{2} ==",r2.exp,sep=""), parse=T, color = "red") +
# 		scale_x_discrete(name = "Number of constituent enhancers")+
# 		scale_y_continuous(name = "Expression - log(RPKM+1)")
#
# 	file2 = "super.exp.png"
# 	png(file2,bg="transparent",units="in",width = 6.25, height= 3.75 ,res=600)
# 	plot(e.plot)
# 	dev.off()
#
#
# 	### typical enh
# 	res = res[res$NumbEnh !=0,]
# 	res$NumbEnh = ifelse(res$NumbEnh > 49, 50, res$NumbEnh)
# 	te = res[res$Group == "Typical",]
# 	r2.exp = round(cor(log10(te$exp+1),te$NumbEnh),2)
# 	#text = paste("corr = ",r2.exp,sep="")
# 	text = paste("r^2 =",r2.exp,sep="")
#
# 	e.plot = ggplot(te, aes(x=factor(NumbEnh), y=log10(exp+1))) +
#   	geom_boxplot(lwd = 0.15, width = 0.4, outlier.alpha=0.5, outlier.size = 0.12)+
# 		stat_summary(fun.y=median, geom="line", aes(group=1), size = 0.1, color = "red", alpha = 1 )+
# 		#geom_smooth(aes(group=1),size=0.4, method = "lm",linetype="dashed", fullrange=TRUE, span = 10) +
# 		theme_bw() +
# 		#theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
# 		theme(plot.title = element_text(size=12, hjust=0.5),
# 			legend.background = element_rect(fill = "transparent",colour = NA),
# 			legend.title=element_text(size=12),
# 			legend.text=element_text(size=11),
# 			legend.position="top",
# 			#legend.position= c(0.70,0.98),
# 			legend.key.size= unit(4,"mm"),
# 			axis.text.x = element_text(size=6),
# 			axis.text.y = element_text(size = 6),
# 			panel.border = element_rect(colour="BLACK",size=0.4),
# 			axis.title.x = element_text(size=10),
# 			axis.title.y = element_text(size=10,angle = 90),
# 			panel.background = element_rect(fill="transparent"),
# 			plot.background = element_rect(fill = "transparent",colour = NA)
# 		)+
# 		ggtitle("Typical-enhancer associated genes")+
# 		annotate("text",x=47,y=4.2,size=3,label=paste("r^{2} ==",r2.exp,sep=""), parse=T, color = "red") +
# 		scale_x_discrete(name = "Number of constituent enhancers")+
# 		scale_y_continuous(name = "Expression - log(RPKM+1)")
#
# 	file3 = "typical.exp.png"
# 	png(file3,bg="transparent",units="in",width = 6.25, height= 3.75 ,res=600)
# 	plot(e.plot)
# 	dev.off()
	
	
	
	##############################################################
	# for tau
	
	tis.exp = data.frame()
	data.exp = data.frame()
	
res = data.frame(Gene= character(0), tau= numeric(0), NumbEnh = numeric(0), Group = character(0))
for(i in 1:length(tissues)){
	tisName = tissues[i]
	
	tisName = ifelse(tisName == "Es-E14","Es.E14",tisName)
	tis.tau = tau[tau$tissue == tisName & tau$Group != "Absent" & tau$Group != "Weak",]
	data.tau = inner_join(tis.tau,mat)
	data.tau = subset(data.tau, select = c("Gene","tau",tisName,"Group"))
	colnames(data.tau) = c("Gene","tau","NumbEnh","Group")
	res = rbind(res,data.tau)
}

table(res$NumbEnh)
	res = res[res$NumbEnh !=0,]
	res$NumbEnh = ifelse(res$NumbEnh > 49, 50, res$NumbEnh)
	r2.tau = round(cor(res$tau,res$NumbEnh,method = "spearman"),2)
	cor.test(res$tau,res$NumbEnh, method = "spearman", exact=F)
	text = paste("r =",r2.tau,sep="")
	pval = "p < 2.2e-16"

	t.plot = ggplot(res, aes(x=factor(NumbEnh), y=tau)) +
  	geom_boxplot(lwd = 0.35, width = 0.4, outlier.alpha=0.5, outlier.size = 0.12)+
		#stat_summary(fun.y=median, geom="line", aes(group=1), size = 0.1, color = "red", alpha = 1 )+
		#geom_smooth(aes(group=1),size=0.4, method = "lm",linetype="dashed", fullrange=TRUE, span = 10) +
		geom_smooth(aes(group=1),size=0.4, method = "lm", fullrange=T, span = 10, color = "red") +
		#geom_jitter(position=position_jitter(width=.05, height=0),size=0.03,alpha=0.3,fill="white", color = "blue", pch=21)+
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
			axis.text.y = element_text(size =5),
			panel.border = element_rect(colour="BLACK",size=0.4),
			axis.title.x = element_text(size=9),
			axis.title.y = element_text(size=9,angle = 90),
			panel.background = element_rect(fill="transparent"),
			plot.background = element_rect(fill = "transparent",colour = NA)
		)+
		#ggtitle(text2)+
		#annotate("text",x=47,y=1.1,size=3,label=paste("r^2 ==",r2.tau,sep=""), parse=T, color = "red") +
		annotate("text",x=47,y=1.04,size=3,label= text, parse=F, color = "red") +
		annotate("text",x=47,y=0.99,size=3,label=pval, parse=T, color = "red") +
		scale_x_discrete(name = "Number of constituent enhancers", breaks = c(1:50), labels = c(1:49,"50+"))+
		scale_y_continuous(name = expression(paste("Tissue-specific expression"," (",tau["exp-frac"],")")))

		file4 = "all.tau.png"
		png(file4,bg="transparent",units="in",width = 6.25, height= 3.75 ,res=600)
		plot(t.plot)
		dev.off()
		
		
		# ### super enh
# 		res = res[res$NumbEnh !=0,]
# 		res$NumbEnh = ifelse(res$NumbEnh > 49, 50, res$NumbEnh)
# 		se = res[res$Group == "Super",]
# 		r2.tau = round(cor(se$tau,se$NumbEnh),2)
# 		text = paste("r^2 =",r2.tau,sep="")
#
# 		t.plot = ggplot(se, aes(x=factor(NumbEnh), y=tau)) +
# 	  	geom_boxplot(lwd = 0.15, width = 0.4, outlier.alpha=0.5, outlier.size = 0.12)+
# 			stat_summary(fun.y=median, geom="line", aes(group=1), size = 0.1, color = "red", alpha = 1 )+
# 			#geom_smooth(aes(group=1),size=0.4, method = "lm",linetype="dashed", fullrange=TRUE, span = 10) +
# 			theme_bw() +
# 			#theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
# 			theme(plot.title = element_text(size=12, hjust=0.5),
# 				legend.background = element_rect(fill = "transparent",colour = NA),
# 				legend.title=element_text(size=12),
# 				legend.text=element_text(size=11),
# 				legend.position="top",
# 				#legend.position= c(0.70,0.98),
# 				legend.key.size= unit(4,"mm"),
# 				axis.text.x = element_text(size=6),
# 				axis.text.y = element_text(size =6),
# 				panel.border = element_rect(colour="BLACK",size=0.4),
# 				axis.title.x = element_text(size=10),
# 				axis.title.y = element_text(size=10,angle = 90),
# 				panel.background = element_rect(fill="transparent"),
# 				plot.background = element_rect(fill = "transparent",colour = NA)
# 			)+
# 			ggtitle("Super-enhancer associated genes")+
# 			annotate("text",x=47,y=0.97,size=3,label=paste("r^{2} ==",r2.tau,sep=""), parse=T, color = "red") +
# 			scale_x_discrete(name = "Number of constituent enhancers")+
# 			scale_y_continuous(name = expression(paste("Tissue-specific expression"," (",tau["exp-frac"],")")))
#
# 			file5 = "se.tau.png"
# 			png(file5,bg="transparent",units="in",width = 6.25, height= 3.75 ,res=600)
# 			plot(t.plot)
# 			dev.off()
#
#
# 			### typical enh
# 			res = res[res$NumbEnh !=0,]
# 			res$NumbEnh = ifelse(res$NumbEnh > 49, 50, res$NumbEnh)
# 			te = res[res$Group == "Typical",]
# 			r2.tau = round(cor(te$tau,te$NumbEnh),2)
# 			text = paste("r^2 =",r2.tau,sep="")
#
# 			t.plot = ggplot(te, aes(x=factor(NumbEnh), y=tau)) +
# 		  	geom_boxplot(lwd = 0.15, width = 0.4, outlier.alpha=0.5, outlier.size = 0.12)+
# 				stat_summary(fun.y=median, geom="line", aes(group=1), size = 0.1, color = "red", alpha = 1 )+
# 				#geom_smooth(aes(group=1),size=0.4, method = "lm",linetype="dashed", fullrange=TRUE, span = 10) +
# 				theme_bw() +
# 				#theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
# 				theme(plot.title = element_text(size=12, hjust=0.5),
# 					legend.background = element_rect(fill = "transparent",colour = NA),
# 					legend.title=element_text(size=12),
# 					legend.text=element_text(size=11),
# 					legend.position="top",
# 					#legend.position= c(0.70,0.98),
# 					legend.key.size= unit(4,"mm"),
# 					axis.text.x = element_text(size=6),
# 					axis.text.y = element_text(size =6),
# 					panel.border = element_rect(colour="BLACK",size=0.4),
# 					axis.title.x = element_text(size=10),
# 					axis.title.y = element_text(size=10,angle = 90),
# 					panel.background = element_rect(fill="transparent"),
# 					plot.background = element_rect(fill = "transparent",colour = NA)
# 				)+
# 				ggtitle("Typical-enhancer associated genes")+
# 				annotate("text",x=47,y=0.97,size=3,label=paste("r^{2} ==",r2.tau,sep=""), parse=T, color = "red") +
# 				scale_x_discrete(name = "Number of constituent enhancers")+
# 				scale_y_continuous(name = expression(paste("Tissue-specific expression"," (",tau["exp-frac"],")")))
#
# 				file5 = "te.tau.png"
# 				png(file5,bg="transparent",units="in",width = 6.25, height= 3.75 ,res=600)
# 				plot(t.plot)
# 				dev.off()


