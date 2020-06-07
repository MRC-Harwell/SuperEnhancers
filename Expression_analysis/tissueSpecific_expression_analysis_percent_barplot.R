#!/bin/Rscript
args<-commandArgs(TRUE)   #1st arument = table , 2nd = output directory path
setwd(args[2])

library("ggplot2")

df = read.table(args[1],header=F)
colnames(df) = c("Group","Tau","Gene","Tissue")

name = sub(".txt","",args[1])
file_1 = paste(name,".high.barplot.png",sep="")
print(file_1)


############# for highly tissue specific ############
res = data.frame() # to store fraction
res_n = data.frame() # to store exact numbers

tissues = unique(df$Tissue)
j=1
for(i in 1:length(tissues)){
	tis = tissues[i]
	data = df[df$Tissue == tissues[i],]

	super = data[data$Group == "Super",]
	tau_s = super[super$Tau >= 0.85,]
	per_s = (nrow(tau_s)/nrow(super))*100
	res[j,"Group"] = "Super"
	res[j,"Fraction"] = per_s
	res[j,"Number"] = nrow(tau_s)
	res[j,"Tissue"] = tis

	typical = data[data$Group == "Typical",]
	tau_t = typical[typical$Tau >= 0.85,]
	per_t = (nrow(tau_t)/nrow(typical))*100
	res[j+1,"Group"] = "Typical"
	res[j+1,"Fraction"] = per_t
	res[j+1,"Number"] = nrow(tau_t)
	res[j+1,"Tissue"] = tis

	weak = data[data$Group == "Weak",]
	tau_w = weak[weak$Tau >= 0.85,]
	per_w = (nrow(tau_w)/nrow(weak))*100
	res[j+2,"Group"] = "Weak"
	res[j+2,"Fraction"] = per_w
	res[j+2,"Number"] = nrow(tau_w)
	res[j+2,"Tissue"] = tis

	j = j + 3
}

cat("\nMean % of genes with tau_fraction >= 0.85\n")
cat("Super: ")
super = mean(res[res$Group == "Super",2])
print(super)
cat("\n")
cat("Typical: ")
typical = mean(res[res$Group == "Typical",2])
print(typical)
cat("\n")
cat("Weak: ")
weak = mean(res[res$Group == "Weak",2])
print(weak)
cat("\n")


# bar plot #########
png(file_1,bg="transparent",units="in",width = 10.25, height= 4.25 ,res=600)
ggplot(res, aes(x=factor(Tissue), y=Fraction)) +
geom_bar(aes(fill = factor(res$Group)), stat = "identity", position = "dodge", alpha=0.9) +
scale_fill_manual(values=c("forestgreen", "orange", "light slate grey"),
						name="",
                         breaks=c("Super", "Typical", "Weak"),
                         labels=c("SEC", "TEC", "WEC")) +
theme_bw() +
geom_hline(aes(yintercept=super), colour="forestgreen", linetype="dashed", lwd = 0.3) +
#annotate("text",x=21.80,y=super+2,size=2,label=c("mean super-enh")) +
geom_hline(aes(yintercept=typical), colour="orange", linetype="dashed", lwd = 0.3) +
#annotate("text",x=21.80,y=typical+2,size=2,label=c("mean typical-enh")) +
geom_hline(aes(yintercept=weak), colour="light slate grey", linetype="dashed", lwd = 0.3) +
#annotate("text",x=21.80,y=weak-1,size=2,label=c("mean weak-enh")) +
ggtitle(expression(paste("Tissue specifcity index ","(",tau["exp-frac"],")"," \u2265 0.85"),sep="")) +
#theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +

theme(plot.title = element_text(size=12, hjust = 0.5, face= "bold", margin = margin(t = 10, b = -17)),
	legend.background = element_rect(fill = "transparent",colour = NA),
	legend.title=element_text(size=11),
	legend.text=element_text(size=11),
	#legend.position="top",
	legend.key.size= unit(4,"mm"),
	legend.position= c(0.96,0.89),
	axis.text.x = element_text(size=11, angle = 90, vjust = 0.3), # 0.8
	axis.text.y = element_text(size = 11),
	panel.border = element_rect(colour="BLACK",size=0.4),
	axis.title.x = element_text(size=12),
	axis.title.y = element_text(size=11,angle = 90),
	panel.background = element_rect(fill="transparent"),
	plot.background = element_rect(fill = "transparent",colour = NA)
)+
scale_x_discrete(name = "Tissues")+
scale_y_continuous(name = "Fraction of target genes with\nhighly tissue-specific expression (%)", limits = c(0,100))
dev.off()

#############################################################
##############################################################
####### For intermediate tissue specificity ################

res2 = data.frame()
tissues = unique(df$Tissue)
j=1
for(i in 1:length(tissues)){
	tis = tissues[i]
	data = df[df$Tissue == tissues[i],]

	super = data[data$Group == "Super",]
	tau_s = super[super$Tau < 0.85 & super$Tau > 0.20,]
	per_s = (nrow(tau_s)/nrow(super))*100
	res2[j,"Group"] = "Super"
	res2[j,"Fraction"] = per_s
	res2[j,"Tissue"] = tis

	typical = data[data$Group == "Typical",]
	tau_t = typical[typical$Tau < 0.85 & typical$Tau > 0.20,]
	per_t = (nrow(tau_t)/nrow(typical))*100
	res2[j+1,"Group"] = "Typical"
	res2[j+1,"Fraction"] = per_t
	res2[j+1,"Tissue"] = tis

	weak = data[data$Group == "Weak",]
	tau_w = weak[weak$Tau < 0.85 & weak$Tau > 0.20,]
	per_w = (nrow(tau_w)/nrow(weak))*100
	res2[j+2,"Group"] = "Weak"
	res2[j+2,"Fraction"] = per_w
	res2[j+2,"Tissue"] = tis

	j = j + 3
}

cat("\nMean % of genes with tau_fraction < 0.85 and >0.20\n")
cat("Super: ")
super = mean(res2[res2$Group == "Super",2])
print(super)
cat("\n")
cat("Typical: ")
typical = mean(res2[res2$Group == "Typical",2])
print(typical)
cat("\n")
cat("Weak: ")
weak = mean(res2[res2$Group == "Weak",2])
print(weak)
cat("\n")

# ### bar plot #########
file_3 = paste(name,".intermediate.barplot.png",sep="")
png(file_3,bg="transparent",units="in",width = 10.25, height= 4.25 ,res=600)
ggplot(res2, aes(x=factor(Tissue), y=Fraction)) +
geom_bar(aes(fill = factor(res2$Group)), stat = "identity", position = "dodge", alpha=0.9) +
scale_fill_manual(values=c("forestgreen", "orange", "light slate grey"),
						name="Genes",
                         breaks=c("Super", "Typical", "Weak"),
                         labels=c("with super enhancers", "with typical enhancers", "with weak enhancers")) +
theme_bw() +
geom_hline(aes(yintercept=super), colour="forestgreen", linetype="dashed", lwd = 0.3) +
#annotate("text",x=21.80,y=super+2,size=2,label=c("mean super-enh")) +
geom_hline(aes(yintercept=typical), colour="orange", linetype="dashed", lwd = 0.3) +
#annotate("text",x=21.80,y=typical+2,size=2,label=c("mean typical-enh")) +
geom_hline(aes(yintercept=weak), colour="light slate grey", linetype="dashed", lwd = 0.3) +
#annotate("text",x=21.80,y=weak+2,size=2,label=c("mean weak-enh")) +
ggtitle(expression(paste("0.20 < ","Tissue specifcity index ","(",tau["exp-frac"],")"," < 0.85"),sep="")) +
#theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
guides(fill=FALSE) +
theme(plot.title = element_text(size=12, hjust = 0.5, face= "bold", margin = margin(t = 10, b = -17)),
	legend.background = element_rect(fill = "transparent",colour = NA),
	legend.title=element_text(size=11),
	legend.text=element_text(size=11),
	#legend.position="top",
	legend.key.size= unit(4,"mm"),
	legend.position= c(0.88,0.85),
	axis.text.x = element_text(size=11, angle = 90, vjust = 0.3), # 0.8
	axis.text.y = element_text(size = 11),
	panel.border = element_rect(colour="BLACK",size=0.4),
	axis.title.x = element_text(size=12),
	axis.title.y = element_text(size=11,angle = 90),
	panel.background = element_rect(fill="transparent"),
	plot.background = element_rect(fill = "transparent",colour = NA)
)+
scale_x_discrete(name = "Tissues")+
scale_y_continuous(name = "Fraction of target genes with\nintermediate tissue-specific expression (%)",limits = c(0,100))
dev.off()

#############################################################
##############################################################
####### For LOW tissue specificity ################

res3 = data.frame()
tissues = unique(df$Tissue)
j=1
for(i in 1:length(tissues)){
	tis = tissues[i]
	data = df[df$Tissue == tissues[i],]

	super = data[data$Group == "Super",]
	tau_s = super[super$Tau <= 0.20,]
	per_s = (nrow(tau_s)/nrow(super))*100
	res3[j,"Group"] = "Super"
	res3[j,"Fraction"] = per_s
	res3[j,"Tissue"] = tis

	typical = data[data$Group == "Typical",]
	tau_t = typical[typical$Tau <= 0.20,]
	per_t = (nrow(tau_t)/nrow(typical))*100
	res3[j+1,"Group"] = "Typical"
	res3[j+1,"Fraction"] = per_t
	res3[j+1,"Tissue"] = tis

	weak = data[data$Group == "Weak",]
	tau_w = weak[weak$Tau <= 0.20,]
	per_w = (nrow(tau_w)/nrow(weak))*100
	res3[j+2,"Group"] = "Weak"
	res3[j+2,"Fraction"] = per_w
	res3[j+2,"Tissue"] = tis

	j = j + 3
}

cat("\nMean % of genes with tau_fraction <= 0.20\n")
cat("Super: ")
super = mean(res3[res3$Group == "Super",2])
print(super)
cat("\n")
cat("Typical: ")
typical = mean(res3[res3$Group == "Typical",2])
print(typical)
cat("\n")
cat("Weak: ")
weak = mean(res3[res3$Group == "Weak",2])
print(weak)
cat("\n")

### bar plot #########
file_4 = paste(name,".low.barplot.png",sep="")
png(file_4,bg="transparent",units="in",width = 10.25, height= 4.25 ,res=600)
ggplot(res3, aes(x=factor(Tissue), y=Fraction)) +
geom_bar(aes(fill = factor(res3$Group)), stat = "identity", position = "dodge", alpha=0.9) +
scale_fill_manual(values=c("forestgreen", "orange", "light slate grey"),
						name="Genes",
                         breaks=c("Super", "Typical", "Weak"),
                         labels=c("with super enhancers", "with typical enhancers", "with weak enhancers")) +
theme_bw() +
geom_hline(aes(yintercept=super), colour="forestgreen", linetype="dashed", lwd = 0.3) +
#annotate("text",x=21.80,y=super+2,size=2,label=c("mean super-enh")) +
geom_hline(aes(yintercept=typical), colour="orange", linetype="dashed", lwd = 0.3) +
#annotate("text",x=21.80,y=typical+2,size=2,label=c("mean typical-enh")) +
geom_hline(aes(yintercept=weak), colour="light slate grey", linetype="dashed", lwd = 0.3) +
#annotate("text",x=21.80,y=weak+2,size=2,label=c("mean weak-enh")) +
ggtitle(expression(paste("Tissue specifcity index ","(",tau["exp-frac"],")"," \u2264 0.20"),sep="")) +
#theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
guides(fill=FALSE) +
theme(plot.title = element_text(size=12, hjust = 0.5, face="bold", margin = margin(t = 10, b = -17)),
	legend.background = element_rect(fill = "transparent",colour = NA),
	legend.title=element_text(size=11),
	legend.text=element_text(size=11),
	#legend.position="top",
	legend.key.size= unit(4,"mm"),
	legend.position= c(0.88,0.85),
	axis.text.x = element_text(size=11, angle = 90, vjust = 0.3, hjust=1), # 0.8
	axis.text.y = element_text(size = 11),
	panel.border = element_rect(colour="BLACK",size=0.4),
	axis.title.x = element_text(size=12),
	axis.title.y = element_text(size=11,angle = 90),
	panel.background = element_rect(fill="transparent"),
	plot.background = element_rect(fill = "transparent",colour = NA)
)+
scale_x_discrete(name = "Tissue")+
scale_y_continuous(name = "Fraction of target genes with\nlow tissue-specific expression (%)",limits = c(0,100))
dev.off()


