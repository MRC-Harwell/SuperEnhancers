#!/bin/Rscript
args<-commandArgs(TRUE)   #1st arument = table , 2nd = output directory path
setwd(args[2])

library("ggplot2")
library("ggsignif")

####### version 2 combining all: gene-tissue pairs  ##############
################# BOXPLOT ########################################

df = read.table(args[1],header=F,sep="\t", stringsAsFactors=F)
colnames(df) = c("Group","RPKM","Gene","Tissue")
df = df[df$Group != "Absent",]
df$Group = ifelse(df$Group == "Super", "SEC",df$Group)
df$Group = ifelse(df$Group == "Typical", "TEC",df$Group)
df$Group = ifelse(df$Group == "Weak", "WEC",df$Group)

df.se = df[df$Group=="SEC",2]
df.te = df[df$Group=="TEC",2]
df.we = df[df$Group=="WEC",2]
man.pval1 = wilcox.test(df.se,df.te)$p.value
man.pval2 = wilcox.test(df.se,df.we)$p.value
label1 = paste("p = ",format.pval(man.pval1,2),sep="")
label2 = paste("p = ",format.pval(man.pval2,2),sep="")

label = "p < 2.2e-16" #custom label in case p-values are extremely low

print(man.pval1)
print(man.pval2)

name = sub(".txt","",args[1])
file_name = paste(name,".summary.png",sep="")
print(file_name)

png(file_name,bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
ggplot(df, aes(x=factor(Group), y=log10(RPKM+1))) +
geom_boxplot(aes(fill = factor(df$Group)), lwd = 0.25, outlier.size = 0, outlier.shape = 1, outlier.alpha = 0.5, alpha=0.7) +
scale_fill_manual(values=c("forestgreen", "orange", "light slate grey")) +
#stat_summary(fun.y=mean, geom="point", shape=5, size=4) +
theme_bw() +
theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
guides(fill=FALSE)+
theme(plot.title = element_text(size=11),
	legend.background = element_rect(fill = "transparent",colour = NA),
	legend.title=element_text(size=9),
	legend.text=element_text(size=9),
	#legend.position="top",
	legend.position= c(0.94,0.85),
	legend.key.size= unit(4,"mm"),
	#axis.text.x = element_text(size=7, angle = 20, hjust = 0.8),
	axis.text.x = element_text(size=11),
	axis.text.y = element_text(size = 11),
	#panel.border = element_rect(colour="BLACK",size=0.4),
	panel.border = element_blank(),
	axis.title.x = element_blank(),
	axis.title.y = element_text(size=12,angle = 90),
	panel.background = element_rect(fill="transparent"),
	plot.background = element_rect(fill = "transparent",colour = NA)
)+
theme(axis.line = element_line(color = "black", size=0.4)) +
scale_x_discrete(name = "\nAssociated genes")+
#scale_y_continuous(name = "log10 expression")+
scale_y_continuous(name = expression(paste(log["10"]," expression")))+
geom_signif(comparisons = list(c("SEC", "WEC")),map_signif_level=TRUE, annotation = label, size=0.4,textsize=3,tip_length=0.02, margin_top = 0.0, color = "black", fontface = 2)+
geom_signif(comparisons = list(c("SEC", "TEC")),map_signif_level=TRUE, annotation = label, size=0.4,textsize=3,tip_length=0.02, margin_top = -0.08, color = "black", fontface = 2)+
geom_signif(comparisons = list(c("TEC", "WEC")),map_signif_level=TRUE, annotation = label, size=0.4,textsize=3,tip_length=0.02, margin_top = -0.16, color = "black", fontface = 2)
dev.off()


# ES
med.diff = median(df.se) - median(df.te)
mad1 = mad(df.se)
mad2 = mad(df.te)
mad.p = sqrt(((mad1)^2 + (mad2)^2)/2)
es1 = round(med.diff/mad.p,2)
#es = ifelse(is.nan(es), 0, es)

med.diff = median(df.se) - median(df.we)
mad1 = mad(df.se)
mad2 = mad(df.we)
mad.p = sqrt(((mad1)^2 + (mad2)^2)/2)
es2 = round(med.diff/mad.p,2)
#es = ifelse(is.nan(es), 0, es)

med.diff = median(df.te) - median(df.we)
mad1 = mad(df.te)
mad2 = mad(df.we)
mad.p = sqrt(((mad1)^2 + (mad2)^2)/2)
es3 = round(med.diff/mad.p,2)

print(es1)
print(es2)
print(es3)
