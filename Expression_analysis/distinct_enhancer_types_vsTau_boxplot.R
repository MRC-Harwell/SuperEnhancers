#!/bin/Rscript
args<-commandArgs(TRUE)   #1st arument = table , 2nd = output directory path
setwd(args[2])

library(ggplot2)

name = sub(".txt","",args[1])
file_name = paste(name,".png",sep="")

df = read.table(args[1], header=F, sep="\t")
colnames(df) = c("Gene","Number","Group","Tau")
df$test = ifelse(df$Number>=4,4,df$Number)

lab = paste("\u2265","4",sep="")

png(file_name,bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
ggplot(df, aes(x=factor(test), y=Tau)) +
#stat_boxplot(aes(fill = factor(df$Group)),geom ='errorbar', lwd = 0.3, width = 0.4, position = "dodge") +
#geom_boxplot(lwd = 0.15, outlier.size = NA, coef = 0, width = 0.5)+
#geom_boxplot(aes(fill = factor(df$Group)),lwd = 0.15, width = 0.4, outlier.alpha=0.1, outlier.size = 0.2, outlier.color="#CCE5FF")+
geom_boxplot(aes(fill = factor(df$Group)), lwd = 0.25, outlier.size = 0, outlier.shape = 1, outlier.alpha = 0.1, alpha=0.7, width = 0.5) +
geom_smooth(aes(group=df$Group, color = df$Group),size=0.4, method = "lm",linetype="dashed", fullrange=TRUE, span = 10) +
scale_color_manual(values=c("forestgreen", "darkorange3"), guide=F)+
theme_bw() +
scale_fill_manual(values=c("forestgreen", "orange"),
						#name="Genes",
                         breaks=c("Super", "Typical"),
                         labels=c("SEC", "TEC")) +
theme_bw() +
theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
guides(fill=guide_legend(ncol=2,title.position = "left"))+
theme(plot.title = element_text(size=10, hjust=0.5),
	legend.background = element_rect(fill = "transparent",colour = NA),
	#legend.title=element_text(size=12),
	legend.title=element_blank(),
	legend.text=element_text(size=11),
	legend.position="top",
	legend.key.width=unit(2,"line"),
	#legend.position= c(0.70,0.98),
	#legend.key.size= unit(4,"mm"),
	axis.text.x = element_text(size=11),
	axis.text.y = element_text(size = 11),
	#panel.border = element_rect(colour="BLACK",size=0.4),
	panel.border = element_blank(),
	axis.title.x = element_text(size=12),
	axis.title.y = element_text(size=12,angle = 90),
	panel.background = element_rect(fill="transparent"),
	plot.background = element_rect(fill = "transparent",colour = NA)
)+
theme(axis.line = element_line(color = "black", size=0.4)) +
scale_x_discrete(name = "# of distinct enhancer tissue-type", breaks = c(1,2,3,4), labels = c(1,2,3,lab))+
scale_y_continuous(name = expression(paste("Tissue-specific expression"," (",tau["exp-frac"],")")))
dev.off()


super = df[df$Group == "Super",]
lm.super = lm(test~Tau,data=super)
print(summary(lm.super))

typ = df[df$Group == "Typical",]
lm.typ = lm(test~Tau,data=typ)
print(summary(lm.typ))
