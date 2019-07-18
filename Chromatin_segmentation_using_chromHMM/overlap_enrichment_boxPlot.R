#!/bin/Rscript
args<-commandArgs(TRUE)   #1st arument = table , 2nd = output directory path
setwd(args[2])

library("ggplot2")

df = read.table(args[1],header=F)
colnames(df) = c("Category","Emission_State","Enrichment")

name = sub(".txt","",args[1])
file_name = paste(name,".png",sep="")
#title = paste(args[1],"( p = ", args[3],")",sep ="")
#print(title)
print(file_name)

png(file_name,bg="transparent",units="in",width = 1.7, height= 2 ,res=600)

ggplot(df, aes(x=Emission_State, y=Enrichment,fill = Emission_State)) + geom_boxplot(outlier.size = 0.1,lwd = 0.15) +

scale_fill_manual(values=c("#FF6666", "#FF0000","#00007F","#FFA500","#333333","#FFFF99")) +
guides(fill=FALSE) +
theme_bw() +
labs(title=df$Category) +
theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
coord_flip() +
#geom_jitter(size = 1, color = "cornflowerblue", alpha = 0.8) +

theme(plot.title = element_text(size=11,hjust=0.5),
	legend.background = element_rect(fill = "transparent",colour = NA),
	axis.text.x = element_text(size=6),
	#axis.text.y = element_text(size = 8),
	axis.text.y = element_blank(),
	panel.border = element_rect(colour="BLACK",size=0.4),
	axis.title.x = element_text(size=7),
	#axis.title.y = element_text(size=9,angle = 90),
	axis.title.y = element_blank(),
	axis.ticks.y = element_blank(),
	panel.background = element_rect(fill="transparent"),
	plot.background = element_rect(fill = "transparent",colour = NA)
)+
scale_x_discrete(limits = rev(levels(df$Emission_State)))
dev.off()

