#!/bin/Rscript
args<-commandArgs(TRUE)   #1st arument = table , 2nd = output directory path

library("ggplot2")
library("grid")
df = read.table(args[1],header=T)
colnames(df) = c("chr","start","end","gene","dist")

breaks = c(-2000,-500,-50,-5,0,5,50,500,2000)
ranges = paste(head(breaks,-1), breaks[-1], sep=" to ")
freq   = hist(df$dist/1000, breaks=breaks, include.lowest=TRUE, plot=FALSE)
data = data.frame(range = ranges, frequency = freq$counts, stringsAsFactors = F)
data$range[1] = "<-500"
data$range[8] = ">500"
data$per = (data$frequency/sum(data$frequency))*100

setwd(args[2])
name = sub(".txt","",args[1])
file_name = paste(name,".png",sep="")
table = paste(name,".resTable.txt",sep="")
print(file_name)
title = "ZT15"

write.table(data,file=table,sep="\t",quote=F,row.names=F,col.names=T)

png(file_name,bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
ggplot(data, aes(x=range,y=per)) +
geom_col(color = "black", fill = "grey", lwd = 0.4,alpha = 0.5, width = 0.7) +
geom_text(label=data$frequency,vjust=-0.3, size= 2.5)+
theme_bw() +
labs(title=title) +
theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
theme(plot.title = element_text(size=14, hjust=0.5, face = "bold"),
        legend.key.size = unit(0.55,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.position= c(0.83,0.9),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=12, angle=45, hjust=1),    #20 before
        axis.text.y = element_text(size=12),
        panel.border = element_rect(colour="BLACK",size=0.4),
        axis.title.x = element_text(size=13, vjust = 0.1),
        axis.title.y = element_text(size=13,angle = 90, vjust = 0.6),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  #scale_y_continuous(name="Region-gene associations (%)",breaks = pretty(data$per, n = 8)) +
  scale_y_continuous(name="Peak-gene associations (%)",breaks = pretty(data$per, n = 7)) +
  scale_x_discrete(name = "Distance to TSS (kb)",
  limits = c("<-500","-500 to -50", "-50 to -5", "-5 to 0", "0 to 5", "5 to 50", "50 to 500", ">500"))
dev.off()
