library(tidyverse)
library(matrixStats)
library(pheatmap)

# SE multiple enhancer tissue-type matrix
se.multiple <- read.table("multiple_cellType_se.genes.txt", header=TRUE, sep="\t") 

# Matrix of number of constituent enhancers per gene-tissue (enhancer usage)
mat = read.table("MATRIX_enhCount_SE.txt", header=TRUE, sep="\t")

mat.se <- dplyr::semi_join(mat, se.multiple, by = "Gene") %>% 
  column_to_rownames("Gene")

write.table(mat.se, "se_multiple_enh_usage.txt", quote=FALSE, row.names=FALSE, sep="\t")

png("eu.se.heat.png",bg="transparent",units="in",width = 7.25, height= 5.75 ,res=600)
p <- pheatmap(mat.se,
         fontsize_col = 9,
         #fontsize_row = 7,
         color = colorRampPalette(c("white","blue","red", "black"))(30),
         cluster_cols = FALSE,
         treeheight_col = 0,
         treeheight_row = 30, 
         show_rownames = FALSE)
plot(p)
dev.off()


mat.se[mat.se ==0] <- NA
mat.se$sd = rowSds(as.matrix(mat.se), na.rm=TRUE)
se.df = mat.se %>% 
  rownames_to_column("gene") %>% 
  dplyr::select(gene, sd) %>% 
  mutate(group = "SEC")




# TE multiple enhancer tissue-type matrix
te.multiple <- read.table("multiple_cellType_te.genes.txt", header=TRUE, sep="\t") 

# Matrix of number of constituent enhancers per gene-tissue (enhancer usage)
mat = read.table("MATRIX_enhCount_TE.txt", header=TRUE, sep="\t")

mat.te <- dplyr::semi_join(mat, te.multiple, by = "Gene") %>% 
  column_to_rownames("Gene")
write.table(mat.te, "te_multiple_enh_usage.txt", quote=FALSE, row.names=FALSE, sep="\t")

png("eu.te.heat.png",bg="transparent",units="in",width = 7.25, height= 5.75 ,res=600)
p <-pheatmap(mat.te,
         fontsize_col = 9,
         #fontsize_row = 7,
         color = colorRampPalette(c("white","blue", "red", "black"))(30),
         cluster_cols = FALSE,
         treeheight_col = 0,
         treeheight_row = 30, 
         show_rownames = FALSE)
plot(p)
dev.off()


mat.te[mat.te ==0] <- NA
mat.te$sd = rowSds(as.matrix(mat.te), na.rm=TRUE)
te.df = mat.te %>% 
  rownames_to_column("gene") %>% 
  dplyr::select(gene, sd) %>% 
  mutate(group = "TEC")



df = bind_rows(se.df, te.df) %>% 
  arrange(desc(sd)) %>% 
  plyr::rename(c("sd" = "enhancer_usage_switch_standard_deviation"))
write.table(df, "Eusage/all_res.txt", quote=FALSE, row.names=FALSE, sep="\t")

top <- df %>% arrange(desc(sd)) %>% dplyr::slice(1:400)
bottom <- df %>% arrange(sd) %>% dplyr::slice(1:400)
write.table(top, "Eusage/top_400.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(bottom, "Eusage/bottom_400.txt", quote=FALSE, row.names=FALSE, sep="\t")


man.pval = format(wilcox.test(se.df$sd,te.df$sd)$p.value,digits=2)
med.diff = median(se.df$sd, na.rm=TRUE) - median(te.df$sd, na.rm=TRUE)
mad1 = mad(se.df$sd, na.rm=TRUE)
mad2 = mad(te.df$sd, na.rm=TRUE)
mad.p = sqrt(((mad1)^2 + (mad2)^2)/2)
es = round(med.diff/mad.p,2)
es = ifelse(is.nan(es), 0, es)
label.es = paste("ES=",es,sep="")
label.p = paste("p=",man.pval,sep="")
label.combine = paste(label.p,"\n",label.es,sep="")


png("Eusage/eus.box.png",bg="transparent",units="in",width = 3.75, height= 4.25 ,res=600)
ggplot(df, aes(x=factor(group), y=enhancer_usage_switch_standard_deviation, color = group)) +
  #geom_boxplot(lwd = 0.4, width = 0.6, outlier.size = 0.5, outlier.alpha=0.8, outlier.color="gray40", outlier.fill=NA, outlier.shape = NA)+
  geom_violin(aes(fill = factor(group)), alpha=0.6, lwd = 0.25, scale = "area") +
  geom_boxplot(width=0.06,lwd = 0.15, outlier.size = 0, outlier.shape = 1, outlier.alpha = 0.5) +
  theme_bw() +
  scale_color_manual(values=c("forestgreen", "orange"),guide=FALSE)+
  scale_fill_manual(values=c("forestgreen", "orange"),guide=FALSE)+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
  theme(plot.title = element_text(size=11, face = "bold", hjust=0.5),
        legend.background = element_rect(fill = "transparent",colour = NA),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        legend.position="top",
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size = 10),
        panel.border = element_rect(colour="BLACK",size=0.4),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  scale_y_continuous(name = "Enhancer usage switch score") +
  geom_signif(comparisons = list(c("SEC", "TEC")),map_signif_level=TRUE, annotation = label.combine, size=0.4,textsize=3,tip_length=0.02, margin_top = -0.08, color = "black", fontface = 2)
dev.off()



