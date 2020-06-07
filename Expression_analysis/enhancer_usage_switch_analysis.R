library(tidyverse)
library(matrixStats)
library(ggpubr)
library(ggrepel)


se = read.table("/tmp_mnt/filer1/bioinformatics/Sid/SE/Eusage/Files/cellTypeSpec_se_MAT.txt", header=TRUE, sep="\t")
se$sum = rowSums(se[,-1])
a = se %>% filter(sum == 1)
write.table(a, "/tmp_mnt/filer1/bioinformatics/Sid/SE/Eusage/Files/one_cellType_se.genes.txt", quote=FALSE, row.names=FALSE, sep="\t")
a = se %>% filter(sum > 1)
write.table(a, "/tmp_mnt/filer1/bioinformatics/Sid/SE/Eusage/Files/multiple_cellType_se.genes.txt", quote=FALSE, row.names=FALSE, sep="\t")

te = read.table("/tmp_mnt/filer1/bioinformatics/Sid/SE/Eusage/Files/cellTypeSpec_te_MAT.txt", header=TRUE, sep="\t")
te$sum = rowSums(te[,-1])
b = te %>% filter(sum == 1)
write.table(b, "/tmp_mnt/filer1/bioinformatics/Sid/SE/Eusage/Files/one_cellType_te.genes.txt", quote=FALSE, row.names=FALSE, sep="\t")
b = te %>% filter(sum > 1)
write.table(b, "/tmp_mnt/filer1/bioinformatics/Sid/SE/Eusage/Files/multiple_cellType_te.genes.txt", quote=FALSE, row.names=FALSE, sep="\t")




se.multiple <- read.table("/tmp_mnt/filer1/bioinformatics/Sid/SE/Eusage/Files/multiple_cellType_se.genes.txt", header=TRUE, sep="\t") 

mat = read.table("/tmp_mnt/filer1/bioinformatics/Sid/SE/Eusage/Files/MATRIX_enhCount_SE.txt", header=TRUE, sep="\t")

mat.se <- dplyr::semi_join(mat, se.multiple, by = "Gene") %>% 
  column_to_rownames("Gene")
#write.table(mat.se, "/tmp_mnt/filer1/bioinformatics/Sid/SE/Eusage/se_multiple_enh_usage.txt", quote=FALSE, row.names=FALSE, sep="\t")

png("/tmp_mnt/filer1/bioinformatics/Sid/SE/Eusage/eu.se.heat.png",bg="transparent",units="in",width = 7.25, height= 5.75 ,res=600)
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
#rowVars(as.matrix(mat.se), na.rm=TRUE)
mat.se$sd = rowSds(as.matrix(mat.se), na.rm=TRUE)
se.df = mat.se %>% 
  rownames_to_column("gene") %>% 
  dplyr::select(gene, sd) %>% 
  mutate(group = "SEC")

se.3qt <- se.df %>% filter(sd >= summary(se.df$sd)[5]) %>% arrange(desc(sd))
se.1qt <- se.df %>% filter(sd <= summary(se.df$sd)[2]) %>% arrange(sd)
write.table(se.3qt, "/tmp_mnt/filer1/bioinformatics/Sid/SE/Eusage/se_3qt.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(se.1qt, "/tmp_mnt/filer1/bioinformatics/Sid/SE/Eusage/se_1qt.txt", quote=FALSE, row.names=FALSE, sep="\t")


te.multiple <- read.table("/tmp_mnt/filer1/bioinformatics/Sid/SE/Eusage/Files/multiple_cellType_te.genes.txt", header=TRUE, sep="\t") 

mat = read.table("/tmp_mnt/filer1/bioinformatics/Sid/SE/Eusage/Files/MATRIX_enhCount_TE.txt", header=TRUE, sep="\t")

mat.te <- dplyr::semi_join(mat, te.multiple, by = "Gene") %>% 
  column_to_rownames("Gene")
#write.table(mat.te, "/tmp_mnt/filer1/bioinformatics/Sid/SE/Eusage/te_multiple_enh_usage.txt", quote=FALSE, row.names=FALSE, sep="\t")

png("/tmp_mnt/filer1/bioinformatics/Sid/SE/Eusage/eu.te.heat.png",bg="transparent",units="in",width = 7.25, height= 5.75 ,res=600)
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
#rowVars(as.matrix(mat.se), na.rm=TRUE)
mat.te$sd = rowSds(as.matrix(mat.te), na.rm=TRUE)
te.df = mat.te %>% 
  rownames_to_column("gene") %>% 
  dplyr::select(gene, sd) %>% 
  mutate(group = "TEC")

te.3qt <- te.df %>% filter(sd >= summary(te.df$sd)[5]) %>% arrange(desc(sd))
te.1qt <- te.df %>% filter(sd <= summary(te.df$sd)[2]) %>% arrange(sd)
write.table(te.3qt, "/tmp_mnt/filer1/bioinformatics/Sid/SE/Eusage/te_3qt.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(te.1qt, "/tmp_mnt/filer1/bioinformatics/Sid/SE/Eusage/te_1qt.txt", quote=FALSE, row.names=FALSE, sep="\t")


df = bind_rows(se.df, te.df) %>% 
  arrange(desc(sd)) %>% 
  plyr::rename(c("sd" = "enhancer_usage_switch_standard_deviation"))
write.table(df, "/tmp_mnt/filer1/bioinformatics/Sid/SE/Eusage/all_res.txt", quote=FALSE, row.names=FALSE, sep="\t")

top <- df %>% arrange(desc(sd)) %>% dplyr::slice(1:400)
bottom <- df %>% arrange(sd) %>% dplyr::slice(1:400)
write.table(top, "/tmp_mnt/filer1/bioinformatics/Sid/SE/Eusage/top_400.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(bottom, "/tmp_mnt/filer1/bioinformatics/Sid/SE/Eusage/bottom_400.txt", quote=FALSE, row.names=FALSE, sep="\t")





ggdensity(df, x = "sd2", y = "..density..",
          add = "median", rug = TRUE,
          color = "group", fill = "group",
          palette = c("#00AFBB", "#404040")
)+
  scale_x_continuous(name = "Number of RBP motif binding per 100 nt", trans="log2")+
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.40,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.position= c(0.85,0.86),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size = 9),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        #panel.border = element_blank(),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )



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




png("/tmp_mnt/filer1/bioinformatics/Sid/SE/Eusage/eus.box.png",bg="transparent",units="in",width = 3.75, height= 4.25 ,res=600)
ggplot(df, aes(x=factor(group), y=enhancer_usage_switch_standard_deviation, color = group)) +
  #geom_boxplot(lwd = 0.4, width = 0.6, outlier.size = 0.5, outlier.alpha=0.8, outlier.color="gray40", outlier.fill=NA, outlier.shape = NA)+
  geom_violin(aes(fill = factor(group)), alpha=0.6, lwd = 0.25, scale = "area") +
  geom_boxplot(width=0.06,lwd = 0.15, outlier.size = 0, outlier.shape = 1, outlier.alpha = 0.5) +
  #geom_jitter(aes(color=group),shape=16, position=position_jitter(0.3), size=0.4, alpha=0.3)+
  theme_bw() +
  # geom_text_repel(
  #   size = 2,
  #   segment.color = "black",
  #   segment.size = 0.1,
  #   direction = "y", 
  #   #hjust = 1, 
  #   nudge_x = 0.35
  # ) +
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



