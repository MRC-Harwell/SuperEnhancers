library(tidyverse)
library(data.table)
library(ggsignif)

#tissues <- (
  #"BAT",
  #"BmarrowDm",
  #	"Bmarrow",
  #"Cerebellum",
  #	"CH12",
  #"Cortex",
  # 	"Esb4",
  #  "Es-E14",
  #"Heart",
  #"Kidney"
  # 	"Limb",
  # 	"Liver",
  # 	"Lung",
  # 	"MEF",
  #	"MEL",
  # 	"OlfactoryBulb",
  #  "Placenta",
  # 	"SmallIntestine",
  #	"Spleen",
  # 	"Testis",
  #	"Thymus",
  #"Wbrain"
#)


analyse <- function(x){
  x.h = dplyr::left_join(x,orths, by = c("V1" = "Mouse_gene_name")) %>% 
    dplyr::select(Human_ID, Human_gene_name) %>% 
    drop_na()
  res <- dplyr::left_join(x.h, eqtl.group, by =c("Human_ID" = "gene_id")) %>% 
    mutate(count = replace_na(count, 0), tissue = tissue)
  return(res)
}

orths <- read.table("/tmp_mnt/filer1/bioinformatics/Sid/Annotation/mouse_human_orthologs.txt", header=TRUE, sep="\t", na.strings = "")

tissue2 = "Heart_Left_Ventricle"
eqtl <- read.table(str_c("/tmp_mnt/filer1/bioinformatics/Sid/Gtex_v8/eQTL/GTEx_Analysis_v8_eQTL_independent/", tissue2, ".v8.independent_eqtls.txt"), header=TRUE, sep="\t")

#eqtl <- fread(str_c("/tmp_mnt/filer1/bioinformatics/Sid/Gtex_v8/eQTL/TissueSpecificAllSnpGene/GTEx_Analysis_v8_eQTL/", tissue2, ".v8.signif_variant_gene_pairs.txt"), header=TRUE, sep="\t")
eqtl$gene_id = str_replace(eqtl$gene_id, "\\.\\d", "")

eqtl.group <- eqtl %>% 
  group_by(gene_id) %>% 
  summarise(count = n() 
            #sum_slope = sum(slope), 
            #avg_slope = mean(slope)
            )


tissue = "Heart"
se = read.table(str_c("/tmp_mnt/filer1/bioinformatics/Sid/SE/", tissue, ".se.genes.txt"), header=FALSE)
te = read.table(str_c("/tmp_mnt/filer1/bioinformatics/Sid/SE/", tissue, ".te.uniq.genes.txt"), header=FALSE)
  
res.se <- analyse(se) %>% mutate(group = "SEC")
res.te <- analyse(te) %>% mutate(group = "TEC")

write.table(res.se, str_c(tissue, ".se.txt"), quote=FALSE, row.names=FALSE, sep="\t")
write.table(res.te, str_c(tissue, ".te.txt"), quote=FALSE, row.names=FALSE, sep="\t")







# ################ plot #############
files <- list.files("/tmp_mnt/filer1/bioinformatics/Sid/SE/Independent_v8", full.names = TRUE)

df <- map_dfr(files, read.table, header=TRUE)


df.se <- df %>% filter(group %in% "SEC")
df.te <- df %>% filter(group %in% "TEC")
man.pval = format(wilcox.test(df.se$count,df.te$count)$p.value,digits=2)
med.diff = median(df.se$count) - median(df.te$count)
mad1 = mad(df.se$count)
mad2 = mad(df.te$count)
mad.p = sqrt(((mad1)^2 + (mad2)^2)/2)
es = round(med.diff/mad.p,2)
es = ifelse(is.nan(es), 0, es)
label.es = paste("ES=",es,sep="")
label.p = paste("p=",man.pval,sep="")
label.combine = paste(label.p,"\n",label.es,sep="")



# png("/tmp_mnt/filer1/bioinformatics/Sid/SE/Variant_gene_pair_v8/eqtl_count.box.png",bg="transparent",units="in",width = 2.20, height= 3.75 ,res=600)
# ggplot(df, aes(x=factor(group), y=count+1, color = group)) +
# geom_boxplot(lwd = 0.4, width = 0.6, outlier.size = 0.5, outlier.alpha=0.8, outlier.color="gray40", outlier.fill=NA, outlier.shape = NA)+
# geom_jitter(aes(color=group),shape=16, position=position_jitter(0.3), size=0.4, alpha=0.3)+
# theme_bw() +
# scale_color_manual(values=c("forestgreen", "orange"),guide=FALSE)+	theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
# theme(plot.title = element_text(size=11, face = "bold", hjust=0.5),
# 		legend.background = element_rect(fill = "transparent",colour = NA),
# 		legend.title=element_text(size=9),
# 		legend.text=element_text(size=8),
# 		legend.position="top",
# 		axis.text.x = element_text(size=12),
# 		axis.text.y = element_text(size = 10),
# 		panel.border = element_rect(colour="BLACK",size=0.4),
# 		axis.title.x = element_blank(),
# 		axis.title.y = element_text(size=12,angle = 90),
# 		panel.background = element_rect(fill="transparent"),
# 		plot.background = element_rect(fill = "transparent",colour = NA)
# 	)+
# scale_y_continuous(name = "# of tissue-specific eQTL associations", trans="log2") +
# geom_signif(comparisons = list(c("SEC", "TEC")),map_signif_level=TRUE, annotation = label.combine, size=0.4,textsize=3,tip_length=0.02, margin_top = -0.08, color = "black", fontface = 2)
# dev.off()


# for independent
png("/tmp_mnt/filer1/bioinformatics/Sid/SE/Independent_v8/eqtl_count.box.png",bg="transparent",units="in",width = 2.20, height= 3.75 ,res=600)
ggplot(df, aes(x=factor(group), y=count, color = group)) +
  geom_boxplot(lwd = 0.4, width = 0.6, outlier.size = 0.5, outlier.alpha=0.8, outlier.color="gray40", outlier.fill=NA, outlier.shape = NA)+
  geom_jitter(aes(color=group),shape=16, position=position_jitter(0.3), size=0.4, alpha=0.3)+
  theme_bw() +
  scale_color_manual(values=c("forestgreen", "orange"),guide=FALSE)+	theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
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
  scale_y_continuous(name = "# conditionally independent eQTL associations") +
  geom_signif(comparisons = list(c("SEC", "TEC")),map_signif_level=TRUE, annotation = label.combine, size=0.4,textsize=3,tip_length=0.02, margin_top = -0.08, color = "black", fontface = 2)
dev.off()


