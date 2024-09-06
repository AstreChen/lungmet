library(Seurat)
library(tidyverse)
library(ggsci)
library(scales)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(Nebulosa)
#source("~/scripts/000.Mylibrary/0.seurat_pipeline_functions.R")
setwd("/raid1/aichen/projects/scLung/00.neutrophil_described/")

neutro <- read_rds("/raid1/aichen/projects/scLung/data/0.our_KP_PyMT_Neu/06.KP_PyMT_Neutrophil_object.rds")



#### compare tissue signatures ###
deg_sc = readxl::read_xlsx("../../DEG_sets_Cell_Hidalgo/DEG_singlecell.xlsx", skip = 1)

## from figure 3 and figure 4, Cell
blood_sigs <- c("Ifitm1","Ifitm3","Retnlg","Oasl2","Wfdc17","Cxcr2","Oas3")

lung_sigs = deg_sc %>% 
  filter(cluster =="Lung") %>% 
  filter(avg_logFC>1.5) %>%
  arrange(-avg_logFC) %>% pull(gene)

df = data.frame(gene = c(lung_sigs, blood_sigs), tissue=c(rep("Lung", length(lung_sigs)), rep("Blood", length(blood_sigs)) ))
write.csv(df, file = "../00.neutrophil_described/Cell_Andres_lung_and_blood_signatures.csv")

neutro <- AddModuleScore(neutro, features = list(lung_sigs), name="LungSig" )
neutro <- AddModuleScore(neutro, features = list(blood_sigs), name="BloodSig" )



p1 = VlnPlot(neutro, group.by = "stage", features = "LungSig1", pt.size = 0) + geom_boxplot( outlier.shape = NA, width=0.3, fill="white" ) +  scale_fill_manual(values = c( brewer.pal(name = "OrRd", n=3)[c(2,3)], brewer.pal(name = "Blues", n=3)[c(3,1)] ))  
p1
p2 = VlnPlot(neutro, group.by = "stage", features = "BloodSig1", pt.size = 0) + geom_boxplot( outlier.shape = NA, width=0.3, fill="white") + scale_fill_manual(values = c( brewer.pal(name = "OrRd", n=3)[c(2,3)], brewer.pal(name = "Blues", n=3)[c(3,1)] )) 
p2

p = p1 / p2
ggsave(p, filename = "/raid1/aichen/projects/scLung/00.neutrophil_described/figures/005.tissue.boxplot.pdf", width = 4, height = 6)

p0 = plot_density(neutro, "LungSig1", method = "wkde")
ggsave(p0, filename = "/raid1/aichen/projects/scLung/00.neutrophil_described/figures/003.a.UMAP.LungSig1.pdf", width = 4, height = 3)



### Analysis of immunosuppressive neutrophil

ref_immun_sup_df = readxl::read_xlsx("/raid1/aichen/projects/scLung/00.neutrophil_described/Gong_RNAseq_immunosuppress_gene.xlsx", skip = 1, col_names = "gene")
head(ref_immun_sup_df)
ref_immun_sup_df$term = "Gong_etal_Immunosuppressive_Neutrophil"
ref_immun_sup_df = ref_immun_sup_df[,c("term","gene")]


neutro <- AddModuleScore(neutro, features = list(immunosup = ref_immun_sup_df$gene), name="Gong_etal_ImmS")
p0= plot_density(neutro, "Gong_etal_ImmS1", method = "wkde")
ggsave(p0, filename = "/raid1/aichen/projects/scLung/00.neutrophil_described/figures/003.a.UMAP_Gongetal_ImmunSup.pdf", width = 4, height = 3)

mat_neu = subset(neutro, cellanno_L3 =="N5")

### 定性比较: 加上显著性

p1 = VlnPlot(mat_neu, features = "Gong_etal_ImmS1", group.by = "Time") +
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA)

df_summary = p1$data %>% group_by(ident) %>%
  summarise(mean = mean(Gong_etal_ImmS1))

kp_data = p1$data %>%
  filter(str_detect(ident,"KP"))

mt_data = p1$data %>%
  filter(str_detect(ident,"MT"))

p.kp = kp_data %>%
  ggplot(aes(x=ident, y=Gong_etal_ImmS1, fill=ident)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA ) +
  theme_classic() +
  scale_fill_manual(values = c( "Grey80", RColorBrewer::brewer.pal(9,"Blues")[c(5,7,9)]) ) +
  NoLegend() +
  ggsignif::geom_signif( test = "t.test",
    comparisons = list(c("KP_Baseline", "KP_W04") ), map_signif_level = F, y_position = 1.55, tip_length = 0, vjust = 0.2
  ) +  
  ggsignif::geom_signif(  test = "t.test",
    comparisons = list(c("KP_Baseline", "KP_W08") ), map_signif_level = F, y_position = 1.65, tip_length = 0, vjust = 0.2
  ) +
    ggsignif::geom_signif(  test = "t.test",
    comparisons = list(c("KP_Baseline", "KP_W12") ), map_signif_level = F, y_position = 1.75, tip_length = 0, vjust = 0.2
  )+
  geom_hline(yintercept = as.numeric(df_summary[1,2]), linetype="dashed", alpha=0.5) +
  ylim(c(-0.4,1.9)) 
p.kp


p.mt = mt_data %>%
  ggplot(aes(x=ident, y=Gong_etal_ImmS1, fill=ident)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA ) +
  theme_classic() +
  scale_fill_manual(values = c( "Grey80", RColorBrewer::brewer.pal(9,"Blues")[c(5,7,9)]) ) +
  NoLegend()  +
  ggsignif::geom_signif(  test = "t.test",
    comparisons = list(c("MT_Baseline", "MT_W08") ), map_signif_level = F, y_position = 1.55, tip_length = 0, vjust = 0.2
  ) +  
  ggsignif::geom_signif(  test = "t.test",
    comparisons = list(c("MT_Baseline", "MT_W12") ), map_signif_level = F, y_position = 1.65, tip_length = 0, vjust = 0.2
  ) +
  ggsignif::geom_signif( test = "t.test",
    comparisons = list(c("MT_Baseline", "MT_W16") ), map_signif_level = F, y_position = 1.75, tip_length = 0, vjust = 0.2
  )+
  ylim(c(-0.4,1.9)) +
  geom_hline(yintercept =  as.numeric(df_summary[5,2]), linetype="dashed", alpha=0.5)

p = p.kp | p.mt
p
ggsave(p, filename = "/raid1/aichen/projects/scLung/00.neutrophil_described/figures/003.c.KP_MT_signature_score.pdf", width=6, height = 4)
write.csv(df_summary, file = "/raid1/aichen/projects/scLung/00.neutrophil_described/figures/003.c.KP_MT_signature_score.data.csv")

# df_summary$mean[c(6,7,8)] / df_summary$mean[5] 
# [1] 1.414918 1.354643 1.367221
# 
# df_summary$mean[c(2,3,4)] / df_summary$mean[1] 
# [1] 1.079396 1.230929 1.053722

### compare neutrophil proportion ###
sample_df = read.csv("../Sample_Table.csv")

cd45_meta = read.csv("../clustering/final_version/02.KP_PyMT_Immune_MetaData.csv", row.names = 1)
cd45_total = table(cd45_meta$Sample) %>% as.data.frame()
neutro_total = table(mat_neu$Sample) %>% as.data.frame()
colnames(neutro_total) = c("Sample", "n_neutro")
colnames(cd45_total) = c("Sample","n_cd45")

plot_df = merge(cd45_total, neutro_total)
plot_df =plot_df %>%
  mutate(ratio = n_neutro/n_cd45) %>%
  inner_join(., y=sample_df , by="Sample") %>%
  arrange(Type, Time_Stage)
plot_df

df_summary <- plot_df %>%
  group_by(., Type,Time_Stage) %>%
  summarise(mean_ratio = mean(ratio), se_ratio = sd(ratio) / sqrt(n()))


p = df_summary %>%
  ggplot(aes(x=Time_Stage, y=mean_ratio, fill=Time_Stage  )) +
  geom_col() +
  geom_errorbar(data = df_summary, aes(ymin=mean_ratio-se_ratio, ymax=mean_ratio + se_ratio ), width = 0.2 ) +
  geom_line(data = df_summary, mapping = aes(x=Time_Stage,y=mean_ratio, group=Type)) +
  facet_wrap(Type ~ . , scales="free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c( "Grey80", RColorBrewer::brewer.pal(9,"Blues")[c(5,7,9)])) +
  xlab("Mature Neutrophil") +
  ylab("Ratio in CD45+")
p
ggsave(p, filename = "/raid1/aichen/projects/scLung/00.neutrophil_described/figures/003.b.Mature_Neu_change.pdf",width = 4, height = 2.5)

write.csv(df_summary,file = "/raid1/aichen/projects/scLung/00.neutrophil_described/figures/003.b.Mature_Neu_change.data.csv")


#### ligand receptor analysis ########
### ligand-receptor analysis of neutrophils:
```{r}

met_obj = read_rds("/raid1/aichen/projects/scLung/analysis/Lib009/clustering/All_metastatic_cells_obj.rds.gz")
met_obj$cellanno_L1.v2 = factor(met_obj$cellanno_L1.v2, levels = rev(c("B","T","NK","DC",
                                                                       "Neutrophil","Monocyte","Macrophage",
                                                                       "Epithelial","Endothelial","Stromal")))
table(met_obj$cellanno_L2)
table(met_obj$Time)


mouse_celltalkdb = read_rds("../cellTalkDB/mouse_lr_pair.rds")
mouse_celltalkdb %>% filter(receptor_gene_symbol %in% plot_receptors)

########## plot neutrophil specific receptor
source("~/scripts/000.Mylibrary/Annotation_gene_set.R")
### GO 中的receptor
gene_receptors = get_genes_from_GO_mouse("GO:0038023")


met_neutro = subset(neutro, Lineage=="fvb")
rep_neu_exp = FetchData(met_neutro, vars = gene_receptors )
flt_rep = colSums(rep_neu_exp > 0)
flt_rep = flt_rep[flt_rep > (ncol(met_neutro) *0.05)] # not final requirement

flt_rep = flt_rep / ncol(met_neutro)
flt_rep = flt_rep[order(flt_rep, decreasing = T)]

## plot receptor
p = DotPlot(met_obj, features = names(flt_rep), group.by = "cellanno_L2") +
  theme( axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p
#ggsave(p,filename = "./output_figures/neutrophil/add_figures/neutro_expressed_receptor_more0.05.pdf",
#       width =  25 , height = 9) # not final results

tmp = p$data
tmp = tmp %>%
  mutate(other = ifelse(id=="Neutrophil","Neutro","Other")) %>%
  group_by(features.plot,other) %>%
  summarise(pct.exp_other = max(pct.exp), avg.exp_other = max(avg.exp))

set1 = pivot_wider(tmp, id_cols = "features.plot", names_from = "other", values_from = "pct.exp_other" ) %>%
  mutate(ratio = Neutro/Other)

#### final filter #
plot_receptors_df = set1 %>% filter(Neutro>=30 ) %>% filter(Other < 20)
## plot ligand:
#other_receptors = set1 %>% filter(Neutro<=40) %>% filter(str_detect(features.plot,"^Il")) %>% filter(ratio>2) %>% filter(Neutro>25)
#plot_receptors_df = rbind(plot_receptors_df, other_receptors) %>%
#  
plot_receptors = plot_receptors_df %>% filter(features.plot != "Itgb2") %>% filter(features.plot!= "Fcer1g") %>% pull(features.plot) %>% as.character()
plot_receptors

#### plot receptor of neutrophil specific high ###
p = DotPlot(met_obj, features = plot_receptors, group.by = "cellanno_L1.v2") +
  theme( axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_distiller(palette = "Reds", direction = 1)
p
ggsave(p, filename = "fig3_neutro_recruitment/20231211-neutro/Selected_NeutroSpec_Receptors.pdf", width = 7, height = 4.5)

################# plot corresponding ligand ###########
to_plot_lr_df = mouse_celltalkdb %>% filter(receptor_gene_symbol %in% plot_receptors) 
to_plot_lr_df = to_plot_lr_df %>%
  mutate(receptor = factor(receptor_gene_symbol, levels = plot_receptors[plot_receptors %in% to_plot_lr_df$receptor_gene_symbol])) %>%
  arrange(receptor)

plot_receptors[plot_receptors %in% to_plot_lr_df$receptor_gene_symbol]
to_plot_lr_df = to_plot_lr_df %>%
  filter(!(ligand_gene_symbol %in% c("Mif","Ppbp","Selplg")))

p = DotPlot(met_obj, features = unique(to_plot_lr_df$ligand_gene_symbol), group.by = "cellanno_L1.v2") +
  theme( axis.text.x = element_text(angle = 45,hjust = 1, vjust = 1)) +
  scale_color_distiller(palette = "Reds", direction = 1)
p
ggsave(p, filename = "fig3_neutro_recruitment/20231211-neutro/Selected_NeutroSpec_Ligands.pdf", width = 7, height = 4.5)


p = DotPlot(met_obj, features = c("Cxcr2", "Cxcl2","Cxcl1","Cxcl3","Cxcl5"), group.by = "cellanno_L1.v2") +
    theme( axis.text.x = element_text(angle = 45,hjust = 1, vjust = 1)) +
  scale_color_distiller(palette = "Reds", direction = 1)
p
ggsave(p, filename = "/raid1/aichen/projects/scLung/00.neutrophil_described/figures/fig3.Cxcr2_cxcl2.pdf", width = 5, height = 3)



### Cxcl2 expression  #######
pymt_neu = subset(neutro, Lineage=="fvb")
p = VlnPlot(pymt_neu, features = "Cxcl2", group.by = "stage") 
p = p$data %>% filter(ident!="4.other") %>%
  ggplot(aes(x=ident, y=Cxcl2)) +
  geom_violin(aes(fill=ident)) +
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA) +
  ggtitle("Cxcl2") +
  ggpubr::theme_classic2() +
  scale_fill_manual(values = c( brewer.pal(name = "OrRd", n=3)[c(2,3)], brewer.pal(name = "Blues", n=3)[c(3,1)] )) +
  NoLegend() 
ggsave(p, filename = "/raid1/aichen/projects/scLung/00.neutrophil_described/figures/fig3.Cxcl2_stage.pdf", width = 4, height = 3)

p = VlnPlot(pymt_neu, features = "Cxcr2", group.by = "stage") 
p = p$data %>% filter(ident!="4.other") %>%
  ggplot(aes(x=ident, y=Cxcr2)) +
  geom_violin(aes(fill=ident)) +
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA) +
  ggtitle("Cxcr2") +
  ggpubr::theme_classic2() +
  scale_fill_manual(values = c( brewer.pal(name = "OrRd", n=3)[c(2,3)], brewer.pal(name = "Blues", n=3)[c(3,1)] )) +
  NoLegend() 
ggsave(p, filename = "/raid1/aichen/projects/scLung/00.neutrophil_described/figures/fig3.Cxcr2_stage.pdf", width = 4, height = 3)

