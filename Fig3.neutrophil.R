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


### UMAP ###
metadf = neutro@meta.data
tb_list = split.data.frame(metadf, f = metadf$Tumor)

names(tb_list) = c("normal", "primary","metastatic")

p0 = ggplot(data=tb_list$normal , aes(x=UMAP1,y=UMAP2))+
  geom_point(data=metadf,aes(x=UMAP1,y=UMAP2), color="grey80", size=0.2) +
  geom_point(aes(color=cellanno_L3), size=0.5) +
  scale_color_manual(
    values = c(
      "#8C564BFF", "#D62728FF", "#2CA02CFF", "#FF7F0EFF" ,"#1F77B4FF" ,"#9467BDFF", "#E377C2FF", 
      "#17BECFFF", "#F57676", "#AEC7E8FF", "#BCBD22FF" 
    )
  ) +
  theme_classic() +
  theme(element_blank(), legend.position = "none") +
  geom_density_2d()


p1=   ggplot(data=tb_list$primary , aes(x=UMAP1,y=UMAP2))+
  geom_point(data=metadf,aes(x=UMAP1,y=UMAP2), color="grey80", size=0.2) +
  geom_point(aes(color=cellanno_L3), size=0.5) +
  scale_color_manual(
    values = c(
      "#8C564BFF", "#D62728FF", "#2CA02CFF", "#FF7F0EFF" ,"#1F77B4FF" ,"#9467BDFF", "#E377C2FF", 
      "#17BECFFF", "#F57676",  "#BCBD22FF" 
    )
  ) +
  theme_classic() +
  theme(element_blank(), legend.position = "none") +
  geom_density_2d()



p2 = ggplot(data=tb_list$metastatic , aes(x=UMAP1,y=UMAP2))+
  geom_point(data=metadf,aes(x=UMAP1,y=UMAP2), color="grey80", size=0.2) +
  geom_point(aes(color=cellanno_L3), size=0.5) +
  scale_color_manual(
    values = c(
      "#8C564BFF", "#D62728FF", "#2CA02CFF", "#FF7F0EFF" ,"#1F77B4FF" ,"#9467BDFF", "#E377C2FF", 
      "#17BECFFF", "#F57676", "#AEC7E8FF", "#BCBD22FF" 
    )
  ) +
  theme_classic() +
  theme(element_blank(), legend.position = "none") +
  geom_density_2d()


p = p0|p1|p2
ggsave(p, filename = "figure/UMAP_compare_3_conditions.pdf", width=7,height = 2.5  )
p



######## plot cell proportion of each condition #########
count_df = table(neutro@meta.data$Tumor, neutro@meta.data$cellanno_L3)
total_num = table(neutro$Tumor)
ratio_df = count_df/ as.vector(total_num) 

ratio_df = data.frame(ratio_df)
ratio_df$Var2= factor(ratio_df$Var2, levels = rev(paste0("N", 1:11)))

p = ratio_df %>%
  ggplot(aes(x=Freq, y=Var1, fill=Var2))+
  geom_bar(width=0.7, stat="identity",position = "stack") +
    scale_fill_manual(
    values = rev(c(
      "#8C564BFF", "#D62728FF", "#2CA02CFF", "#FF7F0EFF" ,"#1F77B4FF" ,"#9467BDFF", "#E377C2FF", 
      "#17BECFFF", "#F57676", "#AEC7E8FF", "#BCBD22FF" 
    ))
  )+
  theme_classic()

write.csv(ratio_df, file = "cell_prop_of_neutrophil.csv")
ggsave(p, filename = "figure/cell_prop_neutrophil_by_condition.pdf", width=6, height = 2.5)

p

## cluster name:
cls_names = c("N1.Ngp", "N2.Mmp8","N3.Retnlg","N4.Steap4",
          "N5.Ptgs2","N6.Isg15","N7.Il1r2","N8.Cd74","N9.Ccl3Cd274","N10.Trps1","N11.mt-Co1")
### get sample table:
sample_tb = data.frame(neutro@meta.data[,c("Sample","Lineage","Type","Tumor","Time")]) %>% unique() %>% as.tibble()
total_num = neutro$Sample %>% table()  %>% data.frame()

mat_cls = table( neutro$Sample, neutro$cellanno_L3) %>% as.matrix()
cell_prop = mat_cls / rowSums(mat_cls)

cell_prop_df = as.data.frame(cell_prop)
colnames(cell_prop_df) = c("Sample","cluster","cell_ratio")
### cell proportion for each sample ###
cell_prop_df = merge(cell_prop_df, sample_tb, by="Sample")

kp_tb = cell_prop_df %>% filter(Lineage=="c57")
tmp1 = kp_tb %>% group_by(cluster,Tumor) %>% summarise(mean_ratio = mean(cell_ratio))
tmp = tmp1 %>% pivot_wider(id_cols = cluster, names_from = Tumor, values_from = mean_ratio)
tmp$fc = tmp$Primary / tmp$Baseline
tmp$log2fc = log2(tmp$fc)

plot_df = tmp

plot_df = plot_df %>%
  mutate(class = ifelse( sign(log2fc)>0, "up","down" )) %>%
  mutate(cluster2 = factor(cls_names[cluster], levels =rev(cls_names)))

p.kp = plot_df %>% 
  ggplot(., aes(x=log2fc, y= cluster2, fill=class)) +
  geom_bar(stat ="identity")+
  theme_classic() +
  scale_fill_manual(values = c("steelblue3", "red2"))  
p.kp

### met 
kp_tb = cell_prop_df %>% filter(Lineage=="fvb")
tmp1 = kp_tb %>% group_by(cluster,Tumor) %>% summarise(mean_ratio = mean(cell_ratio))
tmp = tmp1 %>% pivot_wider(id_cols = cluster, names_from = Tumor, values_from = mean_ratio)
tmp$fc = tmp$Metastatic / tmp$Baseline
tmp$log2fc = log2(tmp$fc)

plot_df = tmp


plot_df = plot_df %>%
  mutate(class = ifelse( sign(log2fc)>0, "up","down" ))  %>%
  mutate(cluster2 = factor(cls_names[cluster], levels = rev(cls_names)))

p.met = plot_df %>% 
  ggplot(., aes(x=log2fc, y= cluster2, fill=class)) +
  geom_bar(stat ="identity")+
  theme_classic() +
  scale_fill_manual(values = c("steelblue3", "red2"))
p.met

p = p.kp | p.met
p
ggsave(p, filename = "KP_vs_MMTV_compare_cell_prop_sample_level.pdf", width = 7, height = 3)


mean_prop = cell_prop_df %>% group_by(cluster,Time) %>% summarise(mean_ratio = mean(cell_ratio))

#### 
tmp = pivot_wider(mean_prop, id_cols = cluster, names_from = Time, values_from = mean_ratio)

cal_logfc  <- function(y, x){
   res = log2(1e-3 + y) - log2(x + 1e-3) 
   res
}

kp_fc = data.frame(cluster = tmp$cluster, 
                   t1 = cal_logfc(tmp$KP_W04, tmp$KP_Baseline),
                   t2 = cal_logfc(tmp$KP_W08, tmp$KP_Baseline),
                   t3 = cal_logfc(tmp$KP_W12, tmp$KP_Baseline))
met_fc = data.frame(cluster = tmp$cluster, 
                   t1 = cal_logfc(tmp$MT_W08, tmp$MT_Baseline),
                   t2 = cal_logfc(tmp$MT_W12, tmp$MT_Baseline),
                   t3 = cal_logfc(tmp$MT_W16, tmp$MT_Baseline))

p1 = kp_fc %>% 
  pivot_longer(2:4, names_to = "time") %>% 
  ggplot() +
  geom_point(aes(x=value, y=cluster, color=time)) +
  xlim(-2,2)

p2 = met_fc %>% 
  pivot_longer(2:4, names_to = "time") %>% 
  ggplot() +
  geom_point(aes(x=value, y=cluster, color=time)) +
  xlim(-2,2)

p1 | p2

plot_single_cls_sample <- function(ptime.kp, ptime.met, test_cls){
  tmp1 = ptime.kp$data %>% filter(cluster==test_cls) %>% mutate(type="KP")
  tmp2 = ptime.met$data %>% filter(cluster==test_cls) %>% mutate(type="MMTV")
  df = rbind(tmp1,tmp2) 
  
  p = df %>%
    ggplot(aes(x=time,y=value, fill=type)) +
    geom_bar(stat="identity") +
    facet_grid(.~type) +
    theme_classic() +
    scale_fill_manual(values = c("orange3","steelblue3"))+
    theme(legend.position = "none") +
    ggtitle(test_cls)+
    ylab(" Change of subset abundance within neutrophils log2FC(tumor/normal) ")
  p
}

p01 = plot_single_cls_sample(p1, p2, test_cls = "N9") +
  ylim(-2,2)

p02 = plot_single_cls_sample(p1, p2, test_cls = "N8") +
  ylim(-1.5, 1.5)

p03= plot_single_cls_sample(p1, p2, test_cls = "N4") +
  ylim(-1.8,1.8)


p= p01|p02| p03
#ggsave(p, filename = "./KP_vs_MMTV_Mean_subset_prop_over_normal_mice.pdf",width=5,height = 2.5)


p04= plot_single_cls_sample(p1, p2, test_cls = "N2") +
  ylim(NA,1)
p04
p = p04| p03| p01
ggsave(p, filename = "./KP_vs_MMTV_Mean_subset_prop_over_normal_mice_Mmp8.pdf",width=5,height = 2.5)



########### Changes across time ############
## Cd74+ neutrophils

#MHC class II protein complex binding
mhc2_gene= get_genes_from_GO_mouse("GO:0023026")

#H2-Ab1: major histocompatibility complex, class II, DQ beta 1
example_gene = c("Cd74", "H2-Ab1","H2-Aa")

FeaturePlot(neutro, example_gene, ncol = 3)
DotPlot(neutro,features = example_gene, group.by = "Tumor")

neutro = AddModuleScore(neutro, features = list(mhc1_gene), name = "MHCI_protein_binding")
neutro = AddModuleScore(neutro, features = list(antigen_present_go), name = "Antigen_presenting")
neutro = AddModuleScore(neutro, features = list(mhc2_gene), name = "MHCII_complex_binding")



p1 = DotPlot(neutro, features = c("Antigen_presenting1", "MHCII_complex_binding1", "MHCI_protein_binding1"), group.by = "cellanno_L3")

p1 = DotPlot(neutro, features = c("Antigen_presenting1", "MHCII_complex_binding1", "MHCI_protein_binding1"), group.by = "Time")
p1


# boxplot of 
VlnPlot(neutro, features = "MHCII_complex_binding1", group.by = "Time", pt.size = 0) + geom_boxplot(width=0.1, fill="white", outlier.shape = NA)

subset(neutro, cellanno_L3=="N8") %>%
  VlnPlot(., features = "MHCII_complex_binding1", group.by = "Time", pt.size = 0) + geom_boxplot(width=0.1, fill="white", outlier.shape = NA)

neu_ag = subset(neutro, cellanno_L3=="N8")

p1 = VlnPlot(neu_ag, features = "MHCII_complex_binding1", group.by = "Time") +
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA)

df_summary = p1$data %>% group_by(ident) %>%
  summarise(mean = mean(MHCII_complex_binding1))

kp_data = p1$data %>%
  filter(str_detect(ident,"KP"))

mt_data = p1$data %>%
  filter(str_detect(ident,"MT"))

p.kp = kp_data %>%
  ggplot(aes(x=ident, y=MHCII_complex_binding1, fill=ident)) +
  geom_boxplot(  outlier.shape = NA ) +
  theme_classic() +
  scale_fill_manual(values = c( "Grey80", RColorBrewer::brewer.pal(9,"Oranges")[c(5,7,9)]) ) +
  NoLegend() +
  #geom_hline(yintercept = as.numeric(df_summary[1,2]), linetype="dashed", alpha=0.5) +
  ylim(c(-0.4,2.7)) 
p.kp

p.mt = mt_data %>%
  ggplot(aes(x=ident, y=MHCII_complex_binding1, fill=ident)) +
  geom_boxplot( outlier.shape = NA ) +
  theme_classic() +
  scale_fill_manual(values = c( "Grey80", RColorBrewer::brewer.pal(9,"Blues")[c(5,7,9)]) ) +
  NoLegend()  +
  #geom_hline(yintercept =  as.numeric(df_summary[5,2]), linetype="dashed", alpha=0.5) +
  ylim(c(-0.4,2.7)) 
p = p.kp | p.mt
p
ggsave(p, filename = "figure/003.KP_MT_Cd74_score.pdf", width=4, height = 2.5)
write.csv(df_summary, file = "figure/003.KP_MT_cd74.data.csv")




######## Analysis of immunosuppressive neutrophil #########

ref_immun_sup_df = readxl::read_xlsx("/raid1/aichen/projects/scLung/00.neutrophil_described/Gong_RNAseq_immunosuppress_gene.xlsx", skip = 1, col_names = "gene")
head(ref_immun_sup_df)
ref_immun_sup_df$term = "Gong_etal_Immunosuppressive_Neutrophil"
ref_immun_sup_df = ref_immun_sup_df[,c("term","gene")]


neutro <- AddModuleScore(neutro, features = list(immunosup = ref_immun_sup_df$gene), name="Gong_etal_ImmS")
p0= plot_density(neutro, "Gong_etal_ImmS1", method = "wkde")
ggsave(p0, filename = "/raid1/aichen/projects/scLung/00.neutrophil_described/figures/003.a.UMAP_Gongetal_ImmunSup.pdf", width = 4, height = 3)

mat_neu = subset(neutro, cellanno_L3 =="N5")

### comparison across time

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



