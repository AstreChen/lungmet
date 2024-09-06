library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(scales)
library(ggsci)



ing_neu <- readr::read_rds("./Integrated_seurat_obj_MNN_v2.rds.gz")
ing_neu@meta.data = ing_neu@meta.data %>%
  mutate(cancer = ifelse(type=="Naive", "None", cancer))
ing_neu$tissue = factor(ing_neu$tissue, levels = c("BM","Blood","Spleen","Lung","Breast","Tumour"))

ing_neu$type = factor(ing_neu$type, levels = c("Naive","Primary","Metastatic"))
p = DimPlot(ing_neu, group.by = "type", label=F) +
  scale_color_manual(values = c(
                                 brewer.pal(n=9, name="Purples")[3], #Naive
                                brewer.pal(n=6, name="Blues")[5],
                               brewer.pal(n=9, name="YlOrBr")[6] # metastasis
                                ))
ggsave(p, filename = "fig4.01.umap_plots//All_cancer_type_umap.pdf", width = 8, height = 5.5)


ing_neu$ClusterName = factor(ing_neu$ClusterName, levels = c( "preNeu1"  ,       "preNeu2" ,        "IMM1"     ,       "IMM2" ,           "MAT1", "MAT2"        ,    "MAT3"   ,         "MAT4"   ,         "MAT5"  ,  "Tissue IMM"   ,   "Tissue MAT",      "Primary TAN",     "Metastatic IMM" ,"Metastatic MAT1", "Metastatic MAT2"          ))

color_code_cls = c(
   brewer.pal(n=4, name="Greys")[c(2,3)], brewer.pal(n=9, name="Purples")[2:6], brewer.pal(n=9, name="BuPu")[2:5], # from preNeu1 to MAT5
  brewer.pal(n=6, name="Blues")[5],  #"Tissue IMM"   ,   "Tissue MAT",      "Primary TAN",  
  brewer.pal(n=9, name="YlOrBr")[c(5,7,9)]
)


p = DimPlot(ing_neu, group.by = "ClusterName", label=F) + NoAxes() +
  scale_color_manual(values = color_code_cls)
p
ggsave(p, filename = "fig4.01.umap_plots/All_ClusterName.pdf", width = 8, height = 5.5)

## 3 tissues
p01 = DimPlot(subset(ing_neu, cancer=="Breast_PyMT"), group.by = "ClusterName" , split.by = "tissue")  + NoAxes() +
  scale_color_manual(values = color_code_cls)
ggsave(p01, filename = "fig4.01.umap_plots/Breast_PyMT_tissues.pdf", width=10, height = 3)


### 4T1 - all tissue
p02 = DimPlot(subset(ing_neu, cancer=="Breast_4T1") %>% subset(., type=="Metastatic"), group.by = "ClusterName" , split.by = "tissue")  + NoAxes() +
  scale_color_manual(values = color_code_cls)
ggsave(p02, filename = "fig4.01.umap_plots/Breast_4T1_tissues.pdf", width=10, height = 3)

### PDAC
p03 = DimPlot(subset(ing_neu, cancer=="PDAC"), group.by = "ClusterName" , split.by = "tissue")  + NoAxes() +
  scale_color_manual(values = color_code_cls)
ggsave(p03, filename = "fig4.01.umap_plots/PDAC_tissues.pdf", width=13, height = 3)


p04 = DimPlot(subset(ing_neu, cancer=="None"), group.by = "ClusterName" , split.by = "tissue")  + NoAxes() +
  scale_color_manual(values = color_code_cls)
ggsave(p04, filename = "fig4.01.umap_plots/Naive_tissues.pdf", width=13, height = 3)


### plot PyMT_Breast
pymt_neu = subset(ing_neu, dataset=="PyMT") %>% subset(cancer !="None")
p.1 = DimPlot(pymt_neu, group.by = "ClusterName" , split.by = "tissue")  + NoAxes() +
  scale_color_manual(values = color_code_cls)
ggsave(p.1, filename = "./fig4.01.umap_plots/pymt_cancer.pdf", width = 7, height = 3)

tmp_neu = subset(ing_neu, dataset =="4T1_Ren") %>% subset(cancer !="None")
p.1 = DimPlot(tmp_neu, group.by = "ClusterName" , split.by = "tissue")  + NoAxes() +
  scale_color_manual(values = color_code_cls)
ggsave(p.1, filename = "./fig4.01.umap_plots/4T1_cancer.pdf", width = 10, height = 3)

### plot KP lung
kp_neu = subset(ing_neu, dataset %in% c("KP","PyMT") ) %>% subset(sample !="PyMT_Breast_M01")
p.1 = DimPlot(kp_neu, group.by = "ClusterName" , split.by = "type")  + NoAxes() +
  scale_color_manual(values = color_code_cls)
ggsave(p.1, filename = "./fig4.01.umap_plots/for_KP_cancer.pdf", width = 10, height = 3)


ing_neu$cancer_tissue = paste0(ing_neu$cancer,".",ing_neu$tissue)


### plot pseudotime ####
to_change = c("preNeu1" = "preNeu",
  "preNeu2" = "preNeu",
  "IMM1" = "IMM",
  "IMM2" = "IMM",
  "MAT1" = "MAT",
  "MAT2" = "MAT",
  "MAT3" = "MAT",
  "MAT4" = "MAT",
  "MAT5" = "MAT",
  "Tissue IMM" =  "Tissue IMM",
  "Tissue MAT" = "Tissue MAT",
  "Primary TAN" = "Primary TAN",
  "Metastatic IMM" = "Metastatic IMM" ,
  "Metastatic MAT1" = "Metastatic MAT1",
  "Metastatic MAT2"=  "Metastatic MAT2")

meta_df = ing_neu@meta.data
meta_df$ClusterNameMarker_coarse =to_change[ as.character( ing_neu$ClusterName)]
meta_df$ClusterNameMarker_coarse = factor(meta_df$ClusterNameMarker_coarse, levels = unique(to_change))
ing_neu@meta.data = meta_df
color_code = c( rev(brewer.pal(n=9, name="YlOrBr")[c(5,7,9)]), ### Met MAT2, MAT1, IMM
                brewer.pal(n=6, name="Blues")[5], rev(brewer.pal(n=9, name="BuPu")[4:5]),
                rev(brewer.pal(n=9, name="Purples")[2:4])) ##


pstime = read.csv("outdir_10x/CytoTRACE_output_ct_pseudotime.csv", row.names = 1)
pstime = pstime[, c("ct_num_exp_genes",  "ct_score" ,"ct_pseudotime")]
ing_neu = AddMetaData(ing_neu, pstime)

p = VlnPlot(ing_neu, group.by = "ClusterNameMarker_coarse", features = "ct_pseudotime", pt.size = 0) +
  scale_fill_manual(values = rev(color_code))

p = p$data %>%
  ggplot() +
  geom_boxplot(aes(x=ct_pseudotime, y=ident, fill=ident ), outlier.shape = NA) +
  theme_classic() +
  scale_fill_manual(values = rev(color_code))

p = RidgePlot(ing_neu, group.by = "ClusterNameMarker_coarse", features = "ct_pseudotime") +     scale_fill_manual(values = rev(color_code))

ggsave(p, filename = "fig4.01.umap_plots/CytoTRACE_ridge_plot.pdf", width = 6, height = 3.5)



########## compare neutrophil abundance in the primary or metastatic tumor sites #########

### get sample table:
sample_tb = data.frame(ing_neu@meta.data[,c("sample","tissue","type","cancer_tissue")]) %>% unique() %>% as.tibble()
total_num = ing_neu$sample %>% table()  %>% data.frame()

mat_cls = table(ing_neu$sample, ing_neu$ClusterName) %>% as.matrix()
cell_prop = mat_cls / rowSums(mat_cls)


### assign met group: 
met_group = sample_tb %>% filter(tissue %in% c("Lung","Breast","Pancreatic")) %>%
  filter(type=="Metastatic")

pri_group = sample_tb %>% filter(tissue %in% c("Lung","Breast","Pancreatic")) %>%
  filter(type=="Primary")

met1 = cell_prop[met_group$sample,] %>% colMeans()
pri1 = cell_prop[pri_group$sample,] %>% colMeans()

####### compare metastasis vs primary #########
## calculate p-value for each cell type
pvalue.vec = lapply(1:ncol(cell_prop), FUN = function(celltype){
  res = t.test(x = cell_prop[met_group$sample, celltype],  y=cell_prop[pri_group$sample, celltype] )
  return(res$p.value)
})

plot_df = data.frame(log2fc = log2(met1/pri1), mean_ratio = 100*(met1 + pri1)/2 )
plot_df$cell_name = rownames(plot_df)
plot_df$cell_name = factor(plot_df$cell_name, levels = levels(ing_neu$ClusterName))
plot_df$pval = as.numeric(pvalue.vec) %>% round(., digits = 4)
plot_df = plot_df %>%
  mutate(signif.label = case_when(pval >= 0.05 ~ "",
                                  pval > 0.01 ~ "*",
                                  pval > 0.001 ~"**",
                                  pval > 0.0001 ~"***",
                                  T ~ "****")) %>%
  mutate(signif.num = ifelse(pval<0.05, paste0("P=",pval), ""))

p = ggplot(plot_df, aes(x=log2fc, y=cell_name, size=mean_ratio, color=log2fc  )) +
  geom_point() +
  geom_text(aes(x=log2fc, y=cell_name, label=signif.num,vjust=0,hjust=-0.2, size=10)) +
  scale_color_distiller(type="seq", palette = "RdBu",direction = -1 ) +
  theme_classic() +
  geom_vline(xintercept = 0,  alpha=0.3)
  
p
ggsave(p, filename = "./fig4.03.composition/1.Metastasis_vs_Primary_composition_Lung_Breast_Pancreas.pdf", width = 5, height = 4)



######## compare DEGs of primary TAN and lung MANs ######
Idents(ing_neu) = c("ClusterName")
deg = FindMarkers(ing_neu, ident.1="Primary TAN", ident.2=c( 'Metastatic IMM','Metastatic MAT1', 'Metastatic MAT2'))
deg$gene = rownames(deg)
deg = deg %>% mutate(mean_pct = (pct.1 + pct.2)/2)
deg = deg %>% filter(mean_pct>0.10)
write.csv(deg, "./04.Integrate_all_neutrophil/DEG_analysis/DEG_pri_vs_met_meanPct_mt0.1.csv")

deg %>% arrange(-avg_log2FC) %>% filter(pct.1>0.3) %>% head(n=50)

sub_ing_neu = subset(ing_neu, ClusterName %in% c("Primary TAN",'Metastatic IMM','Metastatic MAT1', 'Metastatic MAT2'))
sub_ing_neu$group = ifelse(sub_ing_neu$ClusterName =="Primary TAN","primary","metastatic")

avg_exp = AverageExpression(sub_ing_neu, assays = "RNA", group.by = "group", return.seurat = F)
avg_exp

avg_exp = avg_exp$RNA
avg_exp = data.frame(avg_exp)
avg_exp = log1p(avg_exp)


### to label genes
deg %>% arrange(-avg_log2FC) %>% filter(pct.1>0.3) %>% head(n=50) %>% rownames()
deg %>% arrange(avg_log2FC) %>% filter(pct.2>0.3) %>% head(n=50) %>% rownames()

pri_spec = c("Ccl3","Ccl4","Cd274","Spp1","Tnf","Vegfa","Ptgs2","Csf1","Bhlhe40")
met_spec = c("Stfa2","Cstdc5","Stfa3","Mmp8","Anxa1","Wfdc21","Padi4","Prok2","Mapk13")
to_label = avg_exp %>% rownames_to_column(var="gene") %>%
  filter(gene %in% c(pri_spec,met_spec))

to_color_point = deg %>%
  mutate(pct.diff = pct.1 - pct.2) %>%
  filter(abs(avg_log2FC)>2) %>% filter(p_val_adj < 1e-5) %>%
  filter( abs(pct.diff) >= 0.2 )
to_color_point_df = avg_exp[to_color_point$gene,]


library(ggrepel)
p = ggplot(avg_exp) +
  geom_point(aes(x=primary, y=metastatic), color="grey60", alpha=0.8) +
  geom_point(data=to_color_point_df, aes(x=primary, y=metastatic), color="red") +
  geom_abline(slope = 1, linetype="dashed")+
  geom_text_repel(
    data=to_label,
    aes(x=primary, y=metastatic, label=gene),
    color="red",
    box.padding = unit(0.35, "lines"),  
    point.padding = unit(0.5, "lines"), 
    segment.color = 'grey50',           
    segment.size = 0.5,                 
    direction = 'y',                   
    hjust = 0, vjust = 0        ,       
    nudge_x = 0.2, nudge_y =0.6
  ) +
  theme_classic()
p
ggsave(p, filename = "./04.Integrate_all_neutrophil/10.figure4_DEG_primary_TAN_vs_metastatic_Neu.pdf", width = 5, height = 4)


##########################################
##### Find metastasis-specific genes #####
##########################################
# seurat version 5
pseudo_bulk = AggregateExpression(ing_neu, assays="RNA", return.seurat = T, group.by = c("dataset","tissue","sample","type","cancer"))


#### 1. Call Metastatic signatures by pseudobulk in Discovery Cohort ###
#### lung compare ###
obj = subset(pseudo_bulk, tissue=="Lung")
obj$orig.ident = colnames(obj)
Idents(obj) = "type"
bulk.de.naive = FindMarkers(obj, ident.1 = "Metastatic", ident.2 = "Naive",  test.use = "DESeq2", only.pos = T )
bulk.de.prim = FindMarkers(obj, ident.1 = "Metastatic", ident.2 = "Primary",  test.use = "DESeq2", only.pos = T )

deg.naive = bulk.de.naive %>% filter(avg_log2FC>2) %>% filter(p_val_adj<0.05) %>% rownames()
deg.prim = bulk.de.prim %>% filter(avg_log2FC>2) %>% filter(p_val_adj<0.05) %>% row.names()

lung_met_genes = intersect(deg.naive,deg.prim)

get_pct_df = function(obj, features, group.by){
  p = DotPlot(obj, features = lung_met_genes, group.by = "type_cancer")
  pct.df = p$data
  pct.df = pivot_wider(pct.df, id_cols=features.plot, names_from = id, values_from = pct.exp)
  return(pct.df)
}

solid_neu = subset(ing_neu, tissue %in% c( "Lung", "Breast", "Tumour" ) )
solid_neu$type_cancer = paste0(solid_neu$type,"_",solid_neu$cancer)

avg_df = AverageExpression(solid_neu, features = lung_met_genes, group.by = "type_cancer", assays = "RNA")$RNA %>% as.data.frame()
colnames(avg_df) = str_replace(colnames(avg_df), pattern = "-", replacement = "_") %>% str_replace(., pattern = "-", replacement = "_")

pct_df = get_pct_df(solid_neu, features = lung_met_genes, group.by="type_cancer")
pct_df$gene = as.character(pct_df$features.plot)
pct_df$features.plot=NULL 


lung_met_gene_df = avg_df %>%
  rownames_to_column(var = "gene") %>%
  inner_join(., y=pct_df, by="gene",   suffix=c(".avgExp",".pct"))

re_cols = sort(colnames(lung_met_gene_df[,2:ncol(lung_met_gene_df)]))
lung_met_gene_df = lung_met_gene_df[,c("gene", re_cols)]

lung_met_gene_df$ratio_PyMT_Met_vs_brca = lung_met_gene_df$Metastatic_Breast_PyMT.avgExp / lung_met_gene_df$Primary_Breast_PyMT.avgExp

lung_met_gene_df = lung_met_gene_df %>% arrange(-ratio_PyMT_Met_vs_brca) %>%
  relocate(ratio_PyMT_Met_vs_brca, .after = gene) %>%
  arrange(-Metastatic_Breast_PyMT.avgExp)

write.csv(lung_met_gene_df, "./callMetSignature/Metastatic_signatures_lung_tb.csv")

p2 = DotPlot(ing_neu, features = lung_met_gene_df$gene, group.by="cancer_tissue") + RotatedAxis()
tmp = p2$data %>%
  filter(str_detect(id, pattern = "PDAC|None|LUAD|PyMT.Breast")) 
level1_gene = tmp %>%
  group_by(features.plot) %>%
  summarise(max.avg.exp = max(avg.exp), max.pct = max(pct.exp)) %>%
  filter(max.pct < 60) %>%
  pull(features.plot) %>% as.character()


### 2. call DEG in validation cohort #####

b16_neu <- read_rds("../../data/2.our_B16_Neu/B16_Lung_Neutrophil.rds.gz")
neu_AT3 <- read_rds("../../data/Lung_AT3_Ren/Lung_AT3-gcsf_Ren_BrLM_clean.rds.gz") ### 包含三种组织 BM，Blood， Lung
neu_AT3_lung = subset(neu_AT3, tissue=="Lung") %>% JoinLayers()
neu_m39m <- read_rds("../../data/Lung_M39M_Kaplan/Lung_M39M_Kaplan_Neutro.rds.gz") 

Idents(b16_neu) = "Time"
deg_b16 <- FindMarkers(b16_neu, ident.1 = "B16_W04", ident.2 = "B16_Baseline")

Idents(neu_AT3_lung) = "Type"
deg_AT3 <- FindMarkers(neu_AT3_lung, ident.1 = "Metastatic", ident.2 = "Naive")

Idents(neu_m39m) = "type"
deg_m39m <- FindMarkers(neu_m39m, ident.1 = "Tumor_NoGEMy", ident.2 = "NoTumor_NoGEMy")

### list validate gene sets: cutoff 1.0
validate.list = list(B16=deg_b16, AT3=deg_AT3, M39M=deg_m39m) %>%
  lapply(., FUN = function(x){ 
   y = x %>% filter(avg_log2FC>1.0)
   rownames(y)
  })

### calculate the overlap set number for each lung metastasis gene: 
overlap_number_list = list()
for(g in level1_gene){  
  res = sapply(validate.list, FUN = function(to_check_lst){
    tmp = g %in% to_check_lst
    as.numeric(tmp)
  })
  overlap_number_list[[g]] = sum(res)
}
overlap_number_res = unlist(overlap_number_list)
print("All validated genes:")
final_sig_genes = overlap_number_res[overlap_number_res >=2] %>% names()
final_sig_genes

### 3 : identify neutrophil-specific genes  #########
### other immune cells
imm_obj = read_rds("../../../Lib009/clustering/01.KP_PyMT_All_TME_object.rds")
imm_obj = subset(imm_obj, Type=="MMTV") ## only metastasis
imm_obj = subset(imm_obj, cellanno_L1 %in% c("T_NK","B","Myeloid") )

p1 = DotPlot(imm_obj, features = final_sig_genes, group.by="cellanno_L2") + RotatedAxis()
tmp.max = p1$data %>% group_by(features.plot) %>%
  summarise(max.avg.exp = max(avg.exp), max.pct = max(pct.exp)) 

final_gene = p1$data %>% filter(id =="Neutrophil") %>%
  inner_join(tmp.max) %>%
  filter(avg.exp >= max.avg.exp) %>%
  filter(pct.exp >= max.pct) %>%
  filter(avg.exp>=1) %>%
  pull(features.plot) %>% as.character()
final_gene

final_gene_df = data.frame(gene=final_gene)

write.csv(final_gene_df, file = "../final_Neutrophils_metastatic_signature_gene.csv")



final_gene_df = read.csv("../final_Neutrophils_metastatic_signature_gene.csv", row.names = 1)
show_genes =final_gene_df$gene

p = DotPlot(ing_neu,features = show_genes, group.by ="cancer_tissue", cols = "RdBu") + RotatedAxis() 
ggsave(p, filename = "./fig4.02.met_gene/DotPlot_show_genes.pdf", width=6, height = 4.5)
ing_neu= AddModuleScore(ing_neu, features = list(show_genes), name = "MetSig")


### Dotplot for these genes:  lung tissue ###
sub_neu = subset(ing_neu, subset = tissue %in% c("Lung") )
p = DotPlot(sub_neu,features = show_genes, group.by ="cancer_tissue", cols = "RdBu") + RotatedAxis() # 
ggsave(p, filename="./fig4.02.met_gene/Dot_plot_discovery_Lung_tissue.pdf", width = 6, height = 3)


### Dotplot for these genes:  blood  ###
sub_neu = subset(ing_neu, subset = tissue %in% c("Blood") )
sub_neu$sample = factor(sub_neu$sample, levels = rev(c("wt.ctl.pb2","Naive_Blood", "Blood_Rep1","Blood_Rep2","Metastatic_Blood")))
p = DotPlot(sub_neu,features = show_genes, group.by ="sample", cols = "RdBu") + RotatedAxis() # 
ggsave(p, filename="./fig4.02.met_gene/Dot_plot_discovery_blood_tissue.pdf", width = 6, height = 3)

neu_AT3_blood = subset(neu_AT3, tissue=="Blood")
p2 = DotPlot(neu_AT3_blood,features = show_genes, group.by ="type", cols = "RdBu") + RotatedAxis() # 
ggsave(p2, filename="./fig4.02.met_gene/Dot_plot_AT3_blood_tissue.pdf", width = 6, height = 2)

### Featureplot for these genes  ###
p1 = FeaturePlot(ing_neu, features = "MetSig1", split.by = "type", cols = c("Grey90", brewer.pal(n=3, "Reds")[3] ) , min.cutoff = "q10")  
met_neu = subset(ing_neu, type=="Metastatic")
p2 = FeaturePlot(met_neu, features = "MetSig1", split.by = "tissue", cols = c("Grey90", brewer.pal(n=3, "Reds")[3] ), min.cutoff = "q10" )
p2

p = FeaturePlot(ing_neu, features = "MetSig1",  cols = c("Grey90", brewer.pal(n=3, "Reds")[3] ), min.cutoff = "q10" )
ggsave(p, filename = "fig4.02.met_gene/MetSig1_gene_norm_pri_met_merged.pdf", width = 5, height = 3)

ggsave(p1, filename = "fig4.02.met_gene/MetSig1_gene_norm_pri_met.pdf", width = 10, height = 3)
ggsave(p2, filename = "fig4.02.met_gene/MetSig1_MetNeu_tissue.pdf", width = 13, height = 3)
