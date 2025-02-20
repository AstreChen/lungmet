library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(ggpubr)
library(scales)
library(ggsci)
library(igraph)



### compare TME composition ####
meta_tb = read.csv("/raid1/aichen/projects/scLung/01.KP_PyMT_All_TME_MetaData.csv",row.names = 1)
meta_tb = meta_tb %>% filter(cellanno_L2!="Mast") %>%
  filter(CellSorting!="Cd45+")


meta_tb2 = meta_tb %>% 
  mutate(anno = ifelse(cellanno_L1 %in% c("Epithelial",  "Endothelial", "Stromal"), cellanno_L1, cellanno_L2 )) %>%
  mutate(anno = ifelse(cellanno_L2=="DC", cellanno_L3, anno)) %>%
  mutate(anno = ifelse(cellanno_L1=="T_NK", cellanno_L3, anno)) %>%
  mutate(anno = str_replace(anno, pattern ="T_CD4_Teff.Icos", replacement = "T_CD4_Teff" )) %>%
  mutate(anno = str_replace(anno, pattern ="T_CD8_Teff.Gzmk", replacement = "T_CD8_Teff" )) %>%
  mutate(anno = str_replace(anno, pattern ="T_CD4_ISG", replacement = "T_CD4_Teff" )) %>%
  mutate(anno = str_replace(anno, pattern ="T_CD4_Tem.Tox", replacement = "T_CD4_Teff" )) %>%
  mutate(anno = str_replace(anno, pattern ="T_CD8_Texh.Pdcd1", replacement = "T_CD8_Tex" )) 


meta_tb2$CD45_group = ifelse(meta_tb2$cellanno_L1 %in% c("Endothelial","Epithelial","Stromal"), "CD45_neg", "CD45_pos")


## type == PyMT-MMTV #####
meta_df = meta_tb2 %>% filter(Lineage=="fvb")
meta_df$Time_sample = paste0(meta_df$Time, ".", meta_df$Sample)
meta_df_list = split.data.frame(meta_df, meta_df$CD45_group)

#####################################
####### cell module of metastasis ###
#####################################
## get count matrix:
# ----------CD45pos
cell_count_df = meta_df_list$CD45_pos %>%
  group_by(Time_sample, anno) %>%
  summarise(cell_count = n()) %>%
  pivot_wider(data = ., id_cols = Time_sample, names_from = anno, values_from = cell_count, values_fill = 0 )

cell_count_mat = data.frame(cell_count_df[,-1])
rownames(cell_count_mat) = cell_count_df$Time_sample
tmp = table(meta_df_list$CD45_pos$Time_sample) %>% as.data.frame()
rownames(tmp) = tmp$Var1

cell_count_mat = cell_count_mat / tmp$Freq
cell_count_mat_cd45pos = cell_count_mat

#----------CD45 neg
cell_count_df = meta_df_list$CD45_neg %>%
  group_by(Time_sample, anno) %>%
  summarise(cell_count = n()) %>%
  pivot_wider(data = ., id_cols = Time_sample, names_from = anno, values_from = cell_count, values_fill = 0 )

cell_count_mat = data.frame(cell_count_df[,-1], check.names = F)
rownames(cell_count_mat) = cell_count_df$Time_sample
tmp = table(meta_df_list$CD45_neg$Time_sample) %>% data.frame(check.names = F)
rownames(tmp) = tmp$Var1
cell_count_mat = cell_count_mat / tmp$Freq
cell_count_mat_cd45neg = cell_count_mat

cell_count_mat_all = cbind(cell_count_mat_cd45pos, cell_count_mat_cd45neg)

cd45all_sample = tmp %>% 
  filter(!( as.character(Var1) %in% c("MT_W08.we1","MT_W08.we2","MT_W08.we3"))) %>% pull(Var1) %>% as.character()

### remove CD45+only sample
cor_mat = cor(cell_count_mat_all[cd45all_sample,]) #[, anno_df$gep_name])

write.csv(cor_mat, "./007_mmtv_res_corMatrix_cell_module.csv")
write.csv(cell_count_mat_all[cd45all_sample,],"./007_mmtv_res_cell_prop_matrix_all_Cd45_1o1.csv")


################################################
#### Compare cell proportions of  MMTV-PyMT ####
################################################

mean_prop.met = colMeans(cell_count_mat_all[4:12,]) ### Metastasis samples
mean_prop.diff = log2(mean_prop.met) - log2(colMeans(cell_count_mat_all[1:3,])) ### Control samples

df = data.frame(mean_prop.diff)
df$name = factor(rownames(df), levels = rownames(df)[order(df$mean_prop.diff)] )
df$mean_prop.met = mean_prop.met
p = ggplot(df) +
  geom_point(aes(x=mean_prop.diff,y=name, size=mean_prop.met, color=mean_prop.diff)) +
  scale_color_distiller(type="seq", palette = "RdBu",direction = -1 ) +
  scale_size_continuous(limits = c(0,0.5)) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype="dashed",  alpha=0.3) +
  xlab(label="Proportion log2(Metastasis/Naive)")
p
ggsave(p, filename = "./mmtv_celltype_enrichment_plot.pdf", width = 4, height = 6)

### calculate p-value:
pvalue.vec = lapply(1:ncol(cell_count_mat_all), FUN = function(celltype){
  res = t.test(x = cell_count_mat_all[1:3, celltype],
               y=cell_count_mat_all[4:12, celltype] )
  return(res$p.value)
})
names(pvalue.vec) = colnames(cell_count_mat_all)
pvalue.vec = unlist(pvalue.vec)
pvalue.vec[pvalue.vec<0.05]


#############################################
### Calculate correlation of cell types #####
#############################################

p=Heatmap(cor_mat, #row_km = 5, column_km = 5, 
          cluster_rows = T, cluster_columns=T, 
          clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
          row_names_gp = gpar(fontsize=10), column_names_gp = gpar(fontsize=10) )
od =  row_order(p)
p

pdf(file="./mmtv_cellmodule_heatmap.pdf", width=10, height=9)
print(p)
dev.off()


tmp = cor_mat#[od[[3]], od[[3]]]
tmp[tmp==1] =0
tmp[tmp >= 0.5] =1
tmp[tmp< 0.5]=0
adj = tmp

g = graph_from_adjacency_matrix(adj, mode="undirected" , weighted = T)
plot.igraph(g)


vertex_df = data.frame(node = rownames(cor_mat) ) 
rownames(vertex_df) = vertex_df$node


unique(vertex_df$cellanno_L2)



col_vec = c(pal_d3(palette = "category20")(20), "#DE9ED6FF")
names(col_vec) = unique( vertex_df$node)
vertex_df$color = col_vec[vertex_df$node ]

### cellprop:
mean_cell_prop = colMeans(cell_count_mat_all)
vertex_df$cell_prop = mean_cell_prop[vertex_df$node]
vertex_df$cell_prop_plot = ifelse(vertex_df$cell_prop>=0.1, 0.1, vertex_df$cell_prop)

set.seed(1001)
g3 = subgraph(g, vertex_df$node )
vertex_df = vertex_df[names(V(g3)),]
pdf(file = "./MMTV_cell_module_network.pdf", width = 8, height = 8)
p=plot(g, layout=layout_with_fr,
       edge.curved = 0.2,
       vertex.color =  vertex_df$color,
       vertex.size= 40* sqrt(vertex_df$cell_prop_plot), 
       edge.width= E(g3)$weight *2,
       edge.color = "grey",
       label.cex = 6,
       
)
p
print(p)
dev.off()




################################## KP ###############################
## type == KP #####
meta_df = meta_tb2 %>% filter(Lineage=="c57")
meta_df$Time_sample = paste0(meta_df$Time, ".", meta_df$Sample)
meta_df_list = split.data.frame(meta_df, meta_df$CD45_group)

##############################################
####### cell module of primary lung tumors ###
##############################################
## get count matrix:
# ----------CD45pos
cell_count_df = meta_df_list$CD45_pos %>%
  group_by(Time_sample, anno) %>%
  summarise(cell_count = n()) %>%
  pivot_wider(data = ., id_cols = Time_sample, names_from = anno, values_from = cell_count, values_fill = 0 )

cell_count_mat = data.frame(cell_count_df[,-1])
rownames(cell_count_mat) = cell_count_df$Time_sample
tmp = table(meta_df_list$CD45_pos$Time_sample) %>% as.data.frame()
rownames(tmp) = tmp$Var1

cell_count_mat = cell_count_mat / tmp$Freq
cell_count_mat_cd45pos = cell_count_mat

#----------CD45 neg
cell_count_df = meta_df_list$CD45_neg %>%
  group_by(Time_sample, anno) %>%
  summarise(cell_count = n()) %>%
  pivot_wider(data = ., id_cols = Time_sample, names_from = anno, values_from = cell_count, values_fill = 0 )

cell_count_mat = data.frame(cell_count_df[,-1], check.names = F)
rownames(cell_count_mat) = cell_count_df$Time_sample
tmp = table(meta_df_list$CD45_neg$Time_sample) %>% data.frame(check.names = F)
rownames(tmp) = tmp$Var1
cell_count_mat = cell_count_mat / tmp$Freq
cell_count_mat_cd45neg = cell_count_mat

cell_count_mat_all = cbind(cell_count_mat_cd45pos, cell_count_mat_cd45neg)
cell_count_mat_all.kp = cell_count_mat_all
### remove CD45+only sample
cor_mat = cor(cell_count_mat_all) #[, anno_df$gep_name])

#cor_mat = read.csv("./007_KP_res_corMatrix_cell_module.csv", row.names = 1, check.names = F)
write.csv(cor_mat, "./007_KP_res_corMatrix_cell_module.csv")
write.csv(cell_count_mat_all,"./007_KP_res_cell_prop_matrix_all_Cd45_1o1.csv")

#rownames(cor_mat) = anno_df$anno_name
p=Heatmap(cor_mat, #row_km = 5, column_km = 5, 
          cluster_rows = T, cluster_columns=T, 
          clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
          row_names_gp = gpar(fontsize=10), column_names_gp = gpar(fontsize=10) )
od =  row_order(p)
p

pdf(file="./KP_cellmodule_heatmap.pdf", width=10, height=9)
print(p)
dev.off()

####### compare 
cell_count_mat_all = cell_count_mat_all.kp
mean_prop.met = colMeans(cell_count_mat_all[5:13,])
mean_prop.diff = log2(mean_prop.met) - log2(colMeans(cell_count_mat_all[1:4,]))

df = data.frame(mean_prop.diff)
df$name = factor(rownames(df), levels = rownames(df)[order(df$mean_prop.diff)] )
df$mean_prop.met = mean_prop.met
p = ggplot(df) +
  geom_point(aes(x=mean_prop.diff,y=name, size=mean_prop.met, color=mean_prop.diff)) +
  scale_color_distiller(type="seq", palette = "RdBu",direction = -1 ) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype="dashed",  alpha=0.3) +
  xlab(label="Proportion log2(Primary/Naive)")
p
ggsave(p, filename = "./KP_celltype_enrichment_plot.pdf", width = 4, height = 6)


### calculate p-value:
pvalue.vec = lapply(1:ncol(cell_count_mat_all), FUN = function(celltype){
  res = t.test(x = cell_count_mat_all[1:4, celltype],
               y=cell_count_mat_all[5:13, celltype] )
  return(res$p.value)
})
names(pvalue.vec) = colnames(cell_count_mat_all)
pvalue.vec = unlist(pvalue.vec)
pvalue.vec[pvalue.vec<0.05]

tmp = cor_mat#[od[[3]], od[[3]]]
tmp[tmp==1] =0
tmp[tmp >= 0.5] =1
tmp[tmp< 0.5]=0
adj = tmp

g = graph_from_adjacency_matrix(adj, mode="undirected" , weighted = T)
plot.igraph(g)


vertex_df = data.frame(node = rownames(cor_mat) ) 
rownames(vertex_df) = vertex_df$node


col_vec = c(pal_d3(palette = "category20")(20), "#DE9ED6FF")
names(col_vec) = unique( vertex_df$node)
vertex_df$color = col_vec[vertex_df$node ]

### cellprop:
mean_cell_prop = colMeans(cell_count_mat_all)
vertex_df$cell_prop = mean_cell_prop[vertex_df$node]
vertex_df$cell_prop_plot = ifelse(vertex_df$cell_prop>=0.1, 0.1, vertex_df$cell_prop)

set.seed(1001)
g3 = subgraph(g, vertex_df$node )
vertex_df = vertex_df[names(V(g3)),]
pdf(file = "./KP_cell_module_network.pdf", width = 8, height = 8)
p=plot(g, layout=layout_with_fr, edge.curved = 0.2,
       vertex.color =  vertex_df$color,
       vertex.size= 40* sqrt(vertex_df$cell_prop_plot), 
       edge.width= E(g3)$weight *2,
       edge.color = "grey",
       label.cex = 6,
       
)
p
print(p)
dev.off()








