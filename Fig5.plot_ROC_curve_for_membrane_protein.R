library(tidyverse)
library(Seurat)


### load object ###
our_neu = read_rds("/raid1/aichen/projects/scLung/analysis/Neutro_Discuss/data/0.our_KP_PyMT_Neu/06.KP_PyMT_Neutrophil_object.rds")

b16_neu = read_rds("/raid1/aichen/projects/scLung/analysis/Neutro_Discuss/data/2.our_B16_Neu/B16_Lung_Neutrophil.rds.gz")
neu_4t1 = read_rds("/raid1/aichen/projects/scLung/analysis/Neutro_Discuss/data/Lung_4T1_Ren/Lung_4T1_Ren_BrLM_clean.rds.gz")
neu_at3 = read_rds("/raid1/aichen/projects/scLung/analysis/Neutro_Discuss/data/Lung_AT3_Ren/Lung_AT3-gcsf_Ren_BrLM_clean.rds.gz")
neu_m39m = read_rds("/raid1/aichen/projects/scLung/analysis/Neutro_Discuss/data/Lung_M39M_Kaplan/Lung_M39M_Kaplan_Neutro.rds.gz")

### calculate ROC score
library(pROC)


calculate_prop <- function(obj, vars, sample_var, group_var){
  #sample: sampleID
  #group, Met or no Met
  vars = intersect(vars, rownames(obj))
  gene_expr <- FetchData(obj, vars = vars)
  result <- cbind(obj@meta.data, gene_expr)
  result$prop <- NA
  
  sample_count = unique(obj@meta.data[,c(sample_var, group_var)])
  colnames(sample_count)[1] = "Sample"
  rownames(sample_count) = sample_count$Sample
  
  
  prop_matrix = matrix(data=NA, nrow= nrow(sample_count), ncol = length(vars))
  
  for (j in 1:length(sample_count$Sample)){
    subdata = gene_expr[result[,sample_var]==unique(sample_count$Sample)[j], ]
    prop <-   colSums(subdata > 0) / nrow(subdata)
    prop_matrix[j,] = prop
  }
  
  rownames(prop_matrix) = sample_count$Sample
  colnames(prop_matrix) = vars
  
  return(list(mat = prop_matrix, group_df = sample_count) )
}


plot_mem_genes = c("Steap4","S100a9","S100a8","Cd300ld", "Cxcr2")

genes= plot_mem_genes

total_cell = table(neu_4t1_lung@meta.data[,"hash_group"])
sample_count = data.frame(total_cell)
colnames(sample_count) = c("Sample","count")

our_normal = subset(our_neu, Type=="Baseline")
res_normal = calculate_prop(our_normal, vars =genes, sample_var ="Sample", group_var = "Type" )
res_normal$mat
res_normal$group_df$group = ifelse(res_normal$group_df$Type=="Metastatic",1,0)

res_4t1 = calculate_prop(neu_4t1, vars =genes, sample_var ="hash_group", group_var = "Type" )
res_4t1$mat
res_4t1$group_df$group = ifelse(res_4t1$group_df$Type=="Metastatic",1,0)

neu_at3$sample_tissue = paste0(neu_at3$Sample, neu_at3$tissue)
#neu_at3_lung = subset(neu_at3, tissue=="Lung")
res_at3 = calculate_prop(neu_at3, vars =genes, sample_var ="sample_tissue", group_var = "Type" )
res_at3$mat
res_at3$group_df
res_at3$group_df$group = ifelse(res_at3$group_df$Type=="Metastatic",1,0)

res = calculate_prop(b16_neu, vars =genes, sample_var ="orig.ident", group_var = "Tumor" )
res$mat
res$group_df
res$group_df$group = ifelse(res$group_df$Tumor=="Metastatic",1,0)

neu_m39m = subset(neu_m39m, type %in% c("Tumor_NoGEMy","NoTumor_NoGEMy") )
res_m39m = calculate_prop(neu_m39m, vars =genes[genes %in% rownames(neu_m39m)], sample_var ="group", group_var = "type")
res_m39m$mat
res_m39m$group_df
res_m39m$group_df$group = ifelse(str_detect(res_m39m$group_df$type,"NoTumor"), 0, 1)

### calculate ROC

mat = rbind(res_normal$mat, res_4t1$mat, res_at3$mat, res$mat, res_m39m$mat)
group = c( res_normal$group_df$group, res_4t1$group_df$group , res_at3$group_df$group , res$group_df$group,  res_m39m$group_df$group )
table(group)

roc_list = lapply(plot_mem_genes, FUN=function(x){
  roc(group , mat[,x], direction = "<")
})
names(roc_list) = plot_mem_genes
auc_values = sapply(roc_list, auc)
names(auc_values) = plot_mem_genes
auc_values = round(auc_values,digits = 5)
auc_values
roc_list = roc_list[order(-auc_values)]

order_genes = plot_mem_genes[order(-auc_values)]
p2 = ggroc(roc_list, legacy.axes = T) +
  labs(title = "ROC Curves for Metastasis vs NC (B16,AT3,4T1,M3-9-M)",
       x = "False Positive Rate",
       y = "True Positive Rate") +
  geom_abline(slope = 1, linetype="dashed", color="black") +
  theme_classic() +
  scale_color_manual( values = c("#D62728FF" ,
                                 "#1F77B4FF" , "#2CA02CFF" , "#9467BDFF",
                                 "#8C564BFF", "#E377C2FF", "#7F7F7FFF"),
                      labels = c( sprintf("%s, AUC= %0.3f", order_genes[1], auc_values[order_genes[1]]), 
                                  sprintf("%s, AUC= %0.3f", order_genes[2], auc_values[order_genes[2]]), 
                                  sprintf("%s, AUC= %0.3f", order_genes[3], auc_values[order_genes[3]]),
                                  sprintf("%s, AUC= %0.3f", order_genes[4], auc_values[order_genes[4]]),
                                  sprintf("%s, AUC= %0.3f", order_genes[5], auc_values[order_genes[5]]),
                                  sprintf("%s, AUC= %0.3f", order_genes[6], auc_values[order_genes[6]]),
                                  sprintf("%s, AUC= %0.3f", order_genes[7], auc_values[order_genes[7]]),
                                  sprintf("%s, AUC= %0.3f", order_genes[8], auc_values[order_genes[8]])
                      )) 

p2
ggsave(p2, filename = "./03.ROC_curve_validate_data.pdf", width = 6, height = 5)


#########
neu_nc_and_met = subset(our_neu, Type!="KP")

plot_mem_genes = c("Steap4","S100a9","S100a8","Cd300ld", "Cxcr2")

######## metastatic vs normal ######
res_nc= calculate_prop(neu_nc_and_met, vars =plot_mem_genes, sample_var ="Sample", group_var = "Type")
res_nc$mat
res_nc$group_df
res_nc$group_df$group = ifelse(str_detect(res_nc$group_df$Type,"MMTV"), 1, 0)
table(res_nc$group_df$group )

roc_list = lapply(plot_mem_genes, FUN=function(x){
  roc(res_nc$group_df$group , res_nc$mat[,x], direction = "<")
})
names(roc_list) = plot_mem_genes
auc_values = sapply(roc_list, auc)
names(auc_values) = plot_mem_genes
auc_values = round(auc_values,digits = 3)
roc_list = roc_list[order(-auc_values)]

order_genes = plot_mem_genes[order(-auc_values)]
p1 = ggroc(roc_list, legacy.axes = T) +
  labs(title = "ROC Curves for MMTV vs NC",
       x = "False Positive Rate",
       y = "True Positive Rate") +
  geom_abline(slope = 1, linetype="dashed", color="black") +
  theme_classic() +
  scale_color_manual( values = c("#D62728FF" ,
                                 "#1F77B4FF" , "#2CA02CFF" , "#9467BDFF",
                                 "#8C564BFF", "#E377C2FF", "#7F7F7FFF"),
                      labels = c( sprintf("%s, AUC= %0.3f", order_genes[1], auc_values[order_genes[1]]), 
                                  sprintf("%s, AUC= %0.3f", order_genes[2], auc_values[order_genes[2]]), 
                                  sprintf("%s, AUC= %0.3f", order_genes[3], auc_values[order_genes[3]]),
                                  sprintf("%s, AUC= %0.3f", order_genes[4], auc_values[order_genes[4]]),
                                  sprintf("%s, AUC= %0.3f", order_genes[5], auc_values[order_genes[5]]),
                                  sprintf("%s, AUC= %0.3f", order_genes[6], auc_values[order_genes[6]]),
                                  sprintf("%s, AUC= %0.3f", order_genes[7], auc_values[order_genes[7]]),
                                  sprintf("%s, AUC= %0.3f", order_genes[8], auc_values[order_genes[8]])
                      )) 

p1
ggsave(p1, filename = "./02.ROC_curve_MMTV_vs_NC_membrane_protein.pdf", width = 6, height = 5)




p1 = DotPlot(our_neu, rev(plot_mem_genes), group="Type", cols = "RdBu",col.min = -1, col.max = 1) +
  scale_size(range = c(1, 8), breaks = c(0, 25, 50, 75)) + coord_flip() +
  theme(legend.position = "top") + ggtitle("KP,MMTV")
p1

neu_4t1$tumor = ifelse(neu_4t1$Type =="Naive","Baseline","Metastatic")
neu_at3$tumor = ifelse(neu_at3$Type =="Naive","Baseline","Metastatic")
neu_at3_lung = subset(neu_at3, tissue=="Lung")
neu_4t1_lung = subset(neu_4t1, tissue=="Lung")

p2 = DotPlot(neu_4t1_lung, rev(plot_mem_genes), group="tumor",cols = "RdBu", col.min = -0.8, col.max = 0.8) + 
  scale_size(range = c(1, 8), breaks = c(0, 25, 50, 75)) + coord_flip() +
  theme(legend.position = "top") + ggtitle("4T1_lung")
p2
p3 = DotPlot(neu_at3_lung, rev(plot_mem_genes), group="tumor",cols = "RdBu",col.min = -0.8, col.max = 0.8)+ 
  scale_size(range = c(1, 8), breaks = c(0, 25, 50, 75)) + coord_flip() +
  theme(legend.position = "top") + ggtitle("AT3_lung")
p4 = DotPlot(b16_neu, rev(plot_mem_genes), group="Tumor",cols = "RdBu",col.min = -0.8, col.max = 0.8)+
  scale_size(range = c(1, 8), breaks = c(0, 25, 50, 75)) + coord_flip() +
  theme(legend.position = "top") + ggtitle("B16")

neu_m39m$tumor = ifelse(neu_m39m$type=="Tumor_NoGEMy","Metastatic","Baseline") 
p5 = DotPlot(neu_m39m, rev(plot_mem_genes), group="tumor",cols = "RdBu",col.min = -0.8, col.max = 0.8)+ 
  scale_size(range = c(1, 8), breaks = c(0, 25, 50, 75)) + coord_flip() +
  theme(legend.position = "top")  + ggtitle("M39M")

p=p4|p2|p3|p5

ggsave(p1, filename = "./04.dotplot.Discovery_KP_MMTV_membrane_proteins.pdf", width = 2.8, height = 4)

ggsave(p, filename = "./04.dotplot.Validate_membrane_proteins.pdf", width = 10, height = 4)





