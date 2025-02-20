library(tidyverse)
library(Seurat)

### N4.Steap4 marker gene ###
gene_marker = read.csv("/raid1/aichen/projects/scLung/analysis/Neutro_Discuss/analysis/Neutrophil_Met_Up_Detect/N4.Steap4_signature_mouse_human.csv", row.names = 1)

### load object ###
our_neu = read_rds("/raid1/aichen/projects/scLung/analysis/Neutro_Discuss/data/0.our_KP_PyMT_Neu/06.KP_PyMT_Neutrophil_object.rds")

b16_neu = read_rds("/raid1/aichen/projects/scLung/analysis/Neutro_Discuss/data/2.our_B16_Neu/B16_Lung_Neutrophil.rds.gz")
neu_4t1 = read_rds("/raid1/aichen/projects/scLung/analysis/Neutro_Discuss/data/Lung_4T1_Ren/Lung_4T1_Ren_BrLM_clean.rds.gz")
neu_at3 = read_rds("/raid1/aichen/projects/scLung/analysis/Neutro_Discuss/data/Lung_AT3_Ren/Lung_AT3-gcsf_Ren_BrLM_clean.rds.gz")
neu_m39m = read_rds("/raid1/aichen/projects/scLung/analysis/Neutro_Discuss/data/Lung_M39M_Kaplan/Lung_M39M_Kaplan_Neutro.rds.gz")

################ 1. normal KP, and MMTV
### expression percentage ##
obj = our_neu
gene_expr <- FetchData(obj, vars = gene_marker$gene)
result <- cbind(obj@meta.data, gene_expr)
result$prop <- NA

gene = "Steap4"

total_cell = table(result$Sample)
sample_count = data.frame(total_cell)
colnames(sample_count) = c("Sample","count")
sample_count = unique(obj@meta.data[,c("Sample","Time","Type")])
rownames(sample_count) = sample_count$Sample


prop_matrix = matrix(data=NA, nrow= nrow(sample_count), ncol = length(gene_marker$gene))

for (j in 1:length(sample_count$Sample)){
  subdata = result[result$Sample==unique(sample_count$Sample)[j], ]
  prop <-   colSums(subdata[,gene_marker$gene] > 0) / nrow(subdata)
  prop_matrix[j,] = prop
}

rownames(prop_matrix) = sample_count$Sample
colnames(prop_matrix) = gene_marker$gene

prop_df = cbind(sample_count, prop_matrix)

mean_res =aggregate( .~Type, data=prop_df[,-c(1,2)],FUN=mean) 
sd_res = aggregate( .~Type, data=prop_df[,-c(1,2)],FUN=sd)

mean_df = data.frame(t(mean_res[,-1]))
colnames(mean_df) = mean_res$Type

mean_df$met_vs_norm = mean_df$MMTV/mean_df$Baseline
mean_df$met_vs_pri = mean_df$MMTV/mean_df$KP
#mean_df$gene = rownames(mean_df)

mean_df %>% filter(MMTV > 0.1) %>% 
ggplot(., aes(x=met_vs_norm,y=met_vs_pri,label=gene)) +
  geom_point(aes(size=MMTV))+
  ggrepel::geom_text_repel(max.overlaps = 100)

p1 = ggplot(mean_df, aes(x=MMTV,y=Baseline,label=gene)) +
  geom_point()+
  ggrepel::geom_text_repel(max.overlaps = 100) +
  

p2 = ggplot(mean_df, aes(x=MMTV,y=KP,label=gene)) +
  geom_point()+
  ggrepel::geom_text_repel(max.overlaps = 100)


### calculate p-values 
genes = gene_marker$gene

p.vsKP = numeric(length (genes))
p.vsNC = numeric(length(genes))
groups = sample_count$Type

prop_matrix
# Perform Wilcoxon tests for each gene
for (i in seq_along(genes)) {
  # Get group a data
  a_data <- prop_matrix[ groups == "MMTV", genes[i]]
  
  # a vs b comparison
  b_data <- prop_matrix[groups == "KP", genes[i]]
  if (length(a_data) >= 2 && length(b_data) >= 2) {
    p.vsKP[i] <- wilcox.test(a_data, b_data)$p.value
  } else {
    p.vsKP[i] <- NA
  }
  
  # a vs c comparison
  c_data <- prop_matrix[groups == "Baseline", genes[i]]
  if (length(a_data) >= 2 && length(c_data) >= 2) {
    p.vsNC[i] <- wilcox.test(a_data, c_data)$p.value
  } else {
    p.vsNC[i] <- NA
  }
}


result <- data.frame(
  Gene = genes,
  MMTV_vs_KP_pvalue = p.vsKP,
  MMTV_vs_NC_pvalue = p.vsNC,
  stringsAsFactors = FALSE
)

# Optional: Add multiple testing correction
result$MMTV_vs_KP_pvalue_adj <- p.adjust(result$MMTV_vs_KP_pvalue, method = "BH")
result$MMTV_vs_NC_pvalue <- p.adjust(result$MMTV_vs_NC_pvalue, method = "BH")

result = cbind(result, mean_df)

ggplot(result, aes(x=MMTV, y=log2(met_vs_pri), label=Gene)) +
  geom_point() +
  ggrepel::geom_text_repel(max.overlaps = 100) 

result$index_pctxlogFC.vskp = result$MMTV * log2(result$met_vs_pri)
result$index_pctxlogFC.vsnc = result$MMTV * log2(result$met_vs_norm)

source("~/scripts/000.Mylibrary/Annotation_gene_set.R")
membrane_protein = get_genes_from_GO_mouse("GO:0005886")
result$membrane_protein = ifelse(result$Gene %in% membrane_protein, "Membrane","Not Membrane" )


p1 = ggplot(result, aes(x=index_pctxlogFC.vskp, y=index_pctxlogFC.vsnc, label=Gene, color = membrane_protein)) +
  geom_point() +
  ggrepel::geom_text_repel(max.overlaps = 5) +
  xlab(" Proportion weighted log2FC(MMTV / KP) ") +
  ylab(" Proportion weighted log2FC(MMTV / Normal) ") + 
  theme_classic() +
  scale_color_manual(values = c("red3","black")) +
  geom_abline(slope= 3, intercept = 0.12, linetype="dashed", color = "grey" ) +
  geom_abline(slope= 3, intercept = -0.12, linetype="dashed", color = "grey")
p1
ggsave(p1, filename = "./01.point_plot_meanProp_weighted_Met_vs_KP_and_Normal.pdf", width = 6, height = 5)
  



