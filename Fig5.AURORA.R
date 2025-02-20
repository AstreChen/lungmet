library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggpubr)


data_dir = "/raid1/aichen/projects/scLung/analysis/Neutro_Discuss/data/human_data/GSE209998_AURORA_RNAseq/"
exp_data = read_rds(file="/raid1/aichen/projects/scLung/analysis/Neutro_Discuss/data/human_data/GSE209998_AURORA_RNAseq/DESeq_rmBatch_RAP101_AURORA_study.rds.gz")

sample_aur = read.csv( file.path(data_dir, "GSE209998_AUR_129_clinical_data.csv"), row.names = 1)
sample_rap = read.csv(file.path(data_dir,"RAP_study_GSE193103_sample_checkName.csv"), row.names = 1)
sample_df = rbind(sample_rap, sample_aur)



library(GSVA)
library(nichenetr)
## only use Steap4 neutrophils
metUp = read.csv("N4.Steap4_gene_marker.csv", row.names=1)
met_up_genes_human = nichenetr::convert_mouse_to_human_symbols(metUp$gene)
met_up_genes_human = as.character(met_up_genes_human[!is.na(met_up_genes_human)])

gsva_res = gsva(as.matrix(exp_data), gset.idx.list = list(MetNeu=met_up_genes_human))
gsva_res = t(gsva_res)

if(all(sample_df$Name == colnames(exp_data))){
  merged_data <- cbind(sample_df,gsva_res)
}



aur_sup = readxl::read_excel("/raid1/aichen/projects/scLung/analysis/Neutro_Discuss/data/human_data/GSE209998_AURORA_RNAseq/Supplementary_Table.2.xlsx", sheet = "2.AURORA study")
# aur_sup_post_prim = aur_sup %>%
#   filter(`Sample Type`=="Primary") %>%
#   filter(`Tissue primary type treatment` =="Post-treatment") %>%
#   filter(`RNA Seq FreezeSet 123`=="Sequenced")


# Obtain patients with both lung and breast samples
idx_lung = grepl(pattern = "lung|Lung", x=merged_data$Anatomic.Site)
idx_breast = grepl(pattern = "Breast", x=merged_data$Anatomic.Site)

paired_data <- merged_data[(idx_breast | idx_lung) ,  ]

paired_data$Anatomic.Site = ifelse(paired_data$Anatomic.Site == "Breast", "Breast","Lung")

paired_data <- reshape(paired_data, 
                       timevar = "Anatomic.Site", 
                       idvar = "Patient", 
                       direction = "wide", 
                       sep = "_")

# Check if there are matched samples (lung and breast) for each patient
paired_data <- paired_data[complete.cases(paired_data), ]

paired_data_tm= paired_data %>%
  filter(Sample.Type_Lung =="Metastasis") %>%
  filter(!(Patient %in% unique(aur_sup_post_prim$Patient ))  )  

melted_data <- melt(paired_data_tm, id.vars = "Patient", 
                    measure.vars = c("MetNeu_Lung", "MetNeu_Breast"), 
                    variable.name = "Tissue", value.name = "Expression")

melted_data$Tissue <- gsub("MetNeu_", "", melted_data$Tissue)

melted_data = merge(melted_data, unique(paired_data_tm[,c("Patient","Study_Breast")]), by="Patient")

p = ggplot(melted_data %>% dplyr::filter(Study_Breast=="AURORA"), aes(x = Tissue, y =Expression )) +
  geom_boxplot(aes(fill = Tissue), alpha = 0.5, width=0.7) +
  geom_point( size = 2) +
  geom_line(aes(group = Patient), color = "gray", alpha = 0.7) +
  theme_classic() +
  labs(title = "STEAP4Neu marker gene Expression in Lung vs. Breast Samples",
       x = "Tissue", y = "Expression (Batch-Corrected)") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("steelblue", "red3"))

p = p + stat_compare_means(paired = T, method = "t.test")
p

ggsave(p, filename = "./Human_AURORA_compare_preTreatmentbreat_lung_paired_samples.pdf", width = 3.5, height = 5.5)

