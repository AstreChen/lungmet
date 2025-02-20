library(GSVA)
library(survival)
library(survminer)
library(ggplot2)
library(tidyverse)
library(broom)

setwd("/raid1/aichen/projects/scLung/analysis/Neutro_Discuss/analysis/15.steap4_neutrophil_signature/")

tcga_sample_type_df = read_tsv("/raid1/aichen/projects/scLung/public_data/02.TCGA_METAPRISM/TCGA/tcga_from_xena/TCGA_phenotype_denseDataOnlyDownload.tsv.gz")
tcga_surv_df = read.table("/raid1/aichen/projects/scLung/public_data/02.TCGA_METAPRISM/TCGA/Survival_SupplementalTable_S1_20171025_xena_sp", header = T,sep="\t")
rownames(tcga_surv_df) = tcga_surv_df$sample

tcga_subtype = readxl::read_excel("/raid1/aichen/projects/scLung/public_data/02.TCGA_METAPRISM/TCGA/TCGA_subtype_from_PMID29628290.xlsx")

tcga_exp_data = read_tsv("/raid1/aichen/projects/scLung/public_data/02.TCGA_METAPRISM/TCGA/tcga_from_xena/tcga_RSEM_gene_fpkm.gz")
tcga_exp_data_mat = as.matrix(tcga_exp_data[,2:ncol(tcga_exp_data)])
rownames(tcga_exp_data_mat) = tcga_exp_data$sample

##### filter primary sample 
tcga_sample_type_df = tcga_sample_type_df %>%
  filter(sample_type == "Primary Tumor"        )

used_sample = intersect( colnames(tcga_exp_data_mat), tcga_sample_type_df$sample)

tcga_exp_data_mat = tcga_exp_data_mat[,used_sample]
remove(tcga_exp_data)

get_human_metneu_sig <- function(){
  metUp = read.csv("N4.Steap4_genes.csv", row.names=1)
  met_up_genes_human = nichenetr::convert_mouse_to_human_symbols(metUp$gene)
  met_up_genes_human = as.character(met_up_genes_human[!is.na(met_up_genes_human)])
  return(met_up_genes_human)
}

### convert ensembl to gene symbol
library(biomaRt)
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
ensembl_ids <- str_split(rownames(tcga_exp_data_mat), pattern = "\\.",n=2, simplify = T)[,1]
gene_mapping <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), 
                      filters="ensembl_gene_id",
                      values=ensembl_ids, 
                      mart=mart)

gene_mapping <- gene_mapping[!duplicated(gene_mapping$ensembl_gene_id), ]

rownames(tcga_exp_data_mat) = ensembl_ids
rownames(gene_mapping)=gene_mapping$ensembl_gene_id
gene_mapping = gene_mapping[ rownames(tcga_exp_data_mat) ,]

tcga_exp_data_mat2 = as.data.frame(tcga_exp_data_mat)
tcga_exp_data_mat2$Gene <- gene_mapping$hgnc_symbol

tcga_exp_data_mat2 <- tcga_exp_data_mat2[!is.na(tcga_exp_data_mat2$Gene), ]
tcga_exp_data_mat2 = tcga_exp_data_mat2[!duplicated(tcga_exp_data_mat2$Gene), ]
tcga_exp_data_mat2 <- tcga_exp_data_mat2[tcga_exp_data_mat2$Gene!="", ]

rownames(tcga_exp_data_mat2) <- tcga_exp_data_mat2$Gene
tcga_exp_data_mat2$Gene = NULL

tcga_exp_mat = tcga_exp_data_mat2
remove(tcga_exp_data_mat2)
remove(tcga_exp_data_mat)
remove(tcga_exp_data)


###### calculate MetNeu signature score for all samples ##############
met_up_genes_human = get_human_metneu_sig()
gsva_scores <- gsva(as.matrix(tcga_exp_mat), list(MetNeu=met_up_genes_human), method="gsva", kcdf="Gaussian")
gsva_scores = gsva_scores %>% t() 

gsva_score_df = data.frame(sample=colnames(tcga_exp_mat), MetNeu=gsva_scores)
write.csv(gsva_score_df,"TCGA_all_sample_GSVA_score.csv")


######### samples for analysis #############
# 1. group samples
gsva_score_df = read.csv("TCGA_all_sample_GSVA_score.csv",row.names = 1)

tcga_surv_df$label= ifelse( grepl(pattern = "Metastas", x=tcga_surv_df$new_tumor_event_type),"WithMet.", "NoMet"  )

### filter sample
tcga_surv_withmet = tcga_surv_df %>%
  filter(grepl(pattern = "Metastas", x=tcga_surv_df$new_tumor_event_type))

tcga_surv_nomet = tcga_surv_df %>%
  filter(new_tumor_event_type=="")

tcga_flt_sample = rbind(tcga_surv_withmet, tcga_surv_nomet)
tcga_flt_sample$sigScore = gsva_score_df[tcga_flt_sample$sample,"MetNeu"]
tcga_flt_sample = tcga_flt_sample[complete.cases(tcga_flt_sample$sigScore),]   

cancer_types = sort(unique(tcga_flt_sample$cancer.type.abbreviation))


tmp=tcga_flt_sample %>%
  group_by(cancer.type.abbreviation, label) %>%
  summarise(count=n()) %>%
  pivot_wider(., id_cols = cancer.type.abbreviation, names_from = label, values_fill =0 ,values_from = count)

used_cancer_types = cancer_types[tmp$WithMet.>=10]

survival_plots <- list()

# survival time truncated time

for (cancer in used_cancer_types) {
  
  cancer_data <- tcga_flt_sample %>% filter(cancer.type.abbreviation == cancer) %>%
    filter(!is.na(PFI))
  

  cutpoint <- surv_cutpoint(cancer_data, 
                            time = "PFI.time", 
                            event = "PFI", 
                            variables = "sigScore", minprop = 0.2)
  
  print(cutpoint)
  cutpoint$cutpoint[1,1]

  cancer_data$group <- factor(ifelse(cancer_data$sigScore >=   cutpoint$cutpoint[1,1], "High", "Low"))
  
  surv_object <- Surv(cancer_data$PFI.time, cancer_data$PFI)
  
  km_fit <- survfit(surv_object ~ group, data = cancer_data)
  
  plot <- ggsurvplot(km_fit, 
                     data = cancer_data, 
                     pval = TRUE, 
                     title = paste("Survival by Signature Score in", cancer),
                     risk.table = TRUE,
                     palette = c("red", "blue"),
                     legend.title = "Signature Score",
                     legend.labs = c("High", "Low"))
  
  survival_plots[[cancer]] <- plot
}
pfi_survival_plots = survival_plots
for (cancer in names(survival_plots)) {
  pdf(file =paste0("TCGA_alone_compare_pri/TCGA_onlyMet_NoMet_PFI_bestCut_", cancer, ".pdf"),  width = 5, height = 6)
  print(survival_plots[[cancer]])
  dev.off()
}



fig_cancer = c("ACC","SKCM", "PAAD","HNSC")

cancer = "SKCM"
max_truncation_time <- 5000
cancer_data <- tcga_flt_sample %>% filter(cancer.type.abbreviation == cancer) %>%
  filter(!is.na(PFI))
dim(cancer_data)

# Introduce right truncation: truncate all events at 2000 days

cancer_data = cancer_data %>%
  mutate(surv_time_truncated = pmin(PFI.time, max_truncation_time),
         status_truncated =ifelse(PFI.time > max_truncation_time, 0, PFI)  # Truncate censored events beyond 2000 days
  )

cutpoint <- surv_cutpoint(cancer_data, 
                          time = "PFI.time", 
                          event = "PFI", 
                          variables = "sigScore", minprop = 0.2)

print(cutpoint)
cutpoint$cutpoint[1,1]
cancer_data$group <- factor(ifelse(cancer_data$sigScore >=   cutpoint$cutpoint[1,1], "High", "Low"))

surv_object <- Surv(cancer_data$surv_time_truncated, cancer_data$status_truncated)

km_fit <- survfit(surv_object ~ group, data = cancer_data)


plot <- ggsurvplot(km_fit, 
                   data = cancer_data, 
                   pval = TRUE, 
                   title = paste("Survival by Signature Score in", cancer),
                   risk.table = TRUE,
                   palette = c("red", "blue"),
                   legend.title = "Signature Score",
                   legend.labs = c("High", "Low"))

plot
pdf(file =paste0("TCGA_alone_compare_pri/01.TCGA_onlyMet_NoMet_PFI_bestCut_", cancer, ".pdf"),  width = 5, height = 6)
print(plot)
dev.off()


#######################  HR ratio #############################
library(survival)
library(ggplot2)
library(dplyr)
library(broom)


tcga_data = tcga_flt_sample
colnames(tcga_data)[3]="cancer_type"
colnames(tcga_data)[4]="age"
tcga_data = tcga_data %>%
  mutate(stage = case_when(
    ajcc_pathologic_tumor_stage == "Stage 0" ~ "0",
    ajcc_pathologic_tumor_stage %in% c("Stage I", "Stage IA", "Stage IB") ~ "I",
    ajcc_pathologic_tumor_stage %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC") ~ "II",
    ajcc_pathologic_tumor_stage %in% c("Stage IIIA", "Stage IIIB", "Stage IIIC") ~ "III",
    ajcc_pathologic_tumor_stage %in% c("Stage IV", "Stage IVA", "Stage IVB", "Stage IVC") ~ "IV",
    ajcc_pathologic_tumor_stage == "Stage X" ~ "X",
    TRUE ~ NA_character_))

hr_results <- data.frame(cancer_type = character(), HR = numeric(), lower_ci = numeric(), upper_ci = numeric(), p.value=numeric())

for (cancer in used_cancer_types) {
  print(cancer)
  cancer_data <- tcga_data %>% filter(cancer_type == cancer) 
  cancer_data = cancer_data[,c("cancer_type","age","stage","gender","sigScore","PFI.time","PFI","OS.time","OS")]  %>%
    filter(!is.na(PFI.time))
  
  surv_object <- Surv(cancer_data$PFI.time, cancer_data$PFI)
  
  if(cancer %in% c("CESC","OV","PRAD","TGCT","UCEC","UCS")){
    cox_model <- coxph(surv_object ~ sigScore  + age , data = cancer_data)
    
  }else{
    cox_model <- coxph(surv_object ~ sigScore + gender + age , data = cancer_data)
    
  }

  cox_summary <- tidy(cox_model)
  
  sigScore_result <- cox_summary %>%
    filter(term == "sigScore") %>%
    mutate(
      HR = exp(estimate),  # 计算风险比HR
      lower_ci = exp(estimate - 1.96 * std.error),  # 计算置信区间下限
      upper_ci = exp(estimate + 1.96 * std.error) ,  # 计算置信区间上限
      cancer_type = cancer,
      p.value=p.value
    )
  res = sigScore_result[, c("cancer_type", "HR","lower_ci","upper_ci","p.value")]
  
  hr_results <- rbind(hr_results, res)
}

hr_results <- hr_results %>%
  mutate(
    significant = case_when(
      p.value < 0.001 ~ "***",  # p < 0.001
      p.value < 0.01  ~ "**",   # p < 0.01
      p.value < 0.05  ~ "*",    # p < 0.05
      TRUE ~ ""                  # p >= 0.05
    )
  )
hr_results_sorted <- hr_results %>%
  arrange(HR) %>%
  mutate(hr_type = ifelse(HR >1, "Worse", "Better"),
         signif_type = ifelse(p.value<0.05, "S", "NS")) %>%
  mutate(color_type = factor(paste0(hr_type, signif_type), levels = c("BetterS", "BetterNS", "WorseNS", "WorseS" ))  )

p = ggplot(hr_results_sorted, aes(x = HR, y = reorder(cancer_type, HR), xmin = lower_ci, xmax = upper_ci)) +
  geom_point(aes(color = color_type), size = 3) + # 点的颜色根据显著性变化
  geom_errorbarh(aes(color=color_type),height = 0) +  # 添加水平误差条
  scale_x_log10(breaks = c(0.01, 0.1, 1,  10, 100, 1000), labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_color_manual(values = rev(RColorBrewer::brewer.pal(n=4,"RdBu") )) + # 红色表示显著，黑色表示不显著
  geom_text(aes(label = significant), size = 5, hjust = -5 ) +  # 标注显著性符号
  theme_classic() +
  geom_vline(xintercept = 1, linetype="dashed", color="grey",alpha=0.5) +
  theme(axis.text.y = element_text(size = 8)) 
p
ggsave(p, filename = "TCGA_alone_compare_pri/TCGA_PFI_hazard_ratio_for_STEAP4Neu_signature.pdf", width = 5, height = 4)


#### overall survival  #######
os_hr_results <- data.frame(cancer_type = character(), HR = numeric(), lower_ci = numeric(), upper_ci = numeric(), p.value=numeric())

for (cancer in used_cancer_types) {
  print(cancer)
  cancer_data <- tcga_data %>% filter(cancer_type == cancer) 
  cancer_data = cancer_data[,c("cancer_type","age","stage","gender","sigScore","PFI.time","PFI","OS.time","OS")]  %>%
    filter(!is.na(PFI.time))
  
  surv_object <- Surv(cancer_data$OS.time, cancer_data$OS)
  
  if(cancer %in% c("CESC","OV","PRAD","TGCT","UCEC","UCS")){
    cox_model <- coxph(surv_object ~ sigScore  + age , data = cancer_data)
    
  }else{
    cox_model <- coxph(surv_object ~ sigScore + gender + age , data = cancer_data)
    
  }
  
  # 使用 tidy() 获取模型系数结果
  cox_summary <- tidy(cox_model)
  
  # 提取 sigScore 相关的结果
  sigScore_result <- cox_summary %>%
    filter(term == "sigScore") %>%
    mutate(
      HR = exp(estimate),  # 计算风险比HR
      lower_ci = exp(estimate - 1.96 * std.error),  # 计算置信区间下限
      upper_ci = exp(estimate + 1.96 * std.error) ,  # 计算置信区间上限
      cancer_type = cancer,
      p.value=p.value
    )
  res = sigScore_result[, c("cancer_type", "HR","lower_ci","upper_ci","p.value")]
  
  # 显示 sigScore 相关结果
  # 将结果存入hr_results数据框
  os_hr_results <- rbind(os_hr_results, res)
}

# 绘制HR的森林图
os_hr_results <- os_hr_results %>%
  mutate(
    significant = case_when(
      p.value < 0.001 ~ "***",  # p < 0.001
      p.value < 0.01  ~ "**",   # p < 0.01
      p.value < 0.05  ~ "*",    # p < 0.05
      TRUE ~ ""                  # p >= 0.05
    )
  )

hr_results_sorted <- os_hr_results %>%
  arrange(HR) %>%
  mutate(hr_type = ifelse(HR >1, "Worse", "Better"),
         signif_type = ifelse(p.value<0.05, "S", "NS")) %>%
  mutate(color_type = paste0(hr_type, signif_type))

p = ggplot(hr_results_sorted, aes(x = HR, y = reorder(cancer_type, HR), xmin = lower_ci, xmax = upper_ci)) +
  geom_point(aes(color = color_type), size = 3) + # 点的颜色根据显著性变化
  geom_errorbarh(aes(color=color_type),height = 0) +  # 添加水平误差条
  scale_x_log10(breaks = c(0.01, 0.1, 1,  10, 100, 1000), labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_color_manual(values = rev(RColorBrewer::brewer.pal(n=4,"RdBu") )[2:4]) + # 红色表示显著，黑色表示不显著
  geom_text(aes(label = significant), size = 5, hjust = -5 ) +  # 标注显著性符号
  theme_classic() +
  geom_vline(xintercept = 1, linetype="dashed", color="grey",alpha=0.5) +
  theme(axis.text.y = element_text(size = 8)) 
p
ggsave(p, filename = "TCGA_alone_compare_pri/TCGA_OS_hazard_ratio_for_STEAP4Neu_signature.pdf", width = 5, height = 4)

