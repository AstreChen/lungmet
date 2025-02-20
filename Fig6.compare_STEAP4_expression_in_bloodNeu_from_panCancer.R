library(tidyverse)
library(Seurat)

### merge all blood neutrophils from datasets list in Supplementary Table S6
all_blood <- readRDS('/raid1/aichen/projects/scLung/analysis/Evaluate_human/all_blood.rds')

meta_tb = unique(all_blood@meta.data[, c("patient","disease","Cancer.type","Reference","Met","platform")]) %>% unique()
head(meta_tb)
meta_tb = meta_tb %>% filter(disease %in% c("cancer","Metastasis","healthy")) %>%
  filter(platform=="10x")

obj= subset(all_blood, patient %in% meta_tb$patient)
all_blood=obj

gene_expr <- FetchData(all_blood, vars = "STEAP4")
data_combined <- cbind(all_blood@meta.data, gene_expr)
result <- aggregate(STEAP4 ~ cancer+patient, data = data_combined, FUN = mean)

result$prop <- NA
cancer_type <- unique(all_blood$cancer)

for (i in 1:length(result$cancer)){
  subset_data <- subset(data_combined,cancer==cancer_type[i])
  for (j in 1:length(subset_data$patient)){
    subset_patient_data <- subset(subset_data,patient==unique(subset_data$patient)[j])
    prop <- dim(subset_patient_data[subset_patient_data$STEAP4>0, ])[1]/dim(subset_patient_data)[1]
    result$prop[result$patient == unique(subset_data$patient)[j]] <- prop
  }
}

colnames(result) <- c('cancer','patient','gene_expr','prop')
result$cancer <- sub("COAD", "CRC", result$cancer)
result <- result %>% filter(!grepl("^BN-", patient))  ## remove BD 

result1 <- result

counts <- table(result$cancer)
valid_types <- names(counts[counts >= 4])
result <- result[result$cancer %in% valid_types, ]

sorted_data <- result %>%
  dplyr::count(cancer) %>%        
  arrange(desc(n)) %>%           
  left_join(result, by = "cancer") %>% 
  dplyr::select(-n)

sorted_data <- sorted_data %>%
  mutate(
    cancer = factor(cancer,levels = unique(cancer)),  
    patient = factor(patient,levels = unique(patient))  
  )
sorted_data$gene_expr <- as.numeric(sorted_data$gene_expr)
sorted_data$prop <- as.numeric(sorted_data$prop)

order_cc_type = sorted_data %>% filter(cancer!="healthy") %>% dplyr::group_by(cancer) %>% dplyr::summarise(mean_prop = mean(prop)) %>%
  arrange(mean_prop) %>% pull(cancer) %>% as.character()

sorted_data$cancer = factor(sorted_data$cancer,levels = c("healthy", order_cc_type))
hd_line = sorted_data %>% filter(cancer=="healthy") %>% pull(prop) %>% median() *100
p1=ggplot(sorted_data, aes(x=cancer,   
           y=100*prop, group=cancer)) +
  geom_boxplot( aes(color=cancer)) +
  geom_jitter(aes(color=cancer)) +
  scale_color_manual(values = c("grey3","#7288d0","#7d6fbc","#009d55","#e14e55","#e1ac4e","#944b00")) +
  cowplot::theme_cowplot() +
  labs(title = "Boxplot of STEAP4+Neu% by Cancer Type", x = "Cancer Type", y = "STEAP4+ Neutrophils%")+
  geom_hline(yintercept = hd_line, linetype="dashed", alpha=0.5)
my_comparisons <- list(c("healthy",'LC'),c("healthy",'GBC'),c("healthy",'PAAD'),  c("healthy",'CRC'))
pp1 <- p1 + stat_compare_means(comparisons = my_comparisons,size = 3, method = "t.test", label="p.signif")
ggsave(pp1, filename = "./plot00.cancer_type_NeuProp.pdf", width = 5, height = 4)


result <- result1
TNM_info <- read_excel('/raid1/gqh/neutrophil/neutrophil_TNM.xlsx')
colnames(TNM_info) <- c('patient','Met','Tx','TNM')
result_TNM <- merge(result, TNM_info, by = "patient", all = TRUE)
result_TNM[is.na(result_TNM)] <- 'healthy'

result_TNM = result_TNM %>%
  mutate(Met2 = ifelse(Met %in% c("LN","Yes"), "Met", Met )) %>%
  filter(Met !="Undefined") 



p1 = ggplot(result_TNM, aes(x=factor(Met2, levels = c('healthy','No','Met')), y=gene_expr, color=Met2, outlier.shape=NA   )) +
  geom_boxplot(width=0.4, outlier.shape = NA) +
  geom_jitter(aes(color=Met2),width = 0.25,size=2) +
  cowplot::theme_cowplot() +
  scale_color_manual(values = c("grey3","orangered","steelblue")) +
  xlab(label="H=18, No=21, Met=30")
my_comparisons <- list(c("healthy",'No'),c("healthy",'Met'),c("No",'Met'))
p=p1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test") +
  ggtitle("STEAP4 gene expression in human blood neutrophils") +
  ylab("Gene expression")
ggsave(p, filename = "./plot01.panCancer_blood_neutrophil.pdf", width = 4, height = 4)


patient_tb = unique(all_blood@meta.data[,c("patient","dataset","platform")] )
result_TNM = left_join(result_TNM, patient_tb, by="patient")
result_TNM$cancer.y = NULL
result_TNM$treatment= NULL
result_TNM$site="blood"
result_TNM$Met = NULL
write.csv(result_TNM,file="blood_neutrophil_STEAP4_prportion.csv")



p2 = ggplot(result_TNM, aes(x=factor(Met2, levels = c('healthy','No','Met')), y=prop*100, color=Met2, outlier.shape=NA   )) +
  geom_boxplot(width=0.4, outlier.shape = NA) +
  geom_jitter(aes(color=Met2),width = 0.25,size=2) +
  scale_color_manual(values = c("grey3","orangered","steelblue")) +
  xlab(label="H=18, No=21, Met=30") +
  scale_y_continuous(breaks = seq(0, 100, by = 20))+  
  cowplot::theme_cowplot()
my_comparisons <- list(c("healthy",'No'),c("healthy",'Met'),c("No",'Met'))
p=p2 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test", label = "p.signif") +
  ggtitle("STEAP4+(%) in human blood neutrophils") +
  ylab("STEAP4+Neu % in human blood neutrophils")
p
ggsave(p, filename = "./plot01.panCancer_blood_neutrophil_proportion.pdf", width = 4, height = 4)


