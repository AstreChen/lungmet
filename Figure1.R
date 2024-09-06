library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(ggbeeswarm)


load_TCGA_MET <- function(){

  pri_exp = read.table("/raid1/aichen/projects/scLung/public_data/03.lung_GTEX_TCGA_METAPRISM/TCGA_lung_PriT_Kallisto_log2TPM+1_genelevel.tsv", row.names = 1, header = T, check.names = F)
  pri_sample = read.table("/raid1/aichen/projects/scLung/public_data/03.lung_GTEX_TCGA_METAPRISM/TCGA_lung_PriT_sample.tsv", row.names = 1, header = T, check.names = F)
  pri_sample = pri_sample %>% filter(!is.na(Tumor_stage))
  pri_exp = pri_exp[, pri_sample$Tumor_Sample_Barcode]
  if( all(colnames(pri_exp) == pri_sample$Tumor_Sample_Barcode)){
    print("TCGA, checked OK")
  }
  
  
  met_exp = read.table("/raid1/aichen/projects/scLung/public_data/03.lung_GTEX_TCGA_METAPRISM/METAPRISM_lungMet_Kallisto_log2TPM+1_genelevel.tsv", row.names = 1, header = T , check.names = F)
  met_sample = read.table("/raid1/aichen/projects/scLung/public_data/03.lung_GTEX_TCGA_METAPRISM/METAPRISM_lungMet_sample.tsv", row.names = 1, header = T, check.names = F)
  if( all(colnames(met_exp) == met_sample$Sample_Id)){
    print("META-PRISM, checked OK")
  }

  gene_id = intersect( rownames(pri_exp), rownames(met_exp))
  all_data = cbind( pri_exp[gene_id, ], met_exp[gene_id,] )
  sample_df = data.frame(sample_name = colnames(all_data))
  
  tumor_stage = c(
    paste0("1.Pri_", pri_sample$Tumor_stage), 
    rep("2.Met", ncol(met_exp)) 
  )
  
  cancer = c( paste0("1.Pri_", pri_sample$Project ), 
              paste0("2.Met_", met_sample$Cancer_Cohort) )
  
  age = c(pri_sample$Age, floor(met_sample$Age_At_Biopsy))
  sex = c(pri_sample$Gender, met_sample$Sex)
  
  sample_df$tumor_stage = tumor_stage
  sample_df$cancer = cancer
  sample_df$age = age
  sample_df$sex = sex
  
  
  rownames(sample_df) = sample_df$sample_name
  cohort_lung = list(exp = all_data, sample = sample_df)
  return( cohort_lung )
  
}


lung_cohort = load_TCGA_MET()

#####   TME signature
geneset = readxl::read_xlsx(path ="/raid1/aichen/projects/scLung/public_data/02.TCGA_METAPRISM/merged_gene_signatures.xlsx")
geneset_list = split(geneset$Gene, f = geneset$`Gene signature`)

signames = c( "Checkpoint molecules",  "Co-stimulatory ligands" , "Co-stimulatory receptors",
              "CD8 Effector cell traffic" , "CD8 Effector cells",
               "Th1 signature"    , "Th2 signature"  , "Exhausted CD8"   )
geneset_list = geneset_list[signames]

# subsample primary data ###
### prepare data for model ####
zscore_func = function(x){
  (x - mean(x)) / sd(x)
}



###################################################################################
#######  use mean gene expression to calculate signature score              #######
###################################################################################

tmp = lung_cohort$sample %>%
  group_by(cancer) %>%
  summarise(count=n()) %>%
  filter(count>=3)


group_df = lung_cohort$sample %>%
  #filter(cancer %in%  tmp$cancer)
  filter(cancer %in%  c("1.Pri_LUAD","1.Pri_LUSC", paste0("2.Met_", c("ACC","BRCA","COAD","BLCA","SARC","PRAD") )) )
group_df$group = group_df$cancer
exp_data = lung_cohort$exp[, group_df$sample_name]


###  signature x sample, zscored
geneset_list$`T cells` = NULL
geneset_list$`Treg and Th2 traffic` = NULL
gene_sig_mat = lapply(geneset_list, FUN=function(gset){
  mat = lung_cohort$exp[ gset, ]
  return(colMeans(mat, na.rm = T))
}) %>% do.call(rbind,.)

zscore_gene_sig_mat = apply(gene_sig_mat, 1, FUN= zscore_func)

zscore_mat = zscore_gene_sig_mat[rownames(group_df), ]
sample_group = split(rownames(group_df), group_df$group)


mean_df = lapply(sample_group, FUN = function(sample){
  mat = zscore_mat[ sample, ]
  mean_mat = colMeans(mat)
}) %>% do.call(rbind, . ) %>% t() ### signature x cohorts



library(ComplexHeatmap)
library(circlize)
col_used =  colorRamp2(c(-0.3, 0, 0.3), c("#1F77B4FF", "white", "#D62728FF"))
p = Heatmap( t(mean_df), col = col_used,
             cluster_columns = F, cluster_rows = F,
             show_row_names = T, row_names_side = "left",
             name = "mean zscore"
)
p


pdf(file="analysis_bulk/fig1_heatmap_Lung_tumor_compare_expression_score.pdf", width = 4.5, height = 4)
print(p)
dev.off()


######## boxplot of expression ###########
tmp_group = lung_cohort$sample
tmp_group = tmp_group %>%
  mutate(type = ifelse(str_detect(cancer, "1.Pri"), "1.Pri", "2.Met")) 


for (sig in rownames(gene_sig_mat)){
  plot_df = data.frame(sigscore = t(gene_sig_mat)[, sig], sample_name=colnames(gene_sig_mat)  ) %>%
    inner_join(., y= tmp_group, by="sample_name")
  p = plot_df %>%
    ggplot(aes(x=type, y=sigscore, color=type))+
    geom_boxplot(outlier.shape = NA) +
    geom_quasirandom(size=0.5, alpha=1, shape=16, width=0.1)+
    ggpubr::stat_compare_means(,method = "t.test") +
    theme_classic() +
    scale_color_manual(values = c("steelblue","orangered")) +
    ylim(c(-0.1,NA))
  ggsave(plot = p, filename = sprintf("/raid1/aichen/projects/scLung/analysis/Lib009/analysis_bulk/gene_expression_norm_score/%s.pdf", sig),
         width = 3.5, height = 4)
}


########## boxplot of PD1, TIGIT #######
##########
plot_gene_feature <- function(gene_name){
  # gene_name = "CTLA4"
  plot_df = cbind(group_df, t(exp_data)) 
  feature_df = plot_df[,c("group", gene_name)]
  colnames(feature_df)[2] = "gene"
  feature_df = feature_df %>%
    mutate(sample_group = str_replace(group, "1.Pri_", "P.") ) %>%
    mutate(sample_group = str_replace(sample_group, "2.Met_", "M.") ) 
  
  ### calculate p-value
  expr_value = split(feature_df$gene, feature_df$sample_group)
  
  pval = paste0("M.", c("ACC","BRCA","COAD","BLCA","SARC","PRAD") ) %>%
    lapply(., function(x){
      res = t.test(x=expr_value[["P.LUAD"]], y=expr_value[[x]])
      return( round(res$p.value,3) )
    })
  names(pval) =  paste0("M.", c("ACC","BRCA","COAD","BLCA","SARC","PRAD") )
  pval_df = data.frame(x= names(pval), y= as.numeric(pval))
  
  order_grp = feature_df %>%
    group_by(sample_group) %>%
    summarise(med = median(gene)) %>%
    arrange(med) %>%
    filter(sample_group !="P.LUSC") %>%
    filter(sample_group !="P.LUAD") 
  
  feature_df$sample_group = factor(feature_df$sample_group, levels = c("P.LUSC","P.LUAD", order_grp$sample_group ))
  
  my_comparisons=list()
  for(grp in order_grp$sample_group){
    my_comparisons[[grp]] = c("P.LUAD", grp)
  }
  
  # p = ggplot(feature_df) +
  #   geom_boxplot(aes(x=sample_group, y=gene, color=sample_group), outlier.size = 0.1) +
  #   ggpubr::stat_compare_means(method = "t.test",comparisons = my_comparisons) +
  #   theme_classic() +
  #   geom_text(aes(x= x, label=y), y=4.5 , data = pval_df )+
  #   theme(axis.text.x = element_text(angle=45, hjust = 1,vjust = 1)) +
  #   scale_color_manual(values = c( RColorBrewer::brewer.pal(n=4,name="Blues")[c(3,4)], ## pri.LUSC, pri.LUAD
  #                                  rev(RColorBrewer::brewer.pal(n=9,name="RdPu")[2:7] ))) +
  #   ggtitle(gene_name)
  p=ggboxplot(feature_df, x="sample_group", y="gene",color="sample_group",outlier.shape = NA, size = 0.3)+
    ggpubr::stat_compare_means(method = "t.test",comparisons = my_comparisons) +
    geom_quasirandom(aes(x=sample_group, y=gene,color=sample_group), size=1, alpha=0.8, shape=16, width=0.1, data = feature_df)+
    theme_classic()+   
    theme(axis.text.x = element_text(angle=45, hjust = 1,vjust = 1)) +
    scale_color_manual(values = c( RColorBrewer::brewer.pal(n=4,name="Blues")[c(3,4)], ## pri.LUSC, pri.LUAD
                                   rev(RColorBrewer::brewer.pal(n=9,name="RdPu")[2:7] ))) +
    ggtitle(gene_name)
  p
}

to_plot_genes = c("PDCD1","CTLA4","CD27","LAG3")


p.list = lapply( to_plot_genes, plot_gene_feature)
names(p.list) = to_plot_genes

for(i in 1: length(p.list)){
  name = to_plot_genes[i]
  p = p.list[[i]]
  ggsave(p, filename = sprintf("./analysis_bulk/gene_box_plot/Lung_tumor_%s.pdf", name), width = 5, height = 4)
}

