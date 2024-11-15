library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(ggbeeswarm)
library(ggplot2)
library(cowplot)
library(ggpubr)


load_TCGA_MET <- function(){

  pri_exp = read.table("/raid1/aichen/projects/scLung/public_data/03.lung_TCGA_METAPRISM/TCGA_lung_PriT_Kallisto_log2TPM+1_genelevel.tsv", row.names = 1, header = T, check.names = F)
  pri_sample = read.table("/raid1/aichen/projects/scLung/public_data/03.lung_TCGA_METAPRISM/TCGA_lung_PriT_sample.tsv", row.names = 1, header = T, check.names = F)
  pri_sample = pri_sample %>% filter(!is.na(Tumor_stage))
  pri_exp = pri_exp[, pri_sample$Tumor_Sample_Barcode]
  if( all(colnames(pri_exp) == pri_sample$Tumor_Sample_Barcode)){
    print("TCGA, checked OK")
  }
  
  
  met_exp = read.table("/raid1/aichen/projects/scLung/public_data/03.lung_TCGA_METAPRISM/METAPRISM_lungMet_Kallisto_log2TPM+1_genelevel.tsv", row.names = 1, header = T , check.names = F)
  met_sample = read.table("/raid1/aichen/projects/scLung/public_data/03.lung_TCGA_METAPRISM/METAPRISM_lungMet_sample.tsv", row.names = 1, header = T, check.names = F)
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


######### use RUV-III-PRPS correct unwanted variations #########
# library size
# tumor purity
# 


### 
source("./source_ESTIMATE_RUVIIPRPS.R")



### logTPM 
lung_data = load_TCGA_MET()


lung_exp = lung_data$exp
lung_sample = lung_data$sample
lung_sample = lung_sample %>%
  mutate(dataset = ifelse(str_detect(sample_name, "TCGA"), "TCGA", "META-PRISM" ))

lung_luad_met = read.csv("/raid1/aichen/projects/scLung/public_data/02.TCGA_METAPRISM/META-PRISM/META-PRISM_LUAD_in_Lung_Log2TPM.csv", row.names = 1, check.names = F)
lung_luad_met_sample= read.csv("/raid1/aichen/projects/scLung/public_data/02.TCGA_METAPRISM/META-PRISM/META-PRISM_LUAD_in_Lung_sample.csv", row.names = 1)

if( all(rownames(lung_exp)== rownames(lung_luad_met)) ){
  lung_luad_met= lung_luad_met[ rownames(lung_exp), ]
  lung_exp = cbind(lung_exp, lung_luad_met)
  lung_sample2 = data.frame(sample_name = lung_luad_met_sample$Sample_Id, 
                            tumor_stage = "2.Met",
                            cancer = "2.Met_LUAD",
                            age = lung_luad_met_sample$Age_At_Biopsy,
                            sex = lung_luad_met_sample$Sex,
                            dataset = "META-PRISM")
  rownames(lung_sample2)= lung_sample2$sample_name
  lung_sample= rbind(lung_sample, lung_sample2)
}
### filter LUSC #######
idx = which(lung_sample$cancer!="1.Pri_LUSC")
lung_exp = lung_exp[, idx ]
lung_sample = lung_sample[idx,]
lung_sample = lung_sample %>% 
  mutate(cancer_type = ifelse(str_detect(cancer, "LUAD"), "LUAD", "LungMet"))

dim(lung_exp)
tumor.purity = calculate_tumor_purity(lung_exp, output.pre = "TCGA_METAPRISM_+LUAD_tumor_purity")

### prepare input for ruv:
# rawCounts.cancer
lung_exp %>% dim()

lung_exp_tpm = round(2^lung_exp -1)
libsize = colSums(lung_exp_tpm)

cov_df = data.frame(
  sample = lung_sample$sample_name,
  tumor_purity = tumor.purity$TumorPurity, 
  dataset = lung_sample$dataset,
  cancer_type = lung_sample$cancer_type,
  libsize = log2(libsize) 
                    )


### remove outlier #########
####  filter by library size
ls_cutoff = median(cov_df$libsize) - 3* stats::mad(cov_df$libsize)  # Median Absolute Deviation
ls_cutoff
lung_exp = filterSamplesByLibSize(lung_exp, ls_cutoff)
lung_sample = lung_sample[ colnames(lung_exp),]
cov_df = cov_df[colnames(lung_exp),]

### filter genes ###
#To do so, we kept genes that have at least 15 counts in the smallest biological subpopulations within each of the key time intervals in the datasets.
gene_count=5
sample_size= 10 
lung_exp = filterLowExprGenes(lung_exp, gene_count, sample_size)

dim(lung_exp)


####### Select negative control genes:
#' So, we select genes that: 
#' 1) are not highly affected by the cancer type, LUAD and LungMet (major biological variation), 
#' 2) are highly affected by dataset(TCGA and META-PRISM)
#' 3) are highly correlated with tumour purity and 
#' 4) are highly correlated with library size.
#'  We use different cut-offs to selected different sets of NCGs and assessed the performance. We found that the NCGs set 6 (see the R code below) is better than the other sets. Then, we will use that NCGs set for RUV-III normalizations.

gene.Anno.ncg <- data.frame(gene =rownames(lung_exp))

ftest.cancer_type <- computeANOVA(data = lung_exp, sample.info = cov_df, variable = "cancer_type", n.cores = 8, is.log = T )
gene.Anno.ncg$cancer_type = rank(ftest.cancer_type$FValue)

ftest.dataset <- computeANOVA(data = lung_exp[, cov_df$cancer_type=="LUAD" ], sample.info = cov_df[cov_df$cancer_type=="LUAD", ], variable = "dataset", n.cores = 8, is.log = T )
ftest.dataset = ftest.dataset %>% arrange(-FValue)

gene.Anno.ncg$dataset = rank( - ftest.dataset$FValue)



cov_df$dataset_cancer = paste0(cov_df$dataset,"_", cov_df$cancer_type)
corr.purity.per.cal <- lapply(
  unique(cov_df$dataset_cancer),
  function(x){
    idx = cov_df$dataset_cancer ==x
    computeCorr(data=lung_exp[,idx], sample.info=cov_df[idx,], type = "purity", cor.method = "spearman", n.cores = 8, is.log = T)
  }
)
names(corr.purity.per.cal) = unique(cov_df$dataset_cancer)

for(i in names(corr.purity.per.cal) ){
  gene.Anno.ncg[, i] <- abs( corr.purity.per.cal[[i]]$rho)
}

corr.ls <-  computeCorr(data=lung_exp, sample.info=cov_df, type = "librarysize", cor.method = "spearman", n.cores = 8, is.log = T)
gene.Anno.ncg[,"ls_rho"] = abs(corr.ls$rho)

purity.corr <- .5
ls.corr <- .5



nGenes <- c(
  1000,
  2000,
  3000,
  6000
)

ncg.set.dataset =  lapply(nGenes, FUN=function(x){
  gene.Anno.ncg %>%
    filter(cancer_type < x) %>%
    filter(dataset < x ) %>%
    pull(gene)
  }
)

ncg.set.purity =  lapply(nGenes, FUN=function(x){
  gene.Anno.ncg %>%
    filter(cancer_type < x) %>%
    filter(TCGA_LUAD > purity.corr) %>%
    filter(`META-PRISM_LungMet` > purity.corr) %>%
    filter(`META-PRISM_LUAD` > purity.corr) %>%
    pull(gene)
}
)
  
ncg.set.ls <- lapply(nGenes, FUN=function(x){
  gene.Anno.ncg %>%
    filter(cancer_type < x) %>%
    filter(ls_rho > ls.corr) %>%
    pull(gene)
})

all.ncg.sets <- lapply(
  c(1:length(nGenes)),
  function(x){
    ncg.set <- c(
      ncg.set.ls[[x]],
      ncg.set.purity[[x]],
      ncg.set.dataset[[x]]
    ) %>% unique()
  })

names(all.ncg.sets) = nGenes  
write_rds(all.ncg.sets, file="negative_ctrl_gene.rds.gz",compress = "gz" )



## Assessments of negative control genes

ncg.genes= all.ncg.sets$`3000`
pca.res = .pca(lung_exp[ncg.genes, ],nPcs = 10, is.log = T)

cols <- c(
  'cancer_type',
  'dataset'
)
pam50.colors.b <- pam50.colors[2:6]
name.cols <- c(
  'cancer type',
  'dataset'
)
colors <- c(
  'cancer_type',
  'dataset'
)
color.list = list(cancer_type= c("red","blue"), dataset = RColorBrewer::brewer.pal(n=3, name="Pastel1")   )

pp <- lapply(
  c(1:2),
  function(x){
    p <- .scatter.density.pc(
      pcs = pca.res$sing.val$u[, 1:3],
      pc.var = pca.res$variation,
      pcs.no = c(1,2,3),
      group.name = name.cols[x],
      group = cov_df[ , cols[x]],
      color = color.list[colors[x]][[1]],
      strokeSize = .2,
      pointSize = 2,
      strokeColor = 'gray30',
      alpha = .6
    )
    p
  })
do.call(
  gridExtra::grid.arrange,
  c(pp[[1]],
    pp[[2]],
    ncol = 4))

### Create PRPS 
# createPRPS 

lung_ps_list = createPRPS(lung_exp, sample.info = cov_df, librarySize = "libsize", batch = "dataset", biology = "cancer_type", 
                          purity = "tumor_purity", include.ls = T, include.purity = T, 
                          minSamplesPerBatchPS = 10, 
                          minSamplesForPuirtyPS = 10,
                          minSamplesForPurityPerBiology = 60,
                          minSamplesForLibrarySizePerBatch=60,
                          minSamplesForLibrarySizePS = 10
)




### RUV_III.R https://github.com/RMolania/TCGA_PanCancer_UnwantedVariation/blob/master/tcgaCleaneR/R/RUV_III.R
# runRUV_III_PRPS
ruv.data = cbind(
  lung_exp,
  lung_ps_list$ps.ls,
  lung_ps_list$ps.batch,
  lung_ps_list$ps.purity
) %>% t()

## replicate matrix
rep.matrix.ruv <- ruv::replicate.matrix(
  row.names(ruv.data)
) 

## ruv-iii normalization
ruviii.norm = runRUV_III_PRPS(ruv.data = ruv.data, ruv.rep = rep.matrix.ruv, 
                              ncg.set = colnames(ruv.data) %in% all.ncg.sets$`2000`, k = 5)


dim(ruviii.norm)

## compare gene expression


plot_gene_feature <- function(gene_name, group_df, exp_data, group_label="cancer", base_gene="CD8A", ratio=F){
  # gene_name = "CTLA4"
  group_df$group = group_df[,group_label]
  plot_df = cbind(group_df, t(exp_data)) 

  ### calculate p-value
  if(ratio){  
    feature_df = plot_df[,c("group", gene_name, base_gene)]
    colnames(feature_df)[c(2,3)] = c("gene","basegene")
  
    feature_df$ratio = log2(feature_df$gene / feature_df$basegene)

    feature_df$gene=feature_df$ratio
    feature_df = feature_df %>%
      mutate(sample_group = str_replace(group, "1.Pri_", "P.") ) %>%
      mutate(sample_group = str_replace(sample_group, "2.Met_", "M.") ) 
    expr_value = split(feature_df$ratio, feature_df$sample_group)
  }else{
    feature_df = plot_df[,c("group", gene_name)]
    colnames(feature_df)[2] = c("gene") 
    feature_df = feature_df %>%
      mutate(sample_group = str_replace(group, "1.Pri_", "P.") ) %>%
      mutate(sample_group = str_replace(sample_group, "2.Met_", "M.") ) 
    expr_value = split(feature_df$gene, feature_df$sample_group)
  }
    
  feature_df = feature_df %>%
    mutate(sample_group = str_replace(group, "1.Pri_", "P.") ) %>%
    mutate(sample_group = str_replace(sample_group, "2.Met_", "M.") ) 
  used_group = c("P.LUSC","P.LUAD", paste0("M.", c("ACC","BRCA","COAD","BLCA","SARC","PRAD","LUAD") ) )
  idx = feature_df$sample_group %in% used_group
  exp_data= exp_data[,idx]
  feature_df = feature_df[idx,]
  
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
  
  feature_df$sample_group = factor(feature_df$sample_group, levels = c("P.LUAD", order_grp$sample_group ))
  
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
    ggpubr::stat_compare_means(aes(label = after_stat(p.signif)), method = "t.test",comparisons = my_comparisons ) +
    ggbeeswarm::geom_quasirandom(aes(x=sample_group, y=gene,color=sample_group), size=0.5, alpha=0.8, shape=16, width=0.1, data = feature_df)+
    theme_classic()+   
    theme(axis.text.x = element_text(angle=45, hjust = 1,vjust = 1)) +
    scale_color_manual(values = c( RColorBrewer::brewer.pal(n=4,name="Blues")[c(4)], ## pri.LUSC, pri.LUAD
                                   rev(RColorBrewer::brewer.pal(n=9,name="RdPu")[3:8] ))) +
    ylim(c(-0.5,NA)) +
    ggtitle(gene_name)
  p
}

lung_exp_ruv = t(ruviii.norm)[,lung_sample$sample_name ]

p= plot_gene_feature(gene_name="PDCD1",group_df = lung_sample[lung_sample$cancer!="2.Met_LUAD",], 
                     exp_data = lung_exp_ruv[,lung_sample$cancer!="2.Met_LUAD"])
p
ggsave(p, filename = "PDCD1_RUVIIIPRPS_corrected_batch.pdf", width=4, height = 3)

p = plot_gene_feature(gene_name="GZMK",group_df = lung_sample[lung_sample$cancer!="2.Met_LUAD",], 
                      exp_data = lung_exp_ruv[,lung_sample$cancer!="2.Met_LUAD"])
p
ggsave(p, filename = "GZMK_RUVIIIPRPS_corrected_batch.pdf", width=4, height = 3)

p = plot_gene_feature(gene_name="GZMB",group_df = lung_sample[lung_sample$cancer!="2.Met_LUAD",], 
                      exp_data = lung_exp_ruv[,lung_sample$cancer!="2.Met_LUAD"])
p
ggsave(p, filename = "GZMB_RUVIIIPRPS_corrected_batch.pdf", width=4, height = 3)


p = plot_gene_feature(gene_name="PRF1",group_df = lung_sample[lung_sample$cancer!="2.Met_LUAD",], 
                      exp_data = lung_exp_ruv[,lung_sample$cancer!="2.Met_LUAD"])
p
ggsave(p, filename = "PRF1_RUVIIIPRPS_corrected_batch.pdf", width=4, height = 3)


p = plot_gene_feature(gene_name="GNLY",group_df = lung_sample[lung_sample$cancer!="2.Met_LUAD",], 
                      exp_data = lung_exp_ruv[,lung_sample$cancer!="2.Met_LUAD"])
p
ggsave(p, filename = "GNLY_RUVIIIPRPS_corrected_batch.pdf", width=4, height = 3)






#### plot box plot
to_plot_genes1 = c("GNLY","GZMB","GZMK","GZMA","PRF1","CXCL13")

#### plot box plot
to_plot_genes1 = c("CTLA4","TIGIT","LAG3","PDCD1")

#to_plot_genes1 = c("TCF7","IL7R","CCR7")

exp_data = lung_exp_ruv[,lung_sample$cancer!="2.Met_LUAD"]
group_df = lung_sample[lung_sample$cancer!="2.Met_LUAD",]

for (sig in to_plot_genes1){
  plot_df = data.frame(sigscore = exp_data[sig,], type = group_df[, "cancer_type"]  )
  p = plot_df %>%
    ggplot(aes(x=type, y=sigscore, color=type))+
    geom_boxplot(outlier.shape = NA) +
    ggbeeswarm::geom_quasirandom(size=0.5, alpha=1, shape=16, width=0.1)+
    ggpubr::stat_compare_means(aes(label = after_stat(p.signif)),method = "t.test", ) +
    theme_classic() +
    scale_color_manual(values = c("steelblue","orangered")) +
    ylab(label=sig)
  ggsave(plot = p, filename = sprintf("/raid1/aichen/projects/scLung/analysis/Lib009/analysis_bulk/estimate_tumor_purity_RUVIIIPRPS_rmBatch/boxplot_%s.pdf", sig),
         width = 3.5, height = 3)
}


#############   T  cell activation signature ################
geneset = readxl::read_xlsx(path ="/raid1/aichen/projects/scLung/public_data/02.TCGA_METAPRISM/TableS1.xlsx",skip = 1)
geneset_list = split(geneset$Gene, f = geneset$`Gene signature`)

library(GSVA)

### signature x samples
gsva.res = GSVA::gsva(lung_exp_ruv, geneset_list)


group_df = lung_sample %>%
  #filter(cancer %in%  tmp$cancer)
  filter(cancer %in%  c("1.Pri_LUAD", paste0("2.Met_", c("ACC","BRCA","COAD","BLCA","SARC","PRAD") )) )
group_df$group = group_df$cancer
exp_data = gsva.res[, group_df$sample_name]

gsva_mat = gsva.res[, rownames(group_df) ]  %>% t()
sample_group = split(rownames(group_df), group_df$cancer)


mean_df = lapply(sample_group, FUN = function(sample){
  mat = gsva_mat[ sample, ]
  mean_mat = colMeans(mat)
}) %>% do.call(rbind, . ) %>% t() ### signature x cohorts

library(ComplexHeatmap)
library(circlize)
col_used =  colorRamp2(c(-0.3, 0, 0.3), c("#1F77B4FF", "white", "#D62728FF"))
p = Heatmap( t(mean_df), col = col_used,
             cluster_columns = F, cluster_rows = F,
             show_row_names = T, row_names_side = "left",
             name = "GSVA score"
)
p


pdf(file="202411_fig1_heatmap_Lung_tumor_compare_expression_GSVAscore.pdf", width = 4.5, height = 4)
print(p)
dev.off()

######## boxplot of expression ###########

for (sig in rownames(gsva.res)){
  plot_df = data.frame(sigscore = t(gsva.res)[, sig], sample_name=colnames(gsva.res)  ) %>%
    inner_join(., y= lung_sample, by="sample_name") %>%
    filter(cancer != "2.Met_LUAD")
  p = plot_df %>%
    ggplot(aes(x=cancer_type, y=sigscore, color=cancer_type))+
    geom_boxplot(outlier.shape = NA) +
    #geom_jitter(aes(x=cancer_type, y=sigscore), size=0.1) +
    ggbeeswarm::geom_quasirandom(size=0.5, alpha=1, shape=16, width=0.1)+
    ggpubr::stat_compare_means(aes(label = after_stat(p.signif)),method = "t.test") +
    theme_classic() +
    scale_color_manual(values = c("steelblue","orangered")) 
  ggsave(plot = p, filename = sprintf("/raid1/aichen/projects/scLung/analysis/Lib009/analysis_bulk/estimate_tumor_purity_RUVIIIPRPS_rmBatch/SigGSVA_%s.pdf", sig),
         width = 3.5, height = 3)
}

### interferon gamma pathway; and interferon alpha pathway
library(GSVA)
library(GSEABase)
hallmarks_list = getGmt("/raid1/aichen/resources/MsigDB/human_MsigDBv2022.1/msigdb_v2022.1.Hs_files_to_download_locally/msigdb_v2022.1.Hs_GMTs/h.all.v2022.1.Hs.symbols.gmt")

hallmarks.gsva = GSVA::gsva(lung_exp_ruv, hallmarks_list )
plot_gset_name = c("HALLMARK_INTERFERON_GAMMA_RESPONSE" , "HALLMARK_INTERFERON_ALPHA_RESPONSE"   ,"HALLMARK_IL2_STAT5_SIGNALING"    )
for (sig in plot_gset_name){
  plot_df = data.frame(sigscore = t(hallmarks.gsva)[, sig], sample_name=colnames(hallmarks.gsva)  ) %>%
    inner_join(., y= lung_sample, by="sample_name") %>%
    filter(cancer != "2.Met_LUAD")
  p = plot_df %>%
    ggplot(aes(x=cancer_type, y=sigscore, color=cancer_type))+
    geom_boxplot(outlier.shape = NA) +
    #geom_jitter(aes(x=cancer_type, y=sigscore), size=0.1) +
    ggbeeswarm::geom_quasirandom(size=0.5, alpha=1, shape=16, width=0.1)+
    ggpubr::stat_compare_means(aes(label = after_stat(p.signif)),method = "t.test") +
    theme_classic() +
    scale_color_manual(values = c("steelblue","orangered")) 
  ggsave(plot = p, filename = sprintf("/raid1/aichen/projects/scLung/analysis/Lib009/analysis_bulk/estimate_tumor_purity_RUVIIIPRPS_rmBatch/SigGSVA_%s.pdf", sig),
         width = 3.5, height = 3)
}


for (sig in plot_gset_name){
  p=plot_gene_feature(gene_name = sig, group_df = lung_sample[lung_sample$cancer!="2.Met_LUAD",], 
                      exp_data = hallmarks.gsva[, lung_sample$cancer!="2.Met_LUAD"],ratio = F)
  ggsave(plot = p, filename = sprintf("/raid1/aichen/projects/scLung/analysis/Lib009/analysis_bulk/estimate_tumor_purity_RUVIIIPRPS_rmBatch/SigGSVA_%s.pdf", sig),
         width = 4, height = 3)
}



