#' Functions related to RUV-III-PRPS were modified from https://github.com/RMolania/TCGA_PanCancer_UnwantedVariation/
#' 

calculate_tumor_purity <- function(exp_data, output.pre ){
  require(estimate)
  filterCommonGenes_matrix <- function(input.df, output.f, id="GeneSymbol"){
    
    merged.df <- merge(common_genes, input.df, by.x = id, by.y = "row.names")
    rownames(merged.df) <- merged.df$GeneSymbol
    merged.df <- merged.df[, -1:-ncol(common_genes)]
    print(sprintf("Merged dataset includes %d genes (%d mismatched).", 
                  nrow(merged.df), nrow(common_genes) - nrow(merged.df)))
    outputGCT(merged.df, output.f)
  }
  output.f <- paste0(output.pre,'_estimate_gene.gct')
  output.ds <- paste0(output.pre,'_estimate_score.gct')
  
  filterCommonGenes_matrix(exp_data, output.f)
  
  ###计算Score
  estimateScore(input.ds = output.f,output.ds=output.ds, platform="illumina") ## 注意platform
  scores <- read.table(output.ds,skip = 2,header = T)
  rownames(scores) <- scores[,1]
  scores <- t(scores[,3:ncol(scores)])
  TumorPurity  = cos(0.6049872018+0.0001467884 * scores[,3] )
  
  output_df = data.frame(scores, row.names = colnames(exp_data))
  output_df$TumorPurity = TumorPurity
  
  file.remove(output.f)
  file.remove(output.ds)
  return(output_df)
}

# https://github.com/RMolania/TCGA_PanCancer_UnwantedVariation/blob/master/tcgaCleaneR/R/library_size_filter.R
# Removing samples based on library size

#' @title Filter Samples based on Library Size
#'
#' @description This function is a part of the data wrangling functionality of `tcgaCleaneR`.
#' It allows user to handle the bias in \code{SummarizedExperiment} S4 class Cancer Dataset (e.g. TCGA dataset) due to
#' library size by filtering out the samples with sample size greater than the threshold. Using \code{plotLibSize},
#' user can determine the threshold.
#'
#' @param data SummarizedExperiment S4 class Dataset. E.g. TCGA Dataset.
#' @param ls_cutoff numeric: library size threshold
#'
#' @return S4 data object
#' @export
#'
#' @examples
#'
#' filterSamplesByLibSize(data = brca.data, ls_cutoff = 17.5)
#'
filterSamplesByLibSize <- function(data,ls_cutoff){
  raw.count <- round(2^lung_exp -1)
  library_size <- log2(colSums(raw.count) )
  keep.samples <- library_size > ls_cutoff
  brca.se.filtered <- data[ , keep.samples]
  return(brca.se.filtered)
}


filterLowExprGenes <- function(data,gene_count,sample_size){
  raw.count <- round(2^lung_exp -1)
  keep.high <- apply(
    raw.count,
    1,
    function(x) length(x[x>gene_count])>=sample_size
  )
  data.filtered <- data[keep.high , ]
  return(data.filtered)
}

# RUV-III - PRPS((Pseudo replicate of pseudo sample)) generation

#' @title Generate PRPS for RUV-III
#'
#' @description This function is a part of the data analysis functionality of `tcgaCleaneR`. It creates pseudo-replicates
#' of pseudo-samples (PRPS) for unwanted variations like Library size, Batches and Purity in TCGA Pan Cancer Datasets with
#' Cancer biology like Breast Cancer data (BRCA), Lung Cancer (LUAD), Colon Cancer (COAD) and Rectum Cancer (READ). In the
#' function batch refers to the source of Batch Effect variation like Time and Plate which captures variation across
#' biology while factors like Purity captures variation within biology.
#'
#' @param expr.data S4 data object: Cancer Gene expression data
#' @param sample.info S4 data object: Cancer data Sample information
#' @param librarySize character: Library Size variable in input \code{sample.info} data object.
#' @param batch character: Batch effect factors. In current package version batch can take values like 'Year', 'Plate' or both
#' @param biology character: Biology of cancer type. TCGA datasets have biology for only four Cancer types i.e.
#' Lung (LUAD), Breast (BRCA), Rectum (READ) & Colon (COAD). So the function supports only these four datasets for RUV-III
#' and PRPS analysis. Default is 'Subtypes'.
#' @param purity character: Purity variable in input data object
#' @param include.ls logical: Do we need to consider library size in creating pseudo samples
#' @param include.purity logical: Do we need to consider purity in creating pseudo samples
#' @param minSamplesPerBatchPS numeric: Minimum number of samples per batch for creating Pseudo Samples
#' @param minSamplesForPuirtyPS numeric: Minimum number of samples for creating Pseudo Samples for purity.
#' @param minSamplesForPurityPerBiology numeric: Number of samples for purity per biology for creating Pseudo Samples
#' @param minSamplesForLibrarySizePerBatch numeric: Number of samples for library size per batch for creating Pseudo Samples
#' @param minSamplesForLibrarySizePS numeric: Minimum number of samples for creating Pseudo Samples for library size
#'
#' @return A S4 list object with the Pseudo replicate for pseudo samples for different batches, library size and purity.
#' @export
#'
#' @examples
#' \dontrun{
#' createPRPS(expr.data, sample.info, librarySize = 'ls', batch=c('Year', 'Plates'), biology = 'Subtypes',
#' purity='Purity_singscore',include.ls=TRUE, include.purity=TRUE,
#' minSamplesPerBatchPS = 3, minSamplesForPuirtyPS = 3, minSamplesForPurityPerBiology = 12,
#' minSamplesForLibrarySizePerBatch = 6,minSamplesForLibrarySizePS = 3)
#' }
#' 
createPRPS <- function(expr.data, sample.info, librarySize, batch, biology, purity, include.ls, include.purity,
                       minSamplesPerBatchPS, minSamplesForPuirtyPS, minSamplesForPurityPerBiology,
                       minSamplesForLibrarySizePerBatch, minSamplesForLibrarySizePS){
  
  if(include.purity & minSamplesForPuirtyPS > minSamplesForPurityPerBiology){
    stop('error: minSamplesForPuirtyPS can not be smaller than minSamplesForPurityPerBiology')
  } else if(include.purity & minSamplesForPurityPerBiology < 2*minSamplesForPuirtyPS){
    stop('error: minSamplesForPurityPerBiology should be at least two times larger than minSamplesForPuirtyPS')
  } else if(include.purity & minSamplesForLibrarySizePS > minSamplesForLibrarySizePerBatch) {
    stop('error: minSamplesForLibrarySizePerBatch can not be smaller than minSamplesForLibrarySizePS')
  } else if(include.purity & minSamplesForLibrarySizePerBatch < 2*minSamplesForLibrarySizePS){
    stop('error: minSamplesForLibrarySizePerBatch should be at least two times larger than minSamplesForLibrarySizePS')
  }
  
  biology.batch <- c(biology, batch )
  
  ### Biology
  sample.info$biology <- apply(
    sample.info[, biology, drop = FALSE],
    1,
    paste,
    collapse = "-")
  
  ### Biology - Batch
  sample.info$biology.batch <- apply(
    sample.info[, biology.batch],
    1, paste,
    collapse = "-")
  
  ### Create PS per batch
  selected.biology.ps.batch <- unlist(lapply(
    unique(sample.info$biology),
    function(x){
      index <- sample.info$biology == x
      if(sum( table(sample.info$biology.batch[index] ) >= minSamplesPerBatchPS) > 1 ){
        x
      }
    }))
  if(length(selected.biology.ps.batch) > 0){
    message('PRPS are generated for batch effects')
  }else{
    stop('error: there are not enough samples to make pseudo-samples for batch effects, you may want to lower minSamplesPerBatchPS')
  }
  sample.info.ps.batch <- sample.info[sample.info$biology %in% selected.biology.ps.batch , ]
  expr.data.ps.batch <- expr.data[, row.names(sample.info.ps.batch)]
  
  selected.batches <- names(which(table(sample.info.ps.batch$biology.batch) >= minSamplesPerBatchPS))
  ps.batch <- sapply(
    selected.batches,
    function(x) {
      index <- sample.info.ps.batch$biology.batch == x
      Matrix::rowMeans(expr.data.ps.batch[, index])
    })
  
  if(include.ls){
    selected.batches.ls <- names(
      which(table(sample.info$biology.batch) >= minSamplesForLibrarySizePerBatch)
    )
    if(length(selected.batches.ls) > 1){
      message('PRPS are generated for library size effects')
      sample.info <- sample.info[
        with(sample.info,
             order(sample.info[, 'biology.batch'],
                   sample.info[, librarySize])), ]
      expr.data <- expr.data[, row.names(sample.info)]
      ps.ls <- lapply(
        selected.batches.ls,
        function(x){
          index <- sample.info$biology.batch == x
          ls.data <- expr.data[ , index]
          low.ls <- Matrix::rowMeans(ls.data[ , 1:minSamplesForLibrarySizePS])
          high.ls <- rowMeans(ls.data[ , c(ncol(ls.data)-(minSamplesForLibrarySizePS - 1)):ncol(ls.data) ])
          all <- cbind(low.ls, high.ls)
          colnames(all) <- rep(paste(x, 'LS', sep = '-'), 2)
          all
        })
      ps.ls <- do.call(cbind, ps.ls)
      
    }else{
      stop('error: there are not enough samples to make pseudo-samples for library size effects, you may want to lower minSamplesForLibrarySizePerBatch')
    }
  }else if (!include.ls){
    warning('PRPS is not considered for librray size effects')
    ps.ls = list()
  }
  
  if(include.purity ){
    selected.biology.purity <- names(
      which(table(sample.info$biology) >= minSamplesForPurityPerBiology)
    )
    if(length(selected.biology.purity) > 0){
      message('PRPS are generated for purity effects')
      sample.info <- sample.info[
        with(sample.info,
             order(sample.info[, 'biology.batch'],
                   sample.info[, purity])),]
      expr.data <- expr.data[, row.names(sample.info)]
      ps.purity <- lapply(
        selected.biology.purity,
        function(x) {
          index <- sample.info$biology == x
          purity.data <- expr.data[, index]
          low.pur <- rowMeans(purity.data[, 1:minSamplesForPuirtyPS])
          high.pur <- rowMeans(purity.data[, c(ncol(purity.data) - (minSamplesForPuirtyPS - 1)):ncol(purity.data)])
          all <- cbind(low.pur, high.pur)
          colnames(all) <- rep(paste(x, 'purity', sep = '-'), 2)
          all
        })
      ps.purity <- do.call(cbind, ps.purity)
    }else{
      stop('error: there are not enough samples to make pseudo-samples for purity variation, you may want to lower minSamplesForPurityPerBiology')
    }
  } else if (!include.purity){
    warning('PRPS is not considered for purity effects')
    ps.purity = list()
  }
  return(list(ps.purity = ps.purity, ps.ls = ps.ls, ps.batch = ps.batch))
}


### find negative control genes by ANOVA and gene correlation

computeANOVA <- function(data, sample.info, variable, n.cores, is.log=T){
  raw.count <- as.matrix(data)
  sample.info$ls = sample.info$libsize
  
  if(is.log == TRUE){
    raw.count <- raw.count
    average.exp <- rowMeans(raw.count)
  }else{
    raw.count <- log2(raw.count + 1)
    average.exp <- log2(rowMeans(raw.count))
  }
  if (variable == "cancer_type"){
    anova_test <- parallel::mclapply(
      1:nrow(raw.count),
      function(x) {
        MASS::dropterm(lm(raw.count[x , ] ~ sample.info$cancer_type), test = 'F')[c(5:6)]
      }
      , mc.cores = n.cores)
  } else
    if (variable == "dataset"){
      anova_test <- parallel::mclapply(
        1:nrow(raw.count),
        function(x) {
          MASS::dropterm(lm(raw.count[x , ] ~ sample.info$dataset), test = 'F')[c(5:6)]
        }
        , mc.cores = n.cores)
    }
  test.values <- data.frame(
    Genes = row.names(raw.count),
    FValue = round(unlist(lapply(anova_test, function(x) x$`F Value`[2])), digits = 4) ,
    PValue = unlist(lapply(anova_test, function(x) x$`Pr(F)`[2])),
    Adj.PValue = p.adjust(unlist(lapply(anova_test, function(x) x$`Pr(F)`[2])), method = 'BH'),
    Mean = round(average.exp, digits = 2)
  )
  return(test.values)
}


#  Gene Correlation analysis function

#' @title Gene Correlation analysis
#'
#' @description This function is a part of the data analysis functionality of `tcgaCleaneR`.
#' It helps to run correlation analysis between gene expression level and unwanted variation effect because of Library
#' Size or Purity. This function can help to quantify the association between an individual gene’s expression level
#' and library size or tumor purity for a given Cancer type.
#'
#' @param data S4 data object
#' @param is.log logical: Checks if the S4 data has log values. It 'False', it converts data to log scale.
#' @param type character: Variation variable to perform correlation with. type included are 'librarysize', 'purity_HTseq_counts', 'purity_HTseq_FPKM' and 'purity_HTseq_FPKM.UQ'.
#' @param cor.method a character string indicating which correlation coefficient is to be used for the test. One of "pearson", "kendall", or "spearman", can be abbreviated. Default is spearman.
#' @param n.cores The number of cores to use, i.e. at most how many child processes will be run simultaneously. Must be at least one, and parallelization requires at least two cores.
#'
#' @return A S3 data frame. The output contains the correlation test output containing pvalue, adj p-value and
#' Spearman's rank correlation coefficient. Along with the data frame output the function also returns a histogram
#' for Spearman's rank correlation coefficient for easy analysis of the test results.
#' @export
#'
#' @examples
#' \dontrun{
#' df <- computeCorr(data = brca.data,is.log = FALSE,type = "purity",cor.method = 'spearman',n.cores = 1)
#' computeCorr(data = brca.data,is.log = FALSE,type = "librarysize",cor.method = 'pearson',n.cores = 1)
#' }

computeCorr <- function(data,sample.info, type, cor.method, n.cores, is.log=T){
  
  raw.count <- as.matrix(data)
  sample.info$ls = sample.info$libsize
  
  .corr.gene.variable <- function(expr.data, is.log, variable, method, n.cores){
    if(is.log){
      expr.data <- expr.data
    }else{
      expr.data <- log2(expr.data + 1)
    }
    rho <- parallel::mclapply(
      1:nrow(expr.data),
      function(x){
        round(cor.test(x = expr.data[x, ], y = variable, method = method)[[4]], 6)},
      mc.cores = n.cores)
    
    pval <- parallel::mclapply(
      1:nrow(expr.data),
      function(x){
        cor.test(x = expr.data[x, ], y = variable, method = method)[[3]]},
      mc.cores = n.cores)
    
    results <- data.frame(
      genes = row.names(expr.data),
      rho = unlist(rho),
      pvalue = unlist(pval),
      adj.pvalue = p.adjust(unlist(pval), 'BH')
    )
  }
  if(type == "librarysize"){
    genes.df <- .corr.gene.variable(
      expr.data = raw.count,
      is.log = is.log,
      variable = sample.info$ls,
      method = cor.method,
      n.cores = 1
    )
    hist(genes.df$rho)
  } else
    if(type == "purity"){
      genes.df <- .corr.gene.variable(
        expr.data = raw.count,
        is.log = is.log,
        variable = sample.info$tumor_purity,
        method = cor.method,
        n.cores = 1
      )
      hist(genes.df$rho)
    }
  return(genes.df)
}

# RUV-III generation

#' @title Generate RUV-III Data Object
#'
#' @description This function is a part of the data analysis functionality of `tcgaCleaneR`. It captures both the uses PRPS values from library \code{tcgaCleaneR} combined with row counts and run SVD algorithm (\code{runSVD()}) from \code{BiocSingular} on the combined dataset. The function uses RUV-I algorithm from \code{ruv} as a pre-processing step to RUV-III.
#'
#' @param ruv.data S4 data object for RUV-III: A S4 data object with combined data including the row count from original filtered data using assay \code{HTseq_counts}, prps data for batch and library size. This data needs to be further converted to log scale and transposed.
#' @param ruv.rep S4 data matrix for RUV-III: A S4 data object that has been generated using \code{replicate.matrix} functionality from \code{ruv} package. This helps ruv to identify replicate samples.
#' @param ncg.set logical object: Set of Negative Controlled genes.
#' @param k Integer scalar specifying the number of unwanted factors to use. Default is NULL. Currently Value 1 represents the library size, 2 represents purity and 3 is time variation.
#' @param eta Gene-wise (as opposed to sample-wise) covariates. A matrix with n columns for \code{ruv::RUV1}. Default is NULL.
#' @param include.intercept Add an intercept term to eta if it does not include one already for \code{ruv::RUV1}. Default is True.
#' @param average Default is False.
#' @param fullalpha To perform RUV-III calculation. Default is NULL.
#' @param return.info logical: Do you want all the information related to RUV-III object. False gives all information whereas True gives only the
#' @param inputcheck logical: Check the inputs to identify if ruv.data contains missing values or infinite values.
#'
#' @return Based on the return.info we get either a S4 list will all information related to RUV-III object or just the RUV-III result.
#' @export
#'
#' @examples
#' \dontrun{
#' runRUV_III_PRPS(ruv.data = ruv.data, ruv.rep = ruv.rep, ncg.set = ncg.set, k=1, return.info = TRUE)
#' }

runRUV_III_PRPS <- function(ruv.data, ruv.rep, ncg.set, k = NULL, eta = NULL,
                            include.intercept = TRUE, average = FALSE,
                            fullalpha = NULL, return.info = FALSE, inputcheck = TRUE){
  
  if (is.data.frame(ruv.data) ) {
    ruv.data <- data.matrix(ruv.data)
  }
  Y <- ruv.data
  M <- ruv.rep
  
  m <- nrow(ruv.data)
  n <- ncol(ruv.data)
  M <- ruv::replicate.matrix(ruv.rep)
  
  tological <- function(ctl, n) {
    ctl2 <- rep(FALSE, n)
    ctl2[ctl] <- TRUE
    return(ctl2)
  }
  
  ctl <- tological(ncg.set, n)
  if (inputcheck) {
    if (n > m)
      warning(
        "ncol for ruv.data is greater than nrow!  This is not a problem itself, but may
              indicate that you need to transpose your data matrix.
              Please ensure that rows correspond to observations
              (e.g. RNA-Seq assay) and columns correspond to features (e.g. genes).")
    if (sum(is.na(ruv.data)) > 0)
      stop("ruv.data contains missing values.  This is not supported.")
    if (sum(ruv.data == Inf, na.rm = TRUE) + sum(ruv.data == -Inf, na.rm = TRUE) >
        0)
      stop("ruv.data contains infinities.  This is not supported.")
  }
  
  Y <- ruv::RUV1(Y, eta, ctl, include.intercept = include.intercept)
  mu <- colMeans(Y)
  mu_mat <- rep(1, m) %*% t(mu)
  Y_stand <- Y - mu_mat
  if (ncol(M) >= m)
    newY <- Y
  else if (is.null(k)) {
    ycyctinv <- solve(Y[, ctl] %*% t(Y[, ctl]))
    newY <- (M %*% solve(t(M) %*% ycyctinv %*% M) %*% (t(M) %*% ycyctinv)) %*% Y
    fullalpha <- NULL
  } else if (k == 0) {
    newY <- Y
    fullalpha <- NULL
  } else {
    if (is.null(fullalpha) ) {
      Y0 <- ruv::residop(Y, M)
      fullalpha <- t(svd(Y0 %*% t(Y0))$u[, 1:min(m - ncol(M), sum(ctl)), drop = FALSE]) %*% Y
    }
    alpha <- fullalpha[1:min(k, nrow(fullalpha)), , drop = FALSE]
    ac <- alpha[, ctl, drop = FALSE]
    W <- Y_stand[, ctl] %*% t(ac) %*% solve(ac %*% t(ac))
    newY <- Y - W %*% alpha
  }
  if (average)
    newY <- ((1/apply(M, 2, sum)) * t(M)) %*% newY
  if (!return.info) {
    return(newY)
  } else {
    return(list(new.ruv.data = newY, ruv.rep = M, fullalpha = fullalpha,  W =  W))
  }
}


.pca <- function(data, nPcs, is.log=T) {
  if(is.log){
    data <- data
  }else{
    data <- log2(data + 1)
  }
  svd <- BiocSingular::runSVD(
    x = t(data),
    k = nPcs,
    BSPARAM = BiocSingular::bsparam(),
    center = TRUE,
    scale = FALSE
  )
  percent <- svd$d^2/sum(svd$d^2)*100
  percent <-
    sapply(
      seq_along(percent),
      function(i) {round(percent[i], 1)})
  return(list(
    sing.val = svd,
    variation = percent))
}



# Plot PCA function

#' @title PCA Visualization
#'
#' @description This function is a part of the data analysis functionality of `tcgaCleaneR`. It helps to visualize the PCs from \code{get.pca}.
#'
#' @param pca.data list: PCA output from \code{computePCA}.
#' @param data S4 data object
#' @param group character: Color code PCs based on group. groups included are 'Time', 'Tissue', 'Plate', 'TSS', 'Center'
#' @param plot_type character: Plot type
#' @param pcs.no numeric vector: PCs that needs to be plotted
#'
#' @return Density Plot, Box Plot
#' @export
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 geom_density
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom cowplot axis_canvas
#' @importFrom cowplot insert_xaxis_grob
#' @importFrom cowplot insert_yaxis_grob
#'
#' @examples
#' \dontrun{
#' plotPC(pca.data, data = brca.data, group = "Time", plot_type = "DensityPlot", pcs.no = c(1,2,3))
#' plotPC(pca.data = df6, data = brca.data, group = "Plate", plot_type = "BoxPlot", pcs.no = c(1,2,4))
#' }
.scatter.density.pc <- function(
  pcs,
  pc.var,
  pcs.no,
  group.name,
  group,
  color,
  strokeSize,
  pointSize,
  strokeColor,
  alpha){
  colnames(pcs) <- paste0('pc', pcs.no)
  pair.pcs <- utils::combn(ncol(pcs), 2)
  pcs.var <- utils::combn(pcs.no, 2)
  pList <- list()
  for(i in 1:ncol(pair.pcs)){
    if(i == 1){
      x <- pair.pcs[1,i]
      y <- pair.pcs[2,i]
      a <- pcs.var[1,i]
      b <- pcs.var[2,i]
      p <- ggplot(mapping = aes(
        x = pcs[,x],
        y = pcs[,y],
        fill = group)) +
        xlab(paste0('PC', a, ' (', pc.var[a], '%)')) +
        ylab(paste0('PC', b, ' (', pc.var[b], '%)')) +
        geom_point(
          aes(fill = group),
          pch = 21,
          color = strokeColor,
          stroke = strokeSize,
          size = pointSize,
          alpha = alpha) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
        theme(
          legend.position = "right",
          panel.background = element_blank(),
          axis.line = element_line(colour = "black", size = 1.1),
          legend.background = element_blank(),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          legend.key = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14)) +
        guides(fill = guide_legend(override.aes = list(size = 4))) +
        scale_fill_manual(name = group.name, values = color)
      
      le <- ggpubr::get_legend(p)
    }else{
      x <- pair.pcs[1,i]
      y <- pair.pcs[2,i]
      a <- pcs.var[1,i]
      b <- pcs.var[2,i]
      p <- ggplot(mapping = aes(
        x = pcs[,x],
        y = pcs[,y],
        fill = group)) +
        xlab(paste0('PC', a, ' (',pc.var[a],  '%)')) +
        ylab(paste0('PC', b, ' (',pc.var[b], '%)')) +
        geom_point(
          aes(fill = group),
          pch = 21,
          color = strokeColor,
          stroke = strokeSize,
          size = pointSize,
          alpha = alpha) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
        theme(
          panel.background = element_blank(),
          axis.line = element_line(colour = "black", size = 1.1),
          legend.position = "none",
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14)) +
        scale_fill_manual(values = color, name = group.name)
    }
    p <- p + ggplot2::theme(legend.position = "none")
    xdens <- cowplot::axis_canvas(p, axis = "x")+
      ggplot2::geom_density(
        mapping = ggplot2::aes(
          x = pcs[,x],
          fill = group),
        alpha = 0.7,
        size = 0.2
      ) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::scale_fill_manual(values = color)
    
    ydens <- cowplot::axis_canvas(
      p,
      axis = "y",
      coord_flip = TRUE) +
      ggplot2::geom_density(
        mapping = ggplot2::aes(
          x = pcs[,y],
          fill = group),
        alpha = 0.7,
        size = 0.2) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::scale_fill_manual(name = group.name, values = color) +
      ggplot2::coord_flip()
    
    p1 <- insert_xaxis_grob(
      p,
      xdens,
      grid::unit(.2, "null"),
      position = "top"
    )
    p2 <- insert_yaxis_grob(
      p1,
      ydens,
      grid::unit(.2, "null"),
      position = "right"
    )
    pList[[i]] <- ggdraw(p2)
  }
  pList[[i+1]] <- le
  return(pList)
}

plotPC <- function(pca.data, data,sample.info, group, plot_type, pcs.no){

  data.set.names <- names(SummarizedExperiment::assays(data))
  
  if (group == "Time"){
    if(plot_type == "DensityPlot"){
      currentCols <- RColorBrewer::brewer.pal(8, "Dark2")[-5]
      years.colors <- currentCols[1:length(unique(sample.info$Year))]
      pca.plots.time <- lapply(
        data.set.names,#tcga.harmonized,
        function(x){
          pcs <- pca.data[[x]]
          p <- .scatter.density.pc(
            pcs = pcs$sing.val$u[,pcs.no],
            pc.var = pcs$var,
            pcs.no = pcs.no,
            group.name = 'Time',
            group = sample.info$Year,
            color = years.colors,#unique(sample.info$Year),
            strokeSize = .2,
            pointSize = 3,
            strokeColor = 'gray30',
            alpha = .6)
          do.call(
            gridExtra::grid.arrange,
            c(p,
              ncol = 4, top = x)
          )
        })
    } else if (plot_type == "BoxPlot"){
      
      for (i in data.set.names){
        to.plot.pc <- pca.data[[i]]$sing.val$u
        to.plot.pc <- as.data.frame(to.plot.pc)
        to.plot.pc <- to.plot.pc[,pcs.no]
        colnames(to.plot.pc) <- paste0('pc', pcs.no)
        to.plot.pc$variable = sample.info$Year
        to.plot.pc <- to.plot.pc %>%
          tidyr::pivot_longer(-variable, names_to = 'pcs', values_to = 'var') %>%
          data.frame(.)
        p <- ggplot(to.plot.pc, aes(x = variable, y = var)) +
          geom_boxplot()+
          facet_wrap(~pcs) +
          ylab('PC') +
          ggtitle(i) +
          theme(
            panel.background = element_blank(),
            plot.title = element_text(size = 22),
            axis.line = element_line(colour = 'black', size = 1),
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14),
            strip.text.x = element_text(size = 20),
            strip.text = element_text(size = 22))
        print(p)
      }
    }
  }
  else
    if (group == "Tissue"){
      if(plot_type == "DensityPlot"){
        pca.plots.time <- lapply(
          data.set.names,
          function(x){
            pcs <- pca.data[[x]]
            p <- .scatter.density.pc(
              pcs = pcs$sing.val$u[,pcs.no],
              pc.var = pcs$var,
              pcs.no = pcs.no,
              group.name = 'Tissue',
              group = sample.info$Tissues,
              color = c("#252525","#D9D9D9"),#c('red', 'blue'),
              strokeSize = .2,
              pointSize = 3,
              strokeColor = 'gray30',
              alpha = .6)
            do.call(
              gridExtra::grid.arrange,
              c(p,
                ncol = 4,
                top = x)
            )
          })
      } else if (plot_type == "BoxPlot"){
        
        for (i in data.set.names){
          to.plot.pc <- pca.data[[i]]$sing.val$u
          to.plot.pc <- as.data.frame(to.plot.pc)
          to.plot.pc <- to.plot.pc[,pcs.no]
          colnames(to.plot.pc) <- paste0('pc', pcs.no)
          to.plot.pc$variable = sample.info$Tissues
          to.plot.pc <- to.plot.pc %>%
            tidyr::pivot_longer(-variable, names_to = 'pcs', values_to = 'var') %>%
            data.frame(.)
          p <- ggplot(to.plot.pc, aes(x = variable, y = var)) +
            geom_boxplot()+
            facet_wrap(~pcs) +
            ylab('PC') +
            ggtitle(i) +
            theme(
              panel.background = element_blank(),
              plot.title = element_text(size = 22),
              axis.line = element_line(colour = 'black', size = 1),
              axis.title.x = element_text(size = 16),
              axis.title.y = element_text(size = 16),
              axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 14),
              strip.text.x = element_text(size = 20),
              strip.text = element_text(size = 22))
          print(p)
        }
      }
    } else
      if (group == "Plate"){
        if(plot_type == "DensityPlot"){
          pca.plots.time <- lapply(
            data.set.names,
            function(x){
              pcs <- pca.data[[x]]
              p <- .scatter.density.pc(
                pcs = pcs$sing.val$u[,pcs.no],
                pc.var = pcs$var,
                pcs.no = pcs.no,
                group.name = 'Plate',
                group = sample.info$Plates,
                color = unique(factor(sample.info$Plates)),
                strokeSize = .2,
                pointSize = 3,
                strokeColor = 'gray30',
                alpha = .6)
              do.call(
                gridExtra::grid.arrange,
                c(p,
                  ncol = 4,
                  top = x)
              )
            })
        } else if (plot_type == "BoxPlot"){
          
          for (i in data.set.names){
            to.plot.pc <- pca.data[[i]]$sing.val$u
            to.plot.pc <- as.data.frame(to.plot.pc)
            to.plot.pc <- to.plot.pc[,pcs.no]
            colnames(to.plot.pc) <- paste0('pc', pcs.no)
            to.plot.pc$variable = sample.info$Plates
            to.plot.pc <- to.plot.pc %>%
              tidyr::pivot_longer(-variable, names_to = 'pcs', values_to = 'var') %>%
              data.frame(.)
            p <- ggplot(to.plot.pc, aes(x = variable, y = var)) +
              geom_boxplot()+
              facet_wrap(~pcs) +
              ylab('PC') +
              ggtitle(i) +
              theme(
                panel.background = element_blank(),
                plot.title = element_text(size = 22),
                axis.line = element_line(colour = 'black', size = 1),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                axis.text.x = element_text(size = 12),
                axis.text.y = element_text(size = 12),
                legend.text = element_text(size = 12),
                legend.title = element_text(size = 14),
                strip.text.x = element_text(size = 20),
                strip.text = element_text(size = 22))
            print(p)
          }
        }
      } else
        if (group == "TSS"){
          if(plot_type == "DensityPlot"){
            pca.plots.time <- lapply(
              data.set.names,
              function(x){
                pcs <- pca.data[[x]]
                p <- .scatter.density.pc(
                  pcs = pcs$sing.val$u[,pcs.no],
                  pc.var = pcs$var,
                  pcs.no = pcs.no,
                  group.name = 'TSS',
                  group = sample.info$TSS,
                  color = unique(factor(sample.info$TSS)),
                  strokeSize = .2,
                  pointSize = 3,
                  strokeColor = 'gray30',
                  alpha = .6)
                do.call(
                  gridExtra::grid.arrange,
                  c(p,
                    ncol = 4,
                    top = x)
                )
              })
          } else if (plot_type == "BoxPlot"){
            
            for (i in data.set.names){
              to.plot.pc <- pca.data[[i]]$sing.val$u
              to.plot.pc <- as.data.frame(to.plot.pc)
              to.plot.pc <- to.plot.pc[,pcs.no]
              colnames(to.plot.pc) <- paste0('pc', pcs.no)
              to.plot.pc$variable = sample.info$TSS
              to.plot.pc <- to.plot.pc %>%
                tidyr::pivot_longer(-variable, names_to = 'pcs', values_to = 'var') %>%
                data.frame(.)
              p <- ggplot(to.plot.pc, aes(x = variable, y = var)) +
                geom_boxplot()+
                facet_wrap(~pcs) +
                ylab('PC') +
                ggtitle(i) +
                theme(
                  panel.background = element_blank(),
                  plot.title = element_text(size = 22),
                  axis.line = element_line(colour = 'black', size = 1),
                  axis.title.x = element_text(size = 16),
                  axis.title.y = element_text(size = 16),
                  axis.text.x = element_text(size = 12),
                  axis.text.y = element_text(size = 12),
                  legend.text = element_text(size = 12),
                  legend.title = element_text(size = 14),
                  strip.text.x = element_text(size = 20),
                  strip.text = element_text(size = 22))
              print(p)
            }
          }
        } else
          if (group == "Center"){
            if(plot_type == "DensityPlot"){
              pca.plots.time <- lapply(
                data.set.names,
                function(x){
                  pcs <- pca.data[[x]]
                  p <- .scatter.density.pc(
                    pcs = pcs$sing.val$u[,pcs.no],
                    pc.var = pcs$var,
                    pcs.no = pcs.no,
                    group.name = 'Center',
                    group = sample.info$Center,
                    color = unique(sample.info$Center),
                    strokeSize = .2,
                    pointSize = 3,
                    strokeColor = 'gray30',
                    alpha = .6)
                  do.call(
                    gridExtra::grid.arrange,
                    c(p,
                      ncol = 4,
                      top = x)
                  )
                })
            } else if (plot_type == "BoxPlot"){
              
              for (i in data.set.names){
                to.plot.pc <- pca.data[[i]]$sing.val$u
                to.plot.pc <- as.data.frame(to.plot.pc)
                to.plot.pc <- to.plot.pc[,pcs.no]
                colnames(to.plot.pc) <- paste0('pc', pcs.no)
                to.plot.pc$variable = sample.info$Center
                to.plot.pc <- to.plot.pc %>%
                  tidyr::pivot_longer(-variable, names_to = 'pcs', values_to = 'var') %>%
                  data.frame(.)
                p <- ggplot(to.plot.pc, aes(x = variable, y = var)) +
                  geom_boxplot()+
                  facet_wrap(~pcs) +
                  ylab('PC') +
                  ggtitle(i) +
                  theme(
                    panel.background = element_blank(),
                    plot.title = element_text(size = 22),
                    axis.line = element_line(colour = 'black', size = 1),
                    axis.title.x = element_text(size = 16),
                    axis.title.y = element_text(size = 16),
                    axis.text.x = element_text(size = 12),
                    axis.text.y = element_text(size = 12),
                    legend.text = element_text(size = 12),
                    legend.title = element_text(size = 14),
                    strip.text.x = element_text(size = 20),
                    strip.text = element_text(size = 22))
                print(p)
              }
            }
          }
  #return(data.set.names)
}
