#' Normalize RNA-seq matrices
#'
#' This function normalize single cell RNA-seq matrices with several different method
#' @param mat single cell RNA-seq matrices
#' @param method method for normalize. log, quantile, rle, scran(default), cpm, tmm, tpm, upperquantile, downsample, center
#' @param min_size if run scran, minimum size of quick clustering. default=10
#' @return normalized matrices
#' @export
#' @examples
#' NormalizeMatrix(mat, method="scran", min_size=10)
NormalizeMatrix <- function(mat, method="scran", min_size=10){
  mat <- as.matrix(mat)
  switch(method,
         log={
           mat.norm <- log2(mat + 1)
         },
         quantile={
           suppressPackageStartupMessages(library(preprocessCore))
           mat.norm <- normalize.quantiles(log2(mat + 1))
           rownames(mat.norm) <- rownames(mat)
           colnames(mat.norm) <- colnames(mat)
         },
         rle={
           suppressPackageStartupMessages(library(scater))
           suppressPackageStartupMessages(library(scran))
           sce <- SingleCellExperiment(list(counts=mat))
           rm(mat)
           mat.norm <- logcounts(normaliseExprs(
             sce,
             method = "RLE",
             return_log = TRUE,
             return_norm_as_exprs = TRUE
           ))
         },
         scran={
           suppressPackageStartupMessages(library(scater))
           suppressPackageStartupMessages(library(scran))
           sce <- SingleCellExperiment(list(counts=mat))
           rm(mat)
           clusters <- quickCluster(sce, min.size=min_size)
           print(table(clusters))
           sce <- computeSumFactors(sce, cluster=clusters, positive=TRUE)
           print(summary(sizeFactors(sce)))
           sce <- normalize(sce)
           mat.norm <- logcounts(sce)
         },
         cpm={
           suppressPackageStartupMessages(library(scater))
           suppressPackageStartupMessages(library(scran))
           sce <- SingleCellExperiment(list(counts=mat))
           rm(mat)
           mat.norm <- log2(calculateCPM(sce, use.size.factors = FALSE) + 1)
         },
         tmm={
           suppressPackageStartupMessages(library(scater))
           suppressPackageStartupMessages(library(scran))
           sce <- SingleCellExperiment(list(counts=mat))
           rm(mat)
           mat.norm <- logcounts(normaliseExprs(
             sce,
             method = "TMM",
             return_log = TRUE,
             return_norm_as_exprs = TRUE
           ))
         },
         tpm={
           suppressPackageStartupMessages(library(scater))
           suppressPackageStartupMessages(library(scran))
           sce <- SingleCellExperiment(list(counts=mat))
           rm(mat)
           mat.norm <- calculateTPM(
             sce,
             effective_length = 5e04,
             calc_from = "counts"
           )
           mat.norm <- log2(mat.norm/10 + 1)
         },
         upperquantile={
           suppressPackageStartupMessages(library(scater))
           suppressPackageStartupMessages(library(scran))
           sce <- SingleCellExperiment(list(counts=mat))
           rm(mat)
           mat.norm  <- logcounts(normaliseExprs(
             sce,
             method = "upperquartile",
             p = 0.99,
             return_log = TRUE,
             return_norm_as_exprs = TRUE
           ))
         },
         downsample={
           Down_Sample_Matrix <- function (expr_mat) {
             min_lib_size <- min(colSums(expr_mat))
             down_sample <- function(x) {
               prob <- min_lib_size/sum(x)
               return(unlist(lapply(x, function(y) {
                 rbinom(1, y, prob)
               })))
             }
             down_sampled_mat <- apply(expr_mat, 2, down_sample)
             return(down_sampled_mat)
           }
           mat.norm <- log2(Down_Sample_Matrix(mat) + 1)
         },
         center={
           # only centerizing (Tirosh 2016 et al) use normalized value
           mat.norm <- t(apply(mat, 1, scale, center=TRUE, scale=FALSE))
           colnames(mat.norm) <- colnames(mat)
         },
         {
           print('select correct normalization method')
           q()
         }
  )
  mat.norm
}
