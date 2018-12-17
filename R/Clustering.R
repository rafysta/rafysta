#' Clustering by PCA
#'
#' Perform PCA using prcomp function
#' @param DATA input matrices or seurat object. Row is genes, Column is samples
#' @param genes vector of gene list
#' @param pca_max how many principal component will be output (default=30)
#' @return pca object
#' @export
#' @examples
#' Clustering_PCA(DATA, genes=Gene_list)
Clustering_PCA <- function(DATA, genes=NULL, cells=NULL, pca_max = 30){
  Data_type = as.character(class(DATA))
  if(Data_type == "seurat"){
    mat <- DATA@scale.data
  }else if(Data_type != "matrix"){
    mat <- as.matrix(DATA)
  }else{
    mat <- DATA
  }

  if(!is.null(genes)){
    mat <- mat[genes,]
  }
  mat <- t(mat)
  pca <- prcomp(mat, scale=FALSE, tol = 0, rank. = pca_max)
  pca
}

RemoveOutlier <- function(x){
  abs(x - median(x))/mad(x, constant=1)
}


#' remove outlier from PCA result
#'
#' Define outlier from PCA
#' @param pca pca object from Clustering_PCA
#' @param file if set file name, output which cells were removed.
#' @return vector of cell names
#' @examples
#' Filter_PCA(pca, file="output.png")
Filter_PCA <- function(pca, t_outlier=20, file=NULL){
  suppressPackageStartupMessages(library(dplyr))
  D_table <- data.frame(Cell=rownames(pca$x), PC1=pca$x[,1], PC2=pca$x[,2], stringsAsFactors = FALSE)
  mad_PC1 <- RemoveOutlier(D_table$PC1) < 20
  mad_PC2 <- RemoveOutlier(D_table$PC2) < 20
  cell_name <- rownames(pca$x)
  ok_cell <- cell_name[mad_PC1 & mad_PC2]

  ### おかしなcellの名前をプロット
  if(!is.null(file)){
    suppressWarnings(suppressMessages(library(ggplot2)))
    suppressWarnings(suppressMessages(library(cowplot)))
    suppressWarnings(suppressMessages(library(ggrepel)))

    D_pca <- D_table %>% select(PC1, PC2) %>% tidyr::gather(key = "PC", value = "score") %>% mutate(mad=c(mad_PC1, mad_PC2)) %>%
      mutate(name=if_else(mad, "", rep(cell_name,2)))
    p <- ggplot(D_pca, aes(x=PC, y=score, colour=PC, fill=PC,  label =name)) + geom_violin() +
      geom_jitter(colour=if_else(D_pca %>% pull(mad), "gray20", "purple"))+
      geom_text_repel()+
      theme(legend.position="none") + labs(x="", y="PCA score")
    save_plot(file, p, base_height = 5, base_width = 4)
  }
  ok_cell
}



#' plot PCA result
#'
#' Custom PCA output using output from PCA
#' @param pca output from Clustering_PCA function
#' @param file graph image
#' @param legend_file output of legend graph
#' @param title title of graph
#' @param Xcom Specify principal component at X-axis. (Default = 1)
#' @param Ycom Specify principal component at Y-axis. (Default = 2)
#' @param cell_table data frame having cell informations. 1 of column should be "Cell". Other column could be use for specify coloring, shape, size information.
#' @param color_by specify color target column
#' @param shape_by specify shape target column
#' @param size_by specify size target column
#' @param width default = 5
#' @param height default = 5
#' @param pallete color vectors
#' @param option other options for ggplot drawing
#' @return graph file
#' @export
#' @examples
#' plot_PCA(pca, file="output.png")
plot_PCA <- function(pca, file=NULL, legend_file=NULL, title=NULL, Xcom=1, Ycom=2, cell_table=NULL, color_by=NULL, shape_by=NULL,
                     size_by=NULL, width=5, height=5, alpha=0.5, pallete = NULL, option=NULL){
  ### cell_tableはdata.frame
  # Cellというカラムが定義されていること！
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(cowplot))
  D_table <- data.frame(Cell=rownames(pca$x), x=pca$x[,Xcom], y=pca$x[,Ycom], stringsAsFactors = FALSE)
  if(!is.null(cell_table)){
    D_table <- dplyr::left_join(D_table, cell_table, by="Cell", copy=FALSE)
  }

  # 寄与率
  contribution <- pca$sdev^2/sum(pca$sdev^2)*100

  if(is.null(color_by)){
    p <- ggplot(D_table, aes_string(x="x", y="y", colour=color_by, size=size_by, shape=shape_by, stroke=0)) + geom_point(alpha=alpha) +
      labs(x=paste("PC", Xcom, " (", format(contribution[Xcom], digits = 3), "%)", sep=""),
           y=paste("PC", Ycom, " (", format(contribution[Ycom], digits = 3), "%)", sep=""), title=title)
  }else{
    if(length(unique(D_table[,color_by])) < 20){
      if(is.null(pallete)){
        pallete <-c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
      }
      colors <- colorRampPalette(pallete)(length(unique(D_table[,color_by])))
      p <- ggplot(D_table, aes_string(x="x", y="y", colour=factor(D_table[,color_by]), size=size_by, shape=shape_by, stroke=0)) + geom_point(alpha=alpha) +
        scale_color_manual(values=colors, name=color_by) +
        labs(x=paste("PC", Xcom, " (", format(contribution[Xcom], digits = 3), "%)", sep=""),
             y=paste("PC", Ycom, " (", format(contribution[Ycom], digits = 3), "%)", sep=""), title=title)
    }else{
      mid <- median(D_table[,color_by], na.rm = TRUE)
      p <- ggplot(D_table, aes_string(x="x", y="y", colour=color_by, size=size_by, shape=shape_by, stroke=0)) + geom_point(alpha=alpha) +
        scale_color_gradient2(midpoint=mid, low="blue", mid="grey90",　high="red", space ="Lab" )+
        labs(x=paste("PC", Xcom, " (", format(contribution[Xcom], digits = 3), "%)", sep=""),
             y=paste("PC", Ycom, " (", format(contribution[Ycom], digits = 3), "%)", sep=""), title=title)
    }
  }
  if(!is.null(option)){
    p <- p + option
  }
  if(!is.null(legend_file)){
    save_plot(legend_file, get_legend(p))
  }
  if(!is.null(file)){
    save_plot(file, p + theme(legend.position="none"), base_height = height, base_width = width)
  }else{
    p
  }
}


#' perform tSNE clustering
#'
#' Perofrming tSNE clustering using PCA result as input
#' @param DATA pca result or data matrix
#' @param pca_use which principal component will be use for tSNE
#' @param theta 0.0 is original tSNE. 0.5 is quick version
#' @param perplexity perplexity. this value cannot exceed (nrow(mat) -1)/3
#' @param max_iter maximum number of iteration
#' @param seed seed
#' @return tsne result object
#' @export
#' @examples
#' tsne = Clustering_tSNE(pca, pca_use=1:10, perplexity=20)
Clustering_tSNE <- function(DATA, pca_use = 1:10, seed=1, perplexity = 20, theta = 0.0, max_iter=1000){
  Data_type = as.character(class(DATA))
  if(Data_type == "prcomp"){
    mat <- DATA$x[,pca_use]
    rownames(mat) <- rownames(DATA$x)
  }else{
    mat <- DATA
  }

  ### check perplexity maximum
  if(perplexity >= (nrow(mat) -1)/3){
    cat(paste0(perplexity, " is too large\n"))
    q()
  }

  library(Rtsne)
  set.seed(seed) # 再現性の確保
  tsne <- Rtsne(mat, dims = 2, verbose = FALSE, perplexity = perplexity, pca=FALSE,
                theta = theta, max_iter=max_iter, check_duplicates=FALSE, num_threads=2)
  tsne[["Cell"]] <- rownames(mat)
  tsne
}


#' plot tSNE result
#'
#' plotting tSNE results
#' @param tsne result from Clustering_tSNE
#' @param file graph image
#' @param legend_file output of legend graph
#' @param title title of graph
#' @param cell_table data frame having cell informations. 1 of column should be "Cell". Other column could be use for specify coloring, shape, size information.
#' @param color_by specify color target column
#' @param shape_by specify shape target column
#' @param size_by specify size target column
#' @param width default = 5
#' @param height default = 5
#' @param pallete color vectors
#' @param option other options for ggplot drawing
#' @return graph of tSNE result
#' @export
#' @examples
#' plot_tSNE(tsne, file="output.png")
plot_tSNE <- function(tsne, file=NULL, legend_file=NULL, title=NULL,  cell_table=NULL, color_by=NULL, shape_by=NULL, size_by=NULL,
                      width=5, height=5, alpha=0.5, pallete = NULL, option=NULL){
  ### cell_tableはdata.frame
  # Cellというカラムが定義されていること！
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(cowplot))
  D_table <- data.frame(Cell=tsne$Cell, x=tsne$Y[,1], y=tsne$Y[,2], stringsAsFactors = FALSE)
  if(!is.null(cell_table)){
    D_table <- dplyr::left_join(D_table, cell_table, by="Cell", copy=FALSE)
  }
  if(is.null(color_by)){
    p <- ggplot(D_table, aes_string(x="x", y="y", colour=color_by, size=size_by, shape=shape_by, stroke=0)) + geom_point(alpha=alpha) +
      labs(x="tSNE_1", y="tSNE_2", title=title)
  }else{
    if(length(unique(D_table[,color_by])) < 20){
      if(is.null(pallete)){
        pallete <-c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
      }
      colors <- colorRampPalette(pallete)(length(unique(D_table[,color_by])))
      p <- ggplot(D_table, aes_string(x="x", y="y", colour=factor(D_table[,color_by]), size=size_by, shape=shape_by, stroke=0)) + geom_point(alpha=alpha) +
        scale_color_manual(values=colors, name=color_by) +
        labs(x="tSNE_1", y="tSNE_2", title=title)
    }else{
      mid <- median(D_table[,color_by], na.rm = TRUE)
      p <- ggplot(D_table, aes_string(x="x", y="y", colour=color_by, size=size_by, shape=shape_by, stroke=0)) + geom_point(alpha=alpha) +
        scale_color_gradient2(midpoint=mid, low="blue", mid="grey90",　high="red", space ="Lab" )+
        labs(x="tSNE_1", y="tSNE_2", title=title)
    }
  }
  if(!is.null(option)){
    p <- p + option
  }
  if(!is.null(legend_file)){
    save_plot(legend_file, get_legend(p))
  }
  if(!is.null(file)){
    save_plot(file, p + theme(legend.position="none"), base_height = height, base_width = width)
  }else{
    p
  }
}




#' Perform clustering using SIMILR
#'
#' perform clustering using SIMILR method
#' @param DATA output from clustering_PCA
#' @return similr output object
#' @export
#' @examples
#' sim = Clustering_SIMILR(pca, pca_use = 1:10)
Clustering_SIMILR <- function(DATA, pca_use = 1:10, seed=1){
  suppressWarnings(suppressMessages(library(SIMLR)))
  suppressWarnings(suppressMessages(library(igraph)))

  result <- list()

  Data_type = as.character(class(DATA))
  if(Data_type == "prcomp"){
    mat <- DATA$x[,pca_use]
    rownames(mat) <- rownames(DATA$x)
  }else{
    mat <- DATA
  }

  ### 上位3位のclusterを出す
  SIM_estimate <- SIMLR_Estimate_Number_of_Clusters(mat, 2:10, cores.ratio = 0)
  SIM_estimate <- data.frame(cluster_num=2:10, value=SIM_estimate$K1)
  SIM_estimate <- SIM_estimate %>% filter(min_rank(value) < 4) %>% arrange(desc(value)) %>% pull(cluster_num)

  result[["SIM_estimate"]] <- SIM_estimate
  result[["Cell"]] <- rownames(mat)

  i=1
  set.seed(seed)
  for(cluster in SIM_estimate){
    sim <- SIMLR(X = t(mat), c = cluster)
    colnames(sim$ydata) <- paste0("SIMLR", i, "_", c("X", "Y"))
    if(i==1){
      result[["X"]] <- sim$ydata[,1]
      result[["Y"]] <- sim$ydata[,2]
    }else{
      result[["X"]] <- cbind(result[["X"]], sim$ydata[,1])
      result[["Y"]] <- cbind(result[["Y"]], sim$ydata[,2])
    }
    i <- i+1
    # df <- data.frame(sim$ydata)
  }
  result
}


#' make graph of SIMILR result
#'
#' make graph for SIMILR results
#' @param sim Clustering_SIMILR result object
#' @param target which rank of SIMILR reulst will be output (default : best result)
#' @param file graph image
#' @param legend_file output of legend graph
#' @param title title of graph
#' @param Xcom Specify principal component at X-axis. (Default = 1)
#' @param Ycom Specify principal component at Y-axis. (Default = 2)
#' @param cell_table data frame having cell informations. 1 of column should be "Cell". Other column could be use for specify coloring, shape, size information.
#' @param color_by specify color target column
#' @param shape_by specify shape target column
#' @param size_by specify size target column
#' @param width default = 5
#' @param height default = 5
#' @param pallete color vectors
#' @param option other options for ggplot drawing
#' @return graphs
#' @export
#' @examples
#' plot_SIMILR(sim, file="output.png")
plot_SIMILR <- function(sim, file=NULL, target=1, legend_file=NULL, title=NULL,  cell_table=NULL, color_by=NULL, shape_by=NULL, size_by=NULL,
                      width=5, height=5, alpha=0.5, pallete = NULL, option=NULL){
  ### cell_tableはdata.frame
  # Cellというカラムが定義されていること！
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(cowplot))
  D_table <- data.frame(Cell=sim$Cell, x=sim$X[,target], y=sim$Y[,target], stringsAsFactors = FALSE)
  if(!is.null(cell_table)){
    D_table <- dplyr::left_join(D_table, cell_table, by="Cell", copy=FALSE)
  }
  if(is.null(color_by)){
    p <- ggplot(D_table, aes_string(x="x", y="y", colour=color_by, size=size_by, shape=shape_by, stroke=0)) + geom_point(alpha=alpha) +
      labs(x="SIMILR_1", y="SIMILR_2", title=title)
  }else{
    if(length(unique(D_table[,color_by])) < 20){
      if(is.null(pallete)){
        pallete <-c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
      }
      colors <- colorRampPalette(pallete)(length(unique(D_table[,color_by])))
      p <- ggplot(D_table, aes_string(x="x", y="y", colour=factor(D_table[,color_by]), size=size_by, shape=shape_by, stroke=0)) + geom_point(alpha=alpha) +
        scale_color_manual(values=colors, name=color_by) +
        labs(x="SIMILR_1", y="SIMILR_2", title=title)
    }else{
      mid <- median(D_table[,color_by], na.rm = TRUE)
      p <- ggplot(D_table, aes_string(x="x", y="y", colour=color_by, size=size_by, shape=shape_by, stroke=0)) + geom_point(alpha=alpha) +
        scale_color_gradient2(midpoint=mid, low="blue", mid="grey90",　high="red", space ="Lab" )+
        labs(x="SIMILR_1", y="SIMILR_2", title=title)
    }
  }
  if(!is.null(option)){
    p <- p + option
  }
  if(!is.null(legend_file)){
    save_plot(legend_file, get_legend(p))
  }
  if(!is.null(file)){
    save_plot(file, p + theme(legend.position="none"), base_height = height, base_width = width)
  }else{
    p
  }
}

