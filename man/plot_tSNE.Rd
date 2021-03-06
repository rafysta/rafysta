% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Clustering.R
\name{plot_tSNE}
\alias{plot_tSNE}
\title{plot tSNE result}
\usage{
plot_tSNE(tsne, file = NULL, legend_file = NULL, title = NULL,
  cell_table = NULL, color_by = NULL, shape_by = NULL,
  size_by = NULL, width = 5, height = 5, alpha = 0.5,
  pallete = NULL, option = NULL, le_width = NULL, le_height = 4)
}
\arguments{
\item{tsne}{result from Clustering_tSNE}

\item{file}{graph image}

\item{legend_file}{output of legend graph}

\item{title}{title of graph}

\item{cell_table}{data frame having cell informations. 1 of column should be "Cell". Other column could be use for specify coloring, shape, size information.}

\item{color_by}{specify color target column}

\item{shape_by}{specify shape target column}

\item{size_by}{specify size target column}

\item{width}{default = 5}

\item{height}{default = 5}

\item{pallete}{color vectors}

\item{option}{other options for ggplot drawing}
}
\value{
graph of tSNE result
}
\description{
plotting tSNE results
}
\examples{
plot_tSNE(tsne, file="output.png")
}
