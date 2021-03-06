% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DrImpute.R
\name{DrImpute}
\alias{DrImpute}
\title{DrImpute}
\usage{
DrImpute(X, ks = 10:15, dists = c("spearman", "pearson"),
  batch.size = NA, mc.cores = 1)
}
\arguments{
\item{X}{Gene expression matrix (gene by cell).}

\item{ks}{Number of cell clustering groups. Default set to ks = 10:15.}

\item{dists}{Distribution matrices to use. Default is set to c("spearman", "pearson"). "eucleadian" can be added as well.}

\item{batch.size}{Batch size of sampling based MDS approximation; if batch.size is NA, the standard MDS is used.}

\item{mc.cores}{Number of CPU cores (default: 1)}
}
\value{
A matrix object
}
\description{
Imputing dropout events in single-cell RNA-sequencing data.
}
\examples{

library(scDatasets)
library(SummarizedExperiment)
data(usoskin)
X <- assays(usoskin)$count
X <- preprocess(X, min.expressed.gene = 0)
X.log <- log(X + 1)
set.seed(1)
X.imp <- DrImpute(X.log)

}
\references{
Wuming Gong, Il-Youp Kwak, Pruthvi Pota, Kaoko Koyano-Nakagawa and Daniel J. Garry (2017+)
DrImpute: Imputing dropout eveents in single cell RNA sequencing data
}
\author{
Wuming Gong and Il-Youp Kwak
}
