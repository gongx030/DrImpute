\docType{data}
\name{exdata}
\alias{exdata}
\title{Usoskin data}
\usage{
data(exdata)
}
\description{

  This data set is subset from Usoskin et al. (2015). Original data is
  RNA-seq data on 799 cells dissected from the mouse lumbar dorsal root
  ganglion distributed over a total of nine 96-well plates. We randomly
  selected 150 cells from the data. 

  Column names indicate four different cell types, NF, NP, TH, and PEP.
}
\examples{
data(exdata)
exdata <- preprocessSC(exdata)
}
\references{
Usoskin D et al. Unbiased classification of sensory neuron types by large-scale
single-cell RNA sequencing. Nature Neuroscience. Nature Research,2015;18:145-53.
}
