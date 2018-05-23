## DrImpute: Imputing dropout events in single cell RNA sequencing data

Wuming Gong (<gongx030@umn.edu>) and Il-Youp Kwak (<ilyoup.kwak@gmail.com>).

R/DrImpute is an R package for imputing dropout events in single-cell RNA-sequencing data. It improve many statistical tools used for scRNA-seq analysis that do not account for dropout events. 

Details are described [here](http://www.biorxiv.org/content/early/2017/08/28/181479)


## 1. Installation

The recommended installation method for `DrImpute` is using `install_github` command from `devtools` library.  You will first have to have [devtools](https://github.com/hadley/devtools) package installed.

```r
library(devtools)
install_github('gongx030/DrImpute')
```

A number of needed packages are installed in this process.

## 2. Quick Start

We first load the tcm package:
```r
library(DrImpute)
```

We load the scRNA-seq dataset on mouse sensory neurons that have 622 cells from four major types of cells: NP, TH, NF and PEP.  
```r
install_github('gongx030/scDataset')
library(scDataset)
usoskin
```

We extract the read count matrix and the cell labels. The the genes that are expressed in zero or only one cell will be removed by the function `preprocessSC`.  
```r
X <- assays(usoskin)$count
X <- preprocessSC(X)
X.log <- log(X + 1)
class.label <- colData(usoskin)[['Level.1']]
```

Then we run the imputation with the DrImpute:
```r
set.seed(1)
X.imp <- DrImpute(X.log)
```

We can compare the clustering performance (t-SNE followed by k-means) by using the scRNA datasets before and after the imputation:
```r
library(mclust)
library(Rtsne)
# before imputation
adjustedRandIndex(kmeans(Rtsne(t(as.matrix(X.log)))$Y, centers = 4)$cluster, class.label)
# after imputation
adjustedRandIndex(kmeans(Rtsne(t(as.matrix(X.imp)))$Y, centers = 4)$cluster, class.label)
```
