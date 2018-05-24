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
install_github('gongx030/scDatasets')
library(scDatasets)
library(SummarizedExperiment)
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

## 3. Session Information
```r
> sessionInfo()
R version 3.4.3 (2017-11-30)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS release 6.9 (Final)

Matrix products: default
BLAS: /panfs/roc/msisoft/R/3.4.3/lib64/R/lib/libRblas.so
LAPACK: /panfs/roc/msisoft/R/3.4.3/lib64/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices datasets  utils     methods   base

other attached packages:
 [1] Rtsne_0.13                 mclust_5.4                 scDatasets_0.0.3           SummarizedExperiment_1.8.1 DelayedArray_0.4.1         matrixStats_0.53.1         Biobase_2.38.0             GenomicRanges_1.30.3       GenomeInfoDb_1.14.0        IRanges_2.12.0             S4Vectors_0.16.0           BiocGenerics_0.24.0
[13] DrImpute_1.1               BiocInstaller_1.28.0

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.17           XVector_0.18.0         xml2_1.2.0             magrittr_1.5           roxygen2_6.0.1         zlibbioc_1.24.0        devtools_1.13.5        lattice_0.20-35        R6_2.2.2               FNN_1.1                stringr_1.3.0          tools_3.4.3            grid_3.4.3             irlba_2.3.2
[15] withr_2.1.2            commonmark_1.4         digest_0.6.15          Matrix_1.2-12          GenomeInfoDbData_1.0.0 bitops_1.0-6           RCurl_1.95-4.10        memoise_1.1.0          stringi_1.1.7          compiler_3.4.3
```
