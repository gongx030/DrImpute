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
