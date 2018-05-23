## R package, DrImpute

Il-Youp Kwak (<ilyoup.kwak@gmail.com>), with contributions from Wuming Gong.

R/DrImpute is an R package for imputing dropout events in single-cell RNA-sequencing data. It improve many statistical tools used for scRNA-seq analysis that do not account for dropout events. 

Details are described [here](http://www.biorxiv.org/content/early/2017/08/28/181479)


### installation
From `CRAN` :
```S
install.packages("DrImpute")
```

Or, with `devtools`:
```S
library(devtools)
install_github("ikwak2/DrImpute")
```

### License

The R/DrImpute package is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License,
version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but
without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the GNU
General Public License for more details.

A copy of the GNU General Public License, version 3, is available at
<https://www.r-project.org/Licenses/GPL-3>
