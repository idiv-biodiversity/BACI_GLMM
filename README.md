# Overview

This repository contains the R scripts used for performing the analysis in the article:

> Pardini, E. A., Parsons, L. S., Ştefan, V., & Knight, T. M. (2018). GLMM BACI environmental impact analysis shows coastal dune restoration reduces seed predation on an endangered plant. Restoration Ecology, 26(6), 1190-1194.

Link to article on Wiley Online Library [here](http://onlinelibrary.wiley.com/doi/10.1111/rec.12678/full).

The dataset is uploaded with a DOI on zenodo.org, and it can be downloaded from [here](https://zenodo.org/record/4384785). Use the script `get_data.r` to download the data and its metadata.

# How to use this repository

## - Clone or download

You can [download][1] or clone the repository then run the scripts using the *BACI_GLMM.Rproj* file ([R][2] and [R Studio][3] are needed).

For cloning, run this in a terminal (git should be [installed][4]):

```
git clone https://github.com/idiv-biodiversity/BACI_GLMM.git
```

First run the script `get_data.r` to download the data from zenodo.org.

## - Install older versions of R and packages

For increasing reproducibility, one should install the older version of R (3.4.2) and also the older packages.

For installing multiple versions of R on a operating system see these resources:

- https://support.rstudio.com/hc/en-us/articles/360002242413-Multiple-versions-of-R 
- https://support.rstudio.com/hc/en-us/articles/215488098-Installing-multiple-versions-of-R-on-Linux

In the table `package_table.txt` one can find the package names, their versions and their CRAN archive url. These are the packages listed by the command `sessionInfo()` in the script `GLMM_BACI_public.R`.

For example to install the package 'lsmeans', version '2.27-2', one can do:
```r
install.packages(pkgs = "https://cran.r-project.org/src/contrib/Archive/lsmeans/lsmeans_2.27-2.tar.gz",
                 lib = local_lib,
                 type = "source",
                 repos = NULL)
```
You can run into package dependencies issues, in which case pay attention to the error messages and install the needed dependencies. The order of installing dependencies matters.

[1]: https://github.com/idiv-biodiversity/BACI_GLMM/archive/master.zip
[2]: https://www.r-project.org/
[3]: https://www.rstudio.com/products/rstudio/download/
[4]: https://git-scm.com/downloads