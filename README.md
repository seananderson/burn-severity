# Code for 'Incorporating biophysical gradients and uncertainty into burn severity maps in a temperate fire-prone forested region'

[![DOI](https://zenodo.org/badge/117417101.svg)](https://zenodo.org/badge/latestdoi/117417101)

This repository contains code associated with the paper:

Harvey, B.J., R.A. Andrus, S.C. Anderson. 2019. Incorporating biophysical
gradients and uncertainty into burn severity maps in a temperate fire-prone
forested region. In press at Ecosphere. <https://doi.org/10.1002/ecs2.2600>

## Contents

```
|- DESCRIPTION               # project metadata
|- README.md                 # this file
|- LICENSE.md                # specifies the conditions of use and reuse of the code
|- burn-severity.Rproj       # RStudio project file
|
|- data/                     # raw data
|- +- generated/             # data that gets generated by code
|
|- analysis/                 # the main analysis files
|  +- 01-burn-severity.Rmd   # analysis for Q1 portion of paper
|  +- 02-burn-severity.Rmd   # analysis for Q2 portion of paper
|  +- 03-burn-severity.Rmd   # analysis for Q3 portion of paper
|  +- 04-burn-severity.Rmd   # code to generate the map figures
|  +- zoib-functions.R       # helper R functions for the analysis
|  +- zoib1re.stan           # Stan ZOIB model
|  +- oib1re.stan            # Stan OIB model (model without the zero component)
|
|  +- figs/                  # figures get generated here
|
|  +- make-map-data/         # Code for generating predictor data going into the
|                            # map predictions. Note that the various spatial files
|                            # are not included in this repository. They were
|                            # processed with the R file in this repository to generate the
|                            # data that are cached in the `data/generated` folder.
```

## Dependencies

In order to run the included code the following R packages are required:

```r
pkgs <- c("dplyr", "ggplot2", "plyr", "reshape2", "readr", "forcats", "here",
  "purrr", "scales", "rstan", "broom", "viridis", "RColorBrewer", "devtools",
  "doParallel", "cowplot")
install.packages(pkgs)
devtools::install_github("seananderson/ggsidekick")
```

You will also need to have a C++ toolchain setup to build the Stan models. See <https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started>.

The spatial data extraction (for which the raw data is not included here), also
required raster, rgdal, dismo, and sp.

The analysis was run with the following computational environment:

```
─ Session info ──────────────────────────────────────────────────────────
 setting  value
 version  R version 3.5.2 (2018-12-20)
 os       macOS Mojave 10.14.2
 system   x86_64, darwin15.6.0
 language (EN)
 collate  en_CA.UTF-8
 ctype    en_CA.UTF-8

─ Packages ──────────────────────────────────────────────────────────────
 package        version   date        source
 broom          0.5.1     2018-12-05  CRAN (R 3.5.0)
 cowplot        0.9.4     2019-01-08  CRAN (R 3.5.2)   
 doParallel     1.0.14    2018-09-24  CRAN (R 3.5.0)
 dplyr          0.7.8     2018-11-10  CRAN (R 3.5.0)
 forcats        0.3.0     2018-02-19  CRAN (R 3.5.0)
 ggplot2        3.1.0     2018-10-25  CRAN (R 3.5.0)
 here           0.1       2017-05-28  CRAN (R 3.5.0)
 loo            2.0.0     2018-04-11  CRAN (R 3.5.0)
 plyr           1.8.4     2016-06-08  CRAN (R 3.5.0)
 purrr          0.3.0     2019-01-27  CRAN (R 3.5.2)
 RColorBrewer   1.1-2     2014-12-07  CRAN (R 3.5.0)
 Rcpp           1.0.0     2018-11-07  CRAN (R 3.5.0)
 readr          1.3.1     2018-12-21  CRAN (R 3.5.0)
 reshape2       1.4.3     2017-12-11  CRAN (R 3.5.0)
 rstan          2.18.2    2018-11-07  CRAN (R 3.5.0)
 scales         1.0.0     2018-08-09  CRAN (R 3.5.0)
 tidyr          0.8.2     2018-10-28  CRAN (R 3.5.0)
 viridis        0.5.1     2018-03-29  CRAN (R 3.5.0)
```
