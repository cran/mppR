mppR: Multi-Parent Population QTL Analysis
====


## Overview

mppR is an R package to perform QTL analysis of experimental multi-parent populations. The population must be composed of crosses between a set of at least three parents (e.g. factorial design, 'diallel', or nested association mapping). The functions cover data processing, QTL detection, and results visualization.

## Installation

mppR has two different branches: "master" and "mppR_CRAN". The "master" branch allows to perform MPP mixed model QTL detection calling the asreml-R package and function parent_cluster.mppData that call the archived R package clusthaplo for parent clustering. The branch "mppR_CRAN" do not contain the mixed models and the call to clusthaplo.

```
devtools::install_github("vincentgarin/mppR", ref = "master")

```

## Usage

See the two vignettes attached to the package.

# Travis

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/mppR)](https://cran.r-project.org/package=mppR)
