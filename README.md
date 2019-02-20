[![Travis Build Status](https://travis-ci.org/pneuvial/c3co.svg?branch=master)](https://travis-ci.org/pneuvial/c3co)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/pneuvial/c3co?branch=master&svg=true)](https://ci.appveyor.com/project/pneuvial/c3co)
[![Coverage Status](https://img.shields.io/codecov/c/github/pneuvial/c3co/master.svg)](https://codecov.io/github/pneuvial/c3co?branch=master)

# c3co: Inferring cancer cell clonality from copy number data

## Installation

```r
remotes::install_github("pneuvial/c3co", build_vignettes = TRUE)
```

## Model 

This package implements a constraint dictionary learning algorithm to recover subclones across several DNA copy number profiles. The c3co model is described in the figure below. Two heterogeneous tumor samples (green and yellow circles) are composed of a collection of normal cells (gray discs) and two cancer subclones (red triangles and blue squares). One of the cancer subclones is present in both tumor samples. 

![](vignettes/img/features.png)
![](vignettes/img/features2.png)

The corresponding (noiseless) total copy number profiles are displayed in the figure below. They are given by a linear combination of the latent profiles. This Figure is adapted from Nowak *et al.* (Biostatistics, 2011). 

![](vignettes/img/model.png)
![](vignettes/img/model2.png)

The `c3co` method aims at inferring the latent (blue, red and gray) copy-number profiles and the corresponding weights from the observed (yellow and green) copy-number profiles. Importantly, `c3co` takes advantage of the information of parental copy numbers, which are available from allelic ratios (a.k.a. B allele fractions or BAF) as measured by SNP array data or whole exome (WXS) or whole genome sequencing (WGS) data. 

## Documentation

See the package [vignette](vignettes/c3co.Rmd)
