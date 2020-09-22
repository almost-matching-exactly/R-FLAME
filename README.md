[![Build Status](https://travis-ci.com/vittorioorlandi/FLAME.svg?branch=master)](https://travis-ci.com/vittorioorlandi/FLAME)

# Introduction
The FLAME package offers an efficient implementation of the Fast, Large-Scale, Almost Matching Exactly algorithm, described in detail [here](https://arxiv.org/pdf/1707.06315.pdf). FLAME allows for interpretable matching within observational settings in causal inference. It does so by matching units via a learned, weighted Hamming distance that determines which covariates are more important to match on. The package will soon be updated to include other algorithms in the Almost Matching Exactly framework, such as [DAME](https://arxiv.org/pdf/1806.06802.pdf), [FLAME-Networks](https://arxiv.org/pdf/2003.00964.pdf), and [FLAME-IV](http://auai.org/uai2019/proceedings/papers/410.pdf).

# An Example
We can generate some toy data for matching by using the `gen_data` function included in the package:
```
matching_data <- gen_data()
holdout_data <- gen_data()
```

We can then match the units in `matching_data`, by using `holdout_data` to learn a weighted Hamming distance between units by calling the main function in the package:
```
FLAME_out <- FLAME(matching_data, holdout_data)
```

The matched data, which can be accessed via `FLAME_out$data` contains asterisks for the covariates on which units were not matched. The groups of units matched together are stored in `FLAME_out$MGs`, the corresponding conditional average treatment effect (CATE) estimates in `FLAME_out$CATEs`, the covariates values units were matched on in `FLAME_out$matched_on`, and the covariates at every iteration of FLAME by `FLAME_out$matching_covs`. 

For a more in-depth discussion of the package, please see the accompanying vignette.  

# Installation 
FLAME can be installed directly from CRAN via running `install.packages('FLAME')` in R. This will also install some tidyverse dependencies, `gmp` for FLAME's efficient bit-vectors implementation, `mice` for handling missing data, and `glmnet` and `xgboost` for outcome prediction. 

Alternatively, you can install the development version of the package from the author's GitHub, [here](https://github.com/vittorioorlandi/FLAME) via `devtools::install_github('https://github.com/vittorioorlandi/FLAME')`.
