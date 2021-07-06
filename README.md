[![Build Status](https://travis-ci.com/vittorioorlandi/FLAME.svg?branch=master)](https://travis-ci.com/vittorioorlandi/FLAME)

# Introduction
The `FLAME` package offers efficient implementations of the Fast, Large-Scale, Almost Matching Exactly (FLAME) and Dynamic Almost Matching Exactly (DAME) algorithms, described in detail [here](https://arxiv.org/pdf/1707.06315.pdf) and [here](https://arxiv.org/pdf/1806.06802.pdf), respectively. Both algorithms are interpretable matching methods for performing causal inference on observational data with *discrete* covariates. The algorithms work by matching units that share identical values of certain covariates; machine learning on a holdout set is used to learn which covariates are more important for predicting the outcome and prioritize matches on those covariates. `FLAME` provides easy summarization, analysis, and visualization of treatment effect estimates, and features a wide variety of options for how matching is to be performed, allowing for users to make analysis-specific decisions throughout the matching procedure. In the future, the package will be updated to include other algorithms in the Almost Matching Exactly framework, such as a database implementation of FLAME / DAME (allowing for inference on data too large to fit in memory), [FLAME-Networks](https://arxiv.org/pdf/2003.00964.pdf) (allowing for inference on data subject to network interference), [AHB](https://arxiv.org/pdf/2003.01805.pdf) (allowing for inference on mixed discrete and continuous data), [MALTS](https://arxiv.org/pdf/1811.07415.pdf) (allowing for inference on mixed discrete and continuous data), and [FLAME-IV](http://auai.org/uai2019/proceedings/papers/410.pdf) (allowing for inference with instrumental variables).

# An Example
We start by generating some toy data for matching by using the `gen_data` function included in the package:
```
matching_data <- gen_data(n = 500, p = 5)
holdout_data <- gen_data(n = 500, p = 5)
```

We can then match the units in `matching_data`, by calling the `FLAME` function, which matches units according to the FLAME algorithm. The supplied `holdout_data` is used to learn a weighted Hamming distance between units that determines what covariates are more important and therefore what matches should be made. 
```
FLAME_out <- FLAME(matching_data, holdout_data)
```

The output of a call to `FLAME` is of class `ame` and the package includes `print`, `summary`, and `plot` methods to give information on average treatment effect estimates, matched groups formed, and covariates matched on.

Using DAME for matching is analogous and can be done using the `DAME` function. 

For a more in-depth discussion of the package that covers missing data handling, machine learning on the holdout set, matching with replacement, and more, please see the accompanying vignette.  

# Installation 
FLAME can be installed directly from CRAN via running `install.packages('FLAME')` in R. This will also install `gmp` for FLAME's efficient bit-vectors implementation and `glmnet` for outcome prediction. 

Alternatively, you can install the development version of the package from the author's GitHub, [here](https://github.com/vittorioorlandi/FLAME) via `devtools::install_github('https://github.com/vittorioorlandi/FLAME')` or the Almost Matching Exactly GitHub [here](https://github.com/almost-matching-exactly/R-FLAME) via 
`devtools::install_github('https://github.com/almost-matching-exactly/R-FLAME')` (the latter is a mirror of the former). 


