# FLAME 2.1.1
## Minor Changes
1. PE computation now based off new `xgboost::xgboost()` defaults.  
2. Introduced `nnet` as a dependency for compatibility with `mice::mice.impute.polyreg()` for imputing categorical variables. 

# FLAME 2.1.0
## Major Changes
1. The Dynamic Almost Matching Exactly (DAME) algorithm has been implemented. This algorithm has the same features as FLAME and is more computationally intensive as it is _guaranteed_ to find the best set of covariates to match on for each unit. The `n_flame_iters` argument in `DAME` specifies that FLAME should be run for a certain number of iterations to more rapidly eliminate irrelevant covariates, before switching to DAME. See the DAME paper [here](https://arxiv.org/pdf/1806.06802.pdf) for more details. 
2. The `data` entry in the list (now of class `ame`) returned by `FLAME` or `DAME` does not contain `*`'d entries for values that units did not match on. This led to some confusing behavior with `DAME` when `replace = TRUE` and also impeded using the matched data for subsequent analysis (e.g. treatment effect estimation via regression adjustment). 
3. Support was removed for automatically binning continuous covariates, as this is generally not recommended. All variables passed in to a matching algorithm (`DAME` or `FLAME`) will be assumed to be categorical and no longer need to be passed as factors.
4. Whereas users were previously forced to input categorical covariates as factors or numerics with values 0:(n_levels - 1), covariates can now be passed with any values / codings, streamlining the process of running a matching algorithm on one's data and rendering the output more interpretable.
5. Support for fixed weights has been introduced. Users can now pass a vector of covariate importances via the `weights` argument that will be used to determine dropping order. This means that the holdout set is no longer necessary (it is ignored if passed) because the weights implicitly define a dropping order and there is no need to compute PE. Note that this dramatically speeds up the algorithm. 
6. The output of `FLAME` or `DAME` is now an S3 object of class `ame`; print, plot, and summary methods have been introduced.

## Minor Changes
1. The vignette has been updated to reflect all changes. 
2. It is no longer possible to perform more than 1 round of imputation for the matching data. Users that wish to analyze the results of multiple imputations should wrap calls to `FLAME` or `DAME` in a loop. 
3. The values that can be passed to `missing_data` and `missing_holdout` have been given interpretable names based off their function.
4. Functionality for using user-supplied functions to compute PE has been changed. Now, a single function should be passed to `PE_method` that does fitting and prediction and the `user_*` arguments are deprecated and will be removed in a later release. See the documentation and vignette for more details.  
5. The arguments `early_stop_pe` and `early_stop_bf` have been deprecated and will be removed in a later release. 
6. Methods for treatment effect estimation have been changed. The `ATE`, `ATT`, `ATC` functions have been deprecated and will be removed in a later release, as these effect estimates and their variances are computed by `summary.ame`, albeit slightly differently from before; see the documentation for more details. The `CATE` function also estimates CATEs slightly differently and returns estimate variances.
7. CATE estimates are no longer automatically returned; this can be toggled with the `estimate_CATEs` flag, which defaults to `FALSE`. 
8. Small changes in output formats of `MG` and `CATE`. 
9. Support for multi-class outcomes. 
10. Some information about the data and matching process is now returned in the `info` entry of an `ame` object.
11. Eliminated dependence on `xgboost` and `mice`; now only required if using XGBoost to compute PE or if imputing data.
12. Rownames in the matching data are now unchanged and `MG` returns matched groups
with respect to the rownames of the matching data and *not* with respect to the index of the unit in the original data frame. Typically, this will also correspond to the indexing of the data (i.e. the third unit has rowname '3'); however, if a separate holdout set was not passed to the matching algorithm or if the original matching data had rownames other than `1:nrow(data)` then this is not the case.
    - The `index_only` argument of `MG` has therefore been replaced with the more appropriate `id_only` argument. 

# FLAME 2.0.0
Version released to CRAN April 15, 2020. Introduced functionality for matching with FLAME.
