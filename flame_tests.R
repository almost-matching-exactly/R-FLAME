library(MatchIt)
data(lalonde)

lalonde$age %<>% as.factor()
lalonde$educ %<>% as.factor()
lalonde$black %<>% as.factor()
lalonde$hispan %<>% as.factor()
lalonde$married %<>% as.factor()
lalonde$nodegree %<>% as.factor()

lalonde %<>% dplyr::select(-c(re74, re75))

FLAME_out <- FLAME(lalonde, holdout=lalonde, treated_column_name = 'treat', outcome_column_name = 're78')
print(ATE(FLAME_out)) # Can't tell if the Experimental ate in this subsample is ~1700 or ~800.
                      # We're doing well in the second case, not so much in the first.

# Testing verbosity settings
FLAME_out <- FLAME(lalonde, holdout=lalonde, treated_column_name = 'treat', outcome_column_name = 're78',
                   verbose=3)

# Testing replacement
FLAME_out <- FLAME(lalonde, holdout=lalonde, treated_column_name = 'treat', outcome_column_name = 're78',
                   verbose=3, replace=T)
print(ATE(FLAME_out)) # Still pretty good if true ate ~= 800

# Testing PE/BF outputs
FLAME_out <- FLAME(lalonde, holdout=lalonde, treated_column_name = 'treat', outcome_column_name = 're78',
                   verbose=3, replace=T, want_pe=TRUE, want_bf=TRUE)

FLAME_out <- FLAME(lalonde, holdout=lalonde, treated_column_name = 'treat', outcome_column_name = 're78',
                   verbose=3, replace=T, return_pe=TRUE, return_bf=TRUE)
print(sqrt(FLAME_out$PE)/mean(lalonde$re78)) # Is the default PE MSE? Is the interpretation of this
                                             # output as Prediction Error ~= 1.5 average outcome correct?

# Different binning methods
FLAME_out <- FLAME(lalonde, holdout=lalonde, treated_column_name = 'treat', outcome_column_name = 're78',
                   verbose=3, replace=T, binning_method='sturges')
print(ATE(FLAME_out))
FLAME_out <- FLAME(lalonde, holdout=lalonde, treated_column_name = 'treat', outcome_column_name = 're78',
                   verbose=3, replace=T, binning_method='scott')


print(ATE(FLAME_out)) # Both this and the one before are closer to ~1700 interestingly.

FLAME_out <- FLAME(lalonde, holdout=lalonde, treated_column_name = 'treat', outcome_column_name = 're78',
                   verbose=3, replace=T, binning_method='fd')
###### ERROR HERE!!!! #######
# Error in seq.int(rx[1L], rx[2L], length.out = nb) :
# 'length.out' must be a non-negative number
# In addition: Warning messages:
# 1: In bin_continuous_covariates(tmp_df[, 1:n_covs, drop = FALSE], rule = binning_method,  :
# Binning continuous covariates. This is not recommended; users are encouraged
# to use methods specifically designed for continuous covariates.
# 2: In cut.default(X_cont[, i], breaks = n_bins[i], labels = 0:(n_bins[i] -  :
# NAs introduced by coercion to integer range
###########################################
print(ATE(FLAME_out))

# Testing stopping rules
FLAME_out <- FLAME(lalonde, holdout=lalonde, treated_column_name = 'treat', outcome_column_name = 're78',
                   verbose=3, replace=T, binning_method='scott', early_stop_iterations=3)

FLAME_out <- FLAME(lalonde, holdout=lalonde, treated_column_name = 'treat', outcome_column_name = 're78',
                   verbose=3, replace=T, binning_method='scott', early_stop_epsilon=0.001)

FLAME_out <- FLAME(lalonde, holdout=lalonde, treated_column_name = 'treat', outcome_column_name = 're78',
                   verbose=3, replace=T, binning_method='scott', early_stop_treated = 0.1)

FLAME_out <- FLAME(lalonde, holdout=lalonde, treated_column_name = 'treat', outcome_column_name = 're78',
                   verbose=3, replace=T, binning_method='scott', early_stop_control = 0.1)

FLAME_out <- FLAME(lalonde, holdout=lalonde, treated_column_name = 'treat', outcome_column_name = 're78',
                   verbose=3, replace=T, binning_method='scott', early_stop_pe = 2) # Might be better to
# use a standardized pe because users won't know the scale (like here it's in the millions)

FLAME_out <- FLAME(lalonde, holdout=lalonde, treated_column_name = 'treat', outcome_column_name = 're78',
                   verbose=3, replace=T, binning_method='scott', early_stop_bf = 1)

FLAME_out <- FLAME(lalonde, holdout=lalonde, treated_column_name = 'treat', outcome_column_name = 're78',
                   verbose=3, replace=T, binning_method='scott', early_stop_bf = 1, early_stop_pe = Inf)

FLAME_out <- FLAME(lalonde, holdout=lalonde, treated_column_name = 'treat', outcome_column_name = 're78',
                   verbose=3, replace=T, binning_method='scott', early_stop_bf = 0.2,
                   early_stop_iterations = 3, early_stop_epsilon = 0.01)

FLAME_out <- FLAME(lalonde, holdout=lalonde, treated_column_name = 'treat', outcome_column_name = 're78',
                   verbose=3, replace=T, binning_method='scott', early_stop_bf = 0.2,
                   early_stop_iterations = 5, early_stop_epsilon = 0.01, early_stop_control = 0.8)

# Testing PE methods
FLAME_out <- FLAME(lalonde, holdout=lalonde, treated_column_name = 'treat', outcome_column_name = 're78',
                   verbose=3, replace=F, binning_method='scott', PE_method='xgb', return_pe = T)
print(sqrt(FLAME_out$PE))
plot(1:length(FLAME_out$PE), FLAME_out$PE, 'o')

FLAME_out <- FLAME(lalonde, holdout=lalonde, treated_column_name = 'treat', outcome_column_name = 're78',
                   verbose=3, replace=F, binning_method='scott', PE_method='ridge', return_pe = T)

print(sqrt(FLAME_out$PE))
plot(1:length(FLAME_out$PE), FLAME_out$PE, 'o')
# PE is not monotonically increasing in either case (assuming that each value is PE at each iteration).
# This is strange. Is is because of regularization being taken into account in the PE formula?
# Should we have this behavior?

# OLS
custom_lm = function(X, y){
  lm(X ~ y)
}
custom_predict = function(lm_fit, newdata){
  predict(lm_fit, newdata=data.frame(newdata))
}
FLAME_out <- FLAME(lalonde, holdout=lalonde, treated_column_name = 'treat', outcome_column_name = 're78',
                   verbose=3, replace=F, binning_method='scott',
                   user_PE_fit = custom_lm, user_PE_predict=custom_predict, return_pe = T)
print(sqrt(FLAME_out$PE))
plot(1:length(FLAME_out$PE), FLAME_out$PE, 'o')
# Now PE is monotonic (meaning that yes, it was the regularization)
# but decreasing? Shouldn't it be increasing? Are iterations oredered from last to first maybe?

# Works well but I needed a bit of fiddling to get it to work.
# I could do it quickly because I have an idea of what's happening
# under the hood and why the default functions weren't working but you might want to expand more in the
# docs on how to do this. Maybe you could have this exact example here to show what the custom functions
# need to satisfy?

# BART
custom_bart = function(X, y){
  bart(x.train=X, y.train=y, keeptrees = T, verbose=FALSE)
}
custom_predict = function(fit, newdata){
  colMeans(predict(fit, test=newdata))
}
FLAME_out <- FLAME(lalonde, holdout=lalonde, treated_column_name = 'treat', outcome_column_name = 're78',
                   verbose=3, replace=F, binning_method='scott',
                   user_PE_fit = custom_bart, user_PE_predict=custom_predict, return_pe = T)
print(sqrt(FLAME_out$PE))
plot(1:length(FLAME_out$PE), FLAME_out$PE, 'o')
# Now PE is monotonically increasing. This is exactly what I expect this to look like.
print(ATE(FLAME_out)) # Interesting behavior. PE with BART is lowest I get but ATE is completely wrong.

# Testing Missing data
lalonde_mis = lalonde
vvals = rep(FALSE, nrow(lalonde)*ncol(lalonde))
vvals[sample(nrow(lalonde)*ncol(lalonde), 100, F)] = TRUE
lalonde_mis[matrix(vvals, nrow(lalonde), ncol(lalonde))] = NA

FLAME_out <- FLAME(lalonde_mis, holdout=lalonde, treated_column_name = 'treat',
                   outcome_column_name = 're78',
                   verbose=3, replace=F, binning_method='scott')
# Got the correct error

FLAME_out <- FLAME(lalonde_mis, holdout=lalonde, treated_column_name = 'treat',
                   outcome_column_name = 're78',
                   verbose=3, replace=F, binning_method='scott',
                   missing_data=1)


FLAME_out <- FLAME(lalonde_mis, holdout=lalonde, treated_column_name = 'treat',
                   outcome_column_name = 're78',
                   verbose=3, replace=F, binning_method='scott',
                   missing_data=2)

FLAME_out <- FLAME(lalonde_mis, holdout=lalonde, treated_column_name = 'treat',
                   outcome_column_name = 're78',
                   verbose=3, replace=F, binning_method='scott',
                   missing_data=3)

FLAME_out <- FLAME(lalonde, holdout=lalonde_mis, treated_column_name = 'treat',
                   outcome_column_name = 're78',
                   verbose=3, replace=F, binning_method='scott',
                   missing_holdout=1)

FLAME_out <- FLAME(lalonde, holdout=lalonde_mis, treated_column_name = 'treat',
                   outcome_column_name = 're78',
                   verbose=3, replace=F, binning_method='scott',
                   missing_holdout=2)

FLAME_out <- FLAME(lalonde, holdout=lalonde_mis, treated_column_name = 'treat',
                   outcome_column_name = 're78',
                   verbose=3, replace=F, binning_method='scott',
                   missing_holdout=3)
# Got the correct error

FLAME_out <- FLAME(lalonde_mis, holdout=lalonde_mis, treated_column_name = 'treat',
                   outcome_column_name = 're78',
                   verbose=3, replace=F, binning_method='scott',
                   missing_holdout=1, missing_data=1)

FLAME_out <- FLAME(lalonde_mis, holdout=lalonde_mis, treated_column_name = 'treat',
                   outcome_column_name = 're78',
                   verbose=3, replace=F, binning_method='scott',
                   missing_holdout=2, missing_data=2)

FLAME_out <- FLAME(lalonde_mis, holdout=lalonde_mis, treated_column_name = 'treat',
                   outcome_column_name = 're78',
                   verbose=3, replace=F, binning_method='scott',
                   missing_holdout=1, missing_data=3)






