factor_remap <-
  function(df, treated_column_name, outcome_column_name, mapping = NULL) {
  browser()
  cov_inds <- which(!(colnames(df) %in% c(treated_column_name, outcome_column_name)))
  mappings <- list() # Make sure we re-arrange the mappings after we sort!!
  counter <- 1
  for (ind in cov_inds) {
    if (is.factor(df[[ind]])) {
      if (is.null(mapping)) {
        mapping_to_pass <- mapping
      }
      else {
        mapping_to_pass <- mapping[[counter]]
      }
      remapped_col <-
        factor_remap_column(df[[ind]], mapping_to_pass)
      df[[ind]] <- remapped_col$column
      mappings <- c(mappings, list(remapped_col$mapping))

      counter <- counter + 1
    }
  }
  return(list(df = df,
              mappings = mappings))
}

# Check this doesn't break anything if the factor is alright to begin with
factor_remap_column <- function(column, mapping) {
  if (is.null(mapping)) {
    mapping <- levels(column)
    names(mapping) <- as.character(0:(length(mapping) - 1))
  }
  column <- do.call(forcats::fct_recode, c(list(column), mapping))
  return(list(column = column,
              mapping = mapping))
}
