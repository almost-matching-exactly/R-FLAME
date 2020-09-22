factor_remap <-
  function(df,
           treated_column_name = 'treated',
           outcome_column_name = 'outcome',
           mapping = NULL) {

  # Maps the content of mapping to the names of mapping, if supplied

  cov_inds <- which(!(colnames(df) %in% c(treated_column_name, outcome_column_name)))
  mappings <- list()
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

  # Not necessary when mapping original levels to 0:(k - 1),
  # but will come into play if missing_data = 3 to drop the 'fake' levels
  column <- droplevels(column)

  if (is.null(mapping)) {
    mapping <- levels(column)
    names(mapping) <- as.character(0:(length(mapping) - 1))
  }
  levels(column) <- as.list(mapping)
  return(list(column = column,
              mapping = mapping))
}
