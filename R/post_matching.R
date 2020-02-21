mmg_of_unit <- function(unit, flame_obj, output_style = 1) {
  MGs <- flame_obj$MGs
  for (i in 1:length(MGs)) {
    if (unit %in% MGs[[i]]) {
      if (output_style == 1) {
        return(flame_obj$data[MGs[[i]], ])
      }
      return(flame_obj$original_data[MGs[[i]], ])
    }
  }
  stop('Unit not present in data')
}

te_of_unit <- function(unit, flame_obj) {
  MGs <- flame_obj$MGs
  for (i in 1:length(MGs)) {
    if (unit %in% MGs[[i]]) {
      return(flame_obj$CATE[i])
    }
  }
  stop('Unit not present in data')
}

ATE <- function(flame_obj) {
  CATE <- flame_obj$CATE
  MGs <- flame_obj$MGs
  size <- sapply(MGs, length)
  return(sum(size * CATE) / sum(size))
}
