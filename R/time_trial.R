# require(microbenchmark)
#
# aggregate_table <- function(list) {
#   tab = table(as.character(list))
#   tab = unclass(tab)
#   name = names(tab)
#   list_val = as.character(list)
#   return(as.vector(tab[sapply(list_val, function(x) which(name ==  x))]))
# }
#
# opt1 <- function(vec) {
#   # tab <- table(vec)
#   # fast_count(vec, as.numeric(names(tab)), tab)
#   tab = table(as.character(list))
#   tab = unclass(tab)
#   name = names(tab)
#   list_val = as.character(list)
#   return(as.vector(tab[sapply(list_val, function(x) which(name ==  x))]))
# }
#
#
#
#
# bu <- sample(1:100, 5000, replace = T)
#
# microbenchmark(aggregate_table(bu), opt1(bu), check = 'identical')
# microbenchmark({table(bu);as.numeric(names(bu))})
