swap <- function(vec, from, to) {
  tmp <- to[match(vec, from)]
  tmp[is.na(tmp)] <- vec[is.na(tmp)]
  return(tmp)
}
