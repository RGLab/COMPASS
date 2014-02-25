## Clamp a number to a specific range
clamp <- function(x, low, high) {
  if (x < low) return(low)
  else if (x > high) return(high)
  else return(x)
}

## Swap all elements in a vector 'vec' in 'from' to
## corresponding element in 'to'
.swap <- function(vec, from, to) {
  tmp <- to[match(vec, from)]
  tmp[is.na(tmp)] <- vec[is.na(tmp)]
  return(tmp)
}
