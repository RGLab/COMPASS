set.seed(123)
n <- 10 ## number of samples
k <- 3 ## number of markers

## generate some sample data
sid_vec <- paste0("sid_", 1:n) ## sample ids; unique names used to denote samples
iid_vec <- rep_len( paste0("iid_", 1:(n/2) ), n ) ## individual ids

## generate n matrices of 'cell intensities'
data <- replicate(n, {
  nrow <- round(runif(1) * 1E2 + 1000)
  ncol <- k
  vals <- rexp( nrow * ncol, runif(1, 1E-5, 1E-3) )
  vals[ vals < 2000 ] <- 0
  output <- matrix(vals, nrow, ncol)
  output <- output[ apply(output, 1, sum) > 0, ]
  colnames(output) <- paste0("M", 1:k)
  return(output)
})
names(data) <- sid_vec

## make a sample metadata data.frame
meta <- data.frame(
  sid=sid_vec,
  iid=iid_vec,
  trt=rep( c("Control", "Treatment"), each=5 )
)

## generate an example total counts
## recall that cells not expressing any marker are not included
## in the 'data' matrices
counts <- sapply(data, nrow) + round( rnorm(n, 1E3, 1E2) )
counts <- setNames( as.integer(counts), names(counts) )

## insert everything into a COMPASSContainer
CC <- COMPASSContainer(data, counts, meta, "iid", "sid")
