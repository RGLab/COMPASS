set.seed(123)
n <- 20 ## number of samples
k <- 6 ## number of markers

sid_vec <- paste0("sid_", 1:n) ## sample ids; unique names used to denote samples
iid_vec <- rep_len( paste0("iid_", 1:(n/10) ), n ) ## individual ids

data <- replicate(n, {
  nrow <- round(runif(1) * 1E4 + 1000)
  ncol <- k
  vals <- rexp( nrow * ncol, runif(1, 1E-5, 1E-3) )
  vals[ vals < 2000 ] <- 0
  output <- matrix(vals, nrow, ncol)
  output <- output[ apply(output, 1, sum) > 0, ]
  colnames(output) <- paste0("M", 1:k)
  return(output)
})
names(data) <- sid_vec

counts <- sapply(data, nrow) + round( rnorm(n, 1E4, 1E3) )
counts <- setNames( as.integer(counts), names(counts) )

meta <- data.frame(
  sid=sid_vec,
  iid=iid_vec,
  trt=sample( c("Control", "Treatment"), n, TRUE )
)

CC <- COMPASSContainer(
  data=data,
  counts=counts,
  meta=meta,
  individual_id="iid", ## name of individual id vector in metadata
  sample_id="sid" ## name of sample id vector in metadata
)
