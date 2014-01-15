
## ----sim-init------------------------------------------------------------
library(COMPASS)
set.seed(123)
n <- 100 ## number of samples
k <- 6 ## number of markers

sid_vec <- paste0("sid_", 1:n) ## sample ids; unique names used to denote samples
iid_vec <- rep_len( paste0("iid_", 1:(n/10) ), n ) ## individual ids


## ----sim-data------------------------------------------------------------
data <- replicate(n, {
  nrow <- round(runif(1) * 1E4 + 1000)
  ncol <- k
  vals <- rexp( nrow * ncol, 1/1000 )
  vals[ vals < 1000 ] <- 0
  output <- matrix(vals, nrow, ncol)
  colnames(output) <- paste0("C", 1:k)
  return(output)
})
names(data) <- sid_vec


## ----sim-counts----------------------------------------------------------
counts <- sapply(data, nrow) + round( rnorm(n, 1E4, 1E3) )
counts <- setNames( as.integer(counts), names(counts) )


## ----sim-meta------------------------------------------------------------
meta <- data.frame(
  sid=sid_vec,
  iid=iid_vec,
  trt=sample( c("Control", "Treatment"), n, TRUE )
)


## ----sim-CC--------------------------------------------------------------
CC <- COMPASSContainer(
  data=data,
  counts=counts,
  meta=meta,
  individual_id="iid", ## name of individual id vector in metadata
  sample_id="sid", ## name of sample id vector in metadata
  stimulation_id="trt" ## what treatment was applied to each sample?
)


## ----CC-basics-----------------------------------------------------------
CC
summary(CC)


## ----COMPASS-fit---------------------------------------------------------
fit <- COMPASS( CC,
  treatment=trt == "Treatment",
  control=trt == "Control",
  iterations=1000
)


## ----COMPASS-examine-----------------------------------------------------
## Extract the functionality, polyfunctionality scores as described
## within the COMPASS paper
FS <- FunctionalityScore(fit)
PFS <- PolyfunctionalityScore(fit)

## Plot a heatmap of the mean gammas, to visualize differences in expression
## for each category
plot(fit)


## ----citation, echo=FALSE, results='asis'--------------------------------
cite_package("COMPASS", "flowWorkspace", "openCyto")


