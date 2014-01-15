
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
  vals <- rexp( nrow * ncol, runif(1, 1E-5, 1E-3) )
  vals[ vals < 2000 ] <- 0
  output <- matrix(vals, nrow, ncol)
  output <- output[ apply(output, 1, sum) > 0, ]
  colnames(output) <- paste0("M", 1:k)
  return(output)
})
names(data) <- sid_vec
head( data[[1]] )


## ----sim-counts----------------------------------------------------------
counts <- sapply(data, nrow) + round( rnorm(n, 1E4, 1E3) )
counts <- setNames( as.integer(counts), names(counts) )
head(counts)


## ----sim-meta------------------------------------------------------------
meta <- data.frame(
  sid=sid_vec,
  iid=iid_vec,
  trt=sample( c("Control", "Treatment"), n, TRUE )
)
head(meta)


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

## Obtain the posterior difference, posterior log ratio from a COMPASSResult
post <- Posterior(fit)

## Plot a heatmap of the mean gammas, to visualize differences in expression
## for each category
plot(fit)

## Visualize the posterior difference, log difference with a heatmap
plot(fit, measure=PosteriorDiff(fit))
plot(fit, measure=PosteriorLogDiff(fit))



## ----citation, echo=FALSE, results='asis'--------------------------------
citations <- cite_package("COMPASS", "flowWorkspace", "openCyto", "base")
invisible(lapply(citations, function(x) cat(x, "\n\n")))


