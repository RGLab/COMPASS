library(COMPASS)

## simulate some data for COMPASS
set.seed(123)

make.data <- function(n, mean_ncells, sd_ncells, n_markers, 
  mean_intensity, sd_intensity, threshold) {
  
  replicate(n, simplify=FALSE, {
    nCells <- max(0, floor( rnorm(1, mean_ncells, sd_ncells) ))
    tmp <- matrix( 
      rnorm(nCells * n_markers, mean_intensity, sd_intensity), 
      ncol=n_markers 
    )
    tmp[ tmp < threshold ] <- 0
    colnames(tmp) <- paste0("V", 1:n_markers)
    return(tmp)
  })
  
}

n <- 50

unstim <- make.data(n, 50, 100, 6, 3000, 1500, 3500)
stim <- make.data(n, 200, 150, 6, 3500, 1500, 3500)
counts <- floor( rnorm(n*2, 5E4, 1E4) )
names(counts) <- paste0("Sample_", 1:length(counts))

data <- c(unstim, stim)
names(data) <- paste0("Sample_", 1:length(data))

meta <- data.frame(
  Stim=rep( c("Unstimulated", "Stimulated"), each=n ),
  Sample=paste0("Sample_", 1:(n*2)),
  Individual=rep( paste0("Individual_", 1:n), times=2 )
)

CC <- COMPASSContainer(
  data=data,
  meta=meta,
  counts=counts,
  individual_id="Individual",
  sample_id="Sample"
)

discrete <- COMPASS(
  data=CC,
  treatment=Stim == "Stimulated",
  control=Stim == "Unstimulated",
  iterations=1000
)

continuous <- COMPASS(
  data=CC,
  treatment=Stim == "Stimulated",
  control=Stim == "Unstimulated",
  iterations=1000,
  model="continuous"
)