context("test plotCOMPASSResultStack")

library(COMPASS)

generateCompassResultExample <- function(seed, k=6, markernames=NULL) {
  # k specifies number of markers
  set.seed(seed)
  n <- 100 ## number of samples

  sid_vec <- paste0("sid_", 1:n) ## sample ids; unique names used to denote samples
  iid_vec <- rep_len( paste0("iid_", 1:(n/10) ), n ) ## individual ids

  data <- replicate(n, {
    nrow <- round(runif(1) * 1E4 + 1000)
    ncol <- k
    vals <- rexp( nrow * ncol, runif(1, 1E-5, 1E-3) )
    vals[ vals < 2000 ] <- 0
    output <- matrix(vals, nrow, ncol)
    output <- output[ apply(output, 1, sum) > 0, ]
    colnames(output) <- if (is.null(markernames)) {
      paste0("M", 1:k)
    } else {
      markernames
    }
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
    individual_id="iid",
    sample_id="sid"
  )

  fit <- COMPASS( CC,
                  treatment=trt == "Treatment",
                  control=trt == "Control",
                  iterations=100
  )
  fit
}

###################################################

cr1 <- generateCompassResultExample(1)
cr2 <- generateCompassResultExample(2, markernames=paste0("N", 1:6))
cr3 <- generateCompassResultExample(2, k=8)
cr4 <- generateCompassResultExample(2, markernames=c(paste0("M", 4:6), paste0("M", 1:3)))
cr5 <- generateCompassResultExample(2)
cr6 <- generateCompassResultExample(3)

# Take 2 COMPASSResult objects with different markers and put them in a named list
compassResultsDiffMarkers <- list("Seed1" = cr1,
                                   "Seed2" = cr2)

# Take 2 COMPASSResult objects with overlapping markers and put them in a named list
# Seed1 has M1 though M6. Seed2 has M1 through M8.
compassResultsOverlappingMarkers <- list("Seed1" = cr1,
                                  "Seed2" = cr3)

# Same markers but different order
compassResultsSameMarkersDiffOrder <- list("Seed1" = cr1,
                                   "Seed2" = cr4)

test_that("only COMPASSResults with the same markers can be stacked", {
  # Different markers should fail
  expect_error(plotCOMPASSResultStack(compassResultsDiffMarkers,
                                      row_annotation = c("Seed"),
                                      variable = "Seed",
                                      show_rownames = TRUE,
                                      main = "Example Stacked Heatmap of Mean Probability of Response",
                                      fontsize = 14, fontsize_row = 13, fontsize_col = 11),
               "Marker names of all COMPASSResult objects must be the same.")

  # Overlapping markers should also fail
  expect_error(plotCOMPASSResultStack(compassResultsOverlappingMarkers,
                                      row_annotation = c("Seed"),
                                      variable = "Seed",
                                      show_rownames = TRUE,
                                      main = "Example Stacked Heatmap of Mean Probability of Response",
                                      fontsize = 14, fontsize_row = 13, fontsize_col = 11),
               "Marker names of all COMPASSResult objects must be the same.")

  # But same markers with a different order is fine
  expect_error(plotCOMPASSResultStack(compassResultsSameMarkersDiffOrder,
                                      row_annotation = c("Seed"),
                                      variable = "Seed",
                                      show_rownames = TRUE,
                                      main = "Example Stacked Heatmap of Mean Probability of Response",
                                      fontsize = 14, fontsize_row = 13, fontsize_col = 11),
               NA)
})

###################################################

# Test the getCatsAndSubsetNames function which is used by mergeMatricesForPlotCOMPASSResultStack
# Each line below gets turned into a column in the matrix
catsMatrixBefore <- matrix(c(c(1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0),
                       c(0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0),
                       c(0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0),
                       c(0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0),
                       c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0),
                       c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0),
                       c(1, 1, 2, 2, 3, 3, 4, 4, 5, 6, 0)),
                     ncol=7, nrow=11)
colnames(catsMatrixBefore) <- c("M4", "M5", "M6", "M1", "M2", "M3", "Counts")

catsMatrixAfter <- matrix(c(c(0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0),
                            c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0),
                            c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0),
                            c(1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0),
                            c(0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0),
                            c(0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0),
                            c(1, 1, 2, 2, 3, 3, 4, 4, 5, 6, 0)
                            ),
                          ncol=7, nrow=11)
colnames(catsMatrixAfter) <- c("M1", "M2", "M3", "M4", "M5", "M6", "Counts")
subsetNamesAfter <- c("!M1&!M2&!M3&M4&!M5&!M6", "!M1&!M2&!M3&!M4&M5&!M6", "!M1&!M2&!M3&!M4&M5&M6",
                      "M1&!M2&!M3&M4&!M5&!M6", "!M1&!M2&!M3&M4&M5&M6", "M1&!M2&!M3&M4&!M5&M6",
                      "!M1&M2&!M3&M4&M5&M6", "M1&M2&!M3&M4&!M5&M6", "!M1&M2&M3&M4&M5&M6",
                      "M1&M2&M3&M4&M5&M6", "!M1&!M2&!M3&!M4&!M5&!M6")

test_that("getCatsAndSubsetNames function orders catsMatrix columns and creates subset names correctly", {
  getCatsAndSubsetNamesOut <- COMPASS:::getCatsAndSubsetNames(catsMatrixBefore)
  expect_true(identical(getCatsAndSubsetNamesOut$cats, catsMatrixAfter))
  expect_true(identical(getCatsAndSubsetNamesOut$subsets, subsetNamesAfter))
})

###################################################

# This function is used to test whether the mean_gamma values
# from each unmerged COMPASSResult appear in the correct cells of MMerged
# It basically deconstructs MMerged and checks each column individually
checkMeanGammasMergedCorrectly <- function(crList, MMerged) {
  # Keep track of how many rows of MMerged we've checked so far
  checkedRows <- 0
  # We check the contents of MMerged are correct for each COMPASSResult
  ZeroMismatchesFound <- TRUE
  i <- 1
  while (ZeroMismatchesFound & i <= length(crList)) {
    cr <- crList[[i]]
    # The next nrow(cr$fit$mean_gamma) rows of MMerged should contain data from cr$fit$mean_gamma
    finalRow <- checkedRows + nrow(cr$fit$mean_gamma)
    # Grab the subset of MMerged we're concerned about right now
    MMergedSub <- MMerged[(checkedRows + 1):finalRow,]
    checkedRows <- finalRow

    # Corresponding columns of MMergedSub and cr$fit$mean_gamma don't necessarily have the same name,
    # the biggest reason being that MMergedSub column names contain sorted markers (e.g. "M1&M2")
    # whereas the cr$fit$mean_gamma column name markers are not necessarily sorted ("M2&M1" is fine).
    # # # #
    # Therefore we need to construct the *marker sorted* column names of cr$fit$mean_gamma
    # so that we can refer to the corresponding columns of MMergedSub for comparison.
    # We do this using getCatsAndSubsetNames(), defined in plotMeanGammaMulti.R
    # (Note: It is not a great idea to create tests for plotCOMPASSResultStack
    # using a function which is called by plotCOMPASSResultStack)
    catsData <- COMPASS:::getCatsAndSubsetNames(cr$fit$categories)
    # cr$fit$categories <- catsData$cats # Columns might be re-ordered, nothing else
    subsets <- catsData$subsets
    colnames(cr$fit$mean_gamma) <- subsets

    for (colname in colnames(MMergedSub)) {
      # Columns in the merged matrix which appear in cr$fit$mean_gamma
      # should match exactly
      if (colname %in% colnames(cr$fit$mean_gamma)) {
        # These two columns should contain the same contents
        if (all.equal(cr$fit$mean_gamma[, colname], MMergedSub[, colname], check.names = FALSE) != TRUE) {
          ZeroMismatchesFound <- FALSE
        }
      } else {
        # All other columns of the merged matrix should be filled with zeroes
        if (all.equal(rep(0, nrow(cr$fit$mean_gamma)), MmergedSub[, colname], check.names = FALSE) != TRUE) {
          ZeroMismatchesFound <- FALSE
        }
      }
    }
    i <- i + 1
  }
  # Return whether a mismatch was found
  ZeroMismatchesFound
}

test_that("mergeMatricesForPlotCOMPASSResultStack merges mean_gamma matrices correctly", {
  compassResults <- list("Seed1" = cr1,
                         "Seed2" = cr5,
                         "Seed3" = cr6)

  mergedMatrices <- mergeMatricesForPlotCOMPASSResultStack(compassResults,
                                                           row_annotation = c("Seed"),
                                                           variable = "Seed")

  expect_true(checkMeanGammasMergedCorrectly(crList=compassResults, MMerged=mergedMatrices$MMerged))
})
