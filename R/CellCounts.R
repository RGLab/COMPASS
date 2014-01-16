##' Compute Number of Cells Positive for Certain Cytokine Combinations
##' 
##' This function takes a list of cell expression matrices, each matrix \code{i} 
##' of dimension \code{N_i cells} by \code{K} common markers.
##' 
##' @param data A list of matrices. Each matrix \code{i} is
##'   of dimension \code{N_i cells} by \code{K} common markers.
##' @param combinations A list of 'combinations' == integer vectors,
##'   which is used to denote the subsets. See the examples.
##' @export
##' @seealso \code{\link{discrete_combinations}}
##' @examples
##' K <- 6 ## number of markers
##' data <- replicate(10, simplify=FALSE, {
##'   m <- matrix( rnorm(1E4 * K, 2000, 1000 ), ncol=K )
##'   m[m < 2500] <- 0
##'   colnames(m) <- c("IL2", "IL4", "IL6", "Mip1B", "IFNg", "TNFa")
##'   return(m)
##' })
##' names(data) <- sample(letters, 10)
##' 
##' combos <- list(1, 2, 3, 4, 5, 6) ## marginal cell counts
##' cc <- CellCounts(data, combos)
##' f <- function(data) {
##'   do.call(rbind, lapply(data, function(x) apply(x, 2, function(x) sum(x > 0))))
##' }
##' cc2 <- f(data)
##' stopifnot(identical( unname(cc), unname(cc2) ))
##' 
##' ## We can also generate cell counts by expressing various combinations
##' ## of markers (names) in the data.
##' CellCounts(data, "IL2|IL4") ## count cells expressing IL2 or IL4
##' CellCounts(data, "IL2|IL4|IL6") ## count cells expressing IL2, IL4 or IL6
##' CellCounts(data, c("IL2", "IL4", "IL6")) ## counts for each of IL2, IL4, IL6 (marginally)
CellCounts <- function(data, combinations) {
  UseMethod("CellCounts")
}

##' @S3method CellCounts COMPASSContainer
##' @method CellCounts COMPASSContainer
CellCounts.COMPASSContainer <- function(data, combinations) {
  data <- data$data
  NextMethod("CellCounts")
}

.CellCounts_character <- function(data, combinations) {
  
  output <- .Call("COMPASS_CellCounts_character", 
    data, 
    lapply(combinations, function(x) parse(text=x)),
    PACKAGE="COMPASS"
  )
  rownames(output) <- names(data)
  colnames(output) <- unlist(combinations)
  return(output)
}

.CellCounts_numeric <- function(data, combinations) {
  
  cn <- colnames(data[[1]])
  if (is.null(cn) || is.na(cn)) {
    warning("The column names of the matrices in your data are NA or NULL;",
      " they need to be set for the output to have sensible names.")
  }
  
  combinations <- lapply(combinations, function(combo) {
    if (is.character(combo)) {
      splat <- unlist( strsplit(combo, "&", fixed=TRUE) )
      return( sapply(splat, function(y) {
        if (substring(y, 1, 1) == "!") {
          return( match( substring(y, 2, nchar(y)), cn ) )
        } else {
          return( match(y, cn) )
        }
      }))
    } else {
      return(combo)
    }
  })
  
  names(combinations) <- sapply(combinations, function(x) {
    nm <- cn[ abs(x) ]
    nm[ x < 0 ] <- paste0("!", nm[x < 0])
    return( paste(nm[ order(abs(x)) ], collapse="&") )
  })
  
  return( .Call("COMPASS_CellCounts", as.list(data), lapply(combinations, as.integer), PACKAGE="COMPASS") )
}

##' @S3method CellCounts default
##' @method CellCounts default
CellCounts.default <- function(data, combinations) {
  
  if (missing(combinations)) {
    combinations <- BooleanSubsets.default(data)
  }
  
  if (length(unique(sapply(combinations, typeof))) > 1) {
    stop("'combinations' must all be of the same type")
  }
  
  output <- switch( t <- typeof(combinations[[1]]),
    double=.CellCounts_numeric(data, combinations),
    integer=.CellCounts_numeric(data, combinations),
    character=.CellCounts_character(data, combinations),
    stop ("Unexpected value type for combinations (type == ", t, ")")
  )
  
  return(output)

  
}
