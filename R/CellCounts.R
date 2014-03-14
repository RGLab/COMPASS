##' Compute Number of Cells Positive for Certain Cytokine Combinations
##' 
##' Compute the number of cells expressing a particular
##' combination of markers for each sample.
##' 
##' @param data Either a \code{COMPASSContainer}, or a list of matrices. 
##'   Each matrix \code{i} is of dimension \code{N_i} cells (rows) by 
##'   \code{K} common markers (columns).
##' @param combinations A list of 'combinations', used to denote the
##'   subsets of interest. See the examples for usage.
##' @export
##' @seealso \code{\link{Combinations}}
##' @examples
##' set.seed(123)
##' ## generate 10 simulated matrices of flow data
##' K <- 6 ## number of markers
##' data <- replicate(10, simplify=FALSE, {
##'   m <- matrix( rnorm(1E4 * K, 2000, 1000 ), ncol=K )
##'   m[m < 2500] <- 0
##'   colnames(m) <- c("IL2", "IL4", "IL6", "Mip1B", "IFNg", "TNFa")
##'   return(m)
##' })
##' names(data) <- sample(letters, 10)
##' head( data[[1]] )
##' 
##' ## generate counts over all available combinations of markers in data
##' str(CellCounts(data)) ## 64 columns, as all 2^6 combinations expressed
##' 
##' ## generate marginal counts
##' combos <- list(1, 2, 3, 4, 5, 6) ## marginal cell counts
##' cc <- CellCounts(data, combos)
##' 
##' ## a base R way of doing the same thing
##' f <- function(data) {
##'   do.call(rbind, lapply(data, function(x) apply(x, 2, function(x) sum(x > 0))))
##' }
##' cc2 <- f(data)
##' 
##' ## check that they're identical
##' stopifnot(identical( unname(cc), unname(cc2) ))
##' 
##' ## We can also generate cell counts by expressing various combinations
##' ## of markers (names) in the data.
##' 
##' ## count cells expressing IL2 or IL4
##' CellCounts(data, "IL2|IL4")
##' 
##' ## count cells expressing IL2, IL4 or IL6
##' CellCounts(data, "IL2|IL4|IL6")
##' 
##' ## counts for each of IL2, IL4, IL6 (marginally)
##' CellCounts(data, c("IL2", "IL4", "IL6"))
##' 
##' ## counts for cells that are IL2 positive and IL4 negative
##' CellCounts(data, "IL2 & !IL4")
##' 
##' ## expressing the same intent with indices
##' CellCounts(data, list(c(1, -2)))
##' 
##' ## all possible combinations
##' str(CellCounts(data, Combinations(6)))
##' 
##' ## can also call on COMPASSContainers
##' data(COMPASS)
##' CellCounts(CC, "M1&M2")
CellCounts <- function(data, combinations) {
  UseMethod("CellCounts")
}

##' @export
CellCounts.COMPASSContainer <- function(data, combinations) {
  data <- data$data
  NextMethod("CellCounts")
}

.CellCounts_character <- function(data, combinations) {
  
  output <- .Call(C_COMPASS_CellCounts_character, 
    data, 
    lapply(combinations, function(x) parse(text=x))
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
  
  return( .Call(C_COMPASS_CellCounts,
    as.list(data), 
    lapply(combinations, as.integer)
  ) )
}

##' @export
CellCounts.default <- function(data, combinations) {
  
  if (missing(combinations)) {
    combinations <- UniqueCombinations.default(data)
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
