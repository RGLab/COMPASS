context("CellCounts")

n <- 1E2

marginal_counts <- function(x) {
  t(sapply(x, function(xx) {
    apply(xx, 2, function(xxx) sum(xxx > 0))
  }))
}

data <- replicate(10, matrix( rnorm(n * 6, 2500, 500), nrow=n), simplify=FALSE)
data <- lapply(data, function(x) {
  colnames(x) <- LETTERS[1:6]
  x[ x < 1500 ] <- 0
  return (x)
})
combinations <- colnames(data[[1]]) ## [1] "A" "B" "C" "D" "E" "F"

test_that("The integer and character interfaces for CellCounts match up", {

  expect_identical(
    CellCounts(data, combinations),
    CellCounts(data, 1:6)
  )
  
  expect_identical(
    CellCounts(data, list(c(1, 2, 3))),
    CellCounts(data, list("A&B&C"))
  )
  
  y <- "A&B&C"
  expect_identical(
    CellCounts(data, list(c(1, 2, 3))),
    CellCounts(data, list(y))
  )
  
  expect_identical(
    CellCounts(data, 1:6),
    marginal_counts(data)
  )
  
})

test_that("We properly expand with a '*' in the name", {
  
  expect_identical(
    CellCounts(data, "A*B"),
    CellCounts(data, c("A&B", "!A&B", "A&!B", "!A&!B"))
  )
  
  ## TODO: allow combinations of * and [&|]
  expect_error( CellCounts(data, "A*B&C"))
  
  expect_identical(
    CellCounts(data, c("A*B", "C")),
    CellCounts(data, c("A&B", "!A&B", "A&!B", "!A&!B", "C"))
  )
  
})
