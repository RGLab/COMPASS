context("FS, PFS")

test_that("We compute the functionality, polyfunctionality scores correctly for matrices", {

  m <- matrix(1, 5, 8)

  ## Remove null category
  m <- m[, -ncol(m)]

  expect_equal(
    FunctionalityScore(m, 3),
    rep(1, 5)
  )

  degree <- c(1, 1, 1, 2, 2, 2, 3)

  expect_equal(
    PolyfunctionalityScore(m, degree, n = 3),
    rep(1, 5)
  )

  m <- matrix( c(1, 1, 1, 0, 0, 0, 0), 5, 7, byrow=TRUE )
  expect_equal(
    FunctionalityScore(m, n=3),
    rep(3/7, 5)
  )

  expect_equal(
    PolyfunctionalityScore(m, degree, n = 3),
    rep(1/6, 5)
  )

})

test_that("We compute the functionality, polyfunctionality scores correctly for COMPASS results", {

  CRC <- CR
  CRC$fit$mean_gamma <- matrix(1, 5, 8)

  expect_equal(
    FunctionalityScore(CRC, n = 3),
    rep(1, 5)
  )

  expect_equal(
    PolyfunctionalityScore(CRC),
    rep(1, 5)
  )

  CRC$fit$mean_gamma <- matrix( c(1, 1, 1, 0, 0, 0, 0, 0), 5, 8, byrow=TRUE )

  expect_equal(
    FunctionalityScore(CRC),
    rep(3/7, 5)
  )

  expect_equal(
    PolyfunctionalityScore(CRC),
    rep(1/6, 5)
  )

})
