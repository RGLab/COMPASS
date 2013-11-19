library(COMPASS)
library(microbenchmark)

CC <- readRDS("refactor/COMPASS.rds")
data <- CC$data
combinations <- colnames(data[[1]])
CellCounts(data, combinations)
CellCounts(data, list("TNFa|IFNg"))
CellCounts(data, list("TNFa", "IFNg", "TNFa&IFNg"))
stopifnot( identical(
  unname(m1 <- CellCounts(data, combinations)),
  unname(m2 <- CellCounts(data, 1:6))
) )

microbenchmark( times=5,
  CellCounts(data, combinations),
  CellCounts(data, 1:6),
  CellCounts(data, "TNFa|IFNg&IL4")
)
