library(COMPASS)
if (file.exists("refactor/COMPASS.rds")) {
  CC <- readRDS("refactor/COMPASS.rds")
  
  data <- CC
  treatment <- quote(Stim == "92TH023 Env")
  control <- quote(Stim == "negctrl 1")
  model <- "discrete"
  verbose <- TRUE
  
  ## An example call
  set.seed(123)
  gctorture(TRUE)
  discrete1 <- COMPASS(
    data=CC,
    treatment="92TH023 Env",
    control="negctrl 1",
    model="continuous",
    iterations=10,
    replications=8
  )
  gctorture(FALSE)
  
}
