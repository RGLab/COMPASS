context("COMPASS Interface")

source("helper-init.R")

## these should run without problem
tmp <- COMPASS(CC, 
  treatment=trt == "Treatment",
  control=trt == "Control",
  category_filter=NULL,
  iterations=10,
  verbose=FALSE
)

expect_identical( tmp$fit$call$treatment, quote(trt == "Treatment") )

# trt_expr <- quote(trt == "Treatment")
# control_expr <- quote(trt == "Control")
# 
# tmp <- COMPASS(CC, 
#   treatment=trt_expr,
#   control=control_expr,
#   category_filter=NULL,
#   iterations=10,
#   verbose=FALSE
# )
# 
# expect_identical( tmp$fit$call$treatment, quote(trt == "Treatment") )
