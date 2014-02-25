data(COMPASS) ## loads the COMPASSContainer 'CC'
fit <- COMPASS(CC,
  category_filter=NULL,
  treatment=trt == "Treatment",
  control=trt == "Control",
  verbose=FALSE,
  iterations=100 ## set higher for a real analysis
)
