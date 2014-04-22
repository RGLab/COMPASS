CR <- COMPASS(CC,
  treatment=trt == "Treatment",
  control=trt == "Control",
  iterations=1000
)

save(CC, CR, file="data//COMPASS.rda", compress="xz", compression_level=9)
