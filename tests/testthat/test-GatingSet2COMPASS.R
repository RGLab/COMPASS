context("test COMPASSContainerFromGatingSet")
dataDir <- system.file("extdata/gs_manual",package="flowWorkspaceData")
gs <- flowWorkspace::load_gs(dataDir)
gs1 <- flowWorkspace::gs_clone(gs)
gs2 <- flowWorkspace::gs_clone(gs)
flowWorkspace::sampleNames(gs2) <- "sample2.fcs"
gs <- flowWorkspace::merge_list_to_gs(list(gs1, gs2))
cs <- flowWorkspace::gs_cyto_data(gs)
flowWorkspace::cs_lock(cs)#make sure compass code doesn't modify cs/cf object in place
pd <- flowWorkspace::pData(gs)
pd[["name"]] <- rownames(pd)
pd[["PTID"]] <- 1
flowWorkspace::pData(gs) <- pd
test_that("COMPASSContainerFromGatingSet", {
  cc <- COMPASSContainerFromGatingSet(gs, node = "CD8", mp = list("CD8/38+ DR+" = "DR"
                                                            ,"CD8/38+ DR-" = "38")
                                )
  mat <- cc[["data"]][[1]]
  expect_equal(colnames(mat), c("HLA-DR V500", "CD38 APC"))
  expect_equal(nrow(mat), 5785)
  expect_equal(cc[["counts"]][[1]], 14564)

  #wrong subset
  expect_error(cc <- COMPASSContainerFromGatingSet(gs, node = "CD8", mp = list("CD8/38+ DR+" = "DR"
                                                                  ,"CD4/38+ DR-" = "38")
                                      )
               , "not the children")
  #name column to be different from rownames
  flowWorkspace::pData(gs)[["name"]] <- "sample"
  cc <- COMPASSContainerFromGatingSet(gs, node = "CD8", mp = list("CD8/38+ DR+" = "DR"
                                                                  ,"CD8/38+ DR-" = "38")
  )
  mat <- cc[["data"]][[1]]
  expect_equal(colnames(mat), c("HLA-DR V500", "CD38 APC"))
  expect_equal(nrow(mat), 5785)
  expect_equal(cc[["counts"]][[1]], 14564)
})
