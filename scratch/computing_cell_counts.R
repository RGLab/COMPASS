library(COMPASS)

CC_subs <- subset(CC, Stim == "92TH023 Env")
n_s <- CellCounts(CC_subs)

CC_subu <- subset(CC, Stim == "negctrl 1")
n_u <- CellCounts(CC_subu)

## marginal counts
markers <- colnames(CC$data[[1]])
m_s <- CellCounts(CC_subs, list( 2, 4, c(2, 4) ))
