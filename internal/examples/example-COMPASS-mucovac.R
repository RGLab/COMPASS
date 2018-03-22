library(MIMOSA)
library(pheatmap)

CD4 <- readRDS("data/Mucovac_CD4.rds")
unique(CD4$meta$Stim)
table(CD4$meta$Stim)

CD8 <- readRDS("data/Mucovac_CD8.rds")

treatment <- quote(Stim == "ENV-1-MUCO")
control <- quote(Stim == "negctrl")

CD8$meta<-within(CD8$meta,subject2<-factor(factor(PTID):factor(VISITNO)))
CD4$meta<-within(CD4$meta,subject2<-factor(factor(PTID):factor(VISITNO)))
CD8$individual_id<-"subject2"
CD4$individual_id<-"subject2"

#filter some duplicate(?) samples
filter<-na.omit(unique(CD4$meta[Stim%in%"ENV-1-MUCO",nrow(.SD),list(PTID,VISITNO)][V1>1,list(PTID)]))[[1]]
filter<-unique(c(filter,na.omit(CD4$meta[Stim%in%"negctrl",nrow(.SD),list(PTID,VISITNO)][V1>2,PTID])))

## debug
data <- CD4
x <- "ENV-1-MUCO"
treatment <- quote(Stim == x)
control <- quote(Stim == "negctrl")
category_filter <- function(x) colSums(x>4) > 2
verbose <- TRUE
subset <- quote( !PTID %in% filter & !is.na(PTID) )

#all visits
CD4_results <- COMPASS(
  data=CD4,
  treatment=Stim == x,
  control=Stim == "negctrl",
  category_filter=function(x) colSums(x>4) > 2,
  verbose=TRUE,
  iterations=10,
  filter_lowest_frequency=2,
  filter_specific_markers="IL4"
)

saveRDS(CD4_results, file="data/Mucovac_CD4_results.rds")

nc<-ncol(CD4_results$fit$gamma)
M<-apply(CD4_results$fit$gamma,1:2,mean)[,-nc]
rowann<-data.frame(do.call(rbind,strsplit(rownames(M),":")))
colnames(rowann)<-c("PTID","VISIT")
rownames(rowann)<-rownames(M)
cats<-CD4_results$fit$categories[-nc,]
cats<-data.frame(cats)
cats<-data.frame(sapply(cats,factor))[,1:(ncol(cats)-1)]
colnames(M)<-rownames(cats)
o<-order(rowann$VISIT)
pheatmap(M[o,],show_rownames=TRUE,row_annotation=rowann,cluster_rows=FALSE,cluster_cols=FALSE,cytokine_annotation=cats)

## debug
data <- CD8
treatment <- quote(Stim == "ENV-1-MUCO")
control <- quote(Stim == "negctrl")
verbose <- TRUE

filter<-na.omit(unique(CD8$meta[Stim%in%"ENV-1-MUCO",nrow(.SD),list(PTID,VISITNO)][V1>1,list(PTID)]))[[1]]
filter<-unique(c(filter,na.omit(CD8$meta[Stim%in%"negctrl",nrow(.SD),list(PTID,VISITNO)][V1>2,PTID])))

CD8_results <- COMPASS(
  data=CD8,
  treatment=Stim == "ENV-1-MUCO",
  control=Stim == "negctrl",
  verbose=TRUE,
  subset=!PTID%in%filter&!is.na(PTID),
  iterations=10
)
