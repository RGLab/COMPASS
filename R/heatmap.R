## plot heatmap
library(pheatmap)
library(RColorBrewer)

   anno_map<-data.frame(Status = anno$status, row.names=as.numeric(anno$pd)) #  Status: placebo/infected/non-infected # rownames=PTID
   anno_map <-  data.frame(anno_map[order(match(anno$pd, PTID)),], row.names = 1:I)
   colnames(anno_map) = "Status"
   MgammaK1 <- Mgamma[,1:K1] # first K1 categories, excluding Null
   wh0 <- which(colSums(MgammaK1) >0)
   row_ann <- data.frame(d[wh0,1:M])
   colnames(row_ann) <- Cnames
   rownames(row_ann) <- Cate[wh0]
   row_ann <- data.frame(apply(row_ann,2,function(x) as.factor(x)))
   heatmat <- MgammaK1[,wh0] # matrix to be plotted
   colnames(heatmat) <- Cate[wh0] # Categories
   rownames(heatmat) <- 1:I
   
   var1 <- brewer.pal(n=3, "Dark2") # n = number of levels of status
   names(var1) <- levels(anno_map$Status)
   anno_col <- list(Status =var1)
  
   heatmat<-heatmat[order(anno_map[,1]),]

   pal <- brewer.pal(n=9, "Blues")
   pheatmap(heatmat, cluster_cols=FALSE,  cluster_row =FALSE,cytokine_annotation=row_ann, treeheight_row=0, treeheight_col=0, 
            color=colorpanel(length(pal), low=pal[1], high=pal[9]),row_annotation_colors=anno_col, show_rownames=F, row_annotation = anno_map,
            border_color = NA)
