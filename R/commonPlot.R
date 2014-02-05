##' Make a list of COMPASSResults comparable for display as heatamps
##' 
##' This function makes a list of COMPASSResults comparable for display as heatmaps.
##' It has some inefficiencies as it looks at all pairwise comparisons.
##' 
##' @param X An list of \code{COMPASSResult} objects
##' @param subset An \R expression, evaluated within the metadata, used to
##'   determine which individuals should be kept.
##' @param ... Optional arguments passed to \code{pheatmap}.
##' @importFrom clue solve_LSAP
##' @importFrom utils adist
##' @importFrom abind abind
##' @export
##' @examples
##'\dontrun{
##' #X is a list of COMPASSResult objects
##' X<-MakeComparable(X)
##' #These should show the same categories.
##' plot(X[[1]],remove_unexpressed_categories=FALSE)
##' plot(X[[2]],remove_unexpressed_categories=FALSE)
##'}
MakeComparable <- function(X,
                       subset, 
                        ...) {
  
  #check if the type of X is correct
  stopFlag<-FALSE
  if(!inherits(X,"list")){
    stopFlag<-TRUE
  }
  if(!Reduce(all,lapply(X,function(x)inherits(x,"COMPASSResult")))){
    stopFlag<-TRUE
  }
  if(stopFlag){
    stop("X must be a list of COMPASSResult")
  }
  
  subset_expr <- match.call()$subset
  .MakeComparableTwoCOMPASSResults<-function(x,y){
    nc_x <- ncol(x$fit$gamma)
    M_x <- x$fit$mean_gamma[, -nc_x]
    
    nc_y <- ncol(y$fit$gamma)
    M_y <- y$fit$mean_gamma[, -nc_y]
    
    ## make sure the row order is the same
    M_x <- M_x[ order(rownames(M_x)), , drop=FALSE ]
    M_y <- M_y[ order(rownames(M_y)), , drop=FALSE ]
    
    ## augment each matrix with the missing ptids from the other
    
    add_y <- setdiff( rownames(M_x), rownames(M_y) )
    add_x <- setdiff( rownames(M_y), rownames(M_x) )
    .augment<-function(mm,xx){
      if(length(xx)>0){
        rn<-c(rownames(mm),add_y)
        mm<-(rbind(mm,matrix(0,nrow=length(xx),ncol=ncol(mm))))
        rownames(mm)<-rn
      }
      return(mm)
    }
    M_y <- (.augment(M_y,add_y))
    M_x <- (.augment(M_x,add_x))
    
    #update the metadata as well
    meta_x <- data.table(x$data$meta)
    meta_y <- data.table(y$data$meta)
    coln_inter <- intersect(colnames(meta_x),colnames(meta_y))
    setkeyv(meta_x,c(x$data$individual_id))
    setkeyv(meta_y,c(x$data$individual_id))
    meta <- merge(meta_x,meta_y,all=TRUE)
    .consolidateColumns <- function(x){
      cn<-colnames(x)
      xcn<-cn[grepl("\\.x$",cn)]
      ycn<-cn[grepl("\\.y$",cn)]
      match_id<-as.vector(solve_LSAP(adist(xcn,ycn)))
      for(i in 1:length(xcn)){
        imp<-apply(data.frame(x[,get(xcn[i])],x[,get(ycn[match_id[i]])]),1,function(xx)na.omit(unique(xx)))
        if(class(imp)!="list"){
          x[,eval(xcn[i]):=imp]
          x[,eval(ycn[match_id[i]]):=NULL]
          setnames(x,xcn[i],gsub("\\.x$","",xcn)[i])
        }
      }
      #split in to x and y metadata
      xy_cols<-grepl("\\.x$|\\.y$",colnames(x))
      x_cols<-grepl("\\.x$",colnames(x))
      meta_x<-(x[,!xy_cols|x_cols,with=FALSE])
      setnames(meta_x,colnames(meta_x),gsub("\\.x$","",colnames(meta_x)))
      meta_y<-(x[,!xy_cols|!x_cols,with=FALSE])
      setnames(meta_y,colnames(meta_y),gsub("\\.y$","",colnames(meta_y)))
      return(list(meta_x,meta_y))
    }
    meta<-.consolidateColumns(meta)
    meta_x <- meta[[1]]
    meta_y <- meta[[2]]
    
    #order metadata columns consistently
    if(!all(colnames(meta_x)%in%colnames(meta_y))){
      stop("Internal Error: merging metadata columns failed")
    }
    meta_x <- meta_x[,colnames(meta_x)[order(colnames(meta_x))],with=FALSE]
    meta_y <- meta_y[,colnames(meta_x)[order(colnames(meta_x))],with=FALSE]
    
    .setrownames<-function(yy){
      rn<-(yy[,get(x$data$sample_id)])                        
      rn[is.na(rn)]<-"missing"
      rn<-make.unique(rn)
      yy<-data.frame(yy)
      rownames(yy)<-rn
      yy
    }
    meta_x <-.setrownames(meta_x)
    meta_y <-.setrownames(meta_y)
    
    #order rownames of M_x and M_y consistently
    orn<-rownames(M_x)[order(rownames(M_x))]
    M_x<-M_x[orn,]
    M_y<-M_y[orn,]

    
    #Ensure common categories..
    cx <- colnames(M_x)
    cy <- colnames(M_y)
    un<-union(cx,cy)
    cats<-unique(rbind(x$fit$categories,y$fit$categories))
    cats<-cats[order(rowSums(cats)/2,decreasing=FALSE),]
    cats<-cats[c(2:nrow(cats),1),]
    cats_str<-apply(cats[,-ncol(cats)],1,function(x)paste0(x,collapse=""))
    for(i in un[!un%in%cy]){
      M_y<-cbind(M_y,0)
      colnames(M_y)[ncol(M_y)]<-i
    }
    for(i in un[!un%in%cx]){
      M_x<-cbind(M_x,0)
      colnames(M_x)[ncol(M_x)]<-i
    }
    
    #return
    M_x<-(cbind(M_x,0))
    M_y<-cbind(M_y,0)
    colnames(M_x)[ncol(M_x)]<-cats_str[length(cats_str)]
    colnames(M_y)[ncol(M_y)]<-cats_str[length(cats_str)]
    
    M_x<-M_x[,cats_str]
    M_y<-M_y[,cats_str]
    
    nc <- ncol(M_x)
    ncg <- ncol(x$fit$gamma)
    ncd<-nc-ncg
    if(ncd>0){
      for(i in 1:ncd){
        x$fit$gamma<-(abind(x$fit$gamma,matrix(0,nrow=dim(x$fit$gamma)[1],ncol=dim(x$fit$gamma)[3]),along=2))
      }
    }
    nc <- ncol(M_y)
    ncg <- ncol(y$fit$gamma)
    ncd<-nc-ncg
    if(ncd>0){
      for(i in 1:ncd){
        y$fit$gamma<-(abind(y$fit$gamma,matrix(0,nrow=dim(y$fit$gamma)[1],ncol=dim(y$fit$gamma)[3]),along=2))
      }
    }
    
    colnames(M_x)[ncol(M_x)]<-colnames(x$fit$mean_gamma)[nc_x]
    colnames(M_y)[ncol(M_y)]<-colnames(y$fit$mean_gamma)[nc_y]
    x$fit$mean_gamma<-M_x
    y$fit$mean_gamma<-M_y
    x$data$meta<-meta_x
    y$data$meta<-meta_y
    x$fit$categories<-cats
    y$fit$categories<-cats
    return(list(x,y))
  }
  
  
  combo<-combn(length(X),2)
  for(i in 1:ncol(combo)){
    j<-combo[1,i]
    k<-combo[2,i]
    res<-.MakeComparableTwoCOMPASSResults(X[[j]],X[[k]])
    X[[j]]<-res[[1]]
    X[[k]]<-res[[2]]
  }
  return(X)

  
}
