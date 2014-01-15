lo = function(rown, coln, nrow, ncol, cellheight = NA, cellwidth = NA, treeheight_col, treeheight_row, legend, annotation, annotation_colors, annotation_legend, main, fontsize, fontsize_row, fontsize_col,row_annotation,row_annotation_legend,row_annotation_colors, cytokine_annotation, polar=FALSE,...){
  annot_legend_width = unit(0, "npc")
  #cytokine labels
  if(!is.na(cytokine_annotation[[1]][1])){
    cytn<-colnames(cytokine_annotation)
    longest_cytn<-which.max(strwidth(cytn,units='in'))
    gp = list(fontsize = fontsize_row, ...)
    cytn_width = unit(1, "grobwidth", textGrob(cytn[longest_cytn], gp = do.call(gpar, gp))) + unit(10, "bigpts") 
  }else{
    cytn_width = unit(0,"bigpts")
  }
  
  # Get height of colnames and length of rownames
  if(!is.null(coln[1])){
    if(!is.null(row_annotation)[[1]][1]){
      coln<-c(coln,colnames(row_annotation))
      longest_coln = which.max(strwidth(coln, units = 'in'))
      gp = list(fontsize = fontsize_col, ...)
      coln_height = unit(1, "grobheight", textGrob(coln[longest_coln], rot = 90, gp = do.call(gpar, gp))) + unit(5, "bigpts")
    }else{
      longest_coln = which.max(strwidth(coln, units = 'in'))
      gp = list(fontsize = fontsize_col, ...)
      coln_height = unit(1, "grobheight", textGrob(coln[longest_coln], rot = 90, gp = do.call(gpar, gp))) + unit(5, "bigpts")
    }
  }
  else{
    coln_height = unit(5, "bigpts")
  }
  
  if(!is.null(rown[1])){
    longest_rown = which.max(strwidth(rown, units = 'in'))
    gp = list(fontsize = fontsize_row, ...)
    rown_width = unit(1, "grobwidth", textGrob(rown[longest_rown], gp = do.call(gpar, gp))) + unit(10, "bigpts")
  }
  else{
    rown_width = unit(5, "bigpts")
  }
  
  gp = list(fontsize = fontsize, ...)
  # Legend position
  if(!is.na(legend[1])){
    longest_break = which.max(nchar(names(legend)))
    longest_break = unit(1.1, "grobwidth", textGrob(as.character(names(legend))[longest_break], gp = do.call(gpar, gp)))
    title_length = unit(1.1, "grobwidth", textGrob("Scale", gp = gpar(fontface = "bold", ...)))
    legend_width = unit(12, "bigpts") + longest_break * 1.2
    legend_width = max(title_length, legend_width)
    if(polar){
      legend_width<-unit(1.5,"inches")
    }
  }
  else{
    legend_width = unit(0, "bigpts")
  }
  
  # Set main title height
  if(is.na(main)){
    main_height = unit(0, "npc")
  }
  else{
    main_height = unit(1.5, "grobheight", textGrob(main, gp = gpar(fontsize = 1.3 * fontsize, ...)))
  }
  
  # Column annotations
  if(!is.na(annotation[[1]][1])){
    # Column annotation height 
    annot_height = unit(ncol(annotation) * (8 + 2) + 2, "bigpts")
    
    # Width of the correponding legend
    longest_ann = which.max(nchar(as.matrix(annotation)))
    annot_legend_width = unit(1.2, "grobwidth", textGrob(as.matrix(annotation)[longest_ann], gp = gpar(...))) + unit(12, "bigpts")
    if(!annotation_legend&!row_annotation_legend){
      annot_legend_width = unit(0, "npc")
    }
  }
  else{
    annot_height = unit(0, "bigpts")
    annot_legend_width = unit(0, "bigpts")+annot_legend_width
  }
  
  # Row annotations
  if(!is.na(row_annotation[[1]][1])){
    #width of the annotation beside the rows
    row_annotation_width = unit(ncol(row_annotation) * (8 + 2) + 2,"bigpts")
    
    #width of the legend
    
    ## name nchars
    max_nm_nchars <- max(nchar( names(row_annotation) ))
    max_nm_annot <- max(nchar( as.matrix( row_annotation ) ))
    if (max_nm_nchars >= max_nm_annot) {
      annot_legend_width <- unit(1.2, "grobwidth",
        textGrob( names(row_annotation)[ which.max( nchar( names(row_annotation) ) ) ],
          gp=gpar(...))) + unit(12,"bigpts")+annot_legend_width
      if(!row_annotation_legend&!annotation_legend){
        annot_legend_width = unit(0,"npc")
      }
    } else {
      longest_row_annotation = which.max(nchar(as.matrix(row_annotation)))
      annot_legend_width = unit(1.2, "grobwidth",textGrob(as.matrix(row_annotation)[longest_row_annotation],gp=gpar(...))) + unit(12,"bigpts")+annot_legend_width
      if(!row_annotation_legend&!annotation_legend){
        annot_legend_width = unit(0,"npc")
      }
    }
  }else{
    row_annotation_width = unit(0,"bigpts")
    annot_legend_width= unit(0,"bigpts")+annot_legend_width #row legend appears in the same viewport as column legend
  }
  
  # cytokine annotations
  if(!is.na(cytokine_annotation[[1]][1])){
    #height of the annotation below the matrix
    cytokine_height = unit(ncol(cytokine_annotation) * (8 + 2) + 2,"bigpts")
    #no legend
  }else{
    cytokine_height = unit(0,"bigpts")
  }
  
  # Tree height
  treeheight_col = unit(treeheight_col, "bigpts") + unit(5, "bigpts")
  treeheight_row = unit(treeheight_row, "bigpts") + unit(5, "bigpts") + cytn_width
  
  # Set cell sizes
  if(is.na(cellwidth)){
    matwidth = unit(1, "npc") - rown_width - legend_width - treeheight_row - annot_legend_width - row_annotation_width
  }
  else{
    matwidth = unit(cellwidth * ncol, "bigpts")
  }
  
  if(is.na(cellheight)){
    matheight = unit(1, "npc") - main_height - coln_height - treeheight_col - annot_height - cytokine_height
  }
  else{
    matheight = unit(cellheight * nrow, "bigpts")
  }  
  
  
  # Produce layout()
  pushViewport(viewport(layout = 
      grid.layout(nrow = 6, ncol = 6, 
        widths = unit.c(treeheight_row, matwidth, row_annotation_width, rown_width, legend_width, annot_legend_width), 
        heights = unit.c(main_height, treeheight_col, annot_height, matheight, coln_height, cytokine_height)), gp = do.call(gpar, gp)))
  
  # Get cell dimensions
  pushViewport(vplayout(4, 2))
  cellwidth = convertWidth(unit(0:1, "npc"), "bigpts", valueOnly = T)[2] / ncol
  cellheight = convertHeight(unit(0:1, "npc"), "bigpts", valueOnly = T)[2] / nrow
  upViewport()
  
  # Return minimal cell dimension in bigpts to decide if borders are drawn
  mindim = min(cellwidth, cellheight) 
  return(mindim)
}

draw_dendrogram = function(hc, horizontal = T){
  h = hc$height / max(hc$height) / 1.05
  m = hc$merge
  o = hc$order
  n = length(o)
  
  m[m > 0] = n + m[m > 0] 
  m[m < 0] = abs(m[m < 0])
  
  dist = matrix(0, nrow = 2 * n - 1, ncol = 2, dimnames = list(NULL, c("x", "y"))) 
  dist[1:n, 1] = 1 / n / 2 + (1 / n) * (match(1:n, o) - 1)
  
  for(i in 1:nrow(m)){
    dist[n + i, 1] = (dist[m[i, 1], 1] + dist[m[i, 2], 1]) / 2
    dist[n + i, 2] = h[i]
  }
  
  draw_connection = function(x1, x2, y1, y2, y){
    grid.lines(x = c(x1, x1), y = c(y1, y))
    grid.lines(x = c(x2, x2), y = c(y2, y))
    grid.lines(x = c(x1, x2), y = c(y, y))
  }
  
  if(horizontal){
    for(i in 1:nrow(m)){
      draw_connection(dist[m[i, 1], 1], dist[m[i, 2], 1], dist[m[i, 1], 2], dist[m[i, 2], 2], h[i])
    }
  }
  
  else{
    gr = rectGrob()
    pushViewport(viewport(height = unit(1, "grobwidth", gr), width = unit(1, "grobheight", gr), angle = 90))
    dist[, 1] = 1 - dist[, 1] 
    for(i in 1:nrow(m)){
      draw_connection(dist[m[i, 1], 1], dist[m[i, 2], 1], dist[m[i, 1], 2], dist[m[i, 2], 2], h[i])
    }
    upViewport()
  }
}

draw_headerplot = function(order,data,ylabel){
  gr<-rectGrob(y=0.4,height=0.9)
  pushViewport(viewport(height = unit(1, "grobheight", gr), width = unit(1, "grobwidth", gr)))
  pushViewport(dataViewport(1:length(data),data,extension=c(1/length(data)/2,0.05)))
  grid.rect()
  grid.points(x=1:length(data),y=data[order],gp=gpar(cex=0.25,col="red",lwd=3))
  grid.abline(slope=0,gp=gpar(lty=2,lwd=2))
  grid.yaxis()
  gr<-textGrob(ylabel,rot=90)
  grid.text(ylabel,x=unit(-5,units="grobwidth",gr),y=0.5,rot=90)
  popViewport(2)
}

draw_matrix = function(matrix, border_color, fmat, fontsize_number){
  if (!is.matrix(matrix)) {
    matrix <- as.matrix(matrix)
  }
  n = nrow(matrix)
  m = ncol(matrix)
  x = (1:m)/m - 1/2/m
  y = 1 - ((1:n)/n - 1/2/n)
  for(i in 1:m){
    grid.rect(x = x[i], y = y[1:n], width = 1/m, height = 1/n, gp = gpar(fill = matrix[,i], col = border_color))
    if(attr(fmat, "draw")){
      grid.text(x = x[i], y = y[1:n], label = fmat[, i], gp = gpar(col = "grey30", fontsize = fontsize_number))
    }
  }
}

draw_colnames = function(coln, ...){
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = 0.5, hjust = 0, rot = 270, gp = gpar(...))
}

draw_rownames = function(rown, ...){
  n = length(rown)
  y = 1 - ((1:n)/n - 1/2/n)
  grid.text(rown, x = unit(0.04, "npc"), y = y, vjust = 0.5, hjust = 0, gp = gpar(...))	
}

.toCart<-function(r,theta,C=0.5){
  list(x=r*cos(theta)+C,y=r*sin(theta)+C)
}

##' @importFrom plyr ddply
##' @importFrom plyr .
.polarLegend<-function(R=0.25,N=11,C=0.5){
 
  
  polar.grid<-expand.grid(r=seq(0,R,l=N),theta=seq(0,pi,l=N))
  cart.grid<-do.call(cbind,.toCart(polar.grid$r,polar.grid$theta,C=0))/R
  saturation<-ddply(cbind(cart.grid,polar.grid),.(r),function(x)within(x,saturation<-r/R))
  colors<-ddply(saturation,.(theta),function(x)within(x,cbind(A<-theta/pi,B<-1-A)))
  head(colors)
  
  
  pal<-rgb2hsv(colors$A,0,colors$B)
  pal["v",]<-1
  .scale<-function(x){x/max(x)}
  cols<-hsv(pal[1,],colors$saturation,pal[3,],alpha=0.2)
  # plot(cart.grid,col=cols,pch=20,cex=4)
  return(list(colors=cols,position=cart.grid))
}

draw_legend = function(color, breaks, legend, ...){
  height = min(unit(1, "npc"), unit(150, "bigpts"))
  pushViewport(viewport(x = 0, y = unit(1, "npc"), just = c(0, 1), height = height))
  legend_pos = (legend - min(breaks)) / (max(breaks) - min(breaks))
  breaks = (breaks - min(breaks)) / (max(breaks) - min(breaks))
  h = breaks[-1] - breaks[-length(breaks)]
  grid.rect(x = 0, y = breaks[-length(breaks)], width = unit(10, "bigpts"), height = h, hjust = 0, vjust = 0, gp = gpar(fill = color, col = "#FFFFFF00"))
  grid.text(names(legend), x = unit(12, "bigpts"), y = legend_pos, hjust = 0, gp = gpar(...))
  upViewport()
}

convert_cytokine_annotations = function(annotation){
  ck_degree<-apply(annotation,1,function(x)sum(as.numeric(as.character(x))))
  mdegree<-max(ck_degree)
  udegrees<-unique(setdiff(ck_degree,0))
  new_annotation<-matrix("#FFFFFF",ncol=ncol(annotation),nrow=nrow(annotation))
  pal<-brewer.pal("Paired",n=mdegree)
  for(i in udegrees){
    ind<-which(annotation==1,TRUE)[,1]%in%which(ck_degree%in%i)
    for(j in which(ind)){
      new_annotation[which(annotation==1,TRUE)[j,1],which(annotation==1,TRUE)[j,2]]<-pal[i]
    }
  }
  colnames(new_annotation)<-colnames(annotation)
  rownames(new_annotation)<-rownames(annotation)
  new_annotation
}

convert_annotations = function(annotation, annotation_colors){
  new = annotation
  for(i in 1:ncol(annotation)){
    a = annotation[, i]
    b = annotation_colors[[colnames(annotation)[i]]]
    if(is.character(a) | is.factor(a)){
      a = as.character(a)
      if(length(setdiff(a, names(b))) > 0){
        stop(sprintf("Factor levels on variable %s do not match with annotation_colors", colnames(annotation)[i]))
      }
      new[, i] = b[a]
    }
    else{
      a = cut(a, breaks = 100)
      new[, i] = colorRampPalette(b)(100)[a]
    }
  }
  return(as.matrix(new))
}

draw_annotations = function(converted_annotations, border_color){
  n = ncol(converted_annotations)
  m = nrow(converted_annotations)
  x = (1:m)/m - 1/2/m
  y = cumsum(rep(8, n)) - 4 + cumsum(rep(2, n))
  for(i in 1:m){
    grid.rect(x = x[i], unit(y[1:n], "bigpts"), width = 1/m, height = unit(8, "bigpts"), gp = gpar(fill = converted_annotations[i, ], col = border_color))
  }
}

draw_row_annotations = function(converted_annotations, border_color){
  n = ncol(converted_annotations)
  m = nrow(converted_annotations)
  y = rev((1:m)/m - 1/2/m)
  x = cumsum(rep(8, n)) - 4 + cumsum(rep(2, n))
  for(i in 1:m){
    grid.rect(y = y[i], unit(x[1:n], "bigpts"), height = 1/m, width = unit(8, "bigpts"), gp = gpar(fill = converted_annotations[i, ], col = border_color))
  }
}

draw_annotation_legend = function(annotation, annotation_colors, border_color, yoff=c(0,0),...){
  y = unit(1, "npc")
  text_height = unit(1, "grobheight", textGrob("FGH", gp = gpar(...)))
  y<-y-1.5*text_height*(yoff[1]+2*yoff[2])
  for(i in names(annotation_colors)){
    grid.text(i, x = 0, y = y, vjust = 1, hjust = 0, gp = gpar(fontface = "bold", ...))
    y = y - 1.5 * text_height
    if(is.character(annotation[, i]) | is.factor(annotation[, i])){
      for(j in 1:length(annotation_colors[[i]])){
        grid.rect(x = unit(0, "npc"), y = y, hjust = 0, vjust = 1, height = text_height, width = text_height, gp = gpar(col = border_color, fill = annotation_colors[[i]][j]))
        grid.text(names(annotation_colors[[i]])[j], x = text_height * 1.3, y = y, hjust = 0, vjust = 1, gp = gpar(...))
        y = y - 1.5 * text_height
      }
    }
    else{
      yy = y - 4 * text_height + seq(0, 1, 0.02) * 4 * text_height
      h = 4 * text_height * 0.02
      grid.rect(x = unit(0, "npc"), y = yy, hjust = 0, vjust = 1, height = h, width = text_height, gp = gpar(col = "#FFFFFF00", fill = colorRampPalette(annotation_colors[[i]])(50)))
      txt = rev(range(grid.pretty(range(annotation[, i], na.rm = TRUE))))
      yy = y - c(0, 3) * text_height
      grid.text(txt, x = text_height * 1.3, y = yy, hjust = 0, vjust = 1, gp = gpar(...))
      y = y - 4.5 * text_height
    }
    y = y - 1.5 * text_height
  }
}

draw_main = function(text, ...){
  grid.text(text, gp = gpar(fontface = "bold", ...))
}

vplayout = function(x, y){
  return(viewport(layout.pos.row = x, layout.pos.col = y))
}

heatmap_motor = function(matrix, border_color, cellwidth, cellheight, tree_col, tree_row, treeheight_col, treeheight_row, filename, width, height, breaks, color, legend, annotation, annotation_colors, annotation_legend, main, fontsize, fontsize_row, fontsize_col, fmat, fontsize_number, row_annotation, row_annotation_legend, row_annotation_colors, cytokine_annotation, headerplot,polar=polar,...){
  grid.newpage()
  
  
  # Set layout
  mindim = lo(coln = colnames(matrix), rown = rownames(matrix), nrow = nrow(matrix), ncol = ncol(matrix), cellwidth = cellwidth, cellheight = cellheight, treeheight_col = treeheight_col, treeheight_row = treeheight_row, legend = legend, annotation = annotation, annotation_colors = annotation_colors, annotation_legend = annotation_legend, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, row_annotation = row_annotation, row_annotation_legend = row_annotation_legend, row_annotation_colors = row_annotation_colors, cytokine_annotation = cytokine_annotation,headerplot=headerplot,polar=polar,...)
  
  if(!is.na(filename)){
    pushViewport(vplayout(1:5, 1:6)) 
    
    if(is.na(height)){
      height = convertHeight(unit(0:1, "npc"), "inches", valueOnly = T)[2]
    }
    if(is.na(width)){
      width = convertWidth(unit(0:1, "npc"), "inches", valueOnly = T)[2]
    }
    
    # Get file type
    r = regexpr("\\.[a-zA-Z]*$", filename)
    if(r == -1) stop("Improper filename")
    ending = substr(filename, r + 1, r + attr(r, "match.length"))
    
    f = switch(ending,
      pdf = function(x, ...) pdf(x, ...),
      png = function(x, ...) png(x, units = "in", res = 300, ...),
      jpeg = function(x, ...) jpeg(x, units = "in", res = 300, ...),
      jpg = function(x, ...) jpeg(x, units = "in", res = 300, ...),
      tiff = function(x, ...) tiff(x, units = "in", res = 300, compression = "lzw", ...),
      bmp = function(x, ...) bmp(x, units = "in", res = 300, ...),
      stop("File type should be: pdf, png, bmp, jpg, tiff")
    )
    
    # print(sprintf("height:%f width:%f", height, width))
    f(filename, height = height, width = width)
    heatmap_motor(matrix, cellwidth = cellwidth, cellheight = cellheight, border_color = border_color, tree_col = tree_col, tree_row = tree_row, treeheight_col = treeheight_col, treeheight_row = treeheight_row, breaks = breaks, color = color, legend = legend, annotation = annotation, annotation_colors = annotation_colors, annotation_legend = annotation_legend, filename = NA, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, fmat = fmat, fontsize_number =  fontsize_number, row_annotation = row_annotation, row_annotation_legend = row_annotation_legend, cytokine_annotation = cytokine_annotation, ...)
    dev.off()
    upViewport()
    return()
  }
  
  # Omit border color if cell size is too small 
  if(mindim < 3) border_color = NA
  
  # Draw title
  if(!is.na(main)){
    pushViewport(vplayout(1, 2))
    draw_main(main, fontsize = 1.3 * fontsize, ...)
    upViewport()
  }
  
  # Draw tree for the columns
  if(!is.na(tree_col[[1]][1]) & treeheight_col != 0){
    pushViewport(vplayout(2, 2))
    draw_dendrogram(tree_col, horizontal = T)
    upViewport()
  }
  
  #draw headerplot
  if(!is.na(headerplot) & treeheight_col != 0){
    pushViewport(vplayout(2 , 2))
    draw_headerplot(headerplot$order,headerplot$data,headerplot$ylabel)
    upViewport()
  }
  
  # Draw tree for the rows
  if(!is.na(tree_row[[1]][1]) & treeheight_row != 0){
    pushViewport(vplayout(4, 1))
    draw_dendrogram(tree_row, horizontal = F)
    upViewport()
  }
  
  # Draw matrix
  pushViewport(vplayout(4, 2))
  draw_matrix(matrix, border_color, fmat, fontsize_number)
  upViewport()
  
  # Draw colnames
  if(length(colnames(matrix)) != 0){
    pushViewport(vplayout(5, 2))
    pars = list(colnames(matrix), fontsize = fontsize_col, ...)
    do.call(draw_colnames, pars)
    upViewport()
  }
  
  # Draw rownames
  if(length(rownames(matrix)) != 0){
    pushViewport(vplayout(4, 4)) 
    pars = list(rownames(matrix), fontsize = fontsize_row, ...)
    do.call(draw_rownames, pars)
    upViewport()
  }
  
  # Draw annotation tracks
  if(!is.na(annotation[[1]][1])){
    pushViewport(vplayout(3, 2))
    converted_annotation = convert_annotations(annotation, annotation_colors)
    draw_annotations(converted_annotation, border_color)
    upViewport()
  }
  
  #Draw row annotation tracks
  if(!is.na(row_annotation[[1]][1])){
    pushViewport(vplayout(4,3))
    converted_row_annotations = convert_annotations(row_annotation, row_annotation_colors)
    draw_row_annotations(converted_row_annotations, border_color)
    upViewport()
    #label the rows
    pushViewport(vplayout(5,3))
    pars_row_annotations = list(colnames(converted_row_annotations), fontsize = fontsize_col, ...)
    do.call(draw_colnames, pars_row_annotations)
    upViewport()
  }
  
  #Draw cytokine annotation tracks
  if(!is.na(cytokine_annotation[[1]][1])){
    pushViewport(vplayout(6,2))
    
    converted_cytokine_annotations = convert_cytokine_annotations(cytokine_annotation)
    draw_annotations(converted_cytokine_annotations, border_color)
    upViewport()
    
    #label the rows
    pushViewport(vplayout(6,1))
    pars_cytokine_annotations = list(rev(colnames(converted_cytokine_annotations)), fontsize = fontsize_col, ...)
    do.call(draw_rownames, pars_cytokine_annotations)
    upViewport()
  }
  
  # Draw annotation legend
  if(!is.na(annotation[[1]][1]) & annotation_legend){
    if(length(rownames(matrix)) != 0){
      pushViewport(vplayout(4:5, 6)) 
    }
    else{
      pushViewport(vplayout(3:5, 6)) 
    }
    draw_annotation_legend(annotation, annotation_colors, border_color, fontsize = fontsize, ...)
    upViewport()
  }
  
  #  Draw row_annotation legend
  if(!is.na(row_annotation[[1]][1]) & row_annotation_legend){
    if(length(rownames(matrix)) != 0){
      pushViewport(vplayout(4:5, 6)) 
    }
    else{
      pushViewport(vplayout(3:5, 6)) 
    }
    if(!is.na(annotation[[1]][1]) & annotation_legend){
      yoff<-c(sum(unlist(lapply(annotation,function(x)length(unique(x))))),length(annotation))
    }else{
      yoff<-c(0,0)
    }
    draw_annotation_legend(row_annotation, row_annotation_colors, border_color, fontsize = fontsize, yoff=yoff, ...)
    upViewport()
  }
  
  # Draw legend
  if(!is.na(legend[1])){
    length(colnames(matrix))
    if(length(rownames(matrix)) != 0){
      pushViewport(vplayout(4:5, 5)) 
    }
    else{
      pushViewport(vplayout(3:5, 5)) 
    }
    if(!(polar)){
      draw_legend(color, breaks, legend, fontsize = fontsize, ...)
    }else{
      draw_polar_legend(fontsize=fontsize)
    }
    upViewport()
  }
  
  
}
#Todo pass stim condition labels
draw_polar_legend <- function(fontsize=NA,treatmentLabel=c("Condition X","Condition Y")){
  #if we could get the size of the window or root viewport, we could set the size of the legend so
  #that it is scaled to look like a circle in any aspect ratio.
  height =unit(0.25,"npc")
  width = unit(1,"npc")
  pushViewport(viewport(x = unit(0,"npc"), y = unit(0.8, "npc"), just = c(0, 1), width=width,height = height))
  C=0.45 #center
  R=0.25 #radius
  N=51   #number of points
  .aspect<-function(){
    xy<-.toCart(R,seq(0,pi,l=100),C=C)  
    poly<-polygonGrob(x=xy$x, y=xy$y, gp=gpar(fill=NA, col="black",lwd=2),default.units="npc") # outer shell 
    gw<-convertUnit(grobWidth(poly),"cm",valueOnly=TRUE)
    gh<-convertUnit(grobHeight(poly),"cm",valueOnly=TRUE)
    asp.ratio<-gh/gw
  }
  asp.scale<-0.5/.aspect()
  R<-R*asp.scale
  xy<-.toCart(R,seq(0,pi,l=100),C=C)  
  polar_legend <- .polarLegend(R=R,N=N,C=C)
  
  o<-order(polar_legend$position[,"x"],polar_legend$position[,"y"],decreasing=TRUE)
  #for(i in (1:nrow(polar_legend$position))[o]){
    pts<-pointsGrob(x=polar_legend$position[,"x"]*R+C,y=polar_legend$position[,"y"]*R+C,gp=gpar(col=polar_legend$color,cex=0.25),default.units="npc")
  #}
  #for(i in (1:nrow(polar_legend$position))[o]){
  #  grid.points(x=polar_legend$position[i,"x"]*R+C,y=polar_legend$position[i,"y"]*R+C,gp=gpar(col=polar_legend$color[i],cex=0.4),default.units="npc")
  #}

  poly<-polygonGrob(x=xy$x, y=xy$y, gp=gpar(fill=NA, col="black",lwd=2),default.units="npc") # outer shell 
  #get aspect ratio of the polygon
  
  gp2<-pts$gp
  gp2$cex<-0.4
  grid.draw(editGrob(pts,gp=gp2))
  seg<-segmentsGrob(x0=seq(-1,1,l=9)*R+C,y0=rep(0,9)*R+C,x1=seq(-1,1,l=9)*R+C,y1=rep(0-0.1,9)*R+C,default.units="npc",gp=gpar(lwd=2))
  label<-seq(1,0,l=5)
  label<-c(label,rev(label)[-1])
  tx1<-textGrob(label=label,x=seq(-1,1,l=9)*R+C,y=rep(0-0.2,9)*R+C,default.units="npc",rot=90,hjust=1,gp=gpar(cex=0.8))
  tx2<-textGrob(label=treatmentLabel,x=c(-1.1,1.1)*R+C,y=rep(0,2)*R+C,rot=90,just=1,gp=gpar(cex=0.8))
  grid.draw(pts)
  grid.draw(poly)
  grid.draw(seg)
  offset<-unit(1,"grobheight",tx1)
  tx2$y<-tx2$y-offset-unit(10,"bigpts")
  grid.draw(tx1)
  grid.draw(tx2)
  upViewport()
}

generate_breaks = function(x, n, center = F){
  if(center){
    m = max(abs(c(min(x, na.rm = T), max(x, na.rm = T))))
    res = seq(-m, m, length.out = n + 1)
  }
  else{
    res = seq(min(x, na.rm = T), max(x, na.rm = T), length.out = n + 1)
  }
  
  return(res)
}

scale_vec_colours = function(x, col = rainbow(10), breaks = NA){
  return(col[as.numeric(cut(x, breaks = breaks, include.lowest = T))])
}

scale_colours = function(mat, col = rainbow(10), breaks = NA){
  mat = as.matrix(mat)
  return(matrix(scale_vec_colours(as.vector(mat), col = col, breaks = breaks), nrow(mat), ncol(mat), dimnames = list(rownames(mat), colnames(mat))))
}

cluster_mat = function(mat, distance, method){
  if(!(method %in% c("ward", "single", "complete", "average", "mcquitty", "median", "centroid"))){
    stop("clustering method has to one form the list: 'ward', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.")
  }
  if(!(distance[1] %in% c("correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) & class(distance) != "dist"){
    print(!(distance[1] %in% c("correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) | class(distance) != "dist")
    stop("distance has to be a dissimilarity structure as produced by dist or one measure  form the list: 'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'")
  }
  if(distance[1] == "correlation"){
    d = as.dist(1 - cor(t(mat)))
  }
  else{
    if(class(distance) == "dist"){
      d = distance
    }
    else{
      d = dist(mat, method = distance)
    }
  }
  
  return(hclust(d, method = method))
}

scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

scale_mat = function(mat, scale){
  if(!(scale %in% c("none", "row", "column"))){
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  }
  mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
  return(mat)
}

generate_annotation_colours = function(annotation, annotation_colors, drop){
  if(is.na(annotation_colors)[[1]][1]){
    annotation_colors = list()
  }
  count = 0
  for(i in 1:ncol(annotation)){
    if(is.character(annotation[, i]) | is.factor(annotation[, i])){
      if (is.factor(annotation[, i]) & !drop){
        count = count + length(levels(annotation[, i]))
      }
      else{
        count = count + length(unique(annotation[, i]))
      }
    }
  }
  
  factor_colors = hsv((seq(0, 1, length.out = count + 1)[-1] + 
      0.2)%%1, 0.7, 0.95)
  
  set.seed(3453)
  
  for(i in 1:ncol(annotation)){
    if(!(colnames(annotation)[i] %in% names(annotation_colors))){
      if(is.character(annotation[, i]) | is.factor(annotation[, i])){
        n = length(unique(annotation[, i]))
        if (is.factor(annotation[, i]) & !drop){
          n = length(levels(annotation[, i]))
        }
        ind = sample(1:length(factor_colors), n)
        annotation_colors[[colnames(annotation)[i]]] = factor_colors[ind]
        l = levels(as.factor(annotation[, i]))
        l = l[l %in% unique(annotation[, i])]
        if (is.factor(annotation[, i]) & !drop){
          l = levels(annotation[, i])
        }
        names(annotation_colors[[colnames(annotation)[i]]]) = l
        factor_colors = factor_colors[-ind]
      }
      else{
        r = runif(1)
        annotation_colors[[colnames(annotation)[i]]] = hsv(r, c(0.1, 1), 1)
      }
    }
  }
  return(annotation_colors)
}



generate_row_annotation_colours = function(annotation, annotation_colors, drop){
  if(is.na(annotation_colors)[[1]][1]){
    annotation_colors = list()
  }
  count = 0
  for(i in 1:ncol(annotation)){
    if(is.character(annotation[, i]) | is.factor(annotation[, i])){
      if (is.factor(annotation[, i]) & !drop){
        count = count + length(levels(annotation[, i]))
      }
      else{
        count = count + length(unique(annotation[, i]))
      }
    }
  }
  
  factor_colors = hsv((seq(0, 1, length.out = count + 1)[-1] + 
      0.2)%%1, 0.7, 0.95)
  #factor_colors = brewer.pal(name="Paired",n=count)
  
  set.seed(3453)
  
  for(i in 1:ncol(annotation)){
    if(!(colnames(annotation)[i] %in% names(annotation_colors))){
      if(is.character(annotation[, i]) | is.factor(annotation[, i])){
        n = length(unique(annotation[, i]))
        if (is.factor(annotation[, i]) & !drop){
          n = length(levels(annotation[, i]))
        }
        ind = sample(1:length(factor_colors), n)
        annotation_colors[[colnames(annotation)[i]]] = factor_colors[ind]
        l = levels(as.factor(annotation[, i]))
        l = l[l %in% unique(annotation[, i])]
        if (is.factor(annotation[, i]) & !drop){
          l = levels(annotation[, i])
        }
        names(annotation_colors[[colnames(annotation)[i]]]) = l
        factor_colors = factor_colors[-ind]
      }
      else{
        r = runif(1)
        annotation_colors[[colnames(annotation)[i]]] = hsv(r, c(0.1, 1), 1)
      }
    }
  }
  return(annotation_colors)
}

kmeans_pheatmap = function(mat, k = min(nrow(mat), 150), sd_limit = NA, ...){
  
  # Filter data
  if(!is.na(sd_limit)){
    s = apply(mat, 1, sd)
    mat = mat[s > sd_limit, ]	
  }
  
  # Cluster data
  set.seed(1245678)
  km = kmeans(mat, k, iter.max = 100)
  mat2 = km$centers
  
  # Compose rownames
  t = table(km$cluster)
  rownames(mat2) = sprintf("cl%s_size_%d", names(t), t)
  
  # Draw heatmap
  pheatmap(mat2, ...)
}

#' A function to draw clustered heatmaps.
#' 
#' A function to draw clustered heatmaps where one has better control over some graphical 
#' parameters such as cell size, etc. 
#' 
#' The function also allows to aggregate the rows using kmeans clustering. This is 
#' advisable if number of rows is so big that R cannot handle their hierarchical 
#' clustering anymore, roughly more than 1000. Instead of showing all the rows 
#' separately one can cluster the rows in advance and show only the cluster centers. 
#' The number of clusters can be tuned with parameter kmeans_k.
#'
#' @param mat numeric matrix of the values to be plotted.
#' @param color vector of colors used in heatmap.
#' @param kmeans_k the number of kmeans clusters to make, if we want to agggregate the 
#' rows before drawing heatmap. If NA then the rows are not aggregated.
#' @param breaks a sequence of numbers that covers the range of values in mat and is one 
#' element longer than color vector. Used for mapping values to colors. Useful, if needed 
#' to map certain values to certain colors, to certain values. If value is NA then the 
#' breaks are calculated automatically.
#' @param border_color color of cell borders on heatmap, use NA if no border should be 
#' drawn.
#' @param cellwidth individual cell width in points. If left as NA, then the values 
#' depend on the size of plotting window.
#' @param cellheight individual cell height in points. If left as NA, 
#' then the values depend on the size of plotting window.
#' @param scale character indicating if the values should be centered and scaled in 
#' either the row direction or the column direction, or none. Corresponding values are 
#' \code{"row"}, \code{"column"} and \code{"none"}
#' @param cluster_rows boolean values determining if rows should be clustered,
#' @param cluster_cols boolean values determining if columns should be clustered.
#' @param clustering_distance_rows distance measure used in clustering rows. Possible 
#' values are \code{"correlation"} for Pearson correlation and all the distances 
#' supported by \code{\link{dist}}, such as \code{"euclidean"}, etc. If the value is none 
#' of the above it is assumed that a distance matrix is provided.
#' @param clustering_distance_cols distance measure used in clustering columns. Possible 
#' values the same as for clustering_distance_rows.
#' @param clustering_method clustering method used. Accepts the same values as 
#' \code{\link{hclust}}.
#' @param treeheight_row the height of a tree for rows, if these are clustered. 
#' Default value 50 points.
#' @param treeheight_col the height of a tree for columns, if these are clustered. 
#' Default value 50 points.
#' @param legend logical to determine if legend should be drawn or not.
#' @param legend_breaks vector of breakpoints for the legend.
#' @param legend_labels vector of labels for the \code{legend_breaks}.
#' @param annotation data frame that specifies the annotations shown on top of the 
#' columns. Each row defines the features for a specific column. The columns in the data 
#' and rows in the annotation are matched using corresponding row and column names. Note 
#' that color schemes takes into account if variable is continuous or discrete.
#' @param annotation_colors list for specifying annotation track colors manually. It is 
#' possible to define the colors for only some of the features. Check examples for 
#' details.
#' @param annotation_legend boolean value showing if the legend for annotation tracks 
#' should be drawn. 
#' @param drop_levels logical to determine if unused levels are also shown in the legend
#' @param show_rownames boolean specifying if column names are be shown.
#' @param show_colnames boolean specifying if column names are be shown.
#' @param main the title of the plot
#' @param fontsize base fontsize for the plot 
#' @param fontsize_row fontsize for rownames (Default: fontsize) 
#' @param fontsize_col fontsize for colnames (Default: fontsize) 
#' @param display_numbers logical determining if the numeric values are also printed to 
#' the cells. 
#' @param number_format format strings (C printf style) of the numbers shown in cells. 
#' For example "\code{\%.2f}" shows 2 decimal places and "\code{\%.1e}" shows exponential 
#' notation (see more in \code{\link{sprintf}}).    
#' @param fontsize_number fontsize of the numbers displayed in cells
#' @param filename file path where to save the picture. Filetype is decided by 
#' the extension in the path. Currently following formats are supported: png, pdf, tiff,
#'  bmp, jpeg. Even if the plot does not fit into the plotting window, the file size is 
#' calculated so that the plot would fit there, unless specified otherwise.
#' @param width manual option for determining the output file width in inches.
#' @param height manual option for determining the output file height in inches.
#' @param row_annotation data frame that specifies the annotations shown on the 
#' rows. Each row defines the features for a specific row. The rows in the data 
#' and rows in the annotation are matched using corresponding row names.The category labels are
#' given by the data frame column names.
#' @param row_annotation_legend same interpretation as the column parameters.
#' @param row_annotation_colors same interpretation as the column parameters
#' @param cytokine_annotation is a data frame of binary factors with levels 0,1, that defines combinations
#' of the categories for each column. They will be colored by their degree of functionality and ordered by degree of functionality
#' and by amount of expression if column clustering is not done.
#' @param headerplot is a list with two components, order and data. Order tells how to reorder the columns of the matrix. 
#' @param polar Boolean; if \code{TRUE} we draw a polar legend. Primarily for
#'  internal use.
#' Data is some summary statistic over the columns which will be plotted in the header where the column cluster tree usually appears.
#' Cytokine ordering is ignored when the headerplot argument is passed.
#' @param \dots graphical parameters for the text used in plot. Parameters passed to 
#' \code{\link{grid.text}}, see \code{\link{gpar}}. 
#' 
#' @return 
#' Invisibly a list of components 
#' \itemize{
#' 	\item \code{tree_row} the clustering of rows as \code{\link{hclust}} object 
#' 	\item \code{tree_col} the clustering of columns as \code{\link{hclust}} object
#' 	\item \code{kmeans} the kmeans clustering of rows if parameter \code{kmeans_k} was 
#' specified 
#' }
#' 
#' @author  Original version by Raivo Kolde <rkolde@@gmail.com>, with modifications
#'  by Greg Finak <gfinak@@fhcrc.org> and Kevin Ushey <kushey@@fhcrc.org>.
#' @examples
#'  # Generate some data
#' test = matrix(rnorm(200), 20, 10)
#' test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
#' test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
#' test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
#' colnames(test) = paste("Test", 1:10, sep = "")
#' rownames(test) = paste("Gene", 1:20, sep = "")
#' 
#' # Draw heatmaps
#' pheatmap(test)
#' pheatmap(test, kmeans_k = 2)
#' pheatmap(test, scale = "row", clustering_distance_rows = "correlation")
#' pheatmap(test, color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
#' pheatmap(test, cluster_row = FALSE)
#' pheatmap(test, legend = FALSE)
#' pheatmap(test, display_numbers = TRUE)
#' pheatmap(test, display_numbers = TRUE, number_format = "%.1e")
#' pheatmap(test, cluster_row = FALSE, legend_breaks = -1:4, legend_labels = c("0", 
#' "1e-4", "1e-3", "1e-2", "1e-1", "1"))
#' pheatmap(test, cellwidth = 15, cellheight = 12, main = "Example heatmap")
#' #pheatmap(test, cellwidth = 15, cellheight = 12, fontsize = 8, filename = "test.pdf")
#' 
#' 
#' # Generate column annotations
#' annotation = data.frame(Var1 = factor(1:10 %% 2 == 0,
#'                              labels = c("Class1", "Class2")), Var2 = 1:10)
#' annotation$Var1 = factor(annotation$Var1, levels = c("Class1", "Class2", "Class3"))
#' rownames(annotation) = paste("Test", 1:10, sep = "")
#'
#' pheatmap(test, annotation = annotation)
#' pheatmap(test, annotation = annotation, annotation_legend = FALSE)
#' pheatmap(test, annotation = annotation, annotation_legend = FALSE, drop_levels = FALSE)
#' 
#' # Specify colors
#' Var1 = c("navy", "darkgreen")
#' names(Var1) = c("Class1", "Class2")
#' Var2 = c("lightgreen", "navy")
#' 
#' ann_colors = list(Var1 = Var1, Var2 = Var2)
#' 
#' #Specify row annotations
#' row_ann <- data.frame(foo=gl(2,nrow(test)/2),`Bar`=relevel(gl(2,nrow(test)/2),"2"))
#' rownames(row_ann)<-rownames(test)
#' pheatmap(test, annotation = annotation, annotation_legend = FALSE, drop_levels = FALSE,row_annotation = row_ann)
#' 
#' #Using cytokine annotations
#' M<-matrix(rnorm(8*20),ncol=8)
#' row_annotation<-data.frame(A=gl(4,nrow(M)/4),B=gl(4,nrow(M)/4))
#' eg<-expand.grid(factor(c(0,1)),factor(c(0,1)),factor(c(0,1)))
#' colnames(eg)<-c("IFNg","TNFa","IL2")
#' rownames(eg)<-apply(eg,1,function(x)paste0(x,collapse=""))
#' rownames(M)<-1:nrow(M)
#' colnames(M)<-rownames(eg)
#' cytokine_annotation=eg
#' pheatmap(M,annotation=annotation,row_annotation=row_annotation,annotation_legend=TRUE,row_annotation_legend=TRUE,cluster_rows=FALSE,cytokine_annotation=cytokine_annotation,cluster_cols=FALSE)
#' 
#' # Specifying clustering from distance matrix
#' drows = dist(test, method = "minkowski")
#' dcols = dist(t(test), method = "minkowski")
#' pheatmap(test, clustering_distance_rows = drows, clustering_distance_cols = dcols)
#' @importFrom RColorBrewer brewer.pal
#' @export
pheatmap = function(mat, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), kmeans_k = NA, breaks = NA, border_color = "grey60", cellwidth = NA, cellheight = NA, scale = "none", cluster_rows = TRUE, cluster_cols = TRUE, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "complete",  treeheight_row = ifelse(cluster_rows, 50, 0), treeheight_col = ifelse(cluster_cols, 50, 0), legend = TRUE, legend_breaks = NA, legend_labels = NA, annotation = NA, annotation_colors = NA, annotation_legend = TRUE, drop_levels = TRUE, show_rownames = T, show_colnames = T, main = NA, fontsize = 10, fontsize_row = fontsize, fontsize_col = fontsize, display_numbers = F, number_format = "%.2f", fontsize_number = 0.8 * fontsize, filename = NA, width = NA, height = NA, row_annotation = NA, row_annotation_legend = TRUE, row_annotation_colors=NA, cytokine_annotation=NA, headerplot=NA, polar=FALSE,...){
  #Do the arguments even make sense?
  oldwarn<-options("warn")
  options("warn"=-1)
  if(!is.na(headerplot)&(cluster_cols)){
    options("warn"=0)
    warning("columns will NOT be clustered when providing headerplot.")
    options("warn"=-1)
  }
  
  # The headerplot will order the columns by headerplot$order and plot the data in headerplot$data above the heatmap.
  # clustering of columns should be turned off, and ordering of columns by cytokine degree of functionality will be overridden.
  if(!is.na(headerplot)&!cluster_cols){
    #reorder the columns of the matrix according to headerplot
    mat <- mat[,headerplot$order]
    cluster_cols <- FALSE; #Reset this so we don't do clustering when providing a headerplot.
    treeheight_col = 100 #set to default value as if cluster_cols was true
    treeheight_row = 50
  }
  
  #check that the cytokine annotation (if it exists), is in the right form
  if(!is.na(cytokine_annotation[[1]][1])){
    for(i in 1:ncol(cytokine_annotation)){
      if(!(class(cytokine_annotation[[i]])=="factor" & all(levels(cytokine_annotation[[i]])%in%c(0,1)))){
        stop("cytokine annotation must be a data frame with categorical factors containing levels 0,1")
      }else
        show_colnames<-FALSE
    }
  }
  # Preprocess matrix
  mat = as.matrix(mat)
  if(scale != "none"){
    mat = scale_mat(mat, scale)
    if(is.na(breaks)){
      breaks = generate_breaks(mat, length(color), center = T)
    }
  }
  
  # Kmeans
  if(!is.na(kmeans_k)){
    # Cluster data
    km = kmeans(mat, kmeans_k, iter.max = 100)
    mat = km$centers
    
    # Compose rownames
    t = table(km$cluster)
    rownames(mat) = sprintf("cl%s_size_%d", names(t), t)
  }
  else{
    km = NA
  }
  
  # Do clustering
  if(cluster_rows){
    tree_row = cluster_mat(mat, distance = clustering_distance_rows, method = clustering_method)
    mat = mat[tree_row$order, , drop = FALSE]
  }
  else{
    tree_row = NA
    if(is.na(headerplot)){
      treeheight_row = 0
    }
  }
  
  if(cluster_cols){
    tree_col = cluster_mat(t(mat), distance = clustering_distance_cols, method = clustering_method)
    mat = mat[, tree_col$order, drop = FALSE]
  }
  else{
    tree_col = NA
    if(is.na(headerplot))
      treeheight_col = 0
  }
  
  # Format numbers to be displayed in cells 
  if(display_numbers){
    fmat = matrix(sprintf(number_format, mat), nrow = nrow(mat), ncol = ncol(mat))
    attr(fmat, "draw") = TRUE
  }
  else{
    fmat = matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
    attr(fmat, "draw") = FALSE
  }
  
  
  # Colors and scales
  if(!is.na(legend_breaks[1]) & !is.na(legend_labels[1])){
    if(length(legend_breaks) != length(legend_labels)){
      stop("Lengths of legend_breaks and legend_labels must be the same")
    }
  }
  
  
  if(is.na(breaks[1])&is.numeric(mat)){
    breaks = generate_breaks(as.vector(mat), length(color))
  }
  if (legend & is.na(legend_breaks[1])) {
    legend = grid.pretty(range(as.vector(breaks)))
    names(legend) = legend
  }
  else if(legend & !is.na(legend_breaks[1])){
    legend = legend_breaks[legend_breaks >= min(breaks) & legend_breaks <= max(breaks)]
    
    if(!is.na(legend_labels[1])){
      legend_labels = legend_labels[legend_breaks >= min(breaks) & legend_breaks <= max(breaks)]
      names(legend) = legend_labels
    }
    else{
      names(legend) = legend
    }
  }
  else {
    legend = NA
  }
  if(is.numeric(mat)){
    mat = scale_colours(mat, col = color, breaks = breaks)
  }
  
  # Preparing annotation colors
  if(!is.na(annotation[[1]][1])){
    annotation = annotation[colnames(mat), , drop = F]
    annotation_colors = generate_annotation_colours(annotation, annotation_colors, drop = drop_levels)
  }
  
  #Prepare row annotation colors
  if(!is.na(row_annotation[[1]][1])){
    row_annotation = row_annotation[rownames(mat), , drop=F]
    row_annotation_colors = generate_row_annotation_colours(row_annotation,row_annotation_colors, drop = drop_levels)
  }
  
  #prepare cytokine annotations
  if(!is.na(cytokine_annotation[[1]][1])){
    cytokine_annotation = cytokine_annotation[colnames(mat),ncol(cytokine_annotation):1,drop=FALSE]
    
    #order the columns by max functionality 
    cka<-apply(cytokine_annotation,2,function(x)sum(as.numeric(as.character(x))))
    cytokine_annotation = cytokine_annotation[,order(cka,decreasing=TRUE),drop=FALSE]
    
    if(!cluster_cols&is.na(headerplot)){
      ckr<-apply(cytokine_annotation,1,function(x)sum(as.numeric(as.character(x))))
      cytokine_annotation = cytokine_annotation[order(ckr),]
      mat = mat[,rownames(cytokine_annotation)] 
    }#else if(!is.na(headerplot)){
    #reorder the cytokine column annotation according to headerplot order
    #cytokine_annotation<-cytokine_annotation[headerplot$order,]
    #}
  }
  
  if(!show_rownames){
    rownames(mat) = NULL
  }
  
  if(!show_colnames){
    colnames(mat) = NULL
  }
  
  # Draw heatmap
  heatmap_motor(mat, border_color = border_color, cellwidth = cellwidth, cellheight = cellheight, treeheight_col = treeheight_col, treeheight_row = treeheight_row, tree_col = tree_col, tree_row = tree_row, filename = filename, width = width, height = height, breaks = breaks, color = color, legend = legend, annotation = annotation, annotation_colors = annotation_colors, annotation_legend = annotation_legend, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, fmat = fmat, fontsize_number = fontsize_number, row_annotation = row_annotation, row_annotation_legend = row_annotation_legend, row_annotation_colors = row_annotation_colors, cytokine_annotation = cytokine_annotation, headerplot=headerplot, polar=polar,...)
  options("warn"=oldwarn$warn)
  invisible(list(tree_row = tree_row, tree_col = tree_col, kmeans = km))
}

