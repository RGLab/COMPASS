library(Kmisc)
library(data.table.extras)
library(hexbin)
library(scales)
library(gridExtra)
library(ggplot2)
library(shiny)
library(reshape2)
library(data.table)
library(gtools)
library(stringr)
library(RColorBrewer)
library(COMPASS)

DATA <- readRDS("data/data.rds")
sid <- DATA$COMPASS$data$sample_id
iid <- DATA$COMPASS$data$individual_id

## read in the data
## FIXME: only allocate one reference for the data
## need to properly handle subsetting and filtering in the plot commands
counts <- DATA$counts
cell_data <- DATA$cell_data
combo_data <- DATA$COMPASS$fit$mean_gamma
d <- dat <- droplevels(DATA$data[ eval(DATA$stimulated) ])
l <- unique( d$Marker )
orig_markers <- DATA$markers
meta <- DATA$meta
counts <- DATA$counts

CellCounts <- data.frame(
  names(counts),
  counts
)
names(CellCounts) <- c(sid, "counts")

## a dataset for the polyfunctionality score
tmp <- d[, .I[1], by=c(iid)]
pf <- d[ tmp$V1 ]
pf <- melt(pf, measure.vars=c("PolyfunctionalityScore", "FunctionalityScore"),
  variable.name="FunctionalityType",
  value.name="Score")

source("common_functions.R")

## Used for the d3 splom
renderSplom <- function(expr, env=parent.frame(), quoted=FALSE) {
  func <- exprToFunction(expr, env, quoted)
  function() {
    val <- func()
    df <- as.data.frame(val[[1]])
    id <- as.character(val[[2]])
    #str(df)
    #str(id)
    id_ind <- match(id, names(df))
    cbind( df[-id_ind], df[id_ind] )
  }
}

## A function used to filter a dataset based on a vector of
## markers
filter_markers <- function(d, markers, 
  markers_to_marginalize_over, 
  p, order, phenotype, max_combos) {
  
  output <- copy(d)
  
  ## first, only keep markers that have been selected
  for (marker in markers) {
    output <- droplevels(output[ fgrep(marker, Marker), ])
  }
  
  ## then, remove cytokine combos that have an overall proportion < p
  props <- sort(decreasing=TRUE,
    with(output, tapply_(output[[phenotype]], Marker, function(x) mean(x, na.rm=TRUE)))
  )
  
  keep <- names(props)[ props > p ]
  
  ## keep only terms of a certain order
  keep_order <- str_count(keep, "\\+")
  keep <- keep[ keep_order >= order[1] & keep_order <= order[2] ]
  
  ## and only keep the top 'n' of these combinations
  if (!is.na(max_combos) && max_combos > 0 && length(keep) > max_combos) {
    output <- output[ Marker %in% keep[1:max_combos], ]
  } else {
    output <- output[ Marker %in% keep, ]
  }
  
  ## if something broke, return the last output
  if (nrow(output) == 0) {
    cat("No data left after filtering; returning full data\n")
    return(d)
  }
  
  ## and finally, return
  return(output)
  
}

## a function to filter the data, as based on some custom input from
## the user
customFilter <- function(dat, expr) {
  if( missing(expr) || expr == '' ) {
    return(dat)
  } else {
    return( dat[ eval( parse( text=expr ), envir=dat ), ] )
  }
}

## filter a function based on levels of a factor
filter1 <- function(dat, var, levels) {
  if (var == "None") {
    return(dat)
  } else {
    if (is.null(levels) || levels == "") {
      return(dat)
    } else {
      return(dat[ dat[[var]] %in% levels, ])
    }
  }
}

shinyServer( function(input, output, session) {
  
  ## global variables for the data state
  ## NULL variables are updated by other functions
  d <- dat
  m <- NULL
  dof <- NULL
  row_order <- NULL
  col_order <- NULL
  
  ## define getters
  
  getPhenotype <- reactive({
    return( input$phenotype )
  })
  
  getFacet1 <- reactive({
    return( input$facet1 )
  })
  
  getFacet2 <- reactive({
    return( input$facet2 )
  })
  
  getFacet3 <- reactive({
    return( input$facet3 )
  })
  
  getSample <- reactive({
    return( input$sample )
  })
  
  #   getMarker <- reactive({
  #     return( input$marker )
  #   })
  
  getMarkerFilter <- reactive({
    return( input$marker_filter )
  })
  
  getProportionFilter <- reactive({
    return( input$proportion_filter )
  })
  
  getSamplesToKeep <- reactive({
    return( names(overall_sample_prop)[ overall_sample_prop > input$sample_proportion_filter ] )
  })
  
  getPlotType <- reactive({
    return( input$plot_type )
  })
  
  getIndividual <- reactive({
    return( input$individual )
  })
  
  getHeatmapLevel <- reactive({
    return( input$heatmap_level )
  })
  
  getMarkers <- reactive({
    return( input$markers )
  })
  
  getMarkerOrderMin <- reactive({
    return( input$marker_dof_min )
  })
  
  getMarkerOrderMax <- reactive({
    return( input$marker_dof_max )
  })
  
  getMarkerOrder <- reactive({
    return( c(input$marker_dof_min, input$marker_dof_max) )
  })
  
#   getCustomFilter <- reactive({
#     return( input$custom_filter )
#   })
  
  getFilter1 <- reactive({
    return( input$filter1 )
  })
  
  getFilter1Cb <- reactive({
    return( input$filter1_cb )
  })
  
  #   isMarginal <- reactive({
  #     return( as.logical(input$marginal) )
  #   })
  
  getBoxplotOrientation <- reactive({
    return( input$boxplot_by_marker_orientation )
  })
  
  getFlipHeatmap <- reactive({
    return( input$flip_heatmap )
  })
  
  getFlipLinechart <- reactive({
    return( input$flip_linechart )
  })
  
  getBoxplotCoordFlip <- reactive({
    return( input$boxplot_coord_flip )
  })
  
  getBoxplotManualLimits <- reactive({
    return( input$boxplot_manual_limits )
  })
  
  getBoxplotLowerLimit <- reactive({
    return( input$boxplot_lower_limit )
  })
  
  getBoxplotUpperLimit <- reactive({
    return( input$boxplot_upper_limit )
  })
  
  getMaxCombosToShow <- reactive({
    return( as.integer(input$max_combos_to_show) )
  })
  
  getMarkersToMarginalizeOver <- reactive({
    return( input$markers_to_marginalize_over )
  })
  
  ## observers
  observe({
    x <- input$individual
    samples <- meta[[sid]][ meta[[iid]] == x ]
    updateSelectInput(session, "sample",
      label=sid,
      choices=samples
    )
  })
  
  observe({
    x <- input$boxplot_manual_limits
    updateNumericInput(session, "boxplot_lower_limit", value=0)
    updateNumericInput(session, "boxplot_upper_limit", value=0.005)
  })
  
  ## if the phenotype is changed, make sure the 'set limits manually'
  ## thing gets unchecked
  observe({
    x <- input$phenotype
    updateCheckboxInput(session, "boxplot_manual_limits", value=FALSE)
  })
  
  #   observe({
  #     x <- input$heatmap_level
  #     updateSliderInput(session, 
  #       "proportion_filter", 
  #       label=paste("Remove", x, "with proportion < x")
  #     )
  #   })
  
  ## an observer for conditionally updating the filter widget,
  ## depending on the type of variable being used
  observe({
    x <- input$filter1
    if (x != "None") {
      m <- d[[x]]
      if( is.factor(m) || is.character(m) ) {
        
        updateCheckboxGroupInput(
          session, 
          "filter1_cb", 
          choices=levels(d[[x]]),
          selected=levels(d[[x]])
        )
        
      } else {
        
        ## TODO: numeric filter
        
      }
    }
  })
  
  ## wrap data updates in an observer
  observe({
    
    ## The variables we need to observe for updating the data
    markers <- getMarkers()
    marker_filter <- getMarkerFilter()
    marker_dof <- getMarkerOrder()
    phenotype <- getPhenotype()
    filter1 <- getFilter1()
    filter1_cb <- getFilter1Cb()
#     custom_filter <- getCustomFilter()
    max_combos_to_show <- getMaxCombosToShow()
    markers_to_marginalize_over <- getMarkersToMarginalizeOver()
    
    ## filter the markers
    tmp <- filter_markers(
      dat, 
      markers, 
      markers_to_marginalize_over,
      marker_filter, 
      marker_dof, 
      phenotype, 
      max_combos_to_show
    )
    
#     tmp <- customFilter(tmp, custom_filter)
    tmp <- filter1(tmp, filter1, filter1_cb)
    
    ## Marginalize tmp
    if (length(markers_to_marginalize_over)) {
      other_markers <- orig_markers[ !(orig_markers %in% markers_to_marginalize_over) ]
      tmp[, (phenotype) := mean( get(phenotype), na.rm=TRUE), by=c(sid, iid,other_markers)]
      tmp <- tmp[ tmp[, .I[1], by=c(sid, iid, other_markers)]$V1 ]
      for (marker in markers_to_marginalize_over) {
        tmp[, Marker := gsub( paste0(marker, "+"), "", Marker, fixed=TRUE)]
        tmp[, Marker := gsub( paste0(marker, "-"), "", Marker, fixed=TRUE)]
      }
    }
    
    ## Cast the marker data down into a matrix
    form <- as.formula( paste(sid, "~ Marker") )
    m <<- acast(tmp, form, value.var=phenotype)
    
    ## Rows are samples, columns are marker combinations
    
    ## Generate the marker / column annotations
    annot_tmp <- data.frame(do.call(cbind, lapply(orig_markers, function(x) {
      
      ## the positives
      pos <- as.integer(grepl( paste0(x, "\\+"), colnames(m)))
      
      ## the negatives
      neg <- -as.integer(grepl( paste0(x, "-"), colnames(m)))
      
      ## merge them
      output <- rep(0, length(pos))
      output[pos == 1] <- 1
      output[neg == -1] <- 0
      return(output)
      
    })))
    
    rownames(annot_tmp) <- colnames(m)
    colnames(annot_tmp) <- orig_markers
    
    annot_tmp[] <- lapply(annot_tmp, as.factor)
    
    annot <<- annot_tmp
    
    ## Order the columns by degree of functionality, number of
    ## marginals, and then actual mean phenotype.
    dof <<- apply(annot, 1, function(x) {
      sum(x == 1)
    })
    
    nom <<- apply(annot, 1, function(x) {
      sum(x == 0)
    })
    
    mean_pheno <<- apply(m, 2, function(x) {
      mean(x, na.rm=TRUE)
    })
    
    col_order <<- colnames(m)[do.call(order, list(dof, nom, mean_pheno))]
    
    d <<- tmp
    
  })
  
  ## A helper function for rendering views only when the button is clicked
  renderOnUpdateButtonPress <- function(x, env = parent.frame(), quoted = FALSE) {
    force(input$update)
    fun <- exprToFunction(x, env, quoted)
    isolate( fun() )
  }
  
  ## heatmap plot
  output$heatmap <- renderPlot({
    
    renderOnUpdateButtonPress({
      
      markers <- getMarkers()
      marker_filter <- getMarkerFilter()
      marker_dof <- getMarkerOrder()
      phenotype <- getPhenotype()
      filter1 <- getFilter1()
      filter1_cb <- getFilter1Cb()
      #     custom_filter <- getCustomFilter()
      max_combos_to_show <- getMaxCombosToShow()
      markers_to_marginalize_over <- getMarkersToMarginalizeOver()
      
      facet1 <- getFacet1()
      facet2 <- getFacet2()
      facet3 <- getFacet3()
      heatmap_level <- getHeatmapLevel()
      phenotype <- getPhenotype()
      
      ## A quick function used for getting a nice color palette
      heatmap_color <- function(phenotype) {
        switch(phenotype,
          `MeanGamma`=colorRampPalette(brewer.pal(n=9, name="Purples")[2:7])(100),
          #`ModelLogDiff`=colorRampPalette(brewer.pal(n=9, name="Purples")[2:7])(100),
          #`ModelDiff`=colorRampPalette(brewer.pal(n=9, name="Purples")[2:7])(100),
          colorRampPalette(brewer.pal(n=7, name="PuOr"))(100)
        )
      }
      
      heatmap_breaks <- function(m, color, phenotype) {
        if (phenotype %in% c("LogFoldChange", "ModelLogDiff", "ModelDiff")) {
          top <- max( max(m), abs(min(m)) )
          return( seq(-top, top, length.out=length(color)+1) )
        } else if (phenotype %in% c("MeanGamma")) {
          return( seq(0, 1, length.out=length(color)+1) )
        } else {
          stop("Unrecognized phenotype!")
        }
      }
      
      if (facet1 == "Original Ordering") facet1 <- NULL
      if (facet2 == "None") facet2 <- NULL
      if (facet3 == "None") facet3 <- NULL
      
      color <- heatmap_color(phenotype)
      breaks <- heatmap_breaks(m, color, phenotype)
      
      if (any( !is.null( c(facet1, facet2, facet3) ) )) {
        
        row_annot <- d[, c(sid, facet1, facet2, facet3), with=FALSE]
        row_annot <- as.data.frame(row_annot[ row_annot[, .I[1], by=c(sid)]$V1 ])
        row_annot[[sid]] <- as.character(row_annot[[sid]])
        row_annot <- row_annot[ match( rownames(m), row_annot[[sid]] ), ]
        rownames(row_annot) <- row_annot[[sid]]
        row_annot[[sid]] <- NULL
        
        row_order <<- do.call(order, row_annot[ c(facet1, facet2, facet3) ])
        
        pheatmap(m[row_order, col_order, drop=FALSE], 
          color=color,
          breaks=breaks,
          cluster_rows=FALSE,
          cluster_cols=FALSE,
          show_rownames=FALSE,
          show_colnames=FALSE,
          cytokine_annotation=annot[col_order, , drop=FALSE],
          row_annotation=row_annot[row_order, , drop=FALSE],
          main=paste(phenoToLabel(phenotype), "by Sample"),
          order_by_max_functionality=FALSE
        )
        
      } else {
        
        pheatmap(m[, col_order, drop=FALSE], 
          color=color,
          breaks=breaks,
          cluster_rows=FALSE,
          cluster_cols=FALSE,
          show_rownames=FALSE,
          show_colnames=FALSE,
          cytokine_annotation=annot[col_order, , drop=FALSE],
          main=paste(phenoToLabel(phenotype), "by Sample"),
          order_by_max_functionality=FALSE
        )
        
      }
      
    })
    
  })
  
  ## linechart
  output$linechart <- renderPlot({
    
    renderOnUpdateButtonPress({
    
    markers <- getMarkers()
    if (is.null(markers)) {
      markers <- orig_markers
    }
    marker_filter <- getMarkerFilter()
    marker_dof <- getMarkerOrder()
    phenotype <- getPhenotype()
    filter1 <- getFilter1()
    filter1_cb <- getFilter1Cb()
    #     custom_filter <- getCustomFilter()
    max_combos_to_show <- getMaxCombosToShow()
    markers_to_marginalize_over <- getMarkersToMarginalizeOver()
    
    individual <- getIndividual()
    facet1 <- getFacet1()
    facet2 <- getFacet2()
    heatmap_level <- getHeatmapLevel()
    phenotype <- getPhenotype()
    flip_linechart <- getFlipLinechart()
    
    d_sub <- droplevels( d[ d[[iid]] == individual, ] )
    
    ## make Marker a factor now and ensure it has a correct ordering
    ## Ie, it has to match the column order global
    d_sub$Marker <- factor(d_sub$Marker, levels=col_order)
    
    ## We re-order it as well so the colors are properly plotted
    d_sub <- d_sub[ order(d_sub$Marker), ]
    
    p <- ggplot( d_sub, aes_string(x="Marker", y=phenotype, group=sid)) +
      geom_point() +
      geom_line( aes(x=as.integer( factor(Marker) ))) +
      xlab(paste("Individual:", individual)) +
      ggtitle(paste(phenoToLabel(phenotype), "by Marker")) +
      #guides(color=guide_legend(nrow=10)) +
      theme( 
        axis.text.x = element_text(angle=90, hjust=1),
        title = element_text(family="Gill Sans")
      )
    
    grob <- ggplotGrob(p)
    
    ## extract the linechart panel, and the yaxis
    panel <- grob[[1]][[4]]
    yaxis <- grob[[1]][[2]]
    
    ## generate the y-axis label now, so we know how big the cell
    ## should be
    label <- wrap( phenoToLabel(phenotype), width=30 )
    
    layout <- grid.layout(nrow=6, ncol=3,
      heights=unit(
        c(2, 0, 6, 0.5, 2, 1),
        c("lines", "lines", "null", "lines", "null", "lines")
      ),
      widths=unit(
        c(2 * str_count(label, "\n") + 4, 1, 1),
        c("lines", "null", "lines")
      )
    )
    
    ## prepare the plot
    grid.newpage()
    pushViewport( viewport(layout=layout) )
    
    ## push the main chart panel
    pushViewport( viewport(layout.pos.row=3, layout.pos.col=2) )
    grid.draw(panel)
    popViewport()
    
    ## push the cytokine annotations
    border_color <- "grey"
    draw_annotations <- function(d_sub, 
      palette=brewer.pal( length(orig_markers), "Paired" ),
      border_color=NA) {
      palette <- c("white", palette)
      n <- length(orig_markers)
      m <- nrow(d_sub)
      x <- (1:m)/m - 1/2/m
      y <- ((1:n)-0.5)/n
      for (i in 1:m) {
        tmp <- unlist( d_sub[i, orig_markers, with=FALSE] )
        tmp <- tmp * sum(tmp)
        tmp <- tmp + 1
        fill <- palette[tmp]
        if (!length(fill)) fill <- rep("black", n)
        grid.rect( x=x[i], y=y, width=1/m, height=1/(n+1), gp=gpar(fill=fill, col=border_color))
      }
    }
    
    pushViewport( viewport(layout.pos.row=5, layout.pos.col=2) )
    draw_annotations( d_sub )
    popViewport()
    
    ## push the cytokine annotation labels
    pushViewport( viewport(layout.pos.row=5, layout.pos.col=1) )
    grid.draw(grid.text(
      x=0.2,
      y=(1:length(orig_markers)-0.5)/length(orig_markers),
      just="left",
      label=orig_markers,
      gp=gpar(cex=0.7)
    ))
    popViewport()
    
    ## push the yaxis ticks
    pushViewport( viewport(layout.pos.row=3, layout.pos.col=1) )
    grid.draw(yaxis)
    grid.text(label, rot=90)
    popViewport()
    
    ## the title
    pushViewport( viewport(layout.pos.row=1, layout.pos.col=2) )
    grid.draw(grid.text(label=paste(phenoToLabel(phenotype), "by Marker"),
      gp=gpar(fontfamily="Helvetica")))
    popViewport()
    
    ## pop back to home
    popViewport()
    
    ## fill in the annotations
    
    #     if( facet1 != "Original Ordering" ) {
    #       p <- p + aes_string(x="Marker", y=phenotype, group=sid, color=facet1)
    #       if( facet2 != "None" ) {
    #         p <- p + facet_wrap(facet2)
    #       }
    #     }
    #     
    #     ## y label -- for some reason, it has to occur after the facetting code
    #     p <- p +
    #       ylab( phenoToLabel(phenotype) )
    #     
    #     if (flip_linechart) {
    #       p <- p +
    #         coord_flip() +
    #         theme(
    #           title=element_text(family="Gill Sans")
    #         )
    #     }
    
  })
    
  })
  
  output$dofplot <- renderPlot({
    
    renderOnUpdateButtonPress({
    
    markers <- getMarkers()
    marker_filter <- getMarkerFilter()
    marker_dof <- getMarkerOrder()
    phenotype <- getPhenotype()
    filter1 <- getFilter1()
    filter1_cb <- getFilter1Cb()
    #     custom_filter <- getCustomFilter()
    max_combos_to_show <- getMaxCombosToShow()
    markers_to_marginalize_over <- getMarkersToMarginalizeOver()
    
    individual <- getIndividual()
    facet1 <- getFacet1()
    facet2 <- getFacet2()
    heatmap_level <- getHeatmapLevel()
    phenotype <- getPhenotype()
    flip_dofplot <- input$flip_dofplot
    
    ## get the samples belonging to the currently selected individual
    samples <- cell_data[ names(cell_data) %in% d[[sid]] ]
    
    ## computing the degree of functionality
    ## for each cell, we compute the number of markers that were
    ## expressed (expression value > 0). then, we generate 
    ## counts of these values. for example, given 3 cells in a sample like
    
    ## IL1 TNFa INFg
    ## 0 1 0
    ## 0 0 1
    ## 1 1 0
    
    ## we should get output like
    
    ## sample num_markers_expressed count
    ## 1 0 0
    ## 1 1 2
    ## 1 2 1
    ## 1 3 0
    
    ## the zero counts can be left out
    
    ## compute the degree of functionality:
    ## for each cell, the DOF == the number of markers expressed
    samples_dof <- stack( lapply( samples, function(x) {
      rowSums(x > 0)
    }) )
    
    samples_dof <- data.table(samples_dof)
    samples_dof <- samples_dof[, {
      counts <- Kmisc::counts(values)
      list(counts=counts, order=names(counts))
    }, by=ind]
    ## merge in cell counts
    dof <- merge.data.frame( samples_dof, CellCounts, by.x="ind", by.y=sid, all.x=TRUE )
    dof$prop <- dof$counts.x / dof$counts.y
    
    ## merge in meta-info
    dof <- unique( merge( 
      dof, 
      d[ names(d) %in% c(sid, names(meta)) ], 
      by.x="ind", 
      by.y=sid,
      all.x=TRUE
    ) )
    
    p <- ggplot(dof, aes(x=factor(order), y=prop)) +
      theme(
        title=element_text(family="Gill Sans")
      ) +
      xlab("Degree of Functionality") +
      ylab("Proportion of Cells") +
      ggtitle(paste("Degree of Functionality\nIndividual:", individual))
    
    if( facet1 != "Original Ordering" ) {
      if( facet2 == "None" ) {
        p <- p + 
          geom_bar( aes_string(fill=facet1), stat="identity", position="dodge" )
        #facet_grid( paste( facet1, "~ ." ) )
      } else {
        p <- p + 
          geom_bar( aes_string(fill=facet1), stat="identity", position="dodge") +
          facet_grid( paste( facet2, "~ ." ) )
      }
    } else {
      p <- p +
        geom_bar(stat="identity", position="dodge")
      
    }
    
    if (flip_dofplot) {
      p <- p + coord_flip()
    }
    
    plot(p)
    
  })
    
  })
  
  ## boxplot by marker
  output$boxplot_by_marker <- renderPlot({
    
    renderOnUpdateButtonPress({
    
    markers <- getMarkers()
    marker_filter <- getMarkerFilter()
    marker_dof <- getMarkerOrder()
    phenotype <- getPhenotype()
    filter1 <- getFilter1()
    filter1_cb <- getFilter1Cb()
    #     custom_filter <- getCustomFilter()
    max_combos_to_show <- getMaxCombosToShow()
    markers_to_marginalize_over <- getMarkersToMarginalizeOver()
    
    facet1 <- getFacet1()
    facet2 <- getFacet2()
    facet3 <- getFacet3()
    heatmap_level <- getHeatmapLevel()
    boxplot_orientation <- getBoxplotOrientation()
    boxplot_coord_flip <- getBoxplotCoordFlip()
    boxplot_manual_limits <- getBoxplotManualLimits()
    boxplot_lower_limit <- getBoxplotLowerLimit()
    boxplot_upper_limit <- getBoxplotUpperLimit()
    plot_type <- getPlotType()
    
    if( facet1 == "Original Ordering" ) {
      plot( 1 ~ 1, type='n', axes=FALSE, frame.plot=FALSE, ann=FALSE)
      text( x=1, y=1, labels="Select a Facet to view Boxplots of Proportions by Marker")
      return(NULL)
    }
    
    p <- ggplot(d) +
      theme( 
        #axis.text.x = element_text(angle=45, hjust=1),
        title = element_text(family="Gill Sans")
      )
    
    ## what kind of plot to generate?
    if( facet2 == "None" ) {
      p <- p + switch(plot_type,
        boxplot=geom_boxplot( aes_string(x=facet1, y=phenotype), fill="#ABCDEF"),
        histogram=geom_histogram( aes_string(x=phenotype, fill=facet1), alpha=0.5, position="identity" ),
        density=geom_density( aes_string(x=phenotype, fill=facet1), alpha=0.5 )
      )
    } else {
      p <- p + switch(plot_type,
        boxplot=geom_boxplot( aes_string(x=facet2, y=phenotype, fill=facet1)),
        histogram=geom_histogram( aes_string(x=phenotype, fill=facet1), alpha=0.5, position="identity" ),
        density=geom_density( aes_string(x=phenotype, fill=facet1), alpha=0.5 )
      )
    }
    
    if (boxplot_coord_flip) {
      if (facet3 != "None") {
        p <- p +
          ylab( phenoToLabel(phenotype) ) +
          facet_grid(paste("Marker ~", facet3))
      } else {
        p <- p +
          ylab( phenoToLabel(phenotype) ) +
          facet_grid("Marker ~ .")
      }
    } else {
      if (facet3 != "None") {
        p <- p +
          ylab( phenoToLabel(phenotype) ) +
          facet_grid(paste(facet3, "~ Marker"))
      } else {
        p <- p +
          ylab( phenoToLabel(phenotype) ) +
          facet_grid(". ~ Marker")
      }
    }
    
    if (boxplot_orientation == "Horizontal") {
      if (boxplot_manual_limits) {
        p <- p +
          coord_flip(ylim=c(boxplot_lower_limit, boxplot_upper_limit))
      } else {
        p <- p +
          coord_flip()
      }
    } else {
      if (boxplot_manual_limits) {
        p <- p +
          coord_cartesian(ylim=c(boxplot_lower_limit, boxplot_upper_limit))
      } 
    }
    
    plot(p)
    
  })
    
  })
  
  output$stats <- renderDataTable(
    
    options=list(
      iDisplayLength=3,
      bSortClasses=TRUE,
      bAutoWidth=FALSE,
      aoColumnDefs=list(
        list(
          sWidth=c("500px"),
          aTargets=c( list(0) )
        )
      )
    ), {
      
      renderOnUpdateButtonPress({
      
      markers <- getMarkers()
      marker_filter <- getMarkerFilter()
      marker_dof <- getMarkerOrder()
      phenotype <- getPhenotype()
      filter1 <- getFilter1()
      filter1_cb <- getFilter1Cb()
      #     custom_filter <- getCustomFilter()
      max_combos_to_show <- getMaxCombosToShow()
      markers_to_marginalize_over <- getMarkersToMarginalizeOver()
      
      individual <- getIndividual()
      facet1 <- getFacet1()
      facet2 <- getFacet2()
      facet3 <- getFacet3()
      heatmap_level <- getHeatmapLevel()
      
      if (nrow(d) == 0) {
        return("Empty data")
      }
      
      ## construct 'by' from the included facets
      by <- "Marker"
      
      if (facet1 != "Original Ordering") {
        by <- paste(by, facet1, sep=",")
      }
      
      if (facet2 != "None") {
        by <- paste(by, facet2, sep=",")
      }
      
      if (facet3 != "None") {
        by <- paste(by, facet3, sep=",")
      }
      
      m <- d[, list(
        Mean=mean( get(phenotype), na.rm=TRUE ),
        Median=median( get(phenotype), na.rm=TRUE ),
        SD=sd( get(phenotype), na.rm=TRUE ),
        N=.N
      ), keyby=by]
      
      ## Make the values nicer to present
      m$Mean <- sprintf("%.3G", m$Mean)
      m$Median <- sprintf("%.3G", m$Median)
      m$SD <- sprintf("%.3G", m$SD)
      
      return(m)
      
    })
      
    })
  
  #   output$splom <- renderSplom({
  #     
  #     individual <- getIndividual()
  #     facet1 <- getFacet1()
  #     facet2 <- getFacet2()
  #     heatmap_level <- getHeatmapLevel()
  #     phenotype <- getPhenotype()
  #     filter1 <- getFilter1()
  #     filter1_cb <- getFilter1Cb()
  #     custom_filter <- getCustomFilter()
  #     markers_to_marginalize_over <- getMarkersToMarginalizeOver()
  #     max_combos_to_show <- getMaxCombosToShow()
  #     
  #     ## take the selected individual
  #     
  #     samples <- as.character(unique(d[[sid]][ d[[iid]] == individual ]))
  #     dat <- cell_data[ names(cell_data) %in% samples ]
  #     dat <- lapply(dat, as.data.frame)
  #     m <- rbindlistn(dat, TRUE)
  #     setnames(m, ".Names", sid)
  #     list(m, sid)
  #     
  #   })
  
  #   output$polyfunctionality <- renderPlot({
  #     
  #     individual <- getIndividual()
  #     facet1 <- getFacet1()
  #     facet2 <- getFacet2()
  #     
  #     if (facet1 != "Original Ordering") {
  #       
  #       if (facet2 != "None") {
  #         
  #         p <- ggplot(pf, aes_string(y="Score", x=facet2, fill=facet1)) +
  #           geom_boxplot() +
  #           facet_wrap(~ FunctionalityType)
  #         
  #       } else {
  #         
  #         p <- ggplot(pf, aes_string(y="Score", x=factor(1), fill=facet1)) +
  #           geom_boxplot() +
  #           facet_wrap(~ FunctionalityType) +
  #           xlab("")
  #         
  #       }
  #       
  #     } else {
  #       
  #       p <- ggplot(pf, aes_string(x=factor(1), y="Score")) +
  #         geom_boxplot() +
  #         facet_wrap(~ FunctionalityType) +
  #         xlab("")
  #     }
  #     
  #     print(p)
  #     
  #   })
  
})
