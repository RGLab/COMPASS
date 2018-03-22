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

source("common_functions.R")



DATA <- readRDS("data/data.rds")
..sid..  <- DATA$data$sample_id
..iid..  <- DATA$data$individual_id

## global variables for the data state
## NULL variables are updated by other functions
fit <- DATA
data <- DATA$data
meta <- DATA$data$meta
orig_markers <- markers(DATA)

trt <- DATA$fit$call$treatment
if (is.language(trt)) {
  ..stid.. <- as.character(trt[[2]])
} else {
  ..stid.. <- as.character(trt)
}

subsets <- unname(transform_subset_label(colnames(data$n_s)))
dof <- str_count(subsets, "\\+")

evalJS <- function(msg) {
  code <- parse(text="session$sendCustomMessage(type='jsCode', list(value=.__JS_PLACEHOLDER__.))")
  code[[1]][[3]][[2]] <- msg ## don't judge me
  eval(code, envir=parent.frame())
}

## This is JavaScript code for destroying and reconstructing the subset
## widget -- use it when you are updating the subsets available for selection
updateSubsets.js <- paste( readLines("www/js/updateSubsets.js"), collapse="\n" )

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

  getPhenotype <- reactive({
    return( "MeanGamma" )
  })

  getFacet1 <- reactive({
    if (is.null(input$facet1)) return(NULL)
    if (input$facet1 == "None") return(NULL)
    else return( input$facet1 )
  })

  getFacet2 <- reactive({
    if (is.null(input$facet2)) return(NULL)
    if (input$facet2 == "None") return(NULL)
    else return( input$facet2 )
  })

  getFacet3 <- reactive({
    if (is.null(input$facet3)) return(NULL)
    if (input$facet3 == "None") return(NULL)
    else return( input$facet3 )
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

  getSubsets <- reactive({
    return(input$subsets)
  })

#   getCustomFilter <- reactive({
#     return( input$custom_filter )
#   })

  getFilter1 <- function() {
    if (input$filter1 == "None") {
      return(NULL)
    } else {
      return(input$filter1)
    }
  }

  getFilter1Cb <- reactive({
    if (input$filter1 == "None") {
      return(NULL)
    } else {
      return(input$filter1_cb)
    }
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
    return(Inf)
  })

  getMarkersToMarginalizeOver <- reactive({
    return( input$markers_to_marginalize_over )
  })

  ## an observer for conditionally updating the filter widget,
  ## depending on the type of variable being used
  observe({
    x <- input$filter1
    if (!is.null(x)) {
      m <- as.factor(data$meta[[x]])
      updateCheckboxGroupInput(
        session,
        "filter1_cb",
        choices=levels(m),
        selected=levels(m)
      )

      } else {

        ## TODO: numeric filter

      }
  })

  ## an observer for updating the subsets available, based on the markers
  ## wanted + degree of functionality range
  observe({

    marker_dof <- getMarkerOrder()
    markers <- getMarkers()

    ## Contains-selected-marker subsetting
    if (is.null(markers)) {
      keep.ind <- seq_along(subsets)
    } else {
      keep <- Reduce(intersect, lapply(markers, function(x) {
        grep(x, subsets, value=TRUE, fixed=TRUE)
      }))
      keep.ind <- match(keep, subsets)
    }

    ## DOF subsetting
    dof.ind <- which(dof >= marker_dof[1] & dof <= marker_dof[2])

    ## Combine them
    ind <- intersect(keep.ind, dof.ind)

    subsets <- subsets[ind]
    dof <- dof[ind]
    subsets <- subsets[ order(dof, decreasing=TRUE) ]

    subsetsHTML <- tagList(HTML("<select id='subsets' multiple='multiple'>"),
      HTML(
        paste0("<option value='", subsets, "'> ",
          subsets, "</option>")
      ),
      HTML("</select>")
    )
    code <- sprintf("$('#subsets').html(\"%s\")", subsetsHTML )
    code <- gsub("\n", " ", code)
    evalJS( gsub("\n", " ", code) )
    evalJS( updateSubsets.js )

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

      facet1 <- getFacet1()
      facet2 <- getFacet2()
      facet3 <- getFacet3()

      if (length(markers)) {
        must_express <- paste(gsub("+", "", markers, fixed=TRUE), collapse="&")
      } else {
        must_express <- NULL
      }

      if (length(filter1) && length(filter1_cb)) {
        subset_call <- call("%in%", as.symbol(filter1), filter1_cb)
      } else {
        subset_call <- NULL
      }

      plot(fit,
        y=c(facet1, facet2, facet3),
        subset=subset_call,
        minimum_dof=marker_dof[1],
        maximum_dof=marker_dof[2],
        must_express=must_express,
        main=DATA$main
      )

    })

  })

  output$polyfunctionality <- renderPlot({

    renderOnUpdateButtonPress({

      facet1 <- getFacet1()
      facet2 <- getFacet2()
      facet3 <- getFacet3()

      filter1 <- getFilter1()
      filter1_cb <- getFilter1Cb()

      ## Get a melted PolyfunctionalityScore / FunctionalityScore dataset
      df <- data.frame(
        PolyfunctionalityScore=PolyfunctionalityScore(fit),
        FunctionalityScore=FunctionalityScore(fit)
      )
      df$Subject <- rownames(df)

      ## Make sure we merge in the metadata
      meta_dt <- data.table(meta)
      meta_collapsed <- meta_dt[ meta_dt[, .I[1], by=..iid..]$V1 ]

      ## Respect any filtering done
      if (length(filter1) && length(filter1_cb)) {
        subset_call <- call("%in%", as.symbol(filter1), filter1_cb)
        meta_collapsed <- droplevels(meta_collapsed[ eval(subset_call), ])
      }

      df <- merge(df, as.data.frame(meta_collapsed), by.x="Subject", by.y=..iid..)

      pf <- melt(df, value.vars=c("FunctionalityScore", "PolyfunctionaliyScore"),
        variable.name="FunctionalityType",
        value.name="Score"
      )

      pf$FunctionalityType <- factor(pf$FunctionalityType,
        levels=c("FunctionalityScore", "PolyfunctionalityScore"),
        labels=c("Functionality Score", "Polyfunctionality Score")
      )

      if (!is.null(facet3)) {

        p <- ggplot(pf, aes_string(y="Score", x=facet2, fill=facet1)) +
          geom_boxplot(outlier.size = 0) +
          geom_point(position=position_jitterdodge()) +
          facet_grid(paste(facet3, "~", "FunctionalityType"), scales="free_y")

      } else if (!is.null(facet2)) {

        p <- ggplot(pf, aes_string(y="Score", x=facet2, fill=facet1)) +
          geom_boxplot(outlier.size = 0) +
          geom_point(position=position_jitterdodge()) +
          facet_wrap(~ FunctionalityType, scales="free_y")

      } else if (!is.null(facet1)) {

        p <- ggplot(pf, aes_string(y="Score", x=facet1, fill=facet1)) +
          geom_boxplot(outlier.size = 0) +
          geom_point(position=position_jitterdodge()) +
          facet_wrap(~ FunctionalityType, scales="free_y") +
          xlab("")

      } else {

        p <- ggplot(pf, aes_string(x="factor(1)", y="Score")) +
          geom_boxplot(outlier.size = 0) +
          geom_point(position=position_jitter(width=0.1)) +
          facet_wrap(~ FunctionalityType, scales="free_y") +
          xlab("") +
          theme(
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()
          )

      }

      p <- p +
        ggtitle("Functionality and Polyfunctionality Scores") +
        theme(plot.title=element_text(face="bold", size=12, vjust=1))

      print(p)

    })

  })

  output$posterior_plot <- renderPlot({

    renderOnUpdateButtonPress({

      markers <- getMarkers()
      marker_filter <- getMarkerFilter()
      marker_dof <- getMarkerOrder()
      phenotype <- getPhenotype()
      filter1 <- getFilter1()
      filter1_cb <- getFilter1Cb()
      subsets <- getSubsets()

      facet1 <- getFacet1()
      facet2 <- getFacet2()
      facet3 <- getFacet3()

      m <- DATA$fit$mean_gamma
      colnames(m) <- colnames(DATA$data$n_s)

      if (is.null(subsets)) {

        ## Filter by Degree of Functionality
        dof_keep <- intersect(
          which(DATA$fit$categories[, "Counts"] >= marker_dof[1]),
          which(DATA$fit$categories[, "Counts"] <= marker_dof[2])
        )

        ## Filter by markers that must be included
        markers_rex <- gsub("+", "", markers, fixed=TRUE)
        markers_rex <- paste0("(?<!!)", markers_rex)
        filter_keep <- Reduce( intersect, lapply(markers_rex, function(x) {
          grep(x, colnames(DATA$data$n_s), perl=TRUE)
        }))

        ## Only use the top subset
        keep <- max( intersect(dof_keep, filter_keep) )

        m <- m[, keep, drop=FALSE]

      } else {

        colnames(m) <- transform_subset_label(colnames(m))
        m <- m[, subsets, drop=FALSE]

      }

      m <- melt(m)
      names(m) <- c(..iid.., "Subset", "Value")

      ## Respect any filtering done
      if (length(filter1) && length(filter1_cb)) {
        subset_call <- call("%in%", as.symbol(filter1), filter1_cb)
        meta <- as.data.table(meta)
        meta <- droplevels(meta[ eval(subset_call), ])
      }

      ## Make sure we remove the samples from the meta data before merging in
      ## the extra metadata information
      meta <- as.data.table(meta)
      meta_collapsed <- meta[meta[, .I[1], by=..iid..]$V1]
      meta_collapsed <- as.data.frame(meta_collapsed)
      m <- merge(m, meta_collapsed)

      if (is.null(subsets)) {
        m$Subset <- transform_subset_label(m$Subset)
      }

      if (is.null(facet1)) {
        facet1 <- "."
        aes <- aes(x=Value)
      } else {
        aes <- aes_string(x="Value", fill=facet1)
      }

      p <- ggplot(m, aes) +
        geom_histogram() +
        facet_grid( paste(facet1, "~ Subset") ) +
        theme(
          legend.position="none"
        ) +
        xlab("Probability of Response") +
        ylab("Number of Subjects") +
        ggtitle("Histogram of Ag-specific probabilities for selected subsets") +
        theme(plot.title=element_text(face="bold", size=12, vjust=1))

      print(p)

    })

  })

})
