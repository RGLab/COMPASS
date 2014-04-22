library(lattice)
library(shinyGridster)

source("withHelp.R")

zoomButton <- function(inputId) {
  tags$i(class="icon-zoom-in", id=inputId, style="position: absolute; right: 10px; bottom: 10px")
}

splomOutput <- function(outputId) {
  tags$div( id="splom-container",
    tags$div(id=outputId, class="splom")
  )
}

DATA <- readRDS("data/data.rds")
data <- DATA$orig$data
meta <- DATA$data$meta
sid <- DATA$data$sample_id
iid <- DATA$data$individual_id
facet1 <- DATA$facet1
facet2 <- DATA$facet2
facet3 <- DATA$facet3
stimulation <- DATA$stimulation

subsets <- rev( unname(
  transform_subset_label(colnames(DATA$data$n_s)[ -ncol(DATA$data$n_s) ])
) )

if (is.null(stimulation)) {
  stimulations <- DATA$fit$call$treatment[[3]]
} else {
  stimulations <- stimulation
}

facet_vars <- names(meta)
facet_vars <- facet_vars[ !(facet_vars %in% c(sid, iid)) ]

## width, height for gridster + plot elements
width <- 430
height <- 320

## svg output
svgOutput <- function(outputId, width, height) {
  tags$div(
    tag("svg", list(id=outputId, width=width, height=height, class="html-shiny-output"))
  )
}

## Make the HTML representing the facetting variables
make_facets <- function(n) {
  if (n > 3) n <- 3
  width <- floor(100 / n)
  if (n == 1) {
    facets <- selectInput("facet1", label="", choices=c("None", facet_vars), selected=facet1)
  } else {
    facets <- do.call(tagList, lapply(seq_len(n), function(i) {
      tags$div(style=paste0("width: ", width, "%; float: left;"),
        selectInput( paste0("facet", i),
          label=paste("Variable", i),
          choices=c("None", facet_vars),
          selected=get(paste0("facet", i))
        )
      )
    }))
  }
  facets
}

num_meta_vars <- ncol(meta) - 2 ## take away 'sample ID' and 'subject ID'

## ensure that each matrix has the same column names
stopifnot( length( table( table( unlist( lapply( data, names ) ) ) ) ) != 1 )

## markers
markers <- unname( colnames(data[[1]]) )
markers_positive <- paste0( colnames(data[[1]]), "+" )
markers_negative <- paste0( colnames(data[[1]]), "-" )

shinyUI( bootstrapPage(

  includeScript("www/js/d3.js"),
  includeCSS("www/css/styles.css"),
  includeScript("www/js/fancyboxify.js"),

  includeScript("www/jquery-ui/js/jquery-ui-1.10.3.custom.min.js"),
  includeCSS("www/jquery-ui/css/ui-lightness//jquery-ui-1.10.3.custom.min.css"),

  includeScript("www/multiselect/multiselect.js"),
  includeCSS("www/multiselect/multiselect.css"),

  includeScript("www/opentip/opentip.js"),
  includeCSS("www/opentip/opentip.css"),

  includeCSS("www/css/shinySplom.css"),
  includeScript("www/js/shinySplom.js"),

  includeScript("www/js/gridsterExtras.js"),
  includeScript("www/js/tooltips.js"),
  includeScript("scripts.js"),

  ## Allows us to evaluate raw JavaScript code
  tags$head(tags$script(HTML('
      Shiny.addCustomMessageHandler("jsCode",
        function(message) {
          eval(message.value);
        }
      );
    '))),

  singleton( tags$body( style="background-color: #789;" ) ),

  #h1(style="text-align: center; color: white", "Cytokine Visualization"),

  ## Controls exist separate of gridster layout
  tags$div( id="gridster-control-container",

    tags$div(
      id='controls-container',

      HTML("<h3 style='text-align: center;'>ShinyCOMPASS</h3>"),
      HTML("<hr style='margin-top: 0; margin-bottom: 20px;' />"),

      ## multiselect requires the attribute 'multiple' to be set; can't set
      ## this thru regular shiny HTML functions

      h3("Subset Selection"),

      helpIcon("markers-help",
        style="float: right; margin-right: 10px; margin-top: 6px;"
      ),

      h5(style="margin: 6px;", "Marker Selection"),
      HTML("<select id='markers' multiple='multiple'>"),
      HTML(
        paste0("<option value='", markers_positive, "'> ",
        markers, "</option>")
      ),
      HTML("</select>"),

      HTML("<br />"),
      HTML("<br />"),

      ## Visualize subsets...
      helpIcon("subsets-help",
        "float: right; margin-top: 6px; margin-right: 10px;"
      ),

      h5(style="margin: 6px;", "Subsets to Visualize in Histogram"),
      HTML("<select id='subsets' multiple='multiple'>"),
      HTML(
        paste0("<option value='", subsets, "'> ",
          subsets, "</option>")
      ),
      HTML("</select>"),

      HTML("<br />"),
      HTML("<br />"),

      ## overflow: auto keeps div from collapsing to zero height
      ## see: http://stackoverflow.com/questions/218760/how-do-you-keep-parents-of-floated-elements-from-collapsing
      helpIcon("dof-help",
        "float: right; margin-right: 10px; margin-top: 3px;"
      ),
      tags$div(
        h5(style="margin: 6px;", "Range of Degree of Functionality to be Displayed"),
        tags$div( style="width: 50%; float: left;",
          tags$label( `for`="marker_dof_min", "Min"),
          tags$input( id="marker_dof_min", type="number", value="1", min="1", max=ncol( data[[1]] ), step="1" )
        ),
        tags$div( style="width: 50%; float: right;",
          tags$label( `for`="marker_dof_max", "Max"),
          tags$input( id="marker_dof_max", type="number", value="6", min="1", max=ncol( data[[1]] ), step="1" )
        )
      ),

      helpIcon("facets-help",
        style="float: right; margin-right: 10px; margin-top: 12px;"
      ),

      h3("Conditioning Variables"),

      tags$div(
        make_facets(num_meta_vars)
      ),

      h3("Filtering Variables"),

      helpIcon("filter1-help",
        "float: right; margin-right: 10px; margin-top: 3px;"
      ),
      selectizeInput("filter1",
        label="Filter subjects by...",
        choices=c("None", facet_vars)
      ),

      ## this panel will be updated by server.R -- displays available
      ## levels for a factor
      conditionalPanel("input.filter1 != 'None'",
        checkboxGroupInput("filter1_cb", label='', choices='')
      ),
      HTML("<br />"),
      HTML("<hr style='margin-top: 0; margin-bottom: 20px;' />")
    ),
    actionButton("update", "Update View"),
    actionButton("gridster-control-hide", "Show Controls")
  ),

  ## Actual gridster object
  gridster( width=width, height=height,

    gridsterItem(row=1, col=1, sizex=2, sizey=1,
      plotOutput("heatmap", width=width*2, height=height),
      zoomButton("zoom-heatmap")
    ),

    gridsterItem(row=2, col=1, sizex=2, sizey=1,
      plotOutput("polyfunctionality", width=width*2, height=height-20),
      zoomButton("zoom-polyfunctionality")
    ),

    gridsterItem(row=1, col=3, sizex=1, sizey=2,
      h3( style="text-align: center;", "Data Summary"),

      HTML("<hr style='margin-top: 0px;' />"),

      h4("Experiment Description"),
      lapply(DATA$description, function(x) {
        HTML(x, "<br />", "<br />")
      }),

      HTML("<hr style='margin-top: 0px;' />"),

      h4("COMPASS Data Description"),
      p(
        strong("Number of Subjects:"),
        length(unique(DATA$orig$meta[[iid]]))
      ),
#       p(
#         strong("Number of Paired Samples:"),
#         nrow(DATA$data$n_s)
#       ),
      p(
        strong("Number of Markers:"),
        ncol(DATA$orig$data[[1]])
      ),
      p(
        strong("Number of Subsets:"),
        nrow(DATA$data$categories),
        "of",
        2^ncol(DATA$orig$data[[1]]),
        "possible subsets"
      ),
      p(
        strong(
          if (length(stimulations) > 1) "Stimulations" else "Stimulation", " applied:"
        ),
        paste(stimulations, collapse=", ")
      )
    ),

    gridsterItem(row=3, col=1, sizex=2, sizey=1,
      ## custom plot output -- set style manually
      tags$div( style=paste0(
        "width: ", width*2, "px; ",
        "height: ", height, "px; "
      ),
        tags$div( id="posterior_plot", class="shiny-plot-output",
          style=paste0(
            "width: ", width*2, "px; ",
            "height: ", height, "px; "
          )
        ),
        zoomButton("zoom-boxplot")
      )
    )

  )

))
