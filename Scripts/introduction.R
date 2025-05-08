# Scripts/introduction.R

library(shiny)
library(shinyjs)
library(shinyalert)
library(shinydashboard)

# Introduction Module UI Function
introductionUI <- function(id) {
  ns <- NS(id)
  tagList(
    useShinyjs(),
    h1("Welcome to Heraclomics"),
    p("Heraclomics is a comprehensive Shiny app for single-cell RNA-seq data analysis."),
    p("Please choose an option to begin:"),
    actionButton(ns("start_new"), "Start New Analysis", icon = icon("play"), class = "btn-success"),
    actionButton(ns("load_previous"), "Load Previous Analysis", icon = icon("upload"), class = "btn-primary")
  )
}

# Introduction Module Server Function
introductionServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive value to indicate starting a new analysis
    start_new_analysis <- reactiveVal(FALSE)
    
    # Observe when the "Start New Analysis" button is clicked
    observeEvent(input$start_new, {
      # Set the reactive value to TRUE
      start_new_analysis(TRUE)
    })
    
    # Observe when the "Load Previous Analysis" button is clicked
    observeEvent(input$load_previous, {
      # Implement loading of previous analysis if desired
      shinyalert("Feature Not Implemented", "Loading previous analyses is not currently implemented.", type = "warning")
    })
    
    # Return the reactive value
    return(list(start_new_analysis = start_new_analysis))
  })
}
