# app.R
options(future.globals.maxSize = 16 * 1024^3)  # Increase to 2 GiB
library(shiny)
library(shinydashboard)
library(shinyjs)

# Load your modules here
source("Scripts/load_data.R")
source("Scripts/quality_control.R")
source("Scripts/doublet_removal.R")  # Add this line to source doublet_removal.R

ui <- dashboardPage(
  dashboardHeader(title = "Heraclomics"),
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",  # Add an ID for the sidebar menu
      menuItem("Home", tabName = "home"),
      menuItem("Load Data", tabName = "load_data"),
      menuItem("Quality Control", tabName = "quality_control"),
      menuItem("Doublet Removal", tabName = "doublet_removal")  # Add menu item for Doublet Removal
    )
  ),
  dashboardBody(
    useShinyjs(),  # Initialize shinyjs
    tabItems(
      tabItem(tabName = "home",
              h2("Welcome to Heraclomics"),
              p("Please select a dataset to start.")
      ),
      tabItem(tabName = "load_data",
              load_data_ui("load_data")  # Add unique 'id' here
      ),
      tabItem(tabName = "quality_control",
              quality_control_ui("quality_control")  # Add unique 'id' here
      ),
      tabItem(tabName = "doublet_removal",
              doublet_removal_ui("doublet_removal")  # Add unique 'id' here
      )
    )
  )
)

server <- function(input, output, session) {
  # Create a reactiveValues object for the application state
  app_state <- reactiveValues(data = NULL)
  
  # Call the server logic for load_data
  callModule(load_data_server, "load_data", app_state)  # 'id' should match the one in UI
  
  # Call the server logic for quality_control
  callModule(quality_control_server, "quality_control", app_state)  # 'id' should match the one in UI
  
  # Call the server logic for doublet_removal
  callModule(doublet_removal_server, "doublet_removal", app_state)  # 'id' should match the one in UI
}

shinyApp(ui, server)
