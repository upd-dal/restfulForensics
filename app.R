library(shiny)
library(shinyjs)
library(shinydashboard)
library(shinybusy)
library(plotly)
library(dplyr)

options(shiny.maxRequestSize = 5000 * 1024^2)

source_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
suppressWarnings({invisible(lapply(source_files, source))})

suppressWarnings({
   shinyApp(
      ui = app_ui,
      server = app_server
   )
   })
    