navigation_server <- function(input, output, session) {
  observeEvent(input$tophylogentab, {
    updateTabItems(session, "tabs", selected = "PhylogenAnalysis")
  })

  observeEvent(input$tomsatab, {
    updateTabItems(session, "tabs", selected = "MSAtab")
  })

  observeEvent(input$topcatab, {
    updateTabItems(session, "tabs", selected = "PCAtab")
  })

  observeEvent(input$tostructuretab, {
    updateTabItems(session, "tabs", selected = "PopStructure")
  })
}
