references_tab <- function() {
  # REFERENCES ======================================================================================
  tabItem(
    tabName = "AppRef",
    tabBox(
      width = 12,
      tabPanel(
        title = "References",
        div(
          style = "height: 80vh; overflow-y:scroll;",
          uiOutput("referenceTexts")
        )
      )
    )
  )
}
