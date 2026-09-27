restful_header <- function() {
  dashboardHeader(
    title = div(
      tags$span("restful forensics",
        style = "font-family: Carme, sans-serif; font-size: 26px; color: #ffffff; vertical-align: middle; padding-left: 0px;"
      )
    ),
    titleWidth = 300,
    header = tags$li(
      class = "dropdown",
      style = "font-family: 'Carmen', sans-serif; padding: 10px;",
      actionButton("refreshApp", "Refresh App", icon = icon("rotate"))
    )
  )
}
