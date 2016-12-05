library(shiny)

shinyUI(fluidPage(
    titlePanel(title="Le titre : un essai"),
    sidebarLayout(
      sidebarPanel(
        h3("Parameters: INPUT"),
        sliderInput("nb",
                    "Number of elements:",
                    min = 20,
                    max = 500,
                    value = 100,
                    step = 20)
      ),
      mainPanel(
        h3("Results: OUTPUT"), 
        plotOutput("network")
      )
    )
  )
)