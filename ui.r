library(shiny)
require(visNetwork)

shinyUI(fluidPage(
    titlePanel(title="Subsampling in 2D dataset"),
    sidebarLayout(
      sidebarPanel(
        h3("Parameters:"),
        sliderInput("nb",
                    "Number of points:",
                    min = 20,
                    max = 500,
                    value = 100,
                    step = 20),
        radioButtons("distrib", 
                     "Kind of overlapped distributions:", 
                     list("Uniform"=1, "Gaussian"=2), 
                     1),
        sliderInput("maxsize",
                    "Maximum of sub-samples:",
                    min = 1,
                    max = 50,
                    value = 5,
                    step = 1)
      ),
      mainPanel(
        h2("Selection of prototypes in a 2D dataset"), 
        tabsetPanel(type="tab",
                    tabPanel("Plot", plotOutput("distribution") ),
                    tabPanel("Skeleton", plotOutput("skeleton") ),
                    tabPanel("Sampling", plotOutput("sampling"), textOutput("number") ),
                    tabPanel("VisNetwork", visNetworkOutput("network")),
                    tabPanel("Read me", htmlOutput("readme"))
        )
        
      )
    )
  )
)