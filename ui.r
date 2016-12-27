library(shiny)
require(visNetwork)

shinyUI(fluidPage(
    titlePanel(title="Selection of prototypes in 2D dataset with sub-sampling"),
    sidebarLayout(
      sidebarPanel(
        h3("Parameters:"),
        sliderInput("seed",
                    "Random generator seed:",
                    min = 1,
                    max = 100,
                    value = 50,
                    step = 1),
        sliderInput("nb",
                    "Number of sensors:",
                    min = 10,
                    max = 150,
                    value = 30,
                    step = 10),
        radioButtons("distrib", 
                     "Select spatial distributions:", 
                     list("Uniform"=1, "Gaussian"=2), 
                     1),
        sliderInput("maxsize",
                    "Maximum of base stations",
                    min = 1,
                    max = 20,
                    value = 5,
                    step = 1)
      ),
      mainPanel(
        h2("Base Stations of a Wireless Sensor Network"), 
        tabsetPanel(type="tab",
                    tabPanel("Spatial Distribution of Sensors", plotOutput("distribution") ),
                    tabPanel("Skeleton", plotOutput("skeleton") ),
                    tabPanel("Base Stations", plotOutput("sampling"), textOutput("number") ),
                    tabPanel("Network", visNetworkOutput("network") ),
                    tabPanel("Read me", htmlOutput("readme") )
        )
        
      )
    )
  )
)