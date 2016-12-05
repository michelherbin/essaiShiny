library(shiny)
require(mvtnorm)

shinyServer(
  function(input, output){
    nb <- reactive(input$nb)
    output$network <- renderPlot({
      source <- 215
      n1 <- nb()/2
      n2 <- (nb()/5)*2
      n3 <- nb()/10
      set.seed(source)
      data1 <- rmvnorm(n=n1, 
                       mean=c(10, 15), 
                       sigma=matrix(c(3, 0, 0, 3), nrow=2, ncol=2))
      set.seed(source)
      data2 <- rmvnorm(n=n2, 
                       mean=c(15, 8), 
                       sigma=matrix(c(2.5, 0, 0, 2.5), nrow=2, ncol=2))
      set.seed(source)
      data3 <- rmvnorm(n=n3, 
                       mean=c(5, 7), 
                       sigma=matrix(c(2, 0, 0, 2), nrow=2, ncol=2))
      donnees <- rbind(data1, data2, data3)
      label1 <- rep(1, n1)
      label2 <- rep(2, n2)
      label3 <- rep(3, n3)
      label <- c(label1, label2, label3)
      colnames(donnees)<-c("x", "y")
      plot(donnees, 
           xlim=c(0,20), 
           ylim=c(0, 20),
           pch=c(3, 2, 1)[label], col=c("red","blue", "black")[label],
           cex=0.8)
    })
    
  }
)