library(shiny)
require(mvtnorm)
require(igraph)
require(visNetwork)

shinyServer(
  function(input, output){
    nb <- reactive(input$nb)
    maxsample <- reactive(input$maxsize)
    seed <- reactive(input$seed)

    donnees <- reactive({
      n1 <- nb()/2
      n2 <- (nb()/5)*2
      n3 <- nb()-n1-n2
      source <- seed()
      set.seed(source)
      data1 <- NULL
      if(as.numeric(input$distrib) == 1){
        set.seed(source)
        data1 <- cbind(runif(n = n1, min=5, max = 20), 
                       runif(n = n1, min=5, max = 20))
        set.seed(source)
        data2 <- cbind(runif(n = n2, min=0, max = 10), 
                       runif(n = n2, min=0, max = 10))
        set.seed(source)
        data3 <- cbind(runif(n = n3, min=8, max = 15), 
                       runif(n = n3, min=3, max = 12))
      } else {
        set.seed(source)
        data1 <- rmvnorm(n=n1, 
                       mean=c(10, 12), 
                       sigma=matrix(c(4, 0, 0, 4), nrow=2, ncol=2))
        set.seed(source)
        data2 <- rmvnorm(n=n2, 
                       mean=c(13, 8), 
                       sigma=matrix(c(3, 0, 0, 3), nrow=2, ncol=2))
        set.seed(source)
        data3 <- rmvnorm(n=n3, 
                       mean=c(7, 7), 
                       sigma=matrix(c(2.5, 0, 0, 2.5), nrow=2, ncol=2))
      }
      rbind(data1, data2, data3)
    })
    
    
    lien=function(data, n){
      distances <- as.matrix(dist(data, 
                                  diag = TRUE, upper = TRUE, method= "euclidean"))
      pointdevues <- array(rep(0, n*n), c(n, n))
      for(individu in 1:n){
        pointdevues[individu, ] <- rank(distances[individu, ], ties.method= "random")
      }
      entre=function(x, i, j){
        ligne_i <- pointdevues[i,]
        ligne_j <- pointdevues[j,]
        avant_j <-  (1<ligne_i) & (ligne_i<pointdevues[i,j])
        avant_i <-  (1<ligne_j) & (ligne_j<pointdevues[j,i])
        return( !(sum(avant_i & avant_j)) ) 
      }
      lesliens <- matrix(mapply(entre, pointdevues, row(pointdevues), col(pointdevues)), n, n)
      return(matrix(as.numeric(lesliens)-diag(1,n), n, n))
    }
    
    adjacent <- reactive({
      n <- nb()
      data <- donnees()
      lien(data, n)
    })
    
    
    
    orient=function(adj, n){
      gain<-array(rep(0, n*n), c(n, n))
      for(i in 1:n){
        for(j in 1:n){
          if(adj[i, j] != 0){
            gain[i, j]<-sum(adj[,i])
          }
        }
      }
      perte<-array(rep(0, n*n), c(n, n))
      for(i in 1:n){
        for(j in 1:n){
          if(adj[i, j] != 0){
            perte[i, j]<-sum(adj[j,])
          }
        }
      }
      direct<-adj
      for(i in 1:n){
        for(j in 1:n){
          if(gain[i,j] < perte[i,j]){
            direct[j, i]<-0
          }
        }
      }
      # second step
      for(i in 1:n){
        for(j in 1:n){
          if((direct[i, j] != 0) & (gain[i,j] == perte[i,j])){
            coord <- which(direct[,i] != 0, arr.ind=TRUE)
            gainbis <- max(gain[coord,i])
            coord <- which(direct[j,] != 0, arr.ind=TRUE)
            pertebis<- max(perte[j,coord])
            if(gainbis > pertebis){
              direct[i, j]<-0
            } else {
              direct[j, i] <- 0
            }
          } 
        }
      }
      return(matrix(direct, n, n))
    }
    
    gains <- reactive({
      n <- nb()
      data <- donnees()
      maxsize <- as.numeric(input$maxsize)
      sizesample <- n
      selectsample <- c(1:n)
      echantillon <- data[selectsample,]
      
      while(sizesample > maxsize){
        mat.adjacence <- lien(echantillon, sizesample)
        mat.direction <- orient(mat.adjacence, sizesample)
        op<-NULL
        for(i in 1:n){ 
           if( (sum(mat.direction[i,]) == 0) & (sum(mat.direction[,i]) != 0) ){
            op<-c(op, i)
          }
        }
        sizesample <- length(op)
        if(sizesample > 1){
          n <- length(op)
          data <- echantillon
          echantillon <- data[op,]
          choix <- selectsample
          selectsample <- choix[op]
        }
      }
      sample <- NULL
      if(sizesample == 1) {
        sample <- matrix(echantillon[op,], 1, 2)
        choix <- as.vector(selectsample[op])
      } else {
        if(sizesample == 0) {
          sample <- donnees()
          choix <- c(1:n)
        } else {
          sample <- echantillon
          choix <- selectsample
        }
      }
      choix
    })
    
    
    
    output$readme <- renderUI({
      str1 <- h4("Determination of a subsample of a 2D dataset")
      str2 <- paste( "The ", as.character(nb()), " points of a 2D dataset 
              simulate sensors of a Wireless Sensors Network. 
              The sub-samples could determine the locations of the base stations in WSN.")
      str3 <- " "
      str4 <- "Ref.:"
      str4b <- " "
      str5 <- "\"Using Data as Observers: a New Paradigm for Prototypes Selection\" "
      str6 <- "M. Herbin, D. Gillard and L. Hussenet "
      str7 <- "Innovations for Community Services, I4CS2016, Vienna, Austria, june 2016"
      str7b <- " "
      str8 <- "Publication of selected papers in
               \"Communications in Computer and Information Science\" "
      str9 <- "Edts: Gunter Fahrnberger, Gerald Eichler, Christian Erfurth"
      str10 <- "Springer, pp 39-46, 2016, DOI 10.1007  978-3-319-49466-1"
      stra <- " "
      strb <- "Acknowledgment: partially supported by the EC SCOOP project
               INEA/CEF/TRAN/A2014/1042281"
      
      HTML(paste(str1, str2, str3, str4, str4b, str5, str6, 
                 str7, str7b, str8, str9, str10, stra, strb, sep = "<br/>" ))
      })
    
    
    output$distribution <- renderPlot({
      n1 <- nb()/2
      n2 <- (nb()/5)*2
      n3 <- nb()-n1-n2
      label1 <- rep(1, n1)
      label2 <- rep(2, n2)
      label3 <- rep(3, n3)
      label <- c(label1, label2, label3)
      data <- donnees()
      colnames(data) <- c("x", "y")
      plot(data, 
           xlim=c(0,20), 
           ylim=c(0, 20),
           pch=c(3, 2, 1)[label], col=c("navy","darkgreen", "red3")[label],
           label = label,
           cex=c(0.8, 0.8, 1.2)[label])
      title(main=paste("Simulation with n = ", nb()))
    })
    
    
    output$skeleton <- renderPlot({
      n <- nb()
      data <- donnees()
      mat.adjacence <- adjacent()

      colnames(data) <- c("x", "y")
      plot(data, 
           xlim=c(0,20), 
           ylim=c(0, 20),
           pch=1, col="navy",
           cex=2.5)
      title(main=paste("Skeleton of the network"))
      text(data, labels=c(1:n), cex= 0.7)
      
      for(i in 1:(n-1)){
        for(j in i:n){
           if(mat.adjacence[i, j] != 0){
            segments(data[i,1], data[i,2], data[j,1], data[j,2])
           }
        }
      }
    })
    
    
    output$sampling <- renderPlot({
      data <- donnees()
      n <- nb()
      mat.adjacence <- adjacent()
      choix <- gains()
      sample <- matrix(data[choix,], length(choix), 2)
      
      colnames(data) <- c("x", "y")
      plot(data, 
           xlim=c(0,20), 
           ylim=c(0, 20),
           pch=16, col="navy",
           cex=0.8)
      
      for(i in 1:(n-1)){
        for(j in i:n){
          if(mat.adjacence[i, j] != 0){
            segments(data[i,1], data[i,2], data[j,1], data[j,2])
          }
        }
      }
      for(i in c(1:length(choix))){
        points(sample[i,1], sample[i,2], lty="solid", pch=19, cex=1.6)
      }
      text(sample, labels=choix, cex= 0.7, pos = 2)
      title(main=paste("Subsampling: A network with", length(choix), " base stations"))
    })
    
    output$number <- renderText({
      paste( "Maximum number of samples:", maxsample(),".")
    })
    
    
    output$network <- renderVisNetwork({
      n <- nb()
      data <- donnees()

      choix <- gains()
      #sample <- matrix(data[choix,], length(choix), 2)
      sensors <- rep("A", n)
      sensors[choix] <- "B"
      taille <- rep(10, n)
      taille[choix] <- 20
      nodes <- data.frame(id = 1:n, group=sensors, label=1:n, font.size = taille*2,
                          title = paste0("<p><b>", 1:n,"</b><br> Sensor </p>"))
      
      adjacence <- adjacent()
      liens <- adjacence*upper.tri(adjacence, diag=FALSE)
      num <- which(liens != 0)
      origine <- as.vector(((num-1)%%n)+1)
      extremite <- as.vector((trunc((num-1)/n))+1)
      distances <- as.matrix(dist(data, diag = TRUE, upper = TRUE, method= "euclidean"))
      distselect <- sprintf("w = %.2f", distances[num])
      edges <- data.frame(from = origine,
                          to = extremite,
                          label = distselect,
                          color = "lightgreen",
                          width = (4*distances[num]))
      #visNetwork(nodes, edges)  with visIgraphLayout() 
      visNetwork(nodes, edges) %>% 
      visGroups(groupname = "A", 
                color = list(background = "lightblue", 
                             border = "darkblue",
                            highlight = "yellow"),
                shape = "circle") %>% 
      visGroups(groupname = "B",
                color = list(background = "lightsalmon", 
                             border = "darkred",
                             highlight = "yellow"),
                shape = "circle")    
    })

  }
)