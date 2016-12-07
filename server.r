library(shiny)
require(mvtnorm)
require(visNetwork)

shinyServer(
  function(input, output){
    nb <- reactive(input$nb)
    maxsample <- reactive(input$maxsize)
    data <- NULL
    
    donnees <- reactive({
      n1 <- nb()/2
      n2 <- (nb()/5)*2
      n3 <- nb()-n1-n2
      source <- 215
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
    
    
    adjacente <- reactive({
      n <- nb()
      data <- donnees()
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
        # retourne vrai=1 s'il n'y a exclusivement rien entre i et j
        return( !(sum(avant_i & avant_j)) ) 
      }
      # lien logique (true ou false)
      lien <- matrix(mapply(entre, pointdevues, row(pointdevues), col(pointdevues)), n)
      matrix(as.numeric(lien)-diag(1,n), n)
    })
    
    
    
    gains <- reactive({
      n <- nb()
      d <- 2
      data <- donnees()
      op <- c(1:n)
      echantillon <- data
      maxsize <- as.numeric(input$maxsize)
      sizesample <- n
      
      while(sizesample > maxsize){
        
        distances <- as.matrix(dist(echantillon, diag = TRUE, upper = TRUE, method= "euclidean"))
        pointdevues<-array(rep(0, n*n), c(n, n))
        for(individu in 1:n){
          pointdevues[individu, ] <- rank(distances[individu, ], ties.method= "random")
        }
        entre=function(x, i, j){
          ligne_i <- pointdevues[i,]
          ligne_j <- pointdevues[j,]
          avant_j <-  (1<ligne_i) & (ligne_i<pointdevues[i,j])
          avant_i <-  (1<ligne_j) & (ligne_j<pointdevues[j,i])
          # retourne vrai=1 s'il n'y a exclusivement rien entre i et j
          return( !(sum(avant_i & avant_j)) ) 
        }
        lien <- matrix(mapply(entre, pointdevues, row(pointdevues), col(pointdevues)), n)
        mat.adjacence <- as.numeric(lien) - diag(1,n)
      
        gain<-array(rep(0, n*n), c(n, n))
        for(i in 1:n){
          for(j in 1:n){
            if(mat.adjacence[i, j] != 0){
            # nombre de liens strictement entrants sur i pour le lien [i,j]
              gain[i, j]<-sum(mat.adjacence[,i])
            }
          }
        }
        perte<-array(rep(0, n*n), c(n, n))
        for(i in 1:n){
          for(j in 1:n){
            if(mat.adjacence[i, j] != 0){
              # nombre de liens strictement sortants de j pour le lien [i,j]
              perte[i, j]<-sum(mat.adjacence[j,])
            }
          }
        }
        mat.orient<-mat.adjacence
        for(i in 1:n){
          for(j in 1:n){
            # il reste des liens non orientes en cas d'egalite
            if(gain[i,j] < perte[i,j]){
              # orienter vers le debit le plus fort
              # l'autre sens devient nul
              mat.orient[j, i]<-0
            }
          }
        }
        # orientation : 2eme passe
        for(i in 1:n){
          for(j in 1:n){
            # il reste des liens non orientes en cas d'egalite
            if((mat.orient[i, j] != 0) & (gain[i,j] == perte[i,j])){
              # remonter en amont pour l'origine
              # tous les liens orientes ou non, d'extremite i
              coord <- which(mat.orient[,i] != 0, arr.ind=TRUE)
              gainbis <- max(gain[coord,i])
              # descendre en aval pour l'extremite
              # tous les liens d'origine j
              coord <- which(mat.orient[j,] != 0, arr.ind=TRUE)
              pertebis<- max(perte[j,coord])
              if(gainbis > pertebis){
                mat.orient[i, j]<-0
              } else {
                mat.orient[j, i] <- 0
              }
            } 
          }
        }
        op<-NULL
        for(i in 1:n){ 
          # pas de sortant mais au moins un entrant
          if( (sum(mat.orient[i,]) == 0) & (sum(mat.orient[,i]) != 0) ){
            op<-c(op, i)
          }
        }
        sizesample <- length(op)
        if(sizesample > 1){
          n <- length(op)
          data <- echantillon
          echantillon <- data[op,]
        }
      }
      
      sample <- NULL
      if(sizesample == 1) {
        sample <- matrix(echantillon[op,], 1, 2)
      } else {
        if(sizesample == 0) {
          sample <- donnees()
        } else {
          sample <- echantillon
        }
      }
      sample
    })
    
    
    output$readme <- renderUI({
      str1 <- h4("Determination of a subsample of a 2D dataset")
      str2 <- paste( "The ", as.character(nb()), " points of a 2D dataset 
                simulate sensors of a Wireless Sensors Network (WSN). 
                Sub-samples could determine the location of the base stations in the network.")
      str3 <- " "
      str4 <- "Ref.:"
      str5 <- "\"Using Data as Observers: a New Paradigm for Prototypes Selection\" "
      str6 <- "M. Herbin, D. Gillard and L. Hussenet"
      str7 <- "Innovations for Community Services (I4CS),  
                 pp 74-82, Springer, Vienna, Austria, june 2016. "
      HTML(paste(str1, str2, str3, str4, str5, str6, str7, sep="<br/>"))
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
           cex=c(0.8, 0.8, 1.2)[label])
      title(main=paste("Simulation with n = ", nb()))
    })
    
    
    output$skeleton <- renderPlot({
      n <- nb()
      data <- donnees()
      mat.adjacence <- adjacente()

      colnames(data) <- c("x", "y")
      plot(data, 
           xlim=c(0,20), 
           ylim=c(0, 20),
           pch=16, col="navy",
           cex=0.8)
      title(main=paste("Skeleton of the network"))
      
      for(i in 1:(n-1)){
        # non oriente : traitement diagonal superieure
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
      mat.adjacence <- adjacente()
      sample <- gains()
      
      colnames(data) <- c("x", "y")
      plot(data, 
           xlim=c(0,20), 
           ylim=c(0, 20),
           pch=16, col="navy",
           cex=0.8)
      
      for(i in 1:(n-1)){
        # non oriente : traitement diagonal superieure
        for(j in i:n){
          if(mat.adjacence[i, j] != 0){
            segments(data[i,1], data[i,2], data[j,1], data[j,2])
          }
        }
      }
      
      for(i in c(1:dim(sample)[1])){
        points(sample[i,1], sample[i,2], lty="solid", pch=19, cex=1.5)
      }
      
      title(main=paste("Subsampling: A network with", dim(sample)[1], " base stations"))
    })
    
    output$number <- renderText({
      paste( "Maximum number of samples:", maxsample(),".")
    })
    
    
    output$network <- renderVisNetwork({
      n <- nb()
      data <- donnees()
      
      nodes <- data.frame(id = 1:n, shape="circle", label=c(1:n) )
      
      mat.adjacence <- as.vector(adjacente())
      num <- which(mat.adjacence != 0)
      x <- as.vector(((num-1)%%n)+1)
      y <- as.vector((trunc((num-1)/n))+1)
      origine <- x[x>y]
      extremite <- y[x>y]
      distances <- as.vector(dist(data, diag = TRUE, upper = TRUE, method= "euclidean"))
      dist <- distances[num]
      edges <- data.frame(from = origine, 
                          to = extremite)
      visNetwork(nodes, edges)
    })

  }
)