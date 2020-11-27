# *** Analisis Numerico 2020-3 *** 
# -Trabajo realizado por: 
# Daniel Tibaquira Galindo
# Karen Sofia Coral Godoy
# John Gonzalez Martinez
# Stiven Gonzalez Olaya
# -Reto 3. Sistema Epidemologico de COVID en Colombia

library(shiny)
library(shinyjs)
library(deSolve)
ruta = "./Colombia_COVID.csv"
database = read.csv(ruta)
dias=266
print("help")
SeMuestran=TRUE
n = database[266,"TOTAL_CASOS"] + database[266,"TOTAL_MUERTES"] + database[266,"TOTAL_RECUPERADOS"]
print(dias)
print(n)

ui <- fluidPage(
  useShinyjs(),
  titlePanel("COVID-19 Colombia Mar- Nov 2020 "),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("poblacion",
                  "Poblacion:",
                  min = 1,
                  max = n,
                  value = n ),
      sliderInput("infectados",
                  "Infectados:",
                  min = 1,
                  max = n ,
                  value = 1),
      sliderInput("recuperados",
                  "Recuperados:",
                  min = 0,
                  max = n,
                  value = 0),
      sliderInput("transmision",
                  "Tasa transmision:",
                  min = 0,
                  max = 1,
                  value = 0.5),
      sliderInput("recuperacion",
                  "Tasa recuperacion:",
                  min = 0,
                  max = 1,
                  value = 0.1),
      sliderInput("tiempodias",
                  "Tiempo en dias:",
                  min = 1,
                  max = dias,
                  value = dias),
    ),
    
    mainPanel(
      actionButton("botonCalcularSI", "Calcular SI"),
      actionButton("botonCalcularSIR", "Calcular SIR"),
      actionButton("botonReales", "Mostrar Reales"),
      plotOutput("distPlot"),
      plotOutput("distPlot2"),
      plotOutput("distPlotERROR"),
      plotOutput("distPlot2ERROR"),
      plotOutput("distPlot3ERROR"),
    )
  )
)

server <- function(input, output, session) {
  
  
  observeEvent(input$botonCalcularSI, {
    print(class(input$poblacion))
    print(input$poblacion)
    CalcularSI()
  })
  
  observeEvent(input$botonCalcularSIR, {
    print(class(input$poblacion))
    print(input$poblacion)
    calcularSIR()
  })
  observeEvent(input$botonReales, {
    print(class(input$poblacion))
    print(input$poblacion)
    mostrarReales()
  })
  si <- function(times, iniciales, tasa)          #Funcion que calcula el modelo SI
  {
    with(as.list(c(iniciales, tasa)), 
         {
           dS <- - beta*S*I/(S+I)             #Ecuacion diferencial para los susceptibles
           dI <- beta*S*I/(S+I)               #Ecuacion diferencial para los infectados
           return(list(c(dS, dI)))
         })
  }
  
  sir <- function(times, iniciales, tasas)         #Funcion que calcula el modelo SIR
  {
    with(as.list(c(iniciales, tasas)), 
         {
           dS <- - beta*S*I/(S+I+R)           #Ecuacion diferencial para los susceptibles
           dI <- beta*S*I/(S+I+R) - gamma*I   #Ecuacion diferencial para los infectados
           dR <- gamma*I                      #Ecuacion diferencial para los recuperados
           return(list(c(dS, dI, dR)))
         })
  }
  
  CalcularSI <- function(){
    
    iniciales <- c(S = input$poblacion,          #Iniciales contiene en S la poblacion acorde al deslizador  
                   I = input$infectados)         #Iniciales contiene en I los infectados acorde al deslizador
    
    tasa<-c(beta = input$transmision)      #Tasa contiene en Beta la tasa de de transmision acorde al deslizador
    
    output$distPlot3ERROR <- NULL
    
    output$distPlot <- renderPlot({         #Imprime la grafica dada por la solucion del metodo Adams 
      tf = input$tiempodias
      
      
      times <- seq(0, tf, by = 0.1)
      simulacionM1SI.si <- as.data.frame(ode(y=iniciales, times=times, func=si,parms=tasa,method = "adams"))
      attach(simulacionM1SI.si)
      
      plot(times, S, type="l", col="blue", ylim=c(0,sum(iniciales)), xlab="Tiempo (en dias)", ylab="Poblacion",main = "Metodo 1: Adams")
      lines(times, I, type="l", col="red")
      legend(x = "topright", legend=c("Susceptibles", "Infectados"), col=c("blue", "red"), lty=rep(1, 2)) 
    })
    
    output$distPlot2 <- renderPlot({        #Imprime la grafica dada por la solucion del metodo de Euler
      tf = input$tiempodias
      times <- seq(0, tf, by = 0.1)
      simulacionM2SI.si <- as.data.frame(ode(y=iniciales, times=times, func=si,parms=tasa,method = "euler"))
      attach(simulacionM2SI.si)
      plot(times, S, type="l", col="blue", ylim=c(0,sum(iniciales)), xlab="Tiempo (en dias)", ylab="Poblacion", main = "Metodo 2: Euler")
      lines(times, I, type="l", col="red")
      legend(x = "topright", legend=c("Susceptibles", "Infectados"), col=c("blue", "red"), lty=rep(1, 2)) 
      
    })
    
    output$distPlotERROR <- renderPlot({    #Imprime la grafica del error del metodo de Adams contra Euler
      tf = input$tiempodias
      times <- seq(0, tf, by = 0.1)
      simulacionM1SI.si <- as.data.frame(ode(y=iniciales, times=times, func=si,parms=tasa,method = "adams"))
      simulacionM2SI.si <- as.data.frame(ode(y=iniciales, times=times, func=si,parms=tasa,method = "euler"))
      
      i =1
      erroresM1S <- c()
      for (x in simulacionM1SI.si$S) {
        erroresM1S <- c(erroresM1S, ( (abs(simulacionM2SI.si$S[i]-x) )/simulacionM2SI.si$S[i] ) *100)
        i = i +1;
      }
      x <- seq(1,tf)
      plot(simulacionM1SI.si$time , erroresM1S, col="blue", type="l", xlab="Tiempo (dias)", ylab="Error relativo", main = "Error E vs Adams", ylim = c(0,100))
      
    })
    
    output$distPlot2ERROR <- renderPlot({    #Imprime la grafica del error del metodo de Euler contra Adams
      tf = input$tiempodias
      times <- seq(0, tf, by = 0.1)
      simulacionM1SI.si <- as.data.frame(ode(y=iniciales, times=times, func=si,parms=tasa,method = "adams"))
      simulacionM2SI.si <- as.data.frame(ode(y=iniciales, times=times, func=si,parms=tasa,method = "euler"))
      i =1
      erroresM1I <- c()
      for (x in simulacionM1SI.si$I) {
        erroresM1I <- c(erroresM1I, ( (abs(simulacionM2SI.si$I[i]-x) )/simulacionM2SI.si$I[i] ) *100 )
        i = i +1;
      }
      
      x <- seq(1,tf)
      plot(simulacionM1SI.si$time , erroresM1I, col="red", type="l", xlab="Tiempo (en dias)", ylab="Error relativo", main = "Error Adams vs E", ylim = c(0,100))
      
    })
    
  }
  
  calcularSIR <- function(){
    if(SeMuestran){
      SeMuestran=TRUE
      
    }else{
      toggle("distPlotERROR")
      toggle("distPlot2")
      toggle("distPlot2ERROR")
      SeMuestran=TRUE
    }
    
    iniciales <- c(S = input$poblacion,         #Iniciales contiene en S la poblacion acorde al deslizador
                   I = input$infectados,        #Iniciales contiene en I los infectados acorde al deslizador
                   R = input$recuperados)       #Iniciales contiene en R los recuperados acorde al deslizador
    
    tasas<-c(beta = input$transmision,      #Tasa contiene en beta la tasa de transmision acorde al deslizador
             gamma = input$recuperacion)    #Tasa contiene en gamma la tasa de recuperacon acorde al deslizador
    
    output$distPlot <- renderPlot({         #Imprime la grafica dada por la solucion del metodo Adams
      tf = input$tiempodias
      times <- seq(0, tf, by = 0.1)
      simulacionM1SI.sir <- as.data.frame(ode(y=iniciales, times=times, func=sir,parms=tasas,method = "adams"))
      attach(simulacionM1SI.sir)
      
      plot(times, S, type="l", col="blue", ylim=c(0,sum(iniciales)), xlab="Tiempo (en dias)", ylab="Poblacion",main = "Metodo 1: Adams")
      lines(times, I, type="l", col="red")
      lines(times, R, type="l", col="green")
      legend(x = "topright", legend=c("Susceptibles", "Infectados", "Recuperados"), col=c("blue", "red", "green"), lty=rep(1, 2, 3)) 
    })
    
    output$distPlot2 <- renderPlot({        #Imprime la grafica dada por la solucion del metodo de Euler
      tf = input$tiempodias
      times <- seq(0, tf, by = 0.1)
      
      simulacionM2SI.sir <- as.data.frame(ode(y=iniciales, times=times, func=sir, parms=tasas, method = "euler"))
      attach(simulacionM2SI.sir)
      
      plot(times, S, type="l", col="blue", ylim=c(0,sum(iniciales)), xlab="Tiempo (en dias)", ylab="Poblacion", main = "Metodo 2: Euler")
      lines(times, I, type="l", col="red")
      lines(times, R, type="l", col="green")
      legend(x = "topright", legend=c("Susceptibles", "Infectados", "Recuperados"), col=c("blue", "red"), lty=rep(1, 2)) 
      
    })
    
    output$distPlotERROR <- renderPlot({    #Imprime la grafica del error de los infectados del metodo de euler contra la base de datos
      
      tf = input$tiempodias
      times <- seq(0, tf, by = 0.1)
      timesaux <- seq(0, (tf-2), by = 1)
      Steorico = database[["TOTAL_CASOS"]]       
      Iteorico = database[["TOTAL_CASOS"]]
      Rteorico = database[["TOTAL_RECUPERADOS"]]
      Mteorico = database[["TOTAL_MUERTES"]] 
      
      n = database[266,"TOTAL_CASOS"] + database[266,"TOTAL_MUERTES"] + database[266,"TOTAL_RECUPERADOS"]
      
      pobla <- -266:0; 
      pobla[pobla<0] <- n;
      Steorico=pobla-(Rteorico+Iteorico+Mteorico)
      Iteorico=Iteorico[-(tf:length(Iteorico))]
      TeoricoI <- Iteorico
      
      simulacionM2SI.sir <- as.data.frame(ode(y=iniciales, times=timesaux, func=sir,parms=tasas,method = "euler"))
      
      i =1
      erroresM1I <- c()
      for (x in TeoricoI) { 
        erroresM1I <- c(erroresM1I, ( (abs(x-simulacionM2SI.sir$I[i]) )/x ) *100)
        i = i +1;
      }
      
      x <- seq(1,tf)
      plot(simulacionM2SI.sir$time , erroresM1I, col="blue", type="l", xlab="Dias", ylab="Error relativo", main = "Error poblacion Infectada")
      
    })
    
    output$distPlot2ERROR <- renderPlot({   #Imprime la grafica del error de la poblacion susceptible del metodo de Euler contra la base de datos
      tf = input$tiempodias
      times <- seq(0, tf, by = 0.1)
      timesaux <- seq(0, (tf-2), by = 1)
      Steorico = database[["TOTAL_CASOS"]]       
      Iteorico = database[["TOTAL_CASOS"]]
      Rteorico = database[["TOTAL_RECUPERADOS"]]
      Mteorico = database[["TOTAL_MUERTES"]] 
      
      n = database[266,"TOTAL_CASOS"] + database[266,"TOTAL_MUERTES"] + database[266,"TOTAL_RECUPERADOS"]
      
      pobla <- -266:0; 
      pobla[pobla<0] <- n;
      Steorico=pobla-(Rteorico+Iteorico+Mteorico)
      Steorico=Steorico[-(tf:length(Steorico))]
      TeoricoS <- Steorico
      
      simulacionM2SI.sir <- as.data.frame(ode(y=iniciales, times=timesaux, func=sir,parms=tasas,method = "euler"))
      i =1
      erroresM1S <- c()
      for (x in TeoricoS) {
        erroresM1S <- c(erroresM1S, ( (abs(x-simulacionM2SI.sir$S[i]) )/x ) *100 )
        i = i +1;
      }
      x <- seq(1,tf)
      plot(simulacionM2SI.sir$time , erroresM1S, col="red", type="l", xlab="Dias", ylab="Error relativo", main = "Error poblacion Susceptible")
      
    })
    
    output$distPlot3ERROR <- renderPlot({   #Imprime la grafica del error de los recuperados del metodo de Euler contra la base de datos
      tf = input$tiempodias
      times <- seq(0, tf, by = 0.1)
      timesaux <- seq(0, (tf-2), by = 1)
      Steorico = database[["TOTAL_CASOS"]]       
      Iteorico = database[["TOTAL_CASOS"]]
      Rteorico = database[["TOTAL_RECUPERADOS"]]
      Mteorico = database[["TOTAL_MUERTES"]] 
      
      n = database[266,"TOTAL_CASOS"] + database[266,"TOTAL_MUERTES"] + database[266,"TOTAL_RECUPERADOS"]
      
      pobla <- -266:0;   
      pobla[pobla<0] <- n;
      Steorico=pobla-(Rteorico+Iteorico+Mteorico)
      Rteorico=Rteorico[-(tf:length(Rteorico))]
      TeoricoR <- Rteorico
      
      simulacionM2SI.sir <- as.data.frame(ode(y=iniciales, times=timesaux, func=sir,parms=tasas,method = "euler"))
      
      i =1
      erroresM1R <- c() 
      for (x in TeoricoR) {
        erroresM1R <- c(erroresM1R, ( (abs(x-simulacionM2SI.sir$R[i]) )/x ) *100 )
        i = i +1;
      }
      x <- seq(1,tf)
      plot(simulacionM2SI.sir$time , erroresM1R, col="red", type="l", xlab="Dias", ylab="Error relativo", main = "Error poblacion Recuperada")
      
    })
    
  }
  mostrarReales <- function(){
    output$distPlot2 <- NULL
    output$distPlotERROR <- NULL
    output$distPlot2ERROR <- NULL
    output$distPlot3ERROR <- NULL
    
    output$distPlot <- renderPlot({
      tf = input$tiempodias
      times <- seq(0+2, tf, by = 1)
      S = database[["TOTAL_CASOS"]]       
      I = database[["TOTAL_CASOS"]]
      R = database[["TOTAL_RECUPERADOS"]]
      M = database[["TOTAL_MUERTES"]] 
      
      n = database[266,"TOTAL_CASOS"] + database[266,"TOTAL_MUERTES"] + database[266,"TOTAL_RECUPERADOS"]
      
      pobla <- -266:0;   
      pobla[pobla<0] <- n;
      S=pobla-(R+I+M)
      
      S=S[-(tf:length(S))]
      R=R[-(tf:length(R))]
      I=I[-(tf:length(I))]
      M=M[-(tf:length(M))]
      
      plot(times, S, type="l", col="blue", xlab="Tiempo (en dias)", ylab="Poblacion",main = "Datos Reales", ylim=c(0,n))
      lines(times, I, type="l", col="red")
      lines(times, R, type="l", col="green")
      lines(times, M, type="l", col="black")
      legend(x = "topright", legend=c("Susceptibles", "Infectados", "Recuperados", "Muertos"), col=c("blue", "red", "green", "black"), lty=rep(1, 2, 3, 4)) 
    })
    
  }
}

shinyApp(ui = ui, server = server)

#Referencias:
#Tomamos la base de datos de: 
#  
#Obtuvimos informacion de:
#  