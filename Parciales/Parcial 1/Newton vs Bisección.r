
#Parcial 1 Metodo de Newton-Raphson 
#Stiven Gonzalez Olaya
#John Jairo Gonzalez Martinez
#Karen Sofia Coral Godoy
#Daniel Esteban Tibaquira Galindo


library(Rmpfr)

biseccion = function(f, a, b, error){
  
  ai = a
  bi = b
  #Error
  E = (b-a)
  
  #Inicialización de variables
  m<-0
  i <-0
  
  #Listas para graficar
  listaErrorAnt = c()
  listaIteraciones = c()
  listaErrorAct = c()
  
  if( a == 0)
  {
    m2 = 1.63848876953125 
  }
  if( a == 1 )
  {
    m2 = 0.55377197265625
  }
  
  while(E > error)
  {
    
    m<-(a+b)/2
    eAnterior = E
    E<-E/2
    
    if (f(a)*f(m) < 0)
    {
      b = m
    }
    
    if (f(b)*f(m) < 0)
    {
      a = m
    }
    
    if(i>0)
    {
      listaErrorAnt = c(listaErrorAnt , eAnterior)
      listaIteraciones = c(listaIteraciones, i)
      listaErrorAct = c(listaErrorAct, E)
    }
    
    if(i == 0)
    {
      #Par imprimir estado
      cat("---------------------------------------------------------------------------\n")
      cat(formatC( c("a","b","m","Error est.","Error ant."), width = -15, format = "f", flag = " "), "\n")
      cat("---------------------------------------------------------------------------\n")
    }
    
    #Imprimir el estado
    cat(formatC( c(a,b,m,E,eAnterior), digits=7, width = -15, format = "f", flag = " "), "\n")
    
    i = i+1
  }
  
  #graficando h(x) vs x
  #plot.function(f,0,2, pch =19, main = "h(x) en función de x",xlab = " X ",ylab = " h(x) ",col ="black")
  
  
  #graficando relación de error
  #plot.function(f, xlim=c(0,1), ylim=c(0,1),  main = "Relación error",xlab = " Error i ",ylab = " Error i+1 ",col ="white")
  
  #imprime la relación de error -> Convergencia lineal
  points(listaErrorAnt, listaErrorAct, col = "blue")
  lines(listaErrorAnt, listaErrorAct, col = "blue")
  
  lstRaices= c(m,m2)
  lst0 = c(0,0)
  #points(lstRaices,lst0,col="red")
  abline(h = 0, v = 0:2/2, lty = 3, col = "gray")
  cat("Cero de f en [",ai,",",bi,"] las raíces son: ", m , "y ", m2 , "con error <=", E, " El numero de iteraciones es: ", i)
  
  listaErrorAnt
  
  
}



newtonDN = function(f, x0, tol, maxiter, funcionT, equ, xx, yy){
  # Derivada num?rica con diferencia central
  fp = function(x) { h = mpfr("1e-20", 128)
    x1 = mpfr(x,128)
    (f(x1+h) - f(x1-h)) / (2*h)
  }
  k = 0
  #Par imprimir estado
  cat("---------------------------------------------------------------------------\n")
  cat(formatC( c("x_k"," f(x_k)","Error est."), width = -20, format = "f", flag = " "), "\
n")
  cat("---------------------------------------------------------------------------\n")
  Ex = c()
  Ey = c()
  Iter = c()
  Dx = c()
  
  repeat{
    correccion = mpfr(f(x0),128)/mpfr(fp(x0),128)
    #cat("valores:",correccion,f(x0),fp(x0),"\n" )
    x1 = mpfr(x0,128) - correccion
    dx = mpfr(abs(correccion),128)
    
    x1A <- substr(capture.output(x1)[2],5,nchar(capture.output(x1)[2]))
    fxA <- substr(capture.output(f(x1))[2],5,nchar(capture.output(f(x1))[2]))
    dxA <- substr(capture.output(dx)[2],5,nchar(capture.output(dx)[2]))
    
    cat(formatC( c(x1A ,fxA, dxA), digits=22, width = -22, format = "f", flag = " "), "\n")
    
    if (k>0){
      
      En1 = En
      En = substr(capture.output(dx)[2],5,nchar(capture.output(dx)[2]))
      Ey = c(Ey, En1)
      Ex = c(Ex, En)
    }
    x0 = x1
    
    k = k+1
    Iter = c(Iter,k)
    Dx = c(Dx,dxA)
    En =  substr(capture.output(dx)[2],5,nchar(capture.output(dx)[2]))
    # until
    if(dx <= tol || k > maxiter ) break;
  }
  cat("---------------------------------------------------------------------------\n")
  x1 <- substr(capture.output(x1)[2],5,nchar(capture.output(x1)[2]))
  correccion <- substr(capture.output(correccion)[2],5,nchar(capture.output(correccion)[2]))
  fxA <- substr(capture.output(f(x1))[2],5,nchar(capture.output(f(x1))[2])) 
  
  if(k > maxiter){
    cat("Se alcanz? el m?ximo n?mero de iteraciones.\n")
    cat("k = ", k, "Estado: x = ", x1, "Error estimado <= ", correccion)
  } else {
    cat(Iter)
    cat("k = ", k, " x = ", x1, " f(x) = ", fxA, " Error estimado <= ", correccion)
    x1
  }
  Ex

}




## --- Pruebas
fun1 = function(x) {
  x= mpfr(x,128)
  cos(2*x)^2 - x^2}
fun2 = function(x) {
  x= mpfr(x,128)
  x*sin(x)-1}
fun3 = function(x){
  x= mpfr(x,128)
  x^3-2*x^2+(4/3)*x-(8/27)}

options(digits = 22)

obtenerEY = function(E1){
  a = length(E1)
  L = c()
  for(i in 1:a){
    L = c(L, i)
  }
  L
}


#Funcion 1 error -8
toler = 1e-8
E1 = newtonDN(fun1, 1, toler, 500, "", 5, 5)
E2 = biseccion(fun1, 0, 1, toler)
EY1 = obtenerEY(E1)
EY2 = obtenerEY(E2)
print(E1)
print(EY1)

limiX = c(0,110)
limiY = c(1e-36, 1 )

titulo = paste("cos(2*x)^2 - x^2", " tolerancia = ", toler)

plot(EY1, E1, xlim=limiX, ylim=limiY, main=titulo, col="red", type="o", ylab="", xlab="", log="y")
par(new=TRUE)
plot(EY2, E2, xlim=limiX, ylim = limiY,  col="blue", type="o",  ylab="Error", xlab="Iteracion n", log="y")
legend("topright", c("Newton","Bisección"),fill=c("blue","red"))
par(new=FALSE)
plot(1,1)

#Funcion 1 error e-16
toler = 1e-16
E1 = newtonDN(fun1, 1, toler, 500, "", 5, 5)
E2 = biseccion(fun1, 0, 1, toler)
EY1 = obtenerEY(E1)
EY2 = obtenerEY(E2)
print(E1)
print(EY1)

limiX = c(0,110)
limiY = c(1e-36, 1 )

titulo = paste("cos(2*x)^2 - x^2", " tolerancia = ", toler)

plot(EY1, E1, xlim=limiX, ylim=limiY, main=titulo, col="red", type="o", ylab="", xlab="", log="y")
par(new=TRUE)
plot(EY2, E2, xlim=limiX, ylim = limiY,  col="blue", type="o",  ylab="Error", xlab="Iteracion n", log="y")
legend("topright", c("Newton","Bisección"),fill=c("blue","red"))
par(new=FALSE)
plot(1,1)



#Funcion 1 error e-32
toler = 1e-32
E1 = newtonDN(fun1, 1, toler, 500, "", 5, 5)
E2 = biseccion(fun1, 0, 1, toler)
EY1 = obtenerEY(E1)
EY2 = obtenerEY(E2)
print(E1)
print(EY1)

limiX = c(0,110)
limiY = c(1e-36, 1 )

titulo = paste("cos(2*x)^2 - x^2", " tolerancia = ", toler)

plot(EY1, E1, xlim=limiX, ylim=limiY, main=titulo, col="red", type="o", ylab="", xlab="", log="y")
par(new=TRUE)
plot(EY2, E2, xlim=limiX, ylim = limiY,  col="blue", type="o",  ylab="Error", xlab="Iteracion n", log="y")
legend("topright", c("Newton","Bisección"),fill=c("blue","red"))
par(new=FALSE)
plot(1,1)




#Funcion 2 error -8
toler = 1e-8
E1 = newtonDN(fun2, 1, toler, 500, "", 5, 5)
E2 = biseccion(fun2, 0, 1, toler)
EY1 = obtenerEY(E1)
EY2 = obtenerEY(E2)
print(E1)
print(EY1)

limiX = c(0,110)
limiY = c(1e-36, 1 )

titulo = paste("x*sin(x)-1", " tolerancia = ", toler)

plot(EY1, E1, xlim=limiX, ylim=limiY, main=titulo, col="red", type="o", ylab="", xlab="", log="y")
par(new=TRUE)
plot(EY2, E2, xlim=limiX, ylim = limiY,  col="blue", type="o",  ylab="Error", xlab="Iteracion n", log="y")
legend("topright", c("Newton","Bisección"),fill=c("blue","red"))
par(new=FALSE)
plot(1,1)


#Funcion 2 error e-16
toler = 1e-16
E1 = newtonDN(fun2, 1, toler, 500, "", 5, 5)
E2 = biseccion(fun2, 0, 1, toler)
EY1 = obtenerEY(E1)
EY2 = obtenerEY(E2)
print(E1)
print(EY1)

limiX = c(0,110)
limiY = c(1e-36, 1 )

titulo = paste("x*sin(x)-1", " tolerancia = ", toler)

plot(EY1, E1, xlim=limiX, ylim=limiY, main=titulo, col="red", type="o", ylab="", xlab="", log="y")
par(new=TRUE)
plot(EY2, E2, xlim=limiX, ylim = limiY,  col="blue", type="o",  ylab="Error", xlab="Iteracion n", log="y")
legend("topright", c("Newton","Bisección"),fill=c("blue","red"))
par(new=FALSE)
plot(1,1)



#Funcion 2 error e-32
toler = 1e-32
E1 = newtonDN(fun2, 1, toler, 500, "", 5, 5)
E2 = biseccion(fun2, 0, 1, toler)
EY1 = obtenerEY(E1)
EY2 = obtenerEY(E2)
print(E1)
print(EY1)

limiX = c(0,110)
limiY = c(1e-36, 1 )

titulo = paste("x*sin(x)-1)", " tolerancia = ", toler)

plot(EY1, E1, xlim=limiX, ylim=limiY, main=titulo, col="red", type="o", ylab="", xlab="", log="y")
par(new=TRUE)
plot(EY2, E2, xlim=limiX, ylim = limiY,  col="blue", type="o",  ylab="Error", xlab="Iteracion n", log="y")
legend("topright", c("Newton","Bisección"),fill=c("blue","red"))
par(new=FALSE)
plot(1,1)



#Funcion 3 error -8
toler = 1e-8
E1 = newtonDN(fun3, 1, toler, 500, "", 5, 5)
E2 = biseccion(fun3, 0, 1, toler)
EY1 = obtenerEY(E1)
EY2 = obtenerEY(E2)
print(E1)
print(EY1)

limiX = c(0,110)
limiY = c(1e-36, 1 )

titulo = paste("x^3-2*x^2+(4/3)*x-(8/27)", " tolerancia = ", toler)

plot(EY1, E1, xlim=limiX, ylim=limiY, main=titulo, col="red", type="o", ylab="", xlab="", log="y")
par(new=TRUE)
plot(EY2, E2, xlim=limiX, ylim = limiY,  col="blue", type="o",  ylab="Error", xlab="Iteracion n", log="y")
legend("topright", c("Newton","Bisección"),fill=c("blue","red"))
par(new=FALSE)
plot(1,1)


#Funcion 3 error e-16
toler = 1e-16
E1 = newtonDN(fun3, 1, toler, 500, "", 5, 5)
E2 = biseccion(fun3, 0, 1, toler)
EY1 = obtenerEY(E1)
EY2 = obtenerEY(E2)
print(E1)
print(EY1)

limiX = c(0,110)
limiY = c(1e-36, 1 )

titulo = paste("x^3-2*x^2+(4/3)*x-(8/27)", " tolerancia = ", toler)

plot(EY1, E1, xlim=limiX, ylim=limiY, main=titulo, col="red", type="o", ylab="", xlab="", log="y")
par(new=TRUE)
plot(EY2, E2, xlim=limiX, ylim = limiY,  col="blue", type="o",  ylab="Error", xlab="Iteracion n", log="y")
legend("topright", c("Newton","Bisección"),fill=c("blue","red"))
par(new=FALSE)
plot(1,1)



#Funcion 3 error e-32
toler = 1e-32
E1 = newtonDN(fun3, 1, toler, 500, "", 5, 5)
E2 = biseccion(fun3, 0, 1, toler)
EY1 = obtenerEY(E1)
EY2 = obtenerEY(E2)
print(E1)
print(EY1)

limiX = c(0,110)
limiY = c(1e-36, 1 )

titulo = paste("x^3-2*x^2+(4/3)*x-(8/27)", " tolerancia = ", toler)

plot(EY1, E1, xlim=limiX, ylim=limiY, main=titulo, col="red", type="o", ylab="", xlab="", log="y")
par(new=TRUE)
plot(EY2, E2, xlim=limiX, ylim = limiY,  col="blue", type="o",  ylab="Error", xlab="Iteracion n", log="y")
legend("topright", c("Newton","Bisección"),fill=c("blue","red"))






