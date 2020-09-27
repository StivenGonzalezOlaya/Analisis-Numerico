#Taller 1 Punto 4.4  
#Stiven Gonzalez Olaya
#John Jairo Gonzalez Martinez
#Karen Sofia Coral Godoy
#Daniel Esteban Tibaquira Galindo

library(Rmpfr)

steffensen = function(fx, gx, x0, tol){
  error = 1
  k= 0
  Error = c()
  Iter = c()
  while(error >= tol){
    x1 = gx(x0)
    x2 = gx(x1)
    x3 = x0 - (x1-x0)^2/(x2-2*x1+x0)
    x0 = x3
    if (k>0){
      error = abs(x3-x3Ant)/x3Ant  
      Error = c(Error, error)
      Iter = c(Iter, k)
    }
    k = k + 1
    x3Ant = x3
    cat("k = ", k, " ", error,"\n")
    
  }
  titulo = paste("Steffensen x^2-cos(x)\nTolerancia = ", tol)
  plot(y=Error, x=Iter, type="o", xlab="K", ylab="Error K", main=titulo)
  return(x3)
}

fx = function(x){
  x^2-cos(x)
}

gx = function(x){
  cos(x)^(1/2)
}

steffensen(fx, gx, 1, 10e-8)
steffensen(fx, gx, 1, 10e-16)


