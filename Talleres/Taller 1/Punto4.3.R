#Taller 1 punto 4.3 
#Stiven Gonzalez Olaya
#John Jairo Gonzalez Martinez
#Karen Sofia Coral Godoy
#Daniel Esteban Tibaquira Galindo

newtonDN = function(f, x0, tol, maxiter){
  # Derivada numérica con diferencia central
  fp = function(x) { h = 1e-15
  (f(x+h) - f(x-h)) / (2*h)
  }
  Error = c()
  Iteraciones = c()
  k = 0
  #Par imprimir estado
  cat("---------------------------------------------------------------------------\n")
  cat(formatC( c("x_k"," f(x_k)","Error est."), width = -20, format = "f", flag = " "), "\
n")
  cat("---------------------------------------------------------------------------\n")
  repeat{
    correccion = f(x0)/fp(x0)
    x1 = x0 - correccion
    dx = abs(correccion)
    # Imprimir iteraciones
    #cat(formatC( c(x1 ,f(x1), dx), digits=15, width = -15, format = "f", flag = " "), "\n")
    cat(formatC( c(x1 ,f(x1), dx), digits=15, width = -15, format = "f", flag = " "), "\n")
    x0 = x1
    k = k+1
    Iteraciones = c(Iteraciones, k)
    Error = c(Error, dx)
    # until
    if(dx <= tol || k > maxiter ) break;
  }
  cat("---------------------------------------------------------------------------\n")
  if(k > maxiter){
    cat("Se alcanzó el máximo número de iteraciones.\n")
    cat("k = ", k, "Estado: x = ", x1, "Error estimado <= ", correccion)
  } else {
    cat("k = ", k, " x = ", x1, " f(x) = ", f(x1), " Error estimado <= ", correccion) }
  
  plot(y=Error, x=Iteraciones, xlab="K", ylab="Error K", type="o", log="y", main="Error vs K")
  return(x1)
}



## --- Pruebas
f1 = function(x) {
  3* (sin(x))^3 - 1
}

f2 = function(x){
  4*sin(x)*cos(x)
}

f = function(x) {
  f1(x) - f2(x) 
}
options(digits = 15)
x1 = newtonDN(f, 1, 1e-16, 100)

plot(f1, xlim=c(0,2), ylim=c(0,2), col="red", xlab="", ylab="", main="f(t) = g(t)")
par(new=TRUE)
plot(f2, xlim=c(0,2), ylim=c(0,2), col="blue", xlab="t", ylab="f(t)")
par(new=TRUE)
#plot(x=c(x1), y=c(f1(x1)),xlim=c(0,2), ylim=c(0,2), col="green", xlab="t", ylab="f(t)", type="o")
abline(v=c(x1), col="gray")

legend("topleft", c("3sin(t)^3-1","4sin(t)cos(t)"),fill=c("red","blue"))


