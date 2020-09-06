
#Parcial 1 Metodo de Newton-Raphson 
#Stiven Gonzalez Olaya
#John Jairo Gonzalez Martinez
#Karen Sofia Coral Godoy
#Daniel Esteban Tibaquira Galindo

library(Rmpfr)

# Implementación del Metodo de Newton-Raphson 
newtonDN = function(f, x0, tol, maxiter){
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
  #Creación de vectores para graficar
  Iter = c()
  Dx = c()
  
  repeat{
    correccion = mpfr(f(x0),128)/mpfr(fp(x0),128)
    #cat("valores:",correccion,f(x0),fp(x0),"\n" )
    x1 = mpfr(x0,128) - correccion
    dx = mpfr(abs(correccion),128)
    # Se substraen los valores de cada variable
    x1A <- substr(capture.output(x1)[2],5,nchar(capture.output(x1)[2]))
    fxA <- substr(capture.output(f(x1))[2],5,nchar(capture.output(f(x1))[2]))
    dxA <- substr(capture.output(dx)[2],5,nchar(capture.output(dx)[2]))
    
    cat(formatC( c(x1A ,fxA, dxA), digits=22, width = -22, format = "f", flag = " "), "\n")
    
    if (k>1){
      
      En1 = En
      En = substr(capture.output(dx)[2],5,nchar(capture.output(dx)[2]))
      
    }
    x0 = x1
    
    k = k+1
    Iter = c(Iter,k)
    Dx = c(Dx,dxA)
    En =  substr(capture.output(dx)[2],5,nchar(capture.output(dx)[2]))
    # until
    if(dx <= tol || k > maxiter ) break;
  }
  #Impresion de valores resultado final
  cat("---------------------------------------------------------------------------\n")
  x1 <- substr(capture.output(x1)[2],5,nchar(capture.output(x1)[2]))
  correccion <- substr(capture.output(correccion)[2],5,nchar(capture.output(correccion)[2]))
  fxA <- substr(capture.output(f(x1))[2],5,nchar(capture.output(f(x1))[2])) 
  
  if(k > maxiter){
    cat("Se alcanz? el m?ximo n?mero de iteraciones.\n")
    cat("k = ", k, "Estado: x = ", x1, "Error estimado <= ", correccion)
  } else {
    cat("k = ", k, " x = ", x1, " f(x) = ", fxA, " Error estimado <= ", correccion)
    x1
  }
  #Grafica de iteraciones vs tolerancia
  plot(Iter,Dx,type="o",log="y", ylab ="Tolerancia" ,xlab = "N-iteraciones")
  
  
  x1

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

# -----------------------------------
# Resultados funcion 1 para cada raiz
# -----------------------------------

x1f8 = newtonDN(fun1, 1, 10e-8, 500)
x2f8 = newtonDN(fun1, -1, 10e-8, 500)
x1f16 = newtonDN(fun1, 1, 10e-16, 500)
x2f16 = newtonDN(fun1, -1, 10e-16, 500)
x1 = newtonDN(fun1, 1, 10e-32, 500)
x2 = newtonDN(fun1, -1, 10e-32, 500)

#Grafica de la primera función
plot( fun1, -1, 2, main="cos(2*x)^2 - x^2" , ylab ="f(x)" ,xlab = "X")

abline(h=c(0), v=c(0))
points(x1, 0, pch=19)
abline(v=c(x1), col="red")
points(x2, 0, pch=19)
abline(v=c(x2), col="red")

# -----------------------------------
# Resultados funcion 2
# -----------------------------------

x3f8 = newtonDN(fun2, 1, 10e-8, 500)
x3f16 = newtonDN(fun2, 1, 10e-16, 500)
x3 = newtonDN(fun2, 1, 10e-32, 500)

#Grafica de la función
plot( fun2, -1, 2, main="x*sin(x)-1" , ylab ="f(x)" ,xlab = "X")

abline(h=c(0), v=c(0))
points(x3, 0, pch=19)
abline(v=c(x3), col="red")

# -----------------------------------
# Resultados funcion 3
# -----------------------------------
x4f8 = newtonDN(fun3, 1, 10e-8, 500)
x4f16 = newtonDN(fun3, 1, 10e-16, 500)
x4 = newtonDN(fun3, 1, 10e-32, 500)

#Grafica de la función
plot( fun3, -1, 2, main="x^3-2*x^2+(4/3)*x-(8/27)" , ylab ="f(x)" ,xlab = "X")

abline(h=c(0), v=c(0))
points(x4, 0, pch=19)
abline(v=c(x4), col="red")


