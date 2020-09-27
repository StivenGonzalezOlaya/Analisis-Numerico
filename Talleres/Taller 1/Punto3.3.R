#Taller 1 Punto 3.3 Convergencia metodos iterativos
#Stiven Gonzalez Olaya
#John Jairo Gonzalez Martinez
#Karen Sofia Coral Godoy
#Daniel Esteban Tibaquira Galindo

funcO = expression (exp( x ) - x - 1) 
funcD = D(funcO, "x") 

x = 0
aprox = 1
resultado = c(0)

i =1

while ( x != aprox ) {
  
  x = aprox 
  
  corrPoli = eval( funcO )
  corrDeri = eval( funcD ) 
  
  aprox = x - ( corrPoli / corrDeri ) 
  
  resultado[i] =x
  i = i + 1
}

polinomio =function( x )( exp( x ) - x - 1)
multip = function(y,multiplicidad){
  p = polinomio(0)  
  resultado = ( y - p ) ^ multiplicidad
  if( resultado == p ){    
    cat("Cero de multiplicidad ", multiplicidad)  
  }else{    
    cat("No es un cero de multiplicidad", multiplicidad) 
  }
}
multip(0, 2)


curve(polinomio, -2, 4, 100, ylim = c(0,15), main = "Newton",col ="green")
abline(0,0,col="red")

df = data.frame(resultado)
df