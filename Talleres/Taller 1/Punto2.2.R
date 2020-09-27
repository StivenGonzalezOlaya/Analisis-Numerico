#Taller 1 Punto 2.1 Suma de n^{2} numeros naturales
#Stiven Gonzalez Olaya
#John Jairo Gonzalez Martinez
#Karen Sofia Coral Godoy
#Daniel Esteban Tibaquira Galindo
sumaCuadrados= function (n)
{
  suma= c()
  for (z in n)
  {
    cont= 0
    for(i in 1:z)
    {
      cont= cont+ i^2
    }
    cat("\nn=",z, "- suma=", cont)
    suma = c(suma, cont)
    cont =0 
    
  }
  return (suma)
}

val= c(2,4,5,8,12,16,20,24,28,32,36,40)

#Grafica orden de convergencia O(n^{3})
graph= sumaCuadrados(val)
plot(val, graph, main= "Orden de convergencia O(n^3)", col="red", xlab = "n", ylab = "Suma de n^2", type = 'b')