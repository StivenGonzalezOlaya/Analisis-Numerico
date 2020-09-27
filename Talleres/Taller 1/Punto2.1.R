#Taller 1 Punto 2.1 suma de matriz triangular superior o inferior
#Stiven Gonzalez Olaya
#John Jairo Gonzalez Martinez
#Karen Sofia Coral Godoy
#Daniel Esteban Tibaquira Galindo

matTriangular= function (tamano)
{
  suma= c()
  for (z in tamano)
  {
    mat= matrix(1, nrow = z, ncol = z)
    cont= 0
    numop=0
    for(i in 1:z)
    {
      for (j in 1:i)
      {
        cont = cont + mat[i, j]
        numop= numop+1
      }
    }
    cont= cont-1
    numop= numop-1
    suma = c(suma, cont)
    cont =0 
    cat("\n Matriz ",z, "x", z,"-", numop, "sumas")
  }
  return (suma)
}

val= c(2, 4,8,12,16,20,24, 26, 27 ,28,32,36,40)

#Grafica orden de convergencia O(n^{2})
graph= matTriangular(val)
plot(val, graph, main= "Orden de convergencia O(n^2)", col="blue", xlab = "n", ylab = "Número de sumas", type = 'b')