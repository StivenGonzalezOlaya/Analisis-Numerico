#Taller 1 Punto 4.1 Convergencia por Aitken 
#Stiven Gonzalez Olaya
#John Jairo Gonzalez Martinez
#Karen Sofia Coral Godoy
#Daniel Esteban Tibaquira Galindo

convergencia = function(s, ini, fin, val){
  ErrorY = c()
  ErrorX = c()
  cat("----------------------------\n")
  cat("Error n\t\tError n+1\n")
  cat("-----------------------------\n")
  for (i in ini+1:fin){
    errorX = abs(val-s(i-1))  
    errorY = abs(val-s(i))
    ErrorX = c(ErrorX, errorX)
    ErrorY = c(ErrorY, errorY)
    cat(errorX, "\t", errorY, "\n")
  }
  cat("\n\n\n")
  plot(x=ErrorX, y=ErrorY, xlab="Error n", ylab="Error n+1", log="x", main="Convergencia sucesi�n", type="o" )
  text(1e-03, 0.1, "y = -0,6286x� + 0,5542x + 0,0003\nR� = 0,997")
}

aitken = function(s, n) {
  Xn = s(1)
  Xn1 = s(2)
  Xn2 = s(3)
  An = 0
  
  DatAn = c()
  DatSn = c()
  DatX = c()
  cat("------------------------------------------------\n")
  cat("Aitken\t\t\tSucesi�n\n")
  cat("------------------------------------------------\n")
  for (i in 4:(n+3)){
    An = Xn2 - ((Xn2-Xn1)^2/(Xn2-2*Xn1+Xn))
    Xn = s(i-2)
    Xn1 = s(i-1)
    Xn2 = s(i)
    DatAn = c(DatAn, An)
    DatSn = c(DatSn, s(i-3))
    DatX = c(DatX, i-3)
    cat(An, "\t", s(i-3), "\n")
  }
  plot(x=DatX, y=DatSn, ylim=c(0.5,1), xlab="", ylab="", main="Aitken vs sucesi�n", col="blue", type="o");
  par(new=TRUE)
  plot(x=DatX, y=DatAn, ylim=c(0.5,1), xlab="n", ylab="S(n)", col="red", type="o");
  legend("bottomright", c("Aitken An","Sucesi�n Xn"),fill=c("red","blue"))
}

s = function(x){
  cos(1/x)  
}

convergencia(s, 1, 100, 1)
aitken(s, 10)


