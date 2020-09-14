


degs = 5
dats = degs

#Funcion Original
fo = function( x ) {
  exp(sin( x ) - cos( x**2 ))
}

#Derivada de la funcion fo
fd = function( x ){
  exp(sin( x ) - cos( x**2 )) * (cos( x ) + 2 * x * sin( x**2 ))
}

#Grafica de la funcion original
plot(fd, xlim = c(-1, 1), ylim = c(-1, 1), col = 'green')

l1 = c()
l2 = c()
ov = (2 * (2 ^ ( -8 ))) / (dats - 1)

#Remez
i = 1
while( i < dats + 1 ){
  l2[i] = fo( l1[i] )
  i = i + 1
}

l2[ i ]= fd( 2 ^ ( -8 ))

j = 1
while(length( l1 ) < dats){
  if (length( l1 ) == 0 ){
    l1[j]= ( -2 ^ ( -8 ))
    j = j + 1
  }
  
  l1[j] = ( l1[j - 1] + ( ov ))
  j = j + 1
} 



#Auxiliar para evaluar monomio
n = rbind(c( 1, l1[1], ( l1[1] ) ^ 2, ( l1[1] ) ^ 3, ( l1[1] ) ^ 4, ( l1[1] ) ^ 5),
          c( 1, l1[2], ( l1[2] ) ^ 2,( l1[2] ) ^ 3,( l1[2] ) ^ 4,( l1[2] ) ^ 5),
          c( 1, l1[3], ( l1[3] ) ^ 2,( l1[3] ) ^ 3,( l1[3] ) ^ 4,( l1 [3] ) ^ 5),
          c(1, l1[4],( l1[4] ) ^ 2,( l1[4] ) ^ 3,( l1[4] ) ^ 4 ,( l1[4] ) ^ 5 ),
          c(1, l1[5],( l1[5] ) ^ 2,( l1[5] ) ^ 3,( l1[5] )^ 4,( l1[5] ) ^ 5),
          c(0, 1, 2 *( 2 ^ ( -8 )), 3 * ( 2 ^ ( -8 )) ^ 2, 4 * ( 2 ^ ( -8 )) ^ 3, 5 * ( 2 ^ ( -8 )) ^ 4)
)

evalSol=solve(n, l2)
print(evalSol)

#Mononimo a evaluar
monom = function( x ) {
  evalSol[1] + ( evalSol[2] * x ) + ( evalSol[3] * x ^ 2 ) + ( evalSol[4] * x ^ 3 ) + ( evalSol[5] * x ^ 4 ) 
}

#Grafica de Aproximacion de Remez
par(new = TRUE)
plot(monom, xlim = c(-1, 1) , ylim = c(-1, 1), col= 'red', main = 'Remez')




#Calculando errores
EA = 0
ER = 0
valPunt = pi/256
EA = abs((fo(valPunt)- monom(valPunt))*10^-6)
ER  = ((EA / monom(valPunt))*100)*10^-6


cat('Utilizando el punto: ',valPunt, '\n')
cat('Error Relativo: ', ER, '\n')
cat('Error Absoluto:' ,EA, '\n')