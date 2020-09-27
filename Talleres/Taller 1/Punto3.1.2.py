#Taller 1 Punto 3.1.2 Convergencia metodos iterativos
#Stiven Gonzalez Olaya
#John Jairo Gonzalez Martinez
#Karen Sofia Coral Godoy
#Daniel Esteban Tibaquira Galindo
import math
import numpy
from matplotlib import pyplot

def res312(f,a,N,E):
    xn = a
    dp = [0]
    for n in range(1,N+1):
        fxn = f(xn)
        if abs(fxn) < E:
            print('Se encontro la solucion en',n,'iteraciones.')
            return xn
        xn = xn - fxn*(xn-dp[n - 1]/f(xn)-f(dp[n-1]))
        dp.append(xn)
    print('Se excedio el limite de iteraciones, no se encontro respuesta')
    return None

f = lambda x: math.log( x + 2 ) - math.sin( x )
res312(f, -1.5, 5000, 1e-16)