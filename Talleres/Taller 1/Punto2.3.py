#Taller 1 Punto 2.3 Trayectoria de cohete
#Stiven Gonzalez Olaya
#John Jairo Gonzalez Martinez
#Karen Sofia Coral Godoy
#Daniel Esteban Tibaquira Galindo
import math
import numpy
from matplotlib import pyplot

def f2(x):
  return 4.26*x -0.0052*(x**3)
def f(x):
  y=6+2.13*(pow(x,2))
  y-=(0.0013*(pow(x,4)))
  return y

def metodobiseccion(f, x0, x1, e):
  while x1-x0>=e:
    x2=(x0+x1)/2 
    if f(x2)==0:
      return x2
    else:
      if f(x0)*f(x2)>0: 
        x0=x2
      else:
        x1=x2
  return x2

resultados=[]
x=numpy.arange(-50,50,0.1)
pyplot.xlim(0, 50)
pyplot.ylim(0,1000)
pyplot.plot(x, [f(i) for i in x],"orange")
pyplot.plot(x, [f2(i) for i in x],"red")

print(metodobiseccion(f2,25,35,0.0001))
print(f2(28.622207641601562))
print(f(28.622207641601562))
