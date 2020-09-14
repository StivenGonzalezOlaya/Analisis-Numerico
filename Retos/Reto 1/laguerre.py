# Stiven Gonzalez
# John Gonzalez
# Sofia Coral
# Daniel Tibaquira
# Metodo de Laguerre

from random import random
import cmath

#Metodo de deflación de un polinomio
def deflpoly(a,r):
    n = len(a)-1
    b=[0]*n
    b[0] = a[0]
    for i in range(1,n):
        b[i]=a[i] + r*b[i-1]
    return b
#Metodo para evaluar un polinomio dado el valor x, calcula la primera y segunda derivada
def evaluarPoli(a,x):
    n = len(a)-1
    p = a[0]; dp=0.0; ddp = 0.0
    for i in range(0,n):
        ddp = ddp*x + 2.0*dp;
        dp = dp*x +p
        p = p*x + a[i+1]
    return p,dp,ddp
#Metodo para hallar todas las raices de un polinomio, toma como parametros un polinomio y la tolerancia
def polyroots(a,tol):
    n = len(a)-1
    root =[]
    ite= []
    for i in range(0,n):
        x,it = laguerre(a,tol)
        if(abs(x.imag)<tol):
            x = x.real
        root.append(x)
        ite.append(it)
        a= deflpoly(a,x)
    return root,ite
#Metodo de Laguerre, para hallar la raiz de un polinomio, parametros son un polinomio y la tolerancia
def laguerre(a,tol):
    x=-1
    #x=random()*100+1
    n = len(a)-1
    for i in range(0,100):
        p,dp,ddp = evaluarPoli(a,x)
        if (abs(p)<tol):
            return x,i+1
        g = dp/p; h = g*g -ddp/p
        f = cmath.sqrt((n-1)*(n*h -g*g))
        if (abs(g+f)>= abs(g-f)):
            dx=n/(g+f)
        else:
            dx=n/(g-f)
        x = x -dx
        if(abs(dx)<tol):
            return x,i
    print('Ha superado el numero de iteraciones')


#-------------------------------------------------------------
#Pruebas
#-------------------------------------------------------------
pol = [1,-5,-9,155,-250]
pol1 = [1,-2,4/3,-8/27]
print(polyroots(pol,1e-16))
print(polyroots(pol1,1e-16))




