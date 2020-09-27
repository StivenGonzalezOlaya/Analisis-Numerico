#Taller 1 Punto 3.1.1 Convergencia metodos iterativos
#Stiven Gonzalez Olaya
#John Jairo Gonzalez Martinez
#Karen Sofia Coral Godoy
#Daniel Esteban Tibaquira Galindo
import math
import numpy
from matplotlib import pyplot

def res311(f,a,b,N,E):
    if f(a)*f(b) >= 0:
        print("Secant method fails.")
        return None
    a_n = a
    b_n = b
    for n in range(1,N+1):
        m_n = a_n - f(a_n)*(b_n - a_n)/(f(b_n) - f(a_n))
        f_m_n = f(m_n)
        if(abs(m_n - a) < E): 
            break
        if f(a_n)*f_m_n < 0:
            a_n = a_n
            b_n = m_n
        elif f(b_n)*f_m_n < 0:
            a_n = m_n
            b_n = b_n
        elif f_m_n == 0:
            print("Found exact solution.")
            return m_n, n
        else:
            print("Secant method fails.")
            return None
    return a_n - f(a_n)*(b_n - a_n)/(f(b_n) - f(a_n))

f = lambda x: math.log( x + 2 ) - math.sin( x )
res311(f, 0, -1.7, 20,1e-16)