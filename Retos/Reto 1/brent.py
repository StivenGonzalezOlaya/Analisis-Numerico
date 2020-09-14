# Stiven Gonzalez
# John Gonzalez
# Sofia Coral
# Daniel Tibaquira
# Metodo de Brent
import sys

#Funcion de brent que toma por parametros un rango [a,b], machine elipson, tolerancia y función
def Brent(a,b,machep,t,f):
    sa = a
    sb = b
    fa = f(sa)
    fb = f(sb)
    c = sa
    fc = fa
    e = sb-sa
    d = e
    k=0
    while(True):
        if (abs(fc)<abs(fb)):
            sa=sb;sb=c ;c = sa;fa = fb;fb = fc;fc=fa
        tol = 2.0*machep*abs(sb)+t
        m= 0.5*(c-sb)
        if(abs(m)<=tol or fb == 0 ):break
        if (abs(e)<tol or abs(fa)<=abs(fb)):
            e=m;d=e
        else:
            s = fb/fa
            if (sa==c):
                p = 2.0*m*s; q = 1.0 -s
            else:
                q = fa/fc;r = fb/fc
                p = s*(2.0*m*1.0*(q-r)-(sb-sa)*(r-1.0))
                q = (q-1.0)*(r-1.0)*(s-1.0)
            if (0<p):
                q = -q
            else:
                p = -p
            s=e;e=d
            if((2.0*p<3.0*m*q-abs(tol*q))and (p<abs(0.5*s*q))):
                d = p/q
            else:
                e = m; d = e
        sa = sb; fa=fb
        if(tol<abs(d)):
            sb = sb +d
        elif(0<m):
            sb = sb+tol
        else:
            sb = sb-tol
        fb = f(sb)
        if((0.0<fb and 0.0<fc) or (fb<= 0.0 and fc <= 0.0) ):
            c = sa; fc = fa; e = sb-sa; d=e
        k+=1
    
    value = sb
    print(value,k)
    return value


#-----------------------------------------------------------
#Pruebas
#-----------------------------------------------------------
def fun(x):
    return x**3 - 2*x**2 + 4/3*x - 8/27

def fun1(x):
    return x**4 -5*x**3 -9*x**2 + 155*x - 250

x = Brent(0,1,sys.float_info.epsilon,10e-16,fun)
x1 = Brent(-10,0,sys.float_info.epsilon,1e-8,fun1)

print(x)
print(x1)    
        
