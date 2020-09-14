# Stiven Gonzalez
# John Gonzalez
# Sofia Coral
# Daniel Tibaquira
# Metodo de Newton-Horner mejorado con Laguerre

#Libreria que maneja numeros complejos
import cmath

#Metodo de Newton-Horner mejorado
#Pide como parametros, el polinomio, el valor a ser evaluado, maximo numero de iteraciones y la toleracia
def newtonHorner(p,x0,maxIter,tol):
    n = len(p)-1
    k = 0
    x= x0
    while(k<=maxIter ):
        p1,dp,ddp= evaluarPoli(p,x)
        if(abs(dp)<=tol):
            print("Newton's method encountered a slope almost zero.")
        if(p1 ==0):break
        g = dp/p1; h = g*g -ddp/p1
        f = cmath.sqrt((n-1)*(n*h -g*g))
        if (abs(g+f)>= abs(g-f)):
            dx=n/(g+f)
        else:
            dx=n/(g-f)
        x = x -dx
        k+=1
        if(abs(dx)<tol):
            y,q = hornerdefl(p,x)
            return x,k,q
    if(k>maxIter):
        print('Se ha excedido el numero de iteraciones')
    y,q = hornerdefl(p,x)
    return x,k,q

#Metodo de horner para evaluar el polinomio en primera y segunda derivada
def evaluarPoli(a,x):
    n = len(a)-1
    p = a[0]; dp=0.0; ddp = 0.0
    for i in range(0,n):
        ddp = ddp*x + 2.0*dp;
        dp = dp*x +p
        p = p*x + a[i+1]
    return p,dp,ddp

#Se hacer la deflación de horner
def hornerdefl(p,x):
    if(len(p)==0 ): return None
    n = len(p)-1
    q = [0]*(n+1)
    q[0] = p[0]
    for j in range (1,n+1):
        q[j] = p[j] + q[j-1]*x
    return q[-1],q[0:n]


#--------------------------------------------------
#Pruebas
#------------------------------------------------
polinomio =[1,-5,-9,155,-250]
polinomio1 = [1,-2,4/3,-8/27]
x0 = 1
iter=100
tol = 10e-16
k=1
#Se hallan todas las raices del polinomio, no importa si son complejas
while(len(polinomio1)>1 ):
    root,iters,deflate = newtonHorner(polinomio1,x0,iter,tol)
    if(root is not None ):
        print('Valor Raiz',k,':',root, 'Numero iteraciones :',iters)
        polinomio1 = deflate
    else:
        break
    k+=1
