#Taller 1 Punto 1.1 Evaluar Polinomio
#Stiven Gonzalez Olaya
#John Jairo Gonzalez Martinez
#Karen Sofia Coral Godoy
#Daniel Esteban Tibaquira Galindo

def horner(nmaximo, X):
    y = 0.0
    total = 0.0
    for x in range(0, nmaximo):
        if x == 0:
            y += 1
            total = 1*(X)
        else:
            y += 2
            total = total + 1*X*(X**x)
    print('Numero de multiplicaciones totales : ', y)
    print(total+1)


#elevado = int (input("Ingrese el mayor numero exponente : "))
elevado = 50
#X = float(input("Ingrese el valor de X a evaluar : "))
X = 1.00000000001
horner(elevado, X)


def evaluarPoli(a, x):
    n = len(a)-1
    p = a[0]
    dp = 0.0
    ddp = 0.0
    for i in range(0, n):
        ddp = ddp*x + 2.0*dp
        dp = dp*x + p
        p = p*x + a[i+1]
    print("Polinomio :", p, "\n1raDerivada :", dp, "\n2daDerivada :", ddp)
    return p


def verificacion(x):
    y = ((x**51)-1)/(x-1)
    print(y)
    return y


pol = [1]*51
x1 = evaluarPoli(pol, X)
x2 = verificacion(X)

error = abs(x1-x2)
porcentaje = (error/x2) * 100
print(porcentaje)

