# Stiven Gonzalez
# John Gonzalez
# Sofia Coral
# Daniel Tibaquira
# Metodo de Horner para hallar la primera y segunda derivada de un polinomio

def evaluarPoli(a,x):
    n = len(a)-1
    p = a[0]; dp=0.0; ddp = 0.0
    for i in range(0,n):
        ddp = ddp*x + 2.0*dp
        dp = dp*x +p
        p = p*x + a[i+1]
    print("Polinomio :",p,"\n1raDerivada :",dp,"\n2daDerivada :",ddp)

# --------------------------------------------------------------------
# Pruebas del metodo
# ---------------------------------------------------------------------

pol1=[1,-5,-9,155,-250]
pol2=[1,0,-4]
pol3=[2,0,-3,3,-4]
print('Polinomio 1')
evaluarPoli(pol1,1)
evaluarPoli(pol1,2)
print('Polinomio 2')
evaluarPoli(pol2,1)
evaluarPoli(pol2,2)
print('Polinomio 3')
evaluarPoli(pol3,2-3j)
evaluarPoli(pol3,2)


