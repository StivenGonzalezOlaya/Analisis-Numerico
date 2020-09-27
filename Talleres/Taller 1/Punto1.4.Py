#Taller 1 Punto 1.4 Decimal a binario 
#Stiven Gonzalez Olaya
#John Jairo Gonzalez Martinez
#Karen Sofia Coral Godoy
#Daniel Esteban Tibaquira Galindo

def decimal_a_binario(parte_decimal):
    bits_permitidos = 10
    binario = []
    
    for i in range (0 , bits_permitidos):
        parte_decimal = parte_decimal * 2
        
        if parte_decimal >= 1:
            parte_decimal = parte_decimal - 1
            binario.append(1)
        else:
            binario.append(0)
    
    return(binario)



# Parte entera del numero decimal convertido a binario:

def entero_a_binario (parte_entera):
    residuo = 0
    numero_residuos = 0
    binario =[]

    while parte_entera > 0:
        residuo = (parte_entera % 2)
        parte_entera = parte_entera // 2
        numero_residuos = numero_residuos + 1
        binario.append(residuo)
        
    numero_binario = []
    j = 1
    for i in range( 0, numero_residuos):
        numero_binario.append(binario [numero_residuos-j])
        j = j+1
    return (numero_binario)

print("11.25 de base 10 a binario: ", entero_a_binario(11) , "." ,decimal_a_binario (0.25), "\n")
print("2/3 a Base 10 a binario: ", "0" , "." , decimal_a_binario (2/3), "\n")
print("30.6 a Base 10 a binario: ", entero_a_binario(30), "." , decimal_a_binario(0.6), "\n")
print("99.9 a Base 10 a binario: ", entero_a_binario(99), "." , decimal_a_binario(0.9), "\n")