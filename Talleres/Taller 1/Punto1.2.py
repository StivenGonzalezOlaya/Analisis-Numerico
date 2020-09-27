#Taller 1 Punto 1.2 Pi en numero binario 
#Stiven Gonzalez Olaya
#John Jairo Gonzalez Martinez
#Karen Sofia Coral Godoy
#Daniel Esteban Tibaquira Galindo

import math
#Pi en numero binario
print("Pi en decimal es :",math.pi)
print ("Pi en numero binario es: ")


parte_entera_pi = round( math.pi,0)
parte_decimal_pi = math.pi-3


residuo = 0
numero_residuos = 0
binario = []


# Parte entera pi
while (parte_entera_pi > 0):
    residuo = round (parte_entera_pi % 2,0)
    parte_entera_pi = parte_entera_pi//2
    binario.append(residuo)
    numero_residuos = numero_residuos + 1


    
numero_binario = []
j = 1
for i in range(0, numero_residuos):
    numero_binario.append(binario [numero_residuos-j])
    j = j+1


  
  
  
binario_decimal = []
# # Parte decimal de pi
for i in range (0,13):
    parte_decimal_pi = parte_decimal_pi * 2
    if(parte_decimal_pi >= 1):
        parte_decimal_pi = parte_decimal_pi - 1
        binario_decimal.append(1) 
    else:
        binario_decimal.append(0)
    

print ("Pi en numero binario de 15 bits es: " , numero_binario, binario_decimal)