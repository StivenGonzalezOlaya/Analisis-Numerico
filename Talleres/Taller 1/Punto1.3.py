# Parte entera del binario convertido a base 10
def binario_a_entero(parte_entera):
    numero_potencia = 0
    binario = 0
    total = 0
    while parte_entera > 0:
        binario = parte_entera % 10
        total = total + (binario*(2 ** numero_potencia))
        numero_potencia = numero_potencia + 1
        parte_entera = parte_entera // 10
    return total
  

# Parte decimal del binario convertida a base 10

def binarioDecimal(numero):
    total=0
    for i in range (0,len(numero)):
        x=int(numero[i])
        total = total + (x * 2**(-1*(i+1)))    
    return total

print("Primer Numero: ", binario_a_entero (1010101),"\n")

print("Segundo Numero: ",binario_a_entero(1011) + binarioDecimal('101'),"\n")

print("Tercer Numero: ",binario_a_entero(10111) + binarioDecimal('010101010101010101'),"\n")

print("Cuarto Numero: ",binario_a_entero(111) + binarioDecimal('11111111111111111111111'))






