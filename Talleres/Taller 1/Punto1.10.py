#Taller 1 Punto 1.10 Decimal a binario 
#Stiven Gonzalez Olaya
#John Jairo Gonzalez Martinez
#Karen Sofia Coral Godoy
#Daniel Esteban Tibaquira Galindo

def decimal_a_binario(parte_decimal):
    bits_permitidos = 10
    hexa = []
    
    for i in range (0 , bits_permitidos):
        parte_decimal = parte_decimal * 16
        entero = int(parte_decimal)
        

        if parte_decimal-entero != 0:
            parte_decimal = parte_decimal - entero
            
            hexa.append(hex(entero))
        else:
            hexa.append(hex(0))
    
    return(hexa)

print(hex(9) ,'.',decimal_a_binario(0.4))

