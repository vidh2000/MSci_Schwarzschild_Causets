# numero = input()
# cifra = input()

# conto = 0
# for i in numero:
#     if i == cifra:
#         conto += 1
# print(conto)

lista = []
for i in range(5):
    lista.append(int(input()))

def media(lista):
    s = 0
    for n in lista:
        s += n
    return s/len(lista)
mymedia = media(lista)
print(mymedia)