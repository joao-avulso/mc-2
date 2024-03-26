import math
import matplotlib.pyplot as plt
import struct
sqrt_2 = 1.41421356237

#Retorna o valor de x elevado a y
#x,y => x^y
def pow(x: float, y: int):
    if y == 0:
        return 1
    elif y == 1:
        return x
    elif y<0:
        return 1/pow(x,-y)
    elif y%2 == 0:
        x = pow(x, y/2)
        return x*x
    else:
       return x * pow(x, y-1)


#Retorna o valor de log x na base 2 arredondado para baixo, sendo x maior que 0
#Se o valor for menor ou igual a 0, retorna -1.
# x => log2(x) OU x => -1
def logaritm(x: float):
    if(x >= 1):
        logaritmo = -1
        x = int(x)
        while(x>0):
            x = x >> 1
            logaritmo = logaritmo + 1
        return logaritmo


    if(x > 0):
        logaritmo = -1
        while x != 0:
            x = x + x
            if x >= 1:
                x = 0
            else:
                logaritmo = logaritmo -1
        return logaritmo
    return -1


#Retorna o valor de √2^k mais rapidamente, passando apenas o expoente
#k => √2^k
def sqrt_pow(k: int):
    if k < 0:
        return pow(1/sqrt_2, (-k))
    if k == 0:
        return 1
    if k == 1:
        return sqrt_2
    if (k % 2) == 0:
        return pow(2,k/2)
    if (k % 2) == 1:
        return sqrt_2 * pow(2,((k-1)/2))


#Retorna o valor da raiz de x de maneira rápida
#x => √x
def log_taylor_sqrt(x: float):
    k = logaritm(x)
    f = (x / pow(2,k)) - 1
    aux_f = 1 + (f / 2)*(1 - (f/(4 + (2*f))))

    sqrt_2_pow_k = sqrt_pow(k)

    sqrt_x =  sqrt_2_pow_k * aux_f

    return sqrt_x


def aprox_sqrt(A):
    if A == 0:
        return 0
    P = 536870912
    Q = 8388608
    val_x = struct.unpack('<Q', struct.pack('<d', A))[0]
    val_x += Q
    val_x >> 1
    val_x -= P
    return struct.unpack('<d', struct.pack('<Q', val_x))[0]


def newton_raphson_taylor(A, P, list = []):
    P = 5*pow(10,-P)
    erro = math.inf
    x0 = log_taylor_sqrt(A)
    while(erro > P):
        x1 = 0.5*(x0 + (A/x0))
        erro = abs(x1 - x0)
        x0 = x1
        list.append(x0)
    return list

def newton_raphson_aprox(A, P, list = []):
    P = 5*pow(10,-P)
    erro = math.inf
    x0 = aprox_sqrt(A)
    while(erro > P):
        x1 = 0.5*(x0 + (A/x0))
        erro = abs(x1 - x0)
        x0 = x1
        list.append(x0)
    return list


plt.figure().set_size_inches(10, 10)
# Plot 1
x = [i for i in range(1, 30)]
y = [len(newton_raphson_taylor(i, 16, [])) for i in x]
plt.subplot(2, 1, 1)
plt.yticks(range(1,max(y)+1))
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('Iterações')
plt.title('chute inicial com log e taylor')
plt.grid()

# Plot 2
x = [i for i in range(1, 30)]
y = [len(newton_raphson_aprox(i, 16, [])) for i in x]
plt.subplot(2, 1, 2)
plt.plot(x, y, color = 'r')
plt.xlabel('x')
plt.ylabel('Iterações')
plt.title('chute inicial com algoritmo de aproximação')
plt.grid()

plt.suptitle('Newton-Raphson para √x')
plt.tight_layout()
plt.savefig('../images/newton_raphs.png', dpi=300)
