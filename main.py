import math

def calculo_w0_t0(a, b, N):
    w = [0]*N
    t = [0]*N

    for i in range(1, N+1, 1):
        if(i <= N/2):
            w[i-1] = i*(b-a)/(2*N)
            t[i-1] = a + (i*w[i-1])/2
        else:
            w[i-1] = w[math.floor(N/2 + (N/2-i))]
            t[i-1] = (a+b) - t[math.floor(N/2 + (N/2-i))]

    return (w, t)

def g(a,b,N):
    m = 1000;
    h = (b - a)/m

    result = [0]*2*N
    c = 1

    # para cada expoente de 0 até 2N-1
    for exp in range(0, 2*N, 1):
        integral = 0

        # em intervalo simétrico e função par a integral vale 0
        # exp +1 é o valor a ser integrado
        if(b == -a and (exp+1)%2 == 0):
            result[exp] = 0
            continue;

        for i in range(0, m+1, 1):
            x = a + i * h

            #ajuste de coeficientes
            if(i == 0 or i==m):
                c=1
            elif(i%2 == 0):
                c=2
            else: c = 4 

            integral += c * pow(x, exp)
        result[exp] = (h/3) * integral
    return result;

def definir_funcoes(w, t, a, b, N):
    g_values = g(a, b, N)
    funcoes = [0]*(2*N)

    for j in range(1, 2*N + 1):
        for i in range(N):
            funcoes[j - 1] += w[i]*(t[i]**(j-1))
        funcoes[j-1] -= g_values[j-1]

    return funcoes
a = -1
b = 1
N = 4

w, t = calculo_w0_t0(a, b, N)
print("w:", w)
print("t:", t)

funcoes = definir_funcoes(w, t, a, b, N)
print("Funções f_j:", funcoes)