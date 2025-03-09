import math

# regra de simpsom para f(x) = x^exp
def Simpson_1_3(a, b, N):
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
    
a = -1
b = 3
N = 2

w = [0]*N
t = [0]*N

for i in range(1, N+1, 1):
    if(i <= N/2):
        w[i-1] = i*(b-a)/(2*N)
        t[i-1] = a + (i*w[i-1])/2
    else:
        w[i-1] = w[math.floor(N/2 + (N/2-i))]
        t[i-1] = (a+b) - t[math.floor(N/2 + (N/2-i))]

print(w, t)

# calcular as integrais
print(Simpson_1_3(a, b, N))



