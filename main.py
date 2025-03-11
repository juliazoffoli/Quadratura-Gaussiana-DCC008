import math
import sympy as sp
from sympy import symbols
import numpy as np

# Calcula os W0 e T0 iniciais com as condições impostas no enunciado
# Retorna uma tupla com os vetores de W0 e T0
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

# Método de Aproximação das Integrais com uma partição de m = 1000
# Retorna um vetor com as integrais de x elevado aos expoentes de 0 a 2N -1 no intervalo [a , b]
def Simpsom_1_3(a, b, N):
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

# Monta os Sistemas de Equações 
# Retorna um Vetor com as Expressões
def Definir_Funcoes(a, b, N):
    w_sym = sp.symbols(f'w0:{N}')  
    t_sym = sp.symbols(f't0:{N}')  

    # Calculando os valores de g_j
    g_values = Simpsom_1_3(a, b, N)

    funcoes = []
    for j in range(1, 2*N + 1):  
        expr = 0
        for i in range(N):  
            expr += w_sym[i] * (t_sym[i]**(j-1))  # w_i * t_i^(j-1)
        
        # Subtraindo g_j
        expr -= g_values[j-1]
        funcoes.append(expr)

    return funcoes

# Calculo das Derivadas Parciais por Aproximação QuasiNewton
# Retorna uma Tupla de Matrizes de Aproximações para Dw e Dt, em cada linha é a derivada de uma equação   
# e cada coluna é um Wi com i de 0 a N para Dw 
# e cada coluna é um Ti com i de 0 a N para Tw 

def Aproxima_Derivada_QuasiNewton(w, t, N):
    e = 10e-9
    dxW = [[0] * N for _ in range(2 * N)]
    dxT = [[0] * N for _ in range(2 * N)]

    # para a quantidade de equações 2N o expoente de t varia de 0 a 2N-1
    for power in range(0, 2*N):

        # calcula as derivadas de todos os w e t com k indo de 0 a N
        for k in range(0, N):
            fx = 0
            fxW = 0
            fxT = 0

            # para todos os wi e ti identifica o elemento k que estamos derivando e faz os calculos
            for i in range(0, N):
                # calculo do fx normal 
                Fa = w[i] * pow(t[i], power)
                fx += Fa

                # quando estamos derivando o termo k
                if(k == i): 
                    fxW += (w[i] + e) * pow(t[i], power)
                    fxT += w[i] * pow((t[i] + e), power)

                else: 
                    fxW += Fa
                    fxT += Fa

            # matriz com 2N linhas e N colunas que guarda os N elementos derivados
            # ex: w1 w2 e w2 para a primeira linha devem ser sempre =~ 1
            dxW[power][k] = (fxW - fx) / e
            dxT[power][k] = (fxT - fx) / e

    return (dxW, dxT)

# Monta Matriz Jacobiana de Derivadas Parciais a partir dos Vetores de Aproximações 
# Retorna a Matriz M
def Monta_Matriz(dW, dT, N): 
    M = [[0] * 2*N for _ in range(2*N)]
    for linha in range(0, 2*N):
        for k in range(0, N):  
            M[linha][k] = dW[linha][k]
            M[linha][k+N] = dT[linha][k]
    return M

def Metodo_Newton(a, b, N, TOL=1e-8, max_iter=100):
    w, t = calculo_w0_t0(a, b, N)
    
    for num_iter in range(max_iter):
        funcoes = Definir_Funcoes(a, b, N)
        
        # Calculando os valores de f(wk, tk) com os valores atuais de w e t
        f_values = np.array([float(f.subs({f'w{i}': w[i] for i in range(N)} | 
                                       {f't{i}': t[i] for i in range(N)})) for f in funcoes])
        
        dW, dT = Aproxima_Derivada_QuasiNewton(w, t, N)
        J = np.array(Monta_Matriz(dW, dT, N))
        
        # J(wk, tk) * s = -f(wk, tk)   
        s = np.linalg.solve(J, -f_values) #TODO: implementar metodo de solução de sistemas lineares?
        
        w = [w[i] + s[i] for i in range(N)]
        t = [t[i] + s[i + N] for i in range(N)]
        
        # Norma infinito de s para critério de parada
        norma_s = max(abs(s))
        print(f"Iteração {num_iter + 1}: Norma infinito de s = {norma_s}")
        
        # Critério de parada
        if norma_s < TOL:
            print("Convergência alcançada!")
            break
    else:
        print("Número máximo de iterações atingido sem convergência.")
    
    return w, t

# Executando o método iterativo
a = -1
b = 1
N = 2
w_final, t_final = Metodo_Newton(a, b, N)

print("\nValores finais:")
print("w:", [float(wi) for wi in w_final])
print("t:", [float(ti) for ti in t_final])