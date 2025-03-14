import math
import sympy as sp
from sympy import symbols
import numpy as np

# Calcula os W0 e T0 iniciais com as condições impostas no enunciado
# Retorna uma tupla com os vetores de W0 e T0
def calculo_w0_t0(a, b, N):
    w = [0] * N
    t = [0] * N

    mid = N // 2  # Ponto médio (arredondado para baixo)

    for i in range(N):
        if i < mid:
            w[i] = (i + 1) * (b - a) / (2 * N)
            t[i] = a + (i + 1) * w[i] / 2
        elif i == mid and N % 2 != 0:
            # Caso especial para N ímpar: ponto médio
            w[i] = (b - a) / N
            t[i] = (a + b) / 2
        else:
            # Espelha os valores da primeira metade
            mirrored_index = N - i - 1
            w[i] = w[mirrored_index]
            t[i] = (a + b) - t[mirrored_index]

    return w, t

# Método de Aproximação das Integrais com uma partição de m = 1000
# Retorna um vetor com as integrais de x elevado aos expoentes de 0 a 2N -1 no intervalo [a , b]
def Simpson_1_3_Para_Integral(a, b, N):
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

def GeraMatrizIdentidade(n):
    M = []
    for i in range(n):
        linha = [0] * n
        linha[i] = 1
        M.append(linha)
    return M

def Decomposicao_LU(A):
    U = A.copy()
    L = GeraMatrizIdentidade(A.shape[0])

    for linha in range(0, A.shape[0]):
        for proximaLinha in range(linha+1, A.shape[0]):
          # para cada linha seguinte calcula um multiplicador
          m = U[proximaLinha][linha]/U[linha][linha]
          L[proximaLinha][linha] = m

          # para cada elemento restante da linha
          # atualiza com m * linha anterior
          for coluna in range(linha, A.shape[0]):
            U[proximaLinha][coluna] -= m*U[linha][coluna]

    return L, U

def ResolveSistemaLU(L, U, b):
    # Converte L, U e b para arrays do NumPy
    L = np.array(L, dtype=float)
    U = np.array(U, dtype=float)
    b = np.array(b, dtype=float)

    # Ly = b
    y = np.zeros_like(b)
    for i in range(len(b)):
        y[i] = b[i] - np.dot(L[i, :i], y[:i])
    
    # Ux = y
    x = np.zeros_like(b)
    for i in range(len(b) - 1, -1, -1):
        x[i] = (y[i] - np.dot(U[i, i+1:], x[i+1:])) / U[i, i]
    
    return x

# Monta os Sistemas de Equações 
# Retorna um Vetor com as Expressões
def Definir_Funcoes(a, b, N):
    w_sym = sp.symbols(f'w0:{N}')  
    t_sym = sp.symbols(f't0:{N}')  

    # Calculando os valores de g_j
    g_values = Simpson_1_3_Para_Integral(a, b, N)

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
        
        # Decomposição LU da matriz Jacobiana
        L, U = Decomposicao_LU(J)
        
        # Resolvendo o sistema linear J(wk, tk) * s = -f(wk, tk) usando decomposição LU
        s = ResolveSistemaLU(L, U, -f_values)
        
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

def Funcao_Base(a, b, x):
    # return x**5 => função também utilizada para calcularmos integral no relatório
    return math.exp(a * x + b) # função padrão exigida no relatório

def Quadratura_Gaussiana(N, w, t, a, b, func):
    sum = 0
    for i in range(0, N):
        sum += w[i] * func(a, b, t[i])
    return sum

def Simpson_1_3_Para_Sol_Analitica(a, b, func):
    m = 1000;
    h = (b - a)/m

    c = 1
    integral = 0
    for i in range(0, m+1, 1):
        x = a + i * h

        #ajuste de coeficientes
        if(i == 0 or i==m):
            c=1
        elif(i%2 == 0):
            c=2
        else: c = 4 

        integral += c * func(a, b, x)

    return (h/3) * integral

# conjunto de testes que foram usados para f(x) = x^5
# testes = [
#    [-1, 1, 0.0          ],
#    [-3, 3, 0.0          ],
#    [-1, 3, 121.333333   ]
# ]

# [a, b, solucao_analitica ]
testes = [
    [-1, 1, 6.3890561  ],
    [-3, 3, 54251.5963 ],
    [-1, 3, 53.59815   ]
]

# Realiza testes para a função exp(ax + b)
for a, b, solucao in testes:
    solucao_aproximada = Simpson_1_3_Para_Sol_Analitica(a, b, Funcao_Base)
    for N in range(1, 8):
        print(f"\n\n\nCalculando para o intervalo [{a}, {b}] com N={N}")

        # Executando o método iterativo
        w_final, t_final = Metodo_Newton(a, b, N)

        print("\nValores finais:")
        print("w:", [round(float(wi), 8) for wi in w_final])
        print("t:", [round(float(ti),8) for ti in t_final])

        gauss = Quadratura_Gaussiana(N, w_final, t_final, a, b, Funcao_Base)
        erro = (solucao - gauss)/solucao

        print("Solucao Analitica: ", round(float(solucao), 8))
        print("Solucao Aproximada: ", round(float(solucao_aproximada), 8))
        print("Quadratura de Gauss: ", round(float(gauss), 8))
        print("Erro: ", erro)