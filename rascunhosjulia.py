from sympy import symbols, integrate, solve, Eq, diff

def mudanca_variavel(a, b):
    t = symbols('t')
    x = ((b - a) * t) / 2 + (b + a) / 2
    dx_dt = diff(x, t)
    return x, dx_dt

def sol_exata(f, a, b):
    t = symbols('t')
    x, dx_dt = mudanca_variavel(a, b)
    f_substituida = f.subs(t, x) * dx_dt
    return integrate(f_substituida, (t, -1, 1))

def sol_aproximada(f, pesos, pontos):
    t = symbols('t')
    aprox = sum(pesos[i] * f.subs(t, pontos[i]) for i in range(len(pesos)))
    return aprox

def calcula_pesos_e_pontos(n):
    t = symbols('t')
    w = symbols(f'w0:{n}')  # Pesos w0, w1, ..., w(n-1)
    t_points = symbols(f't0:{n}')  # Pontos t0, t1, ..., t(n-1)

    # Sistema de equações
    equacoes = []
    for k in range(2 * n):
        if k == 0:
            # Equação para a integral de 1 (soma dos pesos)
            eq = sum(w[i] for i in range(n)) - 2
        else:
            # Equações para as integrais de t^k
            eq = sum(w[i] * t_points[i]**k for i in range(n)) - integrate(t**k, (t, -1, 1))
        equacoes.append(Eq(eq, 0))

    # Resolver o sistema de equações
    solucao = solve(equacoes, w + t_points, dict=True)

    if not solucao:
        raise ValueError("Não foi possível encontrar uma solução para o sistema de equações.")

    # Extrair os pesos e pontos da solução
    pesos = [solucao[0][wi] for wi in w]
    pontos = [solucao[0][ti] for ti in t_points]

    return pesos, pontos

N = 3  # Número de pontos de integração
pesos, pontos = calcula_pesos_e_pontos(N)
print("Pesos:", pesos)
print("Pontos:", pontos)

# Definindo a função a ser integrada
t = symbols('t')
f = t**2  # Exemplo: f(t) = t^2

# Intervalo de integração
intervalo = (0, 2)

# Calcular a integral exata
resultado_exato = sol_exata(f, *intervalo)
print("Resultado exato:", resultado_exato)

# Calcular a integral aproximada
resultado_aproximado = sol_aproximada(f, pesos, pontos)
print("Resultado aproximado:", resultado_aproximado)