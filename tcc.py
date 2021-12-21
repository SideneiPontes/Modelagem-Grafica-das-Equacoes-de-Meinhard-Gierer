from scipy.integrate import solve_ivp #importando do módulo para simular equações diferenciais o integrador solve_ivp
#from scipy.integrate import odeint #importando do módulo para simular equações diferenciais o integrador odeint
import numpy as np #pacote básico da linguagem Python que permite trabalhar, vetores e matrizes de N dimensões (parece matlab)
import random

# N Conjunto de equações acopladas
# vetor de N componentes, sendo i a posição no vetor e u(i,t) = x(i,t), u(i+N,t) = y(i,t) , ..., u(i+(Ne-1)*N,t)=z(i,t)
# Equações de Meinhardt Gierer [dissertação Fabio]
# dudt[i]=ρx*(u[i]**2)/u[i+N]-μx*u[i]
# dudt[i+N] = ρy*(u[i]**2)-μy*u[i+N]


#________________________________________________Definindo parâmetros__________________________________________________________
#_______Parâmetros da integração______
N=101 #número de sítios
Ne=2 #qtd. de eqs. de 1º ordem em cada sítio, se mexer aqui também de incluir eq. nas condi. de contorno e acoplamento
alpha=1
vizinhotimos=1
tinicial=0
tfinal=3
step=1
Npontos= int(tfinal/step) #quantidade de tempos iterados para obter a solução, que é associado ao número de linhas da matriz solução
tempo = np.linspace(tinicial, tfinal, Npontos)

#______________________________________________________________________________________________

# _______Parâmetros das equações_______
μx = 0.01
μy = 0.02
ρx = 0.01
ρy = 0.02
D0 = 0.005
D1 = 0.2
# _____________________________________
# ________________________________________________Condições Iniciais____________________________________________________________
u0 = np.zeros(Ne * N)  # define um vetor vazio com Ne*N componentes, para as c.i.
for i in range(0, Ne * N, 1):
    #    u0[i]=random.random() # c.i. randômicas no intervalo 0 a 1.
    u0[i] = 0.9  # u[i]=1.0            # u[i]=0.9            u[i]=1          u[i]=1.0        u[i]=1
u0[148] = 0.95  # u[148]=0.90         # u[148]=0.95         u[148]=1.05      u[48]=0.95      u[48]=0.90
u0[149] = 0.99  # u[149]=0.90         # u[149]=0.99         u[149]=1.10      u[49]=0.99      u[49]=0.95
u0[150] = 0.95  # u[150]=0.95         # u[150]=0.95         u[150]=1.05      u[50]=0.95      u[50]=0.90

# imprimindo as condições iniciais
# print("Quantidade de Equações" ,len(u0))
# for iteste in range(0, Ne*N, 1) :
#    print( "indice %d = %f"% (iteste, u0[iteste]) ) #depois do % se aparecer: d=>inteiro e f=float

# ___________________Constante de Normalização e Otimização___________________________
# Otimização para casos de predominancia difusiva
if alpha > 5:
    Nvizinhos = int(vizinhotimos)
else:
    Nvizinhos = int((N - 1) / 2)

# Constante de Normalização
norma = 0
for j in range(1, Nvizinhos + 1, 1):
    norma = norma + 2 / ((j) ** (alpha))


# _____________________________________________________________________________________
# ________________________________________________Definindo Equações____________________________________________________________
# def fun (u,tempo): #definição para Odeint
def fun(t, u):  # definição para solve_ivp
    dudt = np.zeros(Ne * N)  # define um vetor vazio, para as derivadas

    # ________________________Condições de contorno periódicas____________________________
    contor0 = np.zeros(3 * N)  # vetor para condições de contorno peródicas da 1° variável
    contor1 = np.zeros(3 * N)  # vetor para condições de contorno peródicas da 2° variável
    for i0 in range(0, N, 1):
        contor0[i0] = u[i0]
        contor0[i0 + N] = u[i0]
        contor0[i0 + 2 * N] = u[i0]
        contor1[i0] = u[i0 + N]
        contor1[i0 + N] = u[i0 + N]
        contor1[i0 + 2 * N] = u[i0 + N]

    # _____________________________Cálculo do acoplamento lei de potência________________________
    acopl0 = np.zeros(N)  # vetor de acoplamento da 1° variável
    acopl1 = np.zeros(N)  # vetor de acoplamento da 2° variável
    for i1 in range(0, N, 1):
        for i2 in range(1, Nvizinhos + 1, 1):
            acopl0[i1] = (contor0[N + i1 - i2] + contor0[N + i1 + i2]) / (i2 ** alpha) + acopl0[i1]
            acopl1[i1] = (contor1[N + i1 - i2] + contor1[N + i1 + i2]) / (i2 ** alpha) + acopl1[i1]

        # _____________________Entrando com as equações__________________________
    for i in range(0, N, 1):
        dudt[i] = ρx * (u[i] ** 2) / u[i + N] - μx * u[i] + (D0) * (acopl0[i] / norma - u[i])
        dudt[i + N] = ρy * (u[i] ** 2) - μy * u[i + N] + (D1) * (acopl1[i] / norma - u[i + N])
    # _______________________________________________________________________

    return dudt

#___________________________________________________________________________________________________________-


#______________________________Usando o odeint______________________________________________

#sol.y = odeint(fun, u0, tempo, rtol = 1.0e-10, atol = 1.0e-10) # lembrar de mudar o def fun
# Na matriz solução do odeint, as linhas representariam o tempo e as colunas os sítios.
#___________________________________________________________________________________________




#____________________________Usando o solve_ivp_____________________________________________
sol=solve_ivp(fun=lambda t, y: fun(t, y), t_span=[tinicial,tfinal], y0=u0, t_eval=tempo, method = 'RK45', rtol=1.0e-10, atol=1.0e-10) # lembrar de mudar o def fun
# Na matriz solução do solve_ivp, sol.y, as linhas são os sítios e as colunas representam o tempo.

sol.y=np.transpose(sol.y) #Para ficar igual a saída do odeint. As linhas serão os tempos e as colunas os sítios.

#_________________________________________________________________________________________________________________

import matplotlib.pyplot as plt #modulo para geração de gráficos
plt.figure(1)
plt.title("Séries Temporais")
plt.plot(tempo, sol.y[:, 0], 'b', label='X1')
plt.plot(tempo, sol.y[:, 1], 'g', label='X2')
plt.plot(tempo, sol.y[:, 2], 'r', label='X3')
plt.plot(tempo, sol.y[:, 48], 'y', label='X49')
plt.legend(loc='best')
plt.xlabel('t')
plt.axis([0,tfinal, 0, 2])
plt.grid()
plt.show()

#Perfil Espacial
plt.figure(2)
plt.title("Perfil Espacial")
plt.plot(sol.y[Npontos-1,0:N], linestyle='--', marker='o', color='b')#plota os valores da 1º variável da ultima linha do vetor solução, que é o ultimo tempo de integração, assim imprime o perfil espacial da rede.
plt.xlim(0,N)
#plt.ylim(0.85,1.05)
plt.xlabel("k (sítios)")
plt.ylabel("X0_final")
plt.show()

plt.figure(3)
plt.title("Perfil Espacial")
plt.plot(sol.y[Npontos-1,N:2*N+1])#plota os valores da 2º variável da ultima linha do vetor solução, que é o ultimo tempo de integração, assim imprime o perfil espacial da rede.
plt.xlim(0,N)
#plt.ylim(0.85,1.05)
plt.xlabel("k (sítios)")
plt.ylabel("Y0_final")
plt.show()

#___________________________________________________________________________________________________________