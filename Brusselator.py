#!/usr/bin/env python
# coding: utf-8

# # Brusselator - Marlon
# 
# Propagação de ordem exponencial na rede

# In[1]:


#from scipy.integrate import odeint #importando do módulo para simular equações diferenciais o integrador odeint
import numpy as np #pacote básico da linguagem Python que permite trabalhar, vetores e matrizes de N dimensões (parece matlab)
import matplotlib.pyplot as plt #modulo para geração de gráficos
from math import exp
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
#from scipy.integrate import solve_ivp
import pandas as pd
from scipy.integrate import solve_ivp


# In[161]:


#________________________________________________Definindo parâmetros__________________________________________________________

N=101 #número de sítios
Ne=2 #qtd. de eqs. de 1º ordem em cada sítio, se mexer aqui também de incluir eq. nas condi. de contorno e acoplamento
gamma=10
vizinhotimos=5

D0 = 5 
D1 = 8 

a  = 4
b =  24

tinicial=0
tfinal=50
step=0.02
iteracoes = int(tfinal/step)

t = np.linspace(tinicial, tfinal,(iteracoes)+1)

print('iterações: ',iteracoes)


# In[173]:


#________________________________________________Condições Iniciais____________________________________________________________     
u0 = np.zeros(Ne*N) # define um vetor vazio com Ne*N componentes, para as c.i.
for i in range(0, N, 1) : 
    
    u0[i] = a
    u0[i+N] = b/a   
    
u0[48] = a + 1.0
u0[49] = a + 2.0
u0[50] = a + 1.0


#Otimização para casos de predominancia difusiva
if gamma>5:
    Nvizinhos=int(vizinhotimos)
else:
    Nvizinhos=int((N-1)/2)

#Constante de normalização
norma=0
for j in range(1,Nvizinhos+1, 1) : 
    norma=norma+2/(exp(j*gamma))
#___________________________________


# In[174]:


#________________________________________________Definindo Equações____________________________________________________________
def fun(t, u):
    dudt = np.zeros(Ne*N) # define um vetor vazio, para as derivadas
    
    
#________________________Condições de contorno periódicas____________________________    
    contor0= np.zeros(3*N) # vetor para condições de contorno peródicas da 1° variável
    contor1= np.zeros(3*N) # vetor para condições de contorno peródicas da 2° variável
    for i0 in range(0,N,1):
        contor0[i0]=u[i0]   
        contor0[i0+N]=u[i0]
        contor0[i0+2*N]=u[i0]
        contor1[i0]=u[i0+N]   
        contor1[i0+N]=u[i0+N]
        contor1[i0+2*N]=u[i0+N]

#_____________________________Cálculo do acoplamento lei de potência________________________  
    acopl0= np.zeros(N) # vetor de acoplamento da 1° variável
    acopl1= np.zeros(N) # vetor de acoplamento da 2° variável
    
    for i1 in range(0,N,1):
        for i2 in range(1,Nvizinhos+1, 1) :

            acopl0[i1]=(contor0[N+i1-i2]+contor0[N+i1+i2])/(exp(i2*gamma))+acopl0[i1]
            acopl1[i1]=(contor1[N+i1-i2]+contor1[N+i1+i2])/(exp(i2*gamma))+acopl1[i1]        
        
#_____________________Entrando com as equações__________________________
    for i in range(0,N,1): 
        
        dudt[i] = -(b+1)*u[i]+a+(u[i]**2)*u[i+N]+(D0)*(acopl0[i]/norma-u[i])
        dudt[i+N] =  b*u[i]-(u[i]**2)*u[i+N]+(D1)*(acopl1[i]/norma-u[i+N])             
#_______________________________________________________________________

    return dudt
#______________________________________________________________________________________________________________________________


# In[175]:


#sol = odeint(fun, u0, t)
#sol = odeint(fun, u0, t,rtol=1.0e-13,atol=1.0e-13)
#sol = solve_ivp(fun, t, u0,rtol=1.0e-13,atol=1.0e-13)
sol_ivp = solve_ivp(fun=lambda tempo, y: fun(tempo, y), t_span=[min(t),max(t)], y0=u0, t_eval=t,rtol=1.0e-10,atol=1.0e-10,method='LSODA')


# In[176]:


sol=np.transpose(sol_ivp.y) #Para ficar igual a saída do odeint. As linhas serão os tempos e as colunas os sítios.


# In[177]:


#Plotagem de Gráfico#
#-----------------------------------------------------------------------------#
iterada = np.arange(0, iteracoes-1,1) #passo de iterações
sitios = np.arange(0, N,1)
iterada, sitios = np.meshgrid(iterada, sitios)
X0=sol[iterada,sitios]

#-----------------------------------------------------------------------------#
levels = np.arange(0, 16, 2)
norm = cm.colors.Normalize(vmax=abs(X0).max(), vmin=-abs(X0).max())
cmap = cm.gist_rainbow#cor
#-----------------------------------------------------------------------------#

fig = plt.figure(2)
ax = fig.add_subplot(111)
result = ax.contourf(sitios,iterada,X0,levels, norm=norm,cmap=cm.get_cmap(cmap, len(levels) - 1),extend='both')
result2 = ax.contour(result, levels=result.levels[::2], colors='c', norm=norm)

cbar = fig.colorbar(result)
cbar.ax.set_ylabel('X')

#result.cmap.set_under('white')
#result.cmap.set_over('orange')
#-----------------------------------------------------------------------------#

ax.set_xlabel(r'$\bf{Sitios}$', {'color': 'C0', 'fontsize': 20})
ax.set_ylabel(r'$\bf{Tempo}$', {'color': 'C0', 'fontsize': 20})
ax.set_yticklabels((0,10,20,30,40))
#plt.savefig('nomeDaFigura.jpg')
plt.show()


# In[97]:


# Exportando os resultados em arquivos .dat
#______________________________________________________________________________________________________
# 1° Maneira usando o savetxt do numpy do Python
#Escrevendo a matriz solução (para a 1º variável)
passoimp=1
np.savetxt('matriz_solucao2.dat', sol[::passoimp,0:N], fmt='%-.10e', delimiter='   ', newline='\n', header='', footer='', comments='# ', encoding=None)
np.savetxt('matriz_solucao2.dat', sol[:,0:N], fmt='%-.10e', delimiter='   ', newline='\n', header='', footer='', comments='# ', encoding=None)

#Escrevendo as séries temporais para sítios especificados
#np.savetxt('temporal.dat', np.column_stack([t,sol[:,0],sol[:,50]]), fmt='%-.10e', delimiter='   ', newline='\n', header='', footer='', comments='# ', encoding=None)


#_______________________________________________________________________________________________________
# 2° Maneira usando o write do Python
#serie_temporal = open('serie_temporal.dat', 'w')

# Imprimindo séries temporais
#for i in range(0,Npontos,1):
#    tempo=round(i*step,ndigits=3)
#    solu1=round(sol[i,1],ndigits=3)
#    solu2=round(sol[i,2],ndigits=3)
#    serie_temporal.write(str(tempo) +'    '+ str(solu1) +'    '+ str(solu2) +'\n')
#serie_temporal.close() 
#_______________________________________________________________________________________________________


# In[54]:


#Plotando gráficos do arquivo
#2D: Séries temporais
#t, sol0, sol1 = np.loadtxt('temporal.dat', dtype='float', delimiter='   ', unpack=True)
#plt.plot(t,sol0, 'b', label='XK')
#plt.plot(t,sol1, 'r', label='XN')
#plt.xlabel('tempo')
#plt.ylabel('y')
#plt.title('Séries temporais a partir do arquivo dat')
#plt.legend()
#plt.show()


#Exemplo de como plotar 3D a partir de dados externos
#Data3D= np.loadtxt("matriz_solucao1.dat", delimiter='   ')
#pts_imp=int(tfinal/(step*passoimp))
#print(pts_imp,Npontos)
#iterada = np.arange(0, Npontos-1,1000)
#iterada = np.arange(0, pts_imp-1,1)
#sitios = np.arange(0, N,1)
#iterada, sitios = np.meshgrid(iterada, sitios)
#X0=Data3D[iterada,sitios]
#fig = plt.figure()
#ax = Axes3D(fig)
#ax.plot_surface(iterada, sitios, X0, cmap='spring')
#ax.plot_wireframe(iterada, sitios, X0, rstride=10, cstride=10)
#ax.plot_surface(iterada, sitios, X0, rstride=1, cstride=1, cmap=cm.brg, linewidth=0, antialiased=True)
#ax.set_xlabel("t/(step*Qtd. Pontos)")
#ax.set_ylabel("k (sítios)")
#ax.set_zlabel("X0")
#ax.set_xlim([0.,Npontos])
#ax.set_xlim([0.,pts_imp])
#ax.set_ylim([0,N])
#ax.set_zlim([0.85,1.05])
#ax.view_init(elev=45, azim=30)
#plt.show()

