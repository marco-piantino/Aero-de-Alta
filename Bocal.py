from math import *
from sympy import *
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

# cálculo das propriedades isentrópicas
def isentrop_prop(M, gamma):
    #razão de temperaturas
    T0_T = 1 + ((gamma - 1)/2)*M**2
    
    #razão de pressões
    p0_p = (1 + ((gamma - 1)/2)*M**2)**(gamma/(gamma - 1))
    
    #razão de massa específica
    rho0_rho = (1 + ((gamma - 1)/2)*M**2)**(1/(gamma-1))
    
    #razão de áreas
    A_Aestrela = sqrt((1/M**2) * ((2/(gamma + 1))*(1 + ((gamma - 1)/2)*M**2))**((gamma + 1)/(gamma - 1)))
    
    #mi e ni
    if M >= 1:
        k = sqrt((gamma + 1)/(gamma - 1))
        mi = np.rad2deg(float(asin(1/M)))   #ângulo de Mach
        ni = np.rad2deg(float(k*atan((1/k)*sqrt(M**2 - 1)) - atan(sqrt(M**2 - 1)))) #função de Prandtl-Meyer
    else:
        mi = 'Bocal não supersônico'
        ni = 'Bocal não supersônico'
    return p0_p, rho0_rho, T0_T, A_Aestrela, mi, ni

# número de pontos da curva
n = 200

# propriedade do Ar
gamma = 1.4

# Mach desejado no Bocal (Me = mach esperado)
Me = 2

################ Altura da garganta do bocal (m) - analisar melhor
y_throat = 0.05

# Área da garganta do Bocal (m2)
A_throat = pi*y_throat**2

# Definindo ni_e
ni_e = isentrop_prop(Me, gamma)[5]

# Área da seção de testes (A_e = área esperada)
A_e = isentrop_prop(Me, gamma)[3]*A_throat

# Definição de theta 1 de acordo com o que foi proposto pelo NACA1651
theta_1 = float(((isentrop_prop(Me, gamma)[3])**(-2/9)) * (isentrop_prop(Me, gamma)[5])/2)

####################################################################
## Método de Foelsch

# Definição de ni_1 em função de ni_e e theta_1
ni_1 = isentrop_prop(Me, gamma)[5] - theta_1

# Definição do Raio (r_0) para M=1
r_0 = y_throat/np.deg2rad(theta_1)

# Definição de M1 tomando ni_1 pela Tabela A5 do Anderson e interpolando
# ni = 14.27 p/ M=1.58
# ni = 14.86 p/ M=1.60
def f(M, ni=ni_1, gamma=1.4):
    k = np.sqrt((gamma + 1)/(gamma - 1))
    return np.rad2deg(k*np.arctan((1/k)*np.sqrt(M**2 - 1)) - np.arctan(np.sqrt(M**2 - 1))) - ni_1 #função de Prandtl-Meyer

M1 = fsolve(f,1)[0]

# Cálculo de y_1 e x_1 com base nas propriedades isentrópicas para M1
r_1 = r_0*isentrop_prop(M1, gamma)[3]
y_1 = r_1*np.sin(np.deg2rad(theta_1))
x_1 = (3/2)*(y_1 - y_throat)/np.tan(np.deg2rad(theta_1))

# Intervalo de Mach a ser percorrido
n_mach = np.linspace(M1, Me, num=n)

# Estimar as propriedades isentrópicas e calcular os valores de r, l, x_2, y_2, x e y para  cada Mach varrido
x = []
y = []
for M in n_mach:
    r = r_0*isentrop_prop(M, gamma)[3]
    l = M*r*np.deg2rad(isentrop_prop(M, gamma)[5] - ni_1)
    x_2 = -r_1*np.cos(np.deg2rad(theta_1)) + r*np.cos(np.deg2rad(ni_e - isentrop_prop(M, gamma)[5])) + x_1
    y_2 = r*np.sin(np.deg2rad(ni_e  - isentrop_prop(M, gamma)[5]))
    x = x + [x_2 + l*np.cos(np.deg2rad(ni_e - isentrop_prop(M, gamma)[5] + isentrop_prop(M, gamma)[4]))]
    y = y + [y_2 + l*np.sin(np.deg2rad(ni_e - isentrop_prop(M, gamma)[5] + isentrop_prop(M, gamma)[4]))]

plt.plot(np.array(x), np.array(y))
plt.plot(np.array(x), -np.array(y))

x_cubica = np.linspace(-0.20, float(x[0]), num=n)
y_cubica = []
for x in x_cubica:
    y = y_throat + (np.tan(np.deg2rad(theta_1))/x_1)*(x**2)*(1 - x/(3*x_1))
    y_cubica = y_cubica + [y]

plt.plot(np.array(x_cubica), np.array(y_cubica))
plt.plot(np.array(x_cubica), -np.array(y_cubica))