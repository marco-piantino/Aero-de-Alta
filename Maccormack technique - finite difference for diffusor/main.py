from cmath import pi
from ctypes.wintypes import FLOAT
from math import sqrt
from tkinter import N
from matplotlib import pyplot as plt
import numpy as np
import matplotlib 
import DiffuserGeometry


#CONSTANTS
t_max = 3.1
gamma = 1.4
C = 0.5

#INPUTS
N: int = 150 #number of grid points
numTimeSteps: int = 500 #Number of time steps
R: float = 287.05 #constant of air
V_0: float = 2 #stagnation speed
rho_0: float = 1.201 #stagnation density
T_0: float = 300 #stagnation temperature
""" a_0 = sqrt(gamma*R*T_0) #stagnation air speed 
 """
#DIFFUSER GEOMETRY

#converging section
r_0: float = 0.077942132
L_conv: float = 0.069985

#throat
r_t: float = 0.070667
L_t: float = 0.4676

#diverging section
r_f: float = 0.170507
L_div: float = 1.910304

#total length
L = DiffuserGeometry.totalLength(L_conv, L_t, L_div)

##PRE-CALCULUS
#deltas
dx = L/(N-1)
dt = t_max/(numTimeSteps)

timeLength = numTimeSteps + 1
positionLength = N

rho_t=np.zeros([positionLength])
V_t=np.zeros([positionLength])
T_t=np.zeros([positionLength])


t=np.ones([positionLength])


rho_bar=np.ones([positionLength])
V_bar=np.ones([positionLength])
T_bar=np.ones([positionLength])


rho_tbar=np.ones([positionLength])
V_tbar=np.ones([positionLength])
T_tbar=np.ones([positionLength])


rho_avg=np.ones([positionLength])
T_avg=np.ones([positionLength])
V_avg=np.ones([positionLength])

# compliation of flowfield variable at each timestep
rho_r=np.zeros([positionLength,timeLength])
V_r=np.zeros([positionLength,timeLength])
T_r=np.zeros([positionLength,timeLength])

x = np.linspace(0, L, positionLength)

#radius function
radius = np.zeros([positionLength])
i = 0
for i in range(positionLength):
    radius[i] = DiffuserGeometry.radius(L_conv, L_t, L_div, r_0, r_t, r_f, x[i])

#area function
A_aux = pi*radius**2

#throat area
A_t = pi*r_t**2

#dimentionless area
A = A_aux/A_t


#INITIAL CONDITION
rho_e =1/4.170 #densidade na seção de testes
T_e = 1/2.19
V_e = 2
T = T_e + (1 - T_e)*x/L  
rho = rho_e + (1 - rho_e)*x/L
V = V_e + (0.01 - V_e)*x/L 


for time in range(0,numTimeSteps): #time step
    """ V[0] = 2*V[1]-V[2]  # Subsonic Inflow Boundary """
    V[0] = V_e
    V[N-1] = 2*V[N-2]-V[N-3] # Supersonic Outflow Boundary
    rho[N-1] = 2*rho[N-2]-rho[N-3] # Supersonic Outflow Boundary
    T[N-1] = 2*T[N-2]-T[N-3] # Supersonic Outflow Boundary

    V_r[:,time] = V 
    T_r[:,time] = T 
    rho_r[:,time] = rho 

    for i in range(1,N-1):
# Predictor Step
        rho_t[i]=-rho[i]/dx*(V[i+1]-V[i])-rho[i]*V[i]/dx*(np.log1p(A[i+1])-np.log1p(A[i]))-V[i]/dx*(rho[i+1]-rho[i])
        V_t[i]=-V[i]/dx*(V[i+1]-V[i])-1/gamma*((T[i+1]-T[i])/dx+T[i]/rho[i]/dx*(rho[i+1]-rho[i]))
        T_t[i]=-V[i]/dx*(T[i+1]-T[i])-(gamma-1)*T[i]*((V[i+1]-V[i])/dx+V[i]/dx*(np.log1p(A[i+1])-np.log1p(A[i])))

        t[i]=C*dx/(np.sqrt(T[i])+V[i]) # time step

        rho_bar[i]=rho[i]+rho_t[i]*min(t)
        V_bar[i]=V[i]+V_t[i]*min(t)
        T_bar[i]=T[i]+T_t[i]*min(t) 

        rho_tbar[i]=-rho_bar[i]/dx*(V_bar[i]-V_bar[i-1])-rho_bar[i]*V_bar[i]/dx*(np.log1p(A[i])-np.log1p(A[i-1]))-V_bar[i]/dx*(rho_bar[i]-rho_bar[i-1]) #t=t+dt
        V_tbar[i]=-V_bar[i]/dx*(V_bar[i]-V_bar[i-1])-1/gamma*((T_bar[i]-T[i-1])/dx+T_bar[i]/rho_bar[i]/dx*(rho_bar[i]-rho[i-1]))    #t=t+dt
        T_tbar[i]=-V_bar[i]/dx*(T_bar[i]-T_bar[i-1])-(gamma-1)*T_bar[i]*((V[i]-V[i-1])/dx+V_bar[i]/dx*(np.log1p(A[i])-np.log1p(A[i-1])))   #t=t+dt

        rho_avg[i]=0.5*(rho_t[i]+rho_tbar[i])
        T_avg[i]=0.5*(T_t[i]+T_tbar[i])
        V_avg[i]=0.5*(V_t[i]+V_tbar[i])

        rho[i]=rho[i]+rho_avg[i]*min(t)
        V[i]=V[i]+V_avg[i]*min(t)
        T[i]=T[i]+T_avg[i]*min(t)


V_r[:,timeLength-1] = V 
T_r[:,timeLength-1] = T 
rho_r[:,timeLength-1] = rho

print(rho)

print(V)

print(T)

plt.plot(x, rho)
plt.suptitle("Density through the diffuser")
plt.grid()
plt.show()

plt.plot(x, V)
plt.suptitle("Air speed through the diffuser")
plt.grid()
plt.show()

plt.plot(x, T)
plt.suptitle("Temperature through the diffuser")
plt.grid()
plt.show()

plt.plot(V_r[50,:])
plt.suptitle("Air speed on diffuser exit for each time step")
plt.grid()
plt.show()


