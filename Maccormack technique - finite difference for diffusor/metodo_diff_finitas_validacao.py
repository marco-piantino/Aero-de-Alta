from math import *
import numpy as np
import matplotlib.pyplot as plt

# CONSTANTS
gamma = 1.4
C = 0.5         #Courant number

# INPUT
N = 30  #number of grid points
numTimeSteps = 1000 #Number of time steps
R = 287.05  #constant of air
N1 = N + 1  
L = 3

# NOZZLE AREA RATIO AND INITIAL CONDITIONS
dx = L/N
x = 0

A = np.zeros(N1)
rho_t = np.zeros(N1)
T_t = np.zeros(N1)
U_t = np.zeros(N1)
xmflow = np.zeros((numTimeSteps,N1)) 
xr = np.zeros(N1)

for i in np.arange(N1):
    A[i] = 1 + 2.2*(x - 1.5)**2
    rho_t[i] = 1 - 0.3146*x
    T_t[i] = 1 - 0.2314*x
    U_t[i] = (0.1 + 1.09*x)*sqrt(T_t[i])
    xmflow[0,i] = rho_t[i]*U_t[i]*A[i]
    xr[i] = x
    x = x + dx

# CALCULATION OF TIME STEP
delt_y = 1

for i in np.arange(1,N):
    delt_x = dx/(U_t[i] + sqrt(T_t[i]))
    delt_im = min(delt_x,delt_y)
    delt_y = delt_im
    delt_im = C*delt_im
    
time = delt_im

prho = np.zeros((numTimeSteps,N1))
pT = np.zeros((numTimeSteps,N1))
pU = np.zeros((numTimeSteps,N1))
p = np.zeros((numTimeSteps,N1))
x_mach = np.zeros((numTimeSteps,N1))

rho = np.zeros((numTimeSteps,N1))
T = np.zeros((numTimeSteps,N1))
U = np.zeros((numTimeSteps,N1))
P = np.zeros((numTimeSteps,N1))

# SOME ADDITIONAL VALUES TO BE INITIALIZED
prho[0,0] = rho_t[0]
pT[0,0] = T_t[0]
P[0,:] = 1

for i in np.arange(N1):
    rho[0,i] = rho_t[i]
    T[0,i] = T_t[i]
    U[0,i] = U_t[i]
    x_mach[0,i] = U_t[i]/sqrt(T_t[i])

for i in np.arange(1,numTimeSteps):
    
    drho_t = np.zeros((numTimeSteps,N1))
    dU_t = np.zeros((numTimeSteps,N1))
    dT_t = np.zeros((numTimeSteps,N1))
    
    Pdrho = np.zeros((numTimeSteps,N1))
    PdU = np.zeros((numTimeSteps,N1))
    PdT = np.zeros((numTimeSteps,N1))
    
    drho_av = np.zeros((numTimeSteps,N1))
    dU_av = np.zeros((numTimeSteps,N1))
    dT_av = np.zeros((numTimeSteps,N1))
    
    for j in np.arange(1,N):
        # PREDICTED VALUES FOR INTERNAL POINTS
        dxl_AP = (np.log(A[j+1]) - np.log(A[j]))/dx
        dx_UP = (U[i-1,j+1] - U[i-1,j])/dx
        dx_rhoP = (rho[i-1,j+1] - rho[i-1,j])/dx
        dx_TP = (T[i-1,j+1] - T[i-1,j])/dx
        
        drho_t[i,j] = - rho[i-1,j]*dx_UP -rho[i-1,j]*U[i-1,j]*dxl_AP - U[i-1,j]*dx_rhoP
        dU_t[i,j] = -U[i-1,j]*dx_UP - (1/gamma)*(dx_TP + (T[i-1,j]/rho[i-1,j])*dx_rhoP)
        dT_t[i,j] = -U[i-1,j]*dx_TP - (gamma - 1)*(T[i-1,j]*dx_UP + T[i-1,j]*U[i-1,j]*dxl_AP)
        
        #print(drho_t[j])
        prho[i,j] = rho[i-1,j] + delt_im*drho_t[i,j]
        pU[i,j] = U[i-1,j] + delt_im*dU_t[i,j]
        pT[i,j] = T[i-1,j] + delt_im*dT_t[i,j]
        
    # LINEAR EXTRAPOLATION FOR PRHO, PU AND PT
    prho[i,0] = prho[i-1,0]
    prho[i,N] = prho[i-1,N]
    pU[i,0] = 2*pU[i,1] - pU[i,2]
    pU[i,N] = 2*pU[i,N-1] - pU[i,N-2]
    pT[i,0] = pT[i-1,0]
    pT[i,N] = pT[i-1,N]
    
    for j in np.arange(1,N):
        # CORRECTED VALUES FOR INTERNAL POINTS
        dxl_AC = (np.log(A[j]) - np.log(A[j-1]))/dx
        dx_UC = (pU[i,j] - pU[i,j-1])/dx
        dx_rhoC = (prho[i,j] - prho[i,j-1])/dx
        dx_TC = (pT[i,j] - pT[i,j-1])/dx
        
        #print(-prho[i,j])
        Pdrho[i,j] = -prho[i,j]*dx_UC - prho[i,j]*pU[i,j]*dxl_AC - pU[i,j]*dx_rhoC
        #problema aqui
        
        PdU[i,j] = -pU[i,j]*dx_UC - (1/gamma)*(dx_TC + (pT[i,j]/prho[i,j])*dx_rhoC)
        PdT[i,j] = -pU[i,j]*dx_TC - (gamma - 1)*(pT[i,j]*dx_UC + pT[i,j]*pU[i,j]*dxl_AC)
        
        drho_av[i,j] = 0.5*(Pdrho[i,j] + drho_t[i,j])
        dU_av[i,j] = 0.5*(PdU[i,j] + dU_t[i,j])
        dT_av[i,j] = 0.5*(PdT[i,j] + dT_t[i,j])
        
        rho[i,j] = rho[i-1,j] + drho_av[i,j]*delt_im
        U[i,j] = U[i-1,j] + dU_av[i,j]*delt_im
        T[i,j] = T[i-1,j] + dT_av[i,j]*delt_im
        P[i,j] = rho[i,j]*T[i,j]
        x_mach[i,j] = U[i,j]/sqrt(T[i,j])
        
    # EXTRAPOLATION TO END POINTS
    rho[i,0] = rho[i-1,0]
    T[i,0] = T[i-1,0]
    P[i,0] = P[i-1,0]
    U[i,0] = 2*U[i,1] - U[i,2]
    x_mach[i,0] = U[i,0]/sqrt(T[i,0])
    
    rho[i,N] = 2*rho[i,N-1] - rho[i,N-2]
    U[i,N] = 2*U[i,N-1] - U[i,N-2]
    T[i,N] = 2*T[i,N-1] - T[i,N-2]
    P[i,N] = 2*P[i,N-1] - P[i,N-2]
    x_mach[i,N] = U[i,N]/sqrt(T[i,N])
    
    delt_y = 1
    for j in np.arange(1,N):
        delt_x = dx/(U[i,j] + sqrt(T[i,j]))
        delt_im = min((delt_x),(delt_y))
        delt_y = delt_im
        delt_im = C*delt_im
    
    for j in np.arange(N1):
        xmflow[i,j] = rho[i,j]*U[i,j]*A[j]
    
    time = time + delt_im

plt.plot(rho[:,15])
plt.suptitle("Density at i = 15")
plt.grid()
plt.show()

plt.plot(T[:,15])
plt.suptitle("Temperature at i = 15")
plt.grid()
plt.show()

plt.plot(P[:,15])
plt.suptitle("Pressure at i = 15")
plt.grid()
plt.show()

plt.plot(x_mach[:,15])
plt.suptitle("Mach at i = 15")
plt.grid()
plt.show()

plt.plot(x_mach[numTimeSteps - 1,:])
plt.suptitle("Mach through the nozzle")
plt.grid()
plt.show()

