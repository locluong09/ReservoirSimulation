import math as m
import numpy as np
from numpy.linalg import solve

# Assume input data for reservoir
kx,ky = 100,100
delta_x,delta_y = 200,200
c = 1e-5 #psi^-1
phi = 0.2
h = 100 #ft
miu_sc = 0.9 #cp
p_sc = 14.7 #psi
rho_sc = 62.4 #lbm/ft^3
r_w = 0.25 #ft
s = 0

def clc_miu(p1, p2, miu_sc = 0.9):
    return miu_sc + 0.00134*np.log((p1+p2)/2)

def clc_B(p1, p2, c= 1e-5, p_sc = 14.7):
    return np.power(1+c*(p1/2+p2/2-p_sc), -1)

def clc_trans(k_x, k_y, h, delta_x, delta_y, miu, B, beta_c = 1.127e-3):
    Ax = delta_y*h
    return 2*beta_c/(miu*B)/(delta_x/(Ax*k_x) + (delta_x/(Ax*k_y)))

def clc_gamma(delta_x, delta_y, h, phi, delta_t=30, B_0=1, c = 1e-5, alpha_c = 5.615):
    return delta_x*delta_y*h*phi*c/(alpha_c*B_0*delta_t)

def clc_rho(p1,p2, c =1e-5, rho_sc = 62.4):
    return rho_sc*(1+c*(p1/2+p2/2-14.7))

def clc_trans_(k_x, k_x1, h, delta_x, delta_y, rho, miu, B, beta_c = 1.127e-3):
    Ax = delta_y*h
    return 2*beta_c*rho/(144*miu*B)/(delta_x/(Ax*k_x) + (delta_x/(Ax*k_x1)))

def clc_PI(k_x, k_y, h, delta_x, delta_y, miu, B, r_w = 0.25, s = 0, beta_c = 1.127e-3):
    r_e = 0.14*m.sqrt(delta_x**2 + delta_y**2)
    return 2*beta_c*m.pi*m.sqrt(k_x*k_y)*h/(miu*B*(np.log(r_e/r_w) + s))

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~2.1~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("FALSE / TRUE")
print("Because the system of equations can be solved by iteration method, so if we\n\
just solve the pressure distribution by iteration method, the statement should be FALSE.")
print("If we updade the SIP coefficients of single-phase slightly compressible \n\
fluid reservoir within the SIP iterations similar to a gas reservoir, this \n\
pressure distribution will be totally wrong because at time step n+1, all values of transmissibility\n\
miu, B, rho are calculated by values of time step n. So this statement should be TRUE.")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~2.2~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("FALSE")
print("Running the cases below: Reservoir with 50md permeability and\n\
			G=[3100,2990,3100]")

kx1,ky1 = 50,50

gamma = clc_gamma(delta_x, delta_y, h, phi, delta_t=10, B_0=1, c = 1e-5, alpha_c = 5.615)
#print(gamma)
p_init = [3000,3000,3000]
miu1 = clc_miu(p_init[0], p_init[1], miu_sc = 0.9)
miu2 = clc_miu(p_init[1], p_init[2], miu_sc = 0.9)

B1 = clc_B(p_init[0], p_init[1], c= 1e-5, p_sc = 14.7)
B2 = clc_B(p_init[2], p_init[1], c= 1e-5, p_sc = 14.7)

rho1 = clc_rho(p_init[0], p_init[1])
rho2 = clc_rho(p_init[1], p_init[2])

T1 = clc_trans(kx1, ky1, h, delta_x, delta_y, miu1, B1, beta_c = 1.127e-3)
T2 = clc_trans(kx1, ky1, h, delta_x, delta_y, miu2, B2, beta_c = 1.127e-3)

T1_ = clc_trans_(50, 50, h, delta_x, delta_y, rho1, miu1, B1, beta_c = 1.127e-3)
T2_ = clc_trans_(50, 50, h, delta_x, delta_y, rho2, miu1, B1, beta_c = 1.127e-3)

C1 = -(T1 + gamma)
C2 = -(T1+T2+gamma)
C3 = -(T2 + gamma)

G = [3100,2990,3100]
Q1 = -100 +T1_*(G[1]-G[0])
Q2 = T1_*(G[0]-G[1]) + T2_*(G[2] - G[1])
Q3 = T2_*(G[1] - G[2])

#print(C1,C2,C3,T1)
LHS = np.array([[C1,T1,0],[T1,C2,T2],[0,T2,C3]])
RHS = np.array([-gamma*p_init[0]+Q1, -gamma*p_init[1] +Q2, -gamma*p_init[2]+Q3])
x2 = solve(LHS,RHS)
print("Pressure distribution of reservoir",x2)



print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~2.3~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("TRUE")
print("Running 2 cases below: Reservoir 1 with 100md permeability and\n\
			Reservoir 2 with 50md permeability")
kx,ky =100,100
gamma = clc_gamma(delta_x, delta_y, h, phi, delta_t=500, B_0=1, c = 1e-5, alpha_c = 5.615)
omega = clc_PI(kx, ky, h, delta_x, delta_y, miu1, B1, r_w = 0.25, s = 0, beta_c = 1.127e-3)
#p_init = [3000,3000,3000]
#print(omega)
miu = clc_miu(3000, 3000, miu_sc = 0.9)
B = clc_B(3000, 3000, c= 1e-5, p_sc = 14.7)
T1 = clc_trans(kx, ky, h, delta_x, delta_y, miu, B, beta_c = 1.127e-3)
C1 = -(T1 + gamma )
C2 = -(2*T1+gamma)
C3 = -(T1 + gamma )
#print(C1,C2,C3,T1)

LHS = np.array([[C1,T1,0],[T1,C2,T1],[0,T1,C3]])
RHS = np.array([-100-gamma*3000,-gamma*3000,100-gamma*3000])
x = solve(LHS,RHS)
print("Pressure distribution of R1 with greater perm",x)



kx1,ky1 = 50,50
gamma1 = clc_gamma(delta_x, delta_y, h, phi, delta_t=500, B_0=1, c = 1e-5, alpha_c = 5.615)
omega1 = clc_PI(kx1, ky1, h, delta_x, delta_y, miu1, B1, r_w = 0.25, s = 0, beta_c = 1.127e-3)
#p_init = [3000,3000,3000]
#print(omega1)
miu = clc_miu(3000, 3000, miu_sc = 0.9)
B = clc_B(3000, 3000, c= 1e-5, p_sc = 14.7)
T1 = clc_trans(kx1, ky1, h, delta_x, delta_y, miu, B, beta_c = 1.127e-3)
C1 = -(T1 + gamma1)

C2 = -(2*T1+gamma1)
C3 = -(T1 + gamma1)
#print(C1,C2,C3,T1)
LHS = np.array([[C1,T1,0],[T1,C2,T1],[0,T1,C3]])
RHS = np.array([-100-gamma1*3000,-gamma1*3000,100-gamma1*3000])
x1 = solve(LHS,RHS)
print("Pressure distribution of R2 with smaller perm",x1)




