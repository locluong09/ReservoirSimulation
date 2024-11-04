import numpy as np
import math as m

def clc_miu(p1, p2, miu_sc = 0.9):
    return miu_sc + 0.00134*np.log((p1+p2)/2)

def clc_B(p1, p2, c= 1e-5, p_sc = 14.7):
    return np.power(1+c*(p1/2+p2/2-p_sc), -1)

def clc_rho(p1,p2, c =1e-5, rho_sc = 62.4):
    return rho_sc*(1+c*(p1/2+p2/2-14.7))

def clc_trans_x(k_x, k_x1, h, delta_x, delta_y, miu, B, beta_c = 1.127e-3):
    Ax = delta_y*h
    return 2*beta_c/(miu*B)/(delta_x/(Ax*k_x) + (delta_x/(Ax*k_x1)))

def clc_trans_y(k_y, k_y1, h, delta_x, delta_y, miu, B, beta_c = 1.127e-3):
    Ay = delta_x*h
    return 2*beta_c/(miu*B)/(delta_y/(Ay*k_y) + (delta_y/(Ay*k_y1)))

def clc_PI(k_x, k_y, h, delta_x, delta_y, miu, B, r_w = 0.25, s = 0, beta_c = 1.127e-3):
    r_e = 0.14*m.sqrt(delta_x**2 + delta_y**2)
    return 2*beta_c*m.pi*m.sqrt(k_x*k_y)*h/(miu*B*(np.log(r_e/r_w) + s))

def clc_trans_x_(k_x, k_x1, h, delta_x, delta_y, rho, miu, B, beta_c = 1.127e-3):
    Ax = delta_y*h
    return 2*beta_c*rho/(144*miu*B)/(delta_x/(Ax*k_x) + (delta_x/(Ax*k_x1)))

def clc_trans_y_(k_y, k_y1, h, delta_x, delta_y, rho, miu, B, beta_c = 1.127e-3):
    Ay = delta_x*h
    return 2*beta_c*rho/(144*miu*B)/(delta_y/(Ay*k_y) + (delta_y/(Ay*k_y1)))

def clc_gamma(delta_x, delta_y, h, phi, delta_t, B_0, c = 1e-5, alpha_c = 5.615):
    return delta_x*delta_y*h*phi*c/(alpha_c*B_0*delta_t)


'''EXERCISE 1 '''
#Calculate the SIP coefficients C and Q for Block 6 during the time period 0 ≤ t ≤ 30d
#AT TIME STEP 30 DAYS
print("************************")
print("Home work number 1")
print("************************")

k_x, k_y = 50, 50 #md
delta_x, delta_y = 400, 400 #ft
c = 1e-5 #psi^-1
phi = 0.2
h = 100 #ft
miu_sc = 0.9 #cp
p_sc = 14.7 #psi
rho_sc = 62.4 #lbm/ft^3
r_w = 0.25 #ft
s = 0

miu = clc_miu(3000, 3050)
B = clc_B(3000,3050)
rho = clc_rho(3000,3050)
delta_t = 30
B_0 = clc_B(3000,3050)
miu_0 = clc_miu(3000,3050)

E6 = clc_trans_x(k_x, k_x, h, delta_x, delta_y, miu, B)
N6 = clc_trans_y(k_y, k_y, h, delta_x, delta_y, miu, B)
E6_ = clc_trans_x_(k_x, k_x, h, delta_x, delta_y, rho, miu, B, beta_c = 1.127e-3)
N6_ = clc_trans_y_(k_y, k_y, h, delta_x, delta_y, rho, miu, B, beta_c = 1.127e-3)

omega6 = clc_PI(k_x, k_y, h, delta_x, delta_y, miu_0, B_0, r_w = 0.25, s = 0)
gamma6 = clc_gamma(delta_x, delta_y, h, phi, delta_t, B_0=1, c = 1e-5, alpha_c = 5.615)

C6 = -(N6+E6+gamma6 )
Q6 = -500- gamma6*3000+ E6_*(3050-3000) + N6_*(3050-3000)
print("Question a")
print("SIP coefficient of block 6 during time period from 0 to 30 days C6", C6)
print("SIP coefficient of block 6 during time period from 0 to 30 days Q6", Q6)

#AT TIME STEP 50 DAYS
#Recalculate the SIP coefficients C and Q for Block 6 during the time period 30 ≤t ≤ 50d
miu_50 = clc_miu(3800, 4000)
B_50 = clc_B(3800,4000)
rho_50 = clc_rho(3800,4000)
delta_t_50 = 20
#B_0 = clc_B(3800,4000)
B_0_50 = clc_B(3800,3800)
miu_0_50 = clc_miu(3800,3800)


E6 = clc_trans_x(k_x, k_x, h, delta_x, delta_y, miu_50, B_50)
N6 = clc_trans_y(k_y, k_y, h, delta_x, delta_y, miu_50, B_50)
E6_ = clc_trans_x_(k_x, k_x, h, delta_x, delta_y, rho_50, miu_50, B_50, beta_c = 1.127e-3)
N6_ = clc_trans_y_(k_y, k_y, h, delta_x, delta_y, rho_50, miu_50, B_50, beta_c = 1.127e-3)

omega6 = clc_PI(k_x, k_y, h, delta_x, delta_y, miu_0_50, B_0_50, r_w = 0.25, s = 0)
gamma6 = clc_gamma(delta_x, delta_y, h, phi, delta_t_50, B_0_50, c = 1e-5, alpha_c = 5.615)

C6 = -(N6+E6+omega6+gamma6)
Q6 = -omega6*1500 - gamma6*3800+ E6_*(3050-3000) + N6_*(3050-3000)

print("Question b")
print("SIP coefficient of block 6 during time period from 30 to 50 days C6", C6)
print("SIP coefficient of block 6 during time period from 30 to 50 days Q6", Q6)

'''Calculating the coefficients matrix from t=0 to t=30 day'''
#Construct the coefficient matrix during the time period 0 ≤ t ≤ 30d

A = np.ones((3,4))
A[0,0], A[0,2], A[0,3], A[2,0], A[2,3] =0,0,0,0,0
depth = np.zeros((A.shape[0], A.shape[1]))
pressure = np.zeros((A.shape[0], A.shape[1]))
depth[0,1], depth[1,0], depth[1,1], depth[1,2], depth[1,3], depth[2,1], depth[2,2] = 3000, 3000, 3050, 3100, 3150, 3000, 3050
pressure[0,1], pressure[1,0], pressure[1,1], pressure[1,2], pressure[1,3], pressure[2,1], pressure[2,2] = 3000, 3000, 3050, 3100, 3150, 3000, 3050

def boundaries(X):
    mat = np.zeros((X.shape[0]+2, X.shape[1]+2))
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            mat[i+1, j+1] = X[i,j]
    return mat

X, Dep, Pre = boundaries(A), boundaries(depth), boundaries(pressure)
row, col = X.shape[0], X.shape[1]
N,W,E,S = np.zeros((row,col)), np.zeros((row,col)), np.zeros((row,col)), np.zeros((row,col))
N_,W_,E_,S_ = np.zeros((row,col)), np.zeros((row,col)), np.zeros((row,col)), np.zeros((row,col))
C,Q = np.zeros((row,col)), np.zeros((row,col))
gamma = clc_gamma(400,400,100,0.2,30,1,1e-5,5.615)

for i in range(row):
    for j in range(col):
        if X[i,j] > 0:
            if X[i-1, j] > 0:
                miu = clc_miu(Pre[i-1, j], Pre[i,j])
                B = clc_B(Pre[i-1, j], Pre[i,j])
                rho = clc_rho(Pre[i-1, j], Pre[i,j])
                N[i,j] = clc_trans_x(k_x, k_x, h, delta_x, delta_y, miu, B)
                N_[i,j] = clc_trans_x_(k_x, k_x, h, delta_x, delta_y, rho, miu, B, beta_c = 1.127e-3)
            if X[i+1, j] > 0:
                miu = clc_miu(Pre[i+1, j], Pre[i,j])
                B = clc_B(Pre[i+1, j], Pre[i,j])
                rho = clc_rho(Pre[i+1, j], Pre[i,j])
                S[i,j] = clc_trans_x(k_x, k_x, h, delta_x, delta_y, miu, B)
                S_[i,j] = clc_trans_x_(k_x, k_x, h, delta_x, delta_y, rho, miu, B, beta_c = 1.127e-3)
            if X[i,j-1] > 0:
                miu = clc_miu(Pre[i, j-1], Pre[i,j])
                B = clc_B(Pre[i, j-1], Pre[i,j])
                rho = clc_rho(Pre[i, j-1], Pre[i,j])
                W[i,j] = clc_trans_x(k_x, k_x, h, delta_x, delta_y, miu, B)
                W_[i,j] = clc_trans_x_(k_x, k_x, h, delta_x, delta_y, rho, miu, B, beta_c = 1.127e-3)
            if X[i,j+1] > 0:
                miu = clc_miu(Pre[i, j+1], Pre[i,j])
                B = clc_B(Pre[i, j+1], Pre[i,j])
                rho = clc_rho(Pre[i, j+1], Pre[i,j])
                E[i,j] = clc_trans_x(k_x, k_x, h, delta_x, delta_y, miu, B)
                E_[i,j] = clc_trans_x_(k_x, k_x, h, delta_x, delta_y, rho, miu, B, beta_c = 1.127e-3)
            C[i,j] = -(N[i,j]+W[i,j]+E[i,j]+S[i,j] + gamma)
            Q[i,j] = N_[i,j]*(Dep[i-1,j] - Dep[i,j]) + S_[i,j]*(Dep[i+1,j] - Dep[i,j])+\
            W_[i,j]*(Dep[i,j-1] - Dep[i,j]) + E_[i,j]*(Dep[i,j+1] - Dep[i,j]) - gamma*Pre[i,j]      

order = np.zeros((row,col))
count = 0
for i in range(row):
    for j in range(col):
        if X[i,j] > 0:
            count += 1
            order[i,j] = count
            
order = order.astype(int)
LHS = np.zeros((7,7))
for i in range(row):
    for j in range(col):
        if order[i,j] > 0:
            LHS[order[i,j]-1,order[i,j] - 1] = C[i,j]
            if order[i-1,j] > 0:
                LHS[order[i,j] -1, order[i-1,j]-1] = N[i,j]
            if order[i+1,j] > 0:
                LHS[order[i,j] -1, order[i+1,j] -1] = S[i,j]
            if order[i,j-1] > 0:
                LHS[order[i,j] -1, order[i,j-1] -1] = W[i,j]
            if order[i,j+1] > 0:
                LHS[order[i,j] -1, order[i,j+1] -1] = E[i,j]

RHS = []
for i in range(row):
    for j in range(col):
        if order[i,j] > 0:
            RHS.append(Q[i,j])
print("Question c")
print('Coefficient matrix during the time period 0 ≤ t ≤ 30d, LHS is \n',LHS)
print('Coefficient matrix during the time period 0 ≤ t ≤ 30d, RHS is \n',RHS)

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
'''EXERCISE NUMBER 2'''
# Input of homework number 2
print("************************")
print("Homework number 2")
print("************************")

k1, k2 = 200, 100 #md
delta_x, delta_y = 600, 600 #ft
c = 1e-5 #psi^-1
phi = 0.25
h = 100 #ft
miu_sc = 0.9 #cp
p_sc = 14.7 #psi
rho_sc = 62.4 #lbm/ft^3
r_w = 0.25 #ft
s = 0
B = clc_B(3000, 3000)
p_i = 3000 #psi


miu = clc_miu(3000, 3000)
rho = clc_rho(3000, 3000)
delta_t = 30

E1 = clc_trans_x(k1, k2, h, delta_x, delta_y, miu, B)
W2 = clc_trans_y(k1, k2, h, delta_x, delta_y, miu, B)

omega1 = clc_PI(k1, k1, h, delta_x, delta_y, miu, B, r_w = 0.25, s = 0)
gamma1 = clc_gamma(delta_x, delta_y, h, phi, delta_t, B_0 = 1, c = 1e-5, alpha_c = 5.615)
gamma2 = clc_gamma(delta_x, delta_y, h, phi, delta_t, B_0 = 1, c = 1e-5, alpha_c = 5.615)

C1 = - (E1 + gamma1)
C2 = - (W2 + gamma2)

from numpy.linalg import solve
LHS1 = np.array([[C1,E1], [W2,C2]])
RHS1 = np.array([1000-gamma1*3000, -gamma2*3000])
p1_30, p2_30 = solve(LHS1,RHS1)
print("Pressure of block 1  is %d (psi), and block 2 is %d (psi) at time step 30 days" %(p1_30,p2_30))
p1sf_30 = -1000/omega1 + p1_30
print("Sand face pressure at block 1 is %d (psi)" %p1sf_30)

V = delta_x*delta_y*h
B_30 = clc_B(p1_30, p2_30)
#IMBC check
imbc_30 = abs(V*phi/5.615*(1/B_30 - 1/B)/(30*1000) +V*phi/5.615*(1/B_30-1/B)/(30*1000))
print("Incremental Material Balance Check at time step 30 days is", imbc_30)

print("************************")
miu = clc_miu(p1_30, p2_30)
rho = clc_rho(p1_30, p2_30)
delta_t = 30
B_30 = clc_B(p1_30, p2_30)
B_ = clc_B(p1_30, p1_30)


E1 = clc_trans_x(k1, k2, h, delta_x, delta_y, miu, B_30)
W2 = clc_trans_y(k1, k2, h, delta_x, delta_y, miu, B_30)

omega1 = clc_PI(k1, k1, h, delta_x, delta_y, miu, B_, r_w = 0.25, s = 0)
gamma1 = clc_gamma(delta_x, delta_y, h, phi, delta_t, B_0 = 1, c = 1e-5, alpha_c = 5.615)
gamma2 = clc_gamma(delta_x, delta_y, h, phi, delta_t, B_0 = 1, c = 1e-5, alpha_c = 5.615)

C1 = - (E1 + gamma1)
C2 = - (W2 + gamma2)

LHS2 = np.array([[C1,E1], [W2,C2]])
RHS2 = np.array([1000-gamma1*p1_30, -gamma2*p2_30])
p1_60, p2_60 = solve(LHS2,RHS2)
B_60 = clc_B(p1_60,p2_60)
print("Pressure of block 1  is %d (psi), and block 2 is %d (psi) at time step 60 days" %(p1_60,p2_60))
p1sf_60 = -1000/omega1+p1_60
print("Sand face pressure at block 1 is %d (psi)" %p1sf_60)
imbc_60 = abs(V*phi/5.615*(1/B_60 - 1/B_30)/(30*1000) +V*phi/5.615*(1/B_60-1/B_30)/(30*1000))
cmbc_60 = abs(V*phi/5.615*(1/B_60 - 1/B)/(60*1000) +V*phi/5.615*(1/B_60-1/B)/(60*1000))
print("Incremental Material Balance Check at time step 60 days is", imbc_60)
print("Cumulative Material Balance Check at time step 60 days is", cmbc_60)


