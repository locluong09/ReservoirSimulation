import math as m
import numpy as np
from numpy.linalg import solve


#input for calculation
k1,k2 = 200 , 100 #md
delta_x,delta_y = 600, 600 #ft
phi = 0.25 #porosity of each block
h = 100 #ft
rw = 0.25 #ft
s = 0

def calc_miu(p):
	#calculate gas density
	return 5*1e-6*p + 0.008

def calc_Z(p):
	#calculate compressibility factor
	return 5*1e-8*p**2 - 2*1e-4*p +1

def calc_B(p,z):
	#calculate formation volume factor B_g
	return 14.7/520/5.615*z*(120+460)/p

def calc_trans(k1, k2, delta_x, delta_y, miu, B, beta_c = 1.127e-3):
	#calculate transmissibility
    return 2*beta_c/(miu*B)/(delta_x/(h*delta_y*k1) + delta_x/(h*delta_y*k2))

def calc_omega(miu, B, k, h , re, rw = 0.25, beta_c = 1.127e-3):
	#calculate well productivity omega
	re = 0.198*600
	return 2*m.pi*beta_c*k*h/(miu*B*np.log(re/rw))

def calc_gamma(delta_x, delta_y,h,phi = 0.25, Tsc = 520, psc = 14.7, T = 580, t = 1): 
	#calculate gamma
	return delta_x*delta_y*h*phi*Tsc/(psc*T*t)

gamma1 = calc_gamma(delta_x, delta_y, h)
gamma2 = calc_gamma(delta_x, delta_y, h)

#initial pressure of each block at time step 0 day and compressibility Z
p0 = [4000, 4000]
Z0 = [calc_Z(4000), calc_Z(4000)]

##Main loop for solving pressure at time step 20
tol = 1e-6 # convergence tolerance of pressure solve
p1 = 3000 # initial pressure of block 1 at time 20 days.
p2 = 3000 # pressure of block 2 at time step 20 days.
cycle_count = 1
while True:
	ave_p =(p1 + p2)/2
	miu = calc_miu(ave_p)
	Z1 = calc_Z(p1)
	Z2 = calc_Z(p2)
	Z1_ = calc_Z(ave_p)
	Z2_ = calc_Z(ave_p)
	B1 = calc_B(ave_p,Z1_)
	B2 = calc_B(ave_p,Z2_)
	T1, T2 = calc_trans(k1, k2, delta_x, delta_y, miu, B1, beta_c = 1.127e-3),\
			 calc_trans(k1, k2, delta_x, delta_y, miu, B2, beta_c = 1.127e-3)
	print(T1,T2)
	C1 = -(T1 + gamma1/Z1)
	C2 = -(T2 + gamma2/Z2)
	Q1 = gamma1*(-p0[0]/Z0[0])+10e6
	Q2 = gamma2*(-p0[1]/Z0[1])
	LHS = np.array([[C1, T1], [T2, C2]])
	RHS = np.array([Q1,Q2]) 
	res = solve(LHS, RHS)
	error = abs(res[0] - p1) + abs(res[1] - p2)
	if error > tol:
		p1 = res[0]
		p2 = res[1]
		cycle_count += 1
	else:
		print(res)
		break

q2 = 10e6

#Incremental Material Balance Check
Vb1 = delta_x*delta_y*h
Vb2 = delta_x*delta_y*h
numerator =abs(Vb1/5.615*(phi/calc_B(res[0],calc_Z(res[0])) - phi/calc_B(4000,calc_Z(4000)))+\
 			Vb2/5.615*(phi/calc_B(res[1],calc_Z(res[1])) - phi/calc_B(4000,calc_Z(4000)))   )

denominator = q2*1

print("ICBM check",numerator/denominator)
#Residual check
p1,p2 = res[0], res[1]
print("Pressure distribution of block 1 is %d, and block 2 is %d "%(p1,p2))
omega1 = 2*m.pi*1.127e-3*200*100/(calc_miu(p1)*calc_B(p1,calc_Z(p1))*np.log(0.198*600/0.25))
psf1 = p1-10e6/omega1
print("Sand face pressure of well 1 is:",psf1)
