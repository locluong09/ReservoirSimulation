import math as m
import numpy as np

# input parameter
kx,ky,kz = 25,100,20
delta_x, delta_y = 500, 500
h = 100
miu = 1
B  = 1
rw = 0.25

def clc_trans_x(kx1, kx2, h, miu, B, delta_x, delta_y, beta_c = 1.127e-3):
	return 2*beta_c/(miu*B)/(delta_x/(h*delta_y*kx1) + delta_x/(h*delta_y*kx2))

def clc_trans_y(ky1, ky2, h, miu, B, delta_x, delta_y, beta_c = 1.127e-3):
	return 2*beta_c/(miu*B)/(delta_y/(h*delta_x*ky1) + delta_y/(h*delta_x*ky2))


T12 = clc_trans_y(ky,ky,h,miu,B,delta_x,delta_y)
T23 = clc_trans_x(kx,kx,h,miu,B,delta_x,delta_y)

re1 = 0.28*m.sqrt(m.sqrt(ky/kx)*delta_x**2 + m.sqrt(kx/ky)*delta_y**2)/\
 		(pow(ky/kx,0.25) + pow(kx/ky,0.25))
re3 = 0.28*m.sqrt(m.sqrt(kz/kx)*delta_x**2 + m.sqrt(kx/kz)*h**2)/\
 		(pow(kz/kx,0.25) + pow(kx/kz,0.25))

#print(re1, re3)

omega1 = 2*m.pi*m.sqrt(kx*ky)*h*1.127e-3/(miu*B*np.log(re1/rw))
omega3 = 2*m.pi*m.sqrt(kx*kz)*delta_y*1.127e-3/(miu*B*np.log(re3/rw))
#print(2*m.pi*50*100*1.127e-3/np.log(104/rw))
#print(omega1,omega3)

C1 = -(T12 + omega1)
C2 = -(T12+T23)
C3 = - (T23 + omega3)
Q1 = -omega1*3150
Q3 = -omega3*14.7

from numpy.linalg import solve
LHS = np.array([[C1, T12, 0], [T12, C2, T23], [0, T23, C3]])
RHS = np.array([Q1, 0, Q3])
#print(LHS,RHS)
re = solve(LHS, RHS)
print("Pressure distribution is",re)

q1 = -omega1*(re[0] - 3150)
q3 = -omega3*(re[2] - 14.7)
print("Flowrate of well 1 and well 3 is", q1, q3)
print('Material balance check successfully',q1+q3)







