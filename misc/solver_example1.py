import numpy as np
import scipy
import matplotlib.pyplot as plt
import io#for Input Ouput
import math as m
from numpy.linalg import solve
#index = io.StringIO(u"index.txt")
index = np.loadtxt(open("index.txt")) #7x7 matrix
k = np.loadtxt(open("k.txt")) #7x7matrix
isopach = np.loadtxt(open('isopach.txt')) #7x7 matrix
dx = np.loadtxt(open('dx.txt')) #7x7 matrix
dy = np.loadtxt(open('dy.txt')) #7x7 matrix
beta_c = 1.127e-3
B = 1
rho = 62.4
miu = 1
class fluid_prop(object):
    #define class fluid properties
    #all the properties is constant
    def __init__(self,B = 1, rho = 62.4, miu = 1):
        self.B = B
        self.rho = rho
        self.miu = miu

class rock_prop(object):
    #define rock properties
    #include permeability, porosity and thick ness
    def __init__(self, perm, h, dx, dy):
        self.k = k #permeability
        self.h = isopach #thick ness
        self.dx = dx #delta x
        self.dy = dy # delta y

row = index.shape[0]
col = index.shape[1]
count = 0
order = np.zeros((row,col))
for i in range(row):
    for j in range(col):
        if index[i,j] > 0:
            count += 1
            order[i,j] = count

print("Number of cell:", (row-1)*(col-1))
print("Number of active cell", count)

#Define the North, West, East and South of each cell
N = np.zeros((row, col))
S = np.zeros((row, col))
W = np.zeros((row, col))
E = np.zeros((row, col))

#change the zero values of dx to 1
for i in range(7):
    for j in range(7):
        if dx[i,j] == 0:
            dx[i,j] = 1
#change zero values of dy to 1
for i in range(7):
    for j in range(7):
        if dy[i,j] ==0:
            dy[i,j] = 1
#calculation the transmissibility of N,W,N,S
Ax = np.multiply(dy,isopach)
Ay = np.multiply(dx,isopach)
#print(order)
#print(k)
#print(index)
# print(dx)
# print(dy)
#print(order)
for i in range(row):
	for j in range(col):
		if order[i,j] > 0:
			try:
				N[i,j] = 2*beta_c/(miu*B*(dy[i-1,j]/(Ay[i-1,j]*k[i-1,j]) + dy[i,j]/(Ay[i,j]*k[i,j])))

				S[i,j] = 2*beta_c/(miu*B*(dy[i+1,j]/(Ay[i+1,j]*k[i+1,j])) + (dy[i,j]/(Ay[i,j]*k[i,j])))

				W[i,j] = 2*beta_c/(miu*B*(dx[i,j-1]/(Ax[i,j-1]*k[i,j-1])) + (dx[i,j]/(Ax[i,j]*k[i,j])))

				E[i,j] = 2*beta_c/(miu*B*(dx[i,j+1]/(Ax[i,j+1]*k[i,j+1])) + (dx[i,j]/(Ax[i,j]*k[i,j])))
			except ZeroDivisionError:
				N[i,j] = 0
				S[i,j] = 0
				W[i,j] = 0
				E[i,j] = 0
#print(N)
Q = np.zeros(20)
C = -(N+W+E+S)
rw = 0.25
s = -0.5
for i in range(row):
	for j in range(col):
		if order[i,j] ==7:
			re = 0.14*(dx[i,j]**2 + dy[i,j]**2)**0.5
			omega1 = 2*beta_c*m.pi*k[i,j]*isopach[i,j]/(np.log(re/rw) + s)
			C[i,j] = C[i,j] - omega1
			Q[6] = -100*omega1
		if order[i,j] ==10:
			re = 0.14*(dx[i,j]**2 + dy[i,j]**2)**0.5
			omega2 = 2*beta_c*m.pi*k[i,j]*isopach[i,j]/(np.log(re/rw) + s)
			Q[9] = -200

mat = np.zeros((20,20))
for i in range(1,21):
	m = np.where(order ==i)[0][0]
	n = np.where(order ==i)[1][0]		
	mat[i-1,int(order[m,n-1] - 1)] = W[m,n]
	mat[i-1,int(order[m+1,n] - 1)] = S[m,n]
	mat[i-1,int(order[m-1,n] - 1)] = N[m,n]
	mat[i-1,int(order[m,n+1] - 1)] = E[m,n]
	mat[i-1,i-1] = C[m,n]


print(mat)
print(solve(mat,Q))

