from __future__ import division, print_function # to support python2 and python3
import numpy as np
from numpy import meshgrid
from numpy.linalg import solve
import scipy
import matplotlib.pyplot as plt
import math as m
import seaborn as sns
import time


def boundaries(X):
	# Make zero values outside a matrix
	#This function return a matrix with outer boundaries is zero values
	row = X.shape[0]
	col = X.shape[1]
	X_new = np.zeros((row+2, col+2))
	for i in range(row):
		for j in range(col):
			X_new[i+1, j+1] = X[i,j]
	return X_new


def biharmonic_spline_interpolate(x,y,f,X,Y):
	#function of interpolation data
	#same with matlab biharmonic spline interpolate v4
	Z = np.zeros((X.shape[0], X.shape[1]))
	length = len(f)
	G = np.zeros((length, length), dtype = 'float32')
	for i in range(length):
		for j in range(length):
			if i != j:
				Ma = m.sqrt((x[i] - x[j])**2 + (y[i] - y[j])**2)
				if Ma >= 1e-7:
					G[i,j] = Ma**2*(np.log(Ma) -1)

	a = np.asarray(solve(G,f)).reshape(-1,1)
	g = np.zeros((a.shape[0], a.shape[1]))
	for i in range(Z.shape[0]):
		for j in range(Z.shape[1]):
			for k in range(length):
				Ma1 = m.sqrt((X[i,j] - x[k])**2 + (Y[i,j] -y[k])**2)
				if Ma1 >= 1e-7:
					g[k] = Ma1**2*(np.log(Ma1) -1)
				else:
					g[k] = Ma1**2 * (-100)
			Z[i,j] = sum(np.multiply(g,a))

	return Z


def load_data():
	thickness = np.loadtxt(open('new_thickness.txt'))
	top = np.loadtxt(open('Top_new.txt'))
	perm = np.loadtxt(open('k.txt'))
	dx = np.loadtxt(open('x.txt'))
	dy = np.loadtxt(open('y.txt'))
	porosity = np.loadtxt(open('new_por.txt'))

	return thickness, top, perm, dx, dy, porosity


def data_cleaning(thickness, top, perm, dx, dy, porosity):
	'''This fucntuon clean the data and change some properties to
	for further calculation'''
	row, col = perm.shape[0], perm.shape[1]

	#change the perm matrix zones corresponding with their multiplier
	for i in range(row):
		for j in range(col):
			if (i < row//2) and (j < col//2):
				perm[i,j] *= 2.0
			if (i < row//2) and (j >= col//2):
				perm[i,j] *= 0.3
			if (i >= row//2) and (j < col//2):
				perm[i,j] *= 1.0
			if (i >= row//2) and (j >= col//2):
				perm[i,j] *= 3.0


	kx = perm
	ky = perm/2

	for i in range(kx.shape[0]):
		for j in range(kx.shape[1]):
			if kx[i,j] < 0:
				kx[i,j] = 0
			if ky[i,j] < 0:
				ky[i,j] = 0

	#make boundaries
	kx = boundaries(kx)
	ky = boundaries(ky)
	h = thickness
	top = boundaries(top)
	poro = porosity
	Row, Col = kx.shape[0], kx.shape[1]

	#create the dx and dy matrix.
	deltax = np.zeros((Row-2,Col -2))
	deltay = np.zeros((Row -2, Col -2))

	x = np.zeros((1,len(dx) - 1))
	for i in range(Col -2):
		x[0,i] = dx[i+1] - dx[i]

	y = np.zeros((len(dy) - 1,))
	for i in range(Row -2):
		y[i,] = dy[i+1] - dy[i]

	for i in range(deltax.shape[0]):
		deltax[i,:] = x
	for j in range(deltay.shape[1]):
		deltay[:,j] = y
	dx = boundaries(deltax)
	dy = boundaries(deltay)
	return kx, ky, h, top, dx, dy, poro





class props_rock():
	'''This class captures the rock properties of
	reservoir including kx,ky, perm and poro'''
	def __init__(self, kx, ky, poro, top):
		self.kx = kx # permeability in x direction matrix
		self.ky = ky # permeability in y direction matrix
		self.poro = poro # porosity matrix
		self.top = top




class props_fluid():

	'''This class present the fuild properties inside the reservioi
	including formation volume factor of oil, viscosity of oil and oil compresisbility'''
	def __init__(self, cr_o, cr):
		self.cr_o = cr_o # compressibility of oil
		self.cr = cr #compressibility of rock




class props_reservoir():
	def __init__(self, dx, dy, h, p_init):
		self.dx = dx
		self.dy = dy
		self.h = h
		self.p_init = p_init

	def active_block(self, kx):
		Row, Col = kx.shape[0], kx.shape[1]
		count = 0
		order = np.zeros((Row,Col))
		for i in range(Row):
			for j in range(Col):
				if kx[i,j] > 0:
					count += 1
					order[i,j] = count
		order = order.astype(int)
		return order, count


class props_well():
	def __init__(self, well_coor, well_spec, well_ppt, well_state):
		self.well_coor = well_coor
		self.well_spec = well_spec
		self.well_ppt = well_ppt
		self.well_state = well_state
	def PI(self, kx, ky, dx, dy, p, h):
		i = self.well_coor[0]
		j = self.well_coor[1]
		rw = self.well_ppt[0]
		s = self.well_ppt[1]
		re = 0.28*m.sqrt((ky[i,j]/kx[i,j])**0.5*dx[i,j]**2 + (kx[i,j]/ky[i,j])**0.5*dy[i,j]**2)/\
		((ky[i,j]/kx[i,j])**0.25+(kx[i,j]/ky[i,j])**0.25)
		miu_well = clc_miu(p[i,j], p[i,j])
		B_well = clc_B(p[i,j], p[i,j])
		omega = 2*m.pi*1.127e-3*m.sqrt(kx[i,j]*ky[i,j])*h[i,j]/(miu_well*B_well*(np.log(re/rw) +s ))

		return omega


class props_time():
	def __init__(self, time_step = 0, time_interval = 0):
		self.time_step = time_step
		self.time_interval = time_interval


def clc_miu(p1, p2, miu_sc = 0.6):
    return miu_sc + 0.00134*np.log((p1+p2)/2)

def clc_B(p1, p2, c= 3e-5, p_sc = 14.7):
    return np.power(1+c*(p1/2+p2/2-p_sc), -1)

def clc_rho(p1,p2, c =3e-5, rho_sc = 42):
    return rho_sc*(1+c*(p1/2+p2/2-14.7))

def clc_trans_x(k_x, k_x1, h, h1, delta_x, delta_x1, delta_y,delta_y1, miu, B, beta_c = 1.127e-3):
    Ax = delta_y*h
    Ax1 = delta_y1*h1
    return 2*beta_c/(miu*B)/(delta_x/(Ax*k_x) + (delta_x1/(Ax1*k_x1)))

def clc_trans_y(k_y, k_y1, h, h1, delta_x, delta_x1, delta_y, delta_y1, miu, B, beta_c = 1.127e-3):
    Ay = delta_x*h
    Ay1 = delta_x1*h1
    return 2*beta_c/(miu*B)/(delta_y/(Ay*k_y) + (delta_y1/(Ay1*k_y1)))


def clc_trans_x_(k_x, k_x1, h, h1, delta_x, delta_x1, delta_y,delta_y1, rho, miu, B, beta_c = 1.127e-3):
    Ax = delta_y*h
    Ax1 = delta_y1*h1
    return 2*beta_c*rho/(144*miu*B)/(delta_x/(Ax*k_x) + (delta_x1/(Ax1*k_x1)))

def clc_trans_y_(k_y, k_y1, h, h1, delta_x, delta_x1, delta_y, delta_y1, rho, miu, B, beta_c = 1.127e-3):
    Ay = delta_x*h
    Ay1 = delta_x1*h1
    return 2*beta_c*rho/(144*miu*B)/(delta_y/(Ay*k_y) + (delta_y1/(Ay1*k_y1)))

def clc_gamma(delta_x, delta_y, h, phi, delta_t, ct, B_0 = 1, alpha_c = 5.615):
    return delta_x*delta_y*h*phi*ct/(alpha_c*B_0*delta_t)



def construct_T(params):
	kx = params['kx']
	ky = params['ky']
	h = params['h']
	dx = params['dx']
	dy = params['dy']
	por = params['por']
	order = params['order_mat']
	delta_t = params['time']
	Pre = params['pressure']
	cr_o = params['cr_o']
	cr = params['cr']
	ct = cr_o + cr
	Top = params['top']

	Row, Col = order.shape[0], order.shape[1]
	N = np.zeros((Row,Col))
	S = np.zeros((Row,Col))
	W = np.zeros((Row,Col))
	E = np.zeros((Row,Col))
	N_ = np.zeros((Row,Col))
	S_ = np.zeros((Row,Col))
	W_ = np.zeros((Row,Col))
	E_ = np.zeros((Row,Col))
	C,Q = np.zeros((Row,Col)), np.zeros((Row,Col))


	for i in range(Row):
		for j in range(Col):
			if order[i,j] > 0:
				if order[i-1, j] > 0:
					miu = clc_miu(Pre[i-1, j], Pre[i,j])
					B = clc_B(Pre[i-1, j], Pre[i,j])
					rho = clc_rho(Pre[i-1, j], Pre[i,j])
					N[i,j] = clc_trans_x(kx[i-1,j], kx[i,j], h[i-1,j], h[i,j], dx[i-1,j], dx[i,j], dy[i-1,j], dy[i,j], miu, B)
					N_[i,j] = clc_trans_x_(kx[i-1,j], kx[i,j], h[i-1,j], h[i,j], dx[i-1,j], dx[i,j], dy[i-1,j], dy[i,j], rho, miu, B, beta_c = 1.127e-3)
				if order[i+1, j] > 0:
					miu = clc_miu(Pre[i+1, j], Pre[i,j])
					B = clc_B(Pre[i+1, j], Pre[i,j])
					rho = clc_rho(Pre[i+1, j], Pre[i,j])
					S[i,j] = clc_trans_x(kx[i+1,j], kx[i,j], h[i+1,j], h[i,j], dx[i+1,j], dx[i,j], dy[i+1,j], dy[i,j], miu, B)
					S_[i,j] = clc_trans_x_(kx[i+1,j], kx[i,j], h[i+1,j], h[i,j], dx[i+1,j], dx[i,j], dy[i+1,j], dy[i,j], rho, miu, B, beta_c = 1.127e-3)
				if order[i,j-1] > 0:
					miu = clc_miu(Pre[i, j-1], Pre[i,j])
					B = clc_B(Pre[i, j-1], Pre[i,j])
					rho = clc_rho(Pre[i, j-1], Pre[i,j])
					W[i,j] = clc_trans_y(ky[i,j-1], ky[i,j], h[i,j-1], h[i,j], dx[i,j-1], dx[i,j], dy[i,j-1], dy[i,j], miu, B)
					W_[i,j] = clc_trans_y_(ky[i,j-1], ky[i,j], h[i,j-1], h[i,j], dx[i,j-1], dx[i,j], dy[i,j-1], dy[i,j], rho, miu, B, beta_c = 1.127e-3)
				if order[i,j+1] > 0:
					miu = clc_miu(Pre[i, j+1], Pre[i,j])
					B = clc_B(Pre[i, j+1], Pre[i,j])
					rho = clc_rho(Pre[i, j+1], Pre[i,j])
					E[i,j] = clc_trans_y(ky[i,j+1], ky[i,j], h[i,j+1], h[i,j], dx[i,j+1], dx[i,j], dy[i,j+1], dy[i,j], miu, B)
					E_[i,j] = clc_trans_y_(ky[i,j+1], ky[i,j], h[i,j+1], h[i,j], dx[i,j+1], dx[i,j], dy[i,j+1], dy[i,j], rho, miu, B, beta_c = 1.127e-3)
				gamma = clc_gamma(dx[i,j], dy[i,j], h[i,j], por[i,j], delta_t, ct, B_0 = 1, alpha_c = 5.615)
				C[i,j] = -(N[i,j]+W[i,j]+E[i,j]+S[i,j] +gamma)
				Q[i,j] = N_[i,j]*(Top[i-1,j] - Top[i,j]) + S_[i,j]*(Top[i+1,j] - Top[i,j])+W_[i,j]*(Top[i,j-1] - Top[i,j]) + E_[i,j]*(Top[i,j+1] - Top[i,j]) - gamma*Pre[i,j]

	return N,S,W,E,C,Q

def construct_coeff_mat(cnt, N,S,W,E,C,Q, order):
	Row, Col = order.shape[0], order.shape[1]
	LHS = np.zeros((cnt,cnt))
	for i in range(Row):
		for j in range(Col):
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

	RHS = np.zeros((cnt,1))
	for i in range(Row):
		for j in range(Col):
			if order[i,j] > 0:
				RHS[order[i,j]-1] = Q[i,j]

	return LHS, RHS


def run_simulation(props):
	rock = props['rock']
	fluid = props['fluid']
	reservoir = props['reservoir']
	sim_time = props['time']

	wells = props['well']


	kx = rock.kx
	ky = rock.ky
	h = reservoir.h
	dx = reservoir.dx
	dy = reservoir.dy
	por = rock.poro
	order,cnt = reservoir.active_block(kx)
	time_step = sim_time.time_step
	pressure = reservoir.p_init
	cr_o = fluid.cr_o
	cr = fluid.cr
	top = rock.top




	p_reservoir = []
	p_wells = [[] for i in range(len(wells))]
	q_wells = [[] for i in range(len(wells))]
	mbal1 = []
	mbal2 =[]
	por_cmbc = rock.poro
	B_cmbc = np.zeros((order.shape[0], order.shape[1]))
	for i in range(order.shape[0]):
		for j in range(order.shape[1]):
			if order[i,j] > 0:
				B_cmbc[i,j] = (1 + 3e-5*(reservoir.p_init[i,j] - 14.7))**-1

	for t in sim_time.time_interval:
		params = {'kx': kx, 'ky': ky, 'h': h, 'dx' : dx, 'dy':dy, 'por' : por,\
		'order_mat': order, 'time': time_step, 'pressure':pressure, 'cr_o': cr_o, 'cr': cr, 'top': top}
		N,S,W,E,C,Q = construct_T(params)
		# print(N[4,6])
		# print(clc_gamma(dx[4,6], dy[4,6], h[4,6], por[4,6], 30, 3.1e-5, B_0 = 1, alpha_c = 5.615))


		#Productivity index of well will be run by time.
		#Each time step we change PI
		omega = [0 for i in range(len(wells))]

		if t > 720:
			wells[0].well_state = 'S'
			wells[1].well_state = 'S'

		for i, well in enumerate(wells):
			if well.well_state == 'P':
				ii,jj = well.well_coor
				well_spec = well.well_spec
				well_ppt = well.well_ppt
				omega[i] = well.PI(kx, ky, dx, dy, pressure, h )
				if int(well_spec[0]) == 1:
					C[ii,jj] = C[ii,jj] - omega[i] #for central of LHS
					Q[ii,jj] = Q[ii,jj] - omega[i]*well_spec[1] #for RHS
				else: #it mean well_sepc[0] = 2
					Q[ii,jj] = Q[ii,jj] - well_spec[1]


		B_init = np.zeros((order.shape[0], order.shape[1]))
		for i in range(order.shape[0]):
			for j in range(order.shape[1]):
				if order[i,j] > 0:
					B_init[i,j] = (1 + 3e-5*(pressure[i,j] - 14.7))**-1

		LHS,RHS = construct_coeff_mat(cnt,N,S,W,E,C,Q, order)
		x = solve(LHS, RHS)

		#update pressure and porosity
		Pre_update = np.zeros((order.shape[0], order.shape[1]))
		Por_update = np.zeros((order.shape[0], order.shape[1]))
		for i in range(order.shape[0]):
			for j in range(order.shape[1]):
				if order[i,j] > 0:
					Pre_update[i,j] = x[order[i,j]-1,]
					Por_update[i,j] = por[i,j] + 1e-6*(Pre_update[i,j] - pressure[i,j])*por[i,j]

		p_reservoir.append(Pre_update)
		q_total = 0
		for i, well in enumerate(wells):
			iii, jjj = well.well_coor
			well_spec = well.well_spec
			p_wells[i].append(Pre_update[iii,jjj])
			q = omega[i]*(Pre_update[iii,jjj] - well_spec[1])

			q_wells[i].append(q)
			q_total += q
			# print(q_wells[0])

		# print(q_total)
		#update B
		B_update = np.zeros((order.shape[0], order.shape[1]))
		for i in range(order.shape[0]):
			for j in range(order.shape[1]):
				if order[i,j] > 0:
					B_update[i,j] = (1 + 3e-5*(Pre_update[i,j] - 14.7))**-1


		IMBC = 0
		for i in range(order.shape[0]):
			for j in range(order.shape[1]):
				if order[i,j]>0:
					IMBC += dx[i,j]*dy[i,j]*h[i,j]/5.615*(Por_update[i,j]/B_update[i,j] - por[i,j]/B_init[i,j])
		mbal1.append(abs(IMBC / (q_total*time_step)))

		print("IMBC:", abs(IMBC / (q_total*time_step)))

		CMBC = 0
		q_cmbc = sum(sum(np.asarray(q_wells)))
		#print(q_cmbc)
		for i in range(order.shape[0]):
			for j in range(order.shape[1]):
				if order[i,j]>0:
					CMBC += dx[i,j]*dy[i,j]*h[i,j]/5.615*(Por_update[i,j]/B_update[i,j] - por_cmbc[i,j]/B_cmbc[i,j])
		mbal2.append(abs(CMBC / (q_cmbc*time_step)))
		print("CMBC:", abs(CMBC / (q_cmbc*time_step)))
		#update pressure and porosity
		pressure = Pre_update
		por = Por_update


	return p_reservoir, p_wells, q_wells, mbal1, mbal2

def plot_pressure(t, lo, p_well, color, label):
	plt.plot(t, p_well, color = color, markeredgecolor = color, label = label)
	plt.xlabel("Time (days)")
	plt.ylabel("Pressure (psi) ")
	plt.legend(loc = 'best', prop = dict(size =8))
	plt.title("Pressure profile of well %lo" %lo)
	plt.xlim(0, max(t) + 100)
	plt.ylim(min(p_well) - 100, max(p_well))
	plt.grid(True)
	plt.show()


def plot__multiple_pressure(t, p_wells, color, label):

	for i in range(len(p_wells)):
		plt.plot(t, p_wells[i], color = color[i], markeredgecolor = color[i], label = label[i])

	plt.xlabel("Time (days)")
	plt.ylabel("Pressure (psi) ")
	plt.legend(loc = 'best', prop = dict(size =8))
	plt.title("Pressure profile of all wells ")
	# plt.xlim(0, max(t) + 100)
	# plt.ylim(min(p_wells) - 100, max(p_wells))
	plt.grid(True)
	plt.show()


def plot_production(t, lo, q_well, color, label):
	plt.plot(t, q_well, color = color, markeredgecolor = color, label = label)
	plt.xlabel("Time (days) ")
	plt.ylabel("Production (STB/D)")
	plt.legend(loc = 'best', prop = dict(size = 8))
	plt.title("Production profile of well block %lo" %lo)
	plt.xlim(0, max(t) + 100)
	plt.ylim(min(q_well) - 100, max(q_well))
	plt.grid(True)
	plt.show()

def plot__multiple_production(t, q_wells, color, label):

	for i in range(len(q_wells)):
		plt.plot(t, q_wells[i], color = color[i], markeredgecolor = color[i], label = label[i])

	plt.xlabel("Time (days)")
	plt.ylabel("Production (STB/D) ")
	plt.legend(loc = 'best', prop = dict(size =8))
	plt.title("Production profile of all wells ")
	# plt.xlim(0, max(t) + 100)
	# plt.ylim(min(p_wells) - 100, max(p_wells))
	plt.grid(True)
	plt.show()

def plot_MB_check(t, mbal, color, label):
	plt.plot(t, mbal, color = color, markeredgecolor = color, label = label)
	plt.xlabel("Time (days)")
	plt.ylabel(label)
	plt.legend(loc = 'best', prop = dict(size =8))
	plt.title("Incremental material balance check")
	plt.xlim(0, max(t) + 100)
	plt.ylim(0,2)
	plt.grid(True)
	plt.show()

def plot_CMBC_check(t, mbal, color, label):
	plt.plot(t, mbal, color = color, markeredgecolor = color, label = label)
	plt.xlabel("Time (days)")
	plt.ylabel(label)
	plt.legend(loc = 'best', prop = dict(size =8))
	plt.title("Cumulative material balance check")
	plt.xlim(0, max(t) + 100)
	plt.ylim(0,2)
	plt.grid(True)
	plt.show()

def main():
	start = time.time()
	
	thickness, top, perm, dx, dy, porosity = load_data()
	kx, ky, h, top, dx, dy, poro = data_cleaning(thickness, top, perm, dx, dy, porosity)
	p_init = np.zeros((kx.shape[0], kx.shape[1]))
	for i in range(kx.shape[0]):
		for j in range(kx.shape[1]):
			if kx[i,j] > 0:
				p_init[i,j] = 5000


	time_step = 30
	time_interval = np.arange(0, 360*3+30, time_step)
	sim_time = props_time(time_step = time_step, time_interval = time_interval)

	rock = props_rock(kx = kx, ky = ky, poro = poro, top = top)
	fluid = props_fluid(cr_o = 3e-5, cr = 1e-6)
	reservoir = props_reservoir(dx = dx, dy = dy, h = h, p_init = p_init)
	well1 = props_well(well_coor=(4, 6), well_spec = (1,1000), well_ppt = (0.25, 0), well_state = 'P')
	well2 = props_well(well_coor=(9, 6), well_spec = (1,1000), well_ppt = (0.25, 0), well_state = 'P')
	well3 = props_well(well_coor=(4, 20), well_spec = (1,1000), well_ppt = (0.25, 0), well_state = 'P')
	well4 = props_well(well_coor=(8, 10), well_spec = (1,1000), well_ppt = (0.25, 0), well_state = 'P')
	well5 = props_well(well_coor=(12, 13), well_spec = (1,1000), well_ppt = (0.25, 0), well_state = 'P')
	well6 = props_well(well_coor=(6, 14), well_spec = (1,1000), well_ppt = (0.25, 0), well_state = 'P')
	wells = [well1, well2, well3, well4, well5, well6]
	props = {'rock': rock, 'fluid': fluid, 'reservoir': reservoir, 'well': wells, 'time': sim_time}

	p_reservoir, p_wells, q_wells, mbal1, mbal2 = run_simulation(props)

	end = time.time()
	print("Time taken",end - start)

	color = ['b', 'g', 'r', 'c', 'm', 'y',]
	label = ['Well 1', 'Well 2', 'Well 3', 'Well4', 'Well 5', 'Well 6']

	
	for i in range(len(p_reservoir)):
		plt.pcolormesh(p_reservoir[i], cmap = 'RdBu', vmin = 0, vmax = 5000)
		plt.title("Block pressure at time %d" %(i*30) +"days")
		plt.pause(0.1)
	plt.show()

	#plot pressure and production profile
	plot_MB_check(time_interval, mbal1, color = 'r', label = "IMBC")
	plot_CMBC_check(time_interval, mbal2, color = 'r', label = "CMBC")
	plot__multiple_pressure(time_interval, p_wells, color = color, label = label )
	plot__multiple_production(time_interval, q_wells, color = color, label = label )



if __name__ == '__main__':
	main()
