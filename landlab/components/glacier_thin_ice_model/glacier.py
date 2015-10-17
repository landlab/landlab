#! /usr/bin/python

from landlab import Component
import numpy as np
import datetime as datetime
from scipy.sparse import csr_matrix
from scipy.sparse import linalg
import six


class Glacier(Component):

	def __init__(self,grid,dictionary,**kwds):
		'''
		Define physical parameters (here assuming EISMINT-1 values)
		'''
		self.n_GLEN = 3          # Glen's flow law exponent
		self.A_GLEN = 7.5738e-17 #6.05904e-18; Monthly #7.5738e-17 Cuffey & Paterson (4th ed) Glen's law parameter in Pa^{-3} yr^{-1} units (same as A_GLEN=2.4e-24 Pa^{-3} s^{-1})

		self.m_SLIDE = 2        # Sliding law exponent
		self.C_SLIDE = 0    # 1.0e-08;  # 1.0e-06;  # Sliding coefficient in Pa, metre,(Year units)

		self.RHO = 900   # Density (SI units)
		self.g = 9.80    # Gravity (SI units, rho*g has units of Pa)
		self.K_eps = 1.0e-12
		self.OMEGA = 1.5  # 1.6
		super(Glacier,self).__init__(grid)
		self.initialize(dictionary)

	def initialize(self,dictionary,**kwds):
		'''
		Initialize values for calculation:

		S: ice surface ice_elevation
		B: bed elevation
		b_dot: mass of ice added or subtracted from each cell
		dt: time interval
		t_STOP: ending time for modeling
		t: starting time for modeling
		dx: node spacing
		nx: number of columns of nodes
		ny: number of rows of nodes
		N: number of nodes
		'''
		self.S = kwds.pop('S',dictionary['S'])
		self.B = kwds.pop('B',dictionary['B'])
		self.b_dot = kwds.pop('b_dot',dictionary['b_dot'])
		self.dt = kwds.pop('dt',dictionary['dt'])
		self.t_STOP = kwds.pop('t_STOP',dictionary['t_STOP'])
		self.t = kwds.pop('t',dictionary['t'])
		self.dx = kwds.pop('dx',dictionary['dx'])
		self.nx = kwds.pop('nx',dictionary['nx'])
		self.ny = kwds.pop('ny',dictionary['ny'])
		self.N = self.nx * self.ny
		self.setupIndexArrays()

	def step_update(self):
		'''
		calculate S (ice surface elevation) for each timestep

		H_max: maximum ice thickness
		S_max: maximum ice surface elevation
		'''
		self.S, self.t = self.step()
		SB = self.S - self.B
		self.H_max = np.max(SB)
		self.k_H_max = np.argmax(SB)
		self.S_max = np.max(self.S)
		self.k_S_max = np.argmax(self.S)

	def recursive_steps(self):
		'''
		Iterate over each time step to update the ice surface elevation
		'''
		while 1:
			self.step_update()
			# print self.S[0:5]
			self.ALPHA_I = 100*np.sum(self.S > self.B)/float(self.N)

			six.print_('BKS: At t={:8.2f} yr ALPHA_I={:.2f}% and maxima are: H({:d}) = {:f} \
			S({:d})={:f}\n'.format(self.t, self.ALPHA_I, self.k_H_max, self.H_max, self.k_S_max, self.S_max))

			### Stop iterating until the final timestep
			if self.t > self.t_STOP:
				I = np.zeros(self.N)
				I[self.S > self.B] = 1

				# S_map = self.S.reshape(self.ny,self.nx)
				# B_map = self.B.reshape(self.ny,self.nx)
				# I_map = I.reshape(self.ny,self.nx)

				### Note: the difference between python and matlab in matrix orders
				S_map = self.S.reshape(self.nx,self.ny).T
				B_map = self.B.reshape(self.nx,self.ny).T
				I_map = I.reshape(self.nx,self.ny).T
				H_map = S_map - B_map
				self.grid['node']['ice_elevation'] = S_map
				self.grid['node']['B_map'] = B_map
				self.grid['node']['I_map'] = I_map
				self.grid['node']['ice_thickness'] = H_map
				now = datetime.datetime.now().strftime('%H:%M:%S')
				file_str = 'S_map.txt'
				six.print_('main(): Output stored in file "{:s}" at time {:s} \n'.format(file_str,now))
				break

	def step(self):
		'''
		For each timestep, a sparse linear system (Ax = C) need to be solved to update ice surface elevation
		'''

		### update diffusivity for each timestep
		self.diffusion_update()
		D_sum = self.D_IC_jc + self.D_IP_jc + self.D_ic_JC + self.D_ic_JP

		row = np.int64([[self.ic_jc],[self.ic_jc],[self.ic_jc],[self.ic_jc],[self.ic_jc]]).flatten()
		col = np.int64([[self.im_jc],[self.ip_jc],[self.ic_jm],[self.ic_jp],[self.ic_jc]]).flatten()
		val = np.array([[-self.OMEGA * self.D_IC_jc],[-self.OMEGA * self.D_IP_jc],[-self.OMEGA * self.D_ic_JC],[-self.OMEGA * self.D_ic_JP],[1/self.dt + self.OMEGA * D_sum]]).flatten()
		C = (1 - self.OMEGA) * ((self.D_IC_jc * self.S[self.im_jc]) + self.D_IP_jc * self.S[self.ip_jc] + self.D_ic_JC * self.S[self.ic_jm] + self.D_ic_JP * \
			self.S[self.ic_jp]) + (1/self.dt - (1 - self.OMEGA) * D_sum) * self.S[self.ic_jc] + self.b_dot
		C = C.flatten()

		### construct a sparse matrix A
		A = csr_matrix( (val,(row,col)), shape=(self.N, self.N))
		# print 'solving'
		S_out = linalg.spsolve(A,C)
		# print 'solved'

		### ice thickness couldn't be negative, ice surface elevation should not be less than bed elevation
		S_out[S_out < self.B] = self.B[S_out < self.B]

		t_n = self.t + self.dt
		return S_out, t_n

	def diffusion_update(self):
		'''
		calculate diffusivity for each timestep
		'''
		A_tilde = 2 * self.A_GLEN * (self.RHO * self.g) ** self.n_GLEN/(self.n_GLEN + 2)/(self.dx ** 2)
		C_tilde = self.C_SLIDE * (self.RHO * self.g)**self.m_SLIDE/(self.dx**2)
		nm_half = (self.n_GLEN - 1) / 2.0   ### @
		npl = self.n_GLEN + 1
		mm_half = (self.m_SLIDE - 1) / 2.0  ### @
		ml = self.m_SLIDE

		SB = self.S - self.B
		SB[SB<0] = 0
		H = SB

		H_IC_jc = 0.5*(H[self.ic_jc] + H[self.im_jc])
		H_ic_JC = 0.5*(H[self.ic_jc] + H[self.ic_jm])

		H_IC_jc_up = H[self.im_jc]
		H_ic_JC_up = H[self.ic_jm]

		ix = (self.S[self.ic_jc]>self.S[self.im_jc]).reshape(-1)
		H_IC_jc_up[self.S[self.ic_jc]>self.S[self.im_jc]] = H[self.ic_jc[ix]].reshape(-1)

		ix = (self.S[self.ic_jc]>self.S[self.ic_jm]).reshape(-1)
		H_ic_JC_up[self.S[self.ic_jc]>self.S[self.ic_jm]] = H[self.ic_jc[ix]].reshape(-1)

		dS_dx_IC_jc = (self.S[self.ic_jc] - self.S[self.im_jc])/self.dx
		dS_dy_IC_jc = (self.S[self.ic_jp] + self.S[self.im_jp] - self.S[self.ic_jm] - self.S[self.im_jm])/(4*self.dx)
		dS_dx_ic_JC = (self.S[self.ip_jc] + self.S[self.ip_jm] - self.S[self.im_jc] - self.S[self.im_jm])/(4*self.dx)
		dS_dy_ic_JC = (self.S[self.ic_jc] - self.S[self.ic_jm])/self.dx

		S2_IC_jc = np.square(dS_dx_IC_jc) + np.square(dS_dy_IC_jc) + self.K_eps
		S2_ic_JC = np.square(dS_dx_ic_JC) + np.square(dS_dy_ic_JC) + self.K_eps

		if C_tilde == 0:    ### No sliding case
			self.D_IC_jc = A_tilde*H_IC_jc_up*np.power(H_IC_jc,npl)*np.power(S2_IC_jc,nm_half)
			self.D_ic_JC = A_tilde*H_ic_JC_up*np.power(H_ic_JC,npl)*np.power(S2_ic_JC,nm_half)
		elif C_tilde > 0:    ### Sliding case
			self.D_IC_jc = A_tilde*H_IC_jc_up*np.power(H_IC_jc,npl)*np.power(S2_IC_jc,nm_half) \
					+ C_tilde*H_IC_jc_up*np.power(H_IC_jc,ml)*np.power(S2_IC_jc,mm_half)
			self.D_ic_JC = A_tilde*H_ic_JC_up*np.power(H_ic_JC,npl)*np.power(S2_ic_JC,nm_half) \
					+ C_tilde*H_ic_JC_up*np.power(H_ic_JC,ml)*np.power(S2_ic_JC,mm_half)
		else:
			six.print_('diffusion(): C_tilde is undefined or incorrectly defined')

		self.D_IP_jc  = self.D_IC_jc[self.ip_jc]
		self.D_ic_JP  = self.D_ic_JC[self.ic_jp]

	def setupIndexArrays(self):

		ic = np.arange(self.nx)
		ip = np.append(np.array([np.arange(1,self.nx)]),self.nx - 1)
		im = np.append(0,np.array([np.arange(self.nx - 1)]))

		jc = np.arange(self.ny)
		jp = np.append(0,np.array([np.arange(self.ny - 1)]))
		jm = np.append(np.array([np.arange(1,self.ny)]),self.ny - 1)

		ic_jc = np.arange(1, self.N + 1).reshape(self.nx, self.ny)

		self.ip_jc = self.setupArrays(ip,jc,ic_jc) - 1
		self.im_jc = self.setupArrays(im,jc,ic_jc) - 1
		self.ic_jp = self.setupArrays(ic,jp,ic_jc) - 1
		self.ic_jm = self.setupArrays(ic,jm,ic_jc) - 1

		self.im_jm = self.setupArrays(im,jm,ic_jc) - 1
		self.ip_jm = self.setupArrays(ip,jm,ic_jc) - 1
		self.im_jp = self.setupArrays(im,jp,ic_jc) - 1
		self.ip_jp = self.setupArrays(ip,jp,ic_jc) - 1

		self.ic_jc = ic_jc.reshape(-1) - 1


	def setupArrays(self,a, b, ic_jc):
		x,y = np.meshgrid(b,a)
		array = []
		for l in zip(y.ravel(),x.ravel()):
			array.append(ic_jc[l])
		array = np.array(array)
		return array
