
import numpy as np
import scipy.io as io
from landlab.components.glacier_thin_ice_model.glacier import Glacier
from landlab import RasterModelGrid
import matplotlib.pyplot as plt
import matplotlib as mpl

def main():
	'''
	B: bed elevation
	b_dot:
	dx: node spacing (dx = dy)
	nx: number of columns of nodes
	ny: number of rows of nodes
	t_STOP: number of years of simulation
	dt: time step interval, in years
	t: starting time of simulation, default, 0
	'''
	input_file = 'mb4_spin1.mat'
	mat = io.loadmat(input_file)
	B = mat['B']
	b_dot = mat['b_dot']
	dx = mat['dx'][0,0]
	dy = mat['dy'][0,0]
	nx = np.int_(mat['nx'][0,0])
	ny = np.int_(mat['ny'][0,0])

	t_STOP = 500        ### 1000
	dt = 0.08333
	t = 0

	### put input data in a dictionary, and pass the dictionary as arguments
	B,b_dot,S = flatten(B,b_dot)
	dictionary = {'S':S,'B':B,'b_dot':b_dot,'dt':dt,'t_STOP':t_STOP,'t':t,'dx':dx,'nx':nx,'ny':ny}
	grid = RasterModelGrid(nx,ny,dx)
	gla = Glacier(grid,dictionary)
	gla.recursive_steps()

	### save outputs in ascill file
	S_map = gla.grid['node']['ice_elevation'] 	### ice surface elevation matrix
	H_map = gla.grid['node']['ice_thickness']	### ice thickness matrix
	I_map = gla.grid['node']['I_map']			### ice mask matrix

	np.savetxt('S_map.txt',S_map)
	np.savetxt('H_map.txt',H_map)
	np.savetxt('I_map.txt',I_map)

	### plot S_map
	plt.figure(figsize=(8,6))
	plt.imshow(S_map)
	plt.colorbar()
	plt.savefig('S_map_{0}yrs.pdf'.format(t_STOP),dpi=300)

	### plot H_map
	plt.figure(figsize=(8,6))
	plt.imshow(H_map)
	plt.colorbar()
	plt.savefig('H_map_{0}yrs.pdf'.format(t_STOP),dpi=300)

	### plot map of observed and simulated masks of ice
	plot_mask('I_map.txt','obs_map.txt')


def flatten(B,b_dot):
	### flatten two dimensional matrix
	B = B.T.flatten()
	B[np.isnan(B)] = 0
	S = B
	b_dot = b_dot.T.flatten()
	return B,b_dot,S

def plot_mask(ifile_sim,ifile_obs):
	'''
	plot simulated and observed masks of ice

	'''

	# make presence of ice from simulated ice file as 1
	# make presence of ice from observed ice file as 2
	# make presence of ice in overlapping area as 3

	dat_sim = np.genfromtxt(ifile_sim)
	dat_obs = np.genfromtxt(ifile_obs)
	dat_obs[np.where(dat_obs==1)] = 2
	dat_add = dat_sim + dat_obs

	plt.figure(figsize=(10,8))
	# define the colormap
	cmap = plt.cm.jet
	# extract all colors from the .jet map
	cmaplist = [cmap(i) for i in range(cmap.N)]
	# force the first color entry to be grey
	cmaplist[0] = (.5,.5,.5,1.0)
	# create the new map
	cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

	# define the bins and normalize
	bounds = np.linspace(0,4,5)
	norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

	img1 = plt.imshow(dat_sim,cmap=cmap,norm=norm)    ### 1
	img2 = plt.imshow(dat_obs,cmap=cmap,norm=norm)   ### 2
	img3 = plt.imshow(dat_add,cmap=cmap,norm=norm)   ### 3
	cbar = plt.colorbar(img3, ticks=[0.5,1.5, 2.5, 3.5], orientation='vertical')
	cbar.ax.set_yticklabels(['No ice','Simulated Only', 'Observed Only', 'Overlapped'])# horizontal colorbar
	cbar.ax.tick_params(labelsize=12)
	plt.savefig('mask.pdf',dpi=300)

if __name__ == "__main__":
	main()





