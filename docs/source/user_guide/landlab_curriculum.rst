.. _wickert_curriculum:

================================
Landscape Evolution with Landlab
================================

by Andy Wickert, University of Minnesota

The Code
~~~~~~~~

::

   # Python packages, including LandLab; you may safely ignore these.
   ################################################################################
   import numpy as np
   from landlab import RasterModelGrid, HexModelGrid
   from landlab.components import StreamPowerEroder, FlowRouter, \
       PrecipitationDistribution, LinearDiffuser, DepressionFinderAndRouter
   from landlab import imshow_grid
   from copy import deepcopy
   from matplotlib import pyplot as plt
   from matplotlib import rc
   rc('text', usetex=True)
   ################################################################################

   #####################
   # MAKE CHANGES HERE #
   #####################

   # Toggle plotting of topogrpahy at all times
   plotAllTimes = True

   # Change these to explore impacts of different process rates
   uplift_rate = 0.001 # [m/yr]
   K_sp = 1.e-5 # stream power coefficient [meters^(1-2*m_sp)/year]
   K_hs = 0.05 # hillslope coefficient [m^2/year]

   # Change this to change the size of the domain
   ncells_side = 60 # number of raster cells on each side
   dxy  = 200 # side length for raster model cells [m]

   # Choose 1 model grid
   mg = RasterModelGrid((ncells_side, ncells_side), dxy)
   #mg = HexModelGrid(ncells_side, ncells_side, dxy, shape='rect')
   #mg = HexModelGrid(ncells_side, ncells_side/2, dxy, shape='hex')

   # Change this to investigate landscape disequilibrium
   tmax = 1E5

   # More advanced options: coupling time step, exponents for the stream power law
   dt = 20000
   m_sp = 0.5
   n_sp = 1.

   ################################################################################

   if plotAllTimes:
       plt.ion()
       figIter = plt.figure('EvolvingTopography')

   t = np.arange(0, tmax+1, dt)
   title_text_z = 'Topography after '+str(tmax/1000.)+' kyr\n'+ '$K_{sp}$='+str(K_sp) + '; $K_{hs}$='+str(K_hs) + '; $dx$='+str(dxy)
   title_text_SA = 'Slope-area relationship after '+str(tmax/1000.)+' kyr\n'+ '$K_{sp}$='+str(K_sp) + '; $K_{hs}$='+str(K_hs) + '; $dx$='+str(dxy)

   gridlist = []

   # add initial noise to produce convergent flow from the initial conditions
   np.random.seed(0) # so our figures are reproducible
   mg_noise = np.random.rand(mg.number_of_nodes)/1000.

   # set up the input fields
   zr = mg.add_zeros('node', 'topographic__elevation')
   zr += mg_noise

   # Landlab sets fixed elevation boundary conditions by default. This is
   # what we want, so we will not modify these here.

   # instantiate the components:
   frr = FlowRouter(mg) # water__unit_flux_in gets automatically ingested
   spr = StreamPowerEroder(mg, K_sp=K_sp, m_sp=m_sp, n_sp=1, threshold_sp=0,
                           use_Q=None)
   lake = DepressionFinderAndRouter(mg)

   # Hillslopes
   dfn = LinearDiffuser(mg, linear_diffusivity=K_hs)

   for ti in t:
       #last_elev = mg.copy()
       if plotAllTimes:
           plt.clf()
           imshow_grid(mg, 'topographic__elevation', grid_units=('m', 'm'),
                           var_name='Elevation (m)')
           plt.xlabel('$x$ [m]', fontsize=16)
           plt.ylabel('$y$ [m]', fontsize=16)
           plt.title('Topography at ' + str(ti/1000) + ' kyr')
           plt.pause(0.05)
       #zr[mg.core_nodes] += uplift_rate*dt
       # Uncomment the two lines below, and comment out the line above, to create
       # two blocks with different uplift rates
       zr[mg.core_nodes[mg.core_nodes >= 1800]] += uplift_rate*dt / 2.
       zr[mg.core_nodes[mg.core_nodes < 1800]] += uplift_rate*dt
       dfn.run_one_step(dt) # hillslopes always diffusive, even when dry
       frr.run_one_step()
       lake.map_depressions()
       spr.run_one_step(dt, flooded_nodes=lake.lake_at_node)
       print (ti/1000, 'kyr elapsed;', str(100*ti/tmax) + '%')

   plt.ioff()
   #plt.savetxt('landlab_topo.txt', )
   # Do some plotting. First the topography:
   plt.figure('topo')
   imshow_grid(mg, 'topographic__elevation', grid_units=('m', 'm'),
                   var_name='Elevation (m)')
   plt.xlabel('$x$ [m]', fontsize=16)
   plt.ylabel('$y$ [m]', fontsize=16)
   plt.title(title_text_z, fontsize=16)
   plt.tight_layout()

   edge = np.unique(mg.neighbors_at_node[mg.boundary_nodes, :])
   not_edge = np.in1d(mg.nodes.flatten(), edge, assume_unique=True,
                          invert=True)
   plt.figure('S-A')
   plt.loglog(mg.at_node['drainage_area'][not_edge],
              mg.at_node['topographic__steepest_slope'][not_edge], 'x')
   #xlim([1.e3, 1.e7])
   plt.ylabel('Topographic slope', fontsize=16)
   plt.xlabel('Drainage area [m$^2$]', fontsize=16)
   plt.tight_layout()

   plt.show()

The Assignment
~~~~~~~~~~~~~~

[[download pdf \| images/landscape_evolution_assignment_AW.pdf ]]
