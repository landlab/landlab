from six.moves import range
import numpy as np
from matplotlib.pyplot import show, plot, figure, loglog
from landlab import RasterModelGrid, CLOSED_BOUNDARY, imshow_grid_at_node
from landlab.components import FlowRouter, SedDepEroder
from landlab.components import DepressionFinderAndRouter
from landlab.plot import channel_profile as prf

nx = 30
dx = 200.
mg = RasterModelGrid((nx, nx), dx)  # ((10, 3), 200.)
for edge in (mg.nodes_at_left_edge, mg.nodes_at_top_edge,
             mg.nodes_at_right_edge):
    mg.status_at_node[edge] = CLOSED_BOUNDARY

z = mg.add_zeros('node', 'topographic__elevation')
th = mg.add_zeros('node', 'channel_sediment__depth')

m_sp = 0.
n_sp = 0.
m_t = 1.5
n_t = 1.
K_sp = 0.01  # 1.e-5
K_t = 1.e-7  # 1.e-5
dth = 0.1

fr = FlowRouter(mg)
pit = DepressionFinderAndRouter(mg)
sde = SedDepEroder(mg, K_sp=K_sp, m_sp=m_sp, n_sp=n_sp,
                   sed_dependency_type='almost_parabolic',
                   Qc='power_law', m_t=m_t, n_t=n_t, K_t=K_t)

# z[:] = mg.node_y/1000.  # roughness is essential to drive dynamics
z[:] = (mg.node_y + dx*np.random.rand(mg.number_of_nodes))/1000.
th += dth

initz = z.copy()

dt = 1000.  # 1000.
nt = 20000  # 12000
up = 0.0005  # 0.0005
print('m_sp, n_sp, m_t, n_t, K_sp, K_t, dth, dt, nt:')
print(m_sp)
print(n_sp)
print(m_t)
print(n_t)
print(K_sp)
print(K_t)
print(dth)
print(dt)
print(nt)
for i in range(nt):  # 12000
    z[mg.core_nodes] += up*dt
    fr.run_one_step()
    # pit.map_depressions()
    sde.run_one_step(dt)
    if i % 1000 == 0:
        print(i)
        figure(1)
        # plot(z.reshape((nx, nx))[:-1, 1])
        profile_IDs = prf.channel_nodes(
            mg, mg.at_node['topographic__steepest_slope'],
            mg.at_node['drainage_area'],
            mg.at_node['flow__receiver_node'])
        dists_upstr = prf.get_distances_upstream(
            mg, len(mg.at_node['topographic__steepest_slope']),
            profile_IDs, mg.at_node['flow__link_to_receiver_node'])
        plot(dists_upstr[0], z[profile_IDs[0]])
    th += dth

profile_IDs = prf.channel_nodes(mg, mg.at_node['topographic__steepest_slope'],
                                mg.at_node['drainage_area'],
                                mg.at_node['flow__receiver_node'])
dists_upstr = prf.get_distances_upstream(
    mg, len(mg.at_node['topographic__steepest_slope']),
    profile_IDs, mg.at_node['flow__link_to_receiver_node'])

A = mg.at_node['drainage_area']
S = mg.at_node['topographic__steepest_slope']
QsQc = mg.at_node['channel_sediment__relative_flux']

logA = np.log(A[mg.core_nodes])
logS = -np.log(S[mg.core_nodes])

figure(2)
plot(logA, logS, 'x')

# fit trendline
p = np.polyfit(logA, logS, deg=1)
line = np.poly1d(p)
plot(logA, line(logA))
print(p)

figure(3)
imshow_grid_at_node(mg, z)

figure(4)
plot(mg.node_y.reshape((nx, nx))[:-1, 1],
     mg.at_node['channel_sediment__relative_flux'].reshape((nx, nx))[
        :-1, 1:-1], 'y')
plot(dists_upstr[0], QsQc[profile_IDs[0]])
plot(mg.node_y.reshape((nx, nx))[:-1, 1],
     mg.at_node['channel_sediment__relative_flux'].reshape((nx, nx))[
        :-1, 1:-1].min(axis=1))
