import numpy as np
import pickle
from landlab.components.diffusion import LinearDiffuserMultiClass
from landlab.components import ExtendedGravelBedrockEroder
from landlab.components.soil_grading import SoilGrading
from landlab.components import FlowAccumulator
import time
from landlab import imshow_grid
import matplotlib.pyplot as plt

# For Alpine:
# sim_id_to_run = int(sys.argv[1])
# start_new_flag = int(sys.argv[2])

# Read the simulation inputs and description.
# "start_new_flag" allows to continue a stopped simulation.
sim_id_to_run = int(21)
start_new_flag = int(1)

# Load the simulation info
with open(f'./Inputs/sim_num_{sim_id_to_run}.pkl', 'rb') as f:
    loaded_dict = pickle.load(f)


# Init components
grid = loaded_dict['grid']


# Init SoilGrading will throw a warning in case "grid" is
# already associate with "soil__depth" field.
# This is fine.
# SoilGrading will rewrite "soil__depth" to match the weight
# of grains (take into account porosity) as given by the key
# "init_grains_weight".
sg = SoilGrading(grid,
                 meansizes=loaded_dict['grain_sizes'],
                 grains_weight=loaded_dict['init_grains_weight'],
                 phi=loaded_dict['porosity'],
                 soil_density=loaded_dict['rho_sed'])

sg.update_bed_grains_proportions(proportions=loaded_dict['proportions'])

fa = FlowAccumulator(grid, runoff_rate=loaded_dict['runoff_rate'])

eroder = ExtendedGravelBedrockEroder(grid,
                                     intermittency_factor=loaded_dict['Intermittency'],
                                     sediment_porosity=loaded_dict['porosity'],
                                     depth_decay_scale=loaded_dict['depth_decay_scale'],
                                     plucking_coefficient=loaded_dict['plucking_coeff'],
                                     abrasion_coefficients=loaded_dict['sediment_abrasion_coefficients'],
                                     bedrock_abrasion_coefficient=loaded_dict['bedrock_abrasion_coefficient'],
                                     fractions_from_plucking=loaded_dict['proportions'],
                                     rho_sed=loaded_dict['rho_sed'],
                                     rho_water=loaded_dict['rho_water'],
                                     fixed_width_flag=loaded_dict['fixed_width_flag'],
                                     fixed_width_coeff=loaded_dict['fixed_width_coeff'],
                                     fixed_width_expt=loaded_dict['fixed_width_expt'],
                                     mannings_n=loaded_dict['mannings_n'],
                                     tau_star_c_median=loaded_dict['tau_c_median'],
                                     alpha=loaded_dict['alpha'],
                                     tau_c_bedrock=loaded_dict['tau_c_bedrock'],
                                     d_min=loaded_dict['d_min'],
                                     plucking_by_tools_flag=loaded_dict['plucking_by_tools']
                                     )

diffuser = LinearDiffuserMultiClass(grid,
                                    linear_diffusivity_soil=loaded_dict['diffusivity_soil'],
                                    linear_diffusivity_rock=loaded_dict['diffusivity_rock'],
                                    rho_sed=loaded_dict['rho_sed'],
                                    phi=loaded_dict['porosity'])


# Pointers
elev = grid.at_node["topographic__elevation"]
rock_elev = grid.at_node["bedrock__elevation"]
soil_depth = grid.at_node["soil__depth"]

topo_mat = loaded_dict['topo_mat']
soil_mat = loaded_dict['soil_mat']
median_size_mat = loaded_dict['median_size_mat']
grains__weight_mat = loaded_dict['grains__weight_mat']

simulated_year = loaded_dict['simulated_year']
saved_index = loaded_dict['saved_index']
total_time = loaded_dict['total_time']

# this part update the fields in case we want to start the simulation
# from where it stopped.
if start_new_flag ==1:
    elev[:] = rock_elev[:] + soil_depth[:]
    topo_mat[:, saved_index] = grid.at_node['topographic__elevation'][grid.core_nodes]
    soil_mat[:, saved_index] = grid.at_node['soil__depth'][grid.core_nodes]
    median_size_mat[:, saved_index] = grid.at_node['median_size__weight'][grid.core_nodes]
    if np.ndim(grid.at_node['grains__weight']) > 1:
        grains__weight_mat[:, :, saved_index] = grid.at_node['grains__weight'][grid.core_nodes]
    else:
        grains__weight_mat[:, saved_index] = grid.at_node['grains__weight'][grid.core_nodes]

else:
    with open(
            f'./Result_dict/sim_id_{sim_id_to_run}.pkl',
            'rb') as file:
        loaded_object2 = pickle.load(file)

    grid2 = loaded_object2['grid']
    grid.at_node['topographic__elevation'][:] = grid2.at_node['topographic__elevation']
    grid.at_node['soil__depth'][:] = grid2.at_node['soil__depth']
    grid.at_node['bedrock__elevation'][:] = grid2.at_node['bedrock__elevation']
    grid.at_node['median_size__weight'][:] = grid2.at_node['median_size__weight']
    grid.at_node['grains__weight'][:] = grid2.at_node['grains__weight']

    topo_mat[:, :] = loaded_object2['topo_mat']
    soil_mat[:, :] = loaded_object2['soil_mat']
    median_size_mat[:, :] = loaded_object2['median_size_mat']
    grains__weight_mat[:, :] = loaded_object2['grains__weight_mat']

    simulated_year = loaded_object2['simulated_year']
    saved_index = loaded_object2['saved_index']
    total_time = loaded_object2['total_time']

saving_cnt = np.copy(loaded_dict['saving_resolution_years'])+simulated_year

# Within the "max_dt" interval, I assume no change in the median grain size.
# Usually, EGBE reduce "max_dt" by much.
max_dt = 100

# Ok, main loop:
t1 = time.time()
for i in range(simulated_year, total_time, max_dt):

    # Diffuse soil and bedrock
    diffuser.run_one_step(dt=max_dt)

    # Update uplift and accordingly, topography
    rock_elev[grid.core_nodes] += loaded_dict['Uplift rate']*max_dt
    elev[grid.core_nodes] = rock_elev[grid.core_nodes] + grid.at_node['soil__depth'][grid.core_nodes]

    # Run routing and EGBE
    fa.run_one_step()
    eroder.run_one_step(max_dt)

    # Update the median grain size
    # Note: "run_one_step" of SoilGrading is fragmentation which we don't need now.
    sg.update_median_grain_size()


    simulated_year+=max_dt
    #print(i)

    # Here we save the output file according to constant temporal interval/resolution.
    # Keep-in-mind: As we care about grains weight of several classes,
    # the model output file could be heavy.
    if i >= saving_cnt:
        print(time.time()-t1)
        saved_index+=1
        saving_cnt += loaded_dict['saving_resolution_years']

        loaded_dict['saved_index'] = saved_index
        loaded_dict['simulated_year'] = simulated_year
        loaded_dict['grid'] = grid
        topo_mat[:, saved_index] = grid.at_node['topographic__elevation'][grid.core_nodes]
        soil_mat[:, saved_index] = grid.at_node['soil__depth'][grid.core_nodes]
        median_size_mat[:, saved_index] = grid.at_node['median_size__weight'][grid.core_nodes]
        if np.ndim(grid.at_node['grains__weight']) > 1:
            grains__weight_mat[:, :, saved_index] = grid.at_node['grains__weight'][grid.core_nodes]
        else:
            grains__weight_mat[:, saved_index] = grid.at_node['grains__weight'][grid.core_nodes]

        with open(f'./Result_dict/sim_id_{sim_id_to_run}.pkl',
                  'wb') as file:
            pickle.dump(loaded_dict, file, pickle.HIGHEST_PROTOCOL)
        imshow_grid(grid,'topographic__elevation'),plt.show()
        t1 = time.time()



