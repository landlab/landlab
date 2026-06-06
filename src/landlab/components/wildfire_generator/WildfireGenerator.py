#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulation of Wildfires on a 2-D model grid

@author: matheuswanderleydealmeida
"""

import numpy as np
import pandas as pd
import warnings
import random
warnings.filterwarnings("ignore")

from landlab import Component

# =============================================================================
# Defining Classes
# =============================================================================        

class WildfireGenerator(Component):
    """
    Simulate stochastic wildfires of differente sizes and severity, and 
    calculate post-fire vegetation recovery based on aridity. 
    
    Landlab component that calculates that simulates landscape-scale fire 
    activity driven by global climate-vegetation interactions, following
    the theory of Pausas and Paula (2012) and Pausas and Ribeiro (2013). 
    
    See the publication:
        
    de Almeida M., Shobe C.M., Roda-Boluda D.C., Gourbet L., Veraverbeke S., 
    Distelbrink A., Campforts B. (2026) FireLands 1.0: A landscape evolution 
    model for simulating the effects of fires and post-fire erosion on sediment
    dynamics in evolving landscapes. Geosci Model Dev: 
    
    
    Examples
    --------
    
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import WildfireGenerator
    
    >>> np.random.seed(5000)
    
    >>> dt = 1
    >>> mg = RasterModelGrid((5, 5), xy_spacing=1.0)
    >>> z = mg.add_zeros("topographic__elevation", at="node")
    >>> da = mg.add_zeros("drainage_area", at="node")
    >>> fuel = mg.add_zeros("fuel_availability", at="node")
    >>> fuel[:] = np.random.rand(mg.number_of_nodes)*0.75
    
    >>> wg = WildfireGenerator(mg)
    >>> wg.run_one_step(dt)
    >>> len(wg.fire_log['year'])
    >>> 20
    
    
    References
    ----------
    **Required Software Citation(s) Specific to this Component**
    
    de Almeida M., Shobe C.M., Roda-Boluda D.C., Gourbet L., Veraverbeke S., 
    Distelbrink A., Campforts B. (2026) FireLands 1.0: A landscape evolution 
    model for simulating the effects of fires and post-fire erosion on sediment
    dynamics in evolving landscapes. Geosci Model Dev: 
        
        
    **Additional References**

    None Listed

    """

    _name = "WildfireGenerator"

    _unit_agnostic = True
    
    _info = {
        "topographic__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
            },
        
        "drainage_area":{
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m**2",
            "mapping": "node",
            "doc": "Upstream accumulated surface area contributing to the node's discharge",
            },
        
        "fuel_availability":{
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Fuel available above the land surface ",
            }
        
        }
    
    
    
    def __init__(self,
                 grid,
                 potential_fires=100,
                 dt=1,
                 dx=1,
                 vegetation="fuel_availability",
                 riv_min = 2e5,
                 topo="topographic__elevation",
                 alpha = 0.3,
                 aridity=0.5,
                 sev_exponent=0.64):
        self.grid = grid
        self.potential_fires = potential_fires
        self.dt = dt
        self.dx = dx
        self.vegetation = grid.at_node[vegetation]
        self.riv_min = riv_min
        self.rivers = (grid.at_node["drainage_area"] > self.riv_min)
        self.max_vegetation = 1
        self.topo = grid.at_node[topo]
        self.alpha = alpha
        self.aridity = aridity
        self.sev_exponent = sev_exponent
        self.burned = []
        self.last_fire_time = None
        self.fire_log = pd.DataFrame(columns=["year",
                                              "fire_size (km2)",
                                              "center_X",
                                              "center_Y",
                                              "severity_factor",
                                              "aridity",
                                              "changed_nodes",
                                              "pct of vegetation removed"])
        self.regrowth_table = [
            (0.0, 0.1, 5.4, 0.42),  # Trees + shrubs
            (0.1, 0.3, 8.0, 0.29),  # Dense shrubs
            (0.3, 0.5, 12.4, 0.18), # Dense shrubs 
            (0.5, 0.7, 15.0, 0.15), # Sparse shrubs
            (0.7, 1.0, 18.5, 0.12), # Sparse shrubs 
        ]
        
        self.aridity_weights = [
            (0.0, 0.1, 0.5, 0.5), # Trees + shrubs
            (0.1, 0.3, 0.5, 0.5), # Dense shrubs
            (0.3, 0.5, 0.5, 0.5), # Dense shrubs 
            (0.5, 0.7, 0.8, 0.2), # Sparse shrubs
            (0.7, 1.0, 0.9, 0.1), # Sparse shrubs 
        ]
        
    @property
    def fire_sizes(self):
        return self.fire_log["fire_size (km2)"]
    
    @property
    def burned_nodes(self):
        return self.fire_log["changed_nodes"]
    
    @property
    def severity(self):
        return self.fire_log["severity_factor"]
    
    
    def _get_regrowth_params(self):
        for low, high, t, e in self.regrowth_table:
            if low <= self.aridity < high or (self.aridity == 1 and high == 1.0):
                return t, e
        raise ValueError(f"Aridity {self.aridity} is out of range [0,1].")
        
    def get_aridity_weights(self):
        for low, high, f_weight, a_weight in self.aridity_weights:
            if low <= self.aridity < high or (self.aridity == 1 and high == 1.0):
                return f_weight, a_weight
        raise ValueError(f"Aridity {self.aridity} is out of range [0,1].")

        
    def regrow_vegetation(self, t):
        regrowth_time, regrowth_exponent = self._get_regrowth_params()

        delta = 1 - np.exp(-regrowth_exponent / regrowth_time)
        
        # Update the vegetation grid
        self.vegetation[:] = np.minimum(
            self.vegetation + delta,
            self.max_vegetation
            )
    
# method to get the active neighbor nodes (for fire spreading)
    def get_neighbors(self, node):
        
        adj_nodes = self.grid.active_adjacent_nodes_at_node[node]  # D4 approach
        
        return adj_nodes
    
        # adj_nodes = np.concatenate(      # Trying to make it D8
        #     (
        #         self.grid.active_adjacent_nodes_at_node[node],
        #         self.grid.diagonal_adjacent_nodes_at_node[node]
        #         )
        #     )
        # return adj_nodes
        
    
# method to calculate severity based on aridity 
    def calc_severity(self):
        aridity_index = np.random.normal(self.aridity, 0.05)                    # Select a random number based on the 
        aridity_index = np.clip(aridity_index, 0, 1)                            # Aridity cannot be smaller than 0 or higher than 1
        
        severity = np.power(aridity_index, self.sev_exponent)                   # Calculate severity based on the power-law proposed by Grunig et al. (2023)
        severity = np.minimum(severity, 1)                                      # Severity cannot be higer than 1
        severity = np.round(severity, 2)                                        # Rounding severity to 2 decimals
        
        return severity

# method that simulate fire ignition and spreading throughout the grid
    def fire(self, t):
        
        fire_results = []
        
        # Initialize last fire times (if first fire event)
        if self.last_fire_time is None:
            self.last_fire_time = np.full(self.grid.number_of_nodes, -np.inf)

        fire_ignitions = np.random.poisson(self.potential_fires * self.dt)      # define how many fires the model will try to ignite (not all fires will ignite) 
        if fire_ignitions == 0:
            return set()

        all_burned = set()                                                      # Collect all burned nodes from all fires 

        for _ in range(fire_ignitions):                                         # Simulate each ignition event
            
            fuel_weight, aridity_weight = self.get_aridity_weights()
            
            center = random.randint(0, self.grid.number_of_nodes - 1)           # Select a random node to ignite the fire
            if self.rivers[center]:
                continue
            severity_factor = self.calc_severity()                                                                   
            if np.random.rand() > (
                    ((self.vegetation[center]*fuel_weight) +                    # to start the fire: a random number must be smaller than the vegetation value for that node
                     ((self.aridity)*aridity_weight))
                    ):                                                          
                continue                                                        # Skip fire due to low ignition probability

            fire_front = set([center])                                          # Fire starts at the "center" node
            burned = set()                                                      # Track all the nodes affected by this fire

            while fire_front:                                                   # Spread fire iteratively from the front
                new_front = set()
                for node in fire_front:
                    if node in burned:                                          # Skips if the node is already burned by another fire in this timestep
                        continue
                    burned.add(node)
                    for neighbor in self.get_neighbors(node):                   # Looking at the orthogonal neighbors 
                        if neighbor in burned:
                            continue
                        if self.rivers[neighbor]:                                  # stop the fire front if the current fire reaches a firebreak river
                            continue
                        slope = ((self.topo[neighbor] - self.topo[node])/self.dx)  # calculate the slope between the current node and its neighbors 
                        slope_factor = 1 / (1+ np.exp(-self.alpha * slope))        # calculate the slope factor
                        # Calculating the probability spread below 
                        spread_prob = ( 
                            ((self.vegetation[neighbor]*fuel_weight) + 
                            ((self.aridity)*aridity_weight)) * slope_factor
                            )
                        
                        if np.random.rand() < spread_prob:                         # the fire spreads if a random uniform number is smaller than the probability spread
                            new_front.add(neighbor)

                fire_front = new_front

            if not burned:
                continue
            
            changed_nodes = list(burned)

            veg_before = self.vegetation[changed_nodes].copy()
            self.vegetation[changed_nodes] *= (1 - severity_factor)
            veg_after = self.vegetation[changed_nodes].copy()
            veg_change = (1 - np.mean(veg_after/veg_before))*100
            self.last_fire_time[changed_nodes] = t

            fire_area_km2 = len(changed_nodes) * self.dx**2 / 1e6
            s1, s2 = self.grid.x_of_node[center], self.grid.y_of_node[center]

            fire_results.append(
                {"burned_nodes": changed_nodes,
                 "severity_factor": severity_factor})

            self.fire_log = pd.concat([self.fire_log, pd.DataFrame([{           # add all the fire metrics in the DataFrame called "fire_log"
                "year": t,
                "fire_size (km2)": fire_area_km2,
                "center_X": s1,
                "center_Y": s2,
                "severity_factor": severity_factor,
                "aridity": self.aridity,
                "changed_nodes": changed_nodes,
                "pct of vegetation removed": veg_change
            }])], ignore_index=True)

            all_burned.update(burned)

        return fire_results

    def run_one_step(self, dt):         # run one step function to simulate fires and vegetation regrowth 
        self.fire_info = self.fire(dt)
        self.regrow_vegetation(dt)

    
    
    
