"""
Plant integration component of GenVeg - this is the part of GenVeg that handles interactions
between plants and plants and the physical grid. 
"""
from matplotlib.pyplot import grid
#from landlab.components import Radiation
from landlab import Component
import numpy as np
import pandas as pd
from scipy.optimize import fsolve
from .growth import PlantGrowth
import warnings

class GenVeg(Component,PlantGrowth):
    """
    Add Intro Stuff here
    """
    _name = "GenVeg"

    _unit_agnostic = False

    _cite_as = """
    @article{piercygv,
        author = {Piercy, C.D.; Swannack, T.M.; Carrillo, C.C.; Russ, E.R.; Charbonneau, B. M.]
    }
    """
    #Add all variables to be saved or chosen as outputs here
    _info = {
        "plant__total_biomass": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units":"g",
            "mapping":"cell",
            "doc": "Total plant biomass for the plant class at the end of the time step"
        },
    }

    def __init__(
        self,
        grid,
        dt,
        current_day,
        vegparams
    ):

    #save grid object to class
        super().__init__(grid)
    ##Instantiate a plantgrowth object
        PlantGrowth.__init__(self,self._grid,vegparams)
        
        self.dt=dt
        self.current_day=current_day

        #Check to see if there are plants on the grid
        try:
            self.plants_on_grid=self._grid['cell']['vegetation__plant_species']
        except KeyError:
            msg = ("GenVeg requires initial distribution of plant species at-cell field.")
            raise ValueError(msg)   
        #Check to see if grid contains required environmental fields
        try:
            self._air_temp = self._grid['cell']['air__temperature_C'][:].copy()
        except KeyError:
            msg=('GenVeg requires air temperature in Celcius for each time step.')
            raise ValueError(msg)
        
        try:
            self._radiation = self._grid["cell"]["radiation__net_flux"][:].copy()
        except KeyError:
            msg=('GenVeg requires incoming radiation flux for each time step')
            raise ValueError(msg)

        #Count number of plants per cell
        n_plants=np.count_nonzero(np.where(self._grid.at_cell['vegetation__plant_species']=='null',0,1), axis=1)
        _ = self._grid.add_field('vegetation__n_plants',n_plants, at='cell')

        #Calculate initial total biomass per cell
        self._grid.add_empty('vegetation__total_biomass',at='cell')
        biomass_plants=self.plants.groupby(by='cell_index', as_index=False).agg('sum')
        biomass_plants['total_biomass']=biomass_plants['leaf_biomass']+biomass_plants['stem_biomass']+biomass_plants['root_biomass']
        cells_with_plants=biomass_plants.cell_index.tolist()
        self._grid.at_cell['vegetation__total_biomass'][cells_with_plants]=biomass_plants.total_biomass.tolist()


    def run_one_step(self, dt):
        new_plants=self._grow(self.current_day, self.dt)
        n_plants=new_plants.groupby(by='cell_index', as_index=False).agg('count')
        biomass_plants=new_plants.groupby(by='cell_index', as_index=False).agg('sum')
        biomass_plants['total_biomass']=biomass_plants['leaf_biomass']+biomass_plants['stem_biomass']+biomass_plants['root_biomass']
        cells_with_plants=biomass_plants.cell_index.tolist()
        self._grid.at_cell['vegetation__total_biomass'][cells_with_plants]=biomass_plants.total_biomass.tolist()
        #self._grid.at_cell['vegetation__n_plants'][cells_with_plants]=n_plants.counts.tolist()

        
        self.current_day +=1

    def plant_ID(self):
        return self.plants

