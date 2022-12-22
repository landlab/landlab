"""
Plant integration component of GenVeg - this is the part of GenVeg that handles interactions
between plants and plants and the physical grid. 
"""
from logging.config import valid_ident
from matplotlib.pyplot import grid
#from landlab.components import Radiation
from landlab import Component
import numpy as np
import pandas as pd
from scipy.optimize import fsolve
from .growth import PlantGrowth
import warnings
from landlab.data_record import DataRecord

class GenVeg(Component, PlantGrowth):
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
            self._par = self._grid["cell"]["radiation__par_tot"][:].copy()
        except RuntimeWarning:
            msg=('GenVeg requires incoming PAR for each timestep. Empiricial estimation will be used for the run.')
            print(msg)
            #try:
            #    self.__albedo_bare=self._grid['cell']['bare_ground_albedo'][:].copy()
            #except KeyError:
            #   msg=('Empirical estimation of PAR requires bare ground albedo at-cell field.')
            #   raise ValueError(msg)
            #            # From chapter 2 of Teh 2006 (pg. 31; equation 2.13)
            #try:

            #except KeyError:
            #    msg=('Empirical estimation of PAR requires lat/long of grid xy reference.')     
            #try:          
            #    self._clear_sky_index=self._grid['cell']['air__clear_sky_index'][:].copy()
            #except KeyError:
            #    msg=('Empirical estimation of PAR requires a clear sky index value at-cell field')
            #except RuntimeError:
 
            #    raise RuntimeError(msg)
            #else:
            #    self._par_method='empirical_estimation'
        else:
            self._par_method='direct_input'
    
            (_,_latitude)=self._grid.xy_of_reference
            self._lat_rad = np.radians(_latitude)

    #Set initial time variables     
        self.dt=dt
        self.current_day=current_day
        self.start_date=current_day
        self.time_ind=0


    ##Instantiate a plantgrowth object and summarize number of plants and biomass per cell
        #Create empty array to store PlantGrowth objects
        plantspecies=[]
        _ = self._grid.add_zeros('vegetation__total_biomass',at='cell')
        _ = self._grid.add_zeros('vegetation__n_plants',at='cell')
        tot_biomass_all=self._grid.at_cell['vegetation__total_biomass']
        tot_plant_all=self._grid.at_cell['vegetation__n_plants']



        #for each species in the parameters file
        for species in vegparams:
            if species=='null':
                continue
            species_dict=vegparams[species]
            species_obj=PlantGrowth(self._grid,self.dt, self.current_day, self.start_date, species_params=species_dict)
            array_out=species_obj.species_plants()
            plantspecies.append(species_obj)
            #Summarize biomass and number of plants across grid
            tot_bio_species=array_out['root_biomass']+array_out['leaf_biomass']+array_out['stem_biomass']
            cell_biomass=np.bincount(array_out['cell_index'], weights=tot_bio_species, minlength=self._grid.number_of_cells)
            cell_plant_count=np.bincount(array_out['cell_index'], minlength=self._grid.number_of_cells)
            tot_biomass_all=tot_biomass_all+cell_biomass
            tot_plant_all=tot_plant_all+cell_plant_count
        
        self._grid.at_cell['vegetation__total_biomass']=tot_biomass_all
        self._grid.at_cell['vegetation__n_plants']=tot_plant_all

        self.plant_species=plantspecies





    def run_one_step(self):
        tot_biomass_all=np.zeros_like(self._grid.at_cell['vegetation__total_biomass'])
        tot_plant_all=np.zeros_like(self._grid.at_cell['vegetation__n_plants'])
        for species_obj in self.plant_species:
            species_obj._grow(self.current_day)
            array_out=species_obj.species_plants()
            tot_bio_species=array_out['root_biomass']+array_out['leaf_biomass']+array_out['stem_biomass']
            cell_biomass=np.bincount(array_out['cell_index'], weights=tot_bio_species, minlength=self._grid.number_of_cells)
            cell_plant_count=np.bincount(array_out['cell_index'], minlength=self._grid.number_of_cells)
            tot_biomass_all=tot_biomass_all+cell_biomass
            tot_plant_all=tot_plant_all+cell_plant_count
        self._grid.at_cell['vegetation__total_biomass']=tot_biomass_all
        self._grid.at_cell['vegetation__n_plants']=tot_plant_all
        self.current_day +=1


    def view_record_grid(self, ):
        view=self.record_grid.dataset.to_dataframe()
        return view

    
    def save_output(self, save_params=['root_biomass','leaf_biomass','stem_biomass']):
        rel_time=(self.current_day-self.start_date).astype(float)
        for species_obj in self.plant_species:
            species_obj.species_plants()
            species_obj.save_plant_output(rel_time, save_params)
        self.time_ind += 1

    def get_plant_output(self, species='all'):
        if species=='all':
            out_df=pd.DataFrame()
            for species_obj in self.plant_species:
                species_df=species_obj.record_plants.dataset.to_dataframe()
                species_df.reset_index(inplace=True)
                species_df.set_index(['time','vegetation__species','item_id'], inplace=True)
                out_df=pd.concat([out_df,species_df])
        else:
            for species_obj in self.plant_species:
                if species_obj.species_name==species:
                    out_df=species_obj.record_plants.dataset.to_dataframe()
        return out_df