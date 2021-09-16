"""
Growth component of GenVeg - this is the main driver of vegetation growth and
is driven by a photosynthesis model. Vegetation growth depends on the availability
of carbohydrate produced by photosynthetically active plant parts. 
"""
from matplotlib.pyplot import grid
from landlab.components import Radiation
from landlab import Component

class PlantGrowth(Component):
    """
    Add Intro Stuff here
    """
    _name = "PlantGrowth"

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
        plant_ids,
        growthparams={
            'ftype': {'Corn':['1']},
            'ph_type': {'Corn':['C4']},
            'gs_start': {'Corn':[65]},
            'gs_end': {'Corn':[270]},
            'gs_length': {'Corn':[205]},
            'senes': {'Corn':[209]},
            'res_co': {'Corn':[0.15]},
            'glu_req': {'Corn':[55]},
            'le_k': {'Corn':[30]},
            'hi': {'Corn':[0.1]},
            'p_max': {'Corn':[2.5]}
        },
        sizeparams={},
        dispparams={},
        storparams={},
        colparams={},
        mortparams={},
        dt=1,
        current_day=0 
    ):
        """Instantiate PlantGrowth
        Parameters
        ----------
        grid: RasterModelGrid
            A Landlab ModelGrid

        plant_ids: dict, required,
            names of plant communities to be simulated

        growthparams: dict, required,
            dictionary of named growth parameters with below keys.
            A separate dict is stored for each type of plant with 
            the key as the unique plant type ID as a str.

            ftype: dict, int, optional,
                functional type of plant class
            ph_type: dict, string, required,
                photosynthesis type (option are: C3 or C4)
            gs_start: dict, int, required,
                growing season start day of year
            gs_end: dict, int, required,
                growing season end day of year
            gs_length: dict, int, required,
                growing season length in days
            senes: dict, int, required,
                start of senescence period after plant reaches peak biomass
            res_co: dict, float, required,
                respiration coefficient
            glu_req: dict, float, required,
                glucose requirement
            le_k: dict, float, required,
                light extinction coefficient
            hi: dict, float, required,
                something
            p_max: dict, float, required,
                maximum photosyntehtic output
            allocate: dict, float <1, required
                3-element list of proportion of biomass in roots, stems, and leaves
                Values must be between 0-1 and all sum to 1

        dt: int, required,
            time step interval

        current_day: int, required,
            day of simulation
        """
        
        super().__init__(grid)

        #Check to see if grid contains required vegetation fields
        try:
            self._last_veg_biomass = self._grid["cell"]["vegetation__live_biomass"][:].copy()
        except KeyError:
            msg =("GenVeg requires initial vegetation biomass coverage as an at-cell field.")
            raise ValueError(msg)
        
        try:
            self._last_veg_ftype = self._grid["cell"]["vegetation__plant_functional_type"][:].copy()
        except KeyError:
            msg =("GenVeg requires initial vegetation functional type distribution as an at-cell field.")
            raise ValueError(msg)
        
        self.plant_ids=plant_ids
        self.num_ids=len(self.plant_ids)
        self.growth_params=growthparams
        self.dt=dt
        self.current_day=current_day

        #Verify that all named plant classes have defined growth parameters
        self._check_params(self.growth_params)

        #Calculate calculated parameters for speed and error check growth parameter values
        new_item={}
        for item in self.plant_ids:
            end=self.growth_params['gs_end'][item][0]
            begin=self.growth_params['gs_start'][item][0]
            length=end-begin+1
            new_item.update({item:[length]})
            allo=self.growth_params['allocate'][item]
            allosum=0
            for vals in allo:
                if vals <= 1 and vals >=0:
                    allosum += vals
                else:
                    msg='Biomass allocation proportions for plant ID ' + item +' must be between 0 and 1'
                    raise ValueError(msg)
            if allosum != 1.0:
                msg='Biomass allocation for plant ID ' +item +' must sum to 1'
                raise ValueError(msg)
        
        self.growth_params['gs_length']=new_item


    def _check_params(self,check):
        """Check to make sure the number of parameters is consistent
        with the number of vegetation classes in the parameter dictionary
        """
        for items in check:
            num_entries=len(check.get(items).keys())
            if num_entries != self.num_ids:
                msg = ('The number of values in ' + items + ' does not'
                'match the number of plant types')
                raise ValueError(msg)
        
    def run_one_step(self):
        print('Inside growth')

        pass
