# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 15:06:50 2023

@author: LaurentRoberge
"""

import numpy as np

from landlab import Component, NodeStatus
from landlab.utils.return_array import return_array_at_node


class ConcentrationTrackerForSpace(Component):

    """This component tracks the concentration of any user-defined property of
    sediment using a mass balance approach in which the concentration :math:`C`
    is calculated as:

    .. math::
        
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~UPDATE BELOW~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ∂CH/∂t = [-(∂q_x C_x)/∂x - (∂q_y C_y)/∂y] + C_br*H_brw + PH + DH
        
    where :math:`H` is sediment depth, :math:`q_x` and :math:`q_y` are sediment
    fluxed in the x and y directions, :math:`C_br` is concentration in parent 
    bedrock, :math:`H_brw` is the height of bedrock weathered into soil, 
    :math:`P` is the local production rate, :math:`D` is the local decay rate.
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~UPDATE ABOVE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    NOTE: This component requires the sediment__influx and sediment__outflux 
    fields calculated by either the Space or SpaceLargeScaleEroder component.
    This component must be run after every Space or SpaceLargeScaleEroder step
    and before any other flux component. For hillslope sediment tracking see 
    ConcentrationTrackerForDiffusion.

    Examples
    --------
    A 1-D stream channel:
        
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import SpaceLargeScaleEroder
    >>> from landlab.components import ConcentrationTrackerForSpace
    >>> mg = RasterModelGrid((3, 5),xy_spacing=10.)
    >>> mg.set_status_at_node_on_edges(right=4, top=4, left=4, bottom=4)
    >>> mg.status_at_node[5] = mg.BC_NODE_IS_FIXED_VALUE
    >>> c = mg.add_zeros('sediment_property__concentration', at='node')
    >>> h = mg.add_ones("soil__depth", at="node")
    >>> z_br = mg.add_zeros("bedrock__elevation", at="node")
    >>> z = mg.add_zeros("topographic__elevation", at="node")
    >>> c[8] += 1
    >>> z_br += mg.node_x/100
    >>> z += z_br + h
    >>> sp = SpaceLargeScaleEroder(mg)
    >>> ct = ConcentrationTrackerForSpace(mg)
    >>> sp.run_one_step(1.)
    >>> ct.run_one_step(1.)
    >>> np.allclose(mg.at_node["topographic__elevation"][mg.core_nodes],
    ...             np.array([4.11701964, 8.01583689, 11.00247875]))
    True
    >>> np.allclose(mg.at_node["sediment_property__concentration"][mg.core_nodes],
    ...             np.array([0., 0.24839685, 1.]))
    True
    
    Now, a 2-D landscape with stream channels.

    >>> ... EXAMPLE  
    
    
    
    >>> mg = RasterModelGrid((3, 6))
    >>> mg.set_status_at_node_on_edges(right=4, top=4, left=4, bottom=4)
    >>> mg.status_at_node[5] = mg.BC_NODE_IS_FIXED_VALUE
    >>> c = mg.add_zeros('sediment_property__concentration', at='node')
    >>> h = mg.add_zeros("soil__depth", at="node")
    >>> z_br = mg.add_zeros("bedrock__elevation", at="node")
    >>> z = mg.add_zeros("topographic__elevation", at="node")
    >>> c[9] += 1
    >>> h += 1
    >>> z_br += mg.node_x
    >>> z += z_br + h
    >>> sp = SpaceLargeScaleEroder(mg)
    >>> ct = ConcentrationTrackerForSpace(mg)
    >>> sp.run_one_step(1.)
    >>> ct.run_one_step(1.)
    >>> np.allclose(mg.at_node["topographic__elevation"][mg.core_nodes],
    ...             np.array([4.11701964, 8.01583689, 11.00247875]))
    True
    >>> np.allclose(mg.at_node["sediment_property__concentration"][mg.core_nodes],
    ...             np.array([0., 0.24839685, 1.]))
    True
   
    Now, a 2-D landscape with stream channels.

    >>> mg = RasterModelGrid((7, 7),xy_spacing=2.)
    >>> c = mg.add_zeros('sediment_property__concentration', at='node')
    >>> h = mg.add_zeros("soil__depth", at="node")
    >>> z_br = mg.add_zeros("bedrock__elevation", at="node")
    >>> z = mg.add_zeros("topographic__elevation", at="node")
    >>> c[12] += 1
    >>> h += 1
    >>> z_br += 2
    >>> z += z_br + h
    >>> sp = SpaceLargeScaleEroder(mg)
    >>> ct = ConcentrationTrackerForDiffusion(mg)
    >>> sp.run_one_step(1.)
    >>> ct.run_one_step(1.)
    >>> np.allclose(mg.at_node["topographic__elevation"][mg.core_nodes],
    ...             np.array([6.        ,  7.13533528,  6.        ,
    ...                       7.13533528,  8.27067057,  7.13533528,
    ...                       6.        ,  7.13533528,  6.         ]))
    True
    >>> np.allclose(mg.at_node["sediment_property__concentration"][mg.core_nodes],
    ...             np.array([0.        ,  0.38079708,  0.        ,
    ...                       0.38079708,  1.        ,  0.38079708,
    ...                       0.        ,  0.38079708,  0.         ]))
    True
    
    And running one more step.
    
    >>> ddd.run_one_step(1.)
    >>> ct.run_one_step(1.)
    >>> np.allclose(mg.at_node["topographic__elevation"][mg.core_nodes],
    ...             np.array([5.52060315,  6.62473963,  5.52060315,
    ...                       6.62473963,  8.00144598,  6.62473963,
    ...                       5.52060315,  6.62473963,  5.52060315 ]))
    True
    >>> np.allclose(mg.at_node["sediment_property__concentration"][mg.core_nodes],
    ...             np.array([0.09648071,  0.44750673,  0.09648071,
    ...                       0.44750673,  1.        ,  0.44750673,
    ...                       0.09648071,  0.44750673,  0.09648071 ]))
    True
    
    
    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    CITATION

    """

    _name = "ConcentrationTrackerSPACE"

    _unit_agnostic = True

    _cite_as = """
    CITATION
    """

    _info = {
        "soil__depth": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Depth of soil or weathered bedrock",
        },
        "sediment__influx": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m^3/yr",
            "mapping": "node",
            "doc": "flux of sediment into node",
        },
        "sediment__outflux": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m^3/yr",
            "mapping": "node",
            "doc": "flux of sediment out of node",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "sediment_property__concentration": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-/m^3",
            "mapping": "node",
            "doc": "Mass concentration of property per unit volume of sediment",
        },
        "bedrock_property__concentration": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-/m^3",
            "mapping": "node",
            "doc": "Mass concentration of property per unit volume of bedrock",
        },
        "sediment_property_production__rate": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-/m^3/yr",
            "mapping": "node",
            "doc": "Production rate of property per unit volume of sediment per time",
        },
        "sediment_property_decay__rate": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-/m^3/yr",
            "mapping": "node",
            "doc": "Decay rate of property per unit volume of sediment per time",
        },
    }

    def __init__(self, 
                 grid,
                 space_instance,
                 concentration_initial=0, 
                 concentration_in_bedrock=0, 
                 local_production_rate=0,
                 local_decay_rate=0
                 ):
        """
        Parameters
        ----------
        grid: ModelGrid
            Landlab ModelGrid object
        concentration_initial: positive float, array, or field name (optional)
            Initial concentration in soil/sediment, -/m^3
        concentration_in_bedrock: positive float, array, or field name (optional)
            Concentration in bedrock, -/m^3
        local_production_rate: float, array, or field name (optional)
            Rate of local production, -/m^3/yr
        local_decay_rate: float, array, or field name (optional)
            Rate of local decay, -/m^3/yr
        """
        
        super().__init__(grid)
        # Store grid and parameters
        
        self._sp = space_instance

        # use setters for C_init, C_br, P, and D defined below
        self.C_init = concentration_initial
        self.C_br = concentration_in_bedrock
        self.P = local_production_rate
        self.D = local_decay_rate
        
        # get reference to inputs
        self._soil__depth = self._grid.at_node["soil__depth"]
        self._soil__depth_old = self._soil__depth.copy()
        self._Qs_out = self._grid.at_node["sediment__outflux"]
        
        # Define variables used for internal calculations
        self._cell_area = self._grid.dx*self._grid.dy
        self._C_sw = np.zeros(self._grid.number_of_nodes)
        self._QsCsw_in = np.zeros(self._grid.number_of_nodes)
        self._QsCsw_out = np.zeros(self._grid.number_of_nodes)
        self._BED_ero_depo_term = np.zeros(self._grid.number_of_nodes)
        
        # create outputs if necessary and get reference.
        self.initialize_output_fields()
        
        # Define concentration field (if all zeros, then add C_init)
        if not self._grid.at_node["sediment_property__concentration"].any():
            self._grid.at_node["sediment_property__concentration"] += self.C_init
        self._concentration = self._grid.at_node["sediment_property__concentration"]

        if not self._grid.at_node["bedrock_property__concentration"].any():
            self._grid.at_node["bedrock_property__concentration"] += self.C_br
        self.C_br = self._grid.at_node["bedrock_property__concentration"]
        
        if not self._grid.at_node["sediment_property_production__rate"].any():
            self._grid.at_node["sediment_property_production__rate"] += self.P
        self.P = self._grid.at_node["sediment_property_production__rate"]
        
        if not self._grid.at_node["sediment_property_decay__rate"].any():
            self._grid.at_node["sediment_property_decay__rate"] += self.D
        self.D = self._grid.at_node["sediment_property_decay__rate"]
        
        # Check that concentration values are within physical limits
        if isinstance(concentration_initial, np.ndarray):
            if concentration_initial.any() < 0:
                raise ValueError("Concentration cannot be negative.")
        else:
            if concentration_initial < 0:
                raise ValueError("Concentration cannot be negative.")

        if isinstance(concentration_in_bedrock, np.ndarray):
            if concentration_in_bedrock.any() < 0:
                raise ValueError("Concentration in bedrock cannot be negative.")
        else:
            if concentration_in_bedrock < 0:
                raise ValueError("Concentration in bedrock cannot be negative.")
        
    @property
    def C_init(self):
        """Initial concentration in soil/sediment (kg/m^3)."""
        return self._C_init
    
    @property
    def C_br(self):
        """Concentration in bedrock (kg/m^3)."""
        return self._C_br
    
    @property
    def P(self):
        """Rate of local production (kg/m^3/yr)."""
        return self._P
    
    @property
    def D(self):
        """Rate of local decay (kg/m^3/yr)."""
        return self._D

    @C_init.setter
    def C_init(self, new_val):
        self._C_init = return_array_at_node(self._grid, new_val)
        
    @C_br.setter
    def C_br(self, new_val):
        self._C_br = return_array_at_node(self._grid, new_val)
        
    @P.setter
    def P(self, new_val):
        self._P = return_array_at_node(self._grid, new_val)
        
    @D.setter
    def D(self, new_val):
        self._D = return_array_at_node(self._grid, new_val)


    def concentration_watercolumn_and_bed(self, dt):
        """Calculate change in concentration within sediment transported in
        the water column and within sediment on the bed for a time period 'dt'.

        Parameters
        ----------

        dt: float (time)
            The imposed timestep.
        """
        # Define values generated by SPACE/SpaceLargeScaleEroder
        flow_receivers = self._grid.at_node["flow__receiver_node"]
        q = self._grid.at_node["surface_water__discharge"]
        phi = self._sp._phi
        F_f = self._sp._F_f
        v_s = self._sp._v_s
        Es = self._sp.Es
        Er = self._sp.Er
        D_sw = v_s*self._Qs_out/q

        # Calculate WC mass balance terms that don't need downstream iteration
        WC_Es_term = (1-phi)*Es*self._cell_area
        WC_Er_term = (1-F_f)*Er*self._cell_area
        WC_denominator_term = 1 + v_s*self._cell_area/q
        
        # Calculate BED mass balance terms that don't need downstream iteration
        with np.errstate(divide='ignore', invalid='ignore'):
            BED_C_local_term = self._concentration * (self._soil__depth_old/self._soil__depth)
            BED_Production_term = (dt*self._P/2) * (self._soil__depth_old/self._soil__depth + 1)
            BED_Decay_term = (dt*self._D/2) * (self._soil__depth_old/self._soil__depth + 1)
        
        np.nan_to_num(BED_C_local_term[self._soil__depth==0])
        np.nan_to_num(BED_Production_term[self._soil__depth==0])
        np.nan_to_num(BED_Decay_term[self._soil__depth==0])
                
        # Get stack of node ids from top to bottom of channel network
        node_status = self._grid.status_at_node
        stack_flip_ud = np.flipud(self._grid.at_node["flow__upstream_node_order"])
        # Select core nodes where qs >0
        stack_flip_ud_sel = stack_flip_ud[
            (node_status[stack_flip_ud] == NodeStatus.CORE)
            & (q[stack_flip_ud] > 0.0)
            ]
        
        # zero out array values that were updated in the old stack
        self._C_sw[:] = 0
        self._QsCsw_in[:] = 0
        self._BED_ero_depo_term[:] = 0

        # Iterate concentration calc (first BED, then WC) at each node
        for node_id in stack_flip_ud_sel:                
            # Calculate QsCsw_out (i.e., QsCs in the water column)
            self._QsCsw_out[node_id] = (
                (self._QsCsw_in[node_id] 
                 + self._concentration[node_id]*WC_Es_term[node_id]
                 + self.C_br[node_id]*WC_Er_term[node_id]
                 )
                / WC_denominator_term[node_id]
                )
            
            # Send QsCsw_out values to flow receiver nodes
            self._QsCsw_in[flow_receivers[node_id]] += self._QsCsw_out[node_id]

            # Divide QsCsw_out (from above) by Qs_out to get C_sw
            if self._Qs_out[node_id] > 0:
                self._C_sw[node_id] = self._QsCsw_out[node_id] / self._Qs_out[node_id]
            else:
                self._C_sw[node_id] = 0.0

            # Calculate BED erosion/deposition term (requires C_sw from above)
            self._BED_ero_depo_term[node_id] = (
                self._C_sw[node_id] * D_sw[node_id]
                - self._concentration[node_id] * (1-phi)*Es[node_id]
                )
            
            # Calculate BED concentration
            with np.errstate(divide='ignore', invalid='ignore'):
                self._concentration[node_id] = (BED_C_local_term[node_id]
                                                + (dt/self._soil__depth[node_id])
                                                * self._BED_ero_depo_term[node_id]
                                                + BED_Production_term[node_id]
                                                - BED_Decay_term[node_id]
                                                )
            np.nan_to_num(self._concentration, copy=False)

        # Update old soil depth to new value
        self._soil__depth_old = self._soil__depth.copy()
        
        
    def run_one_step(self, dt):
        """

        Parameters
        ----------
        dt: float (time)
            The imposed timestep.
        """

        self.concentration_watercolumn_and_bed(dt)
