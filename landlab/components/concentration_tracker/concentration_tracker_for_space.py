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
    sediment using a mass balance approach in which concentration :math:`C_s`
    is calculated as:

    .. math::
        
                        ∂C_sH/∂t = C_sw*D_sw + C_s*E_s
        
    where :math:`H` is sediment depth, :math:`C_sw` is concentration in 
    sediment suspended in the water column, :math:`D_sw` is volumetric 
    depositional flux of sediment from the water column per unit bed area, and
    :math:`E_s` is volumetric erosional flux of sediment from the bed per unit 
    bed area.
    
    .. note::
    
        This component requires the sediment__influx and sediment__outflux 
        fields calculated by either the :class:`~.Space` or 
        :class:`~.SpaceLargeScaleEroder` component. This component must be run
        after every :class:`~.Space` or :class:`~.SpaceLargeScaleEroder` step
        and before any other flux component. For hillslope sediment tracking,
        see :class:`~.ConcentrationTrackerForDiffusion`.
        
        In-situ production and decay of the material property are handled by
        the ConcentrationTrackerProductionDecay component.

    Examples
    --------
    
    A 1-D stream channel:
        
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import PriorityFloodFlowRouter
    >>> from landlab.components import SpaceLargeScaleEroder
    >>> from landlab.components import ConcentrationTrackerForSpace
    
    >>> mg = RasterModelGrid((3, 5),xy_spacing=10.)
    
    >>> mg.set_status_at_node_on_edges(
    ...     right=NodeStatus.CLOSED,
    ...     top=NodeStatus.CLOSED,
    ...     left=NodeStatus.CLOSED,
    ...     bottom=NodeStatus.CLOSED,
    ... )
    >>> mg.status_at_node[5] = mg.BC_NODE_IS_FIXED_VALUE
    >>> mg.status_at_node.reshape(mg.shape)
    array([[4, 4, 4, 4, 4],
           [1, 0, 0, 0, 4],
           [4, 4, 4, 4, 4]], dtype=uint8)
    
    >>> mg.at_node["sediment_property__concentration"] = [
    ...     [0.0, 0.0, 0.0, 0.0, 0.0],
    ...     [0.0, 0.0, 0.0, 1.0, 0.0],
    ...     [0.0, 0.0, 0.0, 0.0, 0.0],
    ... ]
    >>> mg.at_node["soil__depth"] =  [
    ...     [1.0, 1.0, 1.0, 1.0, 1.0],
    ...     [1.0, 1.0, 1.0, 1.0, 1.0],
    ...     [1.0, 1.0, 1.0, 1.0, 1.0],
    ... ]
    >>> mg.at_node["bedrock__elevation"] = mg.node_x / 100
    >>> mg.at_node["topographic__elevation"] = (
    ...     mg.at_node["soil__depth"] + mg.at_node["bedrock__elevation"]
    ... )
    
    >>> fr = PriorityFloodFlowRouter(mg)
    >>> fr.run_one_step()
    >>> sp = SpaceLargeScaleEroder(mg,phi=0)
    >>> ct = ConcentrationTrackerForSpace(mg,sp)
    
    >>> for i in range(40):
    >>>     fr.run_one_step()
    >>>     sp.run_one_step(10.)
    >>>     ct.run_one_step(10.)

    Erosion has lowered the topography and reduced channel bed sediment depth.
    >>> np.allclose(mg.at_node["topographic__elevation"][mg.core_nodes],
    ...             np.array([1.00292211, 1.00902572, 1.0258774]))
    True
    >>> np.allclose(mg.at_node["soil__depth"][mg.core_nodes],
    ...             np.array([0.90294696, 0.80909071, 0.72601329]))
    True
    
    Some high-concentration sediment has been transported from upstream to be
    deposited on the channel bed further downstream.
    >>> np.allclose(mg.at_node["sediment_property__concentration"][mg.core_nodes],
    ...             np.array([0.04965464, 0.09972316, 0.99991512]))
    True
    
        
    Now, a 2-D landscape with stream channels. All boundaries are closed except
    for Node 0, which is the outlet of the catchment. 

    >>> mg = RasterModelGrid((6, 6),xy_spacing=10.)
    
    >>> mg.set_status_at_node_on_edges(
    ...     right=NodeStatus.CLOSED,
    ...     top=NodeStatus.CLOSED,
    ...     left=NodeStatus.CLOSED,
    ...     bottom=NodeStatus.CLOSED,
    ... )
    >>> mg.status_at_node[0] = mg.BC_NODE_IS_FIXED_VALUE
    >>> mg.status_at_node.reshape(mg.shape)
    array([[4, 4, 4, 4, 4, 4],
           [4, 0, 0, 0, 0, 4],
           [4, 0, 0, 0, 0, 4],
           [4, 0, 0, 0, 0, 4],
           [4, 0, 0, 0, 0, 4],
           [1, 4, 4, 4, 4, 4]], dtype=uint8)
    

    >>> mg.at_node["sediment_property__concentration"] = [
    ...     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ...     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ...     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ...     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ...     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ...     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ... ]
    >>> mg.at_node["soil__depth"] =  [
    ...     [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    ...     [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    ...     [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    ...     [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    ...     [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    ...     [0.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    ... ]
        
    # Add slope and noise to the bedrock to create some topographic structure.
    >>> np.random.seed(5)
    >>> mg.at_node["bedrock__elevation"] = (mg.node_x + mg.node_y) / 1000
    >>> mg.at_node["bedrock__elevation"] += np.random.rand(mg.number_of_nodes) / 100
    >>> mg.at_node["bedrock__elevation"][0] = 0
    
    >>> mg.at_node["topographic__elevation"] = (
    ...     mg.at_node["soil__depth"] + mg.at_node["bedrock__elevation"]
    ...     )
    >>> imshow_grid(mg,"topographic__elevation", cmap=mpl.cm.get_cmap("terrain").copy())
    
    # Instantiate components.
    >>> fr = PriorityFloodFlowRouter(mg)
    >>> fr.run_one_step()
    >>> sp = SpaceLargeScaleEroder(mg,phi=0)
    >>> ct = ConcentrationTrackerForSpace(mg,sp)
    
    >>> imshow_grid(mg,"drainage_area")

    # Run SPACE for 1,000 years to generate a fluvial network.
    >>> for i in range(1000):
    >>>     mg.at_node["bedrock__elevation"][mg.core_nodes] += 0.001
    >>>     mg.at_node["topographic__elevation"] = (
    ...         mg.at_node["soil__depth"] + mg.at_node["bedrock__elevation"]
    ...     )
    >>>     fr.run_one_step()
    >>>     sp.run_one_step(1.)
    >>>     if np.mod(i,100) == 0:
    >>>         mpl.pyplot.figure()
    >>>         imshow_grid(mg,"topographic__elevation", cmap=mpl.cm.get_cmap("terrain").copy())
    
    >>> mpl.pyplot.figure()
    >>> imshow_grid(mg,"drainage_area")

    # Set high concentration at a headwater node to trace sediment downstream.
    >>> mg.at_node["sediment_property__concentration"][22] += 1    
    
    >>> c_22 = np.zeros(100)
    
    >>> for i in range(100):
    >>>     c_22[i] = mg.at_node["sediment_property__concentration"][22]
    >>>     mg.at_node["bedrock__elevation"][mg.core_nodes] += 0.001
    >>>     mg.at_node["topographic__elevation"] = (
    ...         mg.at_node["soil__depth"] + mg.at_node["bedrock__elevation"]
    ...     )
    >>>     fr.run_one_step()
    >>>     sp.run_one_step(1.)
    >>>     ct.run_one_step(1.)
    
    >>> mpl.pyplot.figure()
    >>> mpl.pyplot.plot(c_22)

    WHY DOES CONCENTRATION SKYROCKET ON THE FIRST TIMESTEP???
    WHY DOES CONCENTRATION INCREASE ABOVE THE ORIGINAL INPUT VALUE????
    IT MUST BE SOMETHING TO DO WITH NOT HAVING A PREVIOUS TIMESTEP VALUE FOR
    ONE OF EROSION, DEPOSITION, Qs, ETC......

    >>> mpl.pyplot.figure()
    >>> imshow_grid(mg,"soil__depth")
    
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
    
    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    CITATION

    """

    _name = "ConcentrationTrackerForSpace"

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
    }

    def __init__(self, 
                 grid,
                 space_instance, # NO LONGER NEED TO INPUT SPACE INSTANCE SINCE THE q VALUES ARE NO LONGER ZEROED
                 concentration_initial=0, 
                 concentration_in_bedrock=0, 
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
        """
        
        super().__init__(grid)
        # Store grid and parameters
        
        self._sp = space_instance # INSTEAD OF SPACE INSTANCE, USE VALUES NOW EXPORTED FROM SPACE...

        # use setters for C_init, C_br defined below
        self.C_init = concentration_initial
        self.C_br = concentration_in_bedrock
        
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
        if np.allclose(self._grid.at_node["sediment_property__concentration"], 0.0):
            self._grid.at_node["sediment_property__concentration"] += self.C_init
        self._concentration = self._grid.at_node["sediment_property__concentration"]

        if np.allclose(self._grid.at_node["bedrock_property__concentration"], 0.0):
            self._grid.at_node["bedrock_property__concentration"] += self.C_br
        self.C_br = self._grid.at_node["bedrock_property__concentration"]
                
    @property
    def C_init(self):
        """Initial concentration in soil/sediment (kg/m^3)."""
        return self._C_init
    
    @property
    def C_br(self):
        """Concentration in bedrock (kg/m^3)."""
        return self._C_br

    @C_init.setter
    def C_init(self, new_val):
        if np.any(new_val < 0.0):
            raise ValueError("Concentration in sediment cannot be negative")
        self._C_init = return_array_at_node(self._grid, new_val)
        
    @C_br.setter
    def C_br(self, new_val):
        if np.any(new_val < 0.0):
            raise ValueError("Concentration in bedrock cannot be negative")
        self._C_br = return_array_at_node(self._grid, new_val)


    def calc_concentration_watercolumn_and_bed(self, dt):
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
        phi = self._sp._phi   # USE VALUES INSTEAD OF SPACE INSTANCE
        F_f = self._sp._F_f   # USE VALUES INSTEAD OF SPACE INSTANCE
        v_s = self._sp._v_s   # USE VALUES INSTEAD OF SPACE INSTANCE
        Es = self._sp.Es   # USE VALUES INSTEAD OF SPACE INSTANCE
        Er = self._sp.Er   # USE VALUES INSTEAD OF SPACE INSTANCE
        D_sw = np.zeros(np.shape(q))
        D_sw[q!=0] = v_s*self._Qs_out[q!=0]/q[q!=0]
        
        # Back-calculate soil depth from prior to running SPACE component
        #(REPLACE WITH AN OLD SOIL DEPTH FROM SPACE COMPONENT IF IT HAS ONE)
        self._soil__depth_old = self._soil__depth.copy() + (Es - D_sw)/(1-phi)
        
        # Calculate portions of equation that have soil depth as denominator
        is_soil = self._soil__depth > 0.0
        
        old_depth_over_new = np.divide(self._soil__depth_old, self._soil__depth, where=is_soil)
        old_depth_over_new[~is_soil] = 0.0
        
        dt_over_depth = np.divide(dt, self._soil__depth, where=is_soil)
        dt_over_depth[~is_soil] = 0.0
        
        # Calculate mass balance terms that don't need downstream iteration
        WC_Es_term = Es*self._cell_area/(1-phi)
        WC_Er_term = (1-F_f)*Er*self._cell_area
        WC_denominator_term = np.ones(np.shape(q))
        WC_denominator_term[q!=0] = 1 + v_s*self._cell_area/q[q!=0]        
        BED_C_local_term = self._concentration * old_depth_over_new
                        
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
                self._C_sw[node_id] * D_sw[node_id]/(1-phi)
                - self._concentration[node_id] * Es[node_id]/(1-phi)
                )
            
            # Calculate BED concentration
            self._concentration[node_id] = (BED_C_local_term[node_id]
                                            + dt_over_depth[node_id]
                                            * self._BED_ero_depo_term[node_id]
                                            )
                
            self._concentration[~is_soil] = 0.0      
        
    def run_one_step(self, dt):
        """

        Parameters
        ----------
        dt: float (time)
            The imposed timestep.
        """

        self.calc_concentration_watercolumn_and_bed(dt)
