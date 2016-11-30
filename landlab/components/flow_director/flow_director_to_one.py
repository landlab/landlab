from landlab import FieldError
from landlab.components.flow_director.flow_director import FlowDirector
import numpy

class FlowDirectorToOne(FlowDirector):
    """
    Intermediate level class for calculating flow directions. 
    
    This component is not meant to be used directly in modeling efforts. 
    Instead it has the functionality that all flow direction calculators need
    if they route flow only to one cell, as in D4 or D8 direction finding. 
    It mostly creates model grid fields. 
    
    The primary method of this class, :func:`run_one_step` is not implemented.


    Parameters
    ----------
    grid : ModelGrid
        A grid.
    surface : field name at node or array of length node
        The surface to direct flow across.   
        
    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.flow_director.flow_director_to_one import FlowDirectorToOne
    >>> mg = RasterModelGrid((3,3), spacing=(1, 1))
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> _ = mg.add_field('topographic__elevation', mg.node_x + mg.node_y, at = 'node')
    >>> fd=FlowDirectorToOne(mg, 'topographic__elevation')
    >>> fd.elevs
    array([ 0.,  1.,  2.,  1.,  2.,  3.,  2.,  3.,  4.])
    >>> keys=set(['flow__link_to_receiver_node', 'flow__sink_flag',
    ...           'topographic__steepest_slope', 'topographic__elevation', 
    ...           'flow__receiver_node'])
    >>> len(mg.at_node.keys()-keys)
    0

    """
    
    _name = 'FlowDirectorToOne'

    _input_var_names = ('topographic__elevation',
                        )

    _output_var_names = ('flow__receiver_node',
                         'topographic__steepest_slope',
                         'flow__upstream_node_order',
                         'flow__sink_flag',
                         )

    _var_units = {'topographic__elevation': 'm',
                  'flow__receiver_node': '-',
                  'topographic__steepest_slope': '-',
                  'flow__upstream_node_order': '-',
                  'flow__link_to_receiver_node': '-',
                  'flow__sink_flag': '-',
                  }

    _var_mapping = {'topographic__elevation': 'node',
                    'flow__receiver_node': 'node',
                    'topographic__steepest_slope': 'node',
                    'flow__link_to_receiver_node': 'node',
                    'flow__sink_flag': 'node',
                    }

    _var_doc = {
        'topographic__elevation': 'Land surface topographic elevation',
        'flow__receiver_node':
            'Node array of receivers (node that receives flow from current '
            'node)',
        'topographic__steepest_slope':
            'Node array of steepest *downhill* slopes',
        'surface_water__discharge': 'Discharge of water through each node',
        'flow__link_to_receiver_node':
            'ID of link downstream of each node, which carries the discharge',
        'flow__sink_flag': 'Boolean array, True at local lows',
    }
    

    def __init__(self, grid, surface):
        # run init for the inherited class
        super(FlowDirectorToOne, self).__init__(grid, surface)
        
        # initialize new fields
        try:
            self.receiver = grid.add_zeros('flow__receiver_node', at='node',
                                           dtype=int)
        except FieldError:
            self.receiver = grid.at_node['flow__receiver_node']
        
        try:
            self.steepest_slope = grid.add_zeros(
                'topographic__steepest_slope', at='node', dtype=float)
        except FieldError:
            self.steepest_slope = grid.at_node['topographic__steepest_slope']
        
        try:
            self.links_to_receiver = grid.add_zeros(
                'flow__link_to_receiver_node', at='node', dtype=int)
        except FieldError:
            self.links_to_receiver = grid.at_node[
                'flow__link_to_receiver_node']
        
        grid.add_zeros('flow__sink_flag', at='node', dtype=numpy.int8,
                       noclobber=False)
       
    
    def run_one_step(self):
        raise NotImplementedError('run_one_step()')
    
    
    # set properties. These are the same for all DirectToOne Directors
    @property
    def node_receiving_flow(self):
        return self._grid['node']['flow__receiver_node']

    @property
    def node_steepest_slope(self):
        return self._grid['node']['topographic__steepest_slope']

    @property
    def link_to_flow_receiving_node(self):
        return self._grid['node']['flow__link_to_receiver_node']    
if __name__ == '__main__':
    import doctest
    doctest.testmod()