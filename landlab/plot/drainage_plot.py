"""Plot drainage network.

"""
# KRB, FEB 2017.
import six
from landlab import CORE_NODE, FIXED_VALUE_BOUNDARY, FIXED_GRADIENT_BOUNDARY, CLOSED_BOUNDARY
import matplotlib.pylab as plt
from landlab.plot.imshow import imshow_node_grid
import numpy as np

def drainage_plot(mg, 
                  surface='topographic__elevation', 
                  receivers=None, 
                  proportions=None, 
                  surf_cmap='gray',
                  quiver_cmap='viridis',
                  title = 'Drainage Plot'):
    

    if isinstance(surface, six.string_types):
        colorbar_label = surface
    else:
        colorbar_label = 'topographic_elevation'
    imshow_node_grid(mg, surface, cmap=surf_cmap, colorbar_label=colorbar_label)
    
    if receivers is None:
        try:
            receivers = mg.at_node['flow__receiver_nodes']
            if proportions is None:
                try:
                    proportions = mg.at_node['flow__receiver_proportions']
                except:
                    pass
        except: 
            receivers = np.reshape(mg.at_node['flow__receiver_node'],(mg.number_of_nodes,1))
        
    nreceievers = int(receivers.size/receivers.shape[0])
    
    propColor=plt.get_cmap(quiver_cmap)

    for j in range(nreceievers):
        rec = receivers[:,j]
        is_bad = rec == -1
        
        xdist =  -0.8*(mg.node_x-mg.node_x[rec])
        ydist =  -0.8*(mg.node_y-mg.node_y[rec])
        
        if proportions is None:
           proportions = np.ones_like(receivers, dtype=float)
        
        is_bad[proportions[:,j]==0.]=True
        
        xdist[is_bad] = np.nan
        ydist[is_bad] = np.nan
             
        prop = proportions[:,j]*256.
        lu = np.floor(prop)
        colors = propColor(lu.astype(int))
        
        shape = (mg.number_of_nodes, 1)
        
        plt.quiver(mg.node_x.reshape(shape), mg.node_y.reshape(shape), 
               xdist.reshape(shape), ydist.reshape(shape), 
               color=colors,
               angles='xy',
               scale_units='xy', 
               scale=1,
               zorder=3) 
    
    # Plot differen types of nodes:
    o, = plt.plot(mg.node_x[mg.status_at_node == CORE_NODE], mg.node_y[mg.status_at_node == CORE_NODE], 'b.', label='Core Nodes', zorder=4)
    fg, = plt.plot(mg.node_x[mg.status_at_node == FIXED_VALUE_BOUNDARY], mg.node_y[mg.status_at_node == FIXED_VALUE_BOUNDARY], 'c.', label='Fixed Gradient Nodes', zorder=5)
    fv, = plt.plot(mg.node_x[mg.status_at_node == FIXED_GRADIENT_BOUNDARY], mg.node_y[mg.status_at_node ==FIXED_GRADIENT_BOUNDARY], 'g.', label='Fixed Value Nodes', zorder=6)
    c, = plt.plot(mg.node_x[mg.status_at_node == CLOSED_BOUNDARY], mg.node_y[mg.status_at_node ==CLOSED_BOUNDARY], 'r.', label='Closed Nodes', zorder=7)

    node_id = np.arange(mg.number_of_nodes)
    flow_to_self = receivers[:,0]==node_id
                            
    fts, = plt.plot(mg.node_x[flow_to_self], mg.node_y[flow_to_self], 'kx', markersize=6, label = 'Flows To Self', zorder=8)
    
    ax = plt.gca()
    
    ax.legend(labels = ['Core Nodes', 'Fixed Gradient Nodes', 'Fixed Value Nodes', 'Closed Nodes', 'Flows To Self'],
              handles = [o, fg, fv, c, fts], numpoints=1, loc='center left', bbox_to_anchor=(1.7, 0.5))
    sm = plt.cm.ScalarMappable(cmap=propColor, norm=plt.Normalize(vmin=0, vmax=1))
    sm._A = []
    cx = plt.colorbar(sm)
    cx.set_label('Proportion of Flow')
    plt.title(title)
    plt.show()