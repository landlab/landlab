"""Plot drainage network.

"""
# KRB, FEB 2017.

def drainage_plot(mg, 
                  surface='topographic__elevation', 
                  receivers=None, 
                  proportions=None, 
                  surf_cmap='gray',
                  quiver_cmap='viridis'):
    
    import matplotlib.pylab as plt
    from landlab.plot.imshow import imshow_node_grid
    import numpy as np
    
    imshow_node_grid(mg, surface, cmap=surf_cmap)
    
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
        
        xdist[is_bad] = np.nan
        ydist[is_bad] = np.nan
          
        if proportions is None:
           proportions = np.ones_like(receivers, dtype=float)
       
        prop = proportions[:,j]*256.
        lu = np.floor(prop)
        colors = propColor(lu.astype(int))
        plt.quiver(mg.node_x.reshape(mg.shape), mg.node_y.reshape(mg.shape), 
               xdist.reshape(mg.shape), ydist.reshape(mg.shape), 
               color=colors,
               angles='xy',
               scale_units='xy', 
               scale=1) 
               
    plt.plot(mg.node_x, mg.node_y, 'r.')
    sm = plt.cm.ScalarMappable(cmap=propColor, norm=plt.Normalize(vmin=0, vmax=1))
    sm._A = []
    cx = plt.colorbar(sm)
    cx.set_label('Proportion of Flow')
    plt.show()