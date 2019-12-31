# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 13:08:22 2019

@author: jczuba
"""

# HELP HELP -- link 4 drains to link 2 and it should connect with link 3

# need to determine a specific link on the network
link_number = 4

# def plot_pathway_values(nst_instance)

downstream_link_id = fd.link_to_flow_receiving_node[
    fd.downstream_node_at_link()[link_number]
]

# determine the pathway (set of links) from that link to the outlet

# NEXT NEXT -- collect downstream_link_ids cycling from that link through the
# network to the outlet

# create a 2-d line plot of some specified attribute (y-axis) along that
# pathway (distance downstream; x-axis)

# NEXT NEXT -- for the pathways, index and plot a link attribute value
