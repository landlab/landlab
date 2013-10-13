#! /usr/bin/env python

def map_values_from_link_tail_node_to_link(mg, var_name):
    values_at_nodes = mg.at_node[var_name]
    mg.add_empty('link', var_name)
    values_at_links = mg.at_link[var_name]
    values_at_links[:] = values_at_nodes[mg.node_index_at_link_tail]


def map_values_from_link_head_node_to_link(mg, var_name):
    values_at_nodes = mg.at_node[var_name]
    mg.add_empty('link', var_name)
    values_at_links = mg.at_link[var_name]
    values_at_links[:] = values_at_nodes[mg.node_index_at_link_head]


def map_values_from_link_end_nodes_to_link(mg, var_name):
    values_at_nodes = mg.at_node[var_name]
    mg.add_empty('link', var_name)
    values_at_links = mg.at_link[var_name]
    values_at_links[:] = 0.5 * (values_at_nodes[mg.node_index_at_link_head] +
                                values_at_nodes[mg.node_index_at_link_tail])

def map_values_from_cell_node_to_cell(mg, var_name):
    values_at_nodes = mg.at_node[var_name]
    mg.add_empty('cell', var_name)
    values_at_cells = mg.at_cell[var_name]
    values_at_cells[:] = values_at_nodes[mg.node_index_at_cells]
