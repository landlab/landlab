"""
Created on Thu Jun 30 12:40:39 2016

@author: gtucker
"""

import numpy as np
cimport numpy as np
cimport cython


def update_node_states_cython():
    
	# Remember the previous state of each node so we can detect whether the
	# state has changed
	old_tail_node_state = self.node_state[tail_node]
	old_head_node_state = self.node_state[head_node]

	# Change to the new states
	if self.grid.status_at_node[tail_node] == _CORE:
		self.node_state[tail_node] = self.node_pair[new_link_state][0]
	# landlab.grid.base.CORE_NODE:
	if self.grid.status_at_node[head_node] == _CORE:
		self.node_state[head_node] = self.node_pair[new_link_state][1]

	if _DEBUG:
		print('update_node_states() for', tail_node, 'and', head_node)
		print('  tail_node was', old_tail_node_state,
			  'and is now', self.node_state[tail_node])
		print('  head_node was', old_head_node_state,
			  'and is now', self.node_state[head_node])

	return self.node_state[tail_node] != old_tail_node_state, \
		self.node_state[head_node] != old_head_node_state
