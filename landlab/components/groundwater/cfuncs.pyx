import numpy as np

cimport numpy as np
cimport cython


DTYPE = np.int
ctypedef np.int_t DTYPE_INT_t
DTYPE_FLOAT = np.double
ctypedef np.double_t DTYPE_FLOAT_t

# regularization functions used to deal with numerical demons of seepage
cdef _regularize_G(np.ndarray[DTYPE_FLOAT_t, ndim=1] u, DTYPE_FLOAT_t reg_factor):
	"""Smooths transition of step function with an exponential. 0<=u<=1."""
	return np.exp(-(1 - u) / reg_factor)


cdef _regularize_R(np.ndarray[DTYPE_FLOAT_t, ndim=1] u):
	"""ramp function on u."""
	return u * np.greater_equal(u, 0.)

cdef _update_thickness(DTYPE_FLOAT_t dt,
		np.ndarray[DTYPE_FLOAT_t, ndim=1] h0,
		np.ndarray[DTYPE_FLOAT_t, ndim=1] b,
		np.ndarray[DTYPE_FLOAT_t, ndim=1] f,
		np.ndarray[DTYPE_FLOAT_t, ndim=1] dqdx,
		np.ndarray[DTYPE_FLOAT_t, ndim=1] n,
		DTYPE_FLOAT_t r,
		):
	"""analytical solution for the linearized governing equation."""
	cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] out
	cdef np.ndarray[np.uint8_t, ndim=1] cond = f <= dqdx
	out = b * (
		1
		- r
		* np.log(
			1
			+ np.exp((1 - (h0 + ((f - dqdx) * dt) / n) / b) / r)
			* (1 - np.exp(-(b - h0) / (b * r)))
		)
	)
	out[cond] = (h0 + (1 / n * (f - dqdx)) * dt)[cond]
	return out


def _calc_grad_at_link(np.ndarray[DTYPE_FLOAT_t, ndim=1] values,
		np.ndarray[DTYPE_INT_t, ndim=1] node_head,
		np.ndarray[DTYPE_INT_t, ndim=1] node_tail,
		np.ndarray[DTYPE_FLOAT_t, ndim=1] link_lengths,
		):

	cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] out = np.empty(len(node_head), dtype=float)
	return np.divide(
		values[node_head] - values[node_tail],
		link_lengths,
		out=out,
	)


def _calc_flux_div_at_node(np.ndarray[DTYPE_FLOAT_t, ndim=1] unit_flux,
		np.ndarray[DTYPE_INT_t, ndim=1] node_at_cell,
		np.ndarray[DTYPE_INT_t, ndim=1] link_at_face,
		np.ndarray[DTYPE_INT_t, ndim=2] faces_at_cell,
		np.ndarray[DTYPE_FLOAT_t, ndim=1] area_of_cell,
		np.ndarray[DTYPE_FLOAT_t, ndim=1] length_of_face,
		np.ndarray[signed char, ndim=2] link_dirs_at_node,
		):

	cdef int number_of_dirs = link_dirs_at_node.shape[1]
	cdef int number_of_nodes = link_dirs_at_node.shape[0]
	cdef int c
	cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] buf = np.empty_like(area_of_cell)
	cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] out = np.empty(number_of_nodes, dtype=float)
	cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] total_flux = unit_flux[link_at_face] * length_of_face

	for c in range(number_of_dirs):
		buf -= total_flux[faces_at_cell[:, c]] * link_dirs_at_node[node_at_cell, c]
	out[node_at_cell] = buf / area_of_cell
	return out

def _map_value_at_max_node_to_link(np.ndarray[DTYPE_FLOAT_t, ndim=1] controls,
		np.ndarray[DTYPE_FLOAT_t, ndim=1] values,
		np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_head,
		np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_tail,
		):


	cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] head_control = controls[node_at_link_head]
	cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] tail_control = controls[node_at_link_tail]
	cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] head_vals = values[node_at_link_head]
	cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] tail_vals = values[node_at_link_tail]
	cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] out = np.empty(len(node_at_link_head), dtype=float)

	out[:] = np.where(tail_control > head_control, tail_vals, head_vals)
	return out


def run_substeps(DTYPE_FLOAT_t dt,
		np.ndarray[DTYPE_FLOAT_t, ndim=1] wtable,
		np.ndarray[DTYPE_FLOAT_t, ndim=1] base,
		np.ndarray[DTYPE_FLOAT_t, ndim=1] thickness,
		np.ndarray[DTYPE_FLOAT_t, ndim=1] recharge,
		np.ndarray[DTYPE_FLOAT_t, ndim=1] K,
		np.ndarray[DTYPE_FLOAT_t, ndim=1] n,
		np.ndarray[DTYPE_FLOAT_t, ndim=1] n_link,
		np.ndarray[DTYPE_FLOAT_t, ndim=1] cosa,
		np.ndarray[DTYPE_FLOAT_t, ndim=1] reg_thickness,
		np.ndarray[DTYPE_INT_t, ndim=1] cores,
		np.ndarray[DTYPE_INT_t, ndim=1] node_head,
		np.ndarray[DTYPE_INT_t, ndim=1] node_tail,
		np.ndarray[DTYPE_FLOAT_t, ndim=1] link_length,
		np.ndarray[DTYPE_INT_t, ndim=1] active_links,
		np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_head,
		np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_tail,
		np.ndarray[DTYPE_INT_t, ndim=1] node_at_cell,
		np.ndarray[DTYPE_INT_t, ndim=1] link_at_face,
		np.ndarray[DTYPE_INT_t, ndim=2] faces_at_cell,
		np.ndarray[DTYPE_FLOAT_t, ndim=1] area_of_cell,
		np.ndarray[DTYPE_FLOAT_t, ndim=1] length_of_face,
		np.ndarray[signed char, ndim=2] link_dirs_at_node,
		DTYPE_FLOAT_t vn,
		DTYPE_FLOAT_t cour,
		DTYPE_FLOAT_t r,
		):

	cdef DTYPE_FLOAT_t remaining_time = dt
	cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] hydr_grad = np.zeros_like(wtable)
	cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] vel
	cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] hlink
	cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] q
	cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] dqdx
	cdef DTYPE_FLOAT_t dt_vn
	cdef DTYPE_FLOAT_t dt_courant
	cdef DTYPE_FLOAT_t substep_dt
	cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] qs
	cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] qs_cumulative
	cdef DTYPE_INT_t num_substeps = 0

	while remaining_time > 0.0:

		# Calculate hydraulic gradient
		hydr_grad[active_links] = (
			_calc_grad_at_link(wtable, node_head, node_tail, link_length) * cosa
		)[active_links]

		# Calculate groundwater velocity
		vel = -K * hydr_grad

		# Aquifer thickness at links (upwind)
		hlink = _map_value_at_max_node_to_link(wtable, thickness, node_at_link_head, node_at_link_tail)	* cosa


		# Calculate specific discharge
		q = hlink * vel

		# Groundwater flux divergence
		dqdx = _calc_flux_div_at_node(q,
				node_at_cell,
				link_at_face,
				faces_at_cell,
				area_of_cell,
				length_of_face,
				link_dirs_at_node,
		)

		# calculate criteria for timestep
		dt_vn = vn * np.min(
			np.divide(
				(n_link * link_length ** 2),
				(4 * K * hlink),
				where=hlink > 0.0,
				out=np.ones_like(q) * 1e15,
			)
		)

		dt_courant = cour * np.min(
			np.divide(
				link_length,
				abs(vel / n_link),
				where=vel > 0.0,
				out=np.ones_like(q) * 1e15,
			)
		)
		substep_dt = min([dt_courant, dt_vn, remaining_time])
		# print(np.argmin(np.array([self._dt_courant, self._dt_vn, remaining_time]))) # 0 = courant limited, 1 = vn limited, 2 = not limited

		# update thickness from analytical
		thickness[cores] = _update_thickness(
			substep_dt,
			thickness,
			reg_thickness,
			recharge,
			dqdx,
			n,
			r,
		)[cores]
		thickness[thickness < 0.0] = 0.0

		# Recalculate water surface height
		wtable = base + thickness

		# Calculate surface discharge at nodes
		qs = _regularize_G(
			thickness / reg_thickness, r
		) * _regularize_R(recharge - dqdx)

		# add cumulative sw discharge in substeps
		qs_cumulative += qs * substep_dt

		# calculate the time remaining and advance count of substeps
		remaining_time -= substep_dt
		num_substeps += 1

	return wtable, thickness, q, qs, qs_cumulative, num_substeps
