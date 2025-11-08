from landlab.components.flow_router.ext.single_flow.accumulation import (
    _test_accumulation_c as accumulation,
)


def test_add_to_upstream_ordered_nodes():
    accumulation.test_add_to_upstream_ordered_nodes()


def test_calc_upstream_order_for_nodes_c():
    accumulation.test_calc_upstream_order_for_nodes_c()


def test_calc_upstream_order_for_nodes():
    accumulation.test_calc_upstream_order_for_nodes()


def test_calc_drainage_areas():
    accumulation.test_calc_drainage_areas()
