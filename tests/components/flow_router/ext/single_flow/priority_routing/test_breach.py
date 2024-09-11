from landlab.components.flow_router.ext.single_flow.priority_routing import (
    _test_breach_c as breach,
)


def test_priority_queue():
    breach.test_priority_queue()


def test_init_flow_direction_queues():
    breach.test_init_flow_direction_queues()


def test_set_flooded_and_outlet():
    breach.test_set_flooded_and_outlet()


def test_set_receiver():
    breach.test_set_receiver()


def test_set_donor_properties():
    breach.test_set_donor_properties()


def test_direct_flow_c():
    breach.test_direct_flow_c()


def test_direct_flow():
    breach.test_direct_flow()
