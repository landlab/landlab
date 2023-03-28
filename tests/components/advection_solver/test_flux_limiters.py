from landlab.components.advection.flux_limiters import flux_lim_vanleer


def test_flux_lim_vanleer():
    """Note: same tests apply to all flux limiters."""
    assert flux_lim_vanleer(0.0) == 0.0
    assert flux_lim_vanleer(0.5) >= 0.5
    assert flux_lim_vanleer(0.5) <= 1.0
    assert flux_lim_vanleer(1.0) == 1.0
    assert flux_lim_vanleer(1.5) >= 1.0
    assert flux_lim_vanleer(1.5) <= 1.5
    assert flux_lim_vanleer(2.0) >= 1.0
    assert flux_lim_vanleer(2.0) <= 2.0
    assert flux_lim_vanleer(10.0) >= 1.0
    assert flux_lim_vanleer(10.0) <= 2.0
