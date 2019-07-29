"""Test the SedDepEroder component.

Test the sed dep eroder by turning it over a few times. No attempt has been
made to ensure the solution is stable. Takes a topo already output and runs it
a few more times, to ensure repeatability.
"""
import os
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_equal
from six.moves import range
from matplotlib.pyplot import gca, clf
import pytest

from landlab import RasterModelGrid, VoronoiDelaunayGrid
from landlab import CLOSED_BOUNDARY, ModelParameterDictionary, FieldError
from landlab.components import FlowAccumulator
from landlab.components import SedDepEroder
from landlab.components import FastscapeEroder
from landlab.components import DepressionFinderAndRouter

from landlab.components.stream_power.cfuncs import (
    sed_flux_fn_gen_genhump, sed_flux_fn_gen_lindecl,
    sed_flux_fn_gen_almostparabolic, sed_flux_fn_gen_const,
    get_sed_flux_function_pseudoimplicit_bysedout,
    iterate_sde_downstream
)


def test_flux_fn_const():
    """
    Tests that the const function always returns 1.
    """
    fnval = sed_flux_fn_gen_const(0., np.nan, np.nan, np.nan, np.nan, np.nan)
    assert np.isclose(fnval, 1.)
    fnval = sed_flux_fn_gen_const(1., np.nan, np.nan, np.nan, np.nan, np.nan)
    assert np.isclose(fnval, 1.)
    fnval = sed_flux_fn_gen_const(2., np.nan, np.nan, np.nan, np.nan, np.nan)
    assert np.isclose(fnval, 1.)


def test_flux_fn_lindecl():
    """
    Tests that the linear decline function returns correct values.
    """
    fnval = sed_flux_fn_gen_lindecl(0., np.nan, np.nan, np.nan, np.nan, np.nan)
    assert np.isclose(fnval, 1.)
    fnval = sed_flux_fn_gen_lindecl(1., np.nan, np.nan, np.nan, np.nan, np.nan)
    assert np.isclose(fnval, 0.)
    fnval = sed_flux_fn_gen_lindecl(
        0.5, np.nan, np.nan, np.nan, np.nan, np.nan
    )
    assert np.isclose(fnval, 0.5)
    # observe we permit undefined vals
    fnval = sed_flux_fn_gen_lindecl(2., np.nan, np.nan, np.nan, np.nan, np.nan)
    assert np.isclose(fnval, -1.)


def test_flux_fn_almostpara():
    """
    Tests that the almost parabolic function returns correct values.
    """
    fnval = sed_flux_fn_gen_almostparabolic(
        0., np.nan, np.nan, np.nan, np.nan, np.nan)
    assert np.isclose(fnval, 0.1)
    fnval = sed_flux_fn_gen_almostparabolic(
        0.5, np.nan, np.nan, np.nan, np.nan, np.nan)
    assert np.isclose(fnval, 1.)
    fnval = sed_flux_fn_gen_almostparabolic(
        1., np.nan, np.nan, np.nan, np.nan, np.nan)
    assert np.isclose(fnval, 0.)


def test_flux_fn_genhump():
    """
    Tests that the generalized function returns correct values.
    """
    fnval = sed_flux_fn_gen_genhump(
        0., 13.683, 1.13, 0.00181, 4.24, 1.0000278041373)
    # remember, phi & c are the weird way round
    assert np.isclose(fnval, 0.024766918603659326)
    fnval = sed_flux_fn_gen_genhump(
        1., 13.683, 1.13, 0.00181, 4.24, 1.0000278041373)
    assert np.isclose(fnval, 0.1975013921257844)
    max_val = 0.
    peak_at = np.nan
    for i in np.arange(0., 1.001, 0.001):
        sff = sed_flux_fn_gen_genhump(
            i, 13.683, 1.13, 0.00181, 4.24, 1.0000278041373)
        if sff > max_val:
            max_val = max((sff, max_val))
            peak_at = i
    assert np.isclose(max_val, 1.)
    assert np.isclose(peak_at, 0.264)
    # finally, check the behaviour at high rel sed fluxes
    fnval = sed_flux_fn_gen_genhump(
        1.e6, 13.683, 1.13, 0.00181, 4.24, 1.0000278041373)
    assert np.isclose(fnval, 0.)
    fnval = sed_flux_fn_gen_genhump(
        2., 13.683, 1.13, 0.00181, 4.24, 1.0000278041373)
    assert np.isclose(fnval, 0.006221557384902739)
    # ^note that the function continues to evolve at RSF>1


def test_set_sed_flux_fn_gen():
    """
    This tests the setter for the sff.
    """
    mg = RasterModelGrid((5, 5))
    z = mg.add_zeros('node', 'topographic__elevation')
    pop_me = set([name for name in SedDepEroder.output_var_names])
    pop_me.discard('topographic__elevation')

    fr = FlowAccumulator(mg, flow_director='D8')
    sde = SedDepEroder(mg, K_sp=1.e-4, sed_dependency_type='None',
                       Qc='power_law', K_t=1.e-4)
    sde._sed_flux_fn_gen is sed_flux_fn_gen_const
    for pop in pop_me:
        null = mg.at_node.pop(pop)
    sde = SedDepEroder(mg, K_sp=1.e-4, sed_dependency_type='linear_decline',
                       Qc='power_law', K_t=1.e-4)
    sde._sed_flux_fn_gen is sed_flux_fn_gen_lindecl
    for pop in pop_me:
        null = mg.at_node.pop(pop)
    sde = SedDepEroder(mg, K_sp=1.e-4, sed_dependency_type='almost_parabolic',
                       Qc='power_law', K_t=1.e-4)
    sde._sed_flux_fn_gen is sed_flux_fn_gen_almostparabolic
    for pop in pop_me:
        null = mg.at_node.pop(pop)
    sde = SedDepEroder(
        mg, K_sp=1.e-4, sed_dependency_type='generalized_humped',
        Qc='power_law', phi_hump=4., K_t=1.e-4
    )
    sde._sed_flux_fn_gen is sed_flux_fn_gen_genhump
    # ...& add a test that the auto-normalization is working OK...
    max_val_comp = 0.
    peak_at_comp = np.nan
    for i in np.arange(0., 1.001, 0.001):
        sff = sed_flux_fn_gen_genhump(
            i, 13.683, 1.13, 0.00181, 4., 1./1.0674809986373335)
        compval = sde._sed_flux_fn_gen(
            i, sde.kappa, sde.nu, sde.c, sde.phi, sde.norm
        )
        assert np.isclose(sff, compval)
        if sff > max_val_comp:
            max_val_comp = max((sff, max_val_comp))
            peak_at = i
    assert np.isclose(max_val_comp, 1.)
    assert np.isclose(peak_at, 0.28)


def test_sff_convergence():
    """
    This tests the stable convergence of the sffs.

    get_sed_flux_function_pseudoimplicit(
        sed_in_bydt,
        trans_cap_vol_out_bydt,
        prefactor_for_volume_bydt,
        prefactor_for_dz_bydt,
        sed_flux_fn_gen,
        kappa, nu, c, phi, norm,
        pseudoimplicit_repeats,
        out_array
    )
    """
    humpk = 13.683
    humpnu = 1.13
    humpc = 0.00181
    humpphi = 4.24
    humpnorm = 1.0000278041373
    models = ('linear_decline', )
    funcs = (sed_flux_fn_gen_lindecl, )
    for (sff_style, solver) in zip(models, funcs):
        mg = RasterModelGrid((5, 5))
        z = mg.add_zeros('node', 'topographic__elevation')
        fa = FlowAccumulator(mg)
        out_array = np.empty(4, dtype=float)
        sde = SedDepEroder(mg, sed_dependency_type=sff_style)
        # special case
        get_sed_flux_function_pseudoimplicit_bysedout(1000., 0., 1., 1.,
                                                      sde._sed_flux_fn_gen,
                                                      humpk, humpnu,
                                                      humpc, humpphi,
                                                      humpnorm,
                                                      50, out_array)
        assert np.isclose(out_array[0], 0.)  # dzbydt
        assert np.isclose(out_array[1], 0.)  # vol_pass_rate
        assert np.isclose(out_array[2], 1.)  # rel_sed_flux
        assert np.isclose(out_array[3], 0.)  # error_in_sed_flux_fn

        # very high flux in. Cell area is 10.
        get_sed_flux_function_pseudoimplicit_bysedout(1.e10, 1., 10., 10.,
                                                      sde._sed_flux_fn_gen,
                                                      humpk, humpnu,
                                                      humpc, humpphi,
                                                      humpnorm,
                                                      50, out_array)
        known_fqs = solver(
            1., humpk, humpnu, humpc, humpphi, humpnorm,
        )
        # Note that this test demonstrates that the range of relative sediment
        # flux is limited to 0. <= qs/qc <= 1.
        assert np.isclose(out_array[0], known_fqs)  # dzbydt
        assert np.isclose(out_array[1], 1.)  # vol_pass_rate
        assert np.isclose(out_array[2], 1.)  # rel_sed_flux
        assert np.isclose(out_array[3], 1.)  # error_in_sed_flux_fn

        # zero flux in. Cell area is 10.
        get_sed_flux_function_pseudoimplicit_bysedout(0., 1.e10, 10., 10.,
                                                      sde._sed_flux_fn_gen,
                                                      humpk, humpnu,
                                                      humpc, humpphi,
                                                      humpnorm,
                                                      50, out_array)
        known_fqs = solver(
            0., humpk, humpnu, humpc, humpphi, humpnorm,
        )
        assert np.isclose(out_array[0], known_fqs*1., atol=0.001)  # dzbydt
        assert np.isclose(out_array[1], known_fqs*10., atol=0.001)
        # ^vol_pass_rate
        assert np.isclose(out_array[2], 0., atol=0.001)  # rel_sed_flux
        assert np.less(out_array[3], 0.001)  # error_in_sed_flux_fn

        # very low flux in. Cell area is 100 and erosion rate is 2.
        get_sed_flux_function_pseudoimplicit_bysedout(1., 1.e10, 200., 100.,
                                                      sde._sed_flux_fn_gen,
                                                      humpk, humpnu,
                                                      humpc, humpphi,
                                                      humpnorm,
                                                      50, out_array)
        known_fqs = solver(
            0., humpk, humpnu, humpc, humpphi, humpnorm,
        )
        assert np.isclose(out_array[0], known_fqs*2., atol=0.001)  # dzbydt
        assert np.isclose(out_array[1], known_fqs*200.+1., atol=0.001)
        # ^vol_pass_rate
        assert np.isclose(out_array[2], 0., atol=0.001)  # rel_sed_flux
        assert np.less(out_array[3], 0.001)  # error_in_sed_flux_fn

        # and a test where we produce a middling qs/qc but enough sed to
        # swamp the outlet
        get_sed_flux_function_pseudoimplicit_bysedout(1., 2., 2000., 2.,
                                                      sde._sed_flux_fn_gen,
                                                      humpk, humpnu,
                                                      humpc, humpphi,
                                                      humpnorm,
                                                      50, out_array)
        assert np.isclose(out_array[2], 0.75, atol=0.001)
        # ^mean of 1./2. (in) and 1 (out)
        known_fqs = solver(
            0.75, humpk, humpnu, humpc, humpphi, humpnorm,
        )
        assert np.isclose(out_array[1], 2., atol=0.001)  # output saturates
        assert np.isclose(out_array[0], (2. - 1.)/2., atol=0.001)
        # ^note this is now controlled by export limit, not fqs directly
        assert np.less(out_array[3], 0.001)  # error_in_sed_flux_fn

        # test where erosion can't happen as there's no erosion capacity
        get_sed_flux_function_pseudoimplicit_bysedout(1000., 2000., 0., 1000.,
                                                      sde._sed_flux_fn_gen,
                                                      humpk, humpnu,
                                                      humpc, humpphi,
                                                      humpnorm,
                                                      50, out_array)
        assert np.isclose(out_array[0], 0.)
        assert np.isclose(out_array[1], 1000.)
        assert np.isclose(out_array[2], 0.5)
        assert np.less(out_array[3], 0.001)

    # Now a few specific tests for each form

    # generic case, larger numbers
    mg = RasterModelGrid((5, 5))
    z = mg.add_zeros('node', 'topographic__elevation')
    fa = FlowAccumulator(mg)
    out_array = np.empty(4, dtype=float)
    sde = SedDepEroder(mg, sed_dependency_type='linear_decline')
    get_sed_flux_function_pseudoimplicit_bysedout(
        2500., 5000., 2500., 100000./3.,
        sde._sed_flux_fn_gen,
        humpk, humpnu, humpc, humpphi, humpnorm, 50, out_array
    )
    assert np.isclose(out_array[0], 0.025, atol=0.001)
    assert np.isclose(out_array[1], 3333.333, atol=1e-3)
    assert np.isclose(out_array[2], (1./2. + 2./3.)/2., atol=1.e-3)
    # ^1./2. = 0.5, which gives fqsqc = 0.5, supplies 1 unit sediment (A=2),
    # and so this immediately matches the initial conditions
    assert np.less(out_array[3], 0.001)  # error_in_sed_flux_fn


def test_iteration_dstr():
    """
    This tests the iterate_sde_downstream func.
    We can now test on an arbitrary structure, as we don't need grids!
    This test adopts the simplest func form, but still exercises all the
    crucial parts of the code.
    Recall that the eroder cannot produce sediment that it cannot transport.
    This has significant consequences for the operation of this model re SPL.
    True SPL equivalence for the func=1 case thus only occurs under high
    transport capacity.
    """
    pseudoimplicit_repeats = 50
    funct = sed_flux_fn_gen_const
    cell_areas = np.array([1., 0.1, 0.5, 1., 1., 1.])
    hillsl_sed = np.array([0., 0., 0., 0., 0., 0.])
    porosity = 1.
    # 0 and 1 drain to 2 then 3 then 4 then 5
    upstr_order = np.array([5, 4, 3, 2, 0, 1])
    flow_receiver = np.array([2, 2, 3, 4, 5, 5])
    trans_caps = np.array([0., 1., 1., 1., 0.9, 0.6])
    erosion_prefac_w_S = np.array([1., 2., 1., 1., 1., 1.])
    # output arrays:
    river_volume_flux_into_node = np.zeros(6, dtype=float)
    rel_sed_flux = np.zeros(6, dtype=float)
    is_it_TL = np.zeros(6, dtype=np.int8)
    vol_drop_rate = np.zeros(6, dtype=float)
    dzbydt = np.zeros(6, dtype=float)
    iterate_sde_downstream(
        upstr_order,
        cell_areas,
        hillsl_sed,
        porosity,
        river_volume_flux_into_node,
        trans_caps,
        erosion_prefac_w_S,
        rel_sed_flux,
        is_it_TL,
        vol_drop_rate,
        flow_receiver,
        pseudoimplicit_repeats,
        dzbydt,
        funct,
        0., 0., 0., 0., 0.
    )
    assert np.allclose(dzbydt, np.array([0., -2., -1., -0.3, 0., 0.]))
    assert np.all(np.equal(
        is_it_TL, np.array([1, 0, 0, 0, 1, 1], dtype=np.int8))
    )
    assert np.allclose(river_volume_flux_into_node,
                       np.array([0., 0., 0.2, 0.7, 1., 0.9]))
    assert np.allclose(rel_sed_flux, np.array([1., 0.1, 0.45, 0.85, 1., 1.]))
    assert np.allclose(vol_drop_rate, np.array([0., 0., 0., 0., 0.1, 0.3]))

    # now a very similar test where the capacities are really high
    # i.e., a true SP run.
    # ...also tests porosity is working OK.
    trans_caps.fill(1.e20)
    river_volume_flux_into_node.fill(0.)
    porosity = 2./3.
    iterate_sde_downstream(
        upstr_order,
        cell_areas,
        hillsl_sed,
        porosity,
        river_volume_flux_into_node,
        trans_caps,
        erosion_prefac_w_S,
        rel_sed_flux,
        is_it_TL,
        vol_drop_rate,
        flow_receiver,
        pseudoimplicit_repeats,
        dzbydt,
        funct,
        0., 0., 0., 0., 0.
    )
    assert np.allclose(dzbydt, np.array([-1., -2., -1., -1., -1., -1.]))
    assert np.all(np.equal(is_it_TL, np.int8(0)))
    assert np.allclose(river_volume_flux_into_node,
                       1.5 * np.array([0., 0., 1.2, 1.7, 2.7, 3.7]))
    # note the 1.5 is the effect of the sed porosity
    assert np.allclose(rel_sed_flux, 0., atol=1.e-10)
    assert np.allclose(vol_drop_rate, 0.)


def test_iteration_with_sed_in_channel():
    """
    This tests the iterate_sde_downstream func.
    It uses a lot of sediment in the channel to swamp the channels and halt
    incision.
    """
    pseudoimplicit_repeats = 50
    funct = sed_flux_fn_gen_const
    cell_areas = np.array([1., 0.1, 0.5, 1., 1., 1.])
    hillsl_sed = np.array([10., 10., 10., 10., 10., 10.])
    porosity = 1.
    # 0 and 1 drain to 2 then 3 then 4 then 5
    upstr_order = np.array([5, 4, 3, 2, 0, 1])
    flow_receiver = np.array([2, 2, 3, 4, 5, 5])
    trans_caps = np.array([0., 1., 1., 1., 0.9, 0.6])
    erosion_prefac_w_S = np.array([1., 2., 1., 1., 1., 1.])
    # output arrays:
    river_volume_flux_into_node = np.zeros(6, dtype=float)
    rel_sed_flux = np.zeros(6, dtype=float)
    is_it_TL = np.zeros(6, dtype=np.int8)
    vol_drop_rate = np.zeros(6, dtype=float)
    dzbydt = np.zeros(6, dtype=float)
    iterate_sde_downstream(
        upstr_order,
        cell_areas,
        hillsl_sed,
        porosity,
        river_volume_flux_into_node,
        trans_caps,
        erosion_prefac_w_S,
        rel_sed_flux,
        is_it_TL,
        vol_drop_rate,
        flow_receiver,
        pseudoimplicit_repeats,
        dzbydt,
        funct,
        0., 0., 0., 0., 0.
    )
    assert np.allclose(dzbydt, 0.)
    assert np.all(is_it_TL == 1)
    assert np.allclose(rel_sed_flux, 1.)
    assert np.allclose(river_volume_flux_into_node, np.array(
        [0., 0., 1., 1., 1., 0.9]
    ))
    assert np.allclose(vol_drop_rate, np.array(
        [10., 9., 10., 10., 10.1, 10.3]
    ))


def test_plotting():
    x_pts = np.arange(0., 1.01, 0.01)
    for fname, funct in zip(
        ['generalized_humped', 'None', 'linear_decline', 'almost_parabolic'],
        [sed_flux_fn_gen_genhump, sed_flux_fn_gen_const,
         sed_flux_fn_gen_lindecl, sed_flux_fn_gen_almostparabolic]
    ):
        y_pts = np.empty_like(x_pts)
        for i in range(101):
            y_pts[i] = funct(
                x_pts[i], 13.683, 1.13, 0.00181, 4.24, 1.0000278041373
            )
        print(y_pts)
        mg = RasterModelGrid((5, 5))
        z = mg.add_zeros('node', 'topographic__elevation')
        fa = FlowAccumulator(mg)
        sde = SedDepEroder(mg, sed_dependency_type=fname)
        sde.show_sed_flux_function()
        x_plot, y_plot = gca().lines[0].get_xydata().T
        assert np.allclose(x_plot, x_pts)
        assert np.allclose(y_plot, y_pts)
        clf()


def test_instantiation_trp_laws():
    for bad_Qc in ['MPM', 'Voller_generalized', 'bad_name']:
        mg = RasterModelGrid((5, 5))
        z = mg.add_zeros('node', 'topographic__elevation')
        fa = FlowAccumulator(mg)
        with pytest.raises(NameError):
            sde = SedDepEroder(mg, Qc=bad_Qc)

    mg = RasterModelGrid((5, 5))
    z = mg.add_zeros('node', 'topographic__elevation')
    fa = FlowAccumulator(mg)
    sde = SedDepEroder(mg, Qc='power_law')


def test_instantiation_sed_flux_forms():
    """
    Tests the correct behaviour for the sed_dependency_type term on
    instantiation. This also implicitly but adequately tests the
    set_sed_flux_fn_gen method.
    """
    # Create a fail by supplying a bad term
    mg = RasterModelGrid((5, 5))
    z = mg.add_zeros('node', 'topographic__elevation')
    fa = FlowAccumulator(mg)
    with pytest.raises(NameError):
        sde = SedDepEroder(mg, sed_dependency_type='bad_term')

    for sde_term, sde_fn in zip(
        ('None', 'linear_decline', 'almost_parabolic'),
        (
            sed_flux_fn_gen_const,
            sed_flux_fn_gen_lindecl,
            sed_flux_fn_gen_almostparabolic
        )
    ):
        mg = RasterModelGrid((5, 5))
        z = mg.add_zeros('node', 'topographic__elevation')
        fa = FlowAccumulator(mg)
        sde = SedDepEroder(mg, sed_dependency_type=sde_term)
        assert sde._sed_flux_fn_gen is sde_fn
        assert sde.kappa == 0.
        assert sde.nu == 0.
        assert sde.phi == 0.
        assert sde.c == 0.
        assert sde.norm == 0.

    mg = RasterModelGrid((5, 5))
    z = mg.add_zeros('node', 'topographic__elevation')
    fa = FlowAccumulator(mg)
    sde = SedDepEroder(mg, sed_dependency_type='generalized_humped')
    assert sde._sed_flux_fn_gen is sed_flux_fn_gen_genhump
    assert sde.kappa > 0.
    assert sde.nu > 0.
    assert sde.phi > 0.
    assert sde.c > 0.
    assert sde.norm > 0.


def test_correct_field_input_responses():
    mg = RasterModelGrid((5, 5))
    with pytest.raises(FieldError) as excinfo:
        msg = (
            "In order for the SedDepEroder to work, you must supply a " +
            "topographic__elevation field."
        )
        sde = SedDepEroder(mg)
        assert msg in str(excinfo.value)

    for bad_in_field in (
        "drainage_area",
        "flow__receiver_node",
        "flow__upstream_node_order",
        "topographic__steepest_slope",
        "flow__link_to_receiver_node",
        "flow__sink_flag"
    ):
        mg = RasterModelGrid((5, 5))
        z = mg.add_zeros('node', 'topographic__elevation')
        with pytest.raises(FieldError):
            sde = SedDepEroder(mg)

    # better check the raised error looks right for one of these...
    mg = RasterModelGrid((5, 5))
    z = mg.add_zeros('node', 'topographic__elevation')
    fa = FlowAccumulator(mg)
    catchit = mg.at_node.pop('flow__receiver_node')
    with pytest.raises(FieldError) as excinfo:
        msg = (
            "In order for the SedDepEroder to work, you must " +
            "supply the field flow__receiver_node. You probably need to " +
            "instantiate a FlowAccumulator component *prior* to " +
            "instantiating the SedDepEroder."
        )
        sde = SedDepEroder(mg)
        assert msg in str(excinfo.value)

    # test the channel_sediment__depth field is created and/or bound correctly
    # create the field whole cloth
    mg = RasterModelGrid((5, 5))
    z = mg.add_zeros('node', 'topographic__elevation')
    fa = FlowAccumulator(mg)
    sde = SedDepEroder(mg)
    assert np.allclose(mg.at_node['channel_sediment__depth'], 0.)
    # check bindings
    mg.at_node['channel_sediment__depth'] += 1.
    assert np.allclose(sde._hillslope_sediment, 1.)
    # check binding for an existing field
    mg = RasterModelGrid((5, 5))
    z = mg.add_zeros('node', 'topographic__elevation')
    z += np.random.rand(25) * 1.e-6
    fa = FlowAccumulator(mg)
    d = mg.add_ones('node', 'channel_sediment__depth')
    sde = SedDepEroder(mg)
    assert sde._hillslope_sediment is d
    # check binding is retained through a run cycle
    fa.run_one_step()
    sde.run_one_step(1.e-6)
    assert sde._hillslope_sediment is d


def test_basic_functionality():
    """
    This test runs a single short timestep to match a predictable model
    outcome.
    """
    # Very large capacity; only a small change magnitude
    for sff_type in ['None', 'linear_decline']:
        mg = RasterModelGrid((4, 7))
        closed_nodes = np.array(
            [True,  True,  True,  True,  True,  True,  True,
             True, False, False, False, False, False, False,
             True,  True, False,  True,  True,  True,  True,
             True,  True,  True,  True,  True,  True,  True], dtype=bool
        )
        mg.status_at_node[closed_nodes] = CLOSED_BOUNDARY
        X = np.spacing(4.)
        z_init = np.array(
            [0.,    0.,    0.,    0.,    0.,    0.,    0.,
             0.,    5.,    4.,    3.,    2.,    1.,    0.,
             0.,    0.,  4.+X,    0.,    0.,    0.,    0.,
             0.,    0.,    0.,    0.,    0.,    0.,    0.]
        )
        z = mg.add_field('node', 'topographic__elevation', z_init, copy=True)
        fa = FlowAccumulator(mg)
        sde = SedDepEroder(
            mg, K_sp=1.e-6, K_t=1.e10, m_sp=1., sed_dependency_type=sff_type
        )
        fa.run_one_step()
        sde.run_one_step(1.)
        small_linear_incision = np.array(
            [0.,  0.,  0.,  0.,  0.,  0.,  0.,
             0.,  1.,  3.,  4.,  5.,  6.,  0.,
             0.,  0.,  0.,  0.,  0.,  0.,  0.,
             0.,  0.,  0.,  0.,  0.,  0.,  0.]
        )
        # makes no difference between linear decline & const
        assert np.allclose(
            (z_init - z) * 1.e6, small_linear_incision, atol=1.e-10
        )
        assert np.allclose(
            mg.at_node['channel_sediment__relative_flux'][
                np.array([8, 9, 10, 11, 12])
            ], 0., atol=1.e-10
        )
        # a couple of special cases:
        # the basically flat node gets 1 by definition
        assert np.isclose(
            mg.at_node['channel_sediment__relative_flux'][16], 1.
        )
        # the fixed elev gets 1 (as do the closed nodes)
        assert np.isclose(
            mg.at_node['channel_sediment__relative_flux'][13], 1.
        )
        # now, values aren't important, but we have values in the right places
        assert np.all(np.greater(
            mg.at_node['channel_sediment__volumetric_transport_capacity'][
                mg.core_nodes
            ], 0.
        ))
        assert np.all(np.greater(
            mg.at_node['channel_sediment__volumetric_discharge'][
                mg.at_node['drainage_area'] > 1.
            ], 0.
        ))  # ...because if A==1, there's no sed discharge coming in
        assert np.allclose(
            mg.at_node['channel_sediment__depth'][mg.core_nodes], 0.
        )

    # & we can do the same thing with the aparabolic:
    mg = RasterModelGrid((4, 7))
    closed_nodes = np.array(
        [True,  True,  True,  True,  True,  True,  True,
         True, False, False, False, False, False, False,
         True,  True, False,  True,  True,  True,  True,
         True,  True,  True,  True,  True,  True,  True], dtype=bool
    )
    mg.status_at_node[closed_nodes] = CLOSED_BOUNDARY
    X = np.spacing(4.)
    z_init = np.array(
        [0.,    0.,    0.,    0.,    0.,    0.,    0.,
         0.,    5.,    4.,    3.,    2.,    1.,    0.,
         0.,    0.,  4.+X,    0.,    0.,    0.,    0.,
         0.,    0.,    0.,    0.,    0.,    0.,    0.]
    )
    z = mg.add_field('node', 'topographic__elevation', z_init, copy=True)
    fa = FlowAccumulator(mg)
    sde = SedDepEroder(
        mg, K_sp=1.e-5, K_t=1.e10, m_sp=1.,
        sed_dependency_type='almost_parabolic'
    )  # note 1.e-5 not 1.e-6 now
    fa.run_one_step()
    sde.run_one_step(1.)
    small_linear_incision = np.array(
        [0.,  0.,  0.,  0.,  0.,  0.,  0.,
         0.,  1.,  3.,  4.,  5.,  6.,  0.,
         0.,  0.,  0.,  0.,  0.,  0.,  0.,
         0.,  0.,  0.,  0.,  0.,  0.,  0.]
    )
    assert np.allclose(
        (z_init - z) * 1.e6, small_linear_incision, atol=1.e-10
    )

    # similar end member test - flood the terrain w sed & use lin decl
    mg = RasterModelGrid((4, 7))
    closed_nodes = np.array(
        [True,  True,  True,  True,  True,  True,  True,
         True, False, False, False, False, False, False,
         True,  True, False,  True,  True,  True,  True,
         True,  True,  True,  True,  True,  True,  True], dtype=bool
    )
    mg.status_at_node[closed_nodes] = CLOSED_BOUNDARY
    X = np.spacing(4.)
    z_init = np.array(
        [0.,    0.,    0.,    0.,    0.,    0.,    0.,
         0.,    5.,    4.,    3.,    2.,    1.,    0.,
         0.,    0.,  4.+X,    0.,    0.,    0.,    0.,
         0.,    0.,    0.,    0.,    0.,    0.,    0.]
    )
    z = mg.add_field('node', 'topographic__elevation', z_init, copy=True)
    fa = FlowAccumulator(mg)
    sde = SedDepEroder(
        mg, K_sp=1., K_t=1.e-20, sed_dependency_type='linear_decline'
    )  # note 1.e-5 not 1.e-6 now
    fa.run_one_step()
    sde.run_one_step(1.)
    assert np.allclose(z_init - z, 0., atol=1.e-10)


def test_supplied_sediment():
    """
    This replicates a simple test, but instead of flooding the system with
    eroded sediment, it does it with an external supply.

    Problems here - erosion is occurring at closed nodes, and we aren't
    saturating.
    """
    mg = RasterModelGrid((4, 7), xy_spacing=100.)
    closed_nodes = np.array(
        [True,  True,  True,  True,  True,  True,  True,
         True, False, False, False, False, False, False,
         True,  True, False,  True,  True,  True,  True,
         True,  True,  True,  True,  True,  True,  True], dtype=bool
    )
    mg.status_at_node[closed_nodes] = CLOSED_BOUNDARY
    X = np.spacing(4.)
    z_init = np.array(
        [0.,    0.,    0.,    0.,    0.,    0.,    0.,
         0.,    5.,    4.,    3.,    2.,    1.,    0.,
         0.,    0.,  4.+X,    0.,    0.,    0.,    0.,
         0.,    0.,    0.,    0.,    0.,    0.,    0.]
    )
    z = mg.add_field('node', 'topographic__elevation', z_init, copy=True)
    h = mg.add_ones('node', 'channel_sediment__depth')
    h *= 10.
    fa = FlowAccumulator(mg)
    sde = SedDepEroder(
        mg, K_sp=1.e-4, K_t=1., m_sp=0., n_sp=1., m_t=0., n_t=1.,
        sed_dependency_type='linear_decline'
    )
    fa.run_one_step()
    sde.run_one_step(10.)
    assert np.allclose(z_init - z, 0., atol=1.e-10)
    assert np.allclose(h[np.logical_and(
            mg.at_node['drainage_area'] > 10000.,
            mg.status_at_node == 0
        )], 10.)
    # ^all nodes not at the head are in st st, so keep all sed, and
    assert np.isclose(
        sde.calc_sed_discharge_from_node()[8] * 31557600.,
        1. * 10000.**0. * 0.01**1.
    )
    assert np.isclose(h[8], 10. - 0.01 / mg.area_of_cell[8])
    assert np.allclose(
        mg.at_node['channel_sediment__relative_flux'][mg.core_nodes], 1.
    )


def test_diagonal_route():
    """
    Tests that things make sense on the diagonal.
    """
    mg = RasterModelGrid((5, 5), xy_spacing=(300, 400))
    closed_nodes = np.array(
        [True,  True,  True,  True, False,
         True,  True,  True, False,  True,
         True,  True, False,  True,  True,
         True, False,  True,  True,  True,
         True,  True,  True,  True,  True], dtype=bool
    )
    mg.status_at_node[closed_nodes] = CLOSED_BOUNDARY
    z_init = np.array(
        [0.,  0.,  0.,  0.,  0.,
         0.,  0.,  0.,  5.,  0.,
         0.,  0., 10.,  0.,  0.,
         0., 15.,  0.,  0.,  0.,
         0.,  0.,  0.,  0.,  0.], dtype=float
    )
    z = mg.add_field('node', 'topographic__elevation', z_init, copy=True)
    h = mg.add_ones('node', 'channel_sediment__depth')
    h *= 10.
    fa = FlowAccumulator(mg, flow_director="D8")
    sde = SedDepEroder(
        mg, K_sp=1.e-4, K_t=1., m_sp=0., n_sp=1., m_t=0., n_t=1.,
        sed_dependency_type='linear_decline'
    )
    fa.run_one_step()
    sde.run_one_step(10.)
    assert np.allclose(
        mg.at_node['channel_sediment__volumetric_transport_capacity'][
            mg.core_nodes
        ], 0.01/31557600.)
    assert np.allclose(
        mg.at_node['channel_sediment__volumetric_discharge'][
            mg.core_nodes
        ] * 31557600., np.array([0.01, 0.01, 0.]))  # 8, 12, 16


# # add a regression test :
# def test_arbitrary_grid():
#     """
#     This runs stably but wrongly, though unclear what's up as profiles look OK
#     """
#     from landlab.components import LinearDiffuser
#     mg = RasterModelGrid((50, 50), xy_spacing=100.)
#     for edge in ['top', 'left', 'bottom']:
#         mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
#     z = mg.add_zeros('node', 'topographic__elevation')
#     np.random.seed(50)
#     z += np.random.rand(mg.number_of_nodes) / 1000.
#     z_init = z.copy()
#     dfn = LinearDiffuser(mg, linear_diffusivity=1.e-2)
#     fa = FlowAccumulator(mg, flow_director="D8")
#     pit = DepressionFinderAndRouter(mg, routing="D8")
#     sde = SedDepEroder(mg, K_sp=1.e-4, K_=1.e-4, sed_dependency_type='almost_parabolic')
#     sp = FastscapeEroder(mg, K_sp=1.e-4)
#     th = mg.at_node['channel_sediment__depth']
#     dt = 1000.
#     U = 0.001
#     # burn in a good network
#     for i in range(500):
#         z[mg.core_nodes] += 0.001 * U * dt
#         fa.run_one_step()
#         sp.run_one_step(dt)
#     for i in range(600):
#         z[mg.core_nodes] += U * dt
#         # z_init[:] = z
#         # dfn.run_one_step(dt)
#         # dz = z - z_init  # -ve when lowering
#         # th += dz
#         # np.clip(th, a_max=None, a_min=0., out=th)
#         # z -= th
#         # th[mg.core_nodes] += 0.1 * U * dt  # a st st soil prod of 0.2*U
#         fa.run_one_step()
#         pit.map_depressions()
#         sde.run_one_step(dt, flooded_nodes=pit.lake_at_node)
#         # z += th
#         print(i)
#         if i%25 == 24:
#             figure(i)
#             imshow_grid_at_node(mg, z)
#             figure(i+600)
#             imshow_grid_at_node(mg, 'channel_sediment__depth')
#     # --> Issue here is that top row, in particular top right, is high
#     # relative to rest of grid. Something is screwy here.
#     # Problem goes away if sed_dependency_type='linear_decline'
#     # Check convergence rules where the dstr node is flooded
#     # All this is emerging from the piling up of large amounts of sediment
#     # (10s m) in certain cells, which either diverts the flow (if BR surface
#     # is altered) or destabilises the SFF (topo is altered) in the next
#     # iteration
# 
# 
# # add a regression test :
# def test_arbitrary_grid_with_flood():
#     """
#     This does not run stably, though unclear what's up as profiles look OK
#     """
#     from landlab.components import LinearDiffuser
#     mg = RasterModelGrid((30, 30), xy_spacing=500.)
#     for edge in ['top', 'left', 'bottom']:
#         mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
#     z = mg.add_zeros('node', 'topographic__elevation')
#     np.random.seed(50)
#     z += np.random.rand(mg.number_of_nodes) / 1000.
#     z_init = z.copy()
#     dfn = LinearDiffuser(mg, linear_diffusivity=1.e-2)
#     fa = FlowAccumulator(mg, flow_director="D8")
#     pit = DepressionFinderAndRouter(mg, routing="D8")
#     sde = SedDepEroder(mg, K_sp=1.e-4)
#     th = mg.at_node['channel_sediment__depth']
#     dt = 1000.
#     U = 0.001
#     for i in range(1000):
#         z[mg.core_nodes] += U * dt
#         # z_init[:] = z
#         # dfn.run_one_step(dt)
#         # dz = z - z_init  # -ve when lowering
#         # th += dz
#         # np.clip(th, a_max=None, a_min=0., out=th)
#         # z -= th
#         th[mg.core_nodes] += 0.1 * U * dt  # a st st soil prod of 0.2*U
#         fa.run_one_step()
#         pit.map_depressions()
#         sde.run_one_step(dt, flooded_nodes=pit.lake_at_node)
#         # z += th
#         print(i)
#     # --> Issue here is that top row, in particular top right, is high
#     # relative to rest of grid. Something is screwy here.


def test_flooding():
    mg = RasterModelGrid((3, 7), xy_spacing=100.)
    for edge in ['top', 'left', 'bottom']:
        mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
    z = mg.add_zeros('node', 'topographic__elevation')
    z[mg.core_nodes] = np.array([2.01, 1.999, 1.998, 2., 1.])
    fa = FlowAccumulator(mg, flow_director="D8")
    pit = DepressionFinderAndRouter(mg, routing="D8")
    sde = SedDepEroder(mg, K_sp=1.e-3, K_t=1.e-4)

    fa.run_one_step()
    pit.map_depressions()
    sde.run_one_step(100., flooded_nodes=pit.lake_at_node)
    assert np.allclose(z[mg.core_nodes], np.array(
        [2.00997276, 1.999,  1.998, 1.99502899, 0.95723706]
    ))
    assert np.allclose(
        mg.at_node['channel_sediment__volumetric_transport_capacity'][
            mg.core_nodes
        ],
        np.array([3.47796809e-10, 0., 0., 2.62050580e-07, 3.40770733e-07])
    )
    assert np.allclose(
        mg.at_node['channel_sediment__depth'][mg.core_nodes],
        np.array([0., 2.72372391e-05, 0., 0., 0.])
    )

    fa.run_one_step()
    pit.map_depressions()
    sde.run_one_step(100., flooded_nodes=pit.lake_at_node)
    assert np.allclose(z[mg.core_nodes], np.array(
        [2.00994559, 1.999, 1.99781863, 1.92443714, 0.92132846])
    )
    assert np.allclose(
        mg.at_node['channel_sediment__volumetric_transport_capacity'][
            mg.core_nodes
        ],
        np.array([
            3.46935625e-10,
            9.97358068e-11,
            1.18582125e-08,
            2.52251863e-07,
            3.29771153e-07
        ])
    )
    assert np.allclose(
        mg.at_node['channel_sediment__depth'][mg.core_nodes],
        np.array([0., 2.57859919e-05, 0., 0., 0.])
    )

# def test_large_steps_for_timestepping():
#     pass




############################## Need a test for supplied Voronoi grid









# def test_sed_dep_new_almostpara_fullrun():
#     """
#     This tests only the power_law version of the SDE, using the
#     almost_parabolic form of f(Qs).
#     """
#     mg = RasterModelGrid((10, 5), xy_spacing=200.)
#     for edge in (mg.nodes_at_left_edge, mg.nodes_at_top_edge,
#                  mg.nodes_at_right_edge):
#         mg.status_at_node[edge] = CLOSED_BOUNDARY
#
#     z = mg.add_zeros('node', 'topographic__elevation')
#
#     fr = FlowAccumulator(mg, flow_director='D8')
#     sde = SedDepEroder(mg, K_sp=1.e-4, sed_dependency_type='almost_parabolic',
#                        Qc='power_law', K_t=1.e-4)
#
#     z[:] = mg.node_y/10000.
#     z.reshape((10, 5))[:, 2] *= 2.
#     z.reshape((10, 5))[:, 3] *= 2.1
#
#     initz = z.copy()
#
#     dt = 100.
#     up = 0.05
#
#     for i in range(1):
#         fr.run_one_step()
#         sde.run_one_step(dt)
#
#     assert_array_almost_equal(mg.at_node['channel_sediment__depth'],
#                               np.array([0.,  0.,  0.,  0.,  0.,
#                                         0.,  0.,  0.,  0.,  0.,
#                                         0.,  0.,  0.,  0.,  0.,
#                                         0.,  0.,  0.,  0.,  0.,
#                                         0.,  0.,  0.,  0.,  0.,
#                                         0.,  0.,  0.,  0.,  0.,
#                                         0.,  0.,  0.,  0.,  0.,
#                                         0.,  0.,  0.,  0.,  0.,
#                                         0.,  0.00023431,  0.,  0.,  0.,
#                                         0.,  0.,  0.,  0.,  0.]))
#
#     assert_array_almost_equal(z - initz, np.array(
#         [0.,  0.,  0.,  0.,  0.,
#          0., -0.00076662, -0.0004,     -0.00118507,  0.,
#          0., -0.00068314, -0.00042426, -0.00110741,  0.,
#          0., -0.00065571, -0.0006,     -0.00102346,  0.,
#          0., -0.00056135, -0.0008,     -0.00093111,  0.,
#          0., -0.00045180, -0.0010,     -0.00078882,  0.,
#          0., -0.00031903, -0.0012,     -0.00071743,  0.,
#          0., -0.00015763, -0.0014,     -0.00057541,  0.,
#          0.,  0.00023431, -0.0016,     -0.00042,     0.,
#          0.,  0.,  0.,  0.,  0.]))
#
#     assert_equal(len(sde._error_at_abort), 0)  # good convergence at all nodes
#
#
# def test_sed_dep_new_genhumped_fullrun():
#     """
#     This tests only the power_law version of the SDE, using the
#     generalized_humped form of f(Qs).
#     """
#     mg = RasterModelGrid((10, 5), xy_spacing=200.)
#     for edge in (mg.nodes_at_left_edge, mg.nodes_at_top_edge,
#                  mg.nodes_at_right_edge):
#         mg.status_at_node[edge] = CLOSED_BOUNDARY
#
#     z = mg.add_zeros('node', 'topographic__elevation')
#
#     fr = FlowAccumulator(mg, flow_director='D8')
#     sde = SedDepEroder(mg, K_sp=1.e-4,
#                        sed_dependency_type='generalized_humped',
#                        Qc='power_law', K_t=1.e-4)
#
#     z[:] = mg.node_y/10000.
#     z.reshape((10, 5))[:, 2] *= 2.
#     z.reshape((10, 5))[:, 3] *= 2.1
#
#     initz = z.copy()
#
#     dt = 100.
#     up = 0.05
#
#     for i in range(1):
#         fr.run_one_step()
#         sde.run_one_step(dt)
#
#     ans = np.array([0.,  0.,  0.,  0.,  0.,
#                     0.02,  0.01940896,  0.03969863,  0.04112046,  0.02,
#                     0.04,  0.03948890,  0.07968035,  0.08318039,  0.04,
#                     0.06,  0.05951858,  0.11954794,  0.12524486,  0.06,
#                     0.08,  0.07960747,  0.15939726,  0.16731501,  0.08,
#                     0.10,  0.09969946,  0.19924657,  0.20939250,  0.10,
#                     0.12,  0.11979225,  0.23909589,  0.25147973,  0.12,
#                     0.14,  0.13987971,  0.27894520,  0.29357949,  0.14,
#                     0.16,  0.16023431,  0.31879451,  0.33568356,  0.16,
#                     0.18,  0.18,  0.36,  0.378,  0.18])
#
#     assert_array_almost_equal(z, ans)
#
#
# def test_sed_dep_new_lindecl_fullrun():
#     """
#     This tests only the power_law version of the SDE, using the
#     linear_decline form of f(Qs).
#     """
#     mg = RasterModelGrid((10, 5), xy_spacing=200.)
#     for edge in (
#         mg.nodes_at_left_edge, mg.nodes_at_top_edge, mg.nodes_at_right_edge
#     ):
#         mg.status_at_node[edge] = CLOSED_BOUNDARY
#
#     z = mg.add_zeros("node", "topographic__elevation")
#
#     fr = FlowAccumulator(mg, flow_director='D8')
#     sde = SedDepEroder(
#         mg,
#         K_sp=1.e-4,
#         sed_dependency_type='linear_decline',
#         Qc='power_law',
#         K_t=1.e-4)
#
#     z[:] = mg.node_y/10000.
#     z.reshape((10, 5))[:, 2] *= 2.
#     z.reshape((10, 5))[:, 3] *= 2.1
#
#     initz = z.copy()
#
#     dt = 100.
#     up = 0.05
#
#     for i in range(1):
#         fr.run_one_step()
#         sde.run_one_step(dt)
#
#     ans = np.array([0.,  0.,  0.,  0.,  0.,
#                     0.02,  0.01955879,  0.03980000,  0.04128996,  0.02,
#                     0.04,  0.03961955,  0.07978787,  0.08333878,  0.04,
#                     0.06,  0.05964633,  0.11970000,  0.12539158,  0.06,
#                     0.08,  0.07971138,  0.15960000,  0.16745207,  0.08,
#                     0.10,  0.09978034,  0.19950000,  0.20951474,  0.10,
#                     0.12,  0.11985392,  0.23940000,  0.25158663,  0.12,
#                     0.14,  0.13993079,  0.27930000,  0.29367338,  0.14,
#                     0.16,  0.16023431,  0.31920000,  0.33579000,  0.16,
#                     0.18,  0.18,  0.36,  0.378,  0.18])
#
#     assert_array_almost_equal(z, ans)
#
#
# def test_sed_dep_new_const_fullrun():
#     """
#     This tests only the power_law version of the SDE, using the
#     constant (None) form of f(Qs).
#     """
#     mg = RasterModelGrid((10, 5), xy_spacing=200.)
#     for edge in (
#         mg.nodes_at_left_edge, mg.nodes_at_top_edge, mg.nodes_at_right_edge
#     ):
#         mg.status_at_node[edge] = CLOSED_BOUNDARY
#
#     z = mg.add_zeros('node', 'topographic__elevation')
#
#     fr = FlowAccumulator(mg, flow_director='D8')
#     sde = SedDepEroder(
#         mg,
#         K_sp=1.e-4,
#         sed_dependency_type='None',
#         Qc='power_law',
#         K_t=1.e-4
#     )
#
#     z[:] = mg.node_y/10000.
#     z.reshape((10, 5))[:, 2] *= 2.
#     z.reshape((10, 5))[:, 3] *= 2.1
#
#     initz = z.copy()
#
#     dt = 100.
#     up = 0.05
#
#     for i in range(1):
#         fr.run_one_step()
#         sde.run_one_step(dt)
#
#     ans = np.array([0.,  0.,  0.,  0.,  0.,
#                     0.02,  0.01922540,  0.03960000,  0.04081206,  0.02,
#                     0.04,  0.03927889,  0.07957574,  0.08288878,  0.04,
#                     0.06,  0.05930718,  0.11940000,  0.12497121,  0.06,
#                     0.08,  0.07936754,  0.15920000,  0.16706085,  0.08,
#                     0.10,  0.09943431,  0.19900000,  0.20916000,  0.10,
#                     0.12,  0.11951010,  0.23880000,  0.25127254,  0.12,
#                     0.14,  0.13960000,  0.27860000,  0.29340603,  0.14,
#                     0.16,  0.16023431,  0.31840000,  0.33558000,  0.16,
#                     0.18,  0.18,  0.36,  0.378,  0.18])
#
#     assert_array_almost_equal(z, ans)
#
#
# def test_sed_dep_w_hillslopes():
#     """
#     This tests only the power_law version of the SDE, with a hillslope input.
#     """
#     mg = RasterModelGrid((10, 5), xy_spacing=200.)
#     for edge in (
#         mg.nodes_at_left_edge, mg.nodes_at_top_edge, mg.nodes_at_right_edge
#     ):
#         mg.status_at_node[edge] = CLOSED_BOUNDARY
#
#     z = mg.add_zeros('node', 'topographic__elevation')
#     th = mg.add_zeros('node', 'channel_sediment__depth')
#     th[mg.core_nodes] += 0.001
#
#     fr = FlowAccumulator(mg, flow_director='D8')
#     sde = SedDepEroder(
#         mg,
#         K_sp=1.e-4,
#         sed_dependency_type='almost_parabolic',
#         Qc='power_law',
#         K_t=1.e-4
#     )
#
#     z[:] = mg.node_y/10000.
#     z.reshape((10, 5))[:, 2] *= 2.
#     z.reshape((10, 5))[:, 3] *= 2.1
#
#     initz = z.copy()
#
#     dt = 100.
#     up = 0.05
#
#     for i in range(1):
#         print(i)
#         fr.run_one_step()
#         sde.run_one_step(dt)
#
#     # test binding of field occurs correctly:
#     assert th is mg.at_node['channel_sediment__depth']
#
#     assert_array_almost_equal(
#         mg.at_node['channel_sediment__depth'],
#         np.array([0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
#                   0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
#                   0.00000000e+00,   6.00000000e-04,   0.00000000e+00,
#                   0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
#                   5.75735931e-04,   0.00000000e+00,   0.00000000e+00,
#                   0.00000000e+00,   0.00000000e+00,   4.00000000e-04,
#                   0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
#                   9.28079257e-07,   2.00000000e-04,   0.00000000e+00,
#                   0.00000000e+00,   0.00000000e+00,   4.13904292e-04,
#                   0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
#                   0.00000000e+00,   7.60612309e-04,   0.00000000e+00,
#                   5.55537486e-06,   0.00000000e+00,   0.00000000e+00,
#                   1.16568542e-03,   0.00000000e+00,   2.32060608e-04,
#                   0.00000000e+00,   0.00000000e+00,   1.73431458e-03,
#                   0.00000000e+00,   5.80000000e-04,   0.00000000e+00,
#                   0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
#                   0.00000000e+00,   0.00000000e+00]))
#
#     assert_array_almost_equal(z, np.array(
#         [0.,  0.,  0.,  0.,  0.,
#          0.02,  0.01877957,  0.0396,      0.04055596,  0.02,
#          0.04,  0.03891440,  0.07957574,  0.08263689,  0.04,
#          0.06,  0.05890177,  0.1194,      0.12472345,  0.06,
#          0.08,  0.07900093,  0.1592,      0.16681590,  0.08,
#          0.10,  0.09941390,  0.1990,      0.20891231,  0.10,
#          0.12,  0.11976061,  0.23863333,  0.25100556,  0.12,
#          0.14,  0.14016569,  0.27831429,  0.29323206,  0.14,
#          0.16,  0.16073431,  0.318025,    0.335580,    0.16,
#          0.18,  0.18,  0.36,  0.378,  0.18]))
#
#     # good convergence at all nodes
#     assert_equal(len(sde._error_at_abort), 0)
