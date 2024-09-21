"""
Created on Wed Jul 12 12:25:27 2023

@author: LaurentRoberge
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from landlab import FieldError
from landlab import NodeStatus
from landlab import RasterModelGrid
from landlab.components import ConcentrationTrackerForSpace


### ~~~ TASK ~~~ ###

# THESE TESTS SHOULD ALL FAIL BECAUSE I HAVE NOT YET INCLUDED THE THREE OTHER
# REQUIRED INPUTS: phi, fraction_fines, settling_velocity.

# I NEED TO ADD A TEST FOR THESE INPUTS: phi, fraction_fines, settling_velocity.

class TestInputParameters:
    """Test input errors"""

    def setup_method(self):
        self.mg = RasterModelGrid((3, 3))
        self.mg.add_zeros("soil__depth", at="node")
        self.mg.add_zeros("sediment__outflux", at="node")
        self.mg.add_zeros("bedrock__erosion_flux", at="node")
        self.mg.add_zeros("sediment__erosion_flux", at="node")
        self.mg.add_zeros("sediment__deposition_flux", at="node")
        self.mg.add_zeros("topographic__elevation", at="node")

    @pytest.mark.parametrize(
        "at,required_field",
        [
            ("node", "soil__depth"),
            ("node", "sediment__outflux")
            ("node", "bedrock__erosion_flux")
            ("node", "sediment__erosion_flux")
            ("node", "sediment__deposition_flux")
            ("node", "topographic__elevation"),
        ],
    )
    def test_input_fields_soil(self, at, required_field):
        """
        ConcentrationTrackerForSpace should throw an error when input fields
        are not provided (soil__depth, sediment__outflux, bedrock__erosion_flux,
        sediment__erosion_flux, sediment__deposition_flux, topographic__elevation)
        """
        self.mg[at].pop(required_field)
        with pytest.raises(FieldError):
            ConcentrationTrackerForSpace(self.mg)

    def test_field_instantiation(self):
        """
        ConcentrationTrackerForSpace should instantiate the following fields
        when they do not already exist ('bedrock_property__concentration' and
        'sediment_property__concentration')
        """
        ConcentrationTrackerForSpace(self.mg)

        missing_fields = {
            "bedrock_property__concentration",
            "sediment_property__concentration",
        } - set(self.mg.at_node)
        assert not missing_fields

    @pytest.mark.parametrize(
        "name",
        ["sediment_property__concentration", "bedrock_property__concentration"],
    )
    def test_fields_for_default_input(self, name):
        """Check default input produces correct fields with no pre-existing fields"""
        ConcentrationTrackerForSpace(self.mg)

        assert np.allclose(self.mg.at_node[name], 0.0)

    @pytest.mark.parametrize(
        "name",
        ["sediment_property__concentration", "bedrock_property__concentration"],
    )
    def test_fields_for_default_input_with_preexisting_fields(self, name):
        """Check default input uses correct fields with pre-existing fields"""
        initial = self.mg.add_field(
            name, np.random.uniform(size=self.mg.number_of_nodes), at="node"
        ).copy()

        ConcentrationTrackerForSpace(self.mg)

        assert np.allclose(self.mg.at_node[name], initial)

    @pytest.mark.parametrize(
        "field,keyword",
        [
            ("sediment_property__concentration", "concentration_initial"),
            ("bedrock_property__concentration", "concentration_in_bedrock"),
        ],
    )
    def test_fields_for_user_value_input(self, field, keyword):
        """Check user input of single values produces the correct fields"""
        ConcentrationTrackerForSpace(
            self.mg, **{keyword: (initial_value := np.random.uniform())}
        )

        assert np.allclose(self.mg.at_node[field], initial_value)

    @pytest.mark.parametrize(
        "field,keyword",
        [
            ("sediment_property__concentration", "concentration_initial"),
            ("bedrock_property__concentration", "concentration_in_bedrock"),
        ],
    )
    def test_fields_for_user_array_input(self, field, keyword):
        """Check user input of arrays produces the correct fields"""
        ConcentrationTrackerForSpace(
            self.mg,
            **{
                keyword: (
                    initial_value := np.random.uniform(size=self.mg.number_of_nodes)
                )
            },
        )

        assert np.allclose(self.mg.at_node[field], initial_value)

    @pytest.mark.parametrize(
        "keyword",
        [
            "concentration_initial",
            "concentration_in_bedrock",
            "concentration_from_weathering",
        ],
    )
    def test_properties_concentrations(self, keyword):
        """
        ConcentrationTrackerForSpace should throw an error when input
        concentration values are negative.
        """
        with pytest.raises(ValueError):
            ConcentrationTrackerForSpace(self.mg, **{keyword: -1.0})


class TestAnalytical:
    """Test against analytical solution"""

    def setup_method(self):
        self.mg = RasterModelGrid((3, 5))
        self.mg.axis_units = ("m", "m")
        self.mg.set_status_at_node_on_edges(
            right=NodeStatus.CLOSED,
            top=NodeStatus.CLOSED,
            left=NodeStatus.CLOSED,
            bottom=NodeStatus.CLOSED,
        )
        self.mg.status_at_node[5] = NodeStatus.FIXED_VALUE

        # Grid fields
        self.mg.add_ones("soil__depth", at="node")
        self.mg.add_zeros("bedrock__elevation", at="node")
        self.mg.add_field(
            "topographic__elevation",
            self.mg.at_node["bedrock__elevation"] + self.mg.at_node["soil__depth"],
            at="node",
        )

        # Forced flux rate fields       
        self.mg.add_zeros("sediment__outflux", at="node")
        self.mg.add_zeros("bedrock__erosion_flux", at="node")
        self.mg.add_zeros("sediment__erosion_flux", at="node")
        self.mg.add_zeros("sediment__deposition_flux", at="node")
        
    def test_not_implemented(self):
        """Test that private run_one_step is not implemented"""

        ct = ConcentrationTrackerForSpace(self.mg)
        with pytest.raises(NotImplementedError):
            ct.run_one_step()

### ~~~ TASK ~~~ ###

# UPDATE THESE SO THAT THEY TEST ALL OF THE REQUIRED FLUX FIELDS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    def test_concentration_from_soil_flux(self):
        """
        ConcentrationTrackerForSpace should correctly calculate concentration
        values based on known sediment fluxes.
        """
        # Concentration field for soil
        self.mg.at_node["sediment_property__concentration"] = [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ]

        # Soil flux for middle row of grid is negative 1
        middle_row_link_ids = [9, 10, 11, 12]
        self.mg.at_link["soil__flux"][middle_row_link_ids] -= 1

        # dx is 1 (by default) and soil depth is 1, so soil volume is 1 at each node.
        # Flux of -1 in the middle row should shift all C_sed values left by one node.

        ct = ConcentrationTrackerForDiffusion(self.mg)
        ct.start_tracking()
        ct.stop_tracking(1)

        expected = np.asarray(
            [
                [0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0],
            ]
        ).flatten()
        assert np.allclose(
            self.mg.at_node["sediment_property__concentration"], expected
        )

    def test_concentration_from_weathering_without_conc_w(self):
        """
        ConcentrationTrackerForDiffusion should correctly calculate concentration
        values based on known soil production rates from bedrock weathering when
        concentration_from_weathering input field is not defined.
        """
        # Concentration fields
        self.mg.at_node["sediment_property__concentration"] = [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ]
        self.mg.at_node["bedrock_property__concentration"] = [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ]

        # Soil production rate is 1
        self.mg.at_node["soil_production__rate"].fill(1.0)

        ct = ConcentrationTrackerForDiffusion(self.mg)

        ct.start_tracking()
        # Soil volume is 1 at each node. Soil production rate of 1 doubles volume.
        # This is normally done by the DepthDependentDiffuser. Here, it is forced.
        self.mg.at_node["soil__depth"] += 1.0
        # Node 7: C_sed remains 1 because parent bedrock had conc_br of 1.
        # Node 8: C_sed is halved from 1 to 0.5 because parent bedrock had conc_br = 0.
        ct.stop_tracking(1)

        # Node 7 should have the same concentration as before.
        # Node 8 should have half its previous concentration.
        expected = np.asarray(
            [
                [0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.5, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0],
            ]
        ).flatten()
        assert_allclose(self.mg.at_node["sediment_property__concentration"], expected)

    def test_concentration_from_weathering_with_conc_w(self):
        """
        ConcentrationTrackerForDiffusion should correctly calculate concentration
        values based on known soil production rates from bedrock weathering when
        concentration_from_weathering input field is explicitly defined.
        """
        # Concentration fields
        self.mg.at_node["sediment_property__concentration"] = [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ]
        self.mg.at_node["bedrock_property__concentration"] = [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ]

        # Soil production rate is 1
        self.mg.at_node["soil_production__rate"].fill(1.0)

        ct = ConcentrationTrackerForDiffusion(
            self.mg, concentration_from_weathering=0.0
        )

        ct.start_tracking()
        # Soil volume is 1 at each node. Soil production rate of 1 doubles volume.
        # This is normally done by the DepthDependentDiffuser. Here, it is forced.
        self.mg.at_node["soil__depth"] += 1.0
        # conc_w overrides conc_br values. In this case, no concentration is produced by
        # the weathering process, even at Node 7 where conc_br = 1.

        # Node 7: C_sed is halved from 1 to 0.5 despite parent bedrock with conc_br = 1.
        # Node 8: C_sed is halved from 1 to 0.5 because conc_w = 0.
        ct.stop_tracking(1)

        # Node 7 should have half its previous concentration.
        # Node 8 should have half its previous concentration.
        expected = np.asarray(
            [
                [0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.5, 0.5, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0],
            ]
        ).flatten()
        assert_allclose(self.mg.at_node["sediment_property__concentration"], expected)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

class TestMassBalance:
    """Test mass balance against known scenario"""

    def setup_method(self):
        self.mg = RasterModelGrid((3, 5))
        self.mg.axis_units = ("m", "m")
        self.mg.set_status_at_node_on_edges(
            right=NodeStatus.CLOSED,
            top=NodeStatus.CLOSED,
            left=NodeStatus.CLOSED,
            bottom=NodeStatus.CLOSED,
        )

        self.mg.add_ones("soil__depth", at="node")
        self.mg.add_zeros("bedrock__elevation", at="node")
        self.mg.add_field(
            "topographic__elevation",
            self.mg.at_node["bedrock__elevation"] + self.mg.at_node["soil__depth"],
            at="node",
        )

        # Concentration field for sediment
        self.mg.add_zeros("sediment_property__concentration", at="node")
        self.mg.at_node["sediment_property__concentration"][6] += 1

        # Forced sediment fluxes
        self.mg.add_zeros("sediment__outflux", at="node")
        self.mg.add_zeros("bedrock__erosion_flux", at="node")
        self.mg.add_zeros("sediment__erosion_flux", at="node")
        self.mg.add_zeros("sediment__deposition_flux", at="node")
        # Sediment flux for middle row of grid is negative 1
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ### MAKE SURE FLUXES ARE SET UP CORRECTLY ###
        vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        middle_row_link_ids = [9, 10, 11, 12]
        self.mg.at_link["soil__flux"][middle_row_link_ids] -= 1
        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                    ### MAKE SURE FLUXES ARE SET UP CORRECTLY ###
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    def test_mass_balance_with_all_boundaries_closed(self):
        """
        ConcentrationTrackerForSpace should conserve mass in a scenario with
        all boundaries closed.
        """
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ### MAKE SURE FLUXES ARE SET UP CORRECTLY ###
        vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        # dx is 1 (by default) and soil depth is 1, so soil volume is 1 at each node.
        # Flux of -1 in the middle row should shift all C_sed values left by one node.
        # The boundary is closed, so concentration in Node 6 will remain in the domain.
        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                    ### MAKE SURE FLUXES ARE SET UP CORRECTLY ###
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        total_mass_before = np.sum(
            self.mg.at_node["sediment_property__concentration"]
            * self.mg.at_node["soil__depth"]
            * self.mg.dx**2
        )

        ct = ConcentrationTrackerForSpace(self.mg)
        ct.start_tracking()
        ct.stop_tracking(1)

        total_mass_after = np.sum(
            self.mg.at_node["sediment_property__concentration"]
            * self.mg.at_node["soil__depth"]
            * self.mg.dx**2
        )

        assert_allclose(total_mass_before, total_mass_after)

    def test_mass_balance_with_one_boundary_open(self):
        """
        ConcentrationTrackerForSpace should correctly calculate the mass lost
        from the system over an open boundary.
        """
        self.mg.status_at_node[5] = NodeStatus.FIXED_VALUE

        total_mass_before = np.sum(
            self.mg.at_node["sediment_property__concentration"]
            * self.mg.at_node["soil__depth"]
            * self.mg.dx**2
        )

        ct = ConcentrationTrackerForSpace(self.mg)
        ct.start_tracking()
        ct.stop_tracking(1)

        total_mass_leaving = 1
        total_mass_after = np.sum(
            self.mg.at_node["sediment_property__concentration"]
            * self.mg.at_node["soil__depth"]
            * self.mg.dx**2
        )

        assert_allclose(total_mass_before, total_mass_after + total_mass_leaving)


# %%
class TestFieldCopy:
    """Test that copied field is a copy, but not a reference."""

    def setup_method(self):
        self.mg.add_zeros("soil__depth", at="node")
        self.mg.add_zeros("sediment__outflux", at="node")
        self.mg.add_zeros("bedrock__erosion_flux", at="node")
        self.mg.add_zeros("sediment__erosion_flux", at="node")
        self.mg.add_zeros("sediment__deposition_flux", at="node")
        self.mg.add_zeros("topographic__elevation", at="node")

    def test_copy_is_equal(self):
        """Test that copied values are equal to copied field."""

        ct = ConcentrationTrackerForSpace(self.mg)
        ct._copy_old_soil_depth()

        assert np.allclose(ct._soil__depth_old, self.mg.at_node["soil__depth"])

    def test_copy_is_not_reference(self):
        """Test that copy not a reference."""

        ct = ConcentrationTrackerForSpace(self.mg)
        ct._copy_old_soil_depth()

        self.mg.at_node["soil__depth"] += 1

        assert not np.allclose(ct._soil__depth_old, self.mg.at_node["soil__depth"])
