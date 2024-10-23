"""
Created on Wed Jul 12 12:25:27 2023

@author: LaurentRoberge
"""

import numpy as np
import pytest

from landlab import FieldError
from landlab import NodeStatus
from landlab import RasterModelGrid
from landlab.components import ConcentrationTrackerForSpace


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

    def test_inputs_phi_fraction_fines(self):
        """
        ConcentrationTrackerForSpace should throw an error when phi,
        fraction_fines < 0 or > 1.
        """
        for phi in [-0.2, 1.2]:
            with pytest.raises(ValueError):
                ConcentrationTrackerForSpace(
                    self.mg,
                    phi=phi,
                    fraction_fines=0,
                    settling_velocity=1,
                )
        for fraction_fines in [-0.2, 1.2]:
            with pytest.raises(ValueError):
                ConcentrationTrackerForSpace(
                    self.mg,
                    phi=0,
                    fraction_fines=fraction_fines,
                    settling_velocity=1,
                )

    @pytest.mark.parametrize(
        "at,required_field",
        [
            ("node", "soil__depth"),
            ("node", "sediment__outflux"),
            ("node", "bedrock__erosion_flux"),
            ("node", "sediment__erosion_flux"),
            ("node", "sediment__deposition_flux"),
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
            ConcentrationTrackerForSpace(
                self.mg,
                phi=0,
                fraction_fines=0,
                settling_velocity=1,
            )

    def test_field_instantiation(self):
        """
        ConcentrationTrackerForSpace should instantiate the following fields
        when they do not already exist ('bedrock_property__concentration' and
        'sediment_property__concentration')
        """
        ConcentrationTrackerForSpace(
            self.mg,
            phi=0,
            fraction_fines=0,
            settling_velocity=1,
        )

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
        ConcentrationTrackerForSpace(
            self.mg,
            phi=0,
            fraction_fines=0,
            settling_velocity=1,
        )

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

        ConcentrationTrackerForSpace(
            self.mg,
            phi=0,
            fraction_fines=0,
            settling_velocity=1,
        )

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
            self.mg,
            phi=0,
            fraction_fines=0,
            settling_velocity=1,
            **{keyword: (initial_value := np.random.uniform())},
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
            phi=0,
            fraction_fines=0,
            settling_velocity=1,
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
        ],
    )
    def test_properties_concentrations(self, keyword):
        """
        ConcentrationTrackerForSpace should throw an error when input
        concentration values are negative.
        """
        with pytest.raises(ValueError):
            ConcentrationTrackerForSpace(
                self.mg,
                phi=0,
                fraction_fines=0,
                settling_velocity=1,
                **{keyword: -1.0},
            )


class TestAnalytical:
    """Test against analytical solution"""

    def setup_method(self):
        self.mg = RasterModelGrid((3, 5), xy_spacing=1.0)
        self.mg.axis_units = ("m", "m")
        self.mg.set_status_at_node_on_edges(
            right=NodeStatus.CLOSED,
            top=NodeStatus.CLOSED,
            left=NodeStatus.CLOSED,
            bottom=NodeStatus.CLOSED,
        )
        self.mg.status_at_node[5] = NodeStatus.FIXED_VALUE

        # Grid fields
        self.mg.add_zeros("soil__depth", at="node")
        self.mg.at_node["soil__depth"][:] += 2
        self.mg.add_zeros("bedrock__elevation", at="node")
        self.mg.at_node["bedrock__elevation"][:] += self.mg.node_x
        self.mg.add_field(
            "topographic__elevation",
            self.mg.at_node["bedrock__elevation"] + self.mg.at_node["soil__depth"],
            at="node",
        )
        self.mg.add_zeros("sediment_property__concentration", at="node")
        self.mg.add_zeros("bedrock_property__concentration", at="node")

        # Add forced flow router and flux fields to grid and apply values.
        flow_us_node_order = [5, 6, 7, 8, 3, 2, 1, 0, 4, 9, 10, 11, 12, 13, 14]
        flow_receiver_ids = [
            [0, 1, 2, 3, 4],
            [5, 5, 6, 7, 8],
            [10, 11, 12, 13, 14],
        ]
        discharge = [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ]
        Qs_out = [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ]
        Er = [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ]
        Es = [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ]
        Dsw = [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ]

        self.mg.add_field("flow__receiver_node", flow_receiver_ids, at="node")
        self.mg.add_field("flow__upstream_node_order", flow_us_node_order, at="node")
        self.mg.add_field("surface_water__discharge", discharge, at="node")
        self.mg.add_field("sediment__outflux", Qs_out, at="node")
        self.mg.add_field("bedrock__erosion_flux", Er, at="node")
        self.mg.add_field("sediment__erosion_flux", Es, at="node")
        self.mg.add_field("sediment__deposition_flux", Dsw, at="node")

    def test_not_implemented(self):
        """Test that private run_one_step is not implemented"""

        ct = ConcentrationTrackerForSpace(
            self.mg,
            phi=0,
            fraction_fines=0,
            settling_velocity=1,
        )
        with pytest.raises(NotImplementedError):
            ct.run_one_step()

    @pytest.mark.parametrize(
        "Cs, Cr, expected, expected_sw",
        [
            (0.0, 0.0, 0.0, 0.0),
            (1.0, 0.0, 3.0 / 4.0, 1.0 / 2.0),
            (0.0, 1.0, 1.0 / 4.0, 1.0 / 2.0),
            (1.0, 1.0, 1.0, 1.0),
        ],
    )
    def test_Es_Er_Dsw(self, Cs, Cr, expected, expected_sw):
        self.mg.at_node["sediment_property__concentration"][:] = Cs
        self.mg.at_node["bedrock_property__concentration"][:] = Cr

        ct = ConcentrationTrackerForSpace(
            self.mg,
            phi=0,
            fraction_fines=0,
            settling_velocity=1,
        )

        ct.start_tracking()
        self.mg.at_node["soil__depth"][:] -= self.mg.at_node["sediment__erosion_flux"][
            :
        ]
        self.mg.at_node["soil__depth"][:] += self.mg.at_node[
            "sediment__deposition_flux"
        ][:]
        self.mg.at_node["topographic__elevation"][:] = (
            self.mg.at_node["bedrock__elevation"][:] + self.mg.at_node["soil__depth"][:]
        )
        ct.stop_tracking(1)

        assert np.allclose(ct._C_sw[8], expected_sw)
        assert np.allclose(
            self.mg.at_node["sediment_property__concentration"][8], expected
        )


# %%
class TestFieldCopy:
    """Test that copied field is a copy, but not a reference."""

    def setup_method(self):
        self.mg = RasterModelGrid((3, 3))
        self.mg.add_zeros("soil__depth", at="node")
        self.mg.add_zeros("sediment__outflux", at="node")
        self.mg.add_zeros("bedrock__erosion_flux", at="node")
        self.mg.add_zeros("sediment__erosion_flux", at="node")
        self.mg.add_zeros("sediment__deposition_flux", at="node")
        self.mg.add_zeros("topographic__elevation", at="node")

    def test_copy_is_equal(self):
        """Test that copied values are equal to copied field."""

        ct = ConcentrationTrackerForSpace(
            self.mg,
            phi=0,
            fraction_fines=0,
            settling_velocity=1,
        )
        ct._copy_old_soil_depth()

        assert np.allclose(ct._soil__depth_old, self.mg.at_node["soil__depth"])

    def test_copy_is_not_reference(self):
        """Test that copy not a reference."""

        ct = ConcentrationTrackerForSpace(
            self.mg,
            phi=0,
            fraction_fines=0,
            settling_velocity=1,
        )
        ct._copy_old_soil_depth()

        self.mg.at_node["soil__depth"] += 1

        assert not np.allclose(ct._soil__depth_old, self.mg.at_node["soil__depth"])
