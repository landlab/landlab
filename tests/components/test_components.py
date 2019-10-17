import pytest

from landlab import FieldError, RasterModelGrid
from landlab.components import COMPONENTS

_VALID_LOCS = ("grid", "node", "link", "patch", "corner", "face", "cell")


_REQ_ATTRS = ["doc", "mapping", "dtype", "intent", "optional", "units"]


@pytest.mark.parametrize("Comp", COMPONENTS)
def test_component_metadata(Comp):
    if Comp.name not in (
        "ChannelProfiler",
        "DrainageDensity",
        "gFlex",
        "HackCalculator",
        "Lithology",
        "LithoLayers",
        "SoilMoisture",
        "Vegetation",
    ):
        print(Comp.name)
        grid = RasterModelGrid((10, 10))

        # verify that we can create it
        for name in Comp._info.keys():
            if "in" in Comp._info[name]["intent"]:
                at = Comp.var_loc(name)
                dtype = Comp.var_type(name)
                if at == "grid":
                    grid.at_grid[name] = 0
                else:
                    grid.add_zeros(at, name, dtype=dtype)

        _ = Comp(grid)

        # verify that all output fields are made
        for name in Comp._info.keys():
            if "out" in Comp._info[name]["intent"]:
                if not Comp._info[name]["optional"]:
                    at = Comp._info[name]["mapping"]
                    if name not in grid[at]:
                        raise ValueError(
                            "{component} is missing output variable: {name} at {at}".format(
                                component=Comp._name, name=name, at=at
                            )
                        )

                    field = grid[at][name]
                    dtype = Comp._info[name]["dtype"]

                    try:
                        assert field.dtype == dtype
                    except AssertionError:
                        raise FieldError(
                            "{component} output required variable: {name} at {at} has incorrect dtype. dtype must be {dtype} and is {actual}".format(
                                component=Comp._name,
                                name=name,
                                at=at,
                                dtype=dtype,
                                actual=field.dtype,
                            )
                        )

        # verify all info exist:
        for name in Comp._info.keys():
            info = Comp._info[name].copy()
            at = Comp._info[name]["mapping"]
            for attribute in _REQ_ATTRS:
                if attribute in info:
                    info.pop(attribute)
                else:
                    raise ValueError(
                        "{component} is missing attribute {attribute} about variable: {name} at {at}".format(
                            component=Comp._name, name=name, at=at, attribute=attribute
                        )
                    )

            if len(info) > 0:
                raise ValueError(
                    "{component} has an extra attribute {attribute} about variable: {name} at {at}".format(
                        component=Comp._name, name=name, at=at, attribute=attribute
                    )
                )

            # TODO: Verify that all units are UDUNITS compatible.

            # TODO: Verify that all dtypes are valid.

            # TODO: Verify that all mappings are valid grid locations.
            if Comp._info[name]["mapping"] not in _VALID_LOCS:
                raise ValueError(
                    "{component} mapping for variable: {name} is invalid: {at}".format(
                        component=Comp._name, name=name, at=at, attribute=attribute
                    )
                )
