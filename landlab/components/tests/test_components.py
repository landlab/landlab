import pytest

from landlab import RasterModelGrid
from landlab.components import COMPONENTS

_REQ_ATTRS = ["doc", "mapping", "type", "intent", "optional", "units"]


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
                if at == "grid":
                    grid.at_grid[name] = 0
                else:
                    grid.add_zeros(at, name)

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

        # verify all info exist:
        for name in Comp._info.keys():
            info = Comp._info[name].copy()

            for attribute in _REQ_ATTRS:
                if attribute in info:
                    info.pop(attribute)
                else:
                    raise ValueError(
                        "{component} is missing attribute {attribute} about variable: {name} at {at}".format(
                            component=Comp._name, name=name, at=at, attribute=attribute
                        )
                    )

            # TODO: Verify that all units are UDUNITS compatible.

            # TODO: Verify that all dtypes are valid.


            # TODO: Verify that all mappings are valid grid locations.




            if len(info) > 0:
                raise ValueError(
                    "{component} has an extra attribute {attribute} about variable: {name} at {at}".format(
                        component=Comp._name, name=name, at=at, attribute=attribute
                    )
                )
