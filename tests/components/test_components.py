import pandas as pd
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
        "NetworkSedimentTransporter",
        "Profiler",
        "SoilMoisture",
        "Vegetation",
    ):
        grid = RasterModelGrid((10, 10))

        assert Comp._unit_agnostic in (True, False)

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
                        component=Comp._name, name=name, at=at,
                    )
                )


def test_consistent_doc_names():
    out = []
    for comp in COMPONENTS:
        for name in comp._info:
            temp = {"component": comp.__name__, "field": name}
            for key in comp._info[name].keys():
                temp[key] = comp._info[name][key]
            out.append(temp)
    df = pd.DataFrame(out)

    unique_fields = df.field.unique().astype(str)

    bad_fields = {}
    for field in unique_fields:
        where = df.field == field
        if where.sum() > 1:

            sel = df[where]

            doc_vals = df.doc[where].values.astype(str)

            inconsistent = []
            for i in range(len(doc_vals) - 1):
                if doc_vals[i] != doc_vals[-1]:
                    inconsistent.append(sel.component.values[i])

            if len(inconsistent) > 0:
                bad_fields[field] = inconsistent
    if len(bad_fields) > 0:
        msg = "The following fields have inconsistent documentation:\n"
        for field in bad_fields.keys():
            inconsistent = bad_fields[field]
            msg += "\n" + field + ":\n  " + "\n  ".join(inconsistent)

        raise ValueError(msg)
