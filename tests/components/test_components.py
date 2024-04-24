import numpy as np
import pandas as pd
import pytest

from landlab import FieldError
from landlab import RasterModelGrid
from landlab.components import COMPONENTS
from landlab.components import PriorityFloodFlowRouter

_VALID_LOCS = {"grid", "node", "link", "patch", "corner", "face", "cell"}

_REQUIRED_ATTRS = {"doc", "mapping", "dtype", "intent", "optional", "units"}

_EXCLUDE_COMPONENTS = {
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
    "BedParcelInitializerDischarge",
    "BedParcelInitializerDepth",
    "BedParcelInitializerArea",
    "BedParcelInitializerUserD50",
    "SedimentPulserEachParcel",
    "SedimentPulserAtLinks",
}


try:
    PriorityFloodFlowRouter.load_richdem()
except ModuleNotFoundError:
    _EXCLUDE_COMPONENTS.add("PriorityFloodFlowRouter")


@pytest.mark.parametrize("Comp", COMPONENTS)
def test_component_info_unit_agnostic(Comp):
    """Check for a valid _units_agnostic attribute"""
    assert Comp._unit_agnostic in (True, False)


def _add_input_fields_to_grid(cls, grid):
    for name, meta in cls._info.items():
        if meta["intent"].startswith("in"):
            at = cls.var_loc(name)
            dtype = cls.var_type(name)
            if at == "grid":
                grid.at_grid[name] = np.array(0, dtype=dtype)
            else:
                grid.add_zeros(name, at=at, dtype=dtype)

    return grid


@pytest.mark.parametrize("Comp", COMPONENTS)
def test_component_output_fields(Comp):
    """Check that required output fields exist with correct dtypes and locations"""
    if Comp.name in _EXCLUDE_COMPONENTS:
        pytest.skip("component explicitly excluded")

    component_name = Comp._name
    grid = RasterModelGrid((10, 10))

    _add_input_fields_to_grid(Comp, grid)
    Comp(grid)

    for name, meta in Comp._info.items():
        if meta["intent"].endswith("out") and not meta["optional"]:
            at = meta["mapping"]
            if name not in grid[at]:
                raise ValueError(
                    f"{component_name} is missing output variable: {name} at {at}"
                )

            expected_dtype = meta["dtype"]
            actual_dtype = grid[at][name].dtype

            if actual_dtype != expected_dtype:
                raise FieldError(
                    f"{component_name} output required variable: {name} at {at} has "
                    f"incorrect dtype. dtype must be {expected_dtype} and is "
                    f"{actual_dtype}"
                )


@pytest.mark.parametrize("Comp", COMPONENTS)
def test_component_info_missing_attrs(Comp):
    """Check that in/out fields are not missing attributes"""
    component_name = Comp._name

    for name, meta in Comp._info.items():
        at = meta["mapping"]

        missing = ", ".join(sorted(_REQUIRED_ATTRS - set(meta)))
        if missing:
            raise ValueError(
                f"{component_name} is missing attributes ({missing}) about variable: "
                f"{name} at {at}"
            )


@pytest.mark.parametrize("Comp", COMPONENTS)
def test_component_info_unknown_attrs(Comp):
    """Check that in/out fields have valid attributes"""
    component_name = Comp._name

    for name, meta in Comp._info.items():
        at = meta["mapping"]

        unknown = ", ".join(sorted(set(meta) - _REQUIRED_ATTRS))
        if unknown:
            raise ValueError(
                f"{component_name} has extra attributes ({unknown}) about variable: "
                f"{name} at {at}"
            )


@pytest.mark.parametrize("Comp", COMPONENTS)
def test_component_info_valid_dtype(Comp):
    """Check that fields have a valid numpy dtype"""
    component_name = Comp._name

    for name, meta in Comp._info.items():
        dtype = meta["dtype"]
        try:
            np.dtype(dtype)
        except TypeError as exc:
            raise ValueError(
                f"{component_name} has a bad dtype ({dtype}) for variable: {name}"
            ) from exc


@pytest.mark.parametrize("Comp", COMPONENTS)
def test_component_info_valid_locations(Comp):
    """Check that fields are defined at valid locations"""
    component_name = Comp._name

    # verify all info exist:
    for name, meta in Comp._info.items():
        at = meta["mapping"]

        # TODO: Verify that all units are UDUNITS compatible.

        if at not in _VALID_LOCS:
            raise ValueError(
                f"{component_name} mapping for variable: {name} is invalid: {at}"
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
