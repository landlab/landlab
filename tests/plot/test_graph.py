import pytest

from landlab import RasterModelGrid
from landlab.plot.graph import _parse_locations_as_set
from landlab.plot.graph import plot_graph


def _axes_arrows(ax):
    from matplotlib.patches import FancyArrow

    return [child for child in ax.get_children() if isinstance(child, FancyArrow)]


def test_parse_locations_from_str():
    assert _parse_locations_as_set("node") == {"node"}
    assert _parse_locations_as_set("node,link") == {"node", "link"}
    assert _parse_locations_as_set("node,link,cell,link") == {"node", "link", "cell"}
    assert _parse_locations_as_set("node  , link, patch  ") == {"node", "link", "patch"}


def test_parse_locations_from_iterable():
    assert _parse_locations_as_set(["cell"]) == {"cell"}
    assert _parse_locations_as_set(("patch", "corner")) == {"patch", "corner"}
    assert _parse_locations_as_set(("patch", "corner", "patch")) == {"patch", "corner"}
    assert _parse_locations_as_set((" patch  ", "corner  ")) == {"patch", "corner"}


def test_parse_locations_bad_value():
    with pytest.raises(ValueError, match=r"^unknown location "):
        _parse_locations_as_set("foo")
    with pytest.raises(ValueError, match=r"^unknown locations "):
        _parse_locations_as_set("cells,nodes")


@pytest.mark.parametrize("at", ["node", "link", "patch", "corner", "face", "cell"])
def test_plot_graph(at):
    grid = RasterModelGrid((3, 4))
    ax = plot_graph(grid, at=at)

    assert ax.get_ylim() == (-0.5, 2.5)
    assert ax.get_xlim() == (-0.5, 3.5)

    if at in ("patch", "cell"):
        assert len(ax.patches) == grid.number_of_elements(at)
    elif at in ("link", "face"):
        assert len(_axes_arrows(ax)) == grid.number_of_elements(at)
    else:
        assert len(ax.lines) == grid.number_of_elements(at)


def test_plot_graph_onto_existing():
    grid = RasterModelGrid((3, 4))
    ax = plot_graph(grid, at="node")
    plot_graph(grid, at="cell", axes=ax)

    assert (
        len(ax.patches) + len(ax.lines) == grid.number_of_nodes + grid.number_of_cells
    )
