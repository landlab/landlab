from __future__ import annotations

from collections.abc import Iterable
from collections.abc import Mapping
from collections.abc import Sequence
from enum import Enum
from enum import auto

import numpy as np
from numpy.typing import ArrayLike
from numpy.typing import NDArray

from landlab.graph.structured_quad import DualStructuredQuadGraph
from landlab.grid.base import ModelGrid
from landlab.grid.nodestatus import NodeStatus


class QuadModelGrid(DualStructuredQuadGraph, ModelGrid):
    """Structured quadrilateral grid with arbitrary node spacing.

    This grid represents a logically rectangular array of quadrilateral cells
    whose node coordinates may be uniformly spaced, nonuniformly spaced, or
    fully general (e.g., sheared or curvilinear).

    The grid topology is structured, but geometry is inferred from node
    coordinates.

    Examples
    --------

    >>> from landlab.grid.quad import QuadModelGrid
    >>> grid = QuadModelGrid.from_raster((4, 5))
    >>> grid.shape
    (4, 5)
    >>> len(grid.active_links)
    17

    Set the nodes along the top edge of the grid to be *closed* boundaries.
    This means that any links touching these nodes will be *inactive*.

    >>> grid = QuadModelGrid.from_raster((4, 5), bc={"top": "closed"})
    >>> grid.shape
    (4, 5)
    >>> len(grid.active_links)
    14

    A `RasterModelGrid` can have different node spacings in the *x* and *y*
    directions.

    >>> grid = QuadModelGrid.from_raster((4, 5), xy_spacing=(2, 1))
    >>> grid.xy_spacing
    (2.0, 1.0)
    >>> grid.y_of_node.reshape(grid.shape)
    array([[0.,  0.,  0.,  0.,  0.],
           [1.,  1.,  1.,  1.,  1.],
           [2.,  2.,  2.,  2.,  2.],
           [3.,  3.,  3.,  3.,  3.]])
    >>> grid.x_of_node.reshape(grid.shape)
    array([[0.,  2.,  4.,  6.,  8.],
           [0.,  2.,  4.,  6.,  8.],
           [0.,  2.,  4.,  6.,  8.],
           [0.,  2.,  4.,  6.,  8.]])
    """

    def __init__(
        self,
        x_of_node: ArrayLike,
        y_of_node: ArrayLike,
        *,
        xy_axis_name: Sequence[str] = ("x", "y"),
        xy_axis_units: Sequence[str] = ("m", "m"),
        xy_of_reference: Sequence[float] = (0.0, 0.0),
        bc: str | dict[str, str] = "open",
    ) -> None:
        """Create a quadrilateral grid from 2D node coordinate arrays.

        Parameters
        ----------
        x_of_node, y_of_node : array-like of shape (n_rows, n_cols)
            Node coordinates indexed as ``[row, col]``.
        xy_axis_name : sequence of str, optional
            Names of the spatial axes.
        xy_axis_units : sequence of str, optional
            Units of the spatial axes.
        xy_of_reference : sequence of float, optional
            Reference coordinate for the grid origin.
        bc : str or dict[str, str], optional
            Edge boundary conditions.

        Raises
        ------
        ValueError
            If coordinate arrays have different shapes or are not 2D.
        """
        x_of_node = np.asarray(x_of_node)
        y_of_node = np.asarray(y_of_node)

        if x_of_node.shape != y_of_node.shape:
            raise ValueError("coordinate arrays must have the same shape")
        if x_of_node.ndim != 2:
            raise ValueError("coordinate arrays must be 2D")

        shape: tuple[int, int] = tuple(x_of_node.shape)

        DualStructuredQuadGraph.__init__(
            self,
            (y_of_node.reshape(-1), x_of_node.reshape(-1)),
            shape=shape,
            sort=False,
        )

        ModelGrid.__init__(
            self,
            xy_axis_name=tuple(xy_axis_name),
            xy_axis_units=tuple(xy_axis_units),
            xy_of_reference=tuple(xy_of_reference),
        )

        self._node_status = _make_node_status(
            self.shape, bc_at_edge=_normalize_bc_at_edge(bc)
        ).reshape(-1)

        self._geometry = QuadGridGeometry.UNKNOWN
        self._xy_spacing: tuple[float, float] | None = None

    @classmethod
    def from_rectilinear(cls, x: ArrayLike, y: ArrayLike, **kwds) -> QuadModelGrid:
        """Create a rectilinear grid from 1D coordinate arrays.

        Parameters
        ----------
        x, y : array-like of shape (n_cols,) and (n_rows,)
            Coordinates defining column-wise x positions and row-wise y positions.
        **kwds
            Additional keyword arguments passed to the grid constructor.

        Returns
        -------
        QuadModelGrid
            A rectilinear quadrilateral grid with possibly nonuniform spacing.

        Raises
        ------
        ValueError
            If ``x`` or ``y`` is not one-dimensional.
        """
        if x.ndim != 1 or y.ndim != 1:
            raise ValueError("coordinate arrays must be 1D")

        x_of_node, y_of_node = np.meshgrid(x, y, indexing="xy")
        grid = cls(x_of_node, y_of_node, **kwds)
        grid._geometry = QuadGridGeometry.RECTILINEAR

        return grid

    @classmethod
    def from_raster(
        cls, shape: tuple[int, int], xy_spacing: float | Iterable[float] = 1.0, **kwds
    ) -> QuadModelGrid:
        """Create a raster grid with uniform spacing.

        Parameters
        ----------
        shape : tuple of int
            Number of rows and columns in the grid.
        xy_spacing : float or iterable of float, optional
            Uniform spacing in the x and y directions.
        **kwds
            Additional keyword arguments passed to the grid constructor.

        Returns
        -------
        QuadModelGrid
            A raster-equivalent quadrilateral grid.
        """
        xy_spacing = np.broadcast_to(np.asarray(xy_spacing, dtype=float), (2,))

        x_of_node, y_of_node = np.meshgrid(
            np.arange(shape[1]) * xy_spacing[0],
            np.arange(shape[0]) * xy_spacing[1],
            indexing="xy",
        )
        grid = cls(x_of_node, y_of_node, **kwds)

        grid._xy_spacing = tuple(map(float, xy_spacing))
        grid._geometry = QuadGridGeometry.RASTER

        return grid

    @property
    def geometry(self) -> QuadGridGeometry:
        """Geometric classification of the grid.

        Returns
        -------
        QuadGridGeometry
            The grid's geometry classification.

        Notes
        -----
        The geometry is determined lazily. If not explicitly set during
        construction, it is inferred from node coordinates on first access
        and then cached.
        """
        if self._geometry is QuadGridGeometry.UNKNOWN:
            x = self.x_of_node.reshape(self.shape)
            y = self.y_of_node.reshape(self.shape)

            self._geometry = _classify_quad_grid_geometry(x, y)

            if self._geometry is QuadGridGeometry.RASTER and self._xy_spacing is None:
                self._xy_spacing = (float(x[0, 1] - x[0, 0]), float(y[1, 0] - y[0, 0]))

        return self._geometry

    @property
    def xy_spacing(self) -> tuple[float, float]:
        """Uniform grid spacing in the x and y directions.

        Returns
        -------
        tuple of float
            Spacing in the x and y directions.

        Raises
        ------
        AttributeError
            If the grid is not raster.
        """
        if self._geometry is not QuadGridGeometry.RASTER:
            raise AttributeError(
                "xy_spacing is only defined for raster-equivalent grids"
            )

        return self._xy_spacing


class QuadGridGeometry(Enum):
    """Geometric classifications for structured quadrilateral grids.

    Attributes
    ----------
    UNKNOWN
        Geometry has not yet been classified.
    RASTER
        Rectilinear grid with uniform spacing in both *x* and *y* directions.
    RECTILINEAR
        Axis-aligned grid where *x* varies only with column and *y* varies only
        with row, but spacing may be nonuniform.
    GENERAL
        General quadrilateral grid that is not rectilinear (e.g., sheared
        or curvilinear).
    """

    UNKNOWN = auto()
    RASTER = auto()
    RECTILINEAR = auto()
    GENERAL = auto()


def _classify_quad_grid_geometry(
    x: ArrayLike,
    y: ArrayLike,
    *,
    rtol: float = 1e-12,
    atol: float = 1e-15,
) -> QuadGridGeometry:
    """Classify the geometry of a structured quadrilateral grid from node coordinates.

    A grid is classified as:

    * ``RASTER`` if it is rectilinear with uniform spacing in both directions.
    * ``RECTILINEAR`` if it is axis-aligned but spacing is nonuniform.
    * ``GENERAL`` if either coordinate varies with both row and column.

    Parameters
    ----------
    x, y : array-like of shape (n_rows, n_cols)
        Node coordinates indexed as ``[row, col]``.
    rtol, atol : float, optional
        Relative and absolute tolerances used when comparing floating-point
        values.

    Returns
    -------
    QuadGridGeometry
        The geometric classification of the grid.
    """
    x, y = np.asarray(x), np.asarray(y)

    is_x_rect = np.allclose(x, x[0:1, :], rtol=rtol, atol=atol)
    is_y_rect = np.allclose(y, y[:, 0:1], rtol=rtol, atol=atol)

    if not (is_x_rect and is_y_rect):
        return QuadGridGeometry.GENERAL

    dx = np.diff(x[0, :])
    dy = np.diff(y[:, 0])

    if (
        dx.size > 0
        and dy.size > 0
        and np.allclose(dx, dx[0], rtol=rtol, atol=atol)
        and np.allclose(dy, dy[0], rtol=rtol, atol=atol)
    ):
        return QuadGridGeometry.RASTER

    return QuadGridGeometry.RECTILINEAR


def _make_node_status(
    shape: tuple[int, int],
    bc_at_edge: Mapping[str, NodeStatus],
) -> NDArray:
    """Create a node-status array for a structured grid.

    Parameters
    ----------
    shape : tuple of int
        Grid shape as (n_rows, n_cols).
    bc : mapping[str, NodeStatus], optional
        Boundary condition specification. Keys must be edge names
        ('right', 'top', 'left', 'bottom') and values must be a ``NodeStatus``.

    Returns
    -------
    ndarray
        Array of node statuses with shape ``shape``.
    """
    node_status = np.full(shape, NodeStatus.CORE, dtype=np.uint8)

    for edge, status in bc_at_edge.items():
        match edge:
            case "left":
                node_status[:, 0] = status
            case "right":
                node_status[:, -1] = status
            case "bottom":
                node_status[0, :] = status
            case "top":
                node_status[-1, :] = status
            case _:
                raise KeyError("unknown value for edge")

    return node_status


def _normalize_bc_at_edge(
    bc: str | Mapping[str, str] | NodeStatus | Mapping[str, NodeStatus],
    default_status: NodeStatus = NodeStatus.FIXED_VALUE,
):
    cond_to_status = {
        "open": NodeStatus.FIXED_VALUE,
        "closed": NodeStatus.CLOSED,
    }
    edges = ("right", "top", "left", "bottom")

    def _to_status(bc: str | NodeStatus) -> NodeStatus:
        if isinstance(bc, str):
            return cond_to_status[bc]
        return bc

    if isinstance(bc, str):
        return dict.fromkeys(edges, _to_status(bc))

    if isinstance(bc, NodeStatus):
        return dict.fromkeys(edges, bc)

    if isinstance(bc, Mapping):
        bc_at_edge = {edge: _to_status(status) for edge, status in bc.items()}
        return dict.fromkeys(edges, _to_status(default_status)) | bc_at_edge

    raise TypeError("bc must be either a str, NodeStatus, or a mapping")
