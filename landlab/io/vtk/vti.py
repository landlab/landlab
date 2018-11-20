#! /bin/env python

from landlab.io.vtk.vtktypes import VtkUniformRectilinear
from landlab.io.vtk.vtkxml import (
    VtkCellDataElement,
    VtkExtent,
    VtkGridElement,
    VtkOrigin,
    VtkPieceElement,
    VtkPointDataElement,
    VtkRootElement,
    VtkSpacing,
)
from landlab.io.vtk.writer import VtkWriter


class VtkUniformRectilinearWriter(VtkWriter):
    _vtk_grid_type = VtkUniformRectilinear

    def construct_field_elements(self, field):
        extent = VtkExtent(field.shape[::-1])
        origin = VtkOrigin(field.origin[::-1], field.spacing[::-1])
        spacing = VtkSpacing(field.spacing[::-1])

        element = {
            "VTKFile": VtkRootElement(VtkUniformRectilinear),
            "Grid": VtkGridElement(
                VtkUniformRectilinear,
                WholeExtent=extent,
                Origin=origin,
                Spacing=spacing,
            ),
            "Piece": VtkPieceElement(Extent=extent),
            "PointData": VtkPointDataElement(
                field.at_node, append=self.data, encoding=self.encoding
            ),
            "CellData": VtkCellDataElement(
                field.at_cell, append=self.data, encoding=self.encoding
            ),
        }

        return element
