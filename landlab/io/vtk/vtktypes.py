#! /bin/env python

from landlab.io.vtk.encoders import AsciiEncoder, Base64Encoder, RawEncoder


class VtkEndian(object):
    def __init__(self, endian):
        self.name = endian

    def __str__(self):
        return self.name

    def __repr__(self):
        return "VtkEndian(%s)" % self.name


class VtkType(object):
    def __init__(self, name, size):
        self.bytes = size
        self.name = name

    def __str__(self):
        return self.name

    def __repr__(self):
        return "VtkType(%s, %s)" % (self.name, self.bytes)


class VtkCellType(object):
    def __init__(self, name, id):
        self.id = id
        self.name = name

    def __str__(self):
        return self.name

    def __int__(self):
        return self.id

    def __repr__(self):
        return "VtkCellType(%s, %s)" % (self.name, self.id)


class VtkGridType(object):
    def __init__(self, name):
        self.name = name

    def __str__(self):
        return self.name

    def __repr__(self):
        return "VtkGridType(%s)" % self.name


VtkLittleEndian = VtkEndian("LittleEndian")
VtkBigEndian = VtkEndian("BigEndian")


SYS_TO_VTK_ENDIAN = {"little": VtkLittleEndian, "big": VtkBigEndian}


VtkUInt8 = VtkType("UInt8", 1)
VtkInt32 = VtkType("Int32", 4)
VtkInt64 = VtkType("Int64", 8)
VtkFloat32 = VtkType("Float32", 4)
VtkFloat64 = VtkType("Float64", 8)


NUMPY_TO_VTK_TYPE = {
    "uint8": VtkUInt8,
    "int32": VtkInt32,
    "int64": VtkInt64,
    "float32": VtkFloat32,
    "float64": VtkFloat64,
}

VTK_TO_NUMPY_TYPE = {
    "UInt8": "uint8",
    "Int32": "int32",
    "Int64": "int64",
    "Float32": "float32",
    "Float64": "float64",
}


VtkVertex = VtkCellType("Vertex", 1)
VtkLine = VtkCellType("Line", 3)
VtkTriangle = VtkCellType("Triangle", 5)
VtkQuad = VtkCellType("Quad", 9)
VtkPolygon = VtkCellType("Polygon", 7)


EDGE_COUNT_TO_TYPE = {1: VtkVertex, 2: VtkLine, 3: VtkTriangle, 4: VtkQuad}


VtkUniformRectilinear = VtkGridType("ImageData")
VtkRectilinear = VtkGridType("RectilinearGrid")
VtkStructured = VtkGridType("StructuredGrid")
VtkUnstructured = VtkGridType("UnstructuredGrid")


ENCODERS = {"ascii": AsciiEncoder(), "raw": RawEncoder(), "base64": Base64Encoder()}
