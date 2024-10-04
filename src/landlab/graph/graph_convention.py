import re
from collections import OrderedDict


class GraphConvention:
    """Define a naming convention for graph elements."""

    def __init__(self, node, edge, face, nodes=None, edges=None, faces=None):
        """Define a graph element naming convention.

        Parameters
        ----------
        node : str
            Name to use for "nodes".
        edge : str
            Name to use for "edges".
        face : str
            Name to use for "faces".
        nodes : str, optional
            Plural of "nodes"
        edges : str, optional
            Plural of "edges"
        faces : str, optional
            Plural of "faces"

        Examples
        --------
        >>> from landlab.graph import GraphConvention
        >>> convention = GraphConvention("node", "edge", "face")
        >>> convention.node
        'node'
        >>> convention.edge
        'edge'
        >>> convention.edges
        'edges'

        >>> convention = GraphConvention("node", "link", "patch", faces="patches")
        >>> convention.face
        'patch'
        >>> convention.faces
        'patches'
        """
        self._node = node
        self._edge = edge
        self._face = face
        self._nodes = nodes or node + "s"
        self._edges = edges or edge + "s"
        self._faces = faces or face + "s"

    @property
    def node(self):
        """The name of things that are points."""
        return self._node

    @property
    def nodes(self):
        """The plural name for node."""
        return self._nodes

    @property
    def edge(self):
        """The name of lines that connect two nodes."""
        return self._edge

    @property
    def edges(self):
        """The plural name of edge."""
        return self._edges

    @property
    def face(self):
        """The name of objects made up of a closed set of edges."""
        return self._face

    @property
    def faces(self):
        """The plural name for face."""
        return self._faces


class ConventionConverter:
    """Convert between graph element naming conventions."""

    CONVENTION = {
        "nlp": GraphConvention("node", "link", "patch", faces="patches"),
        "nef": GraphConvention("node", "edge", "face"),
        "cfc": GraphConvention("corner", "face", "cell"),
    }

    def __init__(self, convention="nlp"):
        """

        Parameters
        ----------
        convention : {"nlp", "nef", "cfc"}
            Naming convention to which to conform.
        """
        self.convention = convention

    @property
    def convention(self):
        """Name of the convention to conform to."""
        return self._convention

    @convention.setter
    def convention(self, val):
        if val not in ("nlp", "cfc", "nef"):
            raise ValueError(f"convention not understood ({val})")
        self._convention = val

    @staticmethod
    def _convention_mapper(from_convention, to_convention):
        """Create a mapper to convert names between conventions."""
        elements = ("nodes", "edges", "faces", "node", "edge", "face")

        from_elements = [getattr(from_convention, name) for name in elements]
        to_elements = [getattr(to_convention, name) for name in elements]

        return OrderedDict(
            list(zip(to_elements, from_elements))
            + list(zip(from_elements, to_elements))
        )

    @classmethod
    def _as_convention(cls, name):
        """Get a named convention as a GraphConvention object."""
        try:
            return cls.CONVENTION[name]
        except KeyError as exc:
            raise ValueError(f"convention not understood ({name})") from exc

    def conform(self, name, from_convention):
        """Convert a name to a new convention.

        Parameters
        ----------
        name : str
            A name in the given convention.
        from_convention : {"nlp", "nef", "cfc"}
            Naming convention of the provided string.

        Returns
        -------
        str
            The name converted to the convention.

        Examples
        --------
        >>> from landlab.graph import ConventionConverter

        >>> converter = ConventionConverter("nlp")

        >> converter.conform("nodes_at_edge", "nef")
        'nodes_at_link'

        >>> converter.conform("node_at_cell", "cfc")
        'corner_at_patch'

        >>> converter = ConventionConverter("nef")
        >>> converter.conform("faces_at_cell", "cfc")
        'edges_at_face'
        """
        if from_convention is self.convention:
            return name

        mapping = ConventionConverter._convention_mapper(
            ConventionConverter._as_convention(from_convention),
            ConventionConverter._as_convention(self.convention),
        )

        def rename_words(m):
            key = m.group(0)
            return mapping[key]

        pattern = "|".join(mapping.keys())
        return re.sub(pattern, rename_words, name)
