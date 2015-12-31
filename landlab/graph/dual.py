class DualGraphMixIn(object):
    @property
    def node_at_cell(self):
        return self._node_at_cell

    @property
    def x_of_corner(self):
        return self._dual.x_of_node

    @property
    def y_of_corner(self):
        return self._dual.y_of_node

    @property
    def corners(self):
        return self._dual.nodes

    @property
    def number_of_corners(self):
        return self._dual.number_of_nodes

    @property
    def corners_at_face(self):
        return self._dual.nodes_at_link

    @property
    def corner_at_face_tail(self):
        return self._dual.node_at_link_tail

    @property
    def corner_at_face_head(self):
        return self._dual.node_at_link_head

    @property
    def number_of_faces(self):
        return self._dual.number_of_links

    @property
    def faces_at_cell(self):
        return self._dual.links_at_patch

    @property
    def corners_at_cell(self):
        return self._dual.nodes_at_patch

    @property
    def number_of_cells(self):
        return self._dual.number_of_patches

    @property
    def faces_at_corner(self):
        return self._dual.links_at_node

    @property
    def face_dirs_at_corner(self):
        return self._dual.link_dirs_at_node
