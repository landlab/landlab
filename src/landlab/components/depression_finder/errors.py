class LakeMapperError(RuntimeError):
    pass


class NoOutletError(LakeMapperError):
    def __init__(self, node, lake_size):
        self._node = node
        self._size = lake_size

    def __str__(self):
        return (
            "Unable to find a drainage outlet for a lake."
            " This may be because you have disconnected open nodes, which"
            " sometimes occurs during raster clipping. Consider running"
            " `set_open_nodes_disconnected_from_watershed_to_closed`,"
            " which will remove any isolated open nodes."
            f" Began searching for an outlet at node {self._node} and have"
            f" found {self._size} node{'s' if self._size > 1 else ''} in"
            " the current lake."
        )
