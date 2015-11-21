#! /usr/bin/env python
from .grouped import ModelDataFields


class ModelDataFieldsMixIn(ModelDataFields):
    def empty(self, *args, **kwds):
        if len(args) == 0:
            group = kwds.pop('centering', 'node')
        else:
            group = args[0]

        if self[group].size is None:
            self[group].size = self.number_of_elements(group)
        return ModelDataFields.empty(self, group, **kwds)

    def ones(self, *args, **kwds):
        if len(args) == 0:
            group = kwds.pop('centering', 'node')
        else:
            group = args[0]

        if self[group].size is None:
            self[group].size = self.number_of_elements(group)
        return ModelDataFields.ones(self, group, **kwds)

    def zeros(self, *args, **kwds):
        if len(args) == 0:
            group = kwds.pop('centering', 'node')
        else:
            group = args[0]

        if self[group].size is None:
            self[group].size = self.number_of_elements(group)
        return ModelDataFields.zeros(self, group, **kwds)
