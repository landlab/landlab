#!/usr/bin/env python
"""ComponentMaker is an example of how to make a component.

Congratulations on wanting to make a Landlab component! This file is part of a
set of files meant to make it easier for you to do this.

See also the `Develop your own component <https://github.com/landlab/landlab/wiki/Develop-your-own-component>`_
page on the Landlab Wiki and the `Creating a Component Tutorial <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/master/making_components/making_components.ipynb>`_.

The Creating a Component Tutorial has information about many many parts of your
component. This file is not meant to replace this, but to provide a concrete
working example that includes all of the parts that we expect.

This text is called a docstring and it provides documentation about your
component. You'll use it describe parameters, methods, and attributes that your
component has.

The name of the class should be in CamelCase, and should make sense when used in
the sentence: "A (component-name) is a...".

This docstring is written in ReStructuredText and will be parsed by `Sphinx <http://www.sphinx-doc.org/en/master/index.html>`_ to
make the API documentation. `Here is information <http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_
about how ReStructuredText works.

Recommendations on coding style
-------------------------------
- Please stick to the coding style described by `PEP8 <https://www.python.org/dev/peps/pep-0008/>`_.
- Also read `PEP 257 <https://www.python.org/dev/peps/pep-0257/>`_ which
describes general python docstring conventions.
- Class and function docstrings should follow the `numpydoc conventions <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_.
- See `Standard Naming conventions <https://github.com/landlab/landlab/wiki/Standard-names>`_
for good practice on parameters and variables naming.
- If you have any questions about the process after you have finished reading the
documentation, consider making an `Issue
<https://github.com/landlab/landlab/issues/new/>`_ to ask the development team
for help.

Examples
--------
The next line is a line of code! after three >>>s and a space are lines that
will be evaluated in tests. It is also helpful for users. Use this space to give
an overall example of using your component.

>>> from landlab import RasterModelGrid
>>> from landlab.components import ComponentMaker
>>> mg = RasterModelGrid(3,3)
>>> cm = ComponentMaker(mg, spam=True, eggs=1.0)
>>> cm.spam
True
>>> cm.eggs
1.0
>>> 'component_maker__field' in mg.at_node
True
>>> cm.run_one_step(dt=100)
"""


import numpy as np
from landlab import Component


class ComponentMaker(Component):
    """ComponentMaker is an example of how to make a component.

    If you had more to write, you could write it here. After writing the text,
    we can make two special sections 'Attributes' and 'Methods' in which we list
    the attributes and methods specific to our component. Magically, the first
    line of the docstring for each of these elements will be reproduced next to
    the name of the attribute or method.

    Attributes
    ----------
    spam
    eggs

    Methods
    -------
    run_one_step
    do_a_thing
    """

    _name = 'ComponentMaker'

    _cite_as = """Insert a BibTeX formatted reference here"""

    _input_var_names = (
        'topographic__elevation',
    )

    _output_var_names = (
        'component_maker__field',
    )

    _var_units = {
        'topographic__elevation': 'm',
        'component_maker__field': '-',
    }

    _var_mapping = {
        'topographic__elevation': 'node',
        'component_maker__field': 'node',
    }

    _var_doc = {
        'topographic__elevation':
            'Land surface topographic elevation',
        'component_maker__field':
            'A description of the fields made by the component',
    }

    def __init__(self, grid, spam, eggs=True, **kwds):
        """Initialize NetworkSedimentTransporter component.

        Parameters
        ----------
        grid : RasterModelGrid
            A grid.
        spam : bool
            A variable needed by the component.
        eggs : float, optional
            A variable needed by the component. Default value is 1.0.

        Examples
        --------
        We can put more examples here.

        """
        # In this line, we call the init method of Component, the class from
        # which this component is derived.
        super(ComponentMaker, self).__init__(grid)

        # Then we save an reference to the grid.
        self._grid = grid

        # the rest of initialization goes here.
        # if you have specific types or shapes, you should check that they have
        # been provided correctly. Then you should verify this in the unit test
        # file.
        if isinstance(spam, bool) is False:
            msg = ('Keyword argument to ComponentMaker spam is not of type '
                    'boolean. This is not permitted.')
            raise ValueError(msg)

        if isinstance(eggs, float) is False:
            msg = ('Keyword argument to ComponentMaker eggs is not of type '
                   'float. This is not permitted.')
            raise ValueError(msg)

        self._spam = spam
        self._eggs = eggs

        # You might need to make a new model grid field, its worth checking if
        # exists first though.
        if 'component_maker__field' in self._grid.at_node:
            self.cmf = self._grid.at_node['component_maker__field']
        else:
            self.cmf = self._grid.add_zeros('node', 'component_maker__field')

        # As you develop your component, the best thing to do is to look at
        # the source code of other components and to make an issue if you need
        # help or get stuck.

    @property
    def spam(self):
        """spam."""
        return self._spam

    @property
    def eggs(self):
        """eggs."""
        return self._eggs

    def do_a_thing(self):
        """Do a thing for the ComponentMaker.

        This function sets the value of the field `component_maker__field` to
        random numbers. It is meant as an example of how you can use multiple
        functions within a component.

        If you have Parameter, Examples, or values that this function Returns,
        you can make sections in the docstring for that.
        """
        self.cmf[:] = self._eggs * np.random.rand(self._grid.size('node'))

    def _private_function(self):
        """Sometimes you need a private function.

        A private function starts with an underscore. It is used by your class,
        but it isn't in the documention.

        This private function doesn't do anything.
        """
        pass

    def run_one_step(self, dt):
        """Run ComponentMaker forward in time by dt.

        More text here describing what run_one_step does.

        Parameters
        ----------
        dt : float
            Timestep

        Examples
        --------

        """
        if self._spam:
            self.do_a_thing()
        else:
            self._private_function()
        self.cmf += dt
