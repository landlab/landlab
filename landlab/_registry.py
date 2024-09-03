#! /usr/bin/env python
r"""Registry of landlab components being used.

The landlab registry keeps track of landlab components that have
be instantiated by a user. A user can then get a list of all
the components in use and then print a list of citations for
all of the components they have used.

Examples
--------
>>> from landlab import registry

>>> registry.registered
('landlab',)
>>> print(registry.format_citations())
# Citations
## landlab
    @article{hobley2017creative,
    title={Creative computing with Landlab: an open-source toolkit
        for building, coupling, and exploring two-dimensional
        numerical models of Earth-surface dynamics},
    author={Hobley, Daniel EJ and Adams, Jordan M and Nudurupati,
        Sai Siddhartha and Hutton, Eric WH and Gasparini, Nicole M and
        Istanbulluoglu, Erkan and Tucker, Gregory E},
    journal={Earth Surface Dynamics},
    volume={5},
    number={1},
    pages={21},
    year={2017},
    publisher={Copernicus GmbH}
    }
    @article{barnhart2020short,
    author = {Barnhart, K. R. and Hutton, E. W. H. and Tucker, G. E.
        and Gasparini, N. M. and Istanbulluoglu, E. and Hobley,
        D. E. J. and Lyons, N. J. and Mouchene, M. and Nudurupati,
        S. S. and Adams, J. M. and Bandaragoda, C.},
    title = {Short communication: Landlab v2.0: A software package
        for Earth surface dynamics},
    journal = {Earth Surface Dynamics Discussions},
    volume = {2020},
    year = {2020},
    pages = {1--25},
    url = {https://www.earth-surf-dynam-discuss.net/esurf-2020-12/},
    doi = {10.5194/esurf-2020-12}
    }

When a component contains citation information, and the component has been
instantiated (not just imported) the component citation is also included.

>>> from landlab import RasterModelGrid
>>> from landlab.components import Flexure

>>> grid = RasterModelGrid((4, 5))
>>> _ = grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
>>> _ = grid.add_zeros("lithosphere_surface__elevation_increment", at="node")
>>> flexure = Flexure(grid)
>>> print(registry.format_citations())
# Citations
## landlab
    @article{hobley2017creative,
    title={Creative computing with Landlab: an open-source toolkit
        for building, coupling, and exploring two-dimensional
        numerical models of Earth-surface dynamics},
    author={Hobley, Daniel EJ and Adams, Jordan M and Nudurupati,
        Sai Siddhartha and Hutton, Eric WH and Gasparini, Nicole M and
        Istanbulluoglu, Erkan and Tucker, Gregory E},
    journal={Earth Surface Dynamics},
    volume={5},
    number={1},
    pages={21},
    year={2017},
    publisher={Copernicus GmbH}
    }
    @article{barnhart2020short,
    author = {Barnhart, K. R. and Hutton, E. W. H. and Tucker, G. E.
        and Gasparini, N. M. and Istanbulluoglu, E. and Hobley,
        D. E. J. and Lyons, N. J. and Mouchene, M. and Nudurupati,
        S. S. and Adams, J. M. and Bandaragoda, C.},
    title = {Short communication: Landlab v2.0: A software package
        for Earth surface dynamics},
    journal = {Earth Surface Dynamics Discussions},
    volume = {2020},
    year = {2020},
    pages = {1--25},
    url = {https://www.earth-surf-dynam-discuss.net/esurf-2020-12/},
    doi = {10.5194/esurf-2020-12}
    }
<BLANKLINE>
## Flexure
    @article{hutton2008sedflux,
    title={Sedflux 2.0: An advanced process-response model that
        generates three-dimensional stratigraphy},
    author={Hutton, Eric WH and Syvitski, James PM},
    journal={Computers \& Geosciences},
    volume={34},
    number={10},
    pages={1319--1337},
    year={2008},
    publisher={Pergamon}
    }

Finally, the component's citation information is accessible through an
attribute called ``cite_as``:

>>> print(Flexure.cite_as)
    @article{hutton2008sedflux,
    title={Sedflux 2.0: An advanced process-response model that
        generates three-dimensional stratigraphy},
    author={Hutton, Eric WH and Syvitski, James PM},
    journal={Computers \& Geosciences},
    volume={34},
    number={10},
    pages={1319--1337},
    year={2008},
    publisher={Pergamon}
    }

"""

import os

from . import _info
from .core.messages import indent_and_wrap


class ComponentRegistry:
    """A registry for instantiated landlab components."""

    def __init__(self, objs=None):
        self._registered = []
        if objs is not None:
            try:
                [self.add(obj) for obj in objs]
            except TypeError:
                self.add(objs)

    def add(self, cls):
        """Add a class to the registry.

        Parameters
        ----------
        cls : Component
            A landlab component to register as used.
        """
        if cls not in self._registered:
            self._registered.append(cls)

    @property
    def registered(self):
        """All registered classes.

        Returns
        -------
        tuple
            The components in the registry.

        Examples
        --------
        >>> from landlab._registry import ComponentRegistry
        >>> registry = ComponentRegistry()
        >>> registry.registered
        ()
        >>> class FooBar(object):
        ...     pass
        ...
        >>> registry.add(FooBar)
        >>> registry.registered
        ('FooBar',)
        """
        return tuple(ComponentRegistry.get_name(obj) for obj in self._registered)

    @staticmethod
    def format_citation(obj):
        """Format a single citation.

        Parameters
        ----------
        obj : Component
            A landlab component class or instance.

        Returns
        -------
        str
            The formatted citation, or "None" if there is no citation
            given.

        Examples
        --------
        >>> from landlab._registry import ComponentRegistry
        >>> registry = ComponentRegistry()
        >>> class DoNothingComponent(object):
        ...     pass
        ...
        >>> print(registry.format_citation(DoNothingComponent))
        ## DoNothingComponent
            None

        >>> class SorterAndSearcher(object):
        ...     _cite_as = '''
        ... @book{knuth1998art,
        ... title={The art of computer programming: sorting and searching},
        ... author={Knuth, Donald Ervin},
        ... volume={3},
        ... year={1998},
        ... publisher={Pearson Education}
        ... }'''
        ...
        >>> print(registry.format_citation(SorterAndSearcher))
        ## SorterAndSearcher
            @book{knuth1998art,
            title={The art of computer programming: sorting and searching},
            author={Knuth, Donald Ervin},
            volume={3},
            year={1998},
            publisher={Pearson Education}
            }
        """
        name = ComponentRegistry.get_name(obj)
        header = [f"## {name}"]

        cite_as = ComponentRegistry.get_citations(obj)

        body = []
        for citation in cite_as:
            body.append(indent_and_wrap(citation, indent=" " * 4))

        return os.linesep.join(header + body)

    @staticmethod
    def get_name(obj):
        """Get the display name for an object.

        Examples
        --------
        >>> from landlab._registry import ComponentRegistry
        >>> class MontyPython(object):
        ...     name = "Eric Idle"
        ...
        >>> ComponentRegistry.get_name(MontyPython)
        'Eric Idle'
        >>> class MontyPython(object):
        ...     _name = "Graham Chapman"
        ...
        >>> ComponentRegistry.get_name(MontyPython)
        'Graham Chapman'
        >>> class MontyPython(object):
        ...     pass
        ...
        >>> ComponentRegistry.get_name(MontyPython)
        'MontyPython'

        """
        name = "Unknown"
        for attr in ("name", "_name", "__name__"):
            try:
                name = getattr(obj, attr)
            except AttributeError:
                pass
            else:
                break
        return name

    @staticmethod
    def get_citations(obj):
        """Get a list of citations from an object."""
        citations = "None"
        for attr in ("cite_as", "_cite_as"):
            try:
                citations = getattr(obj, attr)
            except AttributeError:
                pass
            else:
                break
        if isinstance(citations, str):
            citations = [citations]
        return citations

    def format_citations(self):
        """Format citations for all registered components.

        Returns
        -------
        str
            The formatted citations.

        Examples
        --------
        >>> from landlab._registry import ComponentRegistry
        >>> registry = ComponentRegistry()

        >>> class HolyGrailFinder(object):
        ...     _name = "Monty Python"
        ...     _cite_as = [
        ...         '''@book{python2000holy,
        ...         title={the Holy Grail},
        ...         author={Python, Monty and Chapman, Graham and Cleese, John and Gilliam, Terry and Jones, Terry and Idle, Eric and Palin, Michael},
        ...         year={2000},
        ...         publisher={EMI Records}
        ...         }''',
        ...         '''@book{chapman1989complete,
        ...         title={The Complete Monty Python's Flying Circus: All the Words. Volume one},
        ...         author={Chapman, Graham and Python, Monty},
        ...         volume={1},
        ...         year={1989},
        ...         publisher={Pantheon}
        ...         }''',
        ...     ]
        ...
        >>> class Evolution(object):
        ...     _cite_as = '''
        ...         @book{darwin1859origin,
        ...         title={On the origin of species},
        ...         author={Darwin, Charles},
        ...         year={1859},
        ...         publisher={Lulu. com}
        ...         }'''
        ...
        >>> registry.add(HolyGrailFinder)
        >>> registry.add(Evolution)
        >>> print(registry.format_citations())
        # Citations
        ## Monty Python
            @book{python2000holy,
            title={the Holy Grail},
            author={Python, Monty and Chapman, Graham and Cleese, John and
                Gilliam, Terry and Jones, Terry and Idle, Eric and Palin,
                Michael},
            year={2000},
            publisher={EMI Records}
            }
            @book{chapman1989complete,
            title={The Complete Monty Python's Flying Circus: All the Words.
                Volume one},
            author={Chapman, Graham and Python, Monty},
            volume={1},
            year={1989},
            publisher={Pantheon}
            }
        <BLANKLINE>
        ## Evolution
            @book{darwin1859origin,
            title={On the origin of species},
            author={Darwin, Charles},
            year={1859},
            publisher={Lulu. com}
            }
        """  # noqa: B950
        header = ["# Citations"]
        body = []
        for cls in self._registered:
            body.append(self.format_citation(cls))
        return os.linesep.join(header + [(2 * os.linesep).join(body)])

    def __repr__(self):
        return f"ComponentRegistry({repr(self.registered)})"


registry = ComponentRegistry(_info)
