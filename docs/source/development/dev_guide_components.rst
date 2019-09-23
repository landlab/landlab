.. _dev_components:

=================================
Landlab Component Developer Guide
=================================

This section under construction!

For now, see the main section on `Landlab components <https://github.com/landlab/landlab/wiki/Components>`_, or dive right into the `component developer guide <https://github.com/landlab/landlab/wiki/Develop-your-own-component>`_ if you feel like you already know all the basics.
Feel free to contact the development team for more specific advice.

Recommendations and Rules for User-Contributed Components
---------------------------------------------------------

First, thank you for considering contributing a component! The following are
some rules and recommendations.


Erkan suggests: we write that the reason we do this


Rules
`````

- ``super`` is called in the init

- All public attributes are properties
- A component uses all required inputs fields and creates all required output fields.

- Components have complete metadata in ``_info``.

- The component's ``__init__`` method takes a Landlab model grid as the first argument.

- The component has a main method that takes either ``dt`` or nothing. The name is descriptive but does not need to be standardized.
    * Solution is that using setters and getters or making a grid scalar.
    * Some common name patterns include "update", "run_one_step", or "calculate_..."

- A component raises a ``ValueError`` if unused keyword arguments are provided.
- A component raises a ``ValueError`` if a grid type the component does not support is passed.

Recommendations
```````````````

- A component has a ``run_one_step`` method as its main method.
- Arguments and keyword arguments start with lower case letters.
- Keyword arguments have reasonable default values (and the grid is the only argument to ``__init__``.
- Can return or not return information. Some recommendations, Nothing, grid, or a ... or a calculated value. 
