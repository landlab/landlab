.. _landlab_components_page:

=====================
The Component Library
=====================

Landlab offers an ever-growing library of components that aim to describe
individual or closely associated suites of surface processes. Components are
designed to be "plug-and-play" and to interact with each other with the minimum
of technical difficulties. Each component makes use of Landlab grid fields to
enable the sharing of data between the components, and we aim to have a
relatively standardized way of interacting with and using each different one.

Landlab components exist as classes, and can be imported from
``landlab.components``.

To develop your own Landlab component, see
:ref:`this page <dev_contributing>`
and
`this tutorial <https://mybinder.org/v2/gh/landlab/landlab/master?filepath=notebooks/tutorials/making_components/making_components.ipynb>`_.

Component Library Tutorial
----------------------------
For a tutorial introduction to using the component library, see
`here <https://mybinder.org/v2/gh/landlab/landlab/master?filepath=notebooks/tutorials/component_tutorial/component_tutorial.ipynb>`_.

Available Landlab components
----------------------------

For the complete list of Landlab components type the following command in a
command prompt:

``landlab list``

See the :ref:`Components section <api.components>` of the Landlab reference
manual for a list of all Landlab components currently available.

Landlab component classes, their import, and their instantiation
----------------------------------------------------------------

Almost all Landlab components exist as Python classes. This means that to use
them, you must first import the class, then ``instantiate`` a Python instance
of that class, then call a method from the class to run the component. The way
this is done has now been almost totally standardised across Landlab.

A component class is imported from the library as

.. code-block:: python

    from landlab.components import NameOfComponentToImport

e.g., to get the linear diffusion component, we would do:

.. code-block:: python

    from landlab.components import LinearDiffuser

The available components are listed in the
:ref:`Components section <api.components>` of the Landlab Reference Manual.

Component classes always take a copy of the grid as their first argument. They
then take a sequence of additional keyword arguments that set the actual
parameters for the component. This means that the instantiation of a component
looks something like
this:

.. code-block:: python

    dfn = LinearDiffuser(grid, linear_diffusivity=0.01)

These keywords can also be set by passing a Python dictionary, or using a text
input file (see below).

Here, `dfn` is now the `component object`—an "instance" of the component. We
can run it by calling its run method. The component's documentation will
explain how to do this for each individual case, but typically a component will
have a method called `run_one_step`, which can be called like this:

.. code-block:: python

    dt = 100.0  # the timestep
    dfn.run_one_step(dt)

If the component describes a time-varying process, the first argument of
`run_one_step` will be the duration for which to run the component in this
timestep. (If the component is not time sensitive, e.g., the `FlowRouter`,
it won't take `dt`). Some components may also allow/require additional input
parameters to their run method; see individual component documentation for more
details.

Running one of these methods will update the fields held in common by the
single grid object which you linked to all your components during component
instantiation. If you look inside the grid fields having run one of these
methods, you'll see the new fields it has created and populated. The docstrings
for the component should make it clear which fields the component needs to have
in the grid as inputs, and which it modifies and/or creates as outputs.

It should probably be emphasized here to **always read the documentation for
the component you are using**! You can get at this documentation either on this
website, or in a dynamic Python session by getting help for either the imported
class or the instantiated component object. i.e., in this case, any of the
following would work::

.. code-block:: pycon

    >>> help(LinearDiffuser)
    >>> help(dfn)  # LinearDiffuser? or dfn? also works

Quit interactive help in iPython by pressing "q".


.. _input_files:

Inputs to components
--------------------
Landlab components are initialized by passing a copy of the grid, then by
passing additional dynamic Python keyword arguments, almost all of which are
set to default values if a value is not provided. This means all of the ways
that you could call any other Python function using keywords also applies to
our components.

Most simply, components can be initialized by passing only the keyword values
that need to deviate from the defaults. So, for example, the default parameter
values for the `FastscapeEroder` are
`K_sp=None, m_sp=0.5, n_sp=1., threshold_sp=0., rainfall_intensity=1.`. So if
I want to set the `K_sp` to, say, `1.e-6`, but I am happy with these other
parameters, I can simply do:

.. code-block:: python

    fsc = FastscapeEroder(grid, K_sp=1.0e-6)

Because Landlab components make use of Python's native `**kwargs` argument
syntax, we can also pass multiple keywords at once to a component using a
Python dictionary:

.. code-block:: python

    sp_thresholds = grid.add_ones("node", "sp_thresholds")
    myargs = {"K_sp": 1.0e-5, "rainfall_intensity": 0.5, "threshold_sp": sp_thresholds}
    fsc = FastscapeEroder(grid, **myargs)

Note the "magic" `**` decorator that is placed on the dictionary when it is
passed to the component that makes this work. Also note that we can allow the
component default values to continue to set any keywords we still don't want to
supply, and that as long as the component permits it, we can pass in arrays or
field names like this too (see, e.g., `threshold_sp` above). You can have all
of your input parameters for all components in one dictionary if you so wish;
components will ignore any keywords they are passed that they don't recognize.

**Note that Landlab components will raise an error if they are passed
keyword arguments that they do not need.**

Landlab components always want to see a Python dictionary as their input, as
illustrated above. However, Landlab does offer a native file
reader called `load_params` that allows you to create dictionaries to pass to
components from input files. This function recognizes both
`"yaml" <https://yaml.org/spec/1.2/>`_ formatted data files, e.g.,

.. code-block:: yaml

    K_sp: 0.3
    m_sp: 0.5
    n_sp: 1.
    linear_diffusivity: 0.0001

The `load_params` method will figure out which to use by itself, and will do
any necessary typecasting automatically (i.e., floats will be floats, not
strings):

.. code-block:: python

    from landlab import load_params

    my_input_dict = load_params("./mytextinputfile.txt")
    dfn = FastscapeEroder(grid, **my_input_dict)

Component standard properties
-----------------------------

All Landlab components offer a standardized interface. This provides automated information
on the fields, units, etc. that the component works with, creates, and/or modifies. For a
fully compliant component, you will find you can call the following methods and attributes.


+------------------------------------------------------------------------------------------------------+--------------------------------------------------------+
| Property                                                                                             | Description                                            |
+======================================================================================================+========================================================+
| :py:meth:`Component.name <landlab.core.model_component.Component.name>`                              | a string                                               |
+------------------------------------------------------------------------------------------------------+--------------------------------------------------------+
| :py:meth:`Component.input_var_names <landlab.core.model_component.Component.input_var_names>`        | a tuple giving input field names                       |
+------------------------------------------------------------------------------------------------------+--------------------------------------------------------+
| :py:meth:`Component.output_var_names <landlab.core.model_component.Component.output_var_names>`      | a tuple giving output field names                      |
+------------------------------------------------------------------------------------------------------+--------------------------------------------------------+
| :py:meth:`Component.var_loc <landlab.core.model_component.Component.var_loc>`                        | a tuple of (var_name, ['node', 'link', etc])           |
+------------------------------------------------------------------------------------------------------+--------------------------------------------------------+
| :py:meth:`Component.definitions <landlab.core.model_component.Component.definitions>`                | a tuple of pairs of (var_name, short description)      |
+------------------------------------------------------------------------------------------------------+--------------------------------------------------------+
| :py:meth:`Component.units <landlab.core.model_component.Component.units>`                            | a tuple of (var_name, ['m', 'Pa', etc])                |
+------------------------------------------------------------------------------------------------------+--------------------------------------------------------+
| :py:meth:`Component.var_units('field') <landlab.core.model_component.Component.var_units>`           | method to return the unit of 'field'                   |
+------------------------------------------------------------------------------------------------------+--------------------------------------------------------+
| :py:meth:`Component.var_definition('field') <landlab.core.model_component.Component.var_definition>` | method to return a short description of 'field'        |
+------------------------------------------------------------------------------------------------------+--------------------------------------------------------+
| :py:meth:`Component.var_mapping('field') <landlab.core.model_component.Component.var_mapping>`       | method to return the element of 'field' (e.g., 'node') |
+------------------------------------------------------------------------------------------------------+--------------------------------------------------------+
| :py:meth:`Component.var_type('field') <landlab.core.model_component.Component.var_type>`             | method to return dtype of 'field' (e.g., float)        |
+------------------------------------------------------------------------------------------------------+--------------------------------------------------------+
| :py:meth:`Component.var_help('field') <landlab.core.model_component.Component.var_help>`             | a text summary of all of this information for 'field'  |
+------------------------------------------------------------------------------------------------------+--------------------------------------------------------+


See :ref:`the tutorials <tutorials>` for
examples of use cases with one, two, and more coupled components.

You can also get an overview of field usage by all components through Landlab's
command line interface. See
:ref:`here <getting_info_about_fields>`
for more information.

.. _component_standard_names:

Landlab standard naming conventions
-----------------------------------

The Landlab component library attempts to make use of a relatively standardized set of names across
the various components, in order to maximize ease of component coupling. If you're familiar with
the concept of the `CSDMS standard naming conventions
<https://csdms.colorado.edu/wiki/CSDMS_Standard_Names>`_, note that we have tried to strike a balance
between the rigor and uniqueness of those names and a more user-friendly, succinct approach.
Nonetheless, you may recognize the basic style of the names:

	**thing_described__what_is_described**

e.g., *topographic__elevation*, *water_surface__gradient*, *water__volume_flux*

We compile three tables to assist users with the :ref:`Landlab standard names<standard_name_definitions>`.

- First is a list of all names with their definitions.
- Second is a table listing which components use each field.
- Third is a table listing which components provide each field.


See :ref:`here <standard_name_changes>` for a list of changes to the standard
name list associated with the release of Landlab version 1.x (relative to 0.x).


Dealing with nonstandard names
++++++++++++++++++++++++++++++

The large number of developers on Landlab and historical accident have meant that despite our
best efforts you'll inevitably find instances where different components use different names
for the same thing. In these cases, you need to make equivalent two fields in the grid which
have different names so that two components can talk to each other. This is actually easy;
you can just do:

>>> mg.add_field("node", "second_name", mg.at_node["first_name"])

Note that we are making slow progress towards truly standardizing the component library, but
these kind of idiosyncrasies might yet persist for a while!
