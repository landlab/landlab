.. _landlab_components_page:

The Component Library
=====================

Landlab offers an ever-growing library of components that aim to describe individual or closely associated suites of surface processes. Components are designed to be “plug-and-play” and to interact with each other with the minimum of technical difficulties. Each component makes use of Landlab grid fields to enable the sharing of data between the components, and we aim to have a relatively standardized way of interacting with and using each different one.

Landlab components exist as classes, and can be imported from *landlab.components*.

    We emphasize that at the moment most Landlab components are under active development, and it is possible that    changes may occur to things such as input argument format with little or no warning. Please let the Landlab development team know if you’re making heavy use of a component so we can avoid unnecessarily breaking your code in future Landlab releases! Components are presented as-is in this beta release and we don't make any guarantees about the stability or otherwise of the implementations.

To develop your own Landlab component, see `this page <https://github.com/landlab/landlab/wiki/Develop-your-own-component>`_ and `this tutorial <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/master/making_components/making_components.ipynb>`_.

Component Library Tutorial
----------------------------
For a tutorial introduction to using the component library, see `here <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/master/component_tutorial/component_tutorial.ipynb>`_.

Available Landlab components
----------------------------

For the complete list of Landlab components type the following command in a command prompt:

``landlab list``

See the `Components section <http://landlab.readthedocs.io/en/release/#components>`_ of the Landlab reference manual for a list of all Landlab components currently available.


Landlab component classes, their import, and their instantiation
----------------------------------------------------------------

Almost all Landlab components exist as Python classes. This means that to use them, you
must first import the class, then ``instantiate`` a Python instance of that class,
then call a method from the class to run the component. The way this is done has now
been almost totally standardised across Landlab.

A component class is imported from the library as

>>> from landlab.components import [ComponentClass]

e.g., to get the linear diffusion component, we would do::

>>> from landlab.components import LinearDiffuser

The available components are listed
in the `Landlab Reference Manual <http://landlab.readthedocs.io/en/latest/#components>`_.

Component classes always take a copy of the grid as their first argument. They then take a
sequence of additional keyword arguments that set the actual parameters for the component.
This means that the instantiation of a component looks something like this::

>>> dfn = LinearDiffuser(grid, 'linear_diffusivity'=0.01)

These keywords can also be set by passing a Python dictionary, or using a text input file
(see below).

Here, `dfn` is now the `component object`—an "instance" of the component. We can run it
by calling its run method. The component's documentation will explain how to do this for
each individual case, but typically a component will have a method called `run_one_step`,
which can be called like this::

>>> dt = 100.  # the timestep
>>> dfn.run_one_step(dt)

If the component describes a time-varying process, the first argument of `run_one_step`
will be the duration for which to run the component in this timestep. (If the component
is not time sensitive, e.g., the `FlowRouter`, it won't take `dt`). Some components
may also allow/require additional input parameters to their run method; see individual
component documentation for more details.

Running one of these methods will update the fields held in common by the single grid
object which you linked to all your components during component instantiation. If you
look inside the grid fields having run one of these methods, you’ll see the new fields
it has created and populated. The docstrings for the component should make it clear
which fields the component needs to have in the grid as inputs, and which it modifies
and/or creates as outputs.

It should probably be emphasized here to **always read the documentation for the
component you are using**! You can get at this documentation either on this
website, or in a dynamic Python session by getting help for either the imported class
or the instantiated component object. i.e., in this case, any of the following would
work::

>>> help(LinearDiffuser)
>>> help(dfn)
>>> LinearDiffuser?
>>> dfn?

Quit interactive help in iPython by pressing "q".


.. _input_files:

Inputs to components
--------------------

Landlab components are initialized by passing a copy of the grid, then by passing additional
dynamic Python keyword arguments, almost all of which are set to default values if a value
is not provided. This means all of the ways that you could call any other Python function
using keywords also applies to our components.

Most simply, components can be initialized by passing only the keyword values that need to
deviate from the defaults. So, for example, the default parameter values for the
`FastscapeEroder` are `K_sp=None, m_sp=0.5, n_sp=1., threshold_sp=0., rainfall_intensity=1.`.
So if I want to set the `K_sp` to, say, `1.e-6`, but I am happy with these other parameters,
I can simply do::

>>> fsc = FastscapeEroder(grid, K_sp=1.e-6)

Because Landlab components make use of Python's native `**kwargs` argument syntax, we can
also pass multiple keywords at once to a component using a Python dictionary::

>>> sp_thresholds = grid.add_ones('node', 'sp_thresholds')
>>> myargs = {'K_sp': 1.e-5, 'rainfall_intensity': 0.5, 'threshold_sp': sp_thresholds}
>>> fsc = FastscapeEroder(grid, **myargs)

Note the "magic" `**` decorator that is placed on the dictionary when it is passed to the
component that makes this work. Also note that we can allow the component default values to
continue to set any keywords we still don't want to supply, and that as long as the
component permits it, we can pass in arrays or field names like this too (see, e.g.,
`threshold_sp` above). You can have all of your input parameters for all components in
one dictionary if you so wish; components will ignore any keywords they are passed that
they don't recognize.

Landlab components always want to see a Python dictionary as their input, as illustrated
above (passing a string for an actual input file is deprecated functionality in Landlab
version 1). However, Landlab does offer a native file reader called `load_params` that
allows you to create dictionaries to pass to components from input files. This function
recognizes both `"yaml" <http://www.yaml.org/start.html>`_ formatted data files, e.g.,::

    K_sp: 0.3
    m_sp: 0.5
    n_sp: 1.
    linear_diffusivity: 0.0001

...or our own Landlab-native ModelParameterDictionary format::

    K_sp: Additional text can go after the colon to allow annotation
    0.3
    m_sp:
    0.5
    n_sp:
    1.
    linear_diffusivity:
    0.0001

The `load_params` method will figure out which to use by itself, and will do any
necessary typecasting automatically (i.e., floats will be floats, not strings)::

>>> from landlab import load_params
>>> my_input_dict = load_params('./mytextinputfile.txt')
>>> dfn = FastscapeEroder(grid, **my_input_dict)

It is possible to mix passing of parameters by input dictionary and manually specifying
them, like this: `dfn = LinearDiffuser(grid, K_sp=1.e-6, **mydict)`, but we don't
recommend it as it can easily lead to confusion and errors.


Component standard properties
-----------------------------

All Landlab components offer a standardized interface. This provides automated information
on the fields, units, etc. that the component works with, creates, and/or modifies. For a
fully compliant component, you will find you can call:

==================================  ===========================================================
Property                            Description
==================================  ===========================================================
component.name 		                a string
component.input_var_names 	        a tuple giving input field names
component.output_var_names	        a tuple giving output field names
component.var_loc		        a tuple of (var_name, [‘node’, ‘link’, etc])
component.definitions	                a tuple of pairs of (var_name, short description)
component.units                         a tuple of (var_name, [‘m’, ‘Pa’, etc])
component.var_units('field')            method to return the unit of 'field'
component.var_definition('field')       method to return a short description of 'field'
component.var_mapping('field')          method to return the element of 'field' (e.g., 'node')
component.var_type('field')             method to return dtype of 'field' (e.g., float)
component.var_help('field')             a text summary of all of this information for 'field'
==================================  ===========================================================

See `the tutorials <https://github.com/landlab/landlab/wiki/Tutorials>`_ for examples of use cases
with one, two, and more coupled components.

You can also get an overview of field usage by all components through Landlab's command line
interface. See `here <https://github.com/landlab/landlab/wiki/Grid#getting-information-about-fields>`_
for more information.


.. _standard_names:

Landlab standard naming conventions
-----------------------------------

The Landlab component library attempts to make use of a relatively standardized set of names across
the various components, in order to maximize ease of component coupling. If you’re familiar with
the concept of the `CSDMS standard naming conventions
<http://csdms.colorado.edu/wiki/CSDMS_Standard_Names>`_, note that we have tried to strike a balance
between the rigor and uniqueness of those names and a more user-friendly, succinct approach.
Nonetheless, you may recognize the basic style of the names:

	**thing_described__what_is_described**

e.g., *topographic__elevation*, *water_surface__gradient*, *water__volume_flux*

 You can see a list of the names currently in use here: `Landlab Standard Names <https://github.com/landlab/landlab/wiki/Standard-names>`_

See `here <https://github.com/landlab/landlab/wiki/Standard-names#changes-to-standard-names-in-landlab>`_ for a list of recent changes
to the standard name list.


Dealing with nonstandard names
++++++++++++++++++++++++++++++

The large number of developers on Landlab and historical accident have meant that despite our
best efforts you’ll inevitably find instances where different components use different names
for the same thing. In these cases, you need to make equivalent two fields in the grid which
have different names so that two components can talk to each other. This is actually easy;
you can just do:

>>> mg.add_field(‘node’, ‘second_name’, mg.at_node[‘first_name’])

Note that we are making slow progress towards truly standardizing the component library, but
these kind of idiosyncrasies might yet persist for a while!
