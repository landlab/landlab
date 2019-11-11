.. _one_to_two:

Transition from Landlab 1.x to 2.0
==================================

Landlab 2.0 will be released in early 2020. This page is meant to summarize the
major changes in Landlab 2.0 and make it easier for you to transition your
code.

Changes include:

- Drop python 2 support.
- The Landlab grid now inherits from the Landlab graph.
- Deprecation of some grid methods.
- Standardization and deprecation within the component library.
- Some other functions/methods have been removed.
- No more wiki (all the docs are here).

No more Python 2
----------------

TODO: put a link to general 2 to 3 resources


Landlab grid now inherits from the Landlab graph
------------------------------------------------

TODO List of user-facing changes

Deprecation of some grid methods
--------------------------------

TODO List of these, including replacements if any.

Changes to field creation
-------------------------
- no more `noclobber`, instead we have `clobber`
- recommended field init (eg. at="")

Changes to Boundary Condition Flags
-----------------------------------
- no more use of unbound flags AND unbound flags are not importable
  from top level or grid submodule (e.g., `from landlab import CLOSED_NODE`)
  now we use `Grid.BC_NODE_IS_CLOSED`.

  Old name                  New name
  ========================= ====================================
  BAD_INDEX_VALUE           ModelGrid.BAD_INDEX
  CORE_NODE                 ModelGrid.BC_NODE_IS_CORE
  FIXED_VALUE_BOUNDARY      ModelGrid.BC_NODE_IS_FIXED_VALUE
  FIXED_GRADIENT_BOUNDARY   ModelGrid.BC_NODE_IS_FIXED_GRADIENT
  LOOPED_BOUNDARY           ModelGrid.BC_NODE_IS_LOOPED
  CLOSED_BOUNDARY           ModelGrid.BC_NODE_IS_CLOSED
  ACTIVE_LINK               ModelGrid.BC_LINK_IS_ACTIVE
  INACTIVE_LINK             ModelGrid.BC_LINK_IS_INACTIVE
  FIXED_LINK                ModelGrid.BC_LINK_IS_FIXED


Standardization and deprecation within the component library
------------------------------------------------------------

- FlowRouter deprecated. Use FlowAccumulator
- No more flooded_nodes (use ``erode_flooded_nodes=True`` at init)
- Fewer different options for runoff rate, but one, consistent, good option
  (use argument to FlowAccumulator, then use ``surface__water_discharge`` for
  erosion).
- Removal of old deprecated methods such as (``erode``, ``diffuse``...)
- Main method now either take ``dt`` or nothing.

Some other functions/methods have been removed
----------------------------------------------

ModelParameterDictionary
plot.channel_profile (use ChannelProfiler)
