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

No more Python 2 #attnKaty
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

Changes to Boundary Condition Flags
-----------------------------------
- no more use of unbound flags (e.g., `from landlab import CLOSED_NODE`)
  now we use `Grid.BC_NODE_IS_CLOSED`.

Changes to Field initialization
-------------------------------
- PUT HERE IF THIS HAPPENS

Standardization and deprecation within the component library
------------------------------------------------------------

- FlowRouter deprecated. Use FlowAccumulator
- No more flooded_nodes (use ``erode_flooded_nodes=True`` at init)
- Fewer different options for runoff rate, but one, consistent, good option
  (use argument to FlowAccumulator, then use ``surface__water_discharge`` for
  erosion).
- Removal of old deprecated methods such as (``erode``, ``diffuse``...) #attnKaty
- Main method now either take ``dt`` or nothing.

Some other functions/methods have been removed
----------------------------------------------

ModelParameterDictionary
plot.channel_profile (use ChannelProfiler) #attnKaty
