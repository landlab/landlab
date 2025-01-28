# FlowDirAccPf: efficient filling and flow routing



``FlowDirAccPf`` is a ``Landlab`` component that provides an alternative and efficent approach to fill or breach DEMs, calculate flow directions and update flow accumulations. The component is restricted to structured grids and contains a wrapper for the RichDEM python package [@barnes2016parallel,@barnes2017parallel]. [``RichDEM``](https://richdem.readthedocs.io/en/latest/intro.html) is a set of hydrologic analysis tools using parallel processing to process large DEMs and calculate hydrologic properties.

#TODO

FlowDirAccPf is introduced [in this paper]()


## Documentation and installation

Landlab documentation is hosted on this [ReadTheDocs page](https://landlab.csdms.io/),
including instructions to install Landlab. ``FlowDirAccPf`` is installed with
Landlab.

#TODO

``FlowDirAccPf`` documentation is located [here](https://landlab.csdms.io/generated/api/landlab.components.priority_flood_flow_router.priority_flood_flow_router.html).

## FlowDirAccPf tutorial

A ``FlowDirAccPf`` tutorial exists in the form of a Jupyter Notebooks accessible
through the following links:

#TODO adjust links when in main branch. [This](https://github.com/BCampforts/landlab/blob/bc/priority_flood/notebooks/tutorials/flow_direction_and_accumulation/PriorityFlood_realDEMs.ipynb) is a direct link to the notebook.

- [Launch the tutorial](https://mybinder.org/v2/gh/BCampforts/landlab/blob/bc/priority_flood/notebooks/tutorials/PriorityFlood/PriorityFlood_realDEMs.ipynb)
as interactive notebook in your browser, with no need to install software,
launched using Binder.
- [A static version of the same tutorial](https://landlab.csdms.io/tutorials/flow_direction_and_accumulation/PriorityFlood_realDEMs.html)
- All Landlab tutorials can be launched from [this directory](https://mybinder.org/v2/gh/landlab/landlab/release?filepath=notebooks/welcome.ipynb) using Binder.

## Get or give help

[Open an Issue here](https://github.com/landlab/landlab/issues) where we can
respond to your questions, comments, issues, ideas, or any identified bugs
related to Landlab including ``FlowDirAccPf``.
