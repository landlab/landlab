[Landlab](http://landlab.github.io) |
[[About |About]] |
[[Examples |Examples]] |
[[User Guide | User-Guide]] |
[Reference Manual](http://landlab.readthedocs.org/en/latest/#developer-documentation) |
[[Tutorials| Tutorials ]] |
[[FAQs |FAQs]]

[[ Next topic: Intro to Python → | Python,-NumPy,-SciPy,-Cython ]]

## Installation
 - [[ Instructions for a standard install | Installing-Landlab ]]
 - [[ Installing from source code, "developer install" | Installing-Landlab-from-source-code-("developer-install") ]]
 - [[ Help with installing GRASS GIS once Anaconda is on your machine | Installing-GRASS-after-installing-Landlab ]]

## Basics of Python
If you are new to Python or scientific programming, start with an intro to the nuts and bolts of Landlab:

[[Python, NumPy, SciPy, and Cython | Python,-NumPy,-SciPy,-Cython ]]
  - [[ Why Python? | Python,-NumPy,-SciPy,-Cython#cython#why-python ]]
  - [[ Getting to know Python | Python,-NumPy,-SciPy,-Cython#getting-to-know-python ]]
    - [[ If you know MatLab… | Python,-NumPy,-SciPy,-Cython#if-you-know-matlab ]]
  - [[ NumPy, SciPy, and efficient coding style | Python,-NumPy,-SciPy,-Cython#numpy-scipy-and-efficient-coding-style ]]
  - [[ Cython | Python,-NumPy,-SciPy,-Cython#cython ]]

## Landlab's grid
An explanation of Landlab's Grid can be found here:

[[ Landlab's gridding library | grid ]]
  - [[ How a grid is represented | Grid#how-a-grid-is-represented ]]
      - [[ Basic grid elements | Grid#basic-grid-elements ]]
      - [[ Creating a grid | Grid#creating-a-grid ]]
      - [[ Using fields | Grid#adding-data-to-a-landlab-grid-element-using-fields ]]
        - [[ Field initialization | Grid#field-initialization ]]
        - [[ Field creation from existing data | Grid#field-creation-from-existing-data ]]
        - [[ Field access | Grid#field-access ]]
        - [[ Getting information about field usage | Grid#getting-information-about-fields ]]
      - [[ Representing gradients | Grid#representing-gradients-in-a-landlab-grid ]]
      - [[ Other grid elements | Grid#other-grid-elements ]]
    - [[ Managing grid boundaries | Grid#managing-grid-boundaries ]]
      - [[ Boundary condition details and methods | Grid#boundary-condition-details-and-methods ]]
    - [[ Using a different grid type | Grid#using-a-different-grid-type ]]
    - [[ Importing a DEM | Grid#importing-a-dem ]]
    - [[ Ploting and visualization | Grid#plotting-and-visualization ]]
      - [[ Visualizing a grid | Grid#visualizing-a-grid ]]
      - [[ Visualizing transects through your data | Grid#visualizing-transects-through-your-data ]]
      - [[ Visualizing river profiles | Grid#visualizing-river-profiles ]]
      - [[ Making movies | Grid#making-movies ]]

## How to model with Landlab
Build a model and learn how to create and use components and drivers:

[[ Building a Model | Build a Model ]]
- [[ What goes into a Landlab model? | Build-a-Model#what-goes-into-a-landlab-model ]]
- [[ A brief introduction to components | Build-a-Model#a-brief-introduction-to-components ]]
- [[ Implementing a Landlab driver | Build-a-Model#implementing-a-landlab-driver ]]
  - [[ Import the libraries and functions you need | Build-a-Model#1-import-the-libraries-and-functions-you-need ]]
  - [[ Instantiate objects | Build-a-Model#2-instantiate-objects ]]
  - [[ Load/create data fields | Build-a-Model#3-loadcreate-data-in-fields ]]
  - [[ Set boundary conditions | Build-a-Model#4-set-the-boundary-conditions ]]
  - [[ Finalize and handle the data | Build-a-Model#6-finalize-and-handle-the-data ]]
    - [[ Save or export the data | Build-a-Model#save-or-export-the-data ]]
    - [[ Plot the data | Build-a-Model#plot-the-data ]]

[[Components |Components]]
- [[ The Landlab component library | Components#the-landlab-component-library ]]
- [[ Available Landlab components | Components#available-landlab-components ]]
- [[ Input files | Components#input-files ]]
- [[ Implementing a component | Components#implementing-a-component ]]
- [[ Component standard properties | Components#component-standard-properties ]]
- [[ Landlab standard naming conventions | Components#landlab-standard-naming-conventions ]]
  - [[ Dealing with nonstandard names | Components#dealing-with-nonstandard-names ]]

## Cellular automaton functionality
Build a model using Landlab's CellLab-CTS module:

[[ CellLab-CTS User Guide | CellLab-CTS-Users-Manual ]]
- [[ Introduction | CellLab-CTS-Users-Manual#introduction ]]
- [[ Writing a Celllab CTS model | CellLab-CTS-Users-Manual#writing-a-celllab-cts-model ]]
  - [[ What is a Celllab CTS model? | CellLab-CTS-Users-Manual#what-is-a-celllab-cts-model ]]
  - [[ Basic ingredients of Celllab CTS model | CellLab-CTS-Users-Manual#basic-ingredients-of-a-celllab-cts-model ]]
  - [[ Types of Celllab CTS model | CellLab-CTS-Users-Manual#types-of-celllab-cts-model ]]
  - [[ Step 1: Importing Celllab CTS | CellLab-CTS-Users-Manual#step-1-importing-celllab-cts ]]
    - [[ Setting up transitions | CellLab-CTS-Users-Manual#setting-up-transitions ]]
    - [[ Defining parameters | CellLab-CTS-Users-Manual#defining-parameters ]]
  - [[ Step 2: Creating a grid | CellLab-CTS-Users-Manual#step-2-creating-a-grid ]]
  - [[ Step 3: Create a node-state dictionary | CellLab-CTS-Users-Manual#step-3-create-a-node-state-dictionary ]]
  - [[ Step 4: Create the transition list | CellLab-CTS-Users-Manual#step-4-create-the-transition-list ]]
  - [[ Step 5: Create an array with initial node-state values | CellLab-CTS-Users-Manual#step-5-create-an-array-containing-the-initial-node-state-values ]]
  - [[ Step 6: Instantiate a Celllab CTS object | CellLab-CTS-Users-Manual#step-6-instantiate-a-celllab-cts-object ]]
  - [[ Step 7: Set up plotting | CellLab-CTS-Users-Manual#step-7-set-up-plotting ]]
  - [[ Step 8: Run the model | CellLab-CTS-Users-Manual#step-8-run-the-model ]]
  - [[ Step 9: Clean up | CellLab-CTS-Users-Manual#step-9-cleanup ]]
- [[ Reference information | CellLab-CTS-Users-Manual#reference-information ]]
- [[ References | CellLab-CTS-Users-Manual#references ]]

## [[Tutorials |Tutorials]]

Learn more with [[Tutorials |Tutorials]]

## FAQs
Read other frequently asked questions or ask your own:

[[FAQs |FAQs]]
- [[ What is the difference between a cell and a node? | FAQs#What-is-the-difference-between-a-cell-and-a-node ]]
- [[ Why is my node data a 1d array? | FAQs#why-is-my-node-data-a-1d-array-im-using-a-raster ]]
- [[ How do I set the boundary codes for the edges of a grid? | FAQs#how-do-i-set-the-boundary-codes-for-the-edges-of-a-grid ]]
- [[ Can I import Landlab output into ParaView or VisIt? | FAQs#can-i-import-landlab-output-into-paraview-or-visit ]]
- [[ How do I get netCDF output? | FAQs#how-do-i-get-netcdf-output ]]
- [[ How do I assign values from nodes to links? | FAQs#how-do-i-assign-values-from-nodes-to-links ]]
- [[ How do I test whether my grid is regular or irregular? | FAQs#how-do-i-test-whether-my-grid-is-regular-or-irregular ]]
- [[ How do I modify boundary conditions for part of the grid where I know the coordinates? | FAQs#how-do-i-modify-boundary-conditions-for-part-of-the-grid-where-i-know-the-coordinates]]
- [[ How do I keep in touch with Landlab developments? | FAQs#how-do-i-keep-in-touch-with-landlab-developments ]]


## Landlab User Group
To join the Landlab User Group on Slack, send an email request to:
<a href="MAILTO:knuth@colorado.edu?subject=Landlab User Group&body=Please send me an invite to join the Landlab User Group">Join Landlab User Group</a>
