---
title: 'RiverFlowDynamics v1.0: A Landlab component for computing two-dimensional river flow dynamics'

tags:
  - Landlab
  - Python
  - Shallow water equations
  - Saint-Venant equation
  - river

authors:
  - name: Sebastian Bernal
    orcid: 0009-0006-7758-3648
    equal-contrib: true
    affiliation: 1
  - name: Angel Monsalve
    orcid: 0000-0002-7369-1602
    equal-contrib: true
    affiliation: 1
    corresponding: true
  - name: Oscar Link
    orcid: 0000-0002-2188-6504
    equal-contrib: true
    affiliation: 2

affiliations:
 - name: Center for Ecohydraulics Research, Civil and Environmental Engineering, University of Idaho, Boise, ID, USA
   index: 1
 - name: Departamento de Ingeniería Civil, Universidad de Concepción, Concepción, Chile
   index: 2

date: 05 September 2024

bibliography: paper.bib

---
# Summary

Numerical modeling of surface water flow is a critical tool in hydrology, hydraulics, and environmental science. These models play a crucial role in predicting and analyzing flow patterns in rivers, flood plains, and coastal areas, informing decisions in water resource management, flood risk assessment, and ecosystem conservation. This paper introduces a novel two-dimensional flow model, RiverFlowDynamics, developed as a component of the LandLab Python Package [@hobley:2017;@barnhart:2020;@Hutton:2020;@Hutton:2020], designed to simulate the behavior of rivers and streams under various flow conditions over natural and artificial topography

RiverFlowDynamics is founded on the depth-averaged Saint-Venant equations, also known as the shallow water equations [@casulli1990semi;@casulli_semi-implicit_1999]. These equations, derived from the Navier-Stokes equations for incompressible flow, are simplified by integrating over the water depth. This approach assumes that vertical accelerations are negligible compared to horizontal ones, a reasonable approximation for many surface water flows. The governing equations consist of continuity and momentum balance equations in two dimensions, capturing the essential dynamics of free-surface flows.

For the numerical solution of these equations, RiverFlowDynamics employs the finite volume method, chosen for its robustness, capacity to handle complex geometries, and inherent conservation properties [@andersson2011computational;@fletcher2012computational]. The computational domain is discretized into a uniform rectangular grid, with water surface elevation defined at cell centers and velocity components at cell interfaces. This staggered grid arrangement helps in maintaining numerical stability and accuracy. The numerical implementation results in a penta-diagonal, positive-definite system of equations. This system is solved efficiently using the preconditioned conjugate gradient method, ensuring rapid convergence even for large domains. The model's structure allows for easy parallelization, potentially enabling simulations of extensive river networks or large coastal areas.

A key feature of the model is its semi-implicit and semi-Lagrangian representation. The semi-implicit scheme treats water surface elevation and velocity implicitly, while handling advective terms explicitly. This approach allows for larger time steps compared to fully explicit schemes, enhancing computational efficiency. The semi-Lagrangian method for advection involves tracking fluid particles backwards along their flow lines, providing additional stability and accuracy, particularly for flows with strong advective components [@robert1985semi;@robert_stable_1981].

The advection representation uses a semi-Lagrangian scheme on an Eulerian grid, a two-step process where particles are tracked backwards along their flow lines to their starting points [@staniforth1991semi]. This method, combined with the semi-implicit time discretization, allows for a relaxation of the Courant-Friedrichs-Lewy (CFL) condition, typically a limiting factor in explicit schemes [@bates1982multiply].

Source terms in the model primarily account for bottom friction, implemented using the Manning-Chezy formula [@he2017numerical;@brufau2000two]. While the model framework allows for the inclusion of wind stress and Coriolis effects, these are considered negligible in the current implementation, focusing on river and stream applications where these effects are typically less significant.

Flow line tracing, crucial for the semi-Lagrangian advection scheme, employs Pollock's semi-analytical method [@pollock1988semianalytical]. This method assumes linear velocity variation within each grid cell, allowing for an efficient and accurate computation of particle trajectories. This approach is particularly effective in capturing complex flow patterns in natural river systems.

A notable feature of the model is its robust handling of dry/wet cell transitions, crucial for simulating flows over complex topography or in areas with varying water levels. RiverFLowDynamics employs the method of @casulli1992semi, in which the model automatically determines wet and dry cell faces based on local flow conditions, eliminating the need for explicit specification of internal boundaries.

Boundary conditions are implemented to handle both open and closed boundaries. Dirichlet conditions are used for inlet boundaries, specifying flow rates and water depths. For open boundaries where water can freely enter or exit the domain, the model offers both gradient-based and radiation-based conditions. These allow for the realistic simulation of wave propagation and minimize artificial reflections at the domain boundaries.

RiverFlowDynamics, as a full 2D flow model, offers several advantages over simpler flow models, including traditional overland flow models available in Landlab [@adams2017landlab;@de2012improving]. While overland flow models typically focus on shallow sheet flow and often use simplified equations like the kinematic wave approximation, RiverFlowDynamics solves the complete depth-averaged Saint-Venant equations. This approach allows for a more comprehensive representation of complex flow dynamics, including subcritical and supercritical flows, hydraulic jumps, and intricate channel-floodplain interactions. The model's ability to capture these phenomena makes it superior in scenarios involving rapid flood propagation in urban areas, detailed floodplain mapping, or the analysis of complex river morphodynamics. Furthermore, the semi-Lagrangian scheme employed in RiverFlowDynamics provides enhanced stability and accuracy for advection-dominated flows, a critical advantage when modeling high-velocity currents or steep terrain where simpler models might fail. This makes RiverFlowDynamics particularly well-suited for applications in mountainous regions, urban flood modeling, from small-scale stream dynamics to large-scale flood simulations, or any situation where capturing the full range of flow regimes and their transitions is crucial for accurate predictions. The accessibility of this code within the Landlab framework will make it easier for future users to modify and contribute to its continual evolution.

Source code for RiverFlowDynamics is available as part of the Landlab Python package (v2.7.0 and later) and can be found in the RiverFlowDynamics component. The Landlab project maintains a separate repository containing tutorials, including a complete example of RiverFlowDynamics usage in a Jupyter Notebook.

# Statement of need

RiverFlowDynamics is a Python-based 2D flow model developed as a component of the Landlab framework, addressing a critical gap in the modeling of complex river systems and flood dynamics. Prior to RiverFlowDynamics, Landlab lacked a comprehensive 2D flow model capable of handling fully advective-dominated problems, particularly in rivers with complex topographies. This limitation hindered accurate simulations of diverse flow regimes and transitions crucial for advanced hydrological and environmental studies.

RiverFlowDynamics solves the complete depth-averaged Saint-Venant equations, offering a significant advancement over existing Landlab components that typically use simplified equations like the kinematic wave approximation. This approach enables the model to capture complex flow dynamics, including subcritical and supercritical flows, hydraulic jumps, and intricate channel-floodplain interactions. The model's capabilities make it particularly valuable for a wide range of applications, from small-scale stream dynamics to large-scale flood simulations. It is design to be  applicable in scenarios involving rapid flood propagation in urban areas, detailed floodplain mapping, and the analysis of complex river morphodynamics in mountainous regions. By integrating RiverFlowDynamics into the Landlab framework, we provide researchers, students, and practitioners with a powerful, accessible tool for hydraulics modeling. This integration facilitates future modifications and contributions, ensuring the model's continual evolution to meet emerging challenges in river system analysis and flood risk assessment.

# Acknowledgements
Funding for this research was provided by Chilean National Agency for Research and Development – ANID though the programme FONDECYT Iniciación grant 11200949. Landlab is supported by the National Science Foundation (NSF Award Numbers 1147454, 1148305, 1450409, 1450338, and 1450412) and by the Community Surface Dynamics Modeling System (NSF Award Numbers 1226297 and 1831623).

# References
