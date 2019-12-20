.. _time_steps:

Time steps
==========

Anyone who has taken a numerical computing course knows the importance of
time-step size: choose too small a step and your calculation takes forever to
run; choose too big a step and it goes unstable. The details depend, of course,
on the particular equations and the numerical method used to solve them.

Components and time steps
_________________________

Landlab components are responsible for ensuring their own stability. This
means, for example, that if your code asks a component to run for 10 years, the
component needs to make sure 10 years isn't too big a step to handle. If it is,
the component should subdivide the 10 years into smaller steps that are stable
and accurate. For example, if the component determines that 1 year is the
biggest step it can get away with, then it should divide the requested 10-year
run into 10 steps of 1 year each. If the component uses variable-size time
steps, then it should subdivide the 10 years as it needs, and return when it
has used up exactly 10 years.

Even if each component ensures its own stability, it is still possible for a
multi-component model to get into trouble by using "global" time steps that
are too large. The risk arises because during each global time step, you are
effectively de-coupling the components for a small interval of time. For
example, suppose you have one component that calculates solar radiation on a
landscape and another that calculates the resulting evapotranspiration (ET).
When you run the ET component, you assume that the radiation is constant for
the duration of the time step—as if the sun's position in the sky became frozen
for that period of time. That's probably a reasonable assumption if your time
step is 5 minutes long, but obviously problematic if you time step is 12 hours!
The moral of the story: test your model with different global time-step sizes
to identify a time step that is small enough not to impact the solution in a
significant way.

To learn more about time steps, numerical stability, and solution accuracy, a
good source is `Numerical Recipes by Press et al <http://numerical.recipes/>`_.
There are also plenty of textbooks on numerical computing available.
