# Documentation for Better Weatherer (aka ExponentialWeathererIntegrated)


## 1. Why?

Using an exponential weathering function, areas with fast weathering rates or lots of bare bedrock are subject to substantial inaccuaracy in the most common (forward-Eulerian) time stepping implementation, which is a linear extrapolation of present weathering rates through the next timestep. This component uses an analytical integration of the weathering function to exactly determine the weathering total over an arbitrary timestep. 

## 2. Derivation and background
The common exponential weathering function (Anhert, 1970) computes soil production rate that declines as an exponential function of soil thickness, a model that has signficant observational support at least in some geological settings (Heimsath later review). 

Consider that $w_0$ is the maximum soil production rate, ocurring where soil thickness is zero, and that $d^*$ is the characteristic soil production depth. The soil production rate $w$ is given as a function of the soil depth $d$,

$$
    w = w_0 \exp{-\frac{d}{d^*}}   	
$$
	

A common implementation in landscape models is to use the current depth of regolith, apply the above equation to find the erosion rate at each cell, and multiply by a timestep to calculate the amount of soil produced during that timestep. For example, consider the mathematically simple case:

$$
	w=1; \ \ dt=1
$$

The weathering rate $w$ $[L/T]$ multiplied by timestep $dt$ $[T]$ here would imply 1 unit $[L]$ of regolith is produced at the end of that timestep.

Consider, however, that during the timestep, the thickness of soil is increasing, which is decreasing the weathering rate. So the weathering rate at the end of a timestep will be somewhat lower than that at the beginning of the timestep, because of Eq (1). How much lower? Well, a nifty property of exponentials is that they have the same geometry no matter where you are on the curve. What that means is that the final weathering rate in a timestep can be calculated as a fraction of the current weathering rate, regardless of whether the current weathering rate is the maximum $w_0$ precribed for bare bedrock. In other words, there is an exact solution to the amount of regolith produced over any period of time given the current depth and weathering rate. That solution looks like this:

$$
	H_{prod}=d^* log \left[\frac{w*dt}{d^*} + 1\right]
$$

When we compute the thickness of soil production this way, we can see that the linear timestep extrapolation we used above results in about a 30% excess of regolith produced over a single timestep:

```  matlab
>> rstar = 1; rdot_H=1; dt=1;
>> H_prod = rstar .* log(((rdot_H .* dt) ./ rstar) + 1) 
H_prod =
    0.6931
```

## 3. Some nuance

That seems like a dramatic difference, but there is some more nuance here. In the initial scheme, we would end up with a bit more regolith than expected after one time step, which would clamp down the weathering rate for the next timestep, and with small enough timesteps or slow enough weathering rates, or slow erosion rates removing regolith, the thickness would stablize over many timesteps and remain reasonably accurate. However, with large time steps, areas with fast weathering rates or lots of bare bedrock are subject to the most inaccuracy, and these are situations which are commonly of great interest.


## 4. A consideration of density

The above assumes that the regolith being produced is the same as the thickness of bedrock lost; ie that they have equal density. Typically one might compute the thickness of weathered bedrock then multiply by the density ratio to get the additional thickness of regolith. When using $f_{exp}=\rho_{rock}/\rho_{reg}=1.3$, we can see that under the linear timestep extrapolation, the resulting regolith produced is 1.3 units, and it doesn't matter whether this multiplication is done during or after the timestep. However, as the instantaneous thickness feeds back on the production rate, in our integrated solution the expansion factor must be included in the production rate term:

$$
	H_{prod}=d^* log \left[\frac{f_{exp}*w*dt}{d^*} + 1\right]
$$

Using the values above, we see that more regolith is produced over the timestep than in the previous example, because of the expansion factor, but is now 35% less than predicted by simple dt extrapolation corrected for density.

``` matlab
>> f_expand = 1.3
>> H_prod = rstar .* log(((f_expand .* rdot_H .* dt) ./ rstar) + 1)
H_prod =
    0.8329
```

We can confirm that the density ratio must be included in the integral. If we used our integrated result with no density correction (from section 2 above), then post-hoc multiplied by the expansion factor, we overpredict regolith prodiction by about 10%, because the thickness increase caused by expansion isn't included with that of production when integrating the production rates.

``` matlab
>> f_expand * 0.6931
ans =
    0.9010
```

A final note: The production function above integrates assuming we are tracking the amount of regolith produced, expansion included. To determine the final depth of bedrock loss this represents, we may simply divide the final produced depth by the expansion ratio.


---
## 5. Component design
The above discussion informs our component design. In particular, it is necessary to now provide the regolith expansion factor in addition to the soil depth field, maximum production rate, and depth scale. And it must return the integrated depth of soil produced. It may also return the converted depth of bedrock weathered, and the instantaneous weathering rate based on the depth passed in. (This last would be to maintain drop-in compatibility with the existing component). 
```
        Parameters
        ----------
        grid: ModelGrid
            Landlab ModelGrid object
        soil_production__maximum_rate : float
            Maximum weathering rate for bare bedrock
        soil_production__decay_depth : float
            Characteristic weathering depth
        soil_production__expansion_factor : float
            Expansion ratio of regolith (from rel densities of rock and reg) 
```

We are for now keeping the only public attribute the `max_weathering_rate` as in the existing component.
An architectural alternative to integrating the two components directly would be to make this a subclass of the existing exp weatherer with only the new functionality added.
This was dismissed because the integration used explicitly assumes that a particualr weathering function was used for the base rate, and if the original component changes for some reason, it could invalidate the use of the derived component.

Future work should adapt the parameters to be variable across the grid, wheres now they're scalars (but a grid might work... needs testing).

