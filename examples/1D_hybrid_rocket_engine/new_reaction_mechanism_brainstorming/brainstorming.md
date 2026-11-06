Trying to come up with a broadly applicable and fast major-species/equilibrium closure applicable for aerospace

Firstly, I believe that a semi-empirical estimate for the turbulent kinetic energy and specific dissipation rate is required

There's a correlation for turbulent kinetic energy which is:

k = (3 / 2) * (turbulence_intensity * velocity)^2
turbulence_intensity = 0.16 * Re^(-1/8)

Then we can use two empirical correlations for specific dissipation rate using another empirical correlation for turbulent length scale

#NOTE: empirical_turbulence_constant is usually set to 0.09
omega = (k^1/2) / (empirical_turbulence_constant^1/4 * L)

L = 0.07 * Dh

then we can use a function similar to this:
equilibrium_turbulence_limited_relaxation_factor = turbulence_limited_reaction_rate_mix_parameter / turbulence_mixing_timescale

Finally, I think the equation structure should look something like this

NOTE: 
    Notice that the chemical rate has no activation energy or R_gas constant. That's intentional, I think this should be purely empirical, since the activation energy for a mixture slowly approaching equilibrium doesn't seem that physical.

    Perhaps this could be modified slightly where each region in a model has a different rate determining step with the global mixing rate applied to it.

    Another approach would be to have a neural network predict *just* the chemical rate.

    Another approahc would be to have the chemical rate be dependent on pressure, temperature, and three empirical constants

Two options for chemical rate:
1. chemical_rate = (empirical_pre_exponential_factor * exp(empirical_exponential_factor * temperature))
2. chemical_rate = (empirical_pre_exponential_factor * temperature^empirical_exponential_factor * exp(-emprical_activation_temperature / temperature))

I'm not a big fan of the second one
I think adding another fitted parameter to the chemical rate is risky, I think if we wanted to add another empirical constant, we would likely use some pressure or density dependence

mixing_rate = equilibrium_turbulence_limited_relaxation_factor

global_relaxation_rate = 1 / ((1 / chemical_rate) + (1 / mixing_rate))

then we just apply it like so: du.parameter -= global_relaxation_rate * (current_parameter - equilibrium_parameter)

This formulation thus requires these empirical correlations strongly (these are unavoidable):
- turbulence_limited_reaction_rate_mix_parameter
- empirical_pre_exponential_factor
- empirical_exponential_factor

These *should* be able to be separated if an experiment can be run where you run different experiments where you vary turbulence independently, and then temperature independently maybe performing a 2x2 experimental matrix

And can generally rely on these without much issue (hopefully)
- empirical_turbulence_constant (typically 0.09)
- some of the turbulence correlations if there isn't a good one avaliable

