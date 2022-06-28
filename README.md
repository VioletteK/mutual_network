# How does the analysis work ?

Starting at the lign `if 'Vm'` ...



**Definition of the parameters by hand** :
- the *interval* on which we will calculate the MI
- the *number_of_annulus* which is in fact in how many part we divide the disk around the central neuron
- the reference neuron, *ref_neuron*, center of the Annulus
- the ray *r* of the disk around ref_neuron

 The algorithm calculate the distance between each neuron (except the injected ones) and the central ones
 and assign them an annulus depending on the distance. The empty annulus are deleted.

 Then for each time, calculate the mean MI on each annulus.
 After it, it saves the time when the max of the MI occurs for each annulus and plot it (for the raw and filtered with `savgol_filter` data) with a linear regression.
 The slope of the obtained line is the *speed of diffusion*.
