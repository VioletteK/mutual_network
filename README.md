# How does the analysis work ?
There is always two files for each analysis : the one for the injected cells and the one for spontaneous activity.
In the spontaneous ones, the parameters such as :
- *injection_start*
- *ref_neuron*
- ...                 
must be entered by hand. 

## Simulation codes :

**The time var are in ms**

All the analysises start at the lign `if 'Vm'` ...

### `Analysis_gauss` :

**Definition of the parameters by hand** :
- the *interval* on which we will calculate the MI
- the *number_of_annulus* which is in fact in how many part we divide the disk around the central neuron
- the reference neuron, *ref_neuron*, center of the Annulus
- the ray *r* of the disk around ref_neuron
- the *Time_delay* with when it starts, when it ends and the *step*

 The algorithm calculate the distance between each neuron (except the injected ones) and the central ones
 and assign them an annulus depending on the distance. The empty annulus are deleted.

 Then for each time, calculate the mean MI on each annulus.
 After it, it saves the time when the max of the MI occurs for each annulus and plot it (for the raw and filtered with `savgol_filter` data) with a linear regression.
 The slope of the obtained line is the *speed of diffusion*.

### `Analysis_colormap_vm` :
 It plots the colormap of Vm at each time and can also plot the vm of each neuron as a function of time.

### `Analysis_gradient` : 
 It plots the MI matrix with contour based on filtered data set, the norm of each gradient vector and finally, average MI and average Gradient Norm for each cell. 
 
 ### `Analysis_contour_MI_Vm` :
   It plots the MI matrix with contour based on filtered data set.
   
### `Analysis_delay_2cell` :
  It plots the VM of two cells, their MI and their position on the map. 
  
 ### `Analysis_mean_spikes` :
 **Definition of the parameters by hand** :
- the *time_window* on which we will calculate the mean number of spikes
- the *precision* is the size of the square window on which the mean is computed
- the *Time_delay* with when it starts, when it ends and the *step*
It plots a map of the number of spikes for each neuron in the `time_window`, the mean number of spikes in the square window starting frow each neuron as its low left corner and the Vm. 


## Mouses DATA
The data are savec as a 100x100x511 list (sizexsizextime)
