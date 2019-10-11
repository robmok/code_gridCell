# Clustering Spaces

Code for grid / place cells modelling project


## Software and toolboxes

Matlab (2017b used)

Matalb Toolboxes:
- Statistics and Machine Learning Toolbox (bootci) 
- Signal processing toolbox?
- Image processing toolbox (gaussian filter)

Additional tools/functions:
Tools to compute grid scores: gridSCORE_packed (from Roddy Grieves, based on geom2D toolbox: https://uk.mathworks.com/matlabcentral/fileexchange/7844-geom2d)
Plotting - unvariate scatterplots: https://uk.mathworks.com/matlabcentral/fileexchange/37105-plot-spread-points-beeswarm-plot

## Scripts for simulations

covering_map_batch_run.m
- main run script for running clustering algorithm and simulation

covering_map_batch_sim.m
- called by main run script, covering_map_batch_run.m, - simulation script

run_gridnessTestDataPerm.m
- run script for ‘test’ set - plotting cluster activations after learning, with option to perform shuffling/permutation tests for stats on grid score

gridnessTestData_Perm.m
- called by run_gridnessTestDataPerm.m - test set plotting / permutation testing script


run_trapzKfrmSq_covering_map.m
- loads up a square simulation then runs the clustering algorithm (more learning trials) on the trapezoid

covering_map_batch_sim_clusPosIn.m
- called by run_trapzKfrmSq_covering_map.m - loads in an existing set of cluster positions and runs the clustering algorithm 


## Scripts for plotting
covering_map_batch_plot.m
covering_map_batch_plot_testSet_perm.m
catLearn_batch_plot.m
covering_map_plot_simple_wOut_load.m

## FUNCTIONS

bootrm.m - bootstrap confidence intervals for mean and percent/proportion
createTrls.m - creates trials (x-y coordinates) for different spatial environments
nanconv.m - convolution in 1D or 2D ignoring NaNs; used for smoothing - https://uk.mathworks.com/matlabcentral/fileexchange/41961-nanconv

# Functions for univariate scatterplots  - https://uk.mathworks.com/matlabcentral/fileexchange/37105-plot-spread-points-beeswarm-plot
isEven.m
myErrorbar.m 
plotSpread.m
repeatEntries.m

#Computing grid scores - within directory: gridSCORE_packed
ndautoCORR.m - computes an autocorrelation or crosscorrelation across two arrays
gridSCORE.m - takes autocorrelation plot and computes grid score - can compute two grid scores, the one used in this project is the ’allen’ method


Within gridSCORE_packed/gridSCORE_dependencies
- from the geom2d toolbox: https://uk.mathworks.com/matlabcentral/fileexchange/7844-geom2d


Unused but could be used
kmplusInit.m
