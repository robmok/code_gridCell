# Clustering Spaces

Matlab code for grid / place cells modelling project


## Software and toolboxes

Matlab (2017b used)

Matalb Toolboxes:
- Statistics and Machine Learning Toolbox (bootci) 
- Signal processing toolbox?
- Image processing toolbox (gaussian filter)

Additional tools/functions:
Tools to compute grid scores: gridSCORE_packed (from Roddy Grieves, based on geom2D toolbox: https://uk.mathworks.com/matlabcentral/fileexchange/7844-geom2d)
Plotting - unvariate scatterplots: https://uk.mathworks.com/matlabcentral/fileexchange/37105-plot-spread-points-beeswarm-plot


## Set up directories and paths

Set up directories:

Edit directory path at the top of each scripts:
- covering_map_batch_run.m 
- covering_map_plot_simple_wOut_load
- covering_map_batch_plot.m
- catLearn_batch_plot.m
- run_gridnessTestDataPerm.m
- run_trapzKfrmSq_covering_map
- square_splitInHalf_gridness.m

e.g. for covering_map_batch_run.m, set to working directory by editing line 5:
wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';

Directories:

Move code directory (code_gridCell) to the working directory for the code 

Create an empty data directory in the working directory for data to live (data_gridCell)

## Examples 
### Simple example to test things are working

Script: covering_map_batch_run.m (to set up learning phase simulations)
- Run one condition - Edit line 23 to only run one condition to e.g.: clus2run = 20
- Run script (runs 1 iteration)
- Plot using a simplified plotting script: covering_map_plot_simple_wOut_load.m
- run first 3 cells... (ones below require muAll...)

### Full run through of learning plus testing and permutation stats

Learning (or 'training') phase - learning cluster positions
Script: covering_map_batch_run.m
- default, circle. uncomment line 19 or 20 for square environment or category learning
- default, runs cluster conditions from 10 to 30 (set to fewer just to look at a few results)
- if running a large number of iterations, set actOverTime = 0 (this set to 1 gives activation plots over time; takes longer, creates huge files, might crash. In the paper this is set to 200 iterations)
- edit number of iterations on line 70 - nIter = 1, but change to some larger number, e.g. 50 or 100, or 1000 as in the paper
- save data - edit line 69 to: saveDat = 1;

Test phase - fix cluster positions and compute activations and grid scores; default, circle
Script: run_gridnessTestDataPerm.m
- default, circle. uncomment line 19 if previously ran square environment (category learning no need to assess grid scores)
- if you have run the above and saved it, you should be able to run this script and it will give you the activation maps and grid scores)

Note: gA is grid score that is reported in the manuscript (method from Perez-Escobar et al., 2016), gW is a more conservative method (Wills et al., 2012).


Plotting learning phase
Script: covering_map_batch_plot.m

Plotting test phase

Script: covering_map_batch_plot_testSet_perm.m




### Trapezoid learning and stats









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
