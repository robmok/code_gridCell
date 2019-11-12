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

1. Edit working directory path at the top of each of these scripts :
- this working directory (wd) should have directories code_gridCell and data_gridCell

For running the simulations:
- covering_map_batch_run.m 
- covering_map_plot_simple_wOut_load
- run_gridnessTestDataPerm.m
- run_trapzKfrmSq_covering_map
- square_splitInHalf_gridness.m
For plotting:
- covering_map_batch_plot.m
- catLearn_batch_plot.m

e.g. for covering_map_batch_run.m, set to working directory by editing line 5:
wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';

Directories:

2. Move code directory (code_gridCell) to the working directory for the code 

3. Create an empty data directory in the working directory for data to live (data_gridCell)

- OR: download data directory data_gridCell from osf (https://osf.io/2dz3x/)
	

## Examples 
### Simple example to test things are working

Script: covering_map_batch_run.m (to set up learning phase simulations)
- Run one condition - Edit line 23 to only run one condition to e.g.: clus2run = 20
- Run script (runs 1 iteration)
- Plot using a simplified plotting script: covering_map_plot_simple_wOut_load.m
- Run first 3 cells for plotting 
- (Note: commented out cells below require the 'muAll' variable which allows plotting some visualisations over time. To do this, edit script 'covering_map_batch_run' to output muAll - i.e.  comment line 95 and uncomment line 97)

### Full run through of learning plus testing and permutation stats

#### Learning
Learning (or 'training') phase - learning cluster positions
Script: covering_map_batch_run.m
- default (i.e. in the script now): circle (uncomment line 19 or 20 for square environment or category learning)
- default: runs cluster conditions from 10 to 30 (set to fewer just to look at a few results, e.g. 10:12, or [10, 15, 20])
- if running a large number of iterations, set actOverTime = 0 (this set to 1 gives activation plots over time; takes longer, creates huge files, might crash. In the paper this is set to 200 iterations)
- edit number of iterations on line 70 - nIter = 1, but change to some larger number, e.g. 50 or 100, or 1000 as in the paper
- save data - edit line 69 to: saveDat = 1;

Test phase (freeze cluster positions and compute activations and grid scores)
Script: run_gridnessTestDataPerm.m
- default environment/shape is circle. uncomment line 19 if you ran learning in the square environment (category learning no need to assess grid scores)
- if you have run the learning script and saved it, you should be able to run this script and it will give you the activation maps and grid scores)
- permutation testing (on/off)
	- To compute activation maps and grid scores after learning WITHOUT permutation tests, set nIter=1000
	- To compute activation maps, grid scores, and permutation tests, set nIter=200 (else takes long, and not necessary)
(Note: gA is grid score that is reported in the manuscript (method from Perez-Escobar et al., 2016), gW is a more conservative method (from Wills et al., 2012).)

#### Plotting
Plotting learning phase
Script: covering_map_batch_plot.m

- first cell loads in data (default loads in nClus (number of clusters) conditions 10:30 - edit line 38 if ran less/different nClus conditions)

***TO DO
- General: probably need to edit first section for better readability, so that user will only need to edit a couple lines
- need to fix cells 2 and 3 - plotting / stats over time
- cell 4 works fine (univar scatters at the end)
- cell 5 - density plot examples - mostly ok; not sure why the act maps colorscale isn't correct; autocorr maps are ok
- cell 6 - needs muAll; used to generate fig from fig panel 1 (many clusters, spatial case)
***
- NOTE: atm loading in trapz doesn't work. looks like first problem is loading in the file, since the strings/words are mixed up for the trapz compared to circ/sq. fix or just leave 
- think about whether to include orientation, rad, wav measures, or edit so don't even see them. if remove; probabably save a copy for myself



Plotting test phase
Script: covering_map_batch_plot_testSet_perm.m

**TO DO
- General: DEFINITELY need to edit first section for better readability, so that user will only need to edit a couple lines
- loads in ok
- 2nd cell - produces 2 plots. only keep first univar scatter? also, it spits out a bunch of stats, significance things - maybe remove this stuff that is not reported in the paper
- prop grid cell measures - check which one i used. average across all simulations, or across nClus conds? remove the one i did not use



### Trapezoid learning and stats

Learning (or 'training') phase - learning cluster positions in a trapezoid after learning in a square
Script: run_trapzKfrmSq_covering_map.m
- this loads in the data from the data directory (corresponding nClus condition) and runs the learning algorithm in a trapezoid
 
Test phase - fix cluster positions and compute activations and grid scores
Script: run_gridnessTestDataPerm.m
- uncomment line 20 to make dat = 'trapzKfrmSq1'
- this should automatically set doPerm=0, and nIter=1000, since mainly we want activation maps and grid scores after learning in the trapz

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

## Data files for plotting
### SPATIAL CASE:
For each nClus condition, there are 10 files. here are 10 files for the nClus=10 condition.

Main learning results for circle and square environment:
- covering_map_batch_dat_10clus_1000ktrls_eps250_batchSiz200_1000iters_circ_wActNorm_jointTrls_stepSiz_noActOverTime_annEps_101640.mat
- covering_map_batch_dat_10clus_1000ktrls_eps250_batchSiz200_1000iters_square_wActNorm_jointTrls_stepSiz_noActOverTime_annEps_165924.mat

Main test results (clusters fixed - first 2 without permutation test, latter 2 with perm test)
- covering_map_batch_dat_10clus_1000ktrls_eps250_batchSiz200_1000iters_circ_wActNorm_jointTrls_stepSiz_noActOverTime_annEps_trlsTest_noPerm_232126.mat
- covering_map_batch_dat_10clus_1000ktrls_eps250_batchSiz200_1000iters_square_wActNorm_jointTrls_stepSiz_noActOverTime_annEps_trlsTest_noPerm_043400.mat
- covering_map_batch_dat_10clus_1000ktrls_eps250_batchSiz200_200iters_circ_wActNorm_jointTrls_stepSiz_annEps_actNorm_perm_500permsOn200iters_160040.mat
- covering_map_batch_dat_10clus_1000ktrls_eps250_batchSiz200_200iters_square_wActNorm_jointTrls_stepSiz_annEps_actNorm_perm_500permsOn200iters_053554.mat

Main learning results with activations and grid score saved over time (with fewer - 200 - iterations)
covering_map_batch_dat_10clus_1000ktrls_eps250_batchSiz200_200iters_circ_wActNorm_jointTrls_stepSiz_annEps_142357.mat
covering_map_batch_dat_10clus_1000ktrls_eps250_batchSiz200_200iters_square_wActNorm_jointTrls_stepSiz_annEps_200530.mat

Trapezoid:
- covering_map_batch_dat_10clus_250ktrls_eps250_batchSiz200_1000iters_trapzKfrmSq1_wActNorm_epsMuTrapz_25_jointTrls_stepSiz_annEps_201604.mat
- covering_map_batch_dat_10clus_250ktrls_eps250_batchSiz200_1000iters_trapzKfrmSq1_wActNorm_epsMuTrapz_25_jointTrls_stepSiz_annEps_trlsTest_noPerm_trapzKfrmSq1_184626.mat

### CONCEPT CASE:
For the manuscript, I ran 2 simulations for the concept structure example
- covering_map_batch_dat_18clus_50ktrls_eps25_batchSiz10_50iters_catLearn_wActNorm_2cats_stoch0_c0_msExample_145757.mat
- covering_map_batch_dat_20clus_50ktrls_eps25_batchSiz10_50iters_catLearn_wActNorm_2cats_stoch0_c0_msExample_144044.mat






## Functions

bootrm.m - bootstrap confidence intervals for mean and percent/proportion
createTrls.m - creates trials (x-y coordinates) for different spatial environments
nanconv.m - convolution in 1D or 2D ignoring NaNs; used for smoothing - https://uk.mathworks.com/matlabcentral/fileexchange/41961-nanconv

### Functions for univariate scatterplots  - https://uk.mathworks.com/matlabcentral/fileexchange/37105-plot-spread-points-beeswarm-plot
isEven.m
myErrorbar.m 
plotSpread.m
repeatEntries.m

### Computing grid scores - within directory: gridSCORE_packed
ndautoCORR.m - computes an autocorrelation or crosscorrelation across two arrays
gridSCORE.m - takes autocorrelation plot and computes grid score - can compute two grid scores, the one used in this project is the ’allen’ method


Within gridSCORE_packed/gridSCORE_dependencies
- from the geom2d toolbox: https://uk.mathworks.com/matlabcentral/fileexchange/7844-geom2d


Unused but could be used
kmplusInit.m
