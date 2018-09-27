# gridCell_code
# Code for grid / place cells simulation project

Matlab toolbox requirements
- Statistics? Signal processing toolbox?

#########
#SCRIPTS#
#########

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


##################
#PLOTTING SCRIPTS#
##################

covering_map_batch_plot.m
covering_map_batch_plot_testSet_perm.m
catLearn_batch_plot.m
covering_map_plot_simple_wOut_load.m

###########
#FUNCTIONS#
###########

createTrls.m - creates trials (x-y coordinates) for different spatial environments

bootrm.m - bootstrap confidence intervals for mean and percent/proportion


plotSpread.m
isEven.m
myErrorbar.m
repeatEntries.m
nanconv.m


Computing grid scores - within directory: gridSCORE_packed
gridSCORE.m
ndautoCORR.m

Within gridSCORE_packed/gridSCORE_dependencies
- from the geom2d toolbox: https://uk.mathworks.com/matlabcentral/fileexchange/7844-geom2d



Unused but could be used
kmplusInit.m
