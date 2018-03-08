%% run xval over kmeans results

clear all;
% wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
cd(wd);
codeDir = [wd '/code_gridCell'];
saveDir = [wd '/data_gridCell'];
addpath(codeDir); addpath(saveDir);

locRange = [0, 49]; %box dims

%load
kVals = 3:30; 
nKvals = length(kVals);
dat='square'; %square, circ
nKmeans = 1000;
nPoints = 10000; %3k, 5k, 10k
fname = [saveDir, sprintf('/kmeans_nK_%d-%d_%s_nPoints%d_%diters',kVals(1),kVals(end),dat,nPoints,nKmeans)];
load(fname);

%xVal specs
saveDat = 1;
nXvalDataSets = 20; % specify how many datasets to generate
% nDataPtsTest = nPoints;
nDataPtsTestVals = [3000, 5000, 10000]; %run over number of test points (e.g. compare trained of 3k, 5k, and 10k using SSE on 3k pts)
    
%%

%need to edit xVal_clus to take in all clus positions to assess on the same
%xValdatasets


for iDataPtsTest = 1:length(nDataPtsTestVals) 
    %run xval
    xVal_results = xVal_clusConds(muAllkVals, dat,nXvalDataSets, nDataPtsTestVals(iDataPtsTest), locRange, nKmeans);

    if saveDat
        fname = [saveDir sprintf('/kmeans_xVal_nK_%d-%d_%s_%dtrainPts_%dtestPts_%diters_%ddatasets',kVals(1),kVals(end),dat,nPoints, nDataPtsTestVals(iDataPtsTest),nKmeans,nXvalDataSets)];
        save(fname,'xVal_results');
    end
    
end