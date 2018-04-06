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
nKmeans = 200;
nSims = 5; % 5, 10
nIters=nKmeans*nSims;
nPointsVals = [10000, 50000]; % also did 3k, 5k before
% nPointsVals = 5000;
% nPointsVals = 10000;
nPointsVals = 50000;

%xVal specs
saveDat = 1;
nXvalDataSets = 20; % specify how many datasets to generate
nDataPtsTestVals = [10000, 50000]; %run over number of test points (e.g. compare trained of 3k, 5k, and 10k using SSE on 3k pts)
    
%% run

for iNpts = 1:length(nPointsVals)
    fname = [saveDir, sprintf('/kmeans_nK_%d-%d_%s_nPoints%d_%diters_%dsims_merged',kVals(1),kVals(end),dat,nPointsVals(iNpts),nKmeans,nSims)];
    load(fname);
    for iDataPtsTest = 1:length(nDataPtsTestVals)
        %run xval
        fprintf('Running %d dataPts, %d testDataPts\n',nPointsVals(iNpts),nDataPtsTestVals(iDataPtsTest));
        xVal_results = xVal_clusConds(muAllkVals, dat,nXvalDataSets, nDataPtsTestVals(iDataPtsTest), locRange, nIters);
        if saveDat
            fname = [saveDir sprintf('/kmeans_xVal_nK_%d-%d_%s_%dtrainPts_%dtestPts_%diters_%sims_%ddatasets',kVals(1),kVals(end),dat,nPointsVals(iNpts),nDataPtsTestVals(iDataPtsTest),nKmeans,nSims,nXvalDataSets)];
            save(fname,'xVal_results');
        end
    end
end