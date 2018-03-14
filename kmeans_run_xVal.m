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
dat='circ'; %square, circ
nKmeans = 1000;
% nPoints = 5000; %3k, 5k, 10k
nPointsVals = [5000, 10000]; % also did 3000 before

fname = [saveDir, sprintf('/kmeans_nK_%d-%d_%s_nPoints%d_%diters',kVals(1),kVals(end),dat,nPoints,nKmeans)];
load(fname);

%xVal specs
saveDat = 1;
nXvalDataSets = 20; % specify how many datasets to generate
% nDataPtsTest = nPoints;
nDataPtsTestVals = [5000, 10000]; %run over number of test points (e.g. compare trained of 3k, 5k, and 10k using SSE on 3k pts)
    
%%
%need to edit xVal_clus to take in all clus positions to assess on the same
%xValdatasets? 
% but if running all clus conditions and iterations merged, should be fine
% - just takes longer

for iNpts = 1:length(nPointsVals)
    nPoints = nPointsVals(iNpts);
    for iDataPtsTest = 1:length(nDataPtsTestVals)
        %run xval
        xVal_results = xVal_clusConds(muAllkVals, dat,nXvalDataSets, nDataPtsTestVals(iDataPtsTest), locRange, nKmeans);
        
        if saveDat
            %         fname = [saveDir sprintf('/kmeans_xVal_nK_%d-%d_%s_%dtrainPts_%dtestPts_%diters_%ddatasets',kVals(1),kVals(end),dat,nPoints, nDataPtsTestVals(iDataPtsTest),nKmeans,nXvalDataSets)];
            fname = [saveDir sprintf('/kmeans_xVal_nK_%d-%d_%s_%dtrainPts_%dtestPts_%diters_%ddatasets_2',kVals(1),kVals(end),dat,nPoints, nDataPtsTestVals(iDataPtsTest),nKmeans,nXvalDataSets)];
            save(fname,'xVal_results');
        end
    end
end