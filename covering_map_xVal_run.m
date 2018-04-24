clear all;

% wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
cd(wd);

codeDir = [wd '/code_gridCell'];
saveDir = [wd '/data_gridCell'];
addpath(codeDir); addpath(saveDir);
addpath(genpath([codeDir '/gridSCORE_packed']));

% load
nSet        = 22;
gaussSmooth = 1; 
fixBatchSize = 1; %fixed batch size or depend on nClus (for fname)

dat='circ';
dat='square';
annEps=0;
boxSize = 1;
nIter   = 200;
nSets   = 20; % if just clus positions, no need last two

% joined trials
jointTrls = 1;
clus2run  = [8:2:28]; 
clus2run  = [3:25]; %26 
% clus2run  = [3:16, 18, 20:24];  %batchsize400, no 17,19,25, 26, ... 30

% clus2run  = 3:; %26 


epsMuVals = .025;
nTrials   = 1000000;
batchSizeVals = [1000, 400, 100]; 

% clus2run = 18;
batchSizeVals=400;


locRange = [0 49];
nXvalDataSets=20; %20

nDataPtsTestVals = nTrials/100;
iDataPtsTest = 1;



%new
jointTrlsXval = 1; %xval dat also use jointTrls

%no file

%circ
% 17, 19, 26:30, batchsize=400, anneps=0
% 26:30, batchsize=100, anneps=0

%run 17, 19, 26 for first, just 26 for second

%load loop

for iEps = 1:length(epsMuVals)
    epsMuOrig=epsMuVals(iEps);
    epsMuOrig1000=epsMuOrig*1000;
    for iBvals = 1:length(batchSizeVals)
        batchSize=batchSizeVals(iBvals);
        %load up cluster positions for all sim from 3:26/30 clus, then run
        %xVal on it. %NOTE : this does the xVal on the same xVal data
        %within a shape and batchSize/eps condition
        
        %POTENTIAL PROBLEM - different datasets over batchsizes and
        %learning rates; probably ok for purposes here. if need same
        %datasets, might be worth a function for making the 20 datasets,
        %then function for doing the xval.
        
        for iClus2run = 1:length(clus2run)
            nClus = clus2run(iClus2run);
            fprintf('Loading nClus=%d, epsMu=%d, batchSize=%d\n',nClus,epsMuOrig1000,batchSize)
            %load
            fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wAct_jointTrls_stepSiz',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat)];
            if annEps %epsMu is different here
                fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps*_batchSiz%d_%diters_%s_wAct_jointTrls_stepSiz_annEps',nClus,round(nTrials/1000),batchSize,nIter,dat)];
            end
            fname = [saveDir, fname '*']; %finish with directory and * for date/time
            
            f = dir(fname);
            load(f(1).name); %just load first one
            
            %get cluster positions
            for iterI = 1:nIter 
                %what to do with iSet? -atm just use last set. - later: want to show SSE going down over time?
                for iSet = size(densityPlot,3) %  %1:nSets
                    [mu{iClus2run}(:,1,iterI), mu{iClus2run}(:,2,iterI)] = find(densityPlot(:,:,iSet,iterI)); %find cluster positions - note, only last 2 sets
                end
            end
            
            %log gridness to corr below (maybe not nec do here?)
            gA_gAll(:,:,iEps,iBvals,iClus2run,:)   = gA(:,:,1,:);

        end
        xVal_results = xVal_clusConds(mu, dat,nXvalDataSets, nDataPtsTestVals(iDataPtsTest), locRange, nIter,jointTrlsXval);

         fname = [saveDir sprintf('/covering_map_xVal_%s_%d-%dclus_eps%d_batchSiz%d_%diters_annEps%d_%dktestPts_%dxValDatasets',dat,clus2run(1),clus2run(end),epsMuOrig1000,batchSize,nIter,annEps,nDataPtsTestVals(iDataPtsTest)/1000,nXvalDataSets)];
%          save(fname,'xVal_results');



    end
end

%%
for iClus2run = 1:length(clus2run)
    nClus = clus2run(iClus2run);
    
    [r(iClus2run), p(iClus2run)] = corr(gA_gAll(end,:,1,1,iClus2run)',nanmean(xVal_results.tsseXval(:,:,iClus2run),2),'rows','complete','type','spearman');
    
end

[r;p]'


 











