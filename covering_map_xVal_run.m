

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
% dat='square';
annEps=0;
boxSize = 1;
nIter   = 200;
nSets   = 20; % if just clus positions, no need last two

% joined trials
jointTrls = 1;
clus2run  = [8:2:28]; 
epsMuVals = .025;
nTrials   = 1000000;
batchSizeVals = [1000, 400, 100]; 

clus2run = 16;
batchSizeVals=400;


locRange = [0 49];
nXvalDataSets=2; %20

nDataPtsTestVals = nTrials/10;
iDataPtsTest = 1;



%load loop
for iClus2run = 1:length(clus2run) 
    nClus = clus2run(iClus2run);
    for iEps = 1:length(epsMuVals) 
        epsMuOrig=epsMuVals(iEps);
        epsMuOrig1000=epsMuOrig*1000;
        for iBvals = 1:length(batchSizeVals)
            batchSize=batchSizeVals(iBvals);
            fprintf('Loading nClus=%d, epsMu=%d, batchSize=%d\n',nClus,epsMuOrig1000,batchSize)
            
            fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wAct_jointTrls_stepSiz',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat)];
            if annEps %epsMu is different here
                fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps*_batchSiz%d_%diters_%s_wAct_jointTrls_stepSiz_annEps',nClus,round(nTrials/1000),batchSize,nIter,dat)];
            end
            fname = [saveDir, fname '*']; %finish with directory and * for date/time
            
            load(fname);
            
            %need cluster positions
            for iterI = 1:nIter
                for iSet = 1%:nSets
%                     [mu{iClus2run}(:,1,iSet,iterI,iEps,iBvals), mu{iClus2run}(:,2,iSet,iterI,iEps,iBvals)] = find(densityPlot(:,:,iSet,iterI)); %find cluster positions - note, only last 2 sets
                    

                    %what to do with iSet? input to xVal needs to be
                    %nClusx2xiterI, so prob don't put iEPSand iBvals here
                    %either
                    
                    
                    
                    [mu{iClus2run}(:,1,iterI,iEps,iBvals), mu{iClus2run}(:,2,iterI,iEps,iBvals)] = find(densityPlot(:,:,iterI)); %find cluster positions - note, only last 2 sets
                    
                    
                    
                    
                    
                    %NOTE : this does the xVal on the same data - good within dat
                    %type, but different for other shapes. i suppose data
                    %would be different with diff shapes too..
                    
                    %think if this is what i want. probably?
                    
                end
            end
            
            
            xVal_results = xVal_clusConds(mu, dat,nXvalDataSets, nDataPtsTestVals(iDataPtsTest), locRange, nIter);

            
            
        end
    end
end


















