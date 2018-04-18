

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
% batchSizeVals=400;

%load loop
for iClus2run = 1:length(clus2run) 
    nClus = clus2run(iClus2run);
    for iEps = 1:length(epsMuVals) 
        epsMuOrig=epsMuVals(iEps);
        epsMuOrig1000=epsMuOrig*1000;
        for iBvals = 1:length(batchSizeVals)
            fprintf('Loading nClus=%d, epsMu=%d, batchSize=%d\n',nClus,epsMuOrig1000,batchSize)
            
            fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wAct_jointTrls_stepSiz',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat)];
            if annEps %epsMu is different here
                fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps*_batchSiz%d_%diters_%s_wAct_jointTrls_stepSiz_annEps',nClus,round(nTrials/1000),batchSize,nIter,dat)];
            end
            fname = [saveDir, fname '*']; %finish with directory and * for date/time
            
            %need cluster positions
            for iter = 1:nIter
                for iSet = 1:nSets
                    
                    [mu{iClus2Run}(:,1,iSet,iterI,iEps,iBvals), mu{iClus2Run}(:,2,iSet,iterI,iEps,iBvals)] = find(densityPlot(:,:,iSet,iterI)); %find cluster positions - note, only last 2 sets
                    
                    
                    
                    
                    
                end
            end
            
            
        end
    end
end