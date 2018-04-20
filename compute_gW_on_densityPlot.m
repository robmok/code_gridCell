%loads in and computes gW grid measure for those I didn't save gW for

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
boxSize=1;
nIter=200;

% joined trials
jointTrls=1;
% clus2run = [8, 12, 16, 20,24, 28]; 
% clus2run = [8:2:28]; 
clus2run = [3:30]; 
epsMuVals=.025;
nTrials=1000000;
batchSizeVals = [1000, 400, 100]; % 125?
% batchSizeVals=400;
annEps=0;

% dat='trapzKrupic';
% % clus2run = [8:2:28]; 
% clus2run = [12:4:28]; 
% batchSizeVals = 400; %100, 125, 200,400, 1000

rHex=0; %if choose raw 60deg corr values, not gridness

noFile = {}; %log which sims not run
%load loop
for iClus2run = 1:length(clus2run) 
    nClus = clus2run(iClus2run);
    for iEps = 1:length(epsMuVals) 
        epsMuOrig=epsMuVals(iEps);
        epsMuOrig1000=epsMuOrig*1000;
        for iBvals = 1:length(batchSizeVals)
            
            batchSize = batchSizeVals(iBvals);
            fprintf('Loading %s, nClus=%d, epsMu=%d, batchSize=%d\n',dat,nClus,epsMuOrig1000,batchSize)
            fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wAct_jointTrls_stepSiz',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat)];
%             fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wAct_jointTrls_stepSizLR',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat)];

            if boxSize>1
                fname = [fname sprintf('_boxSizex%d',boxSize)];
            end
            
            if annEps %epsMu is different here
                fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps*_batchSiz%d_%diters_%s_wAct_jointTrls_stepSiz_annEps',nClus,round(nTrials/1000),batchSize,nIter,dat)];
            end            
            %finish with directory and * for date/time
            fname = [saveDir, fname '*'];
            
            %edit if want to load more than one file per sim, merge
            f = dir(fname); filesToLoad = cell(1,length(f));
            
            if isempty(f) %if no file, don't load/save - but print a warning
                warning('No file for: %s\n',fname);
                noFile = [noFile, fname];
            elseif ~isempty(f)
                
                for iF = 1%:length(f)
                    filesToLoad{iF} = f(iF).name;
                    load(f(iF).name);
                end
                fprintf('Running: %s\n',f(iF).name);
                gW         = nan(nSet,nIter,9);
                gW_actNorm = nan(nSet,nIter,9);
                for iterI = 1:nIter
                    fprintf('Running Iter %d\n',iterI);
                    for iSet = 1:nSet
                        %get densityplot compute gW and gW_actnorm over sets and iters,
                        densityPlotCentresSm = imgaussfilt(densityPlot(:,:,iSet,iterI),gaussSmooth);
                        [g,gdataA] = gridSCORE(ndautoCORR(densityPlotCentresSm),'wills',0);
                        gW(iSet,iterI,:) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
                        
                        densityPlotCentresSm =  imgaussfilt(densityPlotActNorm(:,:,iSet,iterI),gaussSmooth);
                        [g,gdataA] = gridSCORE(ndautoCORR(densityPlotCentresSm),'wills',0);
                        gW_actNorm(iSet,iterI,:) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
                    end
                end
                
                % save all original variables plus gW and gW_actNorm - use the EXACT same fname!!!
                if ~strcmp(dat(1:4),'trap')
                    save(fname,'densityPlot','densityPlotActNorm','gA','gW','gA_actNorm','gW_actNorm','rSeed','muInit','clusDistB','permPrc','timeTaken');
                end
                
            end
                        
        end
    end
end