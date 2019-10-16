%% Clustering Spaces (with batch update)
clear all;

% Set working directory
wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
% wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';

cd(wd);
codeDir = [wd '/code_gridCell']; %where the code lives
saveDir = [wd '/data_gridCell']; %where to save the output of the simulations
addpath(codeDir); addpath(saveDir);
addpath(genpath([codeDir '/gridSCORE_packed'])); % add path to code for computing grid measures

%%%%%%%%
%SET UP%
%%%%%%%%
%define box / environment -  circ, square, or cat (cat learning; cat = category learning in a 2D feature space)
dat = 'circ';
% dat = 'square';   
% dat = 'catLearn';

% set number of clusters to run (set to one value during testing, e.g. 20 or a few, 20:22)
clus2run = 20; 

%compute activation (densityPlotActNorm) over training time - takes longer
%(only required for inspecting gridness over time)
actOverTime = 1; 

%annealed learning rate: standard = 1
annEps = 1; %1 or 0

% learning rate
if annEps
    epsMuVals = 0.25; %start from a high learning rate, reduce over time
else
    epsMuVals = 0.025; %fixed learning rate
end

% if category learning, specify number of categories (cluster centres) and sigma of the gaussan
catsInfo.nCats = 2;  % 2 categories
sigmaG = [7 0; 0 7]; catsInfo.R = chol(sigmaG); % variance - isotropic
catsInfo.msExample = 1; % set 2 gaussians in opposite sides of the square - example for ms. else set to 0 for random 2 gaussians.

% number of learning/training trials
if ~strcmp(dat(1:3),'cat')
    nTrials = 1000000;
else
    nTrials = 50000; %cat learning - less trials
end

%batch size
nBatches = 5000; % batchSiz=200; % nBatches = 2500; % batchSiz=400
if strcmp(dat(1:3),'cat')
    nBatches = nBatches.*2;
end
batchSizeVals = nTrials./nBatches;
nBvals = length(batchSizeVals);

% environment
nSteps = 50; %to define spacing beween each loc in box
locRange = [0, nSteps-1]; %[-1, 1]; % from locRange(1) to locRange(2)
stepSize = diff(linspace(locRange(1),locRange(2),nSteps)); stepSize=stepSize(1); %smallest diff between locs
jointTrls = 1; % simulating agent moving around environment

% use the same training data (trials) across current sims or gen new data
useSameTrls=0; %set to 0
%%
saveDat = 0; %save simulations

nIter = 1;%1000; %200 for doing activations over time, 1000 for full simulations without activations over time (variable size gets large), 20-50 for catLearn demo

if useSameTrls
    trials=[]; trialsUnique=[];
else
    trials=[]; 
end

tic
for iClus2run = 1:length(clus2run) %nClus conditions to run
    nClus = clus2run(iClus2run);
    for iEps = 1:length(epsMuVals) 
        epsMuOrig=epsMuVals(iEps);
        epsMuOrig1000=epsMuOrig*1000; %for saving
        
        for iBvals = 1:nBvals
            batchSize = batchSizeVals(iBvals);
            fprintf('Running %s, nClus=%d, epsMu=%d, batchSize=%d\n',dat,nClus,epsMuOrig1000, batchSize)
            fname = [saveDir, sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm',nClus,round(nTrials/1000),epsMuOrig1000,round(batchSize),nIter,dat)];
            
            %run
            tic
            if strcmp(dat(1:3),'cat') %save muAll
                [densityPlot,densityPlotActNorm,gA,gW,gA_actNorm,gW_actNorm,muInit,rSeed,muAll,trials] = covering_map_batch_sim(nClus,locRange,catsInfo,epsMuOrig,nTrials,batchSize,nIter,trials,useSameTrls,dat,annEps,jointTrls,actOverTime);
            else
                [densityPlot,densityPlotActNorm,gA,gW,gA_actNorm,gW_actNorm,muInit,rSeed] = covering_map_batch_sim(nClus,locRange,catsInfo,epsMuOrig,nTrials,batchSize,nIter,trials,useSameTrls,dat,annEps,jointTrls,actOverTime);
                %if testing a single sim, could run replace line above with line below to include 'muAll' and 'trials' - so can plot a single sim using: covering_map_plot_simple_wOut_load.m
%                 [densityPlot,densityPlotActNorm,gA,gW,gA_actNorm,gW_actNorm,muInit,rSeed, muAll,trials] = covering_map_batch_sim(nClus,locRange,catsInfo,epsMuOrig,nTrials,batchSize,nIter,trials,useSameTrls,dat,annEps,jointTrls,actOverTime);
            end
            timeTaken=toc;
            
            if saveDat
                if ~strcmp(dat(1:3),'cat')
                    if useSameTrls
                        fname = [fname '_useSameTrls'];
                    end
                    if jointTrls
                        fname = [fname '_jointTrls_stepSiz'];
                    end
                    if ~actOverTime && ~strcmp(dat(1:3),'cat')
                        fname = [fname '_noActOverTime'];
                    end
                    if annEps
                        fname = [fname '_annEps'];
                    end
                else
                    fname = [fname sprintf('_%dcats',catsInfo.nCats)];
                    if catsInfo.msExample
                       fname = [fname '_msExample']; 
                    end
                end
                cTime=datestr(now,'HHMMSS'); fname = sprintf([fname '_%s'],cTime);
                if strcmp(dat(1:4),'trap')
                    save(fname,'densityPlot','gA','gW','rSeed','muInit','timeTaken'); 
                elseif strcmp(dat(1:3),'cat')
                    save(fname,'muAll','densityPlot','densityPlotActNorm','rSeed','muInit','timeTaken');
                else 
                    save(fname,'densityPlot','densityPlotActNorm','gA','gW','gA_actNorm','gW_actNorm','rSeed','muInit','clusDistB','timeTaken'); 
                end
            end
        end
    end
end
toc
