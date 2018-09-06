%% run covering map algorithm with a batch update

clear all;
% close all;

wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
% wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
% wd='/home/robmok/Documents/Grid_cell_model'; %on love01

cd(wd);
codeDir = [wd '/code_gridCell'];
saveDir = [wd '/data_gridCell'];
addpath(codeDir); addpath(saveDir);
addpath(genpath([codeDir '/gridSCORE_packed'])); % ****note edited this - in codeDir now not wd

%define box / environment - random points in a box
dat = 'circ'; % square, circ, rect, or cat (cat learning)cat = category learning in a 2D feature space
% dat = 'square';   
% dat = 'catLearn';

%compute activation and densityPlotActNorm over time? takes longer
actOverTime = 1; 

%annealed learning rate
annEps = 1; %1 or 0

jointTrls = 1;
boxSize = 1; % 1=normal, 2=double size, 3=triple size

% if cat learning specify number of categories (cluster centres) and sigma of the gaussan
catsInfo.nCats=2; %2 categories
sigmaG = [7 0; 0 7]; % variance - isotropic
catsInfo.R=chol(sigmaG);
catsInfo.msExample = 1; % set 2 gaussians in opposite sides of the square - example for ms

% set number of clusters to run
clus2run = 10:30;

% number of learning/training trials
if ~strcmp(dat(1:3),'cat')
    nTrials = 1000000;
%     nTrials = 200000;
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

% parameters
epsMuVals = 0.025;
if annEps
    epsMuVals = 0.25;
end

% use the same training data (trials) across current sims or gen new data
useSameTrls=0;

%box
nSteps = 50; %to define spacing beween each loc in box
locRange = [0, nSteps-1]; %[-1, 1]; % from locRange(1) to locRange(2)
stepSize=diff(linspace(locRange(1),locRange(2),nSteps)); stepSize=stepSize(1); %smallest diff between locs

%%
saveDat=1; %save simulations

nIter=1;%1000; %200 for doing activations over time, 1000 for full simulations without activations over time (variable size gets large), 20-50 for catLearn demo

if useSameTrls
    trials=[]; trialsUnique=[];
else
    trials=[]; %trialsUnique=[];
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
                [densityPlot,densityPlotActNorm,gA,gW,gA_actNorm,gW_actNorm,muInit,rSeed,clusDistB,muAll,trials] = covering_map_batch_sim(nClus,locRange,catsInfo,epsMuOrig,nTrials,batchSize,nIter,trials,useSameTrls,dat,annEps,jointTrls,actOverTime);
            else
                [densityPlot,densityPlotActNorm,gA,gW,gA_actNorm,gW_actNorm,muInit,rSeed,clusDistB] = covering_map_batch_sim(nClus,locRange,catsInfo,epsMuOrig,nTrials,batchSize,nIter,trials,useSameTrls,dat,annEps,jointTrls,actOverTime);
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
                    save(fname,'densityPlot','gA','gW','rSeed','muInit','clusDistB','timeTaken'); 
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
