%% run covering map algorithm with a batch update

clear all;
% close all;

wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
% wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
wd='/home/robmok/Documents/Grid_cell_model'; %on love01

cd(wd);
codeDir = [wd '/code_gridCell'];
saveDir = [wd '/data_gridCell'];
addpath(codeDir); addpath(saveDir);
addpath(genpath([codeDir '/gridSCORE_packed'])); % ****note edited this - in codeDir now not wd

%define box / environment - random points in a box
dat = 'circ'; % square, circ, rect, or cat (cat learning)cat = category learning in a 2D feature space
% dat = 'square'; 
% dat = 'trapzKrupic';

% dat = 'catLearn';

%compute activation and densityPlotActNorm over time? takes longer
actOverTime = 1; 

%annealed learning rate
annEps = 1; %1 or 0

jointTrls = 1;

boxSize = 1; % 1=normal, 2=double size, 3=triple size

% if cat learning specify number of categories (cluster centres) and sigma of the gaussan
catsInfo.nCats=4; %2 categories
% sigmaG = [5 0; 0 5];   % isotropic % sigmaG = [1 .5; .5 2]; R = chol(sigmaG);  % non-isotropic
sigmaG = [3 0; 0 3];
catsInfo.R=chol(sigmaG);

%circ annEps
clus2run = [18, 8, 21, 6, 12]; 
clus2run = [16, 7, 14, 15, 9];  
clus2run = [24, 5, 13, 26, 20];

%love01 - circ annEps, sq annEps
clus2run  = [19,22]; 
%  clus2run = [11, 17];
%  clus2run = [25, 10];
% clus2run  = [3, 23, 4];


%next - batchSize=200

% clus2run=14;

% nTrials
if ~strcmp(dat(1:3),'cat')
    nTrials = 1000000;
else
    nTrials = 50000; %cat learning - need less trials
end

%batch size
fixBatchSize = 1; %fixed, or batchSize depends on mean updates per cluster

if fixBatchSize
    % for 100k trials, joint trials
%     nBatches = [5000, 2500, 10000]; 
%     nBatches = [8000 5000];
%     nBatches = [1000 5000 2500];
    nBatches = 2500;
    if strcmp(dat(1:3),'cat')
        nBatches = nBatches.*2;
    end
    batchSizeVals = nTrials./nBatches;
    nBvals = length(batchSizeVals);
else % define batch size based on average number of updates per cluster
%     avgBatchUpdate = [10, 25, 35, 50]; % 
    avgBatchUpdate = [1, 2, 5]; % avgBatchUpdate = 25;
    nBvals = length(avgBatchUpdate);
end

% parameters
% epsMuVals=[.01, .05, .075, .1, .2, .3];% %learning rate / starting learning rate 
% epsMuVals = 0.075; 
% epsMuVals = [0.05, 0.025]; 
epsMuVals = 0.025;
% epsMuVals = 0.015; 

% use the same training data (trials) across current sims or gen new data
useSameTrls=0;

%box
nSteps = 50; %to define spacing beween each loc in box
locRange = [0, nSteps-1]; %[-1, 1]; % from locRange(1) to locRange(2)
stepSize=diff(linspace(locRange(1),locRange(2),nSteps)); stepSize=stepSize(1); %smallest diff between locs

if boxSize==1.5
    locRange(2)= ceil(locRange(2)*1.5);
elseif boxSize==2 %double
    locRange(2)= locRange(2)*2;
elseif boxSize==3 %triple
    locRange(2)= locRange(2)*3;
end

%stochastic update
stoch = 0; %1 or 0

% c parameter: larger c = less stochastic over trials, smaller c = more stochastic over trials
cVals = [1/nTrials, 2/nTrials, 5/nTrials, 20/nTrials];
% cVals = [1/nTrials, 2/nTrials];
% cVals = [5/nTrials, 20/nTrials];
nC = length(cVals);
if ~stoch
    nC = 1;
else
    nC = length(cVals);
end
% change box shape during learning rectangle
warpBox = 0; %1 or 0
warpType = 'sq2rect';
%%
saveDat=1; %save simulations

nIter=200; %200 for covering map, 20 for cat

if useSameTrls
%     switch dat
%         case 'square'
%         case 'catLearn'
%         otherwise
            trials=[]; trialsUnique=[];
%     end
else
    trials=[]; %trialsUnique=[];
end

tic
for iClus2run = 1:length(clus2run) %nClus conditions to run
    nClus = clus2run(iClus2run);
    for iEps = 1:length(epsMuVals) 
        epsMuOrig=epsMuVals(iEps);
        epsMuOrig1000=epsMuOrig*1000; %for saving
        
        for iStoch = 1:nC %
            c = cVals(iStoch);
            if ~stoch
                c=0;
            end
            c1=round(c.*100000); % for saving file name
        for iBvals = 1:nBvals
            if annEps
%                 epsMuOrig1000 = nBatches(iBvals)/100;%for saving
                epsMuOrig1000 = nBatches*.6; %temp - defined in the function now
            end
            if fixBatchSize
                batchSize = batchSizeVals(iBvals); %fixed batch size
                fprintf('Running %s, nClus=%d, epsMu=%d, c=%d, batchSize=%d\n',dat,nClus,epsMuOrig1000,c1, batchSize)
%                 fname = [saveDir, sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters',nClus,round(nTrials/1000),epsMuOrig1000,round(batchSize),nIter)];
%                 fname = [saveDir, sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wAct',nClus,round(nTrials/1000),epsMuOrig1000,round(batchSize),nIter,dat)];
                %with actNorm computed correctly (no NaNs)
                fname = [saveDir, sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm',nClus,round(nTrials/1000),epsMuOrig1000,round(batchSize),nIter,dat)];
            else % define batch size based on average number of updates per cluster 
%                 batchSize = clus2run(iClus2run).*avgBatchUpdate(iBvals); % batch size depends on average updates per cluster (depends on nClus cond)
%                 fprintf('Running %s, nClus=%d, epsMu=%d, avgBatchUpd=%d; batchSize=%d\n',dat,nClus,epsMuOrig1000,avgBatchUpdate(iBvals),batchSize)
%                 fname = [saveDir, sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_avgBatch%d_batchSiz%d_%diters',nClus,round(nTrials/1000),epsMuOrig1000,round(avgBatchUpdate(iBvals)),round(batchSize),nIter)];
            end
            tic
%             [densityPlot,densityPlotAct,densityPlotActNorm,clusMu,gA,gW,gA_act,gW_act,gA_actNorm,gW_actNorm,rSeed] = covering_map_batch_sim(nClus,locRange,warpType,epsMuOrig,nTrials,batchSize,nIter,warpBox,alpha,trials,useSameTrls,trialsUnique,stochasticType,c,dat,weightEpsSSE);
            
            if strcmp(dat(1:3),'cat') %save muAll
                [densityPlot,densityPlotActNorm,gA,gW,gA_actNorm,gW_actNorm,muInit,rSeed,clusDistB,muAll,trials] = covering_map_batch_sim(nClus,locRange,catsInfo,warpType,epsMuOrig,nTrials,batchSize,nIter,warpBox,trials,useSameTrls,stoch,c,dat,boxSize,annEps,jointTrls,actOverTime);
            else
                [densityPlot,densityPlotActNorm,gA,gW,gA_actNorm,gW_actNorm,muInit,rSeed,clusDistB,muAll, trials] = covering_map_batch_sim(nClus,locRange,catsInfo,warpType,epsMuOrig,nTrials,batchSize,nIter,warpBox,trials,useSameTrls,stoch,c,dat,boxSize,annEps,jointTrls,actOverTime);
            end
            
            timeTaken=toc;
            if saveDat
                if useSameTrls
                    fname = [fname '_useSameTrls'];
                end
                if warpBox
                    fname = [fname '_warpBox'];
                end
                if boxSize>1
                    fname = [fname sprintf('_boxSizex%d',boxSize)];
                end
                if annEps
                    fname = [fname '_annEps'];
                end
                if jointTrls
                    fname = [fname '_jointTrls_stepSiz'];
                end
                if ~actOverTime && ~strcmp(dat(1:3),'cat')
                    fname = [fname '_noActOverTime'];
                end
                
                if strcmp(dat(1:3),'cat')
                    fname = [fname sprintf('_%dcats_stoch%d_c%d',catsInfo.nCats,stoch,c1)];
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
end
toc
