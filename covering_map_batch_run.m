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
% dat = 'trapz1'; % trapz, trapzNorm (without Krupic scaling) trapzSqs,
% dat = 'trapz2';% 
% dat = 'trapz3';
% dat = 'trapzKrupic';  
% dat = 'trapzScaled1';
% dat = 'trapzScaled2';
%  dat = 'trapzScaled3';
% dat = 'trapzNorm';%not scaled, fit into square

% dat = 'catLearn';

boxSize = 1; % 1=normal, 2=double size, 3=triple size

% if cat learning specify number of categories (cluster centres) and sigma of the gaussan
catsInfo.nCats=5; %2 categories
sigmaG = [5 0; 0 5];   % isotropic % sigmaG = [1 .5; .5 2]; R = chol(sigmaG);  % non-isotropic
sigmaG = [3 0; 0 3];
catsInfo.R=chol(sigmaG);

%annealed learning rate
annEps = 0; %1 or 0

%perm testing
doPerm = 0;

jointTrls = 1;

% fewer trials, lower learning rate
% clus2run = [16, 20, 18, 30];
% clus2run = [22, 28, 26, 14];

% %joint trials; 8, 12, 16, 20, 24, 28;    then 6, 10, 14, 18, 22, 26, 28
% clus2run = [26, 8,  10]; 
% clus2run = [16, 28, 14,];  
% clus2run = [12, 20, 22]; 
% clus2run = [28, 24, 5, 6];
% clus2run = [3, 4, 30];

%trapzK love06 more - batchSize=100 / nBatches 8k/10k
%odd, 3 x 2 batchSizes, 6 matlabs
% clus2run = [25, 7,9,15, 29];
% clus2run = [11, 19, 21, 29];
% clus2run = [13, 27, 17, 23];

%now runnig trapzK annEps - 4 batchsizes - 3x nClus, 2x batchSizes
% clus2run = [26, 8,  10, 16, 28, 7, 21, 19, 6,  17]; 
% clus2run = [12, 20, 22, 28, 24, 5,  11,  13,  14]; 
% clus2run = [3,  4,  30, 29, 25, 15, 27, 9,  23];


%run the rest and odd number - 3,4,5,6, 7:2:29
% love01 - circ and sq annEps
% clus2run = [3,  29, 25];
% clus2run = [15, 27,  9];
% clus2run = [21, 19,  5];
% clus2run = [6,  11,  17];
% clus2run = [13, 7,  23];
% clus2run = [4,  30];

% batchSize=200 - circ love01 - 3:26
% clus2run = [20,  8,  10]; 
% clus2run = [3,  7,  22]; 
% clus2run = [26,  5,  6];
% clus2run = [15,  12,  4];
% clus2run = [11,  17, 9];
%  clus2run = [16, 24];  
%  clus2run = [13, 23];
%  clus2run = [21, 19];
%  clus2run = [14, 25];
% clus2run=18; %MISSED

 %square love06 - anneps, 4 sets -batchsize 2500, 10000, 5000
clus2run = [8,  22,  9, 11, 25]; 
clus2run = [12, 18, 16, 7,  17];
 clus2run = [20, 3, 26, 15]; %  5000 - first 2 ran, run rest of this

%next run: batchsize - 5000, 2500, 10000
clus2run = [5,  13, 24, 14, 4];
clus2run = [23, 19, 10, 21, 6];

%next
% clus2run = 27;
% clus2run = 28;
% clus2run = 29;
% clus2run = 30;

%circ not annEps missed:
%batch=400
% clus2run = [17,19];
% clus2run = [26];
%batch=100
% clus2run = 26;

% clus2run = 18; % trapzK missed 18 - 4 batchsizes, 2 matlabs for annEps=1/0

clus2run = 10; 


% nTrials = 5000000; %how many locations in the box / trials 
% nTrials = 2000000;
nTrials = 1000000; %new

% nTrials = 50000; %cat learning - need less trials

%batch size
fixBatchSize = 1; %fixed, or batchSize depends on mean updates per cluster

% 13, 25, 83, 125, 167, 250, 333, 500, 1000, 2000
if fixBatchSize
%     nBatches = [30000, 100000, 200000, 500000, 1250, 2500, 5000, 7500, 10000, 15000, 20000];
% new select batchSizes
%     nBatches = [2500, 20000,5000 50000];
    
    %new - for 100k trials, half nBatches for same batchsize
%     nBatches = [20000, 5000, 50000]./2; %half nBatches
    %joint trials
%     nBatches = [1000, 2500, 10000]; 
%     nBatches = [2500, 10000]; 
    
    nBatches = 2500;
%     nBatches = 10000;

    
    
%     nBatches = [2500, 10000, 5000, 8000]; 
%     nBatches = [2500, 10000, 5000]; 
%     nBatches = [5000, 2500, 10000]; 
%     nBatches = [8000 5000];\
%     nBatches = 5000;


    batchSizeVals = nTrials./nBatches;
    nBvals = length(batchSizeVals);
else % define batch size based on average number of updates per cluster
    avgBatchUpdate = [10, 25, 35, 50]; % 
    avgBatchUpdate = [1, 2, 5]; % avgBatchUpdate = 25;
    nBvals = length(avgBatchUpdate);
%     batchSizePerClus = clus2run.*avgBatchUpdate %just to check
    % nBatches = nTrials./batchSizePerClus; %per clus cond %this is not used..
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
% sTypes = 1;%:1;%  %0, 1 ,2, 3
% % 0. none
% % 1. standard stochastic update - becomes more det over time; becomes
% % basically deterministic at some point
% %  larger c = less stochastic over trials (becomes det quite early on); smaller c = more stochastic over trials (still a bit stochastic by the end)
% % cValsOrig = [2/nTrials, 5/nTrials, 10/nTrials, 20/nTrials]; %removed .1/nTrials and .25/nTrials,  too stochastic. also 3/ntrials, .5/nTrials
% cValsOrig = 5/nTrials;


% change box shape during learning rectangle
warpBox = 0; %1 or 0
warpType = 'sq2rect';


sTypes = 0;%:1;% :3; %0, 1 ,2, 3
stochasticType=0;
c=0;
%%
saveDat=1; %save simulations

nIter=200; %how many iterations (starting points)

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
        for iBvals = 1:nBvals
            if annEps
                epsMuOrig1000 = nBatches(iBvals)/100;%for saving
            end
            if fixBatchSize
                batchSize = batchSizeVals(iBvals); %fixed batch size
                fprintf('Running %s, nClus=%d, epsMu=%d, batchSize=%d\n',dat,nClus,epsMuOrig1000,batchSize)
%                 fname = [saveDir, sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters',nClus,round(nTrials/1000),epsMuOrig1000,round(batchSize),nIter)];
                fname = [saveDir, sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wAct',nClus,round(nTrials/1000),epsMuOrig1000,round(batchSize),nIter,dat)];
            else % define batch size based on average number of updates per cluster 
%                 batchSize = clus2run(iClus2run).*avgBatchUpdate(iBvals); % batch size depends on average updates per cluster (depends on nClus cond)
%                 fprintf('Running %s, nClus=%d, epsMu=%d, avgBatchUpd=%d; batchSize=%d\n',dat,nClus,epsMuOrig1000,avgBatchUpdate(iBvals),batchSize)
%                 fname = [saveDir, sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_avgBatch%d_batchSiz%d_%diters',nClus,round(nTrials/1000),epsMuOrig1000,round(avgBatchUpdate(iBvals)),round(batchSize),nIter)];
            end
            tic
%             [densityPlot,densityPlotAct,densityPlotActNorm,clusMu,gA,gW,gA_act,gW_act,gA_actNorm,gW_actNorm,rSeed] = covering_map_batch_sim(nClus,locRange,warpType,epsMuOrig,nTrials,batchSize,nIter,warpBox,alpha,trials,useSameTrls,trialsUnique,stochasticType,c,dat,weightEpsSSE);
            [densityPlot,densityPlotActNorm,gA,gW,gA_actNorm,gW_actNorm,muInit,rSeed,clusDistB, permPrc,muAll, trials] = covering_map_batch_sim(nClus,locRange,catsInfo,warpType,epsMuOrig,nTrials,batchSize,nIter,warpBox,trials,useSameTrls,stochasticType,c,dat,boxSize,annEps,jointTrls,doPerm);

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
                if doPerm
                    fname = [fname '_doPerm'];
                end
                if jointTrls
                    fname = [fname '_jointTrls_stepSiz'];
                end
                
                if strcmp(dat(1:3),'cat')
                    fname = [fname sprintf('_%dcats',catsInfo.nCats)];
                end

                
                cTime=datestr(now,'HHMMSS'); fname = sprintf([fname '_%s'],cTime);
%                 save(fname,'densityPlot','densityPlotAct','clusMu','gA','gW','gA_act','gW_act','nIter','rSeed','timeTaken'); %added trialsAll for xval - removed, too big.maybe compute at end of each sim? or at each set
%                 save(fname,'densityPlot','densityPlotActNorm','gA','gA_actNorm','rSeed','muInit','clusDistB','timeTaken'); 
                
                if ~strcmp(dat(1:4),'trap')
                    save(fname,'densityPlot','densityPlotActNorm','gA','gW','gA_actNorm','gW_actNorm','rSeed','muInit','clusDistB','permPrc','timeTaken'); 
                else
                    save(fname,'densityPlot','gA','gW','rSeed','muInit','clusDistB','permPrc','timeTaken'); 
                end
            end
        end
        
    end
end
toc
