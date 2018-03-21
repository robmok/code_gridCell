%% run covering map algorithm with a batch update

clear all;
% close all;

wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
% wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
%wd='/home/robmok/Documents/Grid_cell_model'; %on love01

cd(wd);
codeDir = [wd '/code_gridCell'];
saveDir = [wd '/data_gridCell'];
addpath(codeDir); addpath(saveDir);
addpath(genpath([codeDir '/gridSCORE_packed'])); % ****note edited this - in codeDir now not wd

%define box / environment - random points in a box
dat = 'circ'; % rand or cat; rand = uniform points in a box, cat = category learning in a 2D feature space
% dat = 'square'; 
% dat = 'trapz1'; %square rect, trapz, trapzNorm (without Krupic scaling) trapzSqs, or cat (cat learning)
% dat = 'trapz2';% dat = 'trapz3';
% dat = 'trapzKrupic'; % dat = 'trapzKrupic2'; % dat = 'trapzKrupic3';
% dat = 'trapzScaled1';
% dat = 'trapzScaled2';
%  dat = 'trapzScaled3';
% dat = 'trapzNorm';%not scaled, fit into square

boxSize = 1; % 1=normal, 2=double size, 3=triple size

% if cat learning specify number of categories (cluster centres) and sigma of the gaussan
nCats   = 2; %2 categories
sigmaG = [3 0; 0 3]; R = chol(sigmaG);    % isotropic
% sigmaG = [1 .5; .5 2]; R = chol(sigmaG);  % non-isotropic

%run multiple cluster numbers
% clus2run = 12; %20, 30
% % clus2run = [10, 12]; % [11, 14] 
% clus2run  = [20, 24]; 
%    clus2run = [24, 26]; 
%    clus2run = [16, 28]; 
%    clus2run = [18, 22]; 

%odd numbers, and smaller numbers - ran sq, now circ
%love01    
clus2run = [19, 25, 15];
clus2run = [17, 21, 14];
clus2run = [29, 27, 9];

%love06
%clus2run = [5, 12];
% clus2run = [7, 11];
% clus2run = [8, 10];
% clus2run = [3, 4, 6];

%trapz
% clus2run = [18, 24, 26, 28, 16, 30, 20, 22]; %trapzScaled
% clus2run = [24, 26, 30, 28, 20, 22]; %trapzScaled3 x2 size - note 16, 18 actplots dont work
% clus2run = [28 22 14]; %krupic3

% clus2run=20;

% nTrials = 5000000; %how many locations in the box / trials 
nTrials = 2500000; 

%batch size
fixBatchSize = 1; %fixed, or batchSize depends on mean updates per cluster

% 13, 25, 83, 125, 167, 250, 333, 500, 1000, 2000
if fixBatchSize
%     nBatches = [1250, 2500, 5000, 7500, 10000, 15000, 20000];
%     nBatches = [30000, 100000, 200000, 500000, 1250, 2500, 5000, 7500, 10000, 15000, 20000];

% new select batchSizes
    nBatches = [2500, 5000, 20000, 50000];
      nBatches = [2500, 50000];
      nBatches = [5000, 20000];
%     nBatches = 2500;
    batchSizeVals = nTrials./nBatches;
    nBvals = length(batchSizeVals); %length(avgBatchUpdate)
else % define batch size based on average number of updates per cluster
    avgBatchUpdate = [10, 25, 35, 50]; % 
    avgBatchUpdate = [1, 2, 5]; % avgBatchUpdate = 25;
    nBvals = length(avgBatchUpdate);
%     batchSizePerClus = clus2run.*avgBatchUpdate %just to check
    % nBatches = nTrials./batchSizePerClus; %per clus cond %this is not used..
end

% use the same training data (trials) across current sims or gen new data
useSameTrls=0;

%box
nSteps = 50; %to define spacing beween each loc in box
locRange = [0, nSteps-1]; %[-1, 1]; % from locRange(1) to locRange(2)
stepSize=diff(linspace(locRange(1),locRange(2),nSteps)); stepSize=stepSize(1); %smallest diff between locs

if boxSize==2 %double
    locRange(2)= locRange(2)*2;
elseif boxSize==3 %triple
    locRange(2)= locRange(2)*3;
end

% parameters
% epsMuVals=[.01, .05, .075, .1, .2, .3];% %learning rate / starting learning rate 
epsMuVals = 0.075; 

% %tesing 60+clusters
% nTrials = 1000000; 
% nBatches = 1000;
% batchSizeVals = nTrials/nBatches; 
% epsMuVals = 0.025; 

%weight learning rate by SSE 
weightEpsSSE = 0; %1 or 0

% change box shape during learning rectangle
warpBox = 0; %1 or 0
warpType = 'sq2rect';

%mometum-like adaptive learning rate - define alpha (higher = weight
%previous update (direction and magnitude) more; 0 = don't weight previous at all)
alphaVals = 0;
alpha=0;

sTypes = 0;%:1;% :3; %0, 1 ,2, 3
stochasticType=0;
c=0;
% % Create / load in saved test data
% % tile the whole space
% sq=linspace(locRange(1),locRange(2),nSteps);
% allPts=[];
% for i=1:length(sq)
%     for j=1:length(sq)
%         allPts = [allPts; [sq(i), sq(j)]];
%     end
% end
% trials=repmat(allPts,nTrials/length(allPts),1); %note, numel of allPts must be divisble by nTrials atm
% trials=trials(randperm(length(trials)),:);
% % save([saveDir '/randTrialsBox_40k'],'trials');
% trialsUnique=allPts;
% save([saveDir '/randTrialsBox_trialsUnique'],'trialsUnique');
%%
saveDat=1; %save simulations

nIter=200; %how many iterations (starting points)

switch dat
        case 'randUnique'
        %all unique points in box
        load([saveDir '/randTrialsBox_trialsUnique']);
        trials = trialsUnique;
        % does it matter how many points there are if all the same points? e.g.
        % same if just have each trialsUnique twice/x10? - i think not
        % trials = repmat(dataPts,50,1);
    case 'square'
        trials = [randsample(linspace(locRange(1),locRange(2),50),nTrials,'true'); randsample(linspace(locRange(1),locRange(2),50),nTrials,'true')]';
        %for computing sse over trials
%         load([saveDir '/randTrialsBox_trialsUnique']);
        trialsUnique=[];
    case 'cat'
        % draw points from 2 categories (gaussian) from a 2D feature space
        nTrials = floor(nTrials/nCats); % points to sample
        for iCat = 1:nCats 
            mu(iCat,:)=randsample(locRange(1)+10:locRange(2)-10,2,'true'); % �10 so category centres are not on the edge
            datPtsGauss(:,:,iCat) = round(repmat(mu(iCat,:),nTrials,1) + randn(nTrials,2)*R); % key - these are the coordinates of the points
        end
        trials = reshape(datPtsGauss,nTrials,2);
        trials = trials(randperm(length(trials)),:);
        trialsUnique=[];
    otherwise
        trials=[]; trialsUnique=[];
end

tic
for iClus2run = 1:length(clus2run) %nClus conditions to run
    nClus = clus2run(iClus2run);
    for iEps = 1:length(epsMuVals) 
        epsMuOrig=epsMuVals(iEps);
        epsMuOrig1000=epsMuOrig*1000; %for saving
        for iBvals = 1:nBvals
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
            [densityPlot,densityPlotAct,densityPlotActNorm,clusMu,gA,gW,gA_act,gW_act,gA_actNorm,gW_actNorm,rSeed] = covering_map_batch_sim(nClus,locRange,warpType,epsMuOrig,nTrials,batchSize,nIter,warpBox,alpha,trials,useSameTrls,trialsUnique,stochasticType,c,dat,weightEpsSSE);
            timeTaken=toc;
            if saveDat
                if useSameTrls
                    fname = [fname '_useSameTrls'];
                end
                if warpBox
                    fname = [fname '_warpBox'];
                end
%                 if ~strcmp(dat,'square') %atm, only square not append name; later add to above fname
                 fname = [fname sprintf('_%s',dat)]; % for new sims (from adding activations, save name too)
%                 end
                if boxSize>1
                    fname = [fname sprintf('_boxSizex%d',boxSize)];
                end
                cTime=datestr(now,'HHMMSS'); fname = sprintf([fname '_%s'],cTime);
%                 save(fname,'densityPlot','densityPlotAct','clusMu','gA','gW','gA_act','gW_act','nIter','rSeed','timeTaken'); %added trialsAll for xval - removed, too big.maybe compute at end of each sim? or at each set
                save(fname,'densityPlot','densityPlotAct','densityPlotActNorm','clusMu','gA','gW','gA_act','gW_act','gA_actNorm','gW_actNorm','nIter','rSeed','timeTaken'); %added trialsAll for xval - removed, too big.maybe compute at end of each sim? or at each set
            end
        end
        
    end
end
toc
