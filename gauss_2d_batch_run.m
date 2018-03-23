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

dat = 'square'; % rand or cat; rand = uniform points in a box, cat = category learning in a 2D feature space

% if cat learning specify number of categories (cluster centres) and sigma
% of the gaussan
nCats   = 2; %2 categories
sigmaG = [3 0; 0 3]; R = chol(sigmaG);    % isotropic
% sigmaG = [1 .5; .5 2]; R = chol(sigmaG);  % non-isotropic

%run multuple cluster numbers
clus2run = 20; %[10 20]; %20, 30
% nTrials = 40000; %how many locations in the box / trials - 2.5k ; 5k if reset
clus2run = [20, 18];
%clus2run = [16, 28];
%clus2run = [24, 14];

clus2run = 20;

nTrials = 2500000; 

%batch size
fixBatchSize = 1; %fixed, or batchSize depends on mean updates per cluster
% nBatches = [1250, 2500, 5000, 7500, 10000, 15000, 20000];% nBatches = [30000, 100000, 200000, 500000];
nBatches = [1250, 2500, 5000, 10000];
nBatches = 1250;
% nBatches = 1250; %2500
batchSizeVals = nTrials./nBatches;
nBvals = length(batchSizeVals); %length(avgBatchUpdate)
    
%box
nSteps = 50; %to define spacing beween each loc in box
locRange = [0, nSteps-1]; %[-1, 1]; % from locRange(1) to locRange(2)
stepSize=diff(linspace(locRange(1),locRange(2),nSteps)); stepSize=stepSize(1); %smallest diff between locs

% parameters
% learning rate - 
% epsMuVals = [.0075, .01, .02]; 
% epsMuVals = .0075; %atm .01 works well-ish but a bit fast, could be slower; .0075, 005
% epsMuVals = .01;
% with repel
epsMuVals = 0.8;


%weight learning rate by SSE 
weightEpsSSE = 0; %1 or 0

%looks like optimizing the std is also important - stepSize/3.5 is pretty good
% sigmaGaussVals = [stepSize/3, stepSize/3.5, stepSize/4];
sigmaGaussVals = stepSize/3.5;
% sigmaGaussVals = stepSize/3; %with smaller sigmaGauss, need faster learning rate ; e.g. 0.075 was too slow (though ok for sigmaGauss=3.5); 0.01 still slow; 0.05 is gd but fast
% sigmaGaussVals = stepSize/4;


sigmaGaussVals = 5;


%define box / environement - random points in a box
box = 'square'; %square, rect, trapz, trapzSq (trapz and a square box attached)

% change box shape during learning
%rectangle
warpBox = 0; %1 or 0
warpType = 'sq2rect';

%mometum-like adaptive learning rate - define alpha (higher = weight
%previous update (direction and magnitude) more; 0 = don't weight previous at all)
alphaVals = 0;

sTypes = 0;%:1;% :3; %0, 1 ,2, 3
% cValsOrig = [1/nTrials, 5/nTrials, 10/nTrials, 20/nTrials]; %just edited 2/nTrials --> 1/nTrials
% cValsOrig = [2/(nTrials/2), 5/nTrials, 10/(nTrials/2), 20/(nTrials/2)];% if 80k trials...
cValsOrig = 2/nTrials;

%neighbour-weighted update
neigh = 0; %if neigh = 0, no stoch, no alpha
betaVals  = [.3 .5 .7]; %softmax param - higher = less neighbour update; from .25 going toward middle; .3 start OK
if neigh
   sTypes=0;
   alphaVals=0;
   epsMuVals=[.05 .075 .1];
end

%%
saveDat=0; %save simulations

nIter=1; %how many iterations (starting points)

plotGrids = 0; %plot to test? if nIter > 8, then won't plot

switch dat
    case 'square'
        trials = [randsample(linspace(locRange(1),locRange(2),50),nTrials,'true'); randsample(linspace(locRange(1),locRange(2),50),nTrials,'true')]';
        %for computing sse over trials
%         load([saveDir '/randTrialsBox_trialsUnique']);
    case 'cat'
        % draw points from 2 categories (gaussian) from a 2D feature space
        nPoints = floor(nTrials/nCats); % points to sample
        for iCat = 1:nCats
            mu(iCat,:)=randsample(locRange(1)+10:locRange(2)-10,2,'true'); % ±10 so category centres are not on the edge
            datPtsGauss(:,:,iCat) = round(repmat(mu(iCat,:),nPoints,1) + randn(nPoints,2)*R); % key - these are the coordinates of the points
        end
        trials = reshape(datPtsGauss,nTrials,2);
        trials = trials(randperm(length(trials)),:);
        trialsUnique=[];
    case 'circ'
        % Create logical image of a circle
        imageSizeX = nSteps;
        [columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeX);
        centerX = nSteps/2; centerY = nSteps/2;
        radius = nSteps/2-1;
        circIm = (rowsInImage - centerY).^2 ...
            + (columnsInImage - centerX).^2 <= radius.^2;
        circPts=[]; % find circle points in XY coords
        for iX=1:length(circIm)
            yVals = find(circIm(iX,:));
            circPts = [circPts; ones(length(yVals),1)*iX, yVals'];
        end
        trialInd=randi(length(circPts),nTrials,1);
        trials=circPts(trialInd,:);
        %dataPtsTest
        trialIndTest = randi(length(circPts),nTrials,1);
        dataPtsTest  = circPts(trialIndTest,:);
end


tic
if ~neigh %separating neigh and stoch/momentum params
    
    %testing stoch/momentum/eps
    for iClus2run = 1:length(clus2run) %nClus conditions to run
        nClus = clus2run(iClus2run);
%         fprintf('Running clus %d\n',nClus)       
        for iEps = 1:length(epsMuVals)
            epsMuOrig=epsMuVals(iEps);
            epsMuOrig10000=epsMuOrig*10000; %for saving - changed from 1000 to 100000 for slower l rates - changed back now
%             fprintf('Running epsMu %0.4f\n',epsMuOrig)
            for iBvals = 1:nBvals
                batchSize = batchSizeVals(iBvals); %fixed batch size
                for iSigma = 1:length(sigmaGaussVals)
                    sigmaGauss = sigmaGaussVals(iSigma);
                    sigmaGauss100=round(sigmaGauss*100);
%                     fprintf('Running sigmaGauss %0.2f\n',sigmaGauss)
                    for iStype = 1:length(sTypes)
                        stochasticType = sTypes(iStype);
                        if ~stochasticType
                            cVals = 0; else, cVals = cValsOrig;
                        end
                        for iStochastic = 1:length(cVals) %
                            c = cVals(iStochastic);
                            c1m=round(c.*1000000); % for saving file name
%                             fprintf('Running cVal %d\n',c1m)
                            for iAlpha  = 1:length(alphaVals) %
                                alpha = alphaVals(iAlpha);
                                alpha10 = alpha*10; %for saving simulations
%                                 fprintf('Running alphaVal %0.2f\n',alpha);
                                tic
%                                 [actAll,densityPlot,densityPlotAct,clusMu,muAvg,nTrlsUpd,gA,gW,gA_act,gW_act,gA_actNorm,gW_actNorm,muAll] = gauss_2d_batch_sim(nClus,locRange,box,warpType,epsMuOrig,sigmaGauss,nTrials,batchSize,nIter,warpBox,alpha,trials,stochasticType,c,plotGrids,dat,weightEpsSSE);
                                fprintf('Running %s, nClus=%d, epsMu=%d, sigmaGauss=%d, batchSize=%d\n',dat,nClus,epsMuOrig10000,sigmaGauss100,batchSize)
                                [densityPlot,gA,gA_act,gA_actNorm,rSeed,muAll]  = gauss_2d_batch_sim(nClus,locRange,box,warpType,epsMuOrig,sigmaGauss,nTrials,batchSize,nIter,warpBox,alpha,trials,stochasticType,c,plotGrids,dat,weightEpsSSE);
%                                 [actAll, densityPlot,gA,gW,rSeed,muAll] = gauss_2d_batch_sim(nClus,locRange,box,warpType,epsMuOrig,sigmaGauss,nTrials,batchSize,nIter,warpBox,alpha,trials,stochasticType,c,plotGrids,dat,weightEpsSSE);
%                                 fname = [saveDir, sprintf('/gauss_batch_%dclus__%dsigma_%dktrls_eps%d_alpha%d_stype%d_cVal%d_sseW%d_batchSiz%d_%diters',nClus,sigmaGauss100,round(nTrials/1000),epsMuOrig10000,alpha10,stochasticType,c1m,weightEpsSSE,round(batchSize),nIter)];
                                fname = [saveDir, sprintf('/gauss_batch_%dclus_%dsigma_%dktrls_eps%d_batchSiz%d_%diters_%s',nClus,sigmaGauss100,round(nTrials/1000),epsMuOrig10000,round(batchSize),nIter,dat)];                                
                                timeTaken=toc;
                                if saveDat
                                    if warpBox
                                        fname = [fname '_warpBox'];
                                    end
                                    
                                    if 1
                                        fname = [fname '_CE']; %crossentropy - testing out
                                    end
                                    
                                    cTime=datestr(now,'HHMMSS'); fname = sprintf([fname '_%s'],cTime);
%                                     save(fname,'densityPlot','densityPlotAct','clusMu','gA','gW','gA_act','gW_act','gA_actNorm','gW_actNorm','muAvg','nIter','timeTaken');
                                    save(fname,'densityPlot','clusMu','gA','gW','nIter','rSeed','timeTaken');
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
elseif neigh %testing neigh/eps
    
%     for iEps = 1:length(epsMuVals)
%         epsMuOrig=epsMuVals(iEps);
%         epsMuOrig1000=epsMuOrig*1000;
%         for iClus2run = 1:length(clus2run)
%             nClus = clus2run(iClus2run);
%             for iBeta=1:length(betaVals)
%                 beta = betaVals(iBeta);
%                 beta100 = beta*100; %for save filename
%                 tic
%                 [densityPlot,clusMu,muAvg,gA_g,gA_o,gA_wav,gA_rad,gW_g,gW_o,gW_wav,gW_rad]  = covering_map_sim_neigh(nClus,locRange,box,warpType,epsMuOrig,nTrials,nIter,warpBox,trials,beta,dat);
%                 fname = [saveDir, sprintf('/covering_map_dat_gauss_%dclus_%dtrls_eps%d_neigh_beta%d_%diters',nClus,nTrials,epsMuOrig1000,beta100,nIter)];
%                 timeTaken=toc;
%                 if saveDat
%                     if warpBox
%                         fname = [fname '_warpBox'];
%                     end
%                     cTime=datestr(now,'HHMMSS'); fname = sprintf([fname '_%s'],cTime);
% %                     save(fname,'densityPlot','clusMu','gA_g','gA_o','gA_wav','gA_rad','gW_g','gW_o','gW_wav','gW_rad','muAvg','nIter','timeTaken');
%                 end
%             end
%         end
%     end
end

toc
