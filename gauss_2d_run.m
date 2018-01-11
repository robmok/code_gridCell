clear all;
% close all;

% wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
% wd='/home/robmok/Documents/Grid_cell_model'; %on love01

cd(wd);
codeDir = [wd '/code_gridCell'];
saveDir = [wd '/data_gridCell'];
addpath(codeDir); addpath(saveDir);
addpath(genpath([wd '/gridSCORE_packed']));

%run multuple cluster numbers
clus2run = 20; %20, 30
nTrials = 40000; %how many locations in the box / trials - 2.5k ; 5k if reset

%box
nSteps = 50; %to define spacing beween each loc in box
locRange = [0, nSteps-1]; %[-1, 1]; % from locRange(1) to locRange(2)
stepSize=diff(linspace(locRange(1),locRange(2),nSteps)); stepSize=stepSize(1); %smallest diff between locs

% parameters
% epsMuVals=[.05 .075 .1];% %learning rate / starting learning rate 
epsMuVals = .009; %only 0.009/0.01 looks gd now..

%define box / environement - random points in a box
box = 'square'; %square, rect, trapz, trapzSq (trapz and a square box attached)

% change box shape during learning
%rectangle
warpBox = 0; %1 or 0
warpType = 'sq2rect';

%mometum-like adaptive learning rate - define alpha (higher = weight
%previous update (direction and magnitude) more; 0 = don't weight previous at all)
alphaVals = [0, .2, .5, .8];
% alphaVals = [0, .2];
% alphaVals = [.5, .8];
alphaVals=0;

sTypes = 0;%:1;% :3; %0, 1 ,2, 3
cValsOrig = [2/nTrials, 5/nTrials, 10/nTrials, 20/nTrials]; %removed .1/nTrials and .25/nTrials,  too stochastic. also 3/ntrials, .5/nTrials
% cValsOrig = [2/(nTrials/2), 5/nTrials, 10/(nTrials/2), 20/(nTrials/2)];% if 80k trials...
% cValsOrig = 20/nTrials;

%neighbour-weighted update
neigh = 0; %if neigh = 0, no stoch, no alpha
betaVals  = [.3 .5 .7]; %softmax param - higher = less neighbour update; from .25 going toward middle; .3 start OK
if neigh
   sTypes=0;
   alphaVals=0;
   epsMuVals=[.05 .075 .1];
end
% % Create / load in saved test data
% trials = [randsample(linspace(-locRange,locRange,101),nTrials,'true'); randsample(linspace(-locRange,locRange,101),nTrials,'true')]'; % random points in a box
% trials = [randsample(linspace(locRange(1),locRange(2),nSteps),nTrials,'true'); randsample(linspace(locRange(1),locRange(2),nSteps),nTrials,'true')]'; % random points in a box
% save([saveDir '/randTrialsBox_30k'],'trials');

% %new - tile the whole space rather than sample
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
% save([saveDir '/randTrialsBox_80k'],'trials');

%%
saveDat=1; %save simulations

nIter=50; %how many iterations (starting points)

if nTrials==40000
    load([saveDir '/randTrialsBox_40k']); %load in same data with same trial sequence so same for each sim
elseif nTrials==80000
    load([saveDir '/randTrialsBox_80k']);
end

tic
if ~neigh %separating neigh and stoch/momentum params
    
    %testing stoch/momentum/eps
    for iEps = 1:length(epsMuVals)
        epsMuOrig=epsMuVals(iEps);
        epsMuOrig1000=epsMuOrig*1000; %for saving - changed from 1000 to 100000 for slower l rates - changed back now
        for iClus2run = 1:length(clus2run) %nClus conditions to run
            nClus = clus2run(iClus2run);
            for iStype = 1:length(sTypes)
                stochasticType = sTypes(iStype);
                if ~stochasticType
                    cVals = 0;
                else
                    cVals = cValsOrig;
                end
                for iStochastic = 1:length(cVals) %
                    c = cVals(iStochastic);
                    c1m=round(c.*1000000); % for saving file name
                    fprintf('Running cVal %d\n',c1m)
                    for iAlpha  = 1:length(alphaVals) %
                        alpha = alphaVals(iAlpha);
                        alpha10 = alpha*10; %for saving simulations
                        fprintf('Running alphaVal %0.2f\n',alpha);
                        tic
%                         [densityPlot,clusMu,muAvg,nTrlsUpd,gA_g,gA_o,gA_wav,gA_rad,gW_g,gW_o,gW_wav,gW_rad,cParams] = gauss_2d_sim(nClus,locRange,box,warpType,epsMuOrig,nTrials,nIter,warpBox,alpha,trials,stochasticType,c);
%                         [muAll,actAll,densityPlot,densityPlotAct,clusMu,muAvg,nTrlsUpd,gA_g,gA_o,gA_wav,gA_rad,gW_g,gW_o,gW_wav,gW_rad,cParams] = gauss_2d_sim(nClus,locRange,box,warpType,epsMuOrig,nTrials,nIter,warpBox,alpha,trials,stochasticType,c);
                        [muAll,actAll,densityPlot,densityPlotAct,clusMu,muAvg,nTrlsUpd,gA,gW,cParams] = gauss_2d_sim(nClus,locRange,box,warpType,epsMuOrig,nTrials,nIter,warpBox,alpha,trials,stochasticType,c);

                        fname = [saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_alpha%d_stype%d_cVal%d_%diters',nClus,nTrials,epsMuOrig1000,alpha10,stochasticType,c1m,nIter)];
                        timeTaken=toc;
                        if saveDat
                            if warpBox
                                fname = [fname '_warpBox'];
                            end
                            cTime=datestr(now,'HHMMSS'); fname = sprintf([fname '_%s'],cTime);
                            save(fname,'densityPlot','clusMu','gA_g','gA_o','gA_wav','gA_rad','gW_g','gW_o','gW_wav','gW_rad','muAvg','nIter','cParams','timeTaken');
%                             save(fname,'densityPlot','clusMu','gA','gW','muAvg','nIter','cParams','timeTaken');
                        end
                    end
                end
            end
        end
    end
    
elseif neigh %testing neigh/eps
    
    for iEps = 1:length(epsMuVals)
        epsMuOrig=epsMuVals(iEps);
        epsMuOrig1000=epsMuOrig*1000;
        for iClus2run = 1:length(clus2run)
            nClus = clus2run(iClus2run);
            for iBeta=1:length(betaVals)
                beta = betaVals(iBeta);
                beta100 = beta*100; %for save filename
                tic
                [densityPlot,clusMu,muAvg,gA_g,gA_o,gA_wav,gA_rad,gW_g,gW_o,gW_wav,gW_rad]  = covering_map_sim_neigh(nClus,locRange,box,warpType,epsMuOrig,nTrials,nIter,warpBox,trials,beta);
                fname = [saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_neigh_beta%d_%diters',nClus,nTrials,epsMuOrig1000,beta100,nIter)];
                timeTaken=toc;
                if saveDat
                    if warpBox
                        fname = [fname '_warpBox'];
                    end
                    cTime=datestr(now,'HHMMSS'); fname = sprintf([fname '_%s'],cTime);
%                     save(fname,'densityPlot','clusMu','gA_g','gA_o','gA_wav','gA_rad','gW_g','gW_o','gW_wav','gW_rad','muAvg','nIter','timeTaken');
                end
            end
        end
    end
end

toc
