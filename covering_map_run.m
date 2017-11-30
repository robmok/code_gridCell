clear all;
% close all;

% wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
% wd='/home/robmok/Documents/Grid_cell_model'; %on love01

cd(wd);
codeDir = [wd '/code_gridCell'];
saveDir = [wd '/data_gridCell'];
addpath(codeDir); addpath(saveDir);

% clus2run = [20, 40, 60, 80]; %run multuple cluster numbers
clus2run = 20; %20, 30
nTrials = 40000; %how many locations in the box / trials - 2.5k ; 5k if reset

%box
nSteps = 50; %to define spacing beween each loc in box
locRange = [0, nSteps-1]; %[-1, 1]; % from locRange(1) to locRange(2)
stepSize=diff(linspace(locRange(1),locRange(2),nSteps)); stepSize=stepSize(1); %smallest diff between locs

% parameters
epsMuOrig=.075;% %learning rate / starting learning rate %.075
% epsMuOrig=.1;
% epsMuOrig=.05;
% epsMuOrig=.125; %could do this as well
% epsMuOrig = .01;

% for saving simulations - multiple by values to save the files with params
epsMuOrig1000=epsMuOrig*1000;
% deltaEpsMu100 = deltaEpsMu*100;

%define box / environement - random points in a box
box = 'square'; %square, rect, trapz, trapzSq (trapz and a square box attached)

% change box shape during learning
%rectangle
warpBox = 0; %1 or 0
warpType = 'sq2rect';

%mometum-like adaptive learning rate - define alpha (higher = weight
%previous update (direction and magnitude) more; 0 = don't weight previous at all)
% alphaVals = [.2, .5, .8];
alphaVals = .2 ; %.2 .5 .8

sTypes = 1;%:3; %0, 1 ,2, 3
% 0. none
% 1. standard stochastic update - becomes more det over time; becomes
% basically deterministic at some point
% 2. stochastic update with a limit - becomes more det over time, but keeps
%some (very low) stochasticity at the end
% 3. constant stochasticity - keep very low

%  larger c = less stochastic over trials (becomes det quite early on); smaller c = more stochastic over trials (still a bit stochastic by the end)
cVals = [2/nTrials, 3/nTrials, 5/nTrials, 10/nTrials]; %half1
% cVals = [.1/nTrials, .25/nTrials, .5/nTrials, 20/nTrials]; %half2
cVals = 10/nTrials;

% % Create / load in saved test data
% trials = [randsample(linspace(-locRange,locRange,101),nTrials,'true'); randsample(linspace(-locRange,locRange,101),nTrials,'true')]'; % random points in a box
% trials = [randsample(linspace(locRange(1),locRange(2),nSteps),nTrials,'true'); randsample(linspace(locRange(1),locRange(2),nSteps),nTrials,'true')]'; % random points in a box
% save([saveDir '/randTrialsBox_30k'],'trials');

%new - tile the whole space rather than sample
% sq=linspace(locRange(1),locRange(2),nSteps);
% allPts=[];
% for i=1:length(sq)
%     for j=1:length(sq)
%         allPts = [allPts; [sq(i), sq(j)]];
%     end
% end
% trials=repmat(allPts,nTrials/length(allPts),1); %note, numel of allPts must be divisble by nTrials atm
% trials=trials(randperm(length(trials)),:);
% save([saveDir '/randTrialsBox_40k'],'trials');

%%
saveDat=0; %save simulations

load([saveDir '/randTrialsBox_40k']); %load in same data with same trial sequence so same for each sim

nIter=1; %how many iterations (starting points)

tic
for iClus2run = 1:length(clus2run) %nClus conditions to run
    nClus = clus2run(iClus2run);
    for iStype = 1:length(sTypes)
        stochasticType = sTypes(iStype);
        for iStochastic = 1:length(cVals) %
            c = cVals(iStochastic);
            c1m=round(c.*1000000); % for saving file name
            fprintf('Running cVal %d\n',c1m)
            if ~stochasticType, c=0; end
            for iAlpha  = 1:length(alphaVals) %
                alpha = alphaVals(iAlpha);
                alpha10 = alpha*10; %for saving simulations
                fprintf('Running alphaVal %0.2f\n',alpha);
                tic
                [muAll,cParams] = covering_map_sim(nClus,locRange,box,warpType,epsMuOrig,nTrials,nIter,warpBox,alpha,trials,stochasticType,c);
%                 muAll = covering_map_sim_neigh(nClus,locRange,box,warpType,epsMuOrig,nTrials,nIter,warpBox,alpha,trials,stochasticType,c);
                timeTaken=toc;
                if saveDat
                    %                 fname = [saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_alpha%d_%diters',nClus,nTrials,epsMuOrig1000,alpha10,nIter)];
                    fname = [saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_alpha%d_stype%d_cVal%d_%diters',nClus,nTrials,epsMuOrig1000,alpha10,stochasticType,c1m,nIter)];
                    if warpBox
                        fname = [fname '_warpBox'];
                    end
                    cTime=datestr(now,'HHMMSS'); fname = sprintf([fname '_%s'],cTime);
                    save(fname,'muAll','nIter','cParams','timeTaken');
                end
            end
        end
    end
end
toc

% figure; plot(cParams.closestChosen)
% propClosestC = nnz(cParams.closestChosen)/nTrials