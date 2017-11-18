clear all;
% close all;

wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
% wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';

cd(wd);
codeDir = [wd '/code_gridCell'];
saveDir = [wd '/data_gridCell'];
addpath(codeDir); addpath(saveDir);

% nClus   = 20;
% clus2run = [20, 40, 60, 80]; %run multuple cluster numbers
clus2run = [20 40];
nTrials = 30000; %how many locations in the box / trials - 2.5k ; 5k if reset

colgrey = [.5, .5, .5];

%box
nSteps = 500;%1001;%50; %to define spacing beween each loc in box
locRange = [0, nSteps-1]; %[-1, 1]; % from locRange(1) to locRange(2)
stepSize=diff(linspace(locRange(1),locRange(2),nSteps)); stepSize=stepSize(1); %smallest diff between locs

% parameters
epsMuOrig=.075;% %learning rate / starting learning rate %.075
epsMuOrig=.1;
% epsMuOrig=.05;

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
% alphaVals = [0, .1, .2, .3,.4, .5, .6, .7, .8, .9];
alphaVals = [.2, .5, .8];
% alphaVals = [ 0, .1];
% alphaVals = [.2, .3];
% alphaVals = [.4, .5];
% alphaVals = [.6, .7];
% alphaVals = [.8, .9];

stochasticType = 3; %0, 1 ,2, 3
% 0. none
% 1. standard stochastic update - becomes more det over time; becomes
% basically deterministic at some point
% 2. stochastic update with a limit - becomes more det over time, but keeps
%some (very low) stochasticity at the end
% 3. constant stochasticity - keep very low

%  larger c = less stochastic over trials (becomes det quite early on); smaller c = more stochastic over trials (still a bit stochastic by the end)
cVals = [2/nTrials, 3/nTrials, 5/nTrials, 10/nTrials];

% % Create / load in saved test data
% trials = [randsample(linspace(-locRange,locRange,101),nTrials,'true'); randsample(linspace(-locRange,locRange,101),nTrials,'true')]'; % random points in a box
% trials = [randsample(linspace(locRange(1),locRange(2),nSteps),nTrials,'true'); randsample(linspace(locRange(1),locRange(2),nSteps),nTrials,'true')]'; % random points in a box
% save([saveDir '/randTrialsBox_30k'],'trials');

%%
saveDat=1; %save simulations - cluster centres and tsse

load([saveDir '/randTrialsBox_30k']); %load in same data with same trial sequence so same for each sim

nIter=3; %how many iterations (starting points)

% tic
for iStochastic = 1:length(cVals) %
    c = cVals(iStochastic);
    c1m=round(c.*1000000); % for saving file name
    fprintf('Running cVal %d\n',c1m)
    if ~stochasticType, c=[]; end
    tic
    for iAlpha = 1:length(alphaVals) %
        alpha = alphaVals(iAlpha);
        alpha10 = alpha*10; %for saving simulations
            fprintf('Running alphaVal %0.2f\n',alpha);
        for iClus2run = 1:length(clus2run) %nClus conditions to run
            nClus = clus2run(iClus2run);
            %         [muEnd, muAll, tsseTrls,sseTrl,epsMuAll] = covering_map_sim(nClus,locRange,box,warpType,epsMuOrig,nTrials,nIter,warpBox,alpha,trials); %,c);
            [muEnd, muAll, tsseTrls,sseTrl,epsMuAll,cParams] = covering_map_sim(nClus,locRange,box,warpType,epsMuOrig,nTrials,nIter,warpBox,alpha,trials,stochasticType,c);
            
            if saveDat
                % save simulations - might not need all tsseTrls if not viewing
                % them or epsMuAll?
%                 fname = [saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_alpha%d_%diters',nClus,nTrials,epsMuOrig1000,alpha10,nIter)];
                fname = [saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_alpha%d_stype%d_cVal%d_%diters',nClus,nTrials,epsMuOrig1000,alpha10,stochasticType,c1m,nIter)];
                if warpBox
                    fname = [fname '_warpBox'];
                end
                save(fname,'muEnd','muAll','sseTrl','tsseTrls','nIter','cParams');
                %             save(fname,'muEnd','muAll','sseTrl','tsseTrls','nIter','indSSE1','indSSE2');
            end
        end
    end
    toc
end
toc