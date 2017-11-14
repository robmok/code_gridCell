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
clus2run = 40;
nTrials = 30000; %how many locations in the box / trials - 2.5k ; 5k if reset

colgrey = [.5, .5, .5];

%box
nSteps = 50; %to define spacing beween each loc in box
locRange = 1;%1; from -locRange to locRange
stepSize=diff(linspace(-1,1,nSteps)); stepSize=stepSize(1); %smallest diff between locs

% parameters
epsMuOrig=.075;% %learning rate / starting learning rate %.075
epsMuOrig=.1;%
% deltaEpsMu = .96;% %change in learning rate over time (slow down with 'learning')
% deltaEpsMu = .99; % slower decrease in learning rate for expanding (if no
% reset)

% for saving simulations - multiple by values to save the files with params
epsMuOrig1000=epsMuOrig*1000;
% deltaEpsMu100 = deltaEpsMu*100;

%reset learning rate?
% resetEps=0; %0 - no, 1 - once halfway, 2 - twice (quarter way and half way)  % ATM resetting to halfway
% resetMag=2; %1 - back to orig, 2 - 75% of the orig (if resetEps=2, 75%
% the orig then 50% the orig) %- not using this yet

%define box / environement - random points in a box
box = 'square'; %square, rect, trapz, trapzSq (trapz and a square box attached)

% change box shape during learning
%rectangle
warpBox = 0; %1 or 0
warpType = 'sq2rect';

%mometum-like adaptive learning rate - define alpha (higher = weight
%previous update (direction and magnitude) more; 0 = don't weight previous at all)
% alphaVals = [0, .1, .2, .3,.4, .5, .6, .7, .8, .9];
alphaVals = .3;
% alphaVals = [ 0, .1];
% alphaVals = [.2, .3];
% alphaVals = [.4, .5];
% alphaVals = [.6, .7];
% alphaVals = [.8, .9];

%%% things to understand

%stochastic update parameter - %note made 999 hack to test a constant c param
%(not changing over trials)
% cVals = [nTrials/(nTrials*100), nTrials/(nTrials*100*2), nTrials/(nTrials*100*4), 999]; % larger c = less stochastic over trials (becomes det quite early on); smaller c = more stochastic over trials (still a bit stochastic by the end)
% cVals = [nTrials/(nTrials*100), nTrials/(nTrials*100*2)]; 
% cVals = [nTrials/(nTrials*100*4), 999]; 
% cVals = [999]; 
% cVals = [nTrials/(nTrials*100)];

% % Create / load in saved test data

% trials = [randsample(linspace(-locRange,locRange,101),nTrials,'true'); randsample(linspace(-locRange,locRange,101),nTrials,'true')]'; % random points in a box
% save([saveDir '/randTrialsBox_10k'],'trials');

%%
saveDat=1; %save simulations - cluster centres and tsse

load([saveDir '/randTrialsBox_10k']); %load in same data with same trial sequence so same for each sim

nIter=10; %how many iterations (starting points)

% tic
% for iStochastic = 1:length(cVals) %
%     c = cVals(iStochastic);
%     c10k=c*10000; % for saving file name
%     fprintf('Running cVal %.3f\n',c);
tic
for iAlpha = 1:length(alphaVals) %
    alpha = alphaVals(iAlpha);
    alpha10 = alpha*10; %for saving simulations
    %     fprintf('Running alphaVal %0.2f\n',alpha);
    for iClus2run = 1:length(clus2run) %nClus conditions to run
        nClus = clus2run(iClus2run);
        [muEnd, muAll, tsseTrls,sseTrl,epsMuAll] = covering_map_sim(nClus,locRange,box,warpType,epsMuOrig,nTrials,nIter,warpBox,alpha,trials); %,c);
        
        if saveDat,
            % save simulations - might not need all tsseTrls if not viewing
            % them or epsMuAll?
            fname = [saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_alpha%d_%diters',nClus,nTrials,epsMuOrig1000,alpha10,nIter)];
%             fname = [saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_alpha%d_cVal%d_%diters',nClus,nTrials,epsMuOrig1000,alpha10,c10k,nIter)];
            if warpBox,
                fname = [fname '_warpBox'];
            end
            save(fname,'muEnd','muAll','sseTrl','tsseTrls','nIter');
%             save(fname,'muEnd','muAll','sseTrl','tsseTrls','nIter','indSSE1','indSSE2');
        end
    end
end
% toc
% end
toc