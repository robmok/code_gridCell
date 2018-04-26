clear all;

wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
% wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
cd(wd);

codeDir = [wd '/code_gridCell'];
saveDir = [wd '/data_gridCell'];
addpath(codeDir); addpath(saveDir);
addpath(genpath([codeDir '/gridSCORE_packed']));


locRange = [0, 49];


dat='catLearn';
% annEps=0;
boxSize=1;
nIter=20;

% joined trials
jointTrls=0;
clus2run = 2:26; 
epsMuVals=.025;
nTrials=50000;
batchSizeVals= 10;

nCats=2;
stoch=1;
cVals = [2,4, 10, 40];
cVals = 4;

catsInfo.nCats=2; %2 categories
% sigmaG = [5 0; 0 5];   % isotropic % sigmaG = [1 .5; .5 2]; R = chol(sigmaG);  % non-isotropic
sigmaG = [3 0; 0 3];
catsInfo.R=chol(sigmaG);

%load loop
for iClus = 1:length(clus2run)
    nClus = clus2run(iClus);
    for iEps = 1:length(epsMuVals)
        epsMuOrig=epsMuVals(iEps);
        epsMuOrig1000=epsMuOrig*1000;
        for iBvals = 1:length(batchSizeVals)
            for iC = 1:length(cVals)
                fprintf('Loading %s, nClus=%d, epsMu=%d, c=%d, batchSize=%d\n',dat,nClus,epsMuOrig1000,cVals(iC),batchSizeVals(iBvals))

                fname = [saveDir, sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wAct_%dcats_stoch%d_c%d_*',nClus,round(nTrials/1000),epsMuOrig1000,batchSizeVals(iBvals),nIter,dat,nCats,stoch,cVals(iC))];
                f=dir(fname);
                load(f.name);
                muAllClus{iClus}=muAll;

                
                rSeedAll{iClus} = rSeed;
                
            end
        end
    end
end
    
%%
nBatches = size(muAll,3)-1;
trls2Plt = [1, nBatches*.25, nBatches*.5, nBatches*.75, nBatches+1];
% trls2Plt = [nBatches+1];

iterI=1;

for iClus = 1:length(clus2run)
    colors  = distinguishable_colors(clus2run(iClus));    
    figure; hold on;
    for iPlot = 1:length(trls2Plt)
        subplot(1,5,iPlot);
        scatter(muAllClus{iClus}(:,1,trls2Plt(iPlot),iterI),muAllClus{iClus}(:,2,trls2Plt(iPlot),iterI),750,colors,'.');hold on;
        trials = createTrls(dat,nTrials,locRange,1,jointTrls,boxSize,catsInfo,rSeedAll{iClus}(iterI));
        scatter(trials(:,1),trials(:,2),'.');
        xlim(locRange+1);
        ylim(locRange+1);
    end
    
end

%recreate trials - to plot as well
% muInit
% rSeed