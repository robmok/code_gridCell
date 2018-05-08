% load in sq map already run, then run covering map on trapz (smaller)
% shape

clear all;

% wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
wd='/home/robmok/Documents/Grid_cell_model'; %on love01

cd(wd);

codeDir = [wd '/code_gridCell'];
saveDir = [wd '/data_gridCell'];
addpath(codeDir); addpath(saveDir);
addpath(genpath([codeDir '/gridSCORE_packed']));

locRange = [0 49];
dat1 = 'square';
dat2 = 'trapzKfrmSq1'; % run covering map on sq, then assess gridness in trapz

saveDat=1;

%for loading
nIter=200;
epsMuVals=.025;
nTrials=1000000;
% batchSizeVals = [400, 100]; % 125/200?
batchSizeVals=400;

annEps=0;
jointTrls=1;

% for current
nTrials2 = nTrials/2; %/4?
nBatches = 2500;
batchSizeVals = nTrials./nBatches;
nBvals = length(batchSizeVals);
nIter2run = nIter;

% clus2run = 3:26; 

% Split, even
clus2run=[8:4:24, 6]; 
% clus2run=[10:4:26,4];

%testing
% clus2run = 20;
% nIter2run = 1;
%%

for iClus2run = 1:length(clus2run)
    nClus = clus2run(iClus2run);
    for iEps = 1:length(epsMuVals)
        epsMuOrig=epsMuVals(iEps);
        epsMuOrig1000=epsMuOrig*1000;
        for iBvals = 1:length(batchSizeVals)
            batchSize = batchSizeVals(iBvals);
            
            %load in sq
            fname = sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_jointTrls_stepSiz',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat1);
            fname = [saveDir, fname '*']; %finish with directory and * for date/time
            
            f = dir(fname); filesToLoad = cell(1,length(f));
            if isempty(f) %if no file, don't load/save - but print a warning
                warning('No file for: %s\n',fname);
            elseif ~isempty(f)
                for iF = 1%:length(f)
                    %filesToLoad{iF} = f(iF).name;
                    load(f(iF).name);
                end
                
                %get cluster centres
                clusPos = nan(nClus,2,nIter);
                for iterI = 1:nIter
                    [clusPos(:,1,iterI), clusPos(:,2,iterI)] = find(densityPlot(:,:,end,iterI)>0); %find cluster positions - added >0 since counts nans as well
                end
                
                % run - allow starting clus positions (even
                % outside the trapz box)
                
                [densityPlot,densityPlotActNorm,gA,gW,gA_actNorm,gW_actNorm,rSeed,muAll,trials] = covering_map_batch_sim_clusPosIn(clusPos,nClus,locRange,epsMuOrig,nTrials2,batchSize,nIter2run,dat2,annEps,jointTrls);
                
                
                fname = [saveDir, sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm',nClus,round(nTrials2/1000),epsMuOrig1000,round(batchSize),nIter,dat2)];
                
                timeTaken=toc;
                if saveDat
                    if jointTrls
                        fname = [fname '_jointTrls_stepSiz'];
                    end
                    cTime=datestr(now,'HHMMSS'); fname = sprintf([fname '_%s'],cTime);
                    save(fname,'densityPlot','densityPlotActNorm','gA','gW','gA_actNorm','gW_actNorm','rSeed','timeTaken');
                end
            end
        end
    end
end
