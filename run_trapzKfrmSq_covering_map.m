% load in sq map already run, then run covering map on trapz (smaller)
% shape

clear all;

% wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
% wd='/home/robmok/Documents/Grid_cell_model'; %on love01

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

%new running
nIter=1000;


% for current
nTrials2 = nTrials/2; %/4?
nBatches = 2500;
batchSizeVals = nTrials./nBatches;
nBvals = length(batchSizeVals);
nIter2run = nIter;

%new
epsMuTrapz = 0.0025;
% epsMuTrapz = 0.005;
% epsMuTrapz = 0.0015; %to run - slower

% clus2run = 3:26; 

%love06
clus2run=[8:4:20, 23, 10:4:22, 27, 4];
clus2run=[26, 24, 9:4:21,29, 3, 6]; 
clus2run=[7:4:19, 25, 5, 28, 30]; 



%testing
% clus2run = 13;
% nIter2run = 1;
%%

for iClus2run = 1:length(clus2run)
    nClus = clus2run(iClus2run);
    for iEps = 1:length(epsMuVals)
        epsMuOrig=epsMuVals(iEps);
        epsMuOrig1000=epsMuOrig*1000;
        for iBvals = 1:length(batchSizeVals)
            batchSize = batchSizeVals(iBvals);
            fprintf('Running %s, nClus=%d, epsMu=%d, batchSize=%d\n',dat2,nClus,epsMuOrig1000, batchSize)
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
                tic
                [densityPlot,~,gA,gW,~,~,rSeed] = covering_map_batch_sim_clusPosIn(clusPos,nClus,locRange,epsMuTrapz,nTrials2,batchSize,nIter2run,dat2,annEps,jointTrls);
                
                fname = [saveDir, sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_epsMuTrapz10_%d',nClus,round(nTrials2/1000),epsMuOrig1000,round(batchSize),nIter,dat2,epsMuTrapz*10000)];
                
                timeTaken=toc;
                if saveDat
                    if jointTrls
                        fname = [fname '_jointTrls_stepSiz'];
                    end
                    cTime=datestr(now,'HHMMSS'); fname = sprintf([fname '_%s'],cTime);
%                     save(fname,'densityPlot','densityPlotActNorm','gA','gW','gA_actNorm','gW_actNorm','rSeed','timeTaken');
                    save(fname,'densityPlot','gA','gW','rSeed','timeTaken');
                    clear densityPlot gA gW
                end
            end
        end
    end
end
