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
epsMuVals=.025;
nTrials=1000000;
% batchSizeVals=400;
% batchSizeVals=200;

nIter=200;
% nIter=1000;

annEps=1;
if ~annEps
   epsMuVals=.025;
else
    epsMuVals=.25;
end
jointTrls=1;

% trapz sims here
nTrials2 = nTrials/2; %/4?
nBatches = 2500;
nBatches = 5000;%new
batchSizeVals = nTrials./nBatches; %match original batchsize values (nTrials being original nTrials)
nBvals = length(batchSizeVals);
nIter2run = nIter;

%new
epsMuTrapz = 0.0025;
% epsMuTrapz = 0.005; %new 
% epsMuTrapz = 0.003; %new 

%love06
% clus2run = [15, 23,  20, 24, 18, 16, 22, 27, 26, 19]; %half - ;14
% clus2run = [30, 29,  25, 10, 28, 17, 11, 21, 12]; %half - ; 13

% clus2run = [15, 23,  20, 24, 18, 16, 22, 27, 26, 19, 30, 29,  25, 10, 28, 17, 11, 21, 12]; % full - 13, 14 - % 

%10:20 only
% clus2run = [15, 13, 18, 16, 19, 14];
% clus2run = [10, 17, 11, 12, 20];


% rerunning batch400, 0.25 to 0.005, now with correct reduction over time
%rerunning 10:20, but also 21:30
% batchSiz=400, eps250-50 - OK

clus2run = [15, 13, 18, 16, 19, 14,21:2:30];
clus2run = [10, 17, 11, 12, 20, 22:2:30];

% running batchSiz200 0.25 to 00025
% now 0.25 to 0.005
clus2run = [21:2:30 ];
clus2run = 22:2:30;

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
            %finish with directory and * for date/time
            if ~annEps
                fname = [saveDir, fname '*']; %finish with directory and * for date/time
            else
                fname = [saveDir, fname '*annEps*']; %new position of annEps - works now - above no need
            end
            
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
                [densityPlot,~,gA,gW,~,~,rSeed] = covering_map_batch_sim_clusPosIn(clusPos,nClus,locRange,epsMuOrig,epsMuTrapz,nTrials,nTrials2,batchSize,nIter2run,dat2,annEps,jointTrls);
                
                fname = [saveDir, sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_epsMuTrapz_%d',nClus,round(nTrials2/1000),epsMuOrig1000,round(batchSize),nIter,dat2,epsMuTrapz*10000)];
                
                timeTaken=toc;
                if saveDat
                    if jointTrls
                        fname = [fname '_jointTrls_stepSiz'];
                    end
                    if annEps
                        fname = [fname '_annEps'];
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
