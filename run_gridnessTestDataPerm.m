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
nTrialsTest = 100000; %?
dat = 'square';
% dat = 'circ';
% dat = 'trapzKrupic';

%for loading
nTrials = 1000000;
% clus2run = [3:30]; 
% clus2run=[3:10 12:2:26 11:2:25, 27:30]; % 27:30
%split
clus2run=[4:2:10, 12:2:26]; 
clus2run=[3:2:9, 11:2:25]; 

%batchSiz200 - rerunning after fail - love06
clus2run=[14:2:26]; 
clus2run=[13:2:25]; 

%love01 - batchSize=400 - didn't run perms last time; sq, circ
clus2run=[12, 4 16, 8, 24]; 
clus2run=[6 14, 10, 18, 20]; 
clus2run=[11, 5, 15, 9, 19]; 
clus2run=[3, 13, 7, 17, 23]; 
clus2run=[22, 26, 25, 21];

% next, batSiz=100



% % trapz jointTrls=0
% clus2run=[8:2:10, 12:2:26, 4, 6]; 
% clus2run=[7:2:9, 11:2:25,  3, 5]; 

% clus2run = 18;

nIter=200;
epsMuVals=.025;
nTrials=1000000;
% batchSizeVals = [400, 100]; % 125?
batchSizeVals=400;
% batchSizeVals=100;
% batchSizeVals=200;

annEps=0;

jointTrls=1; %for test trials


%running sq, circ, cLus sets x2, batchsize x2 - 8 matlabs - love01 - again

doPerm=1;


% run perm tests on how many iters? takes a bit of time (a couple mins) per
% iter, so with 200 iters plus many conditions, maybe too much (if all the
% perm data are about the same, then just take max, or 95th percentile as
% the threshold value)
nIters2run = 200; 
% nIters2run = 3; 

nPerm = 500;

saveDat=1;

for iClus2run = 1:length(clus2run) 
    nClus = clus2run(iClus2run);
    for iEps = 1%:length(epsMuVals) 
        epsMuOrig=epsMuVals(iEps);
        epsMuOrig1000=epsMuOrig*1000;
        for iBvals = 1%:length(batchSizeVals)
            batchSize = batchSizeVals(iBvals);
            
            if doPerm
                fprintf('Computing test trials gridness and running perm test on %s, nClus=%d, epsMu=%d, batchSize=%d\n',dat,nClus,epsMuOrig1000,batchSize)
            else
                fprintf('Computing test trials gridness on %s, nClus=%d, epsMu=%d, batchSize=%d\n',dat,nClus,epsMuOrig1000,batchSize)
            end
            
            %load
            fname = sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wAct_jointTrls_stepSiz',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat);
            if annEps %epsMu is different here
                fname = sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps*_batchSiz%d_%diters_%s_wAct_annEps_jointTrls_stepSiz',nClus,round(nTrials/1000),batchSize,nIter,dat);
            end
            %finish with directory and * for date/time
            fname = [saveDir, fname '*']; %finish with directory and * for date/time
            
            f = dir(fname); filesToLoad = cell(1,length(f));
            if isempty(f) %if no file, don't load/save - but print a warning
                warning('No file for: %s\n',fname);
            elseif ~isempty(f)
                
                %have to fix overlapping file names
                
                for iF = 1%:length(f)
%                     filesToLoad{iF} = f(iF).name;
                    load(f(iF).name);
                end
                
                %run
                tic
%                 [permPrc_gA, permPrc_gW,densityPlotAct,densityPlotActNorm,gA_act,gA_actNorm,gW_act,gW_actNorm, rSeedTest] = gridnessTestData_Perm(densityPlot,dat,locRange,nClus,nTrialsTest,nPerm,nIters2run);
                [permPrc_gA_act, permPrc_gW_act,permPrc_gA_actNorm, permPrc_gW_actNorm,gA_act,gA_actNorm,gW_act,gW_actNorm, gA_actNormPerm, gW_actNormPerm, densityPlotAct,densityPlotActNorm] = gridnessTestData_Perm(densityPlot,dat,locRange,nClus,nTrialsTest,nPerm,nIters2run,doPerm);
                timeTaken=toc;

                %save
                cTime=datestr(now,'HHMMSS'); 
                if doPerm
                    fname = [fname(1:end-1), sprintf('_perm_%dpermsOn%diters',nPerm,nIters2run)];
                else
                    fname = [fname(1:end-1), '_noPerm'];
                end
                if ~jointTrls
                    fname = [fname '_noJointTrlsTest'];
                end
                
                fname = [fname sprintf('_%s',cTime)];
                

                if saveDat
                    save(fname,'permPrc_gA_act','permPrc_gA_act','permPrc_gA_actNorm','permPrc_gA_actNorm','gA_act','gA_actNorm','gW_act','gW_actNorm','timeTaken');
                end
            end
        end
    end
end

%%
% 
% 
% figure; hist([gA_actNorm(22,:,1,2); gA_actNorm(22,:,1,3)]',25)
% % figure; hist(gA_actNorm(22,:,1,3),25)
% figure; hist(gA_actNorm(22,:,1,2)-gA_actNorm(22,:,1,3),25)
% 
% [h p  c s] = ttest(gA_actNorm(22,:,1,2)-gA_actNorm(22,:,1,3))
