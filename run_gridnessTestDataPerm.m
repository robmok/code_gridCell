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
dat = 'circ';

%for loading
nTrials = 1000000;
jointTrls=1;
% clus2run = [3:30]; 
%clus2run=[3:10 12:2:30 11:2:29];
clus2run=[16:2:26 11:2:25]; %up to 26 only
% clus2run = 14;
nIter=200;
epsMuVals=.025;
nTrials=1000000;
% batchSizeVals = [400, 100]; % 125?
batchSizeVals=400;
% batchSizeVals=100;
annEps=0;

% run perm tests on how many iters? takes a bit of time (a couple mins) per
% iter, so with 200 iters plus many conditions, maybe too much (if all the
% perm data are about the same, then just take max, or 95th percentile as
% the threshold value)
nIters2run = 200; 

nPerm = 500;

for iClus2run = 1:length(clus2run) 
    nClus = clus2run(iClus2run);
    for iEps = 1%:length(epsMuVals) 
        epsMuOrig=epsMuVals(iEps);
        epsMuOrig1000=epsMuOrig*1000;
        for iBvals = 1%:length(batchSizeVals)
            batchSize = batchSizeVals(iBvals);
            
            fprintf('Running perm test on nClus=%d, epsMu=%d, batchSize=%d\n',nClus,epsMuOrig1000,batchSize)

            %load
            fname = sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wAct_jointTrls_stepSiz',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat);
            if annEps %epsMu is different here
                fname = sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps*_batchSiz%d_%diters_%s_wAct_jointTrls_stepSiz_annEps',nClus,round(nTrials/1000),batchSize,nIter,dat);
            end
            %finish with directory and * for date/time
            fname = [saveDir, fname '*']; %finish with directory and * for date/time
            
            f = dir(fname); filesToLoad = cell(1,length(f));
            if isempty(f) %if no file, don't load/save - but print a warning
                warning('No file for: %s\n',fname);
            elseif ~isempty(f)
                for iF = 1%:length(f)
%                     filesToLoad{iF} = f(iF).name;
                    load(f(iF).name);
                end
                
                %run
                tic
%                 [permPrc_gA, permPrc_gW,densityPlotAct,densityPlotActNorm,gA_act,gA_actNorm,gW_act,gW_actNorm, rSeedTest] = gridnessTestData_Perm(densityPlot,dat,locRange,nClus,nTrialsTest,nPerm,nIters2run);
                [permPrc_gA,permPrc_gW,gA_act,gA_actNorm,gW_act,gW_actNorm, gA_actNormPerm, gW_actNormPerm, densityPlotAct,densityPlotActNorm] = gridnessTestData_Perm(densityPlot,dat,locRange,nClus,nTrialsTest,nPerm,nIters2run);
                timeTaken=toc;

                %save
                cTime=datestr(now,'HHMMSS'); 
                fname = [fname(1:end-1), sprintf('_perm_%dpermsOn%diters_%s',nPerm,nIters2run,cTime)];
                save(fname,'permPrc_gA','permPrc_gW','gA_act','gA_actNorm','gW_act','gW_actNorm','timeTaken')
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
