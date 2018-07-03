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
nTrialsTest = 100000; % orig nTrials/10
dat = 'circ';
dat = 'square';

% dat = 'trapzKrupic';

dat = 'trapzKfrmSq1'; % load covering map on sq, then run it on trapz; then assess gridness in trapz
% dat = 'trapzKfrmSq2'; % load covering map on sq, then assess gridness in trapz

saveDat=1;

%for loading
% nTrials = 1000000;
% clus2run = [3:30];
% clus2run=[3:10 12:2:26 11:2:25, 27:30]; % 27:30


% 1k iters, sq/circ annEps, no perm - % rerun
% clus2run = [18,  8, 16,  4, 22, 27, 26,  9, 13, 29,  5, 25, 7, 10, 19, 28, 17, 11, 21, 12,  30, 3]; %
% clus2run = [14, 15, 23,  6, 20, 24, 18,  8, 16,  4, 22, 27, 26,  9, 13, 30, 29,  5, 25, 10,  7, 19, 28, 17, 11, 21, 12,  3]; %all



%%%%%%%
% new fixed perm - 200iters
% annEps
%running on love01, 4 matlabs- circ
clus2run = [3, 15, 23,  26, 20, 18, 4];
% clus2run = [16, 8, 22, 6, 9, 30, 24];
% clus2run = [5,  13, 7, 10,19, 29, 25];
% clus2run = [11, 21, 12, 14, 27, 28, 17]; 


% sq love06 - not started yet - probably just run 2
% clus2run = [14, 15, 23,  30, 6, 20, 16, 22, 26];
% clus2run = [5,  13, 7,  10, 19, 29, 11, 21, 12];
% clus2run = [9, 18, 24,4, 25, 17, 27, 28, 3, 8];




% % next: sq2trapz 1kiters
% clus2run = [14, 15, 23,  6, 20, 24, 18,  8, 16,  4, 22, 27, 26,  9, 13, 30, 29,  5, 25, 10,  7, 19, 28, 17, 11, 21, 12,  3]; %all
% clus2run = [14, 15, 23,  6, 20, 24, 18,  8, 16,  4, 22, 27, 26, 9]; %half
% clus2run = [13, 30, 29,  5, 25, 10,  7, 19, 28, 17, 11, 21, 12,  3]; %half
clus2run = [13, 30, 29,  5, 25, 10,  7, ]; %quarter
clus2run = [19, 28, 17, 11, 21, 12,  3]; %quarter

% clus2run = 11;

%sq 200 iters no annEps no perm
% clus2run = 27:30;

nIter=200;
epsMuVals=.025;
nTrials=1000000;
% batchSizeVals = [400, 100]; % 125?
batchSizeVals=400;
% batchSizeVals=200;
% batchSizeVals=100;

annEps=0;
if annEps
    epsMuVals=1.5; %below multiplies this by 1k
end
jointTrls=1; %for test trials

nIter=1000;

% if trapKfrmSq1
if strcmp(dat,'trapzKfrmSq1')
    nTrials=1000000/2;
    epsMuTrapz10 = 25; %this is 10% of orig learning rate - using this
    % epsMuTrapz10 = 50;
    % epsMuTrapz10 = 15;  %
end


%doPerm or not
if ~strcmp(dat(1:4),'trap')
    doPerm=1;
else
    doPerm=0;
end

%tmp - for getting sq density plots
% doPerm=0;

% run perm tests on how many iters? takes a bit of time (a couple mins) per
% iter, so with 200 iters plus many conditions, maybe too much (if all the
% perm data are about the same, then just take max, or 95th percentile as
% the threshold value)
nIters2run = nIter; %200
nPerm = 500;

% %testing
% nIters2run = 3;
% nPerm=250;

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
            
            %note, annEps slightly different file name
            if ~annEps
                if ~(length(dat)==12)
                    fname = sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_jointTrls_stepSiz',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat);
                else %trapzKfrmSq1 or 2 - i think only 1 now
                    if strcmp(dat(12),'1')
                        fname = sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_epsMuTrapz10_%d_jointTrls_stepSiz',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat,epsMuTrapz10);
                    elseif strcmp(dat(12),'2')
                        fname = sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_jointTrls_stepSiz',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,'square');
                    end
                end
            else
                if ~(length(dat)==12)
                    fname = sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_annEps_jointTrls_stepSiz',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat);
                else %trapzKfrmSq1
                    fname = sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_epsMuTrapz10_%d_jointTrls_stepSiz',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat,epsMuTrapz10);
                    
                end
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
                % if doing trapzKfrmSq, record and run and save square gridness; no perm
                if length(dat)==12
                    if strcmp(dat(12),'2') %only doing withut learning atm
                        [~,~,~,~,gA_act_sq,gA_actNorm_sq,gW_act_sq,gW_actNorm_sq] = gridnessTestData_Perm(densityPlot,'square',locRange,nClus,nTrialsTest,nPerm,nIters2run,0);
                        gA_sq = gA;
                    end
                end
                
                %run
                tic
                %                 [permPrc_gA, permPrc_gW,densityPlotAct,densityPlotActNorm,gA_act,gA_actNorm,gW_act,gW_actNorm, rSeedTest] = gridnessTestData_Perm(densityPlot,dat,locRange,nClus,nTrialsTest,nPerm,nIters2run);
                [permPrc_gA_act, permPrc_gW_act,permPrc_gA_actNorm, permPrc_gW_actNorm,gA_act,gA_actNorm,gW_act,gW_actNorm, gA_actNormPerm, gW_actNormPerm, densityPlotAct,densityPlotActNorm] = gridnessTestData_Perm(densityPlot,dat,locRange,nClus,nTrialsTest,nPerm,nIters2run,doPerm);
                timeTaken=toc;
                
                %save
                cTime=datestr(now,'HHMMSS');
                if doPerm
                    %                     fname = [fname(1:end-1), sprintf('_perm_%dpermsOn%diters',nPerm,nIters2run)];
                    fname = [fname(1:end-1), sprintf('_actNorm_perm_%dpermsOn%diters',nPerm,nIters2run)]; %now correct actNorm; add 'trlsTest'?
                else
                    fname = [fname(1:end-1), '_trlsTest_noPerm'];
                end
                if length(dat)==12
                    fname = [fname sprintf('_%s',dat)];
                end
                if ~jointTrls
                    fname = [fname '_noJointTrlsTest'];
                end
                
                fname = [fname sprintf('_%s',cTime)];
                
                
                if saveDat
                    if length(dat)==12 && strcmp(dat(12),'2')
                        save(fname,'gA_sq','gA_act_sq','gA_actNorm_sq','gW_act_sq','gW_actNorm_sq','gA_act','gA_actNorm','gW_act','gW_actNorm','densityPlotAct','densityPlotActNorm','timeTaken');
                    elseif length(dat)==12 && strcmp(dat(12),'1')
                        save(fname,'gA_act','gA_actNorm','gW_act','gW_actNorm','densityPlotAct','densityPlotActNorm','timeTaken');
                    else
                        %                     save(fname,'permPrc_gA_act','permPrc_gA_act','permPrc_gA_actNorm','permPrc_gA_actNorm','gA_act','gA_actNorm','gW_act','gW_actNorm','timeTaken');
                        %also save some density plots for figures
                        save(fname,'permPrc_gA_act','permPrc_gW_act','permPrc_gA_actNorm','permPrc_gW_actNorm','gA_act','gA_actNorm','gW_act','gW_actNorm','densityPlotAct','densityPlotActNorm','timeTaken');
                    end
                end
                clear densityPlotAct densityPlotActNorm gA_act gA_actNorm gW_act gW_actNorm permPrc_gA_act permPrc_gW_act
            end
        end
    end
end