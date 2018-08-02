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
% dat = 'square';
% dat = 'trapzKfrmSq1'; % load covering map on sq, then run it on trapz; then assess gridness in trapz

saveDat=1;

% nIter=200;
epsMuVals=.025;
nTrials=1000000;
% batchSizeVals=400;
batchSizeVals=200; %new

annEps=1;
if annEps %new
    epsMuVals=.25; %below multiplies this by 1k
end

nIter=200;
nIter=1000;

%%%%%%%
% new fixed perm - 200iters
% annEps 

% sq - not started
clus2run = [16, 22, 24, 14, 15, 27];
% clus2run = [25, 11, 21, 12, 30, 26]; 
% clus2run = [13, 10, 17, 18, 29];

%circ - love06 - not started
% clus2run = [19, 20, 26, 28, 23, 25];


%circ - love01 - perm1-5 - running
% clus2run = [16, 22, 24];
% clus2run = [14, 15, 11];
% clus2run = [27, 21, 12]; 
% clus2run = [13, 10, 17];
% clus2run = [18, 29, 30];

% sq - love01 - perm6 - running
% clus2run = [19, 20, 28, 23];

%%%%
% 1k iters, sq/circ annEps, batchSiz=200, no perm 

%sq
% clus2run = [16, 22, 24, 14, 15, 27, 19, 13, 10, 17, 18];
% clus2run = [11, 21, 12, 30, 25, 28, 20, 29, 23, 26]; 

%circ - love01
% clus2run = [16, 22, 24, 14, 15, 27];
% clus2run = [11, 21, 12, 30, 25]; 
% clus2run = [13, 10, 17, 18, 29];
% clus2run = [19, 20, 26, 28, 23];

%circ - still running - 30, 10
clus2run = 18; % missed
% clus2run = [23, 27]; %then 30, 10

clus2run = [10]; % TO RUN
clus2run = [30]; % TO RUN


%%%%%
% sq2trapz 1kiters - not started
% clus2run = [14, 15, 23,  20, 24, 16,  22, 27, 26,  13, 30, 29,  25, 17, 21, 12, 19, 28, 11, 10, 18]; %all
% clus2run = [14, 15, 23,  20, 24, 18, 16, 22, 27, 26, 19]; %half - 
% clus2run = [13, 30, 29,  25, 10, 28, 17, 11, 21, 12]; %half - 

% clus2run = 11;

jointTrls=1; %for test trials

% if trapKfrmSq1
if strcmp(dat,'trapzKfrmSq1')
    nTrials=1000000/2;
    epsMuTrapz10 = 25; %this is 10% of orig learning rate - using this
%     epsMuTrapz10 = 50; % anEps - 20%
%     epsMuTrapz10 = 30; % anEps -
    % epsMuTrapz10 = 15;  %
end


%doPerm or not
if ~strcmp(dat(1:4),'trap')
    doPerm=1;
else
    doPerm=0;
end

if nIter==1000
    doPerm=0; %if 1000 iters
end
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
            
            if ~(length(dat)==12)
                fname = sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_jointTrls_stepSiz',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat);
            else %trapzKfrmSq1 or 2 - i think only 1 now
                if strcmp(dat(12),'1')
                    fname = sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_epsMuTrapz_%d_jointTrls_stepSiz',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat,epsMuTrapz10);
                elseif strcmp(dat(12),'2')
                    fname = sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_jointTrls_stepSiz',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,'square');
                end
            end

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
%                     fname = [fname(1:end-1), sprintf('_actNorm_perm_%dpermsOn%diters',nPerm,nIters2run)]; %now correct actNorm; add 'trlsTest'?
                    fname = [[f(1).folder '/' f(1).name(1:end-11)], sprintf('_actNorm_perm_%dpermsOn%diters',nPerm,nIters2run)]; %now correct actNorm; add 'trlsTest'?
                else
                    fname = [[f(1).folder '/' f(1).name(1:end-11)], '_trlsTest_noPerm'];
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
