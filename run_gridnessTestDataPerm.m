%% Produce activation maps on test set; do shuffling / permutation test
clear all;

% Set working directory
% wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
% wd='/home/robmok/Documents/Grid_cell_model'; %on love01

cd(wd);
codeDir = [wd '/code_gridCell'];  %where the code lives
saveDir = [wd '/data_gridCell'];  %where to save the output of the simulations
addpath(codeDir); addpath(saveDir);
addpath(genpath([codeDir '/gridSCORE_packed'])); % add path to code for computing grid measures

locRange = [0 49];
nTrialsTest = 100000; % orig nTrials/10

% Set environment
dat = 'circ';
% dat = 'square';
dat = 'trapzKfrmSq1'; % trapz from square

%set
saveDat=1;
nIter = 200;  %for doPerm (500 perms on each iter: 500*200)
nIter = 1000; %noPerm - if just to make activation maps

clus2run = 10:30;
epsMuVals=.025;
nTrials=1000000;
batchSizeVals=200;
annEps=1;
if annEps 
    epsMuVals=.25; %below multiplies this by 1k for file name
end

% if trapKfrmSq1
if strcmp(dat,'trapzKfrmSq1')
    nTrials=1000000/4;
    epsMuTrapz10 = 25;
end
jointTrls=1; %for test trials

%doPerm or not
if ~strcmp(dat(1:4),'trap')
    doPerm=1;
else
    doPerm=0;
end
if nIter==1000
    doPerm=0; %if 1000 iters
end

% number of perm tests
nIters2run = nIter; %200
nPerm = 500;

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
            
            if ~(length(dat)==12) %sq, circ
                fname = sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_jointTrls_stepSiz',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat);
            else %trapzKfrmSq1
                fname = sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_epsMuTrapz_%d_jointTrls_stepSiz',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat,epsMuTrapz10);
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
                %might have problem with overlapping file names if done more than once
                for iF = 1
                    load(f(iF).name);
                end
                                
                %run
                tic
                [permPrc_gA_act, permPrc_gW_act,permPrc_gA_actNorm, permPrc_gW_actNorm,gA_act,gA_actNorm,gW_act,gW_actNorm, gA_actNormPerm, gW_actNormPerm, densityPlotAct,densityPlotActNorm] = gridnessTestData_Perm(densityPlot,dat,locRange,nClus,nTrialsTest,nPerm,nIters2run,doPerm);
                timeTaken=toc;
                
                %save
                cTime=datestr(now,'HHMMSS');
                if doPerm
                    fname = [[f(1).folder '/' f(1).name(1:end-11)], sprintf('_actNorm_perm_%dpermsOn%diters',nPerm,nIters2run)];
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
                        %also save some density plots for figures
                        save(fname,'permPrc_gA_act','permPrc_gW_act','permPrc_gA_actNorm','permPrc_gW_actNorm','gA_act','gA_actNorm','gW_act','gW_actNorm','densityPlotAct','densityPlotActNorm','timeTaken');
                    end
                end
                clear densityPlotAct densityPlotActNorm gA_act gA_actNorm gW_act gW_actNorm permPrc_gA_act permPrc_gW_act
            end
        end
    end
end