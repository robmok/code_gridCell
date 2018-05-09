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
%dat = 'circ';
% dat = 'trapzKrupic';
% dat = 'trapzKfrmSq1'; % run covering map on sq, then assess gridness in trapz

saveDat=1;

%for loading
% nTrials = 1000000;
% clus2run = [3:30]; 
% clus2run=[3:10 12:2:26 11:2:25, 27:30]; % 27:30

%batchSize=100
% clus2run=[4:2:10, 12:2:26]; 
% clus2run=[3:2:9, 11:2:25]; 

% %love06 - batchSize=400 - circ; actNorm
clus2run=[4, 22, 8, 24]; 
clus2run=[25, 10, 7, 20];  
clus2run=[5,  15, 26, 9]; 

% love01
clus2run=[3, 13]; 
% clus2run=[17, 23]; 
% clus2run=[18, 21]; 
% clus2run=[11, 12]; 
% clus2run=[19, 14]; 
% clus2run=[6, 16]; 


% clus2run = 20;

%love01 - sq - batchSize=400;actNorm
%for now run a few on love01
 clus2run=[4, 22]; %8, 24 
%  clus2run=[25, 10]; % 7, 20 
% clus2run=[5,  15];  %26, 9


%old
% clus2run=[4, 22, 16]; 
% % clus2run=[8, 12, 24];
% clus2run=[25, 10, 7];
% clus2run=[20, 14, 6];
% clus2run=[5,  15, 26];
% clus2run=[9, 19, 11];
% clus2run=[3, 13, 18];
% clus2run=[17, 23, 21];

%trapz - trapzKfrmSq
% clus2run=[8:4:24, 6]; 
% clus2run=[10:4:26,4]; 
% clus2run=[3, 7:4:25]; 
% clus2run=[5, 9:4:25]; 

% %split into 4 - running trapzKfrmSq1 love06
% clus2run=[8:8:24]; 
% clus2run=[12, 20,6];    
% clus2run=[10:8:26]; 
% clus2run=[14, 22,4];
% 
% %odd
% clus2run=[7:4:19]; 
% clus2run=[9:4:21]; 
% clus2run=[23,25];
% clus2run = [3,5];

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


%doPerm or not
if ~strcmp(dat(1:4),'trap')
    doPerm=1;
else
    doPerm=0; 
end
% doPerm=0;

% run perm tests on how many iters? takes a bit of time (a couple mins) per
% iter, so with 200 iters plus many conditions, maybe too much (if all the
% perm data are about the same, then just take max, or 95th percentile as
% the threshold value)
nIters2run = 200; 
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
            
            %load
%             fname = sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wAct_jointTrls_stepSiz',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat);
            %new
            if ~(length(dat)==12) 
                fname = sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_jointTrls_stepSiz',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat);
            else %trapzKfrmSq1 or 2
                fname = sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_jointTrls_stepSiz',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,'square');
            end
                       
%             if annEps %epsMu is different here
%                 fname = sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps*_batchSiz%d_%diters_%s_wAct_annEps_jointTrls_stepSiz',nClus,round(nTrials/1000),batchSize,nIter,dat);
%             end

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
                    [~,~,~,~,gA_act_sq,gA_actNorm_sq,gW_act_sq,gW_actNorm_sq] = gridnessTestData_Perm(densityPlot,'square',locRange,nClus,nTrialsTest,nPerm,nIters2run,0);
                    gA_sq = gA;
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
                    fname = [fname(1:end-1), sprintf('_actNorm_perm_%dpermsOn%diters',nPerm,nIters2run)]; %now correct actNorm
                else
                    fname = [fname(1:end-1), '_noPerm'];
                end
                if length(dat)==12
                    fname = [fname sprintf('_%s',dat)];
                end
                if ~jointTrls
                    fname = [fname '_noJointTrlsTest'];
                end
                
                fname = [fname sprintf('_%s',cTime)];
                

                if saveDat 
                    if length(dat)==12
                        save(fname,'gA_sq','gA_act_sq','gA_actNorm_sq','gW_act_sq','gW_actNorm_sq','gA_act','gA_actNorm','gW_act','gW_actNorm','densityPlotAct','densityPlotActNorm','timeTaken');
                    else
                        %                     save(fname,'permPrc_gA_act','permPrc_gA_act','permPrc_gA_actNorm','permPrc_gA_actNorm','gA_act','gA_actNorm','gW_act','gW_actNorm','timeTaken');
                        %also save some density plots for figures
                        save(fname,'permPrc_gA_act','permPrc_gW_act','permPrc_gA_actNorm','permPrc_gW_actNorm','gA_act','gA_actNorm','gW_act','gW_actNorm','densityPlotAct','densityPlotActNorm','timeTaken');
                    end
                end
                clear densityPlotAct densityPlotActNorm
            end
        end
    end
end
