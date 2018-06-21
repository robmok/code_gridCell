%load in square, split in half, assess gridness
clear all;

% wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
cd(wd);

codeDir = [wd '/code_gridCell'];
saveDir = [wd '/data_gridCell'];
% saveDir = [wd '/data_gridCell/toRmv'];
addpath(codeDir); addpath(saveDir);
addpath(genpath([codeDir '/gridSCORE_packed']));

figsDir = [wd '/grid_figs'];

% load
gaussSmooth = 1; 
fixBatchSize = 1; %fixed batch size or depend on nClus (for fname)

dat='square';
boxSize=1;

locRange = [0, 49];

% joined trials
jointTrls=1;
epsMuVals=.025;
nTrials=1000000;
% batchSizeVals = [1000, 400, 100]; 
annEps=0;
nIter=200;

% nIter=1000;

clus2run = [3:26]; 
% clus2run = [3:30]; 

batchSizeVals = 400; %100, 125, 200,400, 1000
% batchSizeVals = 200;

%new - slower learning rate
% epsMuVals=.015;
% batchSizeVals = 100; %100, 125, 200, 400
% clus2run = [12, 16, 24, 28]; %batchSize200 missed 20?

rHex=0; %if choose raw 60deg corr values, not gridness

doPerm = 0;
%load loop
for iClus2run = 1:length(clus2run) 
    nClus = clus2run(iClus2run);
    for iEps = 1:length(epsMuVals) 
        epsMuOrig=epsMuVals(iEps);
        epsMuOrig1000=epsMuOrig*1000;
        for iBvals = 1:length(batchSizeVals)
            batchSize = batchSizeVals(iBvals);
            fprintf('Loading %s, nClus=%d, epsMu=%d, batchSize=%d\n',dat,nClus,epsMuOrig1000,batchSize)
            if doPerm
%                 fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wAct_jointTrls_stepSiz_perm_%dpermsOn%diters',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat,nPerm,nIters2run)];
                fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_jointTrls_stepSiz_actNorm_perm_%dpermsOn%diters',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat,nPerm,nIters2run)];
            else
                if ~(length(dat)==12)
                    fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_jointTrls_stepSiz_trlsTest_noPerm',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat)];
                else
                    if strcmp(dat(12),'1')
                    fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_epsMuTrapz10_%d_jointTrls_stepSiz_trlsTest_noPerm_%s',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat,epsMuTrapz10,dat)];
                    else
                    fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_jointTrls_stepSiz_trlsTest_noPerm_%s',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,'square',dat)];
                    end
                end
            end
                        
%             if annEps % +++







            fname = [saveDir, fname '*'];
            
            %edit if want to load more than one file per sim, merge
            f = dir(fname); filesToLoad = cell(1,length(f));
            for iF = 1%:length(f)
                filesToLoad{iF} = f(iF).name;
                load(f(iF).name);
            end
            
%             %organise gridness values (allen vs willis method)
%             gA_gAll_act(:,iEps,iBvals,iClus2run,:)   = gA_act(:,1,:);
%             gA_oAll_act(:,iEps,iBvals,iClus2run,:)   = gA_act(:,2,:);
%             gA_radAll_act(:,iEps,iBvals,iClus2run,:) = gA_act(:,3,:);
%             gA_wavAll_act(:,iEps,iBvals,iClus2run,:) = gA_act(:,4,:);
%             gW_gAll_act(:,iEps,iBvals,iClus2run,:) = gW_act(:,1,:);
%             gW_oAll_act(:,iEps,iBvals,iClus2run,:) = gW_act(:,2,:);
%             gW_radAll_act(:,iEps,iBvals,iClus2run,:) = gW_act(:,3,:);
%             gW_wavAll_act(:,iEps,iBvals,iClus2run,:) = gW_act(:,4,:);
%             %
%             gA_gAll_actNorm(:,iEps,iBvals,iClus2run,:)   = gA_actNorm(:,1,:);
%             gA_oAll_actNorm(:,iEps,iBvals,iClus2run,:)   = gA_actNorm(:,2,:);
%             gA_radAll_actNorm(:,iEps,iBvals,iClus2run,:) = gA_actNorm(:,3,:);
%             gA_wavAll_actNorm(:,iEps,iBvals,iClus2run,:) = gA_actNorm(:,4,:);
%             gW_gAll_actNorm(:,iEps,iBvals,iClus2run,:) = gW_actNorm(:,1,:);
%             gW_oAll_actNorm(:,iEps,iBvals,iClus2run,:) = gW_actNorm(:,2,:);
%             gW_radAll_actNorm(:,iEps,iBvals,iClus2run,:) = gW_actNorm(:,3,:);
%             gW_wavAll_actNorm(:,iEps,iBvals,iClus2run,:) = gW_actNorm(:,4,:);  
            
            gA_actNorm_sq = nan(nIter,9,3);
            gW_actNorm_sq = nan(nIter,9,3);
                gA_actNorm_sq(:,:,1) = gA_actNorm;
                gW_actNorm_sq(:,:,1) = gW_actNorm;
            for iterI = 1:nIter
%                 densityPlotActAll(:,:,iterI,iEps,iBvals,iClus2run) = densityPlotAct(:,:,iterI);
%                 densityPlotActNormAll(:,:,iterI,iEps,iBvals,iClus2run) = densityPlotActNorm(:,:,iterI);
                
                densityPlotTmp= densityPlotAct(:,:,iterI);
                %left half of box
                aCorrMap = ndautoCORR(densityPlotTmp(:,1:(locRange(2)+1)/2));
                [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
                gA_actNorm_sq(iterI,:,2) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
                [g,gdataA] = gridSCORE(aCorrMap,'wills',0);
                gW_actNorm_sq(iterI,:,2) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
                %right half of box
                aCorrMap = ndautoCORR(densityPlotTmp(:,(locRange(2)+1)/2+1:end));
                [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
                gA_actNorm_sq(iterI,:,3) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
                [g,gdataA] = gridSCORE(aCorrMap,'wills',0);
                gW_actNorm_sq(iterI,:,3) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
            end
            

        end 
    end
end


hist(squeeze(gA_actNorm_sq(:,1,2:3)),25)
hist(squeeze(gA_actNorm_sq(:,1,2)-gA_actNorm_sq(:,1,3)),25)
% mean(squeeze(gA_actNorm_sq(:,1,2)-gA_actNorm_sq(:,1,3)))
