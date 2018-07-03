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

saveDat = 1;

% joined trials
jointTrls=1;
epsMuVals=.025;
nTrials=1000000;
% batchSizeVals = [1000, 400, 100]; 
annEps=0;
nIter=200;

nIter=1000;

clus2run = 3:30; 

% clus2run = 10:30; 
% clus2run = [27:30]; 

batchSizeVals = 400; %100, 125, 200,400, 1000
% batchSizeVals = 200;

rHex=0; %if choose raw 60deg corr values, not gridness

doPerm = 0;


%% load in data

for iClus2run = 1:length(clus2run) 
    nClus = clus2run(iClus2run);
    for iEps = 1:length(epsMuVals) 
        epsMuOrig=epsMuVals(iEps);
        epsMuOrig1000=epsMuOrig*1000;
        for iBvals = 1:length(batchSizeVals)
            batchSize = batchSizeVals(iBvals);
            fprintf('Loading %s, nClus=%d, epsMu=%d, batchSize=%d\n',dat,nClus,epsMuOrig1000,batchSize)

            fname = [sprintf('/covering_map_batch_dat_sqSplitInHalf_gridness_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_jointTrls_stepSiz_trlsTest_noPerm.mat',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat)];
            
            
            % annEps ++
            
            
            
            
            fname = [saveDir, fname];
            load(fname);

            %organise gridness values (allen vs willis method)
            gA_gAll_act(:,iEps,iBvals,iClus2run,:)   = gA_act_sq(:,1,:);
            gA_oAll_act(:,iEps,iBvals,iClus2run,:)   = gA_act_sq(:,2,:);
            gA_radAll_act(:,iEps,iBvals,iClus2run,:) = gA_act_sq(:,3,:);
            gA_wavAll_act(:,iEps,iBvals,iClus2run,:) = gA_act_sq(:,4,:);
%             gW_gAll_act(:,iEps,iBvals,iClus2run,:) = gW_act_sq(:,1,:);
%             gW_oAll_act(:,iEps,iBvals,iClus2run,:) = gW_act_sq(:,2,:);
%             gW_radAll_act(:,iEps,iBvals,iClus2run,:) = gW_act_sq(:,3,:);
%             gW_wavAll_act(:,iEps,iBvals,iClus2run,:) = gW_act_sq(:,4,:);
            
            gA_gAll_actNorm(:,iEps,iBvals,iClus2run,:)   = gA_actNorm_sq(:,1,:);
            gA_oAll_actNorm(:,iEps,iBvals,iClus2run,:)   = gA_actNorm_sq(:,2,:);
            gA_radAll_actNorm(:,iEps,iBvals,iClus2run,:) = gA_actNorm_sq(:,3,:);
            gA_wavAll_actNorm(:,iEps,iBvals,iClus2run,:) = gA_actNorm_sq(:,4,:);
%             gW_gAll_actNorm(:,iEps,iBvals,iClus2run,:) = gW_actNorm_sq(:,1,:);
%             gW_oAll_actNorm(:,iEps,iBvals,iClus2run,:) = gW_actNorm_sq(:,2,:);
%             gW_radAll_actNorm(:,iEps,iBvals,iClus2run,:) = gW_actNorm_sq(:,3,:);
%             gW_wavAll_actNorm(:,iEps,iBvals,iClus2run,:) = gW_actNorm_sq(:,4,:);
            
            for iterI = 1:nIter
               densityPlotActNormAll(:,:,iterI,iEps,iBvals,iClus2run) = densityPlotActNorm(:,:,iterI);
            end
            
        end
    end
end

%% load in, split in half (no need if loading in)
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
            
            gA_act_sq = nan(nIter,9,3);
            gW_act_sq = nan(nIter,9,3);
            gA_act_sq(:,:,1) = gA_act;
            gW_act_sq(:,:,1) = gW_act;
            gA_actNorm_sq = nan(nIter,9,3);
            gW_actNorm_sq = nan(nIter,9,3);
            gA_actNorm_sq(:,:,1) = gA_actNorm;
            gW_actNorm_sq(:,:,1) = gW_actNorm;
            for iterI = 1:nIter
%                 densityPlotActAll(:,:,iterI,iEps,iBvals,iClus2run) = densityPlotAct(:,:,iterI);
%                 densityPlotActNormAll(:,:,iterI,iEps,iBvals,iClus2run) = densityPlotActNorm(:,:,iterI);
                
                %act
                densityPlotTmp= densityPlotAct(:,:,iterI);
                %left half of box
                aCorrMap = ndautoCORR(densityPlotTmp(:,1:(locRange(2)+1)/2));
                [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
                gA_act_sq(iterI,:,2) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
                [g,gdataA] = gridSCORE(aCorrMap,'wills',0);
                gW_act_sq(iterI,:,2) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
                %right half of box
                aCorrMap = ndautoCORR(densityPlotTmp(:,(locRange(2)+1)/2+1:end));
                [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
                gA_act_sq(iterI,:,3) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
                [g,gdataA] = gridSCORE(aCorrMap,'wills',0);
                gW_act_sq(iterI,:,3) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
                
                %actNorm
                densityPlotTmp= densityPlotActNorm(:,:,iterI);
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
            
            %organise gridness values (allen vs willis method)
            gA_gAll_act(:,iEps,iBvals,iClus2run,:)   = gA_act_sq(:,1,:);
            gA_oAll_act(:,iEps,iBvals,iClus2run,:)   = gA_act_sq(:,2,:);
            gA_radAll_act(:,iEps,iBvals,iClus2run,:) = gA_act_sq(:,3,:);
            gA_wavAll_act(:,iEps,iBvals,iClus2run,:) = gA_act_sq(:,4,:);
            gW_gAll_act(:,iEps,iBvals,iClus2run,:) = gW_act_sq(:,1,:);
            gW_oAll_act(:,iEps,iBvals,iClus2run,:) = gW_act_sq(:,2,:);
            gW_radAll_act(:,iEps,iBvals,iClus2run,:) = gW_act_sq(:,3,:);
            gW_wavAll_act(:,iEps,iBvals,iClus2run,:) = gW_act_sq(:,4,:);
            %
            gA_gAll_actNorm(:,iEps,iBvals,iClus2run,:)   = gA_actNorm_sq(:,1,:);
            gA_oAll_actNorm(:,iEps,iBvals,iClus2run,:)   = gA_actNorm_sq(:,2,:);
            gA_radAll_actNorm(:,iEps,iBvals,iClus2run,:) = gA_actNorm_sq(:,3,:);
            gA_wavAll_actNorm(:,iEps,iBvals,iClus2run,:) = gA_actNorm_sq(:,4,:);
            gW_gAll_actNorm(:,iEps,iBvals,iClus2run,:) = gW_actNorm_sq(:,1,:);
            gW_oAll_actNorm(:,iEps,iBvals,iClus2run,:) = gW_actNorm_sq(:,2,:);
            gW_radAll_actNorm(:,iEps,iBvals,iClus2run,:) = gW_actNorm_sq(:,3,:);
            gW_wavAll_actNorm(:,iEps,iBvals,iClus2run,:) = gW_actNorm_sq(:,4,:);  
            
            fname = [fname(1:92) sprintf('_sqSplitInHalf_gridness_%d',nClus) fname(end-92:end-1)]; %92 if 1000iters, 91 if 200 iters..

            %save
            if saveDat
                save(fname,'gA_act_sq', 'gA_actNorm_sq','gW_act_sq', 'gW_actNorm_sq', 'densityPlotAct', 'densityPlotActNorm')
            end
        end 
    end
end

% hist(squeeze(gA_actNorm_sq(:,1,2:3)),25)
% hist(squeeze(gA_actNorm_sq(:,1,2)-gA_actNorm_sq(:,1,3)),25)
% mean(squeeze(gA_actNorm_sq(:,1,2)-gA_actNorm_sq(:,1,3)))


%% sq L-R gridness

savePlots=0;

clusPosAct = 'actNorm'; %'act' or 'actNorm'

gridMsrType = 'a'; % 'a' or 'w' for allen or willis method - a preferred

gridMeasure = 'grid';

switch clusPosAct
%     case 'act'
%         switch gridMsrType
%             case 'a'
%                 gridness    = gA_gAll_act;
%                 orientation = gA_oAll_act;
%                 rad         = gA_radAll_act;
%                 wav         = gA_wavAll_act;
%             case 'w'
%                 gridness    = gW_gAll_act;
%                 orientation = gW_oAll_act;
%                 rad         = gW_radAll_act;
%                 wav         = gW_wavAll_act;
%         end
    case 'actNorm'
        switch gridMsrType
            case 'a'
                gridness    = gA_gAll_actNorm(:,:,:,:,2)-gA_gAll_actNorm(:,:,:,:,3);
            case 'w'
                gridness    = gW_gAll_actNorm(:,:,:,:,2)-gW_gAll_actNorm(:,:,:,:,3);
        end
end

switch gridMeasure
    case 'grid'
        datTmp=gridness;
    case 'angle'
        datTmp=orientation;
    case 'rad'
        datTmp=rad;
    case 'wav'
        datTmp=wav;
end


iEps=1;


clus2plot=(3:30)-2;
clus2plot=(6:30)-2;
clus2plot=(10:30)-2;
clus2plot=(10:26)-2;

% clus2plot=(3:26)-2;

iBatchVals=1;

%fig specs
% xTickLabs = num2cell(clus2run(clus2plot));
xTickLabs = num2cell((clus2plot));
fontSiz=15;

figure; hold on;
    for iEps = 1:length(epsMuVals)
        dat1     = squeeze(datTmp(:,iEps,iBatchVals,clus2plot,1));
%         dat1     = squeeze(datTmp(:,iEps,iBatchVals,:,1));
        barpos  = .25:.5:.5*size(dat1,2);
        colors  = distinguishable_colors(size(dat1,2));
        colgrey = [.6, .6, .6];
        mu      = mean(dat1,1);
        sm      = std(dat1)./sqrt(size(dat1,1));
        ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
        plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
%         errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
%         scatter(barpos,mu,50,colors,'filled','d');
        scatter(barpos,mu,50,colgrey,'filled','d');
        xticklabels(xTickLabs);
        xlim([barpos(1)-.5, barpos(end)+.5]);
        %         ylim([0,1]);
        if strcmp(gridMsrType,'a')
            ylim([-1.5,2]);
        elseif strcmp(gridMsrType,'w')
%             ylim([-1.25,1.4]);
        end
        xlabel('Number of Clusters');
        ylabel('Grid Score');
        title('Left-Right Square box')

        set(gca,'FontSize',fontSiz,'fontname','Arial')
    end

    fname = [figsDir sprintf('/gridness_%s_univarScatters_testSet_nClus%d-%d_eps%d_batchSiz%d_%s_%s','sqL-R',clus2plot(1)+2,clus2plot(end)+2,epsMuVals(iEps)*1000,batchSizeVals(iBatchVals),clusPosAct,gridMsrType)];

if savePlots
   set(gcf,'Renderer','painters');
   print(gcf,'-depsc2',fname)
   saveas(gcf,fname,'png');
end

% mean & bootstrap 95 CIs
nBoot = nIter;
clear ciTest ciTmp ciClusSig
cnter=0;
for iClus=1:length(clus2plot)
    cnter = cnter+1;
    ciTmp = bootci(nBoot,@nanmean,dat1(:,iClus));
    ciTest(:,cnter) = [nanmean(dat1(:,iClus),1); ciTmp]; % mean & CIs
    
    %check CIs sig
    posInd=ciTest(:,cnter)>0;
    negInd=ciTest(:,cnter)<0;
    if all(posInd)
       ciClusSig(cnter)=1;
    end
    if all(negInd)
        ciClusSig(cnter)=-1;
    end
    if ~all(posInd) && ~all(negInd)
        ciClusSig(cnter)=0;
    end
end

ci = bootci(numel(dat1),@nanmean,reshape(dat1,1,numel(dat1))); %CI over ALL runs
fprintf('Sq_L-R, %d to %d clusters: mean=%0.4f; CI=[%0.4f,%0.4f]; %d sig L-R, %d sig R-L, %d n.s.\n',clus2plot(1)+2,clus2plot(end)+2,mean(reshape(dat1,1,numel(dat1))),ci(1),ci(2),nnz(ciClusSig==1),nnz(ciClusSig==-1),nnz(ciClusSig==0));
%% Making figs: density plot examples

savePlots = 0;

doPlot=0; %do plot when computing gridness

fontSiz=15;

clusPosAct = 'actNorm'; %'clus' or 'actNorm'

close all
gridMsrType = 'a';

clus2plot = [7,10,12,25]-2;%[5,8,13,18];
% clus2plot = 10-2;
%
clus2plot = [10,12,18,25]-2;
% clus2plot = [7,10,12,18,25]-2;

% clus2plot = ([7,9,10,11,12,18,20,23,25,28])-2;
clus2plot = ([10,18,20])-2;

% clus2plot = (3:5)-2;


for iClus = clus2plot%:length(clus2run)
for iBvals = 1:length(batchSizeVals)
    for iterI = 1%nIter
        switch clusPosAct
            case 'act'
                densityPlotCentresSm = densityPlotActAll(:,:,iterI,iEps,iBvals,iClus);
                if strcmp(dat(1:2),'ci')
                    zlims=[0,6];
                elseif strcmp(dat(1:2),'sq')
                    zlims=[0,3.5];
                elseif strcmp(dat(1:2),'tr')
                    zlims=[0,11.5];
                end
            case 'actNorm'
                densityPlotCentresSm = densityPlotActNormAll(:,:,iterI,iEps,iBvals,iClus);
                zlims=[0,.08];
                if strcmp(dat(1:2),'ci') || strcmp(dat(1:2),'tr') 
                zlims=[0,.08];
                elseif strcmp(dat(1:2),'sq')
                zlims=[0,.08];
                end
        end

        aCorrMap=ndautoCORR(densityPlotCentresSm); %autocorrelogram        
        [g,gdataA] = gridSCORE(aCorrMap,'allen',doPlot);
%         [g,gdataW] = gridSCORE(aCorrMap,'wills',doPlot);
                
%         %make nans into yellow to edit out
        densityPlotTmp = densityPlotCentresSm;
        densityPlotTmp(isnan(densityPlotCentresSm))=2;
        aCorrMap(isnan(aCorrMap))=.5;
        
        figure; imagesc(densityPlotTmp,zlims); %colorbar;
%         figure; imagesc(densityPlotTmp); colorbar;
        xticks([]); xticklabels({''}); yticks([]); yticklabels({''});
        fname = [figsDir sprintf('/densityPlot_testSet_%s_nClus%d_eps%d_batchSiz%d_%s_%s_iter%d_clusters',dat,iClus+2,epsMuVals(iEps)*1000,batchSizeVals(iBvals),clusPosAct,gridMsrType,iterI)];
        if savePlots
            set(gcf,'Renderer','painters');
            print(gcf,'-depsc2',fname)
            saveas(gcf,fname,'png');
           close all
        end
        figure; imagesc(aCorrMap,[-1.01 1.01]);
%         figure; imagesc(aCorrMap); colorbar
        title(sprintf('Grid Score %.2f',g));
        set(gca,'FontSize',fontSiz,'fontname','Arial')
        xticks([]); xticklabels({''}); yticks([]); yticklabels({''});
        fname = [figsDir sprintf('/densityPlot_testSet_%s_nClus%d_eps%d_batchSiz%d_%s_%s_iter%d_autoCorr',dat,iClus+2,epsMuVals(iEps)*1000,batchSizeVals(iBvals),clusPosAct,gridMsrType,iterI)];
        if savePlots
            set(gcf,'Renderer','painters');
            print(gcf,'-depsc2',fname)
            saveas(gcf,fname,'png');
           close all
        end
        
            aCorrMap = ndautoCORR(densityPlotCentresSm(:,1:(locRange(2)+1)/2));
            [g,gdataA] = gridSCORE(aCorrMap,'allen',doPlot);
            %         [g,gdataW] = gridSCORE(aCorrMap,'wills',doPlot);
            
            aCorrMap(isnan(aCorrMap))=.5;
            figure; imagesc(aCorrMap,[-1.01 1.01]);
            title(sprintf('Grid Score %.2f',g));
            set(gca,'FontSize',fontSiz,'fontname','Arial')
            xticks([]); xticklabels({''}); yticks([]); yticklabels({''});
            fname = [figsDir sprintf('/densityPlot_testSet_%s_L_nClus%d_eps%d_batchSiz%d_%s_%s_iter%d__autoCorr',dat,iClus+2,epsMuVals(iEps)*1000,batchSizeVals(iBvals),clusPosAct,gridMsrType,iterI)];
            if savePlots
                set(gcf,'Renderer','painters');
                print(gcf,'-depsc2',fname)
                saveas(gcf,fname,'png');
                close all
            end
            
            aCorrMap = ndautoCORR(densityPlotCentresSm(:,(locRange(2)+1)/2+1:end));

            [g,gdataA] = gridSCORE(aCorrMap,'allen',doPlot);
            %         [g,gdataW] = gridSCORE(aCorrMap,'wills',doPlot);
            aCorrMap(isnan(aCorrMap))=.5;
            figure; imagesc(aCorrMap,[-1.01 1.01]);
            title(sprintf('Grid Score %.2f',g));
            set(gca,'FontSize',fontSiz,'fontname','Arial')
            xticks([]); xticklabels({''}); yticks([]); yticklabels({''});
            fname = [figsDir sprintf('/densityPlot_testSet_%s_R_nClus%d_eps%d_batchSiz%d_%s_%s_iter%d__autoCorr',dat,iClus+2,epsMuVals(iEps)*1000,batchSizeVals(iBvals),clusPosAct,gridMsrType,iterI)];
            if savePlots
                set(gcf,'Renderer','painters');
                print(gcf,'-depsc2',fname)
                saveas(gcf,fname,'png');
                close all
            end
    end
end
end