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
nSet        = 22;
gaussSmooth = 1; 
fixBatchSize = 1; %fixed batch size or depend on nClus (for fname)

dat='circ';
% dat='square';
compareCircSq = 0; %if compare circ-sq gridness, dat=circ, but load sq too

annEps=0;
boxSize=1;

% joined trials
jointTrls=1;
% clus2run = [8, 12, 16, 20,24, 28]; 
% clus2run = [8:2:28]; 
epsMuVals=.025;
nTrials=1000000;
% batchSizeVals = [1000, 400, 100]; 
% batchSizeVals=400;
nIter=200;

nIter=1000;


%load perm data when nIter=1000- this is with nIter=200. doPerm should=0
loadPerm = 1;

% clus2run = [3:16, 18, 20:26];  %missed 17, 19?

% % dat='trapzKrupic';    
% dat='trapzKfrmSq1';
% %trapzKfrmSq1
% nTrials=1000000/2;
% epsMuTrapz10 = 25; 


% clus2run = [3:26]; 
clus2run = [3:30]; 
% clus2run = [3:26]; 
clus2run = 10:30;

batchSizeVals = 400; %100, 125, 200,400, 1000
% batchSizeVals = 200;

%new - slower learning rate
% epsMuVals=.015;
% batchSizeVals = 100; %100, 125, 200, 400
% clus2run = [12, 16, 24, 28]; %batchSize200 missed 20?



%new annEps
annEps=1;
% epsMuVals = 0.1;
epsMuVals = 0.25;


rHex=0; %if choose raw 60deg corr values, not gridness

%perm
nPerm=500;
nIters2run=nIter;

doPerm = 0; 

if strcmp(dat(1:4),'trap')
    doPerm=0;
end

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
                if ~annEps
                    fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_jointTrls_stepSiz_actNorm_perm_%dpermsOn%diters',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat,nPerm,nIters2run)];
                else
                    %check
                    fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_jointTrls_stepSiz_*annEps_actNorm_perm_%dpermsOn%diters',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat,nPerm,nIters2run)];
                end
            else
                if ~(length(dat)==12)
                    if ~annEps
                        fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_jointTrls_stepSiz_trlsTest_noPerm',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat)];
                    else
                        fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_jointTrls_stepSiz_*annEps_trlsTest_noPerm',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat)];
                    end
                else
                    if strcmp(dat(12),'1')
                        fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_epsMuTrapz10_%d_jointTrls_stepSiz_trlsTest_noPerm_%s',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat,epsMuTrapz10,dat)];
                    else
                        fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_jointTrls_stepSiz_trlsTest_noPerm_%s',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,'square',dat)];
                    end
                end
            end

            %finish with directory and * for date/time
            fname = [saveDir, fname '*']; %finish with directory and * for date/time

            
            %edit if want to load more than one file per sim, merge
            f = dir(fname); filesToLoad = cell(1,length(f));
            for iF = 1%:length(f)
                filesToLoad{iF} = f(iF).name;
                load(f(iF).name);
            end
            
            
            %organise gridness values (allen vs willis method)
            gA_gAll_act(:,iEps,iBvals,iClus2run,:)   = gA_act(:,1,:);
            gA_oAll_act(:,iEps,iBvals,iClus2run,:)   = gA_act(:,2,:);
            gA_radAll_act(:,iEps,iBvals,iClus2run,:) = gA_act(:,3,:);
            gA_wavAll_act(:,iEps,iBvals,iClus2run,:) = gA_act(:,4,:);
            gW_gAll_act(:,iEps,iBvals,iClus2run,:) = gW_act(:,1,:);
            gW_oAll_act(:,iEps,iBvals,iClus2run,:) = gW_act(:,2,:);
            gW_radAll_act(:,iEps,iBvals,iClus2run,:) = gW_act(:,3,:);
            gW_wavAll_act(:,iEps,iBvals,iClus2run,:) = gW_act(:,4,:);
            %
            gA_gAll_actNorm(:,iEps,iBvals,iClus2run,:)   = gA_actNorm(:,1,:);
            gA_oAll_actNorm(:,iEps,iBvals,iClus2run,:)   = gA_actNorm(:,2,:);
            gA_radAll_actNorm(:,iEps,iBvals,iClus2run,:) = gA_actNorm(:,3,:);
            gA_wavAll_actNorm(:,iEps,iBvals,iClus2run,:) = gA_actNorm(:,4,:);
            gW_gAll_actNorm(:,iEps,iBvals,iClus2run,:) = gW_actNorm(:,1,:);
            gW_oAll_actNorm(:,iEps,iBvals,iClus2run,:) = gW_actNorm(:,2,:);
            gW_radAll_actNorm(:,iEps,iBvals,iClus2run,:) = gW_actNorm(:,3,:);
            gW_wavAll_actNorm(:,iEps,iBvals,iClus2run,:) = gW_actNorm(:,4,:);  
            
            
            for iterI = 1:nIter
%                 densityPlotActAll(:,:,iterI,iEps,iBvals,iClus2run) = imgaussfilt(densityPlotAct(:,:,iterI),gaussSmooth);
%                 densityPlotActNormAll(:,:,iterI,iEps,iBvals,iClus2run) = imgaussfilt(densityPlotActNorm(:,:,iterI),gaussSmooth);
                densityPlotActAll(:,:,iterI,iEps,iBvals,iClus2run) = densityPlotAct(:,:,iterI);
                densityPlotActNormAll(:,:,iterI,iEps,iBvals,iClus2run) = densityPlotActNorm(:,:,iterI);
            end
            
            %load perm data when nIter=1000- this is with nIter=200.
            if loadPerm
                nIterOrig = nIter; nIter=200; nIters2run=nIter;
                if ~annEps
                    fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_jointTrls_stepSiz_actNorm_perm_%dpermsOn%diters',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat,nPerm,nIters2run)];
                else
                    fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_jointTrls_stepSiz_annEps_actNorm_perm_%dpermsOn%diters',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat,nPerm,nIters2run)];
                end
                fname = [saveDir, fname '*'];
                f = dir(fname); filesToLoad = cell(1,length(f));
                filesToLoad{1} = f(1).name;
                perm=load(f(1).name);
                nIter = nIterOrig; nIters2run = nIter;
            end
            
            if doPerm || loadPerm
%                 gA_act_permPrc(:,iEps,iBvals,iClus2run,:)   = permPrc_gA_act(:,3,:);
%                 gW_act_permPrc(:,iEps,iBvals,iClus2run,:)   = permPrc_gW_act(:,3,:);
                gA_actNorm_permPrc(:,iEps,iBvals,iClus2run,:)   = perm.permPrc_gA_actNorm(:,3,:);
                gW_actNorm_permPrc(:,iEps,iBvals,iClus2run,:)   = perm.permPrc_gW_actNorm(:,3,:);
            end
            
        end 
    end
end


if strcmp(dat,'trapzKfrmSq1') || compareCircSq %added compare circ vs sq - should work?
    
    datOrig = dat;
    dat = 'square';
    doPerm=0; % sq 200iter noPerm not run yet
    nTrials=1000000;

%load sq gridness to compare with trapz
for iClus2run = 1:length(clus2run) 
    nClus = clus2run(iClus2run);
    for iEps = 1:length(epsMuVals) 
        epsMuOrig=epsMuVals(iEps);
        epsMuOrig1000=epsMuOrig*1000;
        for iBvals = 1:length(batchSizeVals)
            batchSize = batchSizeVals(iBvals);
            fprintf('Loading %s, nClus=%d, epsMu=%d, batchSize=%d\n',dat,nClus,epsMuOrig1000,batchSize)

            if doPerm
                fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_jointTrls_stepSiz_actNorm_perm_%dpermsOn%diters',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat,nPerm,nIters2run)];
            else
                fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_jointTrls_stepSiz_trlsTest_noPerm',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat)];
            end
%             %note, annEps slightly different file name




            %finish with directory and * for date/time
            fname = [saveDir, fname '*'];
            
            %edit if want to load more than one file per sim, merge
            f = dir(fname); filesToLoad = cell(1,length(f));
            for iF = 1%:length(f)
                filesToLoad{iF} = f(iF).name;
                load(f(iF).name);
            end
            
            %organise gridness values (allen vs willis method)
            gA_gAll_actNormSq(:,iEps,iBvals,iClus2run,:)   = gA_actNorm(:,1,:);
            gA_oAll_actNormSq(:,iEps,iBvals,iClus2run,:)   = gA_actNorm(:,2,:);
            gA_radAll_actNormSq(:,iEps,iBvals,iClus2run,:) = gA_actNorm(:,3,:);
            gA_wavAll_actNormSq(:,iEps,iBvals,iClus2run,:) = gA_actNorm(:,4,:);
            gW_gAll_actNormSq(:,iEps,iBvals,iClus2run,:) = gW_actNorm(:,1,:);
            gW_oAll_actNormSq(:,iEps,iBvals,iClus2run,:) = gW_actNorm(:,2,:);
            gW_radAll_actNormSq(:,iEps,iBvals,iClus2run,:) = gW_actNorm(:,3,:);
            gW_wavAll_actNormSq(:,iEps,iBvals,iClus2run,:) = gW_actNorm(:,4,:);  
            
        end 
    end
end
dat = datOrig;
end
%% Making figs - univar scatters 1 - test set

savePlots=0;

clusPosAct = 'actNorm'; %'act' or 'actNorm'

gridMsrType = 'a'; % 'a' or 'w' for allen or willis method - a preferred

gridMeasure = 'grid';

plotFewClus = 0; %plot 3:5 clusters separately

computeCIs = 1; %takes a bit of time

switch clusPosAct
% case 'clus'
%     switch gridMsrType
%         case 'a'
%             gridness    = gA_gAll;
%             orientation = gA_oAll;
%             rad         = gA_radAll;
%             wav         = gA_wavAll;
%         case 'w'
%             gridness    = gW_gAll;
%             orientation = gW_oAll;
%             rad         = gW_radAll;
%             wav         = gW_wavAll;
%     end
    case 'act'
        switch gridMsrType
            case 'a'
                gridness    = gA_gAll_act;
                orientation = gA_oAll_act;
                rad         = gA_radAll_act;
                wav         = gA_wavAll_act;
            case 'w'
                gridness    = gW_gAll_act;
                orientation = gW_oAll_act;
                rad         = gW_radAll_act;
                wav         = gW_wavAll_act;
        end
    case 'actNorm'
        switch gridMsrType
            case 'a'
                gridness    = gA_gAll_actNorm;
                orientation = gA_oAll_actNorm;
                rad         = gA_radAll_actNorm;
                wav         = gA_wavAll_actNorm;
            case 'w'
                gridness    = gW_gAll_actNorm;
                orientation = gW_oAll_actNorm;
                rad         = gW_radAll_actNorm;
                wav         = gW_wavAll_actNorm;
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

% clus2plot=(3:30)-2;
% clus2plot=(6:30)-2;
% clus2plot=(10:30)-2;
% clus2plot=(10:26)-2;

clus2plot=1:length(clus2run);

%trapz
% clus2plot = 1:length(clus2run);

iBatchVals=1;

%fig specs
xTickLabs = num2cell(clus2run(clus2plot));
fontSiz=15;

figure; hold on;
    for iEps = 1:length(epsMuVals)
        %     subplot(2,3,iEps);
        dat1     = squeeze(datTmp(:,iEps,iBatchVals,clus2plot,1));
        barpos  = .25:.5:.5*size(dat1,2);
        colors  = distinguishable_colors(size(dat1,2));
        colgrey = [.6, .6, .6];
        mu      = mean(dat1,1);
        sm      = std(dat1)./sqrt(size(dat1,1));
%         ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
        plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
%         errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
%         scatter(barpos,mu,50,colors,'filled','d');
        scatter(barpos,mu,50,colgrey,'filled','d');
        xticklabels(xTickLabs);
        xlim([barpos(1)-.5, barpos(end)+.5]);
        %         ylim([0,1]);
        if strcmp(gridMsrType,'a')
%             ylim([-.45,1.4]);
        elseif strcmp(gridMsrType,'w')
            ylim([-1.25,1.4]);
        end
        xlabel('Number of Clusters');
        ylabel('Grid Score');
        
%         title(sprintf('%s, %s - eps=%d, batchSize=%d',dat, gridMeasure,epsMuVals(iEps)*1000,batchSizeVals(iBatchVals)))
        if strcmp(dat(1:2),'ci')
            title('Circular Environment')
        elseif strcmp(dat(1:2),'sq')
            title('Square Environment')
        elseif strcmp(dat(1:2),'tr')
            title('Trapezoid')
        end
        set(gca,'FontSize',fontSiz,'fontname','Arial')
    end

    fname = [figsDir sprintf('/gridness_%s_univarScatters_testSet_nClus%d-%d_eps%d_batchSiz%d_%s_%s',dat,clus2run(clus2plot(1)),clus2run(clus2plot(end)),epsMuVals(iEps)*1000,batchSizeVals(iBatchVals),clusPosAct,gridMsrType)];
if annEps
    fname = [fname '_annEps'];
end
if savePlots
   set(gcf,'Renderer','painters');
   print(gcf,'-depsc2',fname)
   saveas(gcf,fname,'png');
   close all
end


%prop grid cells
permPrc = 100; %use top x prctile of the perm distributions - 97.5/99/99.5/100
if ~strcmp(dat(1:4),'trap') && loadPerm
        propGrid = nan(size(dat1,2),1);
        propGridCI = nan(size(dat1,2),2);
        propGridMaxThresh = nan(size(dat1,2),1);
        propGridCIMaxThresh = nan(size(dat1,2),2);
        for i=1:size(dat1,2)
            propGrid(i)     = nnz(dat1(:,i)>squeeze(prctile(gA_actNorm_permPrc(:,iEps,iBvals,clus2plot(i)),permPrc)))./nIter; %cond-wise thresh
            propGridCI(i,:) = bootrm(dat1(:,i),[2.5, 97.5],'percentGr',nIter,prctile(gA_actNorm_permPrc(:,iEps,iBvals,clus2plot(i)),permPrc));
            
            propGridMaxThresh(i)     = nnz(dat1(:,i)>max((prctile(gA_actNorm_permPrc(:,iEps,iBvals,clus2plot),permPrc))))./nIter; %max thresh
            propGridCIMaxThresh(i,:) = bootrm(dat1(:,i),[2.5, 97.5],'percentGr',nIter,max((prctile(gA_actNorm_permPrc(:,iEps,iBvals,clus2plot),permPrc))));
        end
        fprintf('%d to %d clusters: \npropGrid mean across cond means: %.3f\n', clus2run(clus2plot(1)),clus2run(clus2plot(end)), mean(propGrid))
        fprintf('propGrid mean across cond means maxThresh: %.3f\n', mean(propGridMaxThresh))
        
        %prop Grid over all conditions with max threshold
        propGridAll     = nnz(dat1>max((prctile(gA_actNorm_permPrc(:,iEps,iBvals,clus2plot),permPrc))))./nnz(dat1);
        propGridAllCI = bootrm(reshape(dat1,numel(dat1),1),[2.5, 97.5],'percentGr',nIter,max(prctile(squeeze(gA_actNorm_permPrc(:,iEps,iBvals,clus2plot)),permPrc)));
        fprintf('propGrid mean and CI across cond all iters - max thresh: %.3f, CIs [%.3f,%.3f]\n', propGridAll, propGridAllCI)
    
    %plot mean + CI as errorbars
    figure; hold on;
    mu = propGrid;
    ci = propGridCI;
    barpos  = .25:.5:.5*size(dat1,2);
    colors  = distinguishable_colors(size(dat1,2));
    colgrey = [.6, .6, .6];
    scatter(barpos,mu,200,colors,'.');
    errorbar(barpos,mu,mu-ci(:,1),abs(mu-ci(:,2)),'Color',colgrey,'LineStyle','None','LineWidth',1);
    xticks(barpos);
    xticklabels(xTickLabs);
    xlim([barpos(1)-.5, barpos(end)+.5]);
    xlabel('Number of Clusters');
    ylabel('Proportion "Grid Cells"');
    set(gca,'FontSize',fontSiz,'fontname','Arial')

%             title('Circular box')
    fname = [figsDir sprintf('/gridness_%s_propGridCells_permPrc%d_testSet_nClus%d-%d_eps%d_batchSiz%d_%s_%s',dat,permPrc,clus2run(clus2plot(1)),clus2run(clus2plot(end)),epsMuVals(iEps)*1000,batchSizeVals(iBatchVals),clusPosAct,gridMsrType)];
    
    if savePlots
        set(gcf,'Renderer','painters');
        print(gcf,'-depsc2',fname)
        saveas(gcf,fname,'png');
        close all
    end
end


if strcmp(dat(1:4),'trap')
    figure; hold on;
    dat1     = squeeze(datTmp(:,iEps,iBatchVals,clus2plot,2));
    mu      = nanmean(dat1,1);
    sm      = nanstd(dat1)./sqrt(size(dat1,1));
    ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
    plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
    scatter(barpos,mu,50,colgrey,'filled','d');
    xticklabels(xTickLabs);
    xlim([barpos(1)-.5, barpos(end)+.5]);
    ylim([-.5,1.5]);
    title('Left half of Trapezoid box')
    xlabel('Number of Clusters');
    ylabel('Grid Score');
    set(gca,'FontSize',fontSiz,'fontname','Arial')
    
    fname = [figsDir sprintf('/gridness_%s_L_univarScatters_testSet_nClus%d-%d_eps%d_batchSiz%d_%s_%s',dat,clus2run(clus2plot(1)),clus2run(clus2plot(end)),epsMuVals(iEps)*1000,batchSizeVals(iBatchVals),clusPosAct,gridMsrType)];
    
    if savePlots
        set(gcf,'Renderer','painters');
        print(gcf,'-depsc2',fname)
        saveas(gcf,fname,'png');
    end
    
    
    figure; hold on;
    dat1     = squeeze(datTmp(:,iEps,iBatchVals,clus2plot,3));
    mu      = nanmean(dat1,1);
    sm      = nanstd(dat1)./sqrt(size(dat1,1));
    ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
    plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
    scatter(barpos,mu,50,colgrey,'filled','d');
    xticklabels(xTickLabs);
    xlim([barpos(1)-.5, barpos(end)+.5]);
    ylim([-.5,1.5]);
    title('Right half of Trapezoid box')
    xlabel('Number of Clusters');
    ylabel('Grid Score');
    set(gca,'FontSize',fontSiz,'fontname','Arial')
    
    fname = [figsDir sprintf('/gridness_%s_R_univarScatters_testSet_nClus%d-%d_eps%d_batchSiz%d_%s_%s',dat,clus2run(clus2plot(1)),clus2run(clus2plot(end)),epsMuVals(iEps)*1000,batchSizeVals(iBatchVals),clusPosAct,gridMsrType)];
    
    if savePlots
        set(gcf,'Renderer','painters');
        print(gcf,'-depsc2',fname)
        saveas(gcf,fname,'png');
    end
    
    figure; hold on;
    dat1     = squeeze(datTmp(:,iEps,iBatchVals,clus2plot,2)-datTmp(:,iEps,iBatchVals,clus2plot,3));
    mu      = nanmean(dat1,1);
    sm      = nanstd(dat1)./sqrt(size(dat1,1));
    ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
    plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
    scatter(barpos,mu,50,colgrey,'filled','d');
    xticklabels(xTickLabs);
    xlim([barpos(1)-.5, barpos(end)+.5]);
    ylim([-1.5,2]);
    xlabel('Number of Clusters');
    ylabel('Grid Score');
    set(gca,'FontSize',fontSiz,'fontname','Arial')
    
    title('Left-right half of Trapezoid box')
    
    fname = [figsDir sprintf('/gridness_%s_L-R_univarScatters_testSet_nClus%d-%d_eps%d_batchSiz%d_%s_%s',dat,clus2run(clus2plot(1)),clus2run(clus2plot(end)),epsMuVals(iEps)*1000,batchSizeVals(iBatchVals),clusPosAct,gridMsrType)];
    
    if savePlots
        set(gcf,'Renderer','painters');
        print(gcf,'-depsc2',fname)
        saveas(gcf,fname,'png');
        close all
    end
    
end

% mean & bootstrap 95 CIs
if computeCIs
if ~strcmp(dat,'trapzKfrmSq1')
    nBoot = nIter;
    clear ciTest ciTmp
    for iClus=1:length(clus2plot)
        ciTmp = bootci(nBoot,@nanmean,dat1(:,iClus));
        ciTest(iClus,:) = [nanmean(dat1(:,iClus),1); ciTmp]; % mean & CIs
    end
%     ci = bootci(length(clus2run),@nanmean,nanmean(dat1,1)); %CI over mean of nClus conds - prob not as gd
    ci = bootci(numel(dat1),@nanmean,reshape(dat1,1,numel(dat1))); %CI over ALL runs 
    fprintf('%d to %d clusters: mean=%0.4f; CI=[%0.4f,%0.4f]\n',clus2run(clus2plot(1)),clus2run(clus2plot(end)),nanmean(reshape(dat1,1,numel(dat1))),ci(1),ci(2))
else
    %trapz 
    dat1     = squeeze(datTmp(:,iEps,iBatchVals,clus2plot,1));
    ci = bootci(numel(dat1),@nanmean,reshape(dat1,1,numel(dat1))); %CI over ALL runs 
    % individual nClus cond CIs, whole trapz
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
    fprintf('trapz, %d to %d clusters: mean=%0.4f; CI=[%0.4f,%0.4f];  %d sig > 0, out of %d sig conds. %0.2f percent\n',clus2plot(1)+2,clus2plot(end)+2,nanmean(reshape(dat1,1,numel(dat1))),ci(1),ci(2),nnz(ciClusSig==1),nnz(ciClusSig),(nnz(ciClusSig==1)./nnz(ciClusSig))*100);
    dat1     = squeeze(datTmp(:,iEps,iBatchVals,clus2plot,2));
    ci = bootci(numel(dat1),@nanmean,reshape(dat1,1,numel(dat1))); %CI over ALL runs 
    fprintf('trapzL %d to %d clusters: mean=%0.4f; CI=[%0.4f,%0.4f]\n',clus2plot(1)+2,clus2plot(end)+2,nanmean(reshape(dat1,1,numel(dat1))),ci(1),ci(2));
    dat1     = squeeze(datTmp(:,iEps,iBatchVals,clus2plot,3));
    ci = bootci(numel(dat1),@nanmean,reshape(dat1,1,numel(dat1))); %CI over ALL runs 
    fprintf('trapzR %d to %d clusters: mean=%0.4f; CI=[%0.4f,%0.4f]\n',clus2plot(1)+2,clus2plot(end)+2,nanmean(reshape(dat1,1,numel(dat1))),ci(1),ci(2));
    % individual nClus cond CIs, L-R
    dat1     = squeeze(datTmp(:,iEps,iBatchVals,clus2plot,2))-squeeze(datTmp(:,iEps,iBatchVals,clus2plot,3));
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
    fprintf('trapzL-R, %d to %d clusters: mean=%0.4f; CI=[%0.4f,%0.4f];  %d sig > 0, out of %d sig conds. %0.2f percent\n',clus2run(clus2plot(1)),clus2run(clus2plot(end)),nanmean(reshape(dat1,1,numel(dat1))),ci(1),ci(2),nnz(ciClusSig==1),nnz(ciClusSig),(nnz(ciClusSig==1)./nnz(ciClusSig))*100);
end
end
%%

if plotFewClus
%plot 3:5
clus2plot=(3:5)-2;
clus2plot=(3:9)-2;
%fig specs
xTickLabs = num2cell(clus2run(clus2plot));
fontSiz=15;
figure; hold on;
    for iEps = 1:length(epsMuVals)
        %     subplot(2,3,iEps);
        dat1     = squeeze(datTmp(:,iEps,iBatchVals,clus2plot,1));
        barpos  = .25:.5:.5*size(dat1,2);
        colors  = distinguishable_colors(size(dat1,2));
        colgrey = [.5, .5, .5];
        mu      = mean(dat1,1);
        sm      = std(dat1)./sqrt(size(dat1,1));
        ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
        plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
        errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
        scatter(barpos,mu,50,colors,'filled','d');
        xticklabels(xTickLabs);
        xlim([barpos(1)-.5, barpos(end)+.5]);
        %         ylim([0,1]);
        
        if strcmp(gridMsrType,'a')
            ylim([-.45,1.4]);
        elseif strcmp(gridMsrType,'w')
            ylim([-1.25,1.4]);
        end
        xlabel('Number of Clusters');
        ylabel('Grid Score');
        
        title(sprintf('%s, %s - eps=%d, batchSize=%d',dat, gridMeasure,epsMuVals(iEps)*1000,batchSizeVals(iBatchVals)))
        set(gca,'FontSize',fontSiz,'fontname','Arial')
    end

    fname = [figsDir sprintf('/gridness_%s_univarScatters_nClus%d-%d_eps%d_batchSiz%d_%s',dat,clus2run(clus2plot(1)),clus2run(clus2plot(end)),epsMuVals(iEps)*1000,batchSizeVals(iBatchVals),gridMsrType)];
if savePlots
   set(gcf,'Renderer','painters');
   print(gcf,'-depsc2',fname)
   saveas(gcf,fname,'png');
    close all
end
end
%% Making figs: density plot examples

savePlots = 0;

doPlot=0; %do plot when computing gridness

fontSiz=15;

clusPosAct = 'actNorm'; %'clus' or 'actNorm'

gridMsrType = 'a';

% clus2plot = [7,10,12,25]-2;%[5,8,13,18];
% clus2plot = 10-2;
%
clus2plot = [10,12,18,25]-9;
% clus2plot = [7,10,12,18,25]-9;
clus2plot = ([10,12,18,20,23,25,28])-9;
% clus2plot = (3:5)-2;

clus2plot=20-9;

myColorMap = parula;
myColorMap(end,:) = 1;

for iClus = clus2plot%:length(clus2run)
for iBvals = 1:length(batchSizeVals)
    for iterI = 1:10%nIter
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
                zlims=[0,.085];
                if strcmp(dat(1:2),'ci') || strcmp(dat(1:2),'tr') 
                zlims=[0,.085];
                elseif strcmp(dat(1:2),'sq')
                zlims=[0,.085];
                end
        end

        aCorrMap=ndautoCORR(densityPlotCentresSm); %autocorrelogram        
        [g,gdataA] = gridSCORE(aCorrMap,'allen',doPlot);
%         [g,gdataW] = gridSCORE(aCorrMap,'wills',doPlot);
                
%       %make nans into white
        densityPlotTmp = densityPlotCentresSm;
        densityPlotTmp(isnan(densityPlotCentresSm))=1.1;
        aCorrMap(isnan(aCorrMap))=1.1;
        
        figure; imagesc(densityPlotTmp,zlims); colormap(gcf,myColorMap);
%         figure; imagesc(densityPlotTmp); colorbar;
        xticks([]); xticklabels({''}); yticks([]); yticklabels({''});
        fname = [figsDir sprintf('/densityPlot_testSet_%s_nClus%d_eps%d_batchSiz%d_%s_%s_iter%d_clusters',dat,iClus+9,epsMuVals(iEps)*1000,batchSizeVals(iBvals),clusPosAct,gridMsrType,iterI)];
        if annEps
            fname = [fname '_annEps'];
        end
        if savePlots
            set(gcf,'Renderer','painters');
            print(gcf,'-depsc2',fname)
            saveas(gcf,fname,'png');
           close all
        end
        figure; h=imagesc(aCorrMap,[-1.1 1.1]); colormap(gcf,myColorMap);
%         figure; imagesc(aCorrMap); colorbar
        title(sprintf('Grid Score %.2f',g));
        set(gca,'FontSize',fontSiz,'fontname','Arial')
        xticks([]); xticklabels({''}); yticks([]); yticklabels({''});
        fname = [figsDir sprintf('/densityPlot_testSet_%s_nClus%d_eps%d_batchSiz%d_%s_%s_iter%d_autoCorr',dat,iClus+9,epsMuVals(iEps)*1000,batchSizeVals(iBvals),clusPosAct,gridMsrType,iterI)];
        if annEps
            fname = [fname '_annEps'];
        end
        if savePlots
            set(gcf,'Renderer','painters');
            print(gcf,'-depsc2',fname)
            saveas(gcf,fname,'png');
           close all
        end
        
        if strcmp(dat(1:4),'trap')
            h=50; hLeft = 14; hRight = 20;
            aCorrMap = ndautoCORR(densityPlotCentresSm(:,1:hLeft));
            [g,gdataA] = gridSCORE(aCorrMap,'allen',doPlot);
            %         [g,gdataW] = gridSCORE(aCorrMap,'wills',doPlot);
            
            aCorrMap(isnan(aCorrMap))=.5;
            figure; imagesc(aCorrMap,[-1.01 1.01]); colormap(myColorMap);
            title(sprintf('Grid Score %.2f',g));
            set(gca,'FontSize',fontSiz,'fontname','Arial')
            xticks([]); xticklabels({''}); yticks([]); yticklabels({''});
            fname = [figsDir sprintf('/densityPlot_testSet_%s_L_nClus%d_eps%d_batchSiz%d_%s_%s_iter%d_autoCorr',dat,iClus+9,epsMuVals(iEps)*1000,batchSizeVals(iBvals),clusPosAct,gridMsrType,iterI)];
            if annEps
               fname = [fname '_annEps']; 
            end
            if savePlots
                set(gcf,'Renderer','painters');
                print(gcf,'-depsc2',fname)
                saveas(gcf,fname,'png');
                close all
            end
            
            aCorrMap = ndautoCORR(densityPlotCentresSm(:,h-hRight+1:end));
            [g,gdataA] = gridSCORE(aCorrMap,'allen',doPlot);
            %         [g,gdataW] = gridSCORE(aCorrMap,'wills',doPlot);
            aCorrMap(isnan(aCorrMap))=.5;
            figure; imagesc(aCorrMap,[-1.01 1.01]);
            title(sprintf('Grid Score %.2f',g));
            set(gca,'FontSize',fontSiz,'fontname','Arial')
            xticks([]); xticklabels({''}); yticks([]); yticklabels({''});
            fname = [figsDir sprintf('/densityPlot_testSet_%s_R_nClus%d_eps%d_batchSiz%d_%s_%s_iter%d_autoCorr',dat,iClus+9,epsMuVals(iEps)*1000,batchSizeVals(iBvals),clusPosAct,gridMsrType,iterI)];
            if annEps
               fname = [fname '_annEps']; 
            end
            if savePlots
                set(gcf,'Renderer','painters');
                print(gcf,'-depsc2',fname)
                saveas(gcf,fname,'png');
                close all
            end
        end
    end
end
end

%% Thresholds from perm stats
figsDir = [wd '/grid_figs'];

savePlots=0;

clusPosAct = 'actNorm'; %'act' or 'actNorm'

gridMsrType = 'a'; % 'a' or 'w' for allen or willis method - a preferred

gridMeasure = 'grid';

switch clusPosAct
    case 'act'
        switch gridMsrType
            case 'a'
                permThresh = gA_act_permPrc;
            case 'w'
                permThresh = gW_act_permPrc;
        end
    case 'actNorm'
        switch gridMsrType
            case 'a'
                permThresh = gA_actNorm_permPrc;
            case 'w'
                permThresh = gW_actNorm_permPrc;
        end
end

% clus2plot=([3:30])-2;
% clus2plot=([6:26])-2;
% clus2plot=([10:26])-2;
% clus2plot=([10:30])-2;
clus2plot=1:length(clus2run);

iBatchVals=1; %'medium' one

%fig specs
xTickLabs = num2cell(clus2run(clus2plot));
fontSiz=15;

% gA_act_permPrc(:,1,1,:)

dSize=400;

figure; hold on;
    for iEps = 1:length(epsMuVals)
        %     subplot(2,3,iEps);
        dat1     = squeeze(permThresh(:,iEps,iBatchVals,clus2plot,1));
        barpos  = .25:.5:.5*size(dat1,2);
        colors  = distinguishable_colors(size(dat1,2));
        colgrey = [.5, .5, .5];
        mu      = mean(dat1,1);
%         sm      = std(dat1)./sqrt(size(dat1,1));
%         ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
        plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
%         scatter(barpos,mean(dat1),dSize,colors,'d'); %mean
        scatter(barpos,mean(dat1),50,colgrey,'filled','d');

%         scatter(barpos,max(dat1),dSize,colors,'d');  %max
        scatter(barpos,prctile(dat1,100),dSize,colors,'x');  %

%         prctile(dat1,[2.5, 97.5])
        xticklabels(xTickLabs);
        xlim([barpos(1)-.5, barpos(end)+.5]);
%         if strcmp(gridMsrType,'a')
%             ylim([-.45,1.4]);
%         elseif strcmp(gridMsrType,'w')
%             ylim([-1.25,1.4]);
%         end
        xlabel('Number of Clusters');
        ylabel('Grid Score Threshold (95% of permuted distribution)');
        
        
%         title(sprintf('%s, %s - eps=%d, batchSize=%d',dat, gridMeasure,epsMuVals(iEps)*1000,batchSizeVals(iBatchVals)))
        
        set(gca,'FontSize',fontSiz,'fontname','Arial')
    end

    fname = [figsDir sprintf('/gridness_%s_univarScatters_permThresh_nClus%d-%d_eps%d_batchSiz%d_%s_%s',dat,clus2run(clus2plot(1)),clus2run(clus2plot(end)),epsMuVals(iEps)*1000,batchSizeVals(iBatchVals),clusPosAct,gridMsrType)];
    if annEps
       fname = [fname '_annEps']; 
    end
if savePlots
   set(gcf,'Renderer','painters');
   print(gcf,'-depsc2',fname)
   saveas(gcf,fname,'png');
end

   


% squeeze(min(min(permThresh(:,1,1,:))))
% squeeze(mean(mean(permThresh(:,1,1,:))))
% squeeze(max(max(permThresh(:,1,1,:))))
    


%% trapz vs orig sq gridness / circ vs sq

savePlots=1;

clusPosAct = 'actNorm'; %'act' or 'actNorm'

gridMsrType = 'a'; % 'a' or 'w' for allen or willis method - a preferred

gridMeasure = 'grid';

computeCIs=0; 

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
                gridness    = gA_gAll_actNormSq-gA_gAll_actNorm;
            case 'w'
                gridness    = gW_gAll_actNormSq-gW_gAll_actNorm;
        end
end

%flip sign for circ, so plotting circ>sq gridness
if strcmp(dat(1:4),'circ')
    gridness = -gridness;
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
xTickLabs = num2cell(clus2run(clus2plot));
fontSiz=15;

figure; hold on;
    for iEps = 1:length(epsMuVals)
        dat1     = squeeze(datTmp(:,iEps,iBatchVals,clus2plot,1));
        barpos  = .25:.5:.5*size(dat1,2);
        colors  = distinguishable_colors(size(dat1,2));
        colgrey = [.6, .6, .6];
        mu      = mean(dat1,1);
%         sm      = std(dat1)./sqrt(size(dat1,1));
%         ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
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
        title('Square-Trapezoid')

        set(gca,'FontSize',fontSiz,'fontname','Arial')
    end
    fname = [figsDir sprintf('/gridness_%s_univarScatters_testSet_nClus%d-%d_eps%d_batchSiz%d_%s_%s','trapzKfrmSq1_Sq-Trapz',clus2run(clus2plot(1)),clus2run(clus2plot(end)),epsMuVals(iEps)*1000,batchSizeVals(iBatchVals),clusPosAct,gridMsrType)];
if savePlots
   set(gcf,'Renderer','painters');
   print(gcf,'-depsc2',fname)
   saveas(gcf,fname,'png');
end

if computeCIs
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
    if strcmp(dat(1:4),'trap')
        fprintf('Square-Trapezoid box, %d to %d clusters: mean=%0.4f; CI=[%0.4f,%0.4f]; %d sig Square-Trapezoid, %d sig Trapezoid-Square , %0.2f percent Sq>T\n',clus2plot(1)+2,clus2plot(end)+2,mean(reshape(dat1,1,numel(dat1))),ci(1),ci(2),nnz(ciClusSig==1),nnz(ciClusSig==-1),(nnz(ciClusSig==1)./nnz(ciClusSig))*100);
    elseif strcmp(dat(1:4),'circ')
        fprintf('Circ-Square box, %d to %d clusters: mean=%0.4f; CI=[%0.4f,%0.4f]; %d sig Circ-Square, %d sig Square-Circ, %0.2f percent C>Sq\n',clus2plot(1)+2,clus2plot(end)+2,mean(reshape(dat1,1,numel(dat1))),ci(1),ci(2),nnz(ciClusSig==1),nnz(ciClusSig==-1),(nnz(ciClusSig==1)./nnz(ciClusSig))*100);
    end
end


