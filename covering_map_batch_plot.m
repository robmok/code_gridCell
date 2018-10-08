%% Plotting script 1: train set

clear all;

% Set working directory
% wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
cd(wd);

codeDir = [wd '/code_gridCell'];
saveDir = [wd '/data_gridCell'];
figsDir = [wd '/grid_figs'];
addpath(codeDir); addpath(saveDir);
addpath(genpath([codeDir '/gridSCORE_packed']));

%smoothing
gaussSmooth = 1; 
imageFilter = fspecial('gaussian',5,gaussSmooth); %this is default for imgaussfilt

% load
fixBatchSize = 1; %fixed batch size or depend on nClus (for fname)

dat='circ';
% dat='square';

% joined trials
jointTrls=1;

actOverTime = 1;
nSet        = 21;

%1000 iters
nIter=1000;

if nIter == 1000
   actOverTime = 0;
   nSet        = 1;
end

%new - annealed learning rate
clus2run      = 10:30;
epsMuVals     = 0.25;
nTrials       = 1000000;
batchSizeVals = 200;
annEps        =1;

rHex = 0; %if inspect raw 60deg corr values, rather than grid score

%trapzKfrmSq 
dat           = 'trapzKfrmSq1';
nTrials       = 500000;
batchSizeVals = 400;
clus2run      = 10:30;
nIter         = 1000;

% annEps=0;
% epsMuVals = 0.025;

%load loop
for iClus2run = 1:length(clus2run) 
    nClus = clus2run(iClus2run);
    for iEps = 1:length(epsMuVals) 
        epsMuOrig=epsMuVals(iEps);
        epsMuOrig1000=epsMuOrig*1000;
        for iBvals = 1:length(batchSizeVals)
            batchSize = batchSizeVals(iBvals);
            fprintf('Loading %s, nClus=%d, epsMu=%d, batchSize=%d\n',dat,nClus,epsMuOrig1000,batchSize)
            fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_jointTrls_stepSiz',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat)];
            if ~actOverTime
               fname = [fname '_noActOverTime']; 
            end
            %finish with directory and * for date/time
            if ~annEps
                fname = [saveDir, fname '*']; %finish with directory and * for date/time
            else
                fname = [saveDir, fname '*annEps*']; %new position of annEps - works now - above no need
            end
            
            %edit if want to load more than one file per sim, merge
            f = dir(fname); filesToLoad = cell(1,length(f));            
            iF=1;
            filesToLoad{iF} = f(iF).name;
            load(f(iF).name);
            
            for iterI = 1:nIter
                for iSet=1:nSet
                    densityPlotAll(:,:,iSet,iterI,iEps,iBvals,iClus2run) = nanconv(densityPlot(:,:,iSet,iterI),imageFilter, 'nanout'); % smooths ignoring nans
                end
            end
            iSet = 1:nSet; %added this and to ga_gAll (and actNorm) - when iSet = 1
            
            %organise gridness values (allen vs willis method)
            gA_gAll(iSet,:,:,iEps,iBvals,iClus2run,:)   = gA(iSet,:,1,:);
            gA_oAll(iSet,:,:,iEps,iBvals,iClus2run,:)   = gA(iSet,:,2,:);
            gA_radAll(iSet,:,:,iEps,iBvals,iClus2run,:) = gA(iSet,:,3,:);
            gA_wavAll(iSet,:,:,iEps,iBvals,iClus2run,:) = gA(iSet,:,4,:);
            gW_gAll(iSet,:,:,iEps,iBvals,iClus2run,:)   = gW(iSet,:,1,:);
            gW_oAll(iSet,:,:,iEps,iBvals,iClus2run,:)   = gW(iSet,:,2,:);
            gW_radAll(iSet,:,:,iEps,iBvals,iClus2run,:) = gW(iSet,:,3,:);
            gW_wavAll(iSet,:,:,iEps,iBvals,iClus2run,:) = gW(iSet,:,4,:);
            
                        
            if exist('rSeed','var')
                rSeedAll(:,iEps,iBvals,iClus2run) = rSeed;
            end
            
            if exist('gA_actNorm','var') %if have gridness on gauss act
                for iterI = 1:nIter
                    for iSet=1:nSet
                        densityPlotActNormAll(:,:,iSet,iterI,iEps,iBvals,iClus2run) = nanconv(densityPlotActNorm(:,:,iSet,iterI),imageFilter, 'nanout'); %new smooths ignoring nans
                    end
                end
                
                iSet = 1:nSet;
%                 
                gA_gAll_actNorm(iSet,:,:,iEps,iBvals,iClus2run,:)   = gA_actNorm(iSet,:,1,:);
                if rHex % note changes to iSet, for when iSet=1
%                     gA_gAll_actNorm(:,:,iEps,iBvals,iClus2run,:)   = (gA_actNorm(:,:,6,:)+gA_actNorm(:,:,8,:))./2; %rHex 60/120
                    gA_gAll_actNorm(iSet,:,:,iEps,iBvals,iClus2run,:)   = (gA_actNorm(iSet,:,5,:)+gA_actNorm(iSet,:,7,:)+(gA_actNorm(iSet,:,9,:)))./3; %rHex 30,90,150
%                     gA_gAll_actNorm(:,:,iEps,iBvals,iClus2run,:)   = (gA_actNorm(:,:,8,:)); %rHex 
                end
                gA_oAll_actNorm(iSet,:,:,iEps,iBvals,iClus2run,:)   = gA_actNorm(iSet,:,2,:);
                gA_radAll_actNorm(iSet,:,:,iEps,iBvals,iClus2run,:) = gA_actNorm(iSet,:,3,:);
                gA_wavAll_actNorm(iSet,:,:,iEps,iBvals,iClus2run,:) = gA_actNorm(iSet,:,4,:);
                gW_gAll_actNorm(iSet,:,:,iEps,iBvals,iClus2run,:)   = gW_actNorm(iSet,:,1,:);
                gW_oAll_actNorm(iSet,:,:,iEps,iBvals,iClus2run,:)   = gW_actNorm(iSet,:,2,:);
                gW_radAll_actNorm(iSet,:,:,iEps,iBvals,iClus2run,:) = gW_actNorm(iSet,:,3,:);
                gW_wavAll_actNorm(iSet,:,:,iEps,iBvals,iClus2run,:) = gW_actNorm(iSet,:,4,:);

            end            
        end
    end
end
%% plot univar scatters

clusPosAct = 'clus'; %'clus' or 'actNorm'

gridMsrType = 'a'; % 'a' or 'w' for allen or willis method - a preferred

gridMeasure = 'grid';

plotOverSets = 0;

switch clusPosAct
case 'clus'
%     switch gridMsrType
%         case 'a'
            gridness    = gA_gAll;
            orientation = gA_oAll;
            rad         = gA_radAll;
            wav         = gA_wavAll;
%         case 'w'
%             gridness    = gW_gAll;
%             orientation = gW_oAll;
%             rad         = gW_radAll;
%             wav         = gW_wavAll;
%     end
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

iSet=size(datTmp,1);

iEps=1;

%plot univar scatters - over clusters (e.g. one batch size)
if size(gridness,4)==1 %only 1 batchSize
    figure; hold on;
    dat1     = squeeze(datTmp(iSet,:,:,iEps,:,:,1));
    barpos  = .25:.5:.5*size(dat1,2);
    colors  = distinguishable_colors(size(dat1,2));
    colgrey = [.5, .5, .5];
    mu      = mean(dat1,1);
    sm      = std(dat1)./sqrt(size(dat1,1));
    ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
    plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
    errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
    scatter(barpos,mu,100,colors,'d','filled');
    xlim([barpos(1)-.5, barpos(end)+.5]);
    ylim([-.5,1.25]);
    title(sprintf('%s - eps=%d',gridMeasure,epsMuVals(iEps)*1000))
    
    if strcmp(dat(1:4),'trap')
        figure; hold on;
        dat1    = squeeze(datTmp(iSet,:,:,:,:,:,2));
        mu      = nanmean(dat1,1);
        sm      = nanstd(dat1)./sqrt(size(dat1,1));
        ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
        plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
        errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
        scatter(barpos,mu,100,colors,'d','filled');
        xlim([barpos(1)-.5, barpos(end)+.5]);
        ylim([-1.5,1.5]);
        title(sprintf('Left half of box %s - eps=%d',gridMeasure,epsMuVals(iEps)*1000))
        
        figure; hold on;
        dat1    = squeeze(datTmp(iSet,:,:,:,:,:,3));
        mu      = nanmean(dat1,1);
        sm      = nanstd(dat1)./sqrt(size(dat1,1));
        ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
        plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
        errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
        scatter(barpos,mu,100,colors,'d','filled');
        xlim([barpos(1)-.5, barpos(end)+.5]);
        ylim([-1.5,1.5]);
        title(sprintf('Right half of box %s - eps=%d',gridMeasure,epsMuVals(iEps)*1000))
        
        figure; hold on;
        dat1    = squeeze(datTmp(iSet,:,:,:,:,:,2))-squeeze(datTmp(iSet,:,:,:,:,:,3));
        mu      = nanmean(dat1,1);
        sm      = nanstd(dat1)./sqrt(size(dat1,1));
        ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
        plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
        errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
        scatter(barpos,mu,100,colors,'d','filled');
        xlim([barpos(1)-.5, barpos(end)+.5]);
        ylim([-1.5,1.5]);
        title(sprintf('Left-right half of box %s - eps=%d',gridMeasure,epsMuVals(iEps)*1000))
    end
end


% %plot univar scatters - comparing batchSizeVals, with clusters in separate plots
if size(gridness,4)>1 % more than 1 batchSize
figure; hold on;
for iClus2Run = 1:length(clus2run)
    % comparing batch vals, with cluster nums in subplots; can edit for eps
%     figure; hold on;
    subplot(ceil(length(clus2run)/2),2,iClus2Run);
    for iEps = 1:length(epsMuVals)
        %     subplot(2,3,iEps);
        dat1     = squeeze(datTmp(iSet,:,:,iEps,:,iClus2Run,1));
        barpos  = .25:.5:.5*size(dat1,2);
        colors  = distinguishable_colors(size(dat1,2));
        colgrey = [.5, .5, .5];
        mu      = mean(dat1,1);
        sm      = std(dat1)./sqrt(size(dat1,1));
        ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
        plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
        errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
        scatter(barpos,mu,750,colors,'x');
        xlim([barpos(1)-.5, barpos(end)+.5]);
%         ylim([0,1]);
        ylim([-.5,1.5]);
        title(sprintf('%s - eps=%d, nClus=%d',gridMeasure,epsMuVals(iEps)*1000,clus2run(iClus2Run)))
    end
end
end

%plot over sets
if plotOverSets
    for iClus = 1:length(clus2run)
    figure; hold on;
    for iBatchVals = 1:length(batchSizeVals)        
%         subplot(ceil(length(batchSizeVals)/2),2,iBatchVals); hold on;
        dat1     = squeeze(datTmp(:,:,iEps,iBatchVals,iClus,:))';
        barpos  = .25:.5:.5*size(dat1,2);
        colors  = distinguishable_colors(size(dat1,2));
        colgrey = [.5, .5, .5];
        mu      = nanmean(dat1,1);
        sm      = nanstd(dat1)./sqrt(size(dat1,1));
        ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
        plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
        errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
        scatter(barpos,mu,100,colors,'d','filled');
        xlim([barpos(1)-.5, barpos(end)+.5]);
        ylim([-.5,1.25]);
%         title(sprintf('%s - nClus=%d, eps=%d, batchSize=%d',gridMeasure,clus2run(iClus),epsMuVals(iEps)*1000,batchSizeVals(iBatchVals)))
        title(sprintf('nClus=%d, batchSize=%d',clus2run(iClus),batchSizeVals(iBatchVals)))
    end 
    end
end

figure; hold on;
for iClus=1:length(clus2run)
    subplot(ceil(length(clus2run)/2),2,iClus); hold on;
    if strcmp(dat(1:4),'trap')
        dat1     = squeeze(datTmp(:,:,:,:,iClus,2))'-squeeze(datTmp(:,:,:,:,iClus,3))';
        barpos  = .25:.5:.5*size(dat1,2);
        colors  = distinguishable_colors(size(dat1,2));
        mu      = mean(dat1,1);
        sm      = std(dat1)./sqrt(size(dat1,1));
        ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
        plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
        errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
        scatter(barpos,mu,100,colors,'d','filled');
        xlim([barpos(1)-.5, barpos(end)+.5]);
        ylim([-1.25,1.25]);
        title(sprintf('Left-Right half of box %s - eps=%d',gridMeasure,epsMuVals(iEps)*1000))
    end
end
%% check gridness measures

doPlot=0;

iSet=6;
iEps=1;
gaussSmooth=1;
nIter=200;

for iClus2run = 1:length(clus2run)
for iBvals = 1:length(batchSizeVals)
    fprintf(sprintf('Running clus %d batchSizeVals %d\n',clus2run(iClus2run),batchSizeVals(iBvals)));
    for iterI = 1:nIter
        densityPlotCentresSm = nanconv(densityPlotAll(:,:,iSet,iterI,iEps,iBvals,iClus2run),imageFilter, 'nanout');
        if doPlot
            figure; hold on;
            subplot(1,2,1)
            imagesc(densityPlotCentresSm);
            subplot(1,2,2);
        end
        
        aCorrMap=ndautoCORR(densityPlotCentresSm); %autocorrelogram
        [g,gdataA] = gridSCORE(aCorrMap,'allen',doPlot);
%         [g,gdataW] = gridSCORE(aCorrMap,'wills',doPlot);
        
        peaksA(iClus2run,iBvals,iterI,1)=length(gdataA.near_peaks);
%         peaksW(iClus2run,iBvals,iterI,1)=length(gdataW.near_peaks);
        
        if strcmp(dat,'trapz')
            aCorrMap = ndautoCORR(densityPlotCentresSm(:,1:size(densityPlotAll,1)/2));
            [g,gdataA] = gridSCORE(aCorrMap,'allen',doPlot);
            %         [g,gdataW] = gridSCORE(aCorrMap,'wills',doPlot);
            peaksA(iClus2run,iBvals,iterI,2)=length(gdataA.near_peaks);
            gTmpA(iClus2run,iBvals,iterI,2) = gdataA.g_score;
            %             peaksW(iClus2run,iBvals,iterI,2)=length(gdataW.near_peaks);
            
            aCorrMap = ndautoCORR(densityPlotCentresSm(:,size(densityPlotAll,1)/2+1:end));
            [g,gdataA] = gridSCORE(aCorrMap,'allen',doPlot);
            %         [g,gdataW] = gridSCORE(aCorrMap,'wills',doPlot);
            peaksA(iClus2run,iBvals,iterI,3)=length(gdataA.near_peaks);
            gTmpA(iClus2run,iBvals,iterI,3) = gdataA.g_score;
            %             peaksW(iClus2run,iBvals,iterI,3)=length(gdataW.near_peaks);
        end
    end
end
end

nnz(peaksA~=6)
nnz(peaksW~=6)
% figure; scatter(reshape(peaksA(:,:,:,2),1,numel(peaksA(:,:,:,2))), reshape(gTmpA(:,:,:,2),1,numel(gTmpA(:,:,:,2)))); ylim([-.1, 1]);
% figure; scatter(reshape(peaksA(:,:,:,3),1,numel(peaksA(:,:,:,3))), reshape(gTmpA(:,:,:,3),1,numel(gTmpA(:,:,:,3)))); ylim([-.1, 1]);
%% density plots
close all

iSet=21;
iEps=1;
gaussSmooth=1;

if strcmp(dat(1:4),'trap')
subPlt = [2,2];
else
    subPlt = [1,2];
end

%set
iClus2run = 20-9;
iClus2run = 6;
iBvals    = 1;

iters2plot = 1;

fprintf(sprintf('clus %d batchSizeVals %d\n',clus2run(iClus2run),batchSizeVals(iBvals)));
for iterI = iters2plot

    densityPlotCentresSm = densityPlotAll(:,:,iSet,iterI,iEps,iBvals,iClus2run);

    figure; hold on;
    subplot(subPlt(1),subPlt(2),1)
    imagesc(densityPlotCentresSm);
    subplot(subPlt(1),subPlt(2),2)
    
    aCorrMap=ndautoCORR(densityPlotCentresSm); %autocorrelogram
    [g,gdataA] = gridSCORE(aCorrMap,'allen',1);
    %         [g,gdataW] = gridSCORE(aCorrMap,'wills',1);
    
%     peaksA(iClus2run,iBvals,iterI,1)=length(gdataA.near_peaks);
    
    if strcmp(dat(1:4),'trap')
        subplot(subPlt(1),subPlt(2),3);
%         aCorrMap = ndautoCORR(densityPlotCentresSm(:,1:size(densityPlotAll,1)/2));
        aCorrMap = ndautoCORR(densityPlotCentresSm(19:32,1:20));
        [g,gdataA] = gridSCORE(aCorrMap,'allen',1);
        subplot(subPlt(1),subPlt(2),4);
%         aCorrMap = ndautoCORR(densityPlotCentresSm(:,size(densityPlotAll,1)/2+1:end));
        aCorrMap = ndautoCORR(densityPlotCentresSm(19:32,22:end));
        [g,gdataA] = gridSCORE(aCorrMap,'allen',1);
        
    end
    
end
%% Making figs - univar scatters 1 - training set

savePlots=0;

clusPosAct = 'actNorm'; %'clus' or 'actNorm'

gridMsrType = 'a'; % 'a' or 'w' for allen or willis method - a preferred

gridMeasure = 'grid';

switch clusPosAct
case 'clus'
    switch gridMsrType
        case 'a'
            gridness    = gA_gAll;
            orientation = gA_oAll;
            rad         = gA_radAll;
            wav         = gA_wavAll;
        case 'w'
            gridness    = gW_gAll;
            orientation = gW_oAll;
            rad         = gW_radAll;
            wav         = gW_wavAll;
    end
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


iSet=size(datTmp,1); %last or second last (add -1)
iEps=1;

clus2plot = 1:length(clus2run);
% clus2plot=(6:30)-2;
% clus2plot=(10:30)-2;
% clus2plot=(10:26)-2;

iBatchVals=1; 

%fig specs
xTickLabs = num2cell(clus2run(clus2plot));
fontSiz=15;

figure; hold on;
    for iEps = 1:length(epsMuVals)
        %     subplot(2,3,iEps);
        dat1     = squeeze(datTmp(iSet,:,:,iEps,iBatchVals,clus2plot,1));
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
            ylim([-.45,1.4]);
        elseif strcmp(gridMsrType,'w')
            ylim([-1.25,1.4]);
        end
        xlabel('Number of Clusters');
        ylabel('Grid Score');        
%         title(sprintf('%s, %s - eps=%d, batchSize=%d',dat, gridMeasure,epsMuVals(iEps)*1000,batchSizeVals(iBatchVals)))
        if strcmp(dat(1:2),'ci')
        title('Circular box')
        elseif strcmp(dat(1:2),'sq')
        title('Square box')
        end        
        set(gca,'FontSize',fontSiz,'fontname','Arial')
    end
    fname = [figsDir sprintf('/gridness_%s_univarScatters_trainSetEnd_nClus%d-%d_eps%d_nIter%d_batchSiz%d_%s_%s',dat,clus2run(clus2plot(1)),clus2run(clus2plot(end)),epsMuVals(iEps)*1000,nIter,batchSizeVals(iBatchVals),clusPosAct,gridMsrType)];
if savePlots
   set(gcf,'Renderer','painters');
   print(gcf,'-depsc2',fname)
   saveas(gcf,fname,'png');
end

% mean & bootstrap 95 CIs
nBoot = nIter;
clear ciTrain ciTmp
for iClus=1:length(clus2plot)
    ciTmp = bootci(nBoot,@mean,dat1(:,iClus));
    ciTrain(iClus,:) = [mean(dat1(:,iClus),1); ciTmp]; % mean & CIs
end

% ciTrain %print 


%overall mean
ci = bootci(numel(dat1),@nanmean,reshape(dat1,1,numel(dat1))); %CI over ALL runs
fprintf('%d to %d clusters: mean=%0.4f; CI=[%0.4f,%0.4f]\n',clus2run(clus2plot(1)),clus2run(clus2plot(end)),nanmean(reshape(dat1,1,numel(dat1))),ci(1),ci(2));

%%
   
%plot 3:5

clus2plot=(3:5)-2;

%fig specs
xTickLabs = num2cell(clus2run(clus2plot));
fontSiz=15;
figure; hold on;
    for iEps = 1:length(epsMuVals)
        %     subplot(2,3,iEps);
        dat1     = squeeze(datTmp(iSet,:,:,iEps,iBatchVals,clus2plot,1));
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

    fname = [figsDir sprintf('/gridness_%s_univarScatters_nClus%d-%d_eps%d_nIter%d_batchSiz%d_%s',dat,clus2run(clus2plot(1)),clus2run(clus2plot(end)),epsMuVals(iEps)*1000,nIter,batchSizeVals(iBatchVals),gridMsrType)];
if savePlots
   set(gcf,'Renderer','painters');
   print(gcf,'-depsc2',fname)
   saveas(gcf,fname,'png');
end

% for i=1:length(clus2run),
% xx(i)=nnz(dat1(:,i)>.37)/200;
% end
%mean(xx) %prop grid cells, averaged across clusters
    

%% Making figs: plot over sets

savePlots = 0;

plotSubPlots = 0;

fontSiz = 25;

iBatchVals=1;

nTimePts=size(datTmp,1)-1;
colors  = distinguishable_colors(length(clus2run));

if plotSubPlots
%     clus2plot = (6:26)-2;
    clus2plot = 1:length(clus2run);
    figure; hold on;
    for iClus = clus2plot
        %     subplot(2,ceil(length(clus2plot)/2),iClus); hold on;
%         subplot(3,7,iClus-clus2plot(1)+1); hold on;
        subplot(4,6,iClus-clus2plot(1)+1); hold on;
        %         subplot(ceil(length(clus2plot)/2),2,iClus); hold on;
        dat1     = squeeze(datTmp(1:nTimePts,:,iEps,iBatchVals,:,iClus))';
        barpos  = .25:.5:.5*size(dat1,2);
        colgrey = [.5, .5, .5];
        mu      = nanmean(dat1,1);
        sm      = nanstd(dat1)./sqrt(size(dat1,1));
        ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
        plotSpread(dat1,'xValues',barpos,'distributionColors',repmat(colors(iClus,:),nTimePts,1));
        %     errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
        scatter(barpos,mu,10,colgrey,'d','filled');
        %         scatter(barpos,mu,100,repmat(colors(iClus,:),nTimePts,1),'d','filled');
        xlim([barpos(1)-.5, barpos(end)+.5]);
        ylim([-.5,1.25]);
        xticks([]); xticklabels({''});
        %         title(sprintf('nClus=%d, batchSize=%d',clus2plot(iClus),batchSizeVals(iBatchVals)))
%         title(sprintf('%d clusters',clus2plot(iClus-clus2plot(1)+1)+2));
        title(sprintf('%d clusters',clus2run(clus2plot(iClus))));
    end
    fname = [figsDir sprintf('/gridness_%s_univarScatters_overTime_nClus%d-%d_eps%d_batchSiz%d_%s',dat,clus2run(clus2plot(1)),clus2run(clus2plot(end)),epsMuVals(iEps)*1000,batchSizeVals(iBatchVals),gridMsrType)];
    if annEps
        fname = [fname '_annEps'];
    end
    if savePlots
        set(gcf,'Renderer','painters');
        print(gcf,'-depsc2',fname)
        saveas(gcf,fname,'png');
        close all
    end
else
    
    %subset of nClus conds
%     clus2plot = [7,10,12,18,25]-2;
%     clus2plot = [10,12,18,25]-9;
    clus2plot = 1:length(clus2run);
    
%     clus2plot = 18-2;
    for iClus = clus2plot
        figure;
        dat1     = squeeze(datTmp(1:nTimePts,:,:,iEps,iBatchVals,iClus,:))';
        barpos  = .25:.5:.5*size(dat1,2);
        colgrey = [.5, .5, .5];
        mu      = nanmean(dat1,1);
        sm      = nanstd(dat1)./sqrt(size(dat1,1));
        ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
        plotSpread(dat1,'xValues',barpos,'distributionColors',repmat(colors(iClus,:),nTimePts,1));
        %     errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
        scatter(barpos,mu,100,colgrey,'d','filled');
        %         scatter(barpos,mu,100,repmat(colors(iClus,:),nTimePts,1),'d','filled');
        xlim([barpos(1)-.5, barpos(end)+.5]);
        ylim([-.5,1.3]);
        xticks([]); xticklabels({''});
        %         title(sprintf('nClus=%d, batchSize=%d',clus2plot(iClus),batchSizeVals(iBatchVals)))
        title(sprintf('%d clusters',clus2run(iClus)));
        set(gca,'FontSize',fontSiz,'fontname','Arial')

        fname = [figsDir sprintf('/gridness_%s_univarScatters_overTime_singlePlot_nClus%d_eps%d_batchSiz%d_%s',dat,clus2run(iClus),epsMuVals(iEps)*1000,batchSizeVals(iBatchVals),gridMsrType)];
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

%% gridness over time stats

nTimePts = 20;


clus2plot = 1:length(clus2run);
% clus2plot = 4:length(clus2run); %skip 3:5
% clus2plot = 8:length(clus2run); %skip 3:9 %
% clus2plot = 8:length(clus2run)-4; %skip 3:9, 27:30 % turns out last few are positive!

% clear gt_b gt_p
% cnter=0;
% for iClus = clus2plot
%     dat1 = squeeze(datTmp(1:end-1,:,iEps,iBatchVals,:,iClus))';
%     [b,d,s]=glmfit(1:nTimePts,mean(dat1,1)');
%     cnter=cnter+1;
%     gt_b(cnter) = b(2);
%     gt_p(cnter) = s.p(2);
% end
% [h p c s] = ttest(gt_b); %sig 0.0363 for 3:30; 0.031 for 6:30; 0.0068 for 10:30; 0.0079 for 10:26
% 
% ci = bootci(length(gt_b),@mean,gt_b);
% fprintf('%d to %d clusters: mean=%0.4f; CI=[%0.4f,%0.4f]; p=%0.3f, %d sig betas > 0, out of %d sig betas. %0.2f percent\n',clus2plot(1)+2,clus2plot(end)+2,mean(gt_b),ci(1),ci(2),p,nnz(gt_b(gt_p<0.05)>0),nnz(gt_b(gt_p<0.05)),(nnz(gt_b(gt_p<0.05)>0)/numel(gt_b(gt_p<0.05)))*100);


%CORRECT:
% or - get each sim, glmfit, get that beta; use the mean of those betas. calc
% CIs for each nClus cond, report proportions; and over all nClus conds for the main stat

clear gt_b gt_p ciClus ciClusSig
cnter=0;
for iClus = clus2plot
    dat1 = squeeze(datTmp(1:end-1,:,iEps,iBatchVals,:,iClus))';
    cnter=cnter+1;
    for iterI=1:nIter
        [b,d,s]=glmfit(1:nTimePts,dat1(iterI,:)');
        gt_b(cnter,iterI) = b(2);
    end
    ciTmp = bootci(nIter,@mean,gt_b(cnter,:));
    ciClus(:,cnter) = [nanmean(gt_b(cnter,:)); ciTmp]; % mean & CIs    
    
    %check CIs sig
    posInd=ciClus(:,cnter)>0;
    negInd=ciClus(:,cnter)<0;
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
% ci = bootci(length(gt_b),@mean,mean(gt_b,2));
ci = bootci(numel(gt_b),@nanmean,reshape(gt_b,1,numel(gt_b))); %CI over ALL runs

fprintf('%d to %d clusters: mean=%0.4f; CI=[%0.4f,%0.4f]; %d sig betas > 0, out of %d sig betas. %0.2f percent\n',clus2run(clus2plot(1)),clus2run(clus2plot(end)),mean(reshape(gt_b,1,numel(gt_b))),ci(1),ci(2),nnz(ciClusSig==1),nnz(ciClusSig),(nnz(ciClusSig==1)./nnz(ciClusSig))*100);


%TO DO:

% % mean & bootstrap 95 CIs - taken from above. - edit here to compute
% means as well. note: ciClus is 2 x nClus above. below flipped so can
% print to copy and paste

% nBoot = nIter;
% clear ciTrain
% for iClus=clus2plot
%     ciTmp = bootci(nBoot,@mean,dat1(:,iClus));
%     ciTrain(iClus,:) = [mean(dat1(:,iClus),1); ciTmp]; % mean & CIs
% end
% ciTrain

%% Making figs: density plot examples

savePlots = 0;

doPlot=0; %do plot when computing gridness

fontSiz=15;

clusPosAct = 'actNorm'; %'clus' or 'actNorm'

gridMsrType = 'a';

% clus2plot = [5,8,13,18];%[7,10,12,25]
% clus2plot = ([10,11,12,18,20,23,25,28])-2;
clus2plot = ([10,11,12,18,20,23,25,28])-9;
clus2plot = (18)-9;

% iSet=1;
% iSet=10;
iSet=size(densityPlotActNormAll,3);

myColorMap = parula;
myColorMap(end,:) = 1;

%add zlims, aCorr lims
% zlims = [];



for iClus = clus2plot%:length(clus2run)
for iBvals = 1:length(batchSizeVals)
    for iterI = 1%:5%nIter
        
        densityPlotCentresSm = densityPlotActNormAll(:,:,iSet,iterI,iEps,iBvals,iClus);
%         densityPlotCentresSm(isnan(densityPlotCentresSm))=0;
        aCorrMap=ndautoCORR(densityPlotCentresSm); %autocorrelogram        
        
        densityPlotCentresSm(isnan(densityPlotCentresSm))=1.1;
        aCorrMap(isnan(aCorrMap))=1.1;
        
        
        [g,gdataA] = gridSCORE(aCorrMap,'allen',doPlot);
%         [g,gdataW] = gridSCORE(aCorrMap,'wills',doPlot);

% if strcmp(dat,'trapz')
%             aCorrMap = ndautoCORR(densityPlotCentresSm(:,1:size(densityPlotAll,1)/2));
%             [g,gdataA] = gridSCORE(aCorrMap,'allen',doPlot);
%             %         [g,gdataW] = gridSCORE(aCorrMap,'wills',doPlot);
%             aCorrMap = ndautoCORR(densityPlotCentresSm(:,size(densityPlotAll,1)/2+1:end));
%             [g,gdataA] = gridSCORE(aCorrMap,'allen',doPlot);
%             %         [g,gdataW] = gridSCORE(aCorrMap,'wills',doPlot);
%         end
%         
        figure; imagesc(densityPlotCentresSm); colormap(myColorMap);
        xticks([]); xticklabels({''}); yticks([]); yticklabels({''});
        fname = [figsDir sprintf('/densityPlot_%s_nClus%d_eps%d_batchSiz%d_%s_%s_iter%d_set%d_clusters',dat,iClus+2,epsMuVals(iEps)*1000,batchSizeVals(iBvals),clusPosAct,gridMsrType,iterI,iSet)];
        if savePlots
            set(gcf,'Renderer','painters');
            print(gcf,'-depsc2',fname)
            saveas(gcf,fname,'png');
           close all
        end
        figure; imagesc(aCorrMap,[-.7,.7]); colormap(myColorMap);
        title(sprintf('Grid Score %.2f',g));
        set(gca,'FontSize',fontSiz,'fontname','Arial')
        xticks([]); xticklabels({''}); yticks([]); yticklabels({''});
        fname = [figsDir sprintf('/densityPlot_%s_nClus%d_eps%d_batchSiz%d_%s_%s_iter%d_set%d_autoCorr',dat,iClus+2,epsMuVals(iEps)*1000,batchSizeVals(iBvals),clusPosAct,gridMsrType,iterI,iSet)];
        if savePlots
        %     set(gcf,'Renderer','painters');
        %     print(gcf,'-depsc2',fname)
            saveas(gcf,fname,'png');
           close all
        end
    end
end
end
%%
% 
% if strcmp(dat(1:4),'trap')
%     figure; hold on;
%     dat1    = squeeze(datTmp(iSet,:,:,:,:,:,2));
%     mu      = nanmean(dat1,1);
%     sm      = nanstd(dat1)./sqrt(size(dat1,1));
%     ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
%     plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
%     errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
%     scatter(barpos,mu,100,colors,'d','filled');
%     xlim([barpos(1)-.5, barpos(end)+.5]);
%     ylim([-1.5,1.5]);
%     title(sprintf('Left half of box %s - eps=%d',gridMeasure,epsMuVals(iEps)*1000))
%     
%     figure; hold on;
%     dat1    = squeeze(datTmp(iSet,:,:,:,:,:,3));
%     mu      = nanmean(dat1,1);
%     sm      = nanstd(dat1)./sqrt(size(dat1,1));
%     ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
%     plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
%     errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
%     scatter(barpos,mu,100,colors,'d','filled');
%     xlim([barpos(1)-.5, barpos(end)+.5]);
%     ylim([-1.5,1.5]);
%     title(sprintf('Right half of box %s - eps=%d',gridMeasure,epsMuVals(iEps)*1000))
%     
%     figure; hold on;
%     dat1    = squeeze(datTmp(iSet,:,:,:,:,:,2))-squeeze(datTmp(iSet,:,:,:,:,:,3));
%     mu      = nanmean(dat1,1);
%     sm      = nanstd(dat1)./sqrt(size(dat1,1));
%     ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
%     plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
%     errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
%     scatter(barpos,mu,100,colors,'d','filled');
%     xlim([barpos(1)-.5, barpos(end)+.5]);
%     ylim([-1.5,1.5]);
%     title(sprintf('Left-right half of box %s - eps=%d',gridMeasure,epsMuVals(iEps)*1000))
%     
%     
% end
%     
%     
%     
%% Making figs - cluster positions over time, with agent - needs muAll and trials

savePlots=0;

fntSiz=25;

%if plotAgent, need trials as an output argument
plotAgent = 1;

nClus = size(muAll,1);
iterI = 1;
colors = distinguishable_colors(nClus); %function for making distinguishable colors for plotting
colAgent = [.75, .75, .75];

% make cluster positions same over trials in a batch
muTrls = nan(nClus,2,nTrials);
for iBatch = 1:nBatches
    muTrls(:,1,batchSize*(iBatch-1)+1:batchSize*(iBatch-1)+batchSize)=repmat(muAll(:,1,iBatch),1,batchSize);
    muTrls(:,2,batchSize*(iBatch-1)+1:batchSize*(iBatch-1)+batchSize)=repmat(muAll(:,2,iBatch),1,batchSize);
end

figure;
clear h1
fromTrl = 1;
iTrl = 100; %end trial
%agent - have to start from X trials min
plotTrls=fromTrl:iTrl;
h1(iTrl)=plot(trials(plotTrls,1),trials(plotTrls,2),'Color',colAgent); hold on;
%clusters
h2(iTrl)=scatter(squeeze(muTrls(:,1,iTrl)),squeeze(muTrls(:,2,iTrl)),1000,colors,'.'); hold on;
xticks([]); xticklabels({''}); yticks([]); yticklabels({''});
title(sprintf('Trials 1-100 (start)'),'fontname','Arial','fontsize',fntSiz);

fname = [figsDir sprintf('/learnOverTime_clus_wAgent_%s_nClus%d_eps%d_batchSiz%d_trls%d-%d',dat,nClus,epsMuVals(iEps)*1000,batchSizeVals(iBvals),fromTrl,iTrl)];
if savePlots
        set(gcf,'Renderer','painters');
        print(gcf,'-depsc2',fname)
    saveas(gcf,fname,'png');
    close all
end
        

figure;
clear h1
%plot previous trials plotted above
plotTrls=fromTrl:iTrl;
h1(iTrl)=plot(trials(plotTrls,1),trials(plotTrls,2),'Color',colAgent); hold on;

fromTrl = 10000;
iTrl = fromTrl+2000;
plotTrls=fromTrl:iTrl;
h1(iTrl)=plot(trials(plotTrls,1),trials(plotTrls,2),'Color',colAgent); hold on;

%clusters
h2(iTrl)=scatter(squeeze(muTrls(:,1,iTrl)),squeeze(muTrls(:,2,iTrl)),1000,colors,'.'); hold on;
xticks([]); xticklabels({''}); yticks([]); yticklabels({''});
title(sprintf('Trials %d-%d (1%%)',fromTrl,iTrl),'fontname','Arial','fontsize',fntSiz);

fname = [figsDir sprintf('/learnOverTime_clus_wAgent_%s_nClus%d_eps%d_batchSiz%d_trls%dk-%dk',dat,nClus,epsMuVals(iEps)*1000,batchSizeVals(iBvals),fromTrl/1000,iTrl/1000)];
if savePlots
        set(gcf,'Renderer','painters');
        print(gcf,'-depsc2',fname)
    saveas(gcf,fname,'png');
    close all
end

figure;
clear h1
fromTrl = nTrials*.4;
iTrl = fromTrl+7500;
plotTrls=fromTrl:iTrl;
h1(iTrl)=plot(trials(plotTrls,1),trials(plotTrls,2),'Color',colAgent); hold on;
%clusters
h2(iTrl)=scatter(squeeze(muTrls(:,1,iTrl)),squeeze(muTrls(:,2,iTrl)),1000,colors,'.'); hold on;
xticks([]); xticklabels({''}); yticks([]); yticklabels({''});
title(sprintf('Trials %d-%d (40%%)',fromTrl,iTrl),'fontname','Arial','fontsize',fntSiz);
fname = [figsDir sprintf('/learnOverTime_clus_wAgent_%s_nClus%d_eps%d_batchSiz%d_trls%dkplus7500trls',dat,nClus,epsMuVals(iEps)*1000,batchSizeVals(iBvals),fromTrl/1000)];
if savePlots
        set(gcf,'Renderer','painters');
        print(gcf,'-depsc2',fname)
    saveas(gcf,fname,'png');
    close all
end

figure;
clear h1
fromTrl = nTrials*.75;
iTrl = fromTrl+7500;
plotTrls=fromTrl:iTrl;
h1(iTrl)=plot(trials(plotTrls,1),trials(plotTrls,2),'Color',colAgent); hold on;
%clusters
h2(iTrl)=scatter(squeeze(muTrls(:,1,iTrl)),squeeze(muTrls(:,2,iTrl)),1000,colors,'.'); hold on;
xticks([]); xticklabels({''}); yticks([]); yticklabels({''});
title(sprintf('Trials %d-%d (75%%)',fromTrl,iTrl),'fontname','Arial','fontsize',fntSiz);
fname = [figsDir sprintf('/learnOverTime_clus_wAgent_%s_nClus%d_eps%d_batchSiz%d_trls%dkplus7500trls',dat,nClus,epsMuVals(iEps)*1000,batchSizeVals(iBvals),fromTrl/1000)];
if savePlots
        set(gcf,'Renderer','painters');
        print(gcf,'-depsc2',fname)
    saveas(gcf,fname,'png');
    close all
end
%%


figure;
clear h1
for iTrl = 500
    if mod(iTrl,200)==0 %plot centers after x trials
        %agent
        if plotAgent
            
        %have to start from 1000 trials min
        plotTrls=iTrl-499:iTrl;
        if strcmp(dat(1:3),'cat')
            h1(iTrl)=plot(trials(plotTrls,1),trials(plotTrls,2),'.','Color',repmat(colAgentCh(iTrl),1,3)); hold on;
        else
            h1(iTrl)=plot(trials(plotTrls,1),trials(plotTrls,2),'Color',repmat(colAgentCh(iTrl),1,3)); hold on;
        end
        
        if mod(iTrl,5000)==0
            delete(h1);
        end
        end
        
        %clusters
        h2(iTrl)=scatter(squeeze(muTrls(:,1,iTrl)),squeeze(muTrls(:,2,iTrl)),1000,colors,'.'); hold on;
        drawnow;
%         if mod(iTrl,15000)==0
%             delete(h2);
%         end
        
    end
end

%% Inspect properties of only 'grid' cells (> some threshold)

savePlots=0;

clusPosAct = 'clus'; %'clus' or 'actNorm'

gridMsrType = 'a'; % 'a' or 'w' for allen or willis method - a preferred

gridMeasure = 'grid';

switch clusPosAct
case 'clus'
    switch gridMsrType
        case 'a'
            gridness    = gA_gAll;
            orientation = gA_oAll;
            rad         = gA_radAll;
            wav         = gA_wavAll;
        case 'w'
            gridness    = gW_gAll;
            orientation = gW_oAll;
            rad         = gW_radAll;
            wav         = gW_wavAll;
    end
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


iSet=1; %size(datTmp,1); %last or second last (add -1)
iEps=1;

% clus2plot=(3:26)-2;
% clus2plot=(6:26)-2;
clus2plot = 1:length(clus2run);

iBatchVals=1;

% %fig specs
% xTickLabs = num2cell(clus2run(clus2plot));
% fontSiz=15;
% 
% figure; hold on;
%     for iEps = 1:length(epsMuVals)
%         %     subplot(2,3,iEps);
%         dat1     = squeeze(datTmp(iSet,:,:,iEps,iBatchVals,clus2plot,1));
%         barpos  = .25:.5:.5*size(dat1,2);
%         colors  = distinguishable_colors(size(dat1,2));
%         colgrey = [.6, .6, .6];
%         mu      = mean(dat1,1);
%         sm      = std(dat1)./sqrt(size(dat1,1));
%         ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
%         plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
% %         errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
% %         scatter(barpos,mu,50,colors,'filled','d');
%         scatter(barpos,mu,50,colgrey,'filled','d');
%         xticklabels(xTickLabs);
%         xlim([barpos(1)-.5, barpos(end)+.5]);
%         %         ylim([0,1]);
%         if strcmp(gridMsrType,'a')
%             ylim([-.45,1.4]);
%         elseif strcmp(gridMsrType,'w')
%             ylim([-1.25,1.4]);
%         end
%         xlabel('Number of Clusters');
%         ylabel('Grid Score');        
% %         title(sprintf('%s, %s - eps=%d, batchSize=%d',dat, gridMeasure,epsMuVals(iEps)*1000,batchSizeVals(iBatchVals)))
%         if strcmp(dat(1:2),'ci')
%         title('Circular box')
%         elseif strcmp(dat(1:2),'sq')
%         title('Square box')
%         end        
%         set(gca,'FontSize',fontSiz,'fontname','Arial')
%     end



gridInd=dat1>0.35;


gScore=squeeze(gridness);
oScore=squeeze(orientation);
rScore=squeeze(rad);
wScore=squeeze(wav);


figure; hold on;
for iClus = 1:length(clus2run)
    subplot(6,5,iClus);
    hist(oScore(gridInd(:,iClus),iClus),15);
%     hist(rScore(gridInd(:,iClus),iClus),15);
end



    