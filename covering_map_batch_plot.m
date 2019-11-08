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
annEps        = 1;

rHex = 0; %if inspect raw 60deg corr values, rather than grid score
% 
% %trapzKfrmSq 
% dat           = 'trapzKfrmSq1';
% nTrials       = 250000;
% batchSizeVals = 200;
% clus2run      = 10:30;
% nIter         = 1000;
% actOverTime   = 1;

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
                iSet = 1:nSet; %useful for when iSet=1
                 
                gA_gAll_actNorm(iSet,:,:,iEps,iBvals,iClus2run,:)   = gA_actNorm(iSet,:,1,:);
                if rHex
                    gA_gAll_actNorm(:,:,iEps,iBvals,iClus2run,:)   = (gA_actNorm(:,:,6,:)+gA_actNorm(:,:,8,:))./2; %rHex 60/120
%                     gA_gAll_actNorm(iSet,:,:,iEps,iBvals,iClus2run,:)   = (gA_actNorm(iSet,:,5,:)+gA_actNorm(iSet,:,7,:)+(gA_actNorm(iSet,:,9,:)))./3; %rHex 30,90,150
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

%% Making figs: plot univariate scatters over time

savePlots = 0;
plotSubPlots = 0;

fontSiz = 25;
iBatchVals=1;
nTimePts=size(datTmp,1)-1;
colors  = distinguishable_colors(length(clus2run));

if plotSubPlots
    clus2plot = 1:length(clus2run);
    figure; hold on;
    for iClus = clus2plot
        subplot(4,6,iClus-clus2plot(1)+1); hold on;
        dat1     = squeeze(datTmp(1:nTimePts,:,iEps,iBatchVals,:,iClus))';
        barpos  = .25:.5:.5*size(dat1,2);
        colgrey = [.5, .5, .5];
        mu      = nanmean(dat1,1);
        sm      = nanstd(dat1)./sqrt(size(dat1,1));
        plotSpread(dat1,'xValues',barpos,'distributionColors',repmat(colors(iClus,:),nTimePts,1));
        scatter(barpos,mu,10,colgrey,'d','filled');
        xlim([barpos(1)-.5, barpos(end)+.5]);
        ylim([-.5,1.25]);
        xticks([]); xticklabels({''});
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
    
else % plot single plots - but better to plot subset, else too many
    %select subset of nClus conds
    clus2plot = [10,12,18,25]-9; % plot selection
    %     clus2plot = 1:length(clus2run); % plot all
    for iClus = clus2plot
        figure;
        dat1     = squeeze(datTmp(1:nTimePts,:,:,iEps,iBatchVals,iClus,:))';
        barpos  = .25:.5:.5*size(dat1,2);
        colgrey = [.5, .5, .5];
        mu      = nanmean(dat1,1);
        sm      = nanstd(dat1)./sqrt(size(dat1,1));
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

iSet=size(datTmp,1); %last one (final quarter trials)
iEps=1;
clus2plot = 1:length(clus2run);
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
        plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
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

%overall mean
ci = bootci(numel(dat1),@nanmean,reshape(dat1,1,numel(dat1))); %CI over ALL runs
fprintf('%d to %d clusters: mean=%0.4f; CI=[%0.4f,%0.4f]\n',clus2run(clus2plot(1)),clus2run(clus2plot(end)),nanmean(reshape(dat1,1,numel(dat1))),ci(1),ci(2));

%% Making figs: density plot examples

savePlots = 0;

doPlot=0; %do plot when computing gridness

fontSiz=15;
clusPosAct = 'actNorm'; %'clus' or 'actNorm'
gridMsrType = 'a';
clus2plot = ([10,11,12,18,20,23,25,28])-9; %minus 9 because starts from clus=10
% clus2plot = (18)-9;

iSet=size(densityPlotActNormAll,3); %last one
myColorMap = parula;
myColorMap(end,:) = 1; %make bits outside of shape white
for iClus = clus2plot%:length(clus2run)
    for iBvals = 1:length(batchSizeVals)
        for iterI = 1%:5%nIter
            
            densityPlotCentresSm = densityPlotActNormAll(:,:,iSet,iterI,iEps,iBvals,iClus);
            aCorrMap=ndautoCORR(densityPlotCentresSm); %autocorrelogram
            densityPlotCentresSm(isnan(densityPlotCentresSm))=1.1;
            aCorrMap(isnan(aCorrMap))=1.1;
            [g,gdataA] = gridSCORE(aCorrMap,'allen',doPlot);
            %         [g,gdataW] = gridSCORE(aCorrMap,'wills',doPlot);
            
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
                set(gcf,'Renderer','painters');
                print(gcf,'-depsc2',fname)
                saveas(gcf,fname,'png');
                close all
            end
        end
    end
end
%% Making figs - cluster positions over time, with agent - needs muAll and trials

savePlots=0;

%if plotAgent, need trials as an output argument
plotAgent = 1;

fntSiz=25;
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