clear all;

% wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
cd(wd);

codeDir = [wd '/code_gridCell'];
saveDir = [wd '/data_gridCell'];
addpath(codeDir); addpath(saveDir);
addpath(genpath([codeDir '/gridSCORE_packed']));

% load
nSet        = 6;
gaussSmooth = 1; 
fixBatchSize = 1; %fixed batch size or depend on nClus (for fname)

dat='square';
annEps=0;

% clus2run = [7,8,10,12]; 
% nTrials = 100000; 
% nIter=250;
% batchSizeVals = [1, 50, 100, 200, 500];
% epsMuVals=[.01, .05, .075, .1, .2, .3];% learning rate / starting learning rate 

% sims 1 - one single LARGE batchSize, large nTrials, all nClus conds
% clus2run = 3:30;%[7,8,10,12]; 
% nTrials=10000000;
% batchSizeVals=1000;
% nIter=200;
% epsMuVals=.075;

% sims 2 - a smaller val of trials; testing batch sizes (works fine) - also
% have some sims with ntrials = 5000000 (less batchSizeVals)
% clus2run = [10:2:28]; % also have 11 but left out plot
% clus2run = 28; % plotting one nClus cond over sets
nTrials=2500000;
nIter=200;
% batchSizeVals=[13, 25, 83, 125, 167, 250, 333, 500, 1000, 2000];
%wAct
clus2run = [10:2:30];
clus2run = [10, 12 16:2:28]; %no 14 and 30 for sq, no 12, 30 for circ; epsmu=00.75
% batchSizeVals=[1000, 2000];
batchSizeVals = [1000, 500, 125, 50];
epsMuVals=.075;

 %new
%  clus2run = [10, 12, 16:2:26]; %no 14? 28? - errored, running (0.15, circ)
% epsMuVals=.015;
% nTrials=2000000; %new
% batchSizeVals = [800, 400, 100, 40];


% new 3 - fixed batch sizes across clusters
% nTrials=2500000;
% clus2run = [18:2:30]; 
% batchSizeVals=[333, 500, 1000];

% new 4 - batchSizes based on mean updates per clus per batch (avgBatch)
% fixBatchSize = 0;
% clus2run = [18:2:30]; 
% % clus2run = [18 22 28]; % running these on love01
% batchSizeVals=[1, 2, 5, 10, 25, 35, 50]; %avgBatchVals - 1,2,5 new

%new 5 - circ
% fixBatchSize = 1;
% clus2run = [10:2:20, 22:26, 28, 30];
% % clus2run = 14; % plotting one nClus cond over sets
% nTrials=2500000;
% % batchSizeVals=[125, 167, 250, 333, 500, 1000, 2000];
% batchSizeVals=[167, 250, 333, 500, 1000, 2000];
% %wAct
% clus2run = [16 18 20 24];
% batchSizeVals=[1000, 2000];
% dat='circ';

%new 5 - trapz
% fixBatchSize = 1;
% clus2run = [12:2:18 22, 26]; % trapz1-3, krupic
% clus2run = [16, 18, 20, 22, 24, 26, 28, 30]; %new - more clusters (trapzScaled1-3, krupic2,3)
% clus2run = [16, 18, 20, 24, 26, 30]; %krupic3 done
% nTrials=2500000;
% batchSizeVals=1000;
% dat='trapzKrupic3'; %trapz1, trapz2, trapz3, trapzKrupic, trapzNorm (or rename to trapz)

%trapz size x2
% clus2run = [16, 18]; 
% nTrials=2500000;
% batchSizeVals=1000;
% dat='trapzScaled3'; %trapz1, trapz2, trapz3, trapzKrupic, trapzNorm (or rename to trapz)
% boxSize=2;


%new - annealed learning rate
% clus2run = [10:2:30]; 
% epsMuVals=[.2, .05]; %"starting" learning rate - actually these are just numbers, not actualy starting eps
% nTrials=2000000; %new
% % batchSizeVals = [400, 100]; % only 1 batchsize val per learning rate
% batchSizeVals = 0;
% annEps=1;
% covering_map_batch_dat_22clus_2000ktrls_eps200_batchSiz100_200iters_circ_wAct_annEps_132454


%load loop
for iClus2run = 1:length(clus2run) 
    nClus = clus2run(iClus2run);
    for iEps = 1:length(epsMuVals) 
        epsMuOrig=epsMuVals(iEps);
        epsMuOrig1000=epsMuOrig*1000;
        for iBvals = 1:length(batchSizeVals)
            
%             %fixed batch
%             if fixBatchSize
%                 batchSize = batchSizeVals(iBvals);
%                 fprintf('Loading nClus=%d, epsMu=%d, batchSize=%d\n',nClus,epsMuOrig1000,batchSize)
%                 if ~strcmp(dat,'square')
%                     fname = [saveDir, sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s*',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat)];
%                 else
%                     fname = [saveDir, sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters*',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter)];
%                 end
%             else %avgBatch
%                 avgBatch = batchSizeVals(iBvals);
%                 fprintf('Loading nClus=%d, epsMu=%d, avgBatchSize=%d\n',nClus,epsMuOrig1000,avgBatch)
%                 if ~strcmp(dat,'square')
%                     fname = [saveDir, sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_avgBatch%d_batchSiz*_%diters*_%s',nClus,round(nTrials/1000),epsMuOrig1000,avgBatch,nIter,dat)];
%                 else
%                     fname = [saveDir, sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_avgBatch%d_batchSiz*_%diters*',nClus,round(nTrials/1000),epsMuOrig1000,avgBatch,nIter)];
%                 end
%             end
            

            %tmp - wAct and boxSizex2
            batchSize = batchSizeVals(iBvals);
            fprintf('Loading nClus=%d, epsMu=%d, batchSize=%d\n',nClus,epsMuOrig1000,batchSize)
%             fname = [saveDir, sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_wAct_%s*',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat)];
            fname = [saveDir, sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wAct_%s*',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat,dat)];
%             fname = [fname '_boxSizex2*'];
            
            if annEps
                fname = [saveDir, sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz*_%diters_%s_wAct_annEps*',nClus,round(nTrials/1000),epsMuOrig1000,nIter,dat)];
            end


            
            %edit if want to load more than one file per sim, merge
            f = dir(fname); filesToLoad = cell(1,length(f));
            for iF = 1%:length(f)
                filesToLoad{iF} = f(iF).name;
                load(f(iF).name);
            end
            
            for iterI = 1:nIter
                for iSet=1:nSet
                    densityPlotAll(:,:,iSet,iterI,iEps,iBvals,iClus2run) = imgaussfilt(densityPlot(:,:,iSet,iterI),gaussSmooth);
                end
            end
            %organise gridness values (allen vs willis method)
            gA_gAll(:,:,iEps,iBvals,iClus2run,:)   = gA(:,:,1,:);
            gA_oAll(:,:,iEps,iBvals,iClus2run,:)   = gA(:,:,2,:);
            gA_radAll(:,:,iEps,iBvals,iClus2run,:) = gA(:,:,3,:);
            gA_wavAll(:,:,iEps,iBvals,iClus2run,:) = gA(:,:,4,:);
%             gW_gAll(:,:,iEps,iBvals,iClus2run,:) = gW(:,:,1,:);
%             gW_oAll(:,:,iEps,iBvals,iClus2run,:) = gW(:,:,2,:);
%             gW_radAll(:,:,iEps,iBvals,iClus2run,:) = gW(:,:,3,:);
%             gW_wavAll(:,:,iEps,iBvals,iClus2run,:) = gW(:,:,4,:);
            
                        
            if exist('rSeed','var')
                rSeedAll(:,iEps,iBvals,iClus2run) = rSeed;
            end
            
            if exist('gA_actNorm','var') %if have gridness on gauss act
                for iterI = 1:nIter
                    for iSet=1:nSet
%                         densityPlotActAll(:,:,iSet,iterI,iEps,iBvals,iClus2run) = imgaussfilt(densityPlotAct(:,:,iSet,iterI),gaussSmooth);
                        densityPlotActNormAll(:,:,iSet,iterI,iEps,iBvals,iClus2run) = imgaussfilt(densityPlotActNorm(:,:,iSet,iterI),gaussSmooth);
                    end
                end
                
%                 gA_gAll_act(:,:,iEps,iBvals,iClus2run,:)   = gA_act(:,:,1,:);
%                 gA_oAll_act(:,:,iEps,iBvals,iClus2run,:)   = gA_act(:,:,2,:);
%                 gA_radAll_act(:,:,iEps,iBvals,iClus2run,:) = gA_act(:,:,3,:);
%                 gA_wavAll_act(:,:,iEps,iBvals,iClus2run,:) = gA_act(:,:,4,:);
%                 gW_gAll_act(:,:,iEps,iBvals,iClus2run,:) = gW_act(:,:,1,:);
%                 gW_oAll_act(:,:,iEps,iBvals,iClus2run,:) = gW_act(:,:,2,:);
%                 gW_radAll_act(:,:,iEps,iBvals,iClus2run,:) = gW_act(:,:,3,:);
%                 gW_wavAll_act(:,:,iEps,iBvals,iClus2run,:) = gW_act(:,:,4,:);
%                 
                gA_gAll_actNorm(:,:,iEps,iBvals,iClus2run,:)   = gA_actNorm(:,:,1,:);
                gA_oAll_actNorm(:,:,iEps,iBvals,iClus2run,:)   = gA_actNorm(:,:,2,:);
                gA_radAll_actNorm(:,:,iEps,iBvals,iClus2run,:) = gA_actNorm(:,:,3,:);
                gA_wavAll_actNorm(:,:,iEps,iBvals,iClus2run,:) = gA_actNorm(:,:,4,:);
%                 gW_gAll_actNorm(:,:,iEps,iBvals,iClus2run,:) = gW_actNorm(:,:,1,:);
%                 gW_oAll_actNorm(:,:,iEps,iBvals,iClus2run,:) = gW_actNorm(:,:,2,:);
%                 gW_radAll_actNorm(:,:,iEps,iBvals,iClus2run,:) = gW_actNorm(:,:,3,:);
%                 gW_wavAll_actNorm(:,:,iEps,iBvals,iClus2run,:) = gW_actNorm(:,:,4,:);

                gA_gAll_actNorm(:,:,iEps,iBvals,iClus2run,:)   = gA_actNorm(:,:,1,:);
                gA_oAll_actNorm(:,:,iEps,iBvals,iClus2run,:)   = gA_actNorm(:,:,2,:);
                gA_radAll_actNorm(:,:,iEps,iBvals,iClus2run,:) = gA_actNorm(:,:,3,:);
                gA_wavAll_actNorm(:,:,iEps,iBvals,iClus2run,:) = gA_actNorm(:,:,4,:);

            end

            
        end
    end
end
%% plot univar scatters

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
                orientation = gW_oAll_act;
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
    dat1     = squeeze(datTmp(iSet,:,iEps,:,:,1));
    barpos  = .25:.5:.5*size(dat1,2);
    colors  = distinguishable_colors(size(dat1,2));
    colgrey = [.5, .5, .5];
    mu      = mean(dat1,1);
    sm      = std(dat1)./sqrt(size(dat1,1));
    ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
    plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
    errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
    % scatter(barpos,mu,750,colors,'x');
    % scatter(barpos,mu,750,colors,'.');
    scatter(barpos,mu,100,colors,'d','filled');
    xlim([barpos(1)-.5, barpos(end)+.5]);
    ylim([-.5,1.25]);
    title(sprintf('%s - eps=%d',gridMeasure,epsMuVals(iEps)*1000))
    
    if strcmp(dat(1:4),'trap')
        figure; hold on;
        dat1    = squeeze(datTmp(iSet,:,:,:,:,2));
        mu      = nanmean(dat1,1);
        sm      = nanstd(dat1)./sqrt(size(dat1,1));
        ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
        plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
        errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
        scatter(barpos,mu,100,colors,'d','filled');
        xlim([barpos(1)-.5, barpos(end)+.5]);
        ylim([-.5,1.25]);
        title(sprintf('Left half of box %s - eps=%d',gridMeasure,epsMuVals(iEps)*1000))
        
        figure; hold on;
        dat1    = squeeze(datTmp(iSet,:,:,:,:,3));
        mu      = nanmean(dat1,1);
        sm      = nanstd(dat1)./sqrt(size(dat1,1));
        ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
        plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
        errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
        scatter(barpos,mu,100,colors,'d','filled');
        xlim([barpos(1)-.5, barpos(end)+.5]);
        ylim([-.5,1.25]);
        title(sprintf('Right half of box %s - eps=%d',gridMeasure,epsMuVals(iEps)*1000))
        
        figure; hold on;
        dat1    = squeeze(datTmp(iSet,:,:,:,:,2))-squeeze(datTmp(iSet,:,:,:,:,3));
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
    
%     %plot hist
%     figure; pltCount=1;
%     iToPlot = 2:length(clus2run);
%     for iClus2Run = iToPlot
%         subplot(3,3,pltCount)
%         hist(dat1(:,iClus2Run),40); %50
%         xlim([-.5,1.25]);
%         title(sprintf('%d',iClus2Run+2));
%         if mod(pltCount,9)==0 && iClus2Run~=iToPlot(end)
%             pltCount=1; figure;
%         else
%             pltCount = pltCount+1;
%         end
%     end
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
        dat1     = squeeze(datTmp(iSet,:,iEps,:,iClus2Run,1));
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
   
    
    % % comparing learning rate, with batch vals in subplots
    % figure; hold on;
    % for iBvals = 1:length(batchSizeVals)
    %     subplot(2,3,iBvals);
    %     dat1     = squeeze(datTmp(iSet,:,:,iBvals));
    %     barpos  = .25:.5:.5*size(dat1,2);
    %     colors  = distinguishable_colors(size(dat1,2));
    %     colgrey = [.5, .5, .5];
    %     mu      = mean(dat1,1);
    %     sm      = std(dat1)./sqrt(size(dat1,1));
    %     ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
    %     plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
    %     errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
    %     scatter(barpos,mu,750,colors,'x');
    %     xlim([barpos(1)-.5, barpos(end)+.5]);
    %     ylim([-.5,1.25]);
    %     title(sprintf('%s - batchVal=%d',gridMeasure,batchSizeVals(iBvals)))
    % end
end

end





%plot over sets
if 1
   
if size(gridness,3)>1 || size(gridness,4)>1 
    for iClus = 1:length(clus2run)
    figure; hold on;
    for iBatchVals = 1:length(batchSizeVals)        
        subplot(ceil(length(batchSizeVals)/2),2,iBatchVals); hold on;
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
    
%     if strcmp(dat,'trapz') || strcmp(dat,'trapzNorm')
%         figure; hold on;
%         dat1     = squeeze(datTmp(iSet,:,:,:,:,2));
%         mu      = mean(dat1,1);
%         sm      = std(dat1)./sqrt(size(dat1,1));
%         ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
%         plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
%         errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
%         scatter(barpos,mu,100,colors,'d','filled');
%         xlim([barpos(1)-.5, barpos(end)+.5]);
%         ylim([-.5,1.25]);
%         title(sprintf('Left half of box %s - eps=%d',gridMeasure,epsMuVals(iEps)*1000))
%         
%         figure; hold on;
%         dat1     = squeeze(datTmp(iSet,:,:,:,:,3));
%         mu      = mean(dat1,1);
%         sm      = std(dat1)./sqrt(size(dat1,1));
%         ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
%         plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
%         errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
%         scatter(barpos,mu,100,colors,'d','filled');
%         xlim([barpos(1)-.5, barpos(end)+.5]);
%         ylim([-.5,1.25]);
%         title(sprintf('Right half of box %s - eps=%d',gridMeasure,epsMuVals(iEps)*1000))
%     end
    end
end
end








%testing one sim only

% dat1=squeeze(gA(iSet,:,1));
% figure;
% barpos  = .25:.5:.5*size(dat1,1);
% colors  = distinguishable_colors(size(dat1,1));
% colgrey = [.5, .5, .5];
% mu      = mean(dat1);
% sm      = std(dat1)./sqrt(size(dat1,2));
% ci      = sm.*tinv(.025,size(dat1,2)-1); %compute conf intervals
% plotSpread(dat1','xValues',barpos);
% errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
% scatter(barpos,mu,750,'x');
% xlim([barpos(1)-.5, barpos(end)+.5]);
% ylim([-.5,1.25]);


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
        densityPlotCentresSm = imgaussfilt(densityPlotAll(:,:,iSet,iterI,iEps,iBvals,iClus2run),gaussSmooth);
        
        if doPlot
            figure; hold on;
            subplot(1,2,1)
            imagesc(densityPlotCentresSm);
            % figure;
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

iSet=6;
iEps=1;
gaussSmooth=1;


if strcmp(dat(1:4),'trap')
subPlt = [2,2];
else
    subPlt = [1,2];
end

%set
iClus2run = 2;
iBvals    = 1;

iters2plot = 40:45;

fprintf(sprintf('clus %d batchSizeVals %d\n',clus2run(iClus2run),batchSizeVals(iBvals)));
for iterI = iters2plot
%     densityPlotCentresSm = imgaussfilt(densityPlotAll(:,:,iSet,iterI,iEps,iBvals,iClus2run),gaussSmooth);
%     densityPlotCentresSm = imgaussfilt(densityPlotActAll(:,:,iSet,iterI,iEps,iBvals,iClus2run),gaussSmooth);
    densityPlotCentresSm = imgaussfilt(densityPlotActNormAll(:,:,iSet,iterI,iEps,iBvals,iClus2run),gaussSmooth);
    
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
        aCorrMap = ndautoCORR(densityPlotCentresSm(:,1:size(densityPlotAll,1)/2));
        [g,gdataA] = gridSCORE(aCorrMap,'allen',1);
        subplot(subPlt(1),subPlt(2),4);
        aCorrMap = ndautoCORR(densityPlotCentresSm(:,size(densityPlotAll,1)/2+1:end));
        [g,gdataA] = gridSCORE(aCorrMap,'allen',1);
        
    end
    
end

%% regen trials and activation?
% 
% 
% %set
% iClus2run = 3;
% iBvals    = 3;
% 
% iters2plot = 40:45;
% 
% locRange = [0 49];
% 
% % fprintf(sprintf('clus %d batchSizeVals %d\n',clus2run(iClus2run),batchSizeVals(iBvals)));
% for iterI = iters2plot
% 
%     rng(rSeedAll(iterI,iEps,iBvals,iClus2run));
%     
%     trials = [randsample(linspace(locRange(1),locRange(2),50),nTrials,'true'); randsample(linspace(locRange(1),locRange(2),50),nTrials,'true')]';
%     
%     
%     %need mu starting points
%     
%     
%     %for each trial, compute closest cluster and save activation + location
%     
%     
%     %need also to update the cluster means over batches...
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
% end
% 
% 
%         
% 
