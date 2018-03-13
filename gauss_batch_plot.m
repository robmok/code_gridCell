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


dat='circ';
nIter=200;
epsMuVals = [.0075, .01, .02]; 


clus2run = [16, 20, 24]; %14 18, 28 still running
nTrials=2500000;
batchSizeVals=[250, 500, 1000, 2000];


%sigmaGauss - temp
stepSize = 1;
sigmaGaussVals = round([stepSize/3, stepSize/3.5].*100);
sigmaGauss100  = sigmaGaussVals(2);


%load loop
for iClus2run = 1:length(clus2run) 
    nClus = clus2run(iClus2run);
    for iEps = 1:length(epsMuVals) 
        epsMuOrig=epsMuVals(iEps);
        epsMuOrig10000=epsMuOrig*10000;
        for iBvals = 1:length(batchSizeVals)
            
            %fixed batch
            if fixBatchSize
                batchSize = batchSizeVals(iBvals);
                fprintf('Loading nClus=%d, epsMu=%d, batchSize=%d\n',nClus,epsMuOrig10000,batchSize)
                fname = [saveDir, sprintf('/gauss_batch_%dclus_%dsigma_%dktrls_eps%d_batchSiz%d_%diters_%s*',nClus,sigmaGauss100,round(nTrials/1000),epsMuOrig10000,batchSize,nIter,dat)];
                
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
                gW_gAll(:,:,iEps,iBvals,iClus2run,:) = gW(:,:,1,:);
                gW_oAll(:,:,iEps,iBvals,iClus2run,:) = gW(:,:,2,:);
                gW_radAll(:,:,iEps,iBvals,iClus2run,:) = gW(:,:,3,:);
                gW_wavAll(:,:,iEps,iBvals,iClus2run,:) = gW(:,:,4,:);

                
            end
        end
    end
end

%% plot univar scatters

gridMsrType = 'a'; % 'a' or 'w' for allen or willis method

gridMeasure = 'grid';

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


iSet=6;

%plot univar scatters - over clusters (e.g. one batch size)
if size(gridness,4)==1 %only 1 batchSize
    figure; hold on;
    dat1     = squeeze(datTmp(iSet,:,:,:,:,1));
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
    
    if strcmp(dat(1:5),'trapz')
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


%plot univar scatters - comparing batchSizeVals, with learning rates in sep
%plts, with clusters in separate figs
if size(gridness,4)>1 % more than 1 batchSize
    for iClus2Run = 1:length(clus2run)
        figure; hold on;
        for iEps = 1:length(epsMuVals)
            % comparing batch vals, with cluster nums in subplots; can edit for eps
            %             subplot(ceil(length(epsMuVals)/2),2,iEps);
            subplot(1,3,iEps);
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
            ylim([-.5,1.25]);
            title(sprintf('%s - eps=%d, nClus=%d',gridMeasure,epsMuVals(iEps)*10000,clus2run(iClus2Run)))
            
        end
    end
    
end


%plot univar scatters - comparing batchSizeVals, with clusters in separate plots
% figure; hold on;
% if size(gridness,4)>1 % more than 1 batchSize
%     for iEps = 1:length(epsMuVals)
%         for iClus2Run = 1:length(clus2run)
%             % comparing batch vals, with cluster nums in subplots
%             subplot(ceil(length(clus2run)/2),2,iClus2Run);
%             %     subplot(2,3,iEps);
%             dat1     = squeeze(datTmp(iSet,:,iEps,:,iClus2Run,1));
%             barpos  = .25:.5:.5*size(dat1,2);
%             colors  = distinguishable_colors(size(dat1,2));
%             colgrey = [.5, .5, .5];
%             mu      = mean(dat1,1);
%             sm      = std(dat1)./sqrt(size(dat1,1));
%             ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
%             plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
%             errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
%             scatter(barpos,mu,750,colors,'x');
%             xlim([barpos(1)-.5, barpos(end)+.5]);
%             %         ylim([0,1]);
%             ylim([-.5,1.25]);
%             title(sprintf('%s - eps=%d, nClus=%d',gridMeasure,epsMuVals(iEps)*10000,clus2run(iClus2Run)))
%             
%             
%             
%             
%             % % comparing learning rate, with batch vals in subplots
%             % figure; hold on;
%             % for iBvals = 1:length(batchSizeVals)
%             %     subplot(2,3,iBvals);
%             %     dat1     = squeeze(datTmp(iSet,:,:,iBvals));
%             %     barpos  = .25:.5:.5*size(dat1,2);
%             %     colors  = distinguishable_colors(size(dat1,2));
%             %     colgrey = [.5, .5, .5];
%             %     mu      = mean(dat1,1);
%             %     sm      = std(dat1)./sqrt(size(dat1,1));
%             %     ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
%             %     plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
%             %     errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
%             %     scatter(barpos,mu,750,colors,'x');
%             %     xlim([barpos(1)-.5, barpos(end)+.5]);
%             %     ylim([-.5,1.25]);
%             %     title(sprintf('%s - batchVal=%d',gridMeasure,batchSizeVals(iBvals)))
%             % end
%         end
%     end
%     
% end





%plot over sets
if 0
if size(gridness,4)>1 
    figure; hold on;
    for iBatchVals = 1:length(batchSizeVals)        
        subplot(ceil(length(batchSizeVals)/2),2,iBatchVals); hold on;
        dat1     = squeeze(datTmp(:,:,:,iBatchVals,:,:))';
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
        title(sprintf('%s - eps=%d',gridMeasure,epsMuVals(iEps)*1000))
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
iBvals    = 3;

iters2plot = 40:60;

% fprintf(sprintf('clus %d batchSizeVals %d\n',clus2run(iClus2run),batchSizeVals(iBvals)));
for iterI = iters2plot
    densityPlotCentresSm = imgaussfilt(densityPlotAll(:,:,iSet,iterI,iEps,iBvals,iClus2run),gaussSmooth);
    
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

%% plot gauss activations (batch)
% 
% densityPlotAct        = zeros(50,50);
% densityPlotAct1       = zeros(50,50); %one cluster
% figure; hold on;
% xlim([0,50]); ylim([0,50]);
% 
% 
% iterI = 1; 
% %seed to regenerate trials
% rng(rSeed(iterI));
% 
% locRange = [0, 49];
% spacing  = linspace(locRange(1),locRange(2),locRange(2)+1);
% trials   = [randsample(locRange(1):diff( v(1:2)):locRange(2)*2,nTrialsTest,'true'); randsample(spacing,nTrialsTest,'true')]';
% 
% zlims=[15 40];
% for iTrl = nTrials*.5:nTrials
%     
%     densityPlotAct(trials(iTrl,1)+1, trials(iTrl,2)+1) = densityPlotAct(trials(iTrl,1)+1, trials(iTrl,2)+1)+ sum(actAll(:,iTrl)).^2; %squaring might make it look better..
%     if mod(iTrl,2000)==0 %plot centers after x trials
% %     imagesc(densityPlotAct,zlims); drawnow;
%     imagesc(imgaussfilt(densityPlotAct,1)); drawnow;
%     end
%     
% %     densityPlotAct1(trials(iTrl,1)+1, trials(iTrl,2)+1) = densityPlotAct1(trials(iTrl,1)+1, trials(iTrl,2)+1)+ (actAll(1,iTrl)).^2;
% %     if mod(iTrl,2000)==0 %plot centers after x trials
% %         imagesc(imgaussfilt(densityPlotAct1,1)); drawnow;
% %     end
%     
%     
% end
%     %%
% zlims=[10 20];
% 
% figure;
% imagesc(densityPlotAct,zlims);
% figure;
% imagesc(imgaussfilt(densityPlotAct,1),zlims)
% % 
% % zlims1=[.5 1.5];
% % figure;
% % imagesc(densityPlotAct1,zlims1);
% % figure;
% % imagesc(imgaussfilt(densityPlotAct1,1),zlims1)
% %     
% %%
% aCorrMap=ndautoCORR(imgaussfilt(densityPlotAct,1)); %autocorrelogram
% [g,gdataA] = gridSCORE(aCorrMap,'allen',1);