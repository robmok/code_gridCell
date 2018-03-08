%% kmeans plot

% clear all

% wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
cd(wd);

codeDir = [wd '/code_gridCell'];
saveDir = [wd '/data_gridCell'];
addpath(codeDir); addpath(saveDir);


dat = 'circ'; %'square', 'circ'

kVals   = 3:30; 
nKmeans = 1000;
nPoints = 3000; %3k, 5k, 10k
nKvals  = length(kVals);

nDataPtsTest = 3000; %3k, 5k, 10k
nXvalDataSets = 20;

%load in data
fname=[saveDir sprintf('/kmeans_nK_3-30_%s_nPoints%d_1000iters',dat,nPoints)];
load(fname);

% load in xVal
fname = [saveDir sprintf('/kmeans_xVal_nK_%d-%d_%s_%dtrainPts_%dtestPts_%diters_%ddatasets',kVals(1),kVals(end),dat,nPoints, nDataPtsTest,nKmeans,nXvalDataSets)];
load(fname);

tsseXval=xVal_results.tsseXval;
sseSprdSdXval=xVal_results.sseSprdSdXval;
sseSprdVarXval=xVal_results.sseSprdVarXval;

%comparing between nPoints
% figure;
% plot(squeeze((mean(mean(tsseXval,2),1))));
% hold on;
%% plot hist and density  plots

gridMsr = 'a'; % 'a' or 'w' for allen or willis method

switch gridMsr
    case 'a'
        g = gA;
    case 'w'
        g = gW;
end

%plot univar scatters
figure;
dat1     = squeeze(g(:,1,:));
barpos  = .25:.5:.5*size(dat1,2);
colors  = distinguishable_colors(size(dat1,2));
colgrey = [.5, .5, .5];
mu      = mean(dat1,1);
sm      = std(dat1)./sqrt(size(dat1,1));
ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
scatter(barpos,mu,750,colors,'x','LineWidth',1);
xlim([barpos(1)-.5, barpos(end)+.5]);
ylim([-.5,1.25]);

%plot some grids

% densityPlotCentres(:,:,kMeansIter,iKvals) = sum(densityPlotClus,3);
% densityPlotCentresSm = imgaussfilt(densityPlotCentres(:,:,kMeansIter,iKvals),gaussSmooth);
% aCorrMap=ndautoCORR(densityPlotCentresSm); %autocorrelogram
% [~,gdataA] = gridSCORE(aCorrMap,'allen',0);
% [~,gdataW] = gridSCORE(aCorrMap,'wills',0);

%% Crossvalidation - SSE on the test dataset
% note: tsseXval is already ordered from low to high SSE on train set

% iToPlot = 1:28;
% % SSE on test data,  showing lowest SSE on training data
% figure; hold on;
% iPlt=0;
% for iKvals = iToPlot %1:nKvals
%     iPlt=iPlt+1;
%     subplot(6,5,iPlt); hold on;
%     nK = kVals(iKvals);
%     plot(tsseXval(:,:,iKvals)','Color',[.75 .75 .75]);
%     plot(tsseXval(indSSE2(iKvals,:)==1,:,iKvals)','Color',[.2 .2 .2]); 
% %     title(sprintf('SSE on Test data showing lowest SSE on training data (k=%d) (nTest=%d)',nK,nXvalDataSets));
% end

% cross-validation - plot SSE over datasets
figure; hold on;
iPlt=0;
for iKvals = iToPlot
    iPlt=iPlt+1;
    subplot(6,5,iPlt); hold on;
    nK = kVals(iKvals);
    plot(tsseXval(:,:,iKvals)','Color',[.75 .75 .75]);
%     plot(tsseXval(min(nanmean(tsseXval(:,:,iKvals),2))==nanmean(tsseXval(:,:,iKvals),2),:,iKvals),'Color',[.5 .5 .5],'LineWidth',2); %lowest mean SSE
    plot(tsseXval(indSSE2(iKvals,:)==1,:,iKvals)','Color',[.2 .2 .2],'LineWidth',2); %lowest mean SSE from training data
%     title(sprintf('SSE over all Test data sets,  showing lowest SSE on training data (k=%d) (nTest=%d)(datPts=%dk)',nK,nXvalDataSets,nPoints/1000));

end

%% Scatters - corr with gridness

%which measure of error; SSE or SSE per cluster
err=tssekVals;
%need to add / savev spread out measure for the kmeans orig?

xValErr = tsseXval;
% xValErr = sseSprdSdXval; %works well too
%sseSprdVarXval - gives same corr values


% corr of gridness and SSE: original data
iToPlot = [1:28];
figure; hold on;
iPlt=0;
for iKvals = iToPlot
    iPlt=iPlt+1;
    subplot(6,5,iPlt)
    [r p] = corr(g(:,1,iKvals),err(iKvals,:)','type','spearman');
    r1(iKvals)=r;
    p1(iKvals)=p;
    p1corr(iKvals)=p*nKvals;
    scatter(err(iKvals,:)',g(:,1,iKvals),'.');
    title(sprintf('r = %.2f, p (crcted) = %.5f',r,p*nKvals))
end
[r1; p1corr]'

% corr of gridness and SSE: cross-validated data
iToPlot = [1:28];
figure; hold on;
iPlt=0;
for iKvals = iToPlot
    iPlt=iPlt+1;
    subplot(6,5,iPlt)
    [r p] = corr(g(:,1,iKvals),nanmean(xValErr(:,:,iKvals),2),'type','spearman');
    r1(iKvals)=r;
    p1(iKvals)=p;
    p1corr(iKvals)=p*nKvals;
    scatter(nanmean(xValErr(:,:,iKvals),2),g(:,1,iKvals),'.');
    title(sprintf('r = %.2f, p (crcted) = %.5f',r,p*nKvals))
end

[r1; p1corr]'


% are 30/60 degree orientations --> better gridness? nk=3 suggests so.
% others not sure; most orientations at 0, 90 and 40/50.. (3 peaks in the
% hist). just need more iters to get more 30/60 degree grids then will get
% better gridness scores?
% figure; hold on;
% iPlt=0;
% for iKvals = iToPlot
%     iPlt=iPlt+1;
%     subplot(6,5,iPlt)
%     % hist(g(:,2,iKvals),50);
%     scatter(g(:,1,iKvals),g(:,2,iKvals),'.')
% end




% % correlate over all nClus conditions?
% clus2use=2:28;
% [r p] = corr(reshape(squeeze(g(:,1,clus2use)),numel(g(:,1,clus2use)),1),reshape(tssekVals(clus2use,:)',numel(tssekVals(clus2use,:)),1),'type','spearman')
% figure; scatter(reshape(tssekVals(clus2use,:)',numel(tssekVals(clus2use,:)),1),reshape(squeeze(g(:,1,clus2use)),numel(g(:,1,clus2use)),1),'.')
% %multiply by nK 
% [r p] = corr(reshape(squeeze(g(:,1,clus2use)),numel(g(:,1,clus2use)),1),reshape(tssekVals(clus2use,:)'.*kVals(clus2use),numel(tssekVals(clus2use,:)),1),'type','spearman')
% figure; scatter(reshape(tssekVals(clus2use,:)'.*kVals(clus2use),numel(tssekVals(clus2use,:)),1),reshape(squeeze(g(:,1,clus2use)),numel(g(:,1,clus2use)),1),'.')

%% Density plots for top / bottom 3

iToPlot = 12;%[4, 22, 28];

gaussSmooth             = 1;
nSets                   = 6; %top and bottom 3 SSE
locRange                = [0, 49];
spacing                 = linspace(locRange(1),locRange(2),locRange(2)+1); 
densityPlotCentresBest  = zeros(length(spacing),length(spacing),nSets,nKvals);

%sort top/bottom 3
bestWorst = [1,2,3,nKmeans-3,nKmeans-2,nKmeans-1];
muBest = cell(1,4);
for iKvals = 1:nKvals
    for iterI=1:length(bestWorst)
        muBest{iKvals}(:,:,iterI) = muAllkVals{iKvals}(:,:,indSSE2(iKvals,:)==bestWorst(iterI));
    end
end

for iKvals = iToPlot%:nKvals %prob don't want to plot all at once...
    
    for iSet=1:nSets
        nK = kVals(iKvals);
        densityPlotClus = zeros(length(spacing),length(spacing),nK,nSets);
        figure; hold on;

        for iClus=1:nK
            clusTmp  = squeeze(round(muBest{iKvals}(iClus,:,iSet)))';
            for iTrlUpd=1:size(clusTmp,2)
                densityPlotClus(clusTmp(1,iTrlUpd),clusTmp(2,iTrlUpd),iClus,iSet) = densityPlotClus(clusTmp(1,iTrlUpd),clusTmp(2,iTrlUpd),iClus,iSet)+1;
            end
        end
        
        %make combined (grid cell) plot, smooth
        densityPlotCentresBest(:,:,iSet,iKvals) = sum(densityPlotClus(:,:,:,iSet),3);
        densityPlotCentresBestSm = imgaussfilt(densityPlotCentresBest(:,:,iSet,iKvals),gaussSmooth);
        
        subplot(1,3,1);
        imagesc(densityPlotCentresBestSm);
        aCorrMap=ndautoCORR(densityPlotCentresBestSm); %autocorrelogram
        subplot(1,3,2);
        imagesc(aCorrMap,[-.45 .45]);
        subplot(1,3,3);
        [g,gdataA] = gridSCORE(aCorrMap,'allen',1);
        % [g,gdataW] = gridSCORE(aCorrMap,'wills',1);
        
    end
end
%% plot cluster centres

saveplots=0;

iToPlot = 25:28;

colgrey = [.5, .5, .5];

for iKvals = iToPlot %nKvals - don't plot all if too many loaded up
    nK=kVals(iKvals);
    colors = distinguishable_colors(nK); %function for making distinguishable colors for plotting
    muAll = muBest{iKvals};
    
    %%%%%%%%%
    %edit this bit? currently just taking from above (computed muBest from
    %above)
    %%%%%%%%
    
    
    figure;
    for i = 1:6
        subplot(2,3,i); hold on;
        
        voronoi(muAll(:,1,i),muAll(:,2,i),'k');
        
        for iClus = 1:nK
            plot(muAll(iClus,1,i),muAll(iClus,2,i),'.','Color',colors(iClus,:),'MarkerSize',25); hold on; %plot cluster final point
        end
        xlim([locRange(1),locRange(2)]); ylim([locRange(1),locRange(2)]);
        hold on;
        if i==2
            title(sprintf('Lowest 3 and Highest 3 SSE cluster locations (k=%d)',nK));
        end
    end
    % if saveplots
    %     fname=[wd, sprintf('/kmeans_clus_k_%d_locs_top_bottom_3_datPts%dk',nK,nPoints/1000)];
    %     saveas(gcf,fname,'png');
    % end
end