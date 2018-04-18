%% kmeans plot

% clear all

% wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
cd(wd);

codeDir = [wd '/code_gridCell'];
saveDir = [wd '/data_gridCell'];
addpath(codeDir); addpath(saveDir);
addpath(genpath([codeDir '/gridSCORE_packed']));

dat = 'circ'; %'square', 'circ'

kVals   = 3:30; 
% nKmeans = 1000;
nKmeans = 200;
nSims = 5; % nIters=nKmeans*nSims;
% nSims = 10;
nPoints = 10000; %5k, 10k
nKvals  = length(kVals);

nDataPtsTest = 10000; %5k, 10k, 50k
nXvalDataSets = 20;

%load in data
fname = [saveDir, sprintf('/kmeans_nK_%d-%d_%s_nPoints%d_%diters_%dsims_merged',kVals(1),kVals(end),dat,nPoints,nKmeans,nSims)];
load(fname);

% load in xVal
fname = [saveDir sprintf('/kmeans_xVal_nK_%d-%d_%s_%dtrainPts_%dtestPts_%diters_%sims_%ddatasets',kVals(1),kVals(end),dat,nPoints,nDataPtsTest,nKmeans,nSims,nXvalDataSets)];
load(fname);

tsseXval=xVal_results.tsseXval;
sseSprdSdXval=xVal_results.sseSprdSdXval;
sseSprdVarXval=xVal_results.sseSprdVarXval;

%comparing between nPoints
% figure;
% plot(squeeze((mean(mean(tsseXval,2),1))));
% hold on;
%% plot hist and density  plots

figsDir = [wd '/grid_figs'];

savePlots = 0;

gridMsr = 'a'; % 'a' or 'w' for allen or willis method

switch gridMsr
    case 'a'
        g = gA;
    case 'w'
        g = gW;
end

%fig specs
xTickLabs = num2cell(kVals);
fontSiz=13;

%plot univar scatters
figure;
dat1     = squeeze(g(:,1,:));
barpos  = .25:.5:.5*size(dat1,2);
colors  = distinguishable_colors(size(dat1,2));
colgrey = [.65, .65, .65];
mu      = nanmean(dat1,1);
% sm      = std(dat1)./sqrt(size(dat1,1));
% ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
% errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
% scatter(barpos,mu,750,colors,'x','LineWidth',1);
scatter(barpos,mu,50,colgrey,'filled','d');
xticklabels(xTickLabs);

xlim([barpos(1)-.5, barpos(end)+.5]);

if strcmp(gridMsr,'a')
    ylim([-.5,1.35]);
elseif strcmp(gridMsr,'w')
    ylim([-1.35,1.35]);
end

xlabel('Number of Clusters');
ylabel('Gridness score');
title(sprintf('kmeans %s',dat))
set(gca,'FontSize',fontSiz,'fontname','Arial')


fname = [figsDir sprintf('/kmeans_%s_gridness_univarScatters_nK%d-%d_%s',dat,kVals(1),kVals(end),gridMsr)];
if savePlots
   set(gcf,'Renderer','painters');
   print(gcf,'-depsc2',fname)
   saveas(gcf,fname,'png');
    
end
    

%plot some grids

% densityPlotCentres(:,:,kMeansIter,iKvals) = sum(densityPlotClus,3);
% densityPlotCentresSm = imgaussfilt(densityPlotCentres(:,:,kMeansIter,iKvals),gaussSmooth);
% aCorrMap=ndautoCORR(densityPlotCentresSm); %autocorrelogram
% [~,gdataA] = gridSCORE(aCorrMap,'allen',0);
% [~,gdataW] = gridSCORE(aCorrMap,'wills',0);

%% Crossvalidation - SSE on the test dataset
% note: tsseXval is already ordered from low to high SSE on train set

iToPlot = 1:28;
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
% iToPlot = [1:28];
% figure; hold on;
% iPlt=0;
% for iKvals = iToPlot
%     iPlt=iPlt+1;
%     subplot(6,5,iPlt)
%     [r p] = corr(g(:,1,iKvals),err(iKvals,:)','type','spearman');
%     r1(iKvals)=r;
%     p1(iKvals)=p;
%     p1corr(iKvals)=p*nKvals;
%     if p1corr(iKvals)>1, p1corr(iKvals)=1; end
%     scatter(err(iKvals,:)',g(:,1,iKvals),'.');
%     title(sprintf('nK=%d, r = %.2f, p (crcted) = %.3f',kVals(iKvals),r,p1corr(iKvals)))
% end
% [r1; p1corr]'

% corr of gridness and SSE: cross-validated data

iToPlot = [1:24]; % excluding more clusters (3:26 here)
figure; hold on;
iPlt=0;
for iKvals = iToPlot
    iPlt=iPlt+1;
    subplot(6,4,iPlt)
    [r p] = corr(g(:,1,iKvals),nanmean(xValErr(:,:,iKvals),2),'type','spearman');
    rx1(iKvals)=r;
    px1(iKvals)=p;
    px1corr(iKvals)=p*(nKvals-4); %-4 if only reporting 26
    if px1corr(iKvals)>1, px1corr(iKvals)=1; end
    h = scatter(nanmean(xValErr(:,:,iKvals),2),g(:,1,iKvals),'.'); hold on;
    [b, d, s] = glmfit(nanmean(xValErr(:,:,iKvals),2),g(:,1,iKvals));
    xMinMax=[min(nanmean(xValErr(:,:,iKvals),2)), max(nanmean(xValErr(:,:,iKvals),2))];
    yhat = b(1)+b(2)*xMinMax;
    hLine=plot(xMinMax,yhat,'Color',colgrey+.1','LineWidth',.5);
    ylim([-.5, 1.5]);
    if px1corr(iKvals) <0.0001
        plabels = sprintf('rho = %.2f\npCorr < 0.001',r);
    else
        plabels = sprintf('rho = %.2f\npCorr = %.3f',r,px1corr(iKvals));
    end
    hleg = legend(h,plabels,'Location','NorthEast');
    title(sprintf('%d clusters',kVals(iKvals),r,px1corr(iKvals)))
end

[rx1(iToPlot); px1corr(iToPlot)]'


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
%% Density plots for top / bottom 3

% making nicer imagesc images - circ: grey out the bits that are not a
% circle. find the points and those outside of them should be white
% ([1,1,1])

% circ - might want to cut out the circ before doing autoCorrelogram?

%think about autocorr plots for gridness

iToPlot = 10;%[4, 22, 28];

gaussSmooth             = 1;
nSets                   = 6; %top and bottom 3 SSE
locRange                = [0, 49];
spacing                 = linspace(locRange(1),locRange(2),locRange(2)+1); 
densityPlotCentresBest  = zeros(length(spacing),length(spacing),nSets,nKvals);

%re-sort sse after combining tsseXval values-atm just one nK value at a time
[y, indSSEall] = sort(mean(tsseXval(:,:,iKvals),2));
[y, indSSEall2] = sort(indSSEall);

bestWorst = [1,2,3,nKmeans*nSims-2,nKmeans*nSims-1,nKmeans*nSims]; 
muBest = cell(1,4);
for iKvals = 1:nKvals
    for iterI=1:length(bestWorst)
        muBest{iKvals}(:,:,iterI) = muAllkVals{iKvals}(:,:,indSSEall2==bestWorst(iterI));
    end
end

for iKvals = iToPlot%:nKvals %prob don't want to plot all at once...
    
    for iSet=1%:nSets
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
        imagesc(aCorrMap,[-.3 .3]);
        subplot(1,3,3);
        [g,gdataA] = gridSCORE(aCorrMap,'allen',1);
        % [g,gdataW] = gridSCORE(aCorrMap,'wills',1);
        
    end
end
%% plot cluster centres

% for kmeans, thinking to use this; nice colours. but for circ - need draw
% a circ to show its a circ box, trapz for trapz, etc.



% plot best/worse 3. or best and worse for each nClus cond (since so many);
% or just plot a selection


saveplots=0;

locRange = [0, 49];

iToPlot = [5, 8, 15, 17, 19, 24, 26, 28];

colgrey = [.5, .5, .5];

for iKvals = iToPlot %nKvals - don't plot all if too many loaded up
    nK=kVals(iKvals);
    colors = distinguishable_colors(nK); %function for making distinguishable colors for plotting
    muAll = muAllkVals{iKvals};
    
    %%%%%%%%%
    %edit this bit? currently just taking from above (computed muBest from
    %above)
    %%%%%%%%
    
    
    figure;
    for i = 1:6 %currently not best/worst 6
        subplot(2,3,i); hold on;
        
        voronoi(muAll(:,1,i),muAll(:,2,i),'k');
        
        for iClus = 1:nK
            plot(muAll(iClus,1,i),muAll(iClus,2,i),'.','Color',colors(iClus,:),'MarkerSize',25); hold on; %plot cluster final point
        end
        xlim([locRange(1),locRange(2)]); ylim([locRange(1),locRange(2)]);
        hold on;
%         if i==2
%             title(sprintf('Lowest 3 and Highest 3 SSE cluster locations (k=%d)',nK));
%         end
    end
    % if saveplots
    %     fname=[wd, sprintf('/kmeans_clus_k_%d_locs_top_bottom_3_datPts%dk',nK,nPoints/1000)];
    %     saveas(gcf,fname,'png');
    % end
end