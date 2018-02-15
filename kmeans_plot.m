%% kmeans plot

clear all

% wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
cd(wd);

codeDir = [wd '/code_gridCell'];
saveDir = [wd '/data_gridCell'];
addpath(codeDir); addpath(saveDir);


dat = 'randUnique'; %'rand' points in a box, randUnique (all unique points in box), or 'cat'

kVals = 3:30; 
nKmeans = 1000;
nPoints = 10000; 
nKvals = length(kVals);

nXvalDataSets = 20;

%load in data
switch dat
    case 'randUnique'
        fname1 = [saveDir, sprintf('/kmeans_nK_%d-%d_uniquePts_%diters',kVals(1),kVals(end),nKmeans)];
        fname2 = [saveDir, sprintf('/kmeans_nK_%d-%d_uniquePts_xVal_%ddatasets_%diters',kVals(1),kVals(end),nXvalDataSets,nKmeans)];
    case 'rand'
        fname1 = [saveDir, sprintf('/kmeans_nK_%d-%d_randPts_%diters',kVals(1),kVals(end),nKmeans)];
        fname2 = [saveDir, sprintf('/kmeans_nK_%d-%d_randPts_xVal_%ddatasets_%diters',kVals(1),kVals(end),nXvalDataSets,nKmeans)];
end
load(fname1);
load(fname2);

%merge
% fname=[saveDir '/kmeans_nK_3-17_randPts_1000iters'];
% d1=load(fname);
% fname=[saveDir '/kmeans_nK_18-25_randPts_1000iters'];
% d2=load(fname);
% fname=[saveDir '/kmeans_nK_26-30_randPts_1000iters'];
% d3=load(fname);
% densityPlotCentres  = cat(4,d1.densityPlotCentres,d2.densityPlotCentres,d3.densityPlotCentres);
% gA                  = cat(3,d1.gA,d2.gA,d3.gA);
% gW                  = cat(3,d1.gW,d2.gW,d3.gW);
% indSSE1             = cat(1,d1.indSSE1,d2.indSSE1,d3.indSSE1);
% indSSE2             = cat(1,d1.indSSE2,d2.indSSE2,d3.indSSE2);
% kVals               = cat(2,d1.kVals,d2.kVals,d3.kVals);
% muAllkVals          = cat(2,d1.muAllkVals,d2.muAllkVals,d3.muAllkVals);
% tssekVals           = cat(1,d1.tssekVals,d2.tssekVals,d3.tssekVals);
% if 0
%     fname=[saveDir '/kmeans_nK_3-30_randPts_1000iters'];
%     save(fname,'muAllkVals','tssekVals', 'gA','gW','densityPlotCentres','indSSE1','indSSE2','kVals')
% end
% 
% fname=[saveDir '/kmeans_nK_3-17_randPts_xVal_20datasets_1000iters'];
% d1=load(fname);
% fname=[saveDir '/kmeans_nK_18-25_randPts_xVal_20datasets_1000iters'];
% d2=load(fname);
% fname=[saveDir '/kmeans_nK_26-30_randPts_xVal_20datasets_1000iters'];
% d3=load(fname);
% tsseXval = cat(3,d1.tsseXval,d2.tsseXval,d3.tsseXval);
% kVals    = cat(2,d1.kVals,d2.kVals,d3.kVals);
% if 0
%    fname=[saveDir '/kmeans_nK_3-30_randPts_xVal_20datasets_1000iters'];
%    save(fname,'tsseXval','kVals')
% end

%% plot hist and density  plots

gridMsr = 'a'; % 'a' or 'w' for allen or willis method

switch gridMsr
    case 'a'
        g = gA;
    case 'w'
        g = gW;
end

%plot hist
figure; pltCount=1;
for iK = 1:nKvals    
    subplot(3,2,pltCount)
    hist(squeeze(g(:,1,iK)),50);
    xlim([-.5,1.25]);
    
    pltCount = pltCount+1;
    if mod(iK,6)==0
    pltCount=1; figure;
    end
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
scatter(barpos,mu,750,colors,'x');
xlim([barpos(1)-.5, barpos(end)+.5]);
ylim([-.5,1.25]);

%plot some grids

% densityPlotCentres(:,:,kMeansIter,iKvals) = sum(densityPlotClus,3);
% densityPlotCentresSm = imgaussfilt(densityPlotCentres(:,:,kMeansIter,iKvals),gaussSmooth);
% aCorrMap=ndautoCORR(densityPlotCentresSm); %autocorrelogram
% [~,gdataA] = gridSCORE(aCorrMap,'allen',0);
% [~,gdataW] = gridSCORE(aCorrMap,'wills',0);

%% Crossvalidation - SSE on the test dataset, corr with gridness
% note: tsseXval is already ordered from low to high SSE on train set

iToPlot = 25:38; %1:28 - but plotting all will create too many plots atm - make subplots




% SSE on test data,  showing lowest SSE on training data
for iKvals = iToPlot %1:nKvals
    nK = kVals(iKvals);
    figure;
    plot(tsseXval(:,:,iKvals)','Color',[.75 .75 .75]); hold on;
    plot(tsseXval(1,:,iKvals)','Color',[.2 .2 .2]);
    title(sprintf('SSE on Test data showing lowest SSE on training data (k=%d) (nTest=%d)',nK,nXvalDataSets));
    % if saveplots
    %     fname=[wd, sprintf('/kmeans_clus_k_%d_sse_test_sortedByTrain_datPts%dk',nK,nPoints/1000)];
    %     saveas(gcf,fname,'png');
    % end
end

% cross-validation - plot SSE over datasets
for iKvals = iToPlot
    nK = kVals(iKvals);
    figure; hold on;
    plot(tsseXval(:,:,iKvals)','Color',[.75 .75 .75]);
    plot(tsseXval(min(nanmean(tsseXval(:,:,iKvals),2))==nanmean(tsseXval(:,:,iKvals),2),:,iKvals),'Color',[.5 .5 .5],'LineWidth',2); %lowest mean SSE
    plot(tsseXval(1,:,iKvals)','Color',[.2 .2 .2],'LineWidth',2); %lowest mean SSE from training data
    title(sprintf('SSE over all Test data sets,  showing lowest SSE on training data (k=%d) (nTest=%d)(datPts=%dk)',nK,nXvalDataSets,nPoints/1000));
    % if saveplots
    %     fname=[wd, sprintf('/kmeans_clus_k_%d_sse_test_overDatasets_datPts%dk',nK,nPoints/1000)];
    % %     fname = [fname '_limDatPtsTo85'];
    %     saveas(gcf,fname,'png');
    % end
end

%scatter plots - corr of gridness and SSE: original data
for iKvals = iToPlot
figure;
[r, p] = corr(g(:,1,iKvals),(tssekVals(iKvals,:)'),'type','pearson');
% [r p] = corr(g(:,1,iKvals),(tssekVals(iKvals,:)'),'type','spearman');
scatter(g(:,1,iKvals),(tssekVals(iKvals,:)'),'.');
title(sprintf('r = %.2f, p = %.5f',r,p))
end

% corr of gridness and SSE: cross-validated data
for iKvals = iToPlot
figure;
[r, p] = corr(g(:,1,iKvals),nanmean(tsseXval(indSSE1(iKvals,:),:,iKvals),2),'type','pearson');
% [r p] = corr(g(:,1,iKvals),mean(tsseXval(indSSE1(iKvals,:),:,iKvals),2),'type','spearman');
scatter(g(:,1,iKvals),nanmean(tsseXval(indSSE1(iKvals,:),:,iKvals),2),'.');
title(sprintf('r = %.2f, p = %.5f',r,p))
end





% are 30/60 degree orientations --> better gridness? nk=3 suggests so.
% others not sure; most orientations at 0, 90 and 40/50.. (3 peaks in the
% hist). just need more iters to get more 30/60 degree grids then will get
% better gridness scores?
for iKvals = iToPlot
figure;

% hist(g(:,2,iKvals),50);
% scatter(g(:,1,iKvals),g(:,2,iKvals))

end
%% Density plots for top / bottom 3

iToPlot = 15:18;

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