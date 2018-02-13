%kmeans plot

%load....

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
    k = kVals(iK);
    
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

% gridness corr with tsse
% figure;
% scatter(g(:,1,iK),tssekVals(iK,:))
% [r p] = corr(g(:,1),tssekVals(iK,:)')
% [r p] = corr(g(:,1),tssekVals(iK,:)','type','spearman')


%plot some grids


densityPlotCentres(:,:,kMeansIter,iKvals) = sum(densityPlotClus,3);
densityPlotCentresSm = imgaussfilt(densityPlotCentres(:,:,kMeansIter,iKvals),gaussSmooth);
aCorrMap=ndautoCORR(densityPlotCentresSm); %autocorrelogram
[g,gdataA] = gridSCORE(aCorrMap,'allen',0);
[g,gdataW] = gridSCORE(aCorrMap,'wills',0);







%% plots

saveplots=0;

colgrey = [.5, .5, .5];

% iToPlot=[1,2,3,nKmeans-2,nKmeans-1,nKmeans];

iK = 3;
nK=kVals(iK);
colors = distinguishable_colors(nK); %function for making distinguishable colors for plotting
muAll = muAllkVals{iK};

figure;
for i = 1:6
    subplot(2,3,i); hold on;

    voronoi(muAll(:,1,i),muAll(:,2,i),'k');

    for iClus = 1:nK
        plot(muAll(iClus,1,i),muAll(iClus,2,i),'.','Color',colors(iClus,:),'MarkerSize',25); hold on; %plot cluster final point
%         plot(muAll(iClus,1,(indSSE2==iToPlot(i))),muAll(iClus,2,(indSSE2==iToPlot(i))),'.','Color',colors(iClus,:),'MarkerSize',25); hold on; %plot cluster final point
    end
    xlim([locRange(1),locRange(2)]); ylim([locRange(1),locRange(2)]);
    hold on;
    if i==2
        title(sprintf('Lowest 3 and Highest 3 SSE cluster locations (k=%d)',nK));
    end
end
if saveplots
    fname=[wd, sprintf('/kmeans_clus_k_%d_locs_top_bottom_3_datPts%dk',nK,nPoints/1000)];
    saveas(gcf,fname,'png');
end

%%
% SSE on the test dataset, ranked from the training data set
% note: tsseXval is already ordered from low to high SSE on train set

% SSE sorted on training data
% figure; plot((tsse(indSSE))); title(sprintf('SSE over diff initializations sorted (training data) (k=%d)',nK));
% if saveplots
%     fname=[wd, sprintf('/kmeans_clus_k_%d_sse_train_datPts%dk',nK,nPoints/1000)];
% %     fname = [fname '_limDatPtsTo85'];
%     saveas(gcf,fname,'png');
% end

%on test data
for iKvals = 1:nKvals
figure;
plot(tsseXval(:,:,iKvals)','Color',[.75 .75 .75]); hold on;
plot(tsseXval(1,:,iKvals)','Color',[.2 .2 .2]);
title(sprintf('SSE on Test data sorted by training data (k=%d) (nTest=%d)',nK,nXvalDataSets));
% if saveplots
%     fname=[wd, sprintf('/kmeans_clus_k_%d_sse_test_sortedByTrain_datPts%dk',nK,nPoints/1000)];
%     saveas(gcf,fname,'png');
% end
end


% cross-validation - plot SSE over datasets
for iKvals = 1:nKvals
    figure; hold on;
    plot(tsseXval(:,:,iKvals)','Color',[.75 .75 .75]);
    plot(tsseXval(min(mean(tsseXval(:,:,iKvals),2))==mean(tsseXval(:,:,iKvals),2),:,iKvals),'Color',[.5 .5 .5],'LineWidth',2); %lowest mean SSE
    plot(tsseXval(1,:,iKvals)','Color',[.2 .2 .2],'LineWidth',2); %lowest mean SSE from training data
    title(sprintf('SSE over all Test data sets (k=%d) (nTest=%d)(datPts=%dk)',nK,nXvalDataSets,nPoints/1000));
    % if saveplots
    %     fname=[wd, sprintf('/kmeans_clus_k_%d_sse_test_overDatasets_datPts%dk',nK,nPoints/1000)];
    % %     fname = [fname '_limDatPtsTo85'];
    %     saveas(gcf,fname,'png');
    % end
end

% figure;hist([tsseHex;tsseXval(min(mean(tsseXval,2))==mean(tsseXval,2),:)]',20);


%scatters - gridness vs original data
for iKvals = 1:nKvals
figure;
[r p] = corr(g(:,1,iKvals),(tssekVals(iKvals,:)'),'type','pearson');
% [r p] = corr(g(:,1,iKvals),(tssekVals(iKvals,:)'),'type','spearman');
scatter(g(:,1,iKvals),(tssekVals(iKvals,:)'),'.');
title(sprintf('r = %.2f, p = %.5f',r,p))
end


%scatters - gridness vs cross-validated data

%%%%%
% - this doesn't work yet; tsseXval is sorted by sse in training set; need
% reverse the ind from indSSE2, or sth.. save the indSSE1?

%%%%%%% - to CHECK %%%%%
% - saved indSSE1 now; put into as an index below, check it works
%%%%%

for iKvals = 1:nKvals
figure;
[r p] = corr(g(:,1,iKvals),mean(tsseXval(indSSE1,:,iKvals),2),'type','pearson');
% [r p] = corr(g(:,1,iKvals),mean(tsseXval(indSSE1,:,iKvals),2),'type','spearman');
scatter(g(:,1,iKvals),mean(tsseXval(indSSE1,:,iKvals),2),'.');
title(sprintf('r = %.2f, p = %.5f',r,p))
end
%%
gaussSmooth=1;

nSets=6; %top and bottom 3 SSE
spacing=linspace(locRange(1),locRange(2),locRange(2)+1); 
densityPlotClus      = zeros(length(spacing),length(spacing),nK,nSets);

for iSet=1:nSets
figure; hold on;
for iClus=1:nK
    clusTmp  = squeeze(round(muBest(iClus,:,iSet)))';
    for iTrlUpd=1:size(clusTmp,2)
        densityPlotClus(clusTmp(1,iTrlUpd),clusTmp(2,iTrlUpd),iClus,iSet) = densityPlotClus(clusTmp(1,iTrlUpd),clusTmp(2,iTrlUpd),iClus,iSet)+1;
    end
end

%make combined (grid cell) plot, smooth
densityPlotCentres(:,:,iSet) = sum(densityPlotClus(:,:,:,iSet),3);
densityPlotCentresSm = imgaussfilt(densityPlotCentres(:,:,iSet),gaussSmooth);

subplot(1,3,1);
imagesc(densityPlotCentresSm);


aCorrMap=ndautoCORR(densityPlotCentresSm); %autocorrelogram
subplot(1,3,2);
imagesc(aCorrMap,[-.45 .45]);

subplot(1,3,3);
[g,gdataA] = gridSCORE(aCorrMap,'allen',1);
% [g,gdataW] = gridSCORE(aCorrMap,'wills',1);
    
end

