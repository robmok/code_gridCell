%% plot centres 

nClus = size(muAvg,1);
colors = distinguishable_colors(nClus); %function for making distinguishable colors for plotting

iterI=1; 

for iSet=1:4 %size(muAvg,3) %plot - diff averaging over nTrials
    figure; hold on;
    scatter(muAvg(:,1,iSet,iterI),muAvg(:,2,iSet,iterI),20e+2,colors,'.');
    xlim(locRange); ylim(locRange);
%     voronoi(muAvg(:,1,iSet,iterI),muAvg(:,2,iSet,iterI),'k');
end
% 
figure; hold on;
scatter(muAll(:,1,end),muAll(:,2,end),20e+2,colors,'.');
xlim(locRange); ylim(locRange);
% voronoi(muAll(:,1,end),muAll(:,2,end),'k');
%% density map, autocorrelogram

iterI=1; 
gaussSmooth=1;

for iSet=1:size(muAvg,3) %plot - diff averaging over nTrials
    densityPlotSm = imgaussfilt(densityPlot(:,:,iSet,iterI),gaussSmooth);
    aCorrMap=ndautoCORR(densityPlotSm); %autocorrelogram

    figure; hold on;
    subplot(1,2,1); imagesc(densityPlotSm);
    subplot(1,2,2); imagesc(aCorrMap,[-.5 .5]);
end

%% density map, unsmoothed, smoothed

iterI=1; 
gaussSmooth=1;

for iSet=1:4%size(muAvg,3) 

   figure; hold on;
    subplot(1,2,1); imagesc(densityPlot(:,:,iSet,iterI));
    densityPlotSm = imgaussfilt(densityPlot(:,:,iSet,iterI),gaussSmooth);
    subplot(1,2,2); imagesc(densityPlotSm);
end

%%
gaussSmooth=1;

nSets=11;
spacing=linspace(locRange(1),locRange(2),locRange(2)+1); 
% densityPlotClus      = zeros(length(spacing),length(spacing),nClus,nSets,nIter);
densityPlotClus      = zeros(length(spacing),length(spacing),nClus,nSets);

for iSet=1:4%:nSets

for iClus=1:nClus
    clusTmp  = squeeze(round(muAvg(iClus,:,iSet)))';
    for iTrlUpd=1:size(clusTmp,2)
%         densityPlotClus(clusTmp(1,iTrlUpd),clusTmp(2,iTrlUpd),iClus,iSet,iterI) = densityPlotClus(clusTmp(1,iTrlUpd),clusTmp(2,iTrlUpd),iClus,iSet, iterI)+1;
        densityPlotClus(clusTmp(1,iTrlUpd),clusTmp(2,iTrlUpd),iClus,iSet) = densityPlotClus(clusTmp(1,iTrlUpd),clusTmp(2,iTrlUpd),iClus,iSet)+1;
    end
end

%make combined (grid cell) plot, smooth
densityPlotCentres(:,:,iSet) = sum(densityPlotClus(:,:,:,iSet),3);
densityPlotCentresSm = imgaussfilt(densityPlotCentres(:,:,iSet),gaussSmooth);

figure;
imagesc(densityPlotCentresSm);

aCorrMap=ndautoCORR(densityPlotCentresSm); %autocorrelogram
figure;
imagesc(aCorrMap,[-.45 .45]);

figure;
[g,gdataA] = gridSCORE(aCorrMap,'allen',1);
% [g,gdataW] = gridSCORE(aCorrMap,'wills',1);
    
end


%% mu - plot each trial; could compute gridness on each trial?

gaussSmooth=1;
% 
% spacing=linspace(locRange(1),locRange(2),locRange(2)+1); 
% densityPlotClus      = zeros(length(spacing),length(spacing),nClus,nSets);
% 
% 
% i = 30000; % could compute the gridness over trials..
% 
% for iClus=1:nClus
%     clusTmp  = squeeze(round(muAll(iClus,:,i)))';
%     for iTrlUpd=1:size(clusTmp,2)
%         densityPlotClus(clusTmp(1,iTrlUpd),clusTmp(2,iTrlUpd),iClus,iSet) = densityPlotClus(clusTmp(1,iTrlUpd),clusTmp(2,iTrlUpd),iClus,iSet)+1;
%     end
% end
% 
% %make combined (grid cell) plot, smooth
% densityPlotCentres(:,:,iSet) = sum(densityPlotClus(:,:,:,iSet),3);
% densityPlotCentresSm = imgaussfilt(densityPlotCentres(:,:,iSet),gaussSmooth);
% 
% figure;
% imagesc(densityPlotCentresSm);
% 
% aCorrMap=ndautoCORR(densityPlotCentresSm); %autocorrelogram
% figure;
% imagesc(aCorrMap,[-.45 .45]);
% 
% figure;
% [g,gdataA] = gridSCORE(aCorrMap,'allen',1);
% % [g,gdataW] = gridSCORE(aCorrMap,'wills',1);
    

%% plot over time (need muAll as output arg)

savePlots=0;

iterI = 1; %1 to 3, lowest sse, 4:6, highest sse

nClus = size(muAll,1);
colors = distinguishable_colors(nClus); %function for making distinguishable colors for plotting

% figure;   
figure('units','normalized','outerposition',[0 0 1 1]);
iPlot = 1; subplot(3,4,iPlot); hold on;%subplot
for iTrl = 1:nTrials
    if mod(iTrl,4000)==0
%         voronoi(muAll(:,1,iTrl,toPlot),muAll(:,2,iTrl,toPlot),'k'); %plot at the END before starting new subplot
        iPlot=iPlot+1;
        subplot(3,4,iPlot); hold on;
%         voronoi(muAll(:,1,iTrl,iterI),muAll(:,2,iTrl,iterI),'k'); %plot at the START showing the previous final positions
    end
    xlim(locRange); ylim(locRange);

    if mod(iTrl,500)==0 %plot centers after x trials
        for i=1:nClus
            plot(squeeze(muAll(i,1,iTrl,iterI)),squeeze(muAll(i,2,iTrl,iterI)),'.','Color',colors(i,:),'MarkerSize',20); hold on;
        end
        drawnow;
    end
end
% voronoi(muAll(:,1,iTrl,iterI),muAll(:,2,iTrl,iterI),'k'); %final one - if plotting at END above
% xlim([min(trials(:,1))-.1,max(trials(:,1))+.1]); ylim([min(trials(:,2))-.1,max(trials(:,2))+.1]);

if savePlots
    fname=[saveDir, sprintf('/covering_map_%d_clus_locs_plotOverTime_top_%d_trls_eps%d_alpha%d_%diters',nClus,nTrials,epsMuOrig100,alpha10,nIter)];
    if warpBox
        fname = [fname sprintf('_warpBox_resetEps%d',resetEps)];
    end
    saveas(gcf,fname,'png');
end

%% over time - one plot

nClus = size(muAll,1);

iterI = 1;
colors = distinguishable_colors(nClus); %function for making distinguishable colors for plotting

figure;
% figure('units','normalized','outerposition',[0 0 1 1]);
for iTrl = 1000:nTrials
    if mod(iTrl,500)==0
%         iPlot=iPlot+1;
%         voronoi(muAll(:,1,iTrl,iterI),muAll(:,2,iTrl,iterI),'k')
    end
%     xlim(locRange); ylim(locRange);

    if mod(iTrl,200)==0 %plot centers after x trials
        for i=1:nClus
            plot(squeeze(muAll(i,1,iTrl,iterI)),squeeze(muAll(i,2,iTrl,iterI)),'.','Color',colors(i,:),'MarkerSize',10); hold on; %make marker size bigger - larger/smoother firing field!
        end
        drawnow;
    end
end
