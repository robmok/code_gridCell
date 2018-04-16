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

%% gridness, autocorrelogram

gaussSmooth=1;
for iSet=1:20+2 %size(densityPlotActTNorm,3)
    
for iterI = 1%:3%:10 
    
% densityPlotCentresSm = imgaussfilt(densityPlot(:,:,iSet,iterI),gaussSmooth);
% densityPlotCentresSm = imgaussfilt(densityPlotAct(:,:,iSet,iterI),gaussSmooth);
densityPlotCentresSm = imgaussfilt(densityPlotActNorm(:,:,iSet,iterI),gaussSmooth);

figure; hold on;
subplot(1,2,1)
imagesc(densityPlotCentresSm);
% clos
aCorrMap=ndautoCORR(densityPlotCentresSm); %autocorrelogram
% subplot(1,3,2)
% imagesc(aCorrMap,[-.45 .45]);

% figure;
subplot(1,2,2)

[g,gdataA] = gridSCORE(aCorrMap,'allen',1);
% [g,gdataW] = gridSCORE(aCorrMap,'wills',1);
end
end

%% gridness - left vs right half

spacing=linspace(locRange(1),locRange(2),locRange(2)+1)+1;

%spacing split into 2 equal areas - trapzKrupic only for now]

%krupic
spacingTrapz = spacing(14:37);
%boxSize=2
spacingTrapz = spacing(14:37+length(14:37));

%trapzScaled1
% spacingTrapz = spacing(10:41);

a = length(spacingTrapz);
b = length(spacing);
h = length(spacing);
halfArea = (((a+b)/2)*h)/2;
c = sqrt(((a^2)+(b^2))/2);
% hLeft  = round(halfArea/((b+c)/2)); %bigger side
hLeft  = floor(halfArea/((b+c)/2)); %bigger side
hRight = ceil(halfArea/((b+c)/2))+1; %smaller size - to equalize area more

gaussSmooth=1;
for iSet=21:22
for iterI = 1%:10
    
    densityPlotCentresSm = imgaussfilt(densityPlotActNorm(:,:,iSet,iterI),gaussSmooth);
    
    figure; hold on;
    subplot(1,3,1)
    imagesc(densityPlotCentresSm);
    
    %left half of box
    subplot(1,3,2)
%     aCorrMap = ndautoCORR(densityPlotCentresSm(:,spacing(1):spacing(ceil(length(spacing)/2))));
    aCorrMap = ndautoCORR(densityPlotCentresSm(:,1:hLeft));
    [g,gdataA] = gridSCORE(aCorrMap,'allen',1);
    % [g,gdataW] = gridSCORE(aCorrMap,'wills',1);
    
    %right half of box
    subplot(1,3,3)
%     aCorrMap = ndautoCORR(densityPlotCentresSm(:,spacing(ceil(length(spacing)/2))+1:spacing(end)));
    aCorrMap = ndautoCORR(densityPlotCentresSm(:,hRight:end));
    [g,gdataA] = gridSCORE(aCorrMap,'allen',1);
    % [g,gdataW] = gridSCORE(aCorrMap,'wills',1);
end
end

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

%if plotAgent, need trials as an output argument
plotAgent = 1;


nClus = size(muAll,1);

iterI = 1;
colors = distinguishable_colors(nClus); %function for making distinguishable colors for plotting
colAgent = [.75, .75, .75];

% colAgent gets darker over time? also clear the lines?
colAgentCh = fliplr(linspace(0,.9,nTrials));

%if plotting agent over batches - need make cluster positions same over
%trials in a batch
% if plotAgent
muTrls = nan(nClus,2,nTrials);
for iBatch = 1:nBatches
    muTrls(:,1,batchSize*(iBatch-1)+1:batchSize*(iBatch-1)+batchSize)=repmat(muAll(:,1,iBatch),1,batchSize);
    muTrls(:,2,batchSize*(iBatch-1)+1:batchSize*(iBatch-1)+batchSize)=repmat(muAll(:,2,iBatch),1,batchSize);
end
% end

figure;
% figure('units','normalized','outerposition',[0 0 1 1]);
for iTrl = 100000:nTrials
    if mod(iTrl,10)==0 %plot centers after x trials
        %agent
        if plotAgent
%         scatter(trials(iTrl,1),trials(iTrl,2),1000,colAgent,'.');
%         plot(trials(1:iTrl,1),trials(1:iTrl,2),'Color',colAgent); 
%         scatter1 = scatter(trials(iTrl,1),trials(iTrl,2),'MarkerFaceColor',colAgent,'MarkerEdgeColor',colAgent);alpha(scatter1,.2)
%         plot(trials(1:iTrl,1),trials(1:iTrl,2),'Color',repmat(colAgentCh(iTrl),1,3)); % make it darker over time 

        %have to start from 1000 trials min
        plotTrls=iTrl-999:iTrl;
        h1(iTrl)=plot(trials(plotTrls,1),trials(plotTrls,2),'Color',repmat(colAgentCh(iTrl),1,3));
        if mod(iTrl,5000)==0
            delete(h1);
        end
        end
        
        %clusters
%         scatter(squeeze(muAll(:,1,iTrl,iterI)),squeeze(muAll(:,2,iTrl,iterI)),200,colors,'.'); hold on;
        scatter(squeeze(muTrls(:,1,iTrl)),squeeze(muTrls(:,2,iTrl)),2000,colors,'.'); hold on;
        drawnow;
    end
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
    