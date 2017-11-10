%% Plot (with/without load in data)

% save plots?
savePlots=0;

%set to load in data (note divide by value in variable)
if 1
%     epsMuOrig100 = 75;
    epsMuOrig1000 = 75;
%     deltaEpsMu100 = 96; %96, 98
    nClus = 40;
    nTrials = 30000;
    nIter=5;
%     nIter=20; %for warpBox, for now
    warpBox=0;
    alpha10=8;
%     fname=[saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_deltaMu%d_%diters',nClus,nTrials,epsMuOrig10,deltaEpsMu100,nIter)];
%     fname=[saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_alpha%d_%diters',nClus,nTrials,epsMuOrig1000 ,alpha10,nIter)];
    fname=[saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_alpha%d_cVal%d_%diters',nClus,nTrials,epsMuOrig1000,alpha10,c10k,nIter)];
    
    if warpBox,
        fname = [fname '_warpBox'];
    end
    load(fname);
end



colors = distinguishable_colors(nClus); %function for making distinguishable colors for plotting

figure; plot(tsseTrls);
title(sprintf('SSE for each iteration (nClus=%d)',nClus));
if savePlots,
%     fname=[saveDir, sprintf('/covering_map_%d_clus_tsseTrls_%d_trls_eps%d_alpha%d_%diters',nClus,nTrials,epsMuOrig100,alpha10,nIter)];
    fname=[saveDir, sprintf('/covering_map_%d_clus_tsseTrls_%d_trls_eps%d_alpha%d_cVal%d_%diters',nClus,nTrials,epsMuOrig100,alpha10,c10k,nIter)];
    if warpBox,
        fname = [fname '_warpBox'];
    end
    saveas(gcf,fname,'png');
end

% plot change in learning rate for each cluster:
% figure;
% plot(epsMuAll(:,1),'.');

% plot lowest and highest SSE
iToPlot=size(muEnd,3);
figure;
for i = 1:iToPlot
    subplot(2,3,i); hold on;
%     plot(dataPtsTest(:,1),dataPtsTest(:,2),'.','Color',colgrey,'MarkerSize',2); hold on;
    voronoi(squeeze(muEnd(:,1,i)),squeeze(muEnd(:,2,i)),'k')
    for iClus = 1:nClus
        plot(mean(muEnd(iClus,1,(i))),mean(muEnd(iClus,2,(i))),'.','Color',colors(iClus,:),'MarkerSize',25); hold on; %plot cluster final point
    end
    if warpBox,
        xlim([-1.1,2.1]); ylim([-1.1,1.1]);
    else
        xlim([-1.1,1.1]); ylim([-1.1,1.1]);
    end 
    hold on;
    if i==2,
        title(sprintf('Lowest 3 and Highest 3 SSE cluster locations (nClus=%d)',nClus));
    end
end
if savePlots,
%     fname=[saveDir, sprintf('/covering_map_%d_clus_locs_top_bottom_3_%d_trls_eps%d_alpha%d_%diters',nClus,nTrials,epsMuOrig100,alpha10,nIter)];
    fname=[saveDir, sprintf('/covering_map_%d_clus_locs_top_bottom_3_%d_trls_eps%d_alpha%d_cVal%d_%diters',nClus,nTrials,epsMuOrig100,alpha10,c10k,nIter)];
    if warpBox,
        fname = [fname '_warpBox'];
    end
    saveas(gcf,fname,'png');
end

%% plot cluster centers over time - subplots

savePlots=0;

toPlot = 1; %1 to 3, lowest sse, 4:6, highest sse

colors = distinguishable_colors(nClus); %function for making distinguishable colors for plotting

% figure;   
figure('units','normalized','outerposition',[0 0 1 1]);
iPlot = 1; subplot(3,4,iPlot); hold on;%subplot
for iTrl = 1:nTrials
    if mod(iTrl,1000)==0
%         voronoi(muAll(:,1,iTrl,toPlot),muAll(:,2,iTrl,toPlot),'k'); %plot at the END before starting new subplot
        iPlot=iPlot+1;
        subplot(3,4,iPlot); hold on;
        voronoi(muAll(:,1,iTrl,toPlot),muAll(:,2,iTrl,toPlot),'k'); %plot at the START showing the previous final positions

    end
    if warpBox,
        xlim([-1.1,2.1]); ylim([-1.1,1.1]);
    else
        xlim([-1.1,1.1]); ylim([-1.1,1.1]);
    end
    if mod(iTrl,50)==0, %plot centers after x trials
        for i=1:nClus
            plot(squeeze(muAll(i,1,iTrl,toPlot)),squeeze(muAll(i,2,iTrl,toPlot)),'.','Color',colors(i,:),'MarkerSize',20); hold on;
        end
        drawnow;
    end
end
voronoi(muAll(:,1,iTrl,toPlot),muAll(:,2,iTrl,toPlot),'k'); %final one - if plotting at END above
% xlim([min(trials(:,1))-.1,max(trials(:,1))+.1]); ylim([min(trials(:,2))-.1,max(trials(:,2))+.1]);

if savePlots,
    fname=[saveDir, sprintf('/covering_map_%d_clus_locs_plotOverTime_top_%d_trls_eps%d_alpha%d_%diters',nClus,nTrials,epsMuOrig100,alpha10,nIter)];
    if warpBox
        fname = [fname sprintf('_warpBox_resetEps%d',resetEps)];
    end
    saveas(gcf,fname,'png');
end

%% over time - one plot

toPlot = 1;
colors = distinguishable_colors(nClus); %function for making distinguishable colors for plotting

figure;
% figure('units','normalized','outerposition',[0 0 1 1]);
for iTrl = 1:nTrials
    if mod(iTrl,250)==0
%         iPlot=iPlot+1;
%         voronoi(muAll(:,1,iTrl,toPlot),muAll(:,2,iTrl,toPlot),'k')
    end
    if warpBox,
        xlim([-1.1,2.1]); ylim([-1.1,1.1]);
    else
        xlim([-1.1,1.1]); ylim([-1.1,1.1]);
    end
    if mod(iTrl,50)==0, %plot centers after x trials
        for i=1:nClus
            plot(squeeze(muAll(i,1,iTrl,toPlot)),squeeze(muAll(i,2,iTrl,toPlot)),'.','Color',colors(i,:),'MarkerSize',10); hold on; %make marker size bigger - larger/smoother firing field!
        end
        drawnow;
    end
end

%% plot final centres 
toPlot=1;
colors = distinguishable_colors(nClus); %function for making distinguishable colors for plotting

% for i = 1:nIter
%     figure;
% %     plot(trials(:,1),trials(:,2),'Color',colgrey); hold on;
%     for iClus = 1:nClus
%         plot(squeeze(muEnd(iClus,1,i)),squeeze(muEnd(iClus,2,i)),'.','Color',colors(iClus,:),'MarkerSize',25); hold on; %plot cluster final point
%     end
%     if warpBox,
%         xlim([-1.1,2.1]); ylim([-1.1,1.1]);
%     else
%         xlim([-1.1,1.1]); ylim([-1.1,1.1]);
%     end
% end
% voronoi(squeeze(muEnd(:,1,i)),squeeze(muEnd(:,2,i)))

%plot over average of final trials
nTrlsToPlot = [3000, 2000, 1000, 500, 200, 100];

nTrlsToPlot = 1000;

for iToPlot=1:length(nTrlsToPlot)
    for i = 1%:nIter
        figure;
        for iClus = 1:nClus
            plot(mean(squeeze(muAll(iClus,1,nTrials-nTrlsToPlot(iToPlot):end,toPlot))),squeeze(mean(muAll(iClus,2,nTrials-nTrlsToPlot(iToPlot):end,toPlot))),'.','Color',colors(i,:),'MarkerSize',20); hold on;
            
        end
        if warpBox,
            xlim([-1.1,2.1]); ylim([-1.1,1.1]);
        else
            xlim([-1.1,1.1]); ylim([-1.1,1.1]);
        end
    end
    voronoi(squeeze(mean(muAll(:,1,nTrials-nTrlsToPlot(iToPlot):end,toPlot),3)),squeeze(mean(muAll(:,2,nTrials-nTrlsToPlot(iToPlot):end,toPlot),3)),'k');
end

%%

avgTrlSSE = mean(sseTrl,1);

for iClus=1:nClus,
    devAvgSSE(iClus,:)=(sseTrl(iClus,:)-avgTrlSSE);
end

stdAcrossClus=std(devAvgSSE);
varAcrossClus=var(devAvgSSE);

figure; plot(devAvgSSE');
figure; plot(varAcrossClus)
figure; plot(stdAcrossClus)


