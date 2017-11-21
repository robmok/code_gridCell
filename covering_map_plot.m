%% Plot (with/without load in data)
clear all;

wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
cd(wd);
saveDir = [wd '/data_gridCell'];

addpath(genpath([wd '/gridSCORE_packed']));

nSteps = 500;
locRange = [0, nSteps-1]; 

% save plots?
% savePlots=0;

%set to load in data (note divide by value in variable)
if 1
    epsMuOrig1000 = 75;
    nClus = 40;
    nTrials = 30000;
    nIter=3;
    warpBox=0;
    alpha10=8; 

    
    %atm have only 1 c value run for alpha=2/8, but all stypes
    %also have 2 and 3 for alpha=5, can check those out too
    stochasticType=1;
    cVals = round([2/nTrials, 3/nTrials, 5/nTrials, 10/nTrials].*1000000);
    c=cVals(1);
    
%     fname=[saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_deltaMu%d_%diters',nClus,nTrials,epsMuOrig10,deltaEpsMu100,nIter)];
%     fname=[saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_alpha%d_%diters',nClus,nTrials,epsMuOrig1000 ,alpha10,nIter)];
    fname = [saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_alpha%d_stype%d_cVal%d_%diters',nClus,nTrials,epsMuOrig1000,alpha10,stochasticType,c,nIter)];
    
    if warpBox
        fname = [fname '_warpBox'];
    end
    load(fname);
end

%%

savePlots=0;

colors = distinguishable_colors(nClus); %function for making distinguishable colors for plotting

figure; plot(tsseTrls);
title(sprintf('SSE for each iteration (nClus=%d)',nClus));
if savePlots
%     fname=[saveDir, sprintf('/covering_map_%d_clus_tsseTrls_%d_trls_eps%d_alpha%d_%diters',nClus,nTrials,epsMuOrig100,alpha10,nIter)];
    fname=[saveDir, sprintf('/covering_map_%d_clus_tsseTrls_%d_trls_eps%d_alpha%d_cVal%d_%diters',nClus,nTrials,epsMuOrig100,alpha10,c10k,nIter)];
    if warpBox
        fname = [fname '_warpBox'];
    end
    saveas(gcf,fname,'png');
end

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
    xlim(locRange); ylim(locRange);

    hold on;
    if i==2
        title(sprintf('Lowest 3 and Highest 3 SSE cluster locations (nClus=%d)',nClus));
    end
end
if savePlots
%     fname=[saveDir, sprintf('/covering_map_%d_clus_locs_top_bottom_3_%d_trls_eps%d_alpha%d_%diters',nClus,nTrials,epsMuOrig100,alpha10,nIter)];
    fname=[saveDir, sprintf('/covering_map_%d_clus_locs_top_bottom_3_%d_trls_eps%d_alpha%d_cVal%d_%diters',nClus,nTrials,epsMuOrig100,alpha10,c10k,nIter)];
    if warpBox
        fname = [fname '_warpBox'];
    end
    saveas(gcf,fname,'png');
end

%% plot cluster centers over time - subplots

savePlots=0;

iterI = 1; %1 to 3, lowest sse, 4:6, highest sse

colors = distinguishable_colors(nClus); %function for making distinguishable colors for plotting

% figure;   
figure('units','normalized','outerposition',[0 0 1 1]);
iPlot = 1; subplot(3,4,iPlot); hold on;%subplot
for iTrl = 1:nTrials
    if mod(iTrl,1000)==0
%         voronoi(muAll(:,1,iTrl,toPlot),muAll(:,2,iTrl,toPlot),'k'); %plot at the END before starting new subplot
        iPlot=iPlot+1;
        subplot(3,4,iPlot); hold on;
        voronoi(muAll(:,1,iTrl,iterI),muAll(:,2,iTrl,iterI),'k'); %plot at the START showing the previous final positions

    end
    xlim(locRange); ylim(locRange);

    if mod(iTrl,50)==0 %plot centers after x trials
        for i=1:nClus
            plot(squeeze(muAll(i,1,iTrl,iterI)),squeeze(muAll(i,2,iTrl,iterI)),'.','Color',colors(i,:),'MarkerSize',20); hold on;
        end
        drawnow;
    end
end
voronoi(muAll(:,1,iTrl,iterI),muAll(:,2,iTrl,iterI),'k'); %final one - if plotting at END above
% xlim([min(trials(:,1))-.1,max(trials(:,1))+.1]); ylim([min(trials(:,2))-.1,max(trials(:,2))+.1]);

if savePlots
    fname=[saveDir, sprintf('/covering_map_%d_clus_locs_plotOverTime_top_%d_trls_eps%d_alpha%d_%diters',nClus,nTrials,epsMuOrig100,alpha10,nIter)];
    if warpBox
        fname = [fname sprintf('_warpBox_resetEps%d',resetEps)];
    end
    saveas(gcf,fname,'png');
end

%% over time - one plot

iterI = 1;
colors = distinguishable_colors(nClus); %function for making distinguishable colors for plotting

figure;
% figure('units','normalized','outerposition',[0 0 1 1]);
for iTrl = 1:nTrials
    if mod(iTrl,250)==0
%         iPlot=iPlot+1;
%         voronoi(muAll(:,1,iTrl,iterI),muAll(:,2,iTrl,iterI),'k')
    end
    xlim(locRange); ylim(locRange);

    if mod(iTrl,50)==0, %plot centers after x trials
        for i=1:nClus
            plot(squeeze(muAll(i,1,iTrl,iterI)),squeeze(muAll(i,2,iTrl,iterI)),'.','Color',colors(i,:),'MarkerSize',10); hold on; %make marker size bigger - larger/smoother firing field!
        end
        drawnow;
    end
end

%% plot final centres 
colors = distinguishable_colors(nClus); %function for making distinguishable colors for plotting

iterI=1;

%plot over average of final trials
% nTrlsToPlot = [3000, 2000, 1000, 500, 200, 100];
nTrlsToPlot = 10000;

for iToPlot=1:length(nTrlsToPlot)
    figure;
    for iClus = 1:nClus
        plot(mean(squeeze(muAll(iClus,1,nTrials-nTrlsToPlot(iToPlot):nTrials,iterI))),squeeze(mean(muAll(iClus,2,nTrials-nTrlsToPlot(iToPlot):nTrials,iterI))),'.','Color',colors(iClus,:),'MarkerSize',20); hold on;
    end
    xlim(locRange); ylim(locRange);
%     voronoi(squeeze(mean(muAll(:,1,nTrials-nTrlsToPlot(iToPlot):end,toPlot),3)),squeeze(mean(muAll(:,2,nTrials-nTrlsToPlot(iToPlot):end,toPlot),3)),'k');
end


%  spreadout-ness measure
for iterI=1:nIter
    avgTrlSSE(:,iterI) = squeeze(mean(sseTrl(:,:,iterI),1));
    
    for iClus=1:nClus
        devAvgSSE(iClus,:,iterI)=sseTrl(iClus,:,iterI)-avgTrlSSE(:,iterI)';
    end
end
stdAcrossClus=squeeze(std(devAvgSSE));
varAcrossClus=squeeze(var(devAvgSSE));

% figure; plot(devAvgSSE);
figure; plot(varAcrossClus)
figure; plot(stdAcrossClus)

%% Plot / check stochastic update vals

%    figure; plot(cParams.closestChosen);
%    propClosestC = nnz(cParams.closestChosen)/nTrials
%    figure; plot(cParams.betaAll);


% stochasticType=1;
% 
% cVals = round([3.25/nTrials, 3.5/nTrials, 4/nTrials].*1000000);
% for iTestVal = 2:length(cVals)
% 
%     epsMuOrig1000 = 75;
%     nClus = 40;
%     nTrials = 30000;
%     nIter=10;
%     warpBox=0;
%     alpha10=5;
% 
%     %atm have only 1 c value run for alpha=2/8, but all stypes
%     %also have 2 and 3 for alpha=5, can check those out too
%     
%     c =cVals(iTestVal);
%     fname = [saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_alpha%d_stype%d_cVal%d_%diters',nClus,nTrials,epsMuOrig1000,alpha10,stochasticType,c,nIter)];
%     
%     load(fname);
%     
%     figure; plot(cParams.closestPr); %probably most useful    
%     
%     propClosestC(iTestVal) = nnz(cParams.closestPr==1)/nTrials; % proportion of times the closesst cluster was chosen 
% end
%%
% %plot clusters at certain times of learning?
% snapShots = [1500, 10000, 15000, 20000, 25000]; %first one is at learning, others are already more stable but want to check how stable and compare across parameters used
% trls2Plot = [];
% for iShots = 1:length(snapShots)
%     trls2Plot = [trls2Plot, snapShots(iShots)-500, snapShots(iShots), snapShots(iShots)+500];
% end
% 
% figure; hold on;
% for iShots=1:length(trls2Plot)
%     iTrl=trls2Plot(iShots);
%     subplot(length(snapShots),3,iShots);
%     for i=1:nClus
%         plot(squeeze(muAll(i,1,iTrl,toPlot)),squeeze(muAll(i,2,iTrl,toPlot)),'.','Color',colors(i,:),'MarkerSize',10); hold on; %make marker size bigger - larger/smoother firing field!
%     end
% end
