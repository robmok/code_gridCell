clear all;
% close all;

% wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';

cd(wd);
codeDir = [wd '/code_gridCell'];
saveDir = [wd '/data_gridCell'];
addpath(codeDir); addpath(saveDir);

% nClus   = 20;
clus2run = [20, 40, 60, 80]; %run multuple cluster numbers
clus2run = 40;
nTrials = 2500; %how many locations in the box / trials - 2.5k ; 5k if reset

colgrey = [.5, .5, .5];

%box
nSteps = 50; %to define spacing beween each loc in box
locRange = 1;%1; from -locRange to locRange
stepSize=diff(linspace(-1,1,nSteps)); stepSize=stepSize(1); %smallest diff between locs

% parameters
epsMuOrig=.6;% %learning rate / starting learning rate %.6
deltaEpsMu = .96;% %change in learning rate over time (slow down with 'learning')
% deltaEpsMu = .99; % slower decrease in learning rate for expanding (if no
% reset)

% for saving simulations - multiple by values to save the files with params
epsMuOrig10=epsMuOrig*10;
deltaEpsMu100 = deltaEpsMu*100;

% testing with 60 clusters, 5k trials
% starting learning rate = .6
% deltaEpsMu:
    %.99   - SSE ~ 60 - fluctuating more at end, not stabilising
    %.98   - SSE ~.58
    %.96/7 - SSE ~ 56-7
    %.95   - SSE ~56-7 (decrease in epsmu curve looks gd, exponential, half near zero; 
    %.94   - SSE ~ 56 (curve gd, steeper, most near zero at end)
    %.935  - SSE ~ 59-60 (curve gd, steeper, maybe too steep? most near zero at end)

    
%diff number of clusters
% 30 clusters
% deltaEpsMu:
    %.99  - SSE ~120-124 - curve not near 0
    %.985 - SSE 120
    %.98  - SSE 116/117 - curve now steep but not all 0 at end
    %.97  - SSE ~ 117 - curve pretty steep
    %.95  - SSE ~ 114 - curve too steep, at 0 
    
    % although curves are different when less clusters, clusters look hexagonal anyway - less
    % clusters easier?
%40 cluster, .96; SSE~90; .97; SSE~88-91; .98; SSE~90
    

%different starting learning rate
%30 clusters
%.99, SSE = 124 (too shallow,), .98; SSE~118-123, .96; SSE~115 (steeper), .95, SSE~ 120 (v steep)

%reset learning rate?
resetEps=0; %0 - no, 1 - once halfway, 2 - twice (quarter way and half way)  % ATM resetting to halfway
% resetMag=2; %1 - back to orig, 2 - 75% of the orig (if resetEps=2, 75%
% the orig then 50% the orig) %- not using this yet

%define box / environement - random points in a box
box = 'square';

% change box shape during learning
%rectangle
warpBox = 0; %1 or 0
warpType = 'sq2rect';

%%
saveDat=0; %save simulations - cluster centres and tsse

nIter=6; %how many iterations (starting points)
tic
for iClus2run = 1:length(clus2run)    
    nClus = clus2run(iClus2run);
    [muEnd, muAll, tsseTrls, epsMuAll] = covering_map_sim(nClus,locRange,box,warpType,epsMuOrig,deltaEpsMu,nTrials,nIter,warpBox,resetEps);

    if saveDat,
        % save simulations - might not need all tsseTrls if not viewing
        % them or epsMuAll?
        fname = [saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_deltaMu%d_%diters',nClus,nTrials,epsMuOrig10,deltaEpsMu100,nIter)];
        if warpBox,
            fname = [fname '_warpBox'];
        end
        if resetEps,
            fname = [fname '_resetHalfwayHalfEps'];
        end
        save(fname,'muEnd','muAll', 'tsseTrls', 'epsMuAll','nIter');
    end
end
toc
%% Plot (with/without load in data)

% save plots?
savePlots=0;

%set to load in data (note divide by value in variable)
if 1
    epsMuOrig10 = 6;
    deltaEpsMu100 = 96; %96, 98
    nClus = 40;
    nTrials = 2500;
    nIter=500;
%     nIter=20; %for warpBox, for now
    warpBox=0;
    fname=[saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_deltaMu%d_%diters',nClus,nTrials,epsMuOrig10,deltaEpsMu100,nIter)];
    if warpBox,
        fname = [fname '_warpBox'];
    end
    load(fname);
end

colors = distinguishable_colors(nClus); %function for making distinguishable colors for plotting

figure; plot(tsseTrls);
title(sprintf('SSE for each iteration (nClus=%d)',nClus));
if savePlots,
    fname=[saveDir, sprintf('/covering_map_%d_clus_tsseTrls_%d_trls_eps%d_deltaMu%d_%diters',nClus,nTrials,epsMuOrig10,deltaEpsMu100,nIter)];
    if warpBox,
        fname = [fname '_warpBox'];
    end
    saveas(gcf,fname,'png');
end

% plot change in learning rate for each cluster:
figure;
plot(epsMuAll(:,1),'.');

% plot lowest and highest SSE
iToPlot=[1,2,3,nIter-2,nIter-1,nIter];
figure;
for i = 1:6
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
    fname=[saveDir, sprintf('/covering_map_%d_clus_locs_top_bottom_3_%d_trls_eps%d_deltaMu%d_%diters',nClus,nTrials,epsMuOrig10,deltaEpsMu100,nIter)];
    if warpBox,
        fname = [fname '_warpBox'];
    end
    saveas(gcf,fname,'png');
end

%% plot cluster centers over time - subplots

savePlots=0;

toPlot = 1; %1 to 3, lowest sse, 4:6, highest sse

% figure;   
figure('units','normalized','outerposition',[0 0 1 1]);
iPlot = 1; subplot(3,4,iPlot); hold on;%subplot
for iTrl = 1:nTrials
    if mod(iTrl,250)==0
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
    fname=[saveDir, sprintf('/covering_map_%d_clus_locs_plotOverTime_top_%d_trls_eps%d_deltaMu%d_%diters',nClus,nTrials,epsMuOrig10,deltaEpsMu100,nIter)];
    if warpBox
        fname = [fname sprintf('_warpBox_resetEps%d',resetEps)];
    end
    saveas(gcf,fname,'png');
end

%% over time - one plot

figure('units','normalized','outerposition',[0 0 1 1]);
for iTrl = 1:nTrials
    if mod(iTrl,250)==0
%         iPlot=iPlot+1;
        voronoi(muAll(:,1,iTrl,toPlot),muAll(:,2,iTrl,toPlot),'k')
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


%% plot final centres 
for i = 1:nIter
    figure;
    for iClus = 1:nClus
        plot(squeeze(muEnd(iClus,1,i)),squeeze(muEnd(iClus,2,i)),'.','Color',colors(iClus,:),'MarkerSize',25); hold on; %plot cluster final point
    end
    xlim([min(trials(:,1))-.1,max(trials(:,1))+.1]); ylim([min(trials(:,2))-.1,max(trials(:,2))+.1]);
end
voronoi(squeeze(muEnd(:,1,i)),squeeze(muEnd(:,2,i)))

