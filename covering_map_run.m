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
nTrials = 20000; %how many locations in the box / trials - 2.5k ; 5k if reset

colgrey = [.5, .5, .5];

%box
nSteps = 50; %to define spacing beween each loc in box
locRange = 1;%1; from -locRange to locRange
stepSize=diff(linspace(-1,1,nSteps)); stepSize=stepSize(1); %smallest diff between locs

% parameters
epsMuOrig=.075;% %learning rate / starting learning rate %.075
% epsMuOrig=.1;%
% deltaEpsMu = .96;% %change in learning rate over time (slow down with 'learning')
% deltaEpsMu = .99; % slower decrease in learning rate for expanding (if no
% reset)

% for saving simulations - multiple by values to save the files with params
epsMuOrig1000=epsMuOrig*1000;
% deltaEpsMu100 = deltaEpsMu*100;

% testing with 60 clusters, 5k trials
% starting learning rate = .6x



    
%diff number of clusters
% 30 clusters


%different starting learning rate
%30 clusters



%reset learning rate?
% resetEps=0; %0 - no, 1 - once halfway, 2 - twice (quarter way and half way)  % ATM resetting to halfway
% resetMag=2; %1 - back to orig, 2 - 75% of the orig (if resetEps=2, 75%
% the orig then 50% the orig) %- not using this yet

%define box / environement - random points in a box
box = 'square'; %square, rect, trapz, trapzSq (trapz and a square box attached)

% change box shape during learning
%rectangle
warpBox = 0; %1 or 0
warpType = 'sq2rect';

%mometum-like adaptive learning rate - define alpha (higher = weight
%previous update (direction and magnitude) more; 0 = don't weight previous at all)
alpha = 0.8;
alpha10 = alpha*10; %for saving simulations - multiple by values to save the files with params



%%% things to understand
% - plot changes in the updates with a higher alpha (should slow down, but
% not sure

% - how does the first trial affect it? (now setting all to 0.5 or sth) -
% if start at 0 should be slow first, then fast then slow, but when i set
% it high, still like that? setting it high, at 1.5 still doesnt change it
% much.... - 0.9 ad 0.2 both get good solutions, maybe doesnt matter as
% long as it takes prev trial into acc?



%%
saveDat=0; %save simulations - cluster centres and tsse

nIter=1; %how many iterations (starting points)
tic
for iClus2run = 1:length(clus2run)    
    nClus = clus2run(iClus2run);
    [muEnd, muAll, tsseTrls,sseTrl,epsMuAll,deltaMu,clusUpdAll] = covering_map_sim(nClus,locRange,box,warpType,epsMuOrig,nTrials,nIter,warpBox,alpha);

    if saveDat,
        % save simulations - might not need all tsseTrls if not viewing
        % them or epsMuAll?
        fname = [saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_alpha%d_%diters',nClus,nTrials,epsMuOrig1000,alpha10,nIter)];
        if warpBox,
            fname = [fname '_warpBox'];
        end
%         if resetEps,
%             fname = [fname '_resetHalfwayHalfEps'];
%         end
        save(fname,'muEnd','muAll','sseTrl','tsseTrls','nIter');
    end
end
toc



%% Plot (with/without load in data)

% save plots?
savePlots=0;

%set to load in data (note divide by value in variable)
if 0
    epsMuOrig100 = 75;
%     deltaEpsMu100 = 96; %96, 98
    nClus = 40;
    nTrials = 30000;
    nIter=1;
%     nIter=20; %for warpBox, for now
    warpBox=0;
    alpha10=8;
%     fname=[saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_deltaMu%d_%diters',nClus,nTrials,epsMuOrig10,deltaEpsMu100,nIter)];
    fname=[saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_alpha%d_%diters',nClus,nTrials,epsMuOrig100,alpha10,nIter)];
    
    if warpBox,
        fname = [fname '_warpBox'];
    end
    load(fname);
end

colors = distinguishable_colors(nClus); %function for making distinguishable colors for plotting

figure; plot(tsseTrls);
title(sprintf('SSE for each iteration (nClus=%d)',nClus));
if savePlots,
    fname=[saveDir, sprintf('/covering_map_%d_clus_tsseTrls_%d_trls_eps%d_alpha%d_%diters',nClus,nTrials,epsMuOrig100,alpha10,nIter)];
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
    fname=[saveDir, sprintf('/covering_map_%d_clus_locs_top_bottom_3_%d_trls_eps%d_alpha%d_%diters',nClus,nTrials,epsMuOrig100,alpha10,nIter)];
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
%         drawnow;
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
            plot(squeeze(muAll(i,1,iTrl,toPlot)),squeeze(muAll(i,2,iTrl,toPlot)),'.','Color',colors(i,:),'MarkerSize',40); hold on; %make marker size bigger - larger/smoother firing field!
        end
%         drawnow;
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






