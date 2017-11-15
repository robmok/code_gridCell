%% Plot (with/without load in data)
clear all;

wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
% wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
cd(wd);
saveDir = [wd '/data_gridCell'];

% save plots?
savePlots=0;

%set to load in data (note divide by value in variable)
if 0
%     epsMuOrig100 = 75;
    epsMuOrig1000 = 100;
%     deltaEpsMu100 = 96; %96, 98
    nClus = 40;
    nTrials = 30000;
%     nIter=5;
    nIter=20;
    warpBox=0;
    alpha10=9; %0,2, 5, 9, 95 (95 looks like it goes really bad)
    
%     c10k=9990000; %25, 50, 100, 9990000
    
    
%     fname=[saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_deltaMu%d_%diters',nClus,nTrials,epsMuOrig10,deltaEpsMu100,nIter)];
    fname=[saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_alpha%d_%diters',nClus,nTrials,epsMuOrig1000 ,alpha10,nIter)];
%     fname=[saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_alpha%d_cVal%d_%diters',nClus,nTrials,epsMuOrig1000,alpha10,c10k,nIter)];
    
    if warpBox,
        fname = [fname '_warpBox'];
    end
    load(fname);
end

%%

savePlots=0;

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

iter = 1; %1 to 3, lowest sse, 4:6, highest sse

colors = distinguishable_colors(nClus); %function for making distinguishable colors for plotting

% figure;   
figure('units','normalized','outerposition',[0 0 1 1]);
iPlot = 1; subplot(3,4,iPlot); hold on;%subplot
for iTrl = 1:nTrials
    if mod(iTrl,1000)==0
%         voronoi(muAll(:,1,iTrl,toPlot),muAll(:,2,iTrl,toPlot),'k'); %plot at the END before starting new subplot
        iPlot=iPlot+1;
        subplot(3,4,iPlot); hold on;
        voronoi(muAll(:,1,iTrl,iter),muAll(:,2,iTrl,iter),'k'); %plot at the START showing the previous final positions

    end
    xlim(locRange); ylim(locRange);

    if mod(iTrl,50)==0, %plot centers after x trials
        for i=1:nClus
            plot(squeeze(muAll(i,1,iTrl,iter)),squeeze(muAll(i,2,iTrl,iter)),'.','Color',colors(i,:),'MarkerSize',20); hold on;
        end
        drawnow;
    end
end
voronoi(muAll(:,1,iTrl,iter),muAll(:,2,iTrl,iter),'k'); %final one - if plotting at END above
% xlim([min(trials(:,1))-.1,max(trials(:,1))+.1]); ylim([min(trials(:,2))-.1,max(trials(:,2))+.1]);

if savePlots,
    fname=[saveDir, sprintf('/covering_map_%d_clus_locs_plotOverTime_top_%d_trls_eps%d_alpha%d_%diters',nClus,nTrials,epsMuOrig100,alpha10,nIter)];
    if warpBox
        fname = [fname sprintf('_warpBox_resetEps%d',resetEps)];
    end
    saveas(gcf,fname,'png');
end

%% over time - one plot

iter = 1;
colors = distinguishable_colors(nClus); %function for making distinguishable colors for plotting

figure;
% figure('units','normalized','outerposition',[0 0 1 1]);
for iTrl = 1:nTrials
    if mod(iTrl,250)==0
%         iPlot=iPlot+1;
%         voronoi(muAll(:,1,iTrl,toPlot),muAll(:,2,iTrl,toPlot),'k')
    end
    xlim(locRange); ylim(locRange);

    if mod(iTrl,50)==0, %plot centers after x trials
        for i=1:nClus
            plot(squeeze(muAll(i,1,iTrl,iter)),squeeze(muAll(i,2,iTrl,iter)),'.','Color',colors(i,:),'MarkerSize',10); hold on; %make marker size bigger - larger/smoother firing field!
        end
        drawnow;
    end
end

%% plot final centres 
iter=1;
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

nTrlsToPlot = 100;

for iToPlot=1:length(nTrlsToPlot)
    figure;
    for iClus = 1:nClus
        plot(mean(squeeze(muAll(iClus,1,nTrials-nTrlsToPlot(iToPlot):end,iter))),squeeze(mean(muAll(iClus,2,nTrials-nTrlsToPlot(iToPlot):end,iter))),'.','Color',colors(iClus,:),'MarkerSize',20); hold on;
    end
    xlim(locRange); ylim(locRange);

%     voronoi(squeeze(mean(muAll(:,1,nTrials-nTrlsToPlot(iToPlot):end,toPlot),3)),squeeze(mean(muAll(:,2,nTrials-nTrlsToPlot(iToPlot):end,toPlot),3)),'k');
end

%

for iter=1:nIter
    avgTrlSSE(:,iter) = squeeze(mean(sseTrl(:,:,iter),1));
    
    for iClus=1:nClus,
        devAvgSSE(iClus,:,iter)=sseTrl(iClus,:,iter)-avgTrlSSE(:,iter)';
    end
end

stdAcrossClus=squeeze(std(devAvgSSE));
varAcrossClus=squeeze(var(devAvgSSE));

% figure; plot(devAvgSSE);
figure; plot(varAcrossClus)
figure; plot(stdAcrossClus)

%% comparing parameters

% at the moment I didn't save indSSE1 and 2 - so not immediately sure which
% tsseTrls and sseTrls are the top3bottom 3. now i'm saving them.

% more importantly, i didn't save all the cluster positions - i will do
% this now (either downsampled or averaging over the last x trials). Also
% order the clusters by SSE - after averaging cluster positions over 10k
% trials


%for now i will compute SSE on new datasets on cluster positions for top3
%and bottom3 only. this probably is good enough for now




%load in different versions

%%%%
% loop and plot different versions? or loop and load in variables with diff
% names?
%%%%

% if 1
%     epsMuOrig1000 = 75; %75, 100
%     nClus = 40;
%     colors = distinguishable_colors(nClus); %function for making distinguishable colors for plotting
%     nTrials = 30000;
%     nIter=20;
%     warpBox=0;
%     alpha10=5; %0,2, 5, 9, 95 (95 looks like it goes really bad)
% %     c10k=9990000; %25, 50, 100, 9990000
%     
%     fname=[saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_alpha%d_%diters',nClus,nTrials,epsMuOrig1000 ,alpha10,nIter)];
% %     fname=[saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_alpha%d_cVal%d_%diters',nClus,nTrials,epsMuOrig1000,alpha10,c10k,nIter)];
% %     if warpBox,
% %         fname = [fname '_warpBox'];
% %     end
%     load(fname);
% end





%make test data - keep same for comparison
locRange = 1;
% nTrialsTest = nTrials;
nTrialsTest = 50000;
dataPtsTest = [randsample(linspace(-locRange,locRange,101),nTrialsTest,'true'); randsample(linspace(-locRange,locRange,101),nTrialsTest,'true')]'; % random points in a box

%%
    
%define which parmeters to test (if not testing, keep at 1 value)
%%%%%%%
%how to test several; or should test one at a time?
%%%%%%%
% alphaVals = [0,2,5,6,7,8,9];
alphaVals = [0,2,5,9];
% alphaVals = 9;
% cVals = [nTrials/(nTrials*100), nTrials/(nTrials*100*2), nTrials/(nTrials*100*4), 999]; % larger c = less stochastic over trials (becomes det quite early on); smaller c = more stochastic over trials (still a bit stochastic by the end)
% cVals = cVals*10000;
% testVals = cVals;


%may want to test epsMu too later?


testVals = alphaVals;

nClus = 40; colors = distinguishable_colors(nClus);
nTrials = 30000;
nIter = 20;
%warpBox=0;

nIterSaved = 6;  %edit all nIterSaved to nIter later!


clusMu=nan(nClus,2,nIterSaved,length(testVals));
sse=nan(1,nClus);
tsse=nan(nIterSaved,length(testVals));
devAvgSSE=nan(nClus,nIterSaved,length(testVals));
for iDataSet = 1:length(testVals)
    
    %set params, load in data
    epsMuOrig1000 = 75; %75, 100


    alpha10 = alphaVals(iDataSet);
%     alpha10 = alphaVals;
    
    fname=[saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_alpha%d_%diters',nClus,nTrials,epsMuOrig1000,alpha10,nIter)];
%     c10k = cVals(iDataSet);
%     fname=[saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_alpha%d_cVal%d_%diters',nClus,nTrials,epsMuOrig1000,alpha10,c10k,nIter)];
    load(fname);
    
    iter = 1:nIterSaved; %which of top3/bottom3 to use. all?
    for iterI = 1:nIterSaved,
        
        nTrlsToAvg = 10000; % also check 5k, 15k, 20k - could make a loop to check all later
        
        clusMu(:,:,iterI,iDataSet)=[mean(muAll(:,1,nTrials-nTrlsToAvg:end,iter(iterI)),3),mean(muAll(:,2,nTrials-nTrlsToAvg:end,iter(iterI)),3)];
        
        % compute sse with respect to all test data points
        distTrl=[];
        for iClus = 1:size(clusMu,1),
            distTrl(:,iClus)=sum([clusMu(iClus,1,iterI,iDataSet)-dataPtsTest(:,1), clusMu(iClus,2,iterI,iDataSet)-dataPtsTest(:,2)].^2,2);
        end
        [indValsTrl, indTmp]=min(distTrl,[],2); % find which clusters are points closest to
        for iClus = 1:size(clusMu,1),
            sse(iClus)=sum(sum([clusMu(iClus,1,iterI,iDataSet)-dataPtsTest(indTmp==iClus,1), clusMu(iClus,2,iterI,iDataSet)-dataPtsTest(indTmp==iClus,2)].^2,2)); %distance from each cluster from training set to datapoints closest to that cluster
        end
        tsse(iterI,iDataSet)=sum(sse);
        
        %compute 'spreaded-ness' - variance of SE across clusters is a measure
        %of this, assuing uniform data points
        devAvgSSE(:,iterI,iDataSet)   = sse-mean(sse);
        stdAcrossClus(iterI,iDataSet) = std(devAvgSSE(:,iterI,iDataSet));
        varAcrossClus(iterI,iDataSet) = var(devAvgSSE(:,iterI,iDataSet));
    end
end






% figure; plot(tsse);
% figure; plot(mean(tsse,1));
% figure; plot(mean(tsse(1:3,:),1));
% figure; plot(mean(tsse(4:6,:),1));
% 
% 
% figure; plot(stdAcrossClus);
% figure; plot(mean(stdAcrossClus,1));
% figure; plot(mean(stdAcrossClus(1:3,:),1));
% figure; plot(mean(stdAcrossClus(4:6,:),1));


% plot cluster centres, print their sse and stdAcrossClus


for iterI=1:nIter
    % figure('units','normalized','outerposition',[0 0 1 1]); hold on; iPlot=0;
    figure; hold on; iPlot=0;
    for iDataSet = 1:length(testVals)
        
        iPlot = iPlot+1;
        
        subplot(2,2,iPlot);
        
        scatter(clusMu(:,1,iterI,iDataSet),clusMu(:,2,iterI,iDataSet),20e+2,colors,'.')
        
        title(sprintf('%d alpha, %.1f tsse, %.2f. spreadedness', testVals(iDataSet), tsse(iterI,iDataSet), stdAcrossClus(iterI,iDataSet)))
%         title(sprintf('%d cVal, %.1f tsse, %.2f. spreadedness', testVals(iDataSet), tsse(iterI,iDataSet), stdAcrossClus(iterI,iDataSet)))
    end
end

%% make logical map of cluster positions, density plots

nTrlsToAvg = 5000;

spacing = linspace(locRange(1),locRange(2),locRange(2)+1);

for iter=1, 
%     iter=1;%just get 1 iter, 1 of the param settings for now

clus = round(muAll(:,:,nTrials-nTrlsToAvg+1:end,iter)); 
% nDecimals = 3; f = 10.^nDecimals; %round to 3 dec places
% clus = round(f.*clusTmp)./f;
% clus = (clus*1000)+1000; %add 1000 to make it into the range of 0-2001


densityPlot(:,:,iter) = zeros(length(spacing),length(spacing));

% clus(:,:,1)
% clus(1,:,1)

for iTrls=1:nTrlsToAvg,
    for i=1:nClus
        densityPlot(clus(i,1,iTrls),clus(i,2,iTrls),iter)=1; % works, but better way / faster to vectorise?
    end
end

densityPlot(:,:,iter) = imgaussfilt(densityPlot(:,:,iter),10);

figure;
imagesc(densityPlot(:,:,iter));
% imagesc(densityPlot(:,:,iter)./mean(reshape(densityPlot(:,:,iter),1,numel(densityPlot(:,:,iter)))));

end

% figure;
% imagesc(densityPlot(:,:,iter));



% also should think about the best way to 'add' points to the map....
% logging each point as a 1 (or +1) seems very slow...





%%
% %plot clusters at certain times of learning
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







% plot these values over their timepoints above
% iter = 1;
% tsseTrls(iTrl,iter)
% stdAcrossClus(iTrl)



%atm - tsseTrls and sseTrl are not orders by lowest SSE... I suppose this
%is not hard to check. BUT if want to plot top3/bottom 3, should have this
%info at hand...

% notable that i only save best/worst 3 from the LAST trial; might be worth
% to save the AVERAGE SSE over last 1000 trials to make it more robust? or
% average cluster positions then compute the sSE
    % - need to do this in the sim script
    % - also might be worth checking ranking by the spreaded-outness
    % (stdAcrosClus)

% code it in sim.m too, but for the ones already run, recompute
tsseIter = mean(tsseTrls(end-10000:end,:)); %average over last x trials 
[indVal, indSSE1] = sort(tsseIter);
[y, indSSE2] = sort(indSSE1); %use indSSE2 to sort for plotting

figure;
plot(tsseTrls,'color',[.5 .5 .5]); hold on;
plot(tsseTrls(:,indSSE2==1))

figure;
plot(stdAcrossClus,'color',[.5 .5 .5]); hold on;
plot(stdAcrossClus(:,indSSE2==1)); %interesting - lowest spread here is diff to lowest SSE. prob corr, but
%still


%sort according to spread now
tsseVarIter = mean(stdAcrossClus(end-10000:end,:)); %average over last x trials 
[indVal, indSSEvar1] = sort(tsseVarIter);
[y, indSSEvar2] = sort(indSSEvar1); %use indSSE2 to sort for plotting

figure;
plot(tsseTrls,'color',[.5 .5 .5]); hold on;
plot(tsseTrls(:,indSSEvar2==1))

figure;
plot(stdAcrossClus,'color',[.5 .5 .5]); hold on;
plot(stdAcrossClus(:,indSSEvar2==1));


% some that do well at the last 10k were quite bad before that... need to
% think about this a bit. Also the initial learning rate: need to look at
% alpha = 1, but maybe also need to sim less than .75 and more than 1?












%plot over average of final trials - taken from above
% % need to add a loop for iterations? or plot more than just one?

% nTrlsToPlot = [1, 100, 500, 1000];
% figure; hold on;
% for iToPlot=1:length(nTrlsToPlot)
%     subplot(2,2,iToPlot);
%     for iClus = 1:nClus
%         plot(mean(squeeze(muAll(iClus,1,nTrials-nTrlsToPlot(iToPlot):end,toPlot))),squeeze(mean(muAll(iClus,2,nTrials-nTrlsToPlot(iToPlot):end,toPlot))),'.','Color',colors(i,:),'MarkerSize',20); hold on;
%     end
%     title(sprintf('Avg over %d trials',nTrlsToPlot(iToPlot)));
%     if warpBox,
%         xlim([-1.1,2.1]); ylim([-1.1,1.1]);
%     else
%         xlim([-1.1,1.1]); ylim([-1.1,1.1]);
%     end
% %     voronoi(squeeze(mean(muAll(:,1,nTrials-nTrlsToPlot(iToPlot):end,toPlot),3)),squeeze(mean(muAll(:,2,nTrials-nTrlsToPlot(iToPlot):end,toPlot),3)),'k');
% end
% 













