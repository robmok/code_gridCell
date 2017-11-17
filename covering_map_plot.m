%% Plot (with/without load in data)
clear all;

% wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
cd(wd);
saveDir = [wd '/data_gridCell'];

nSteps = 500;
locRange = [0, nSteps-1]; 

% save plots?
% savePlots=0;

%set to load in data (note divide by value in variable)
if 0
    epsMuOrig1000 = 75;
    nClus = 40;
    nTrials = 30000;
    nIter=10;
    warpBox=0;
    alpha10=8; 

    
    %atm have only 1 c value run for alpha=2/8, but all stypes
    %also have 2 and 3 for alpha=5, can check those out too
    stochasticType=3;
    cVals = [2/nTrials, 3/nTrials, 5/nTrials, 10/nTrials];
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
        voronoi(muAll(:,1,iTrl,iterI),muAll(:,2,iTrl,iterI),'k')
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

%% comparing parameters

%make test data - keep same for comparison
nSteps = 500;
locRange = [0, nSteps-1]; 
nTrialsTest = nTrials;
% nTrialsTest = 50000;
dataPtsTest = [randsample(linspace(locRange(1),locRange(2),nSteps),nTrialsTest,'true'); randsample(linspace(locRange(1),locRange(2),nSteps),nTrialsTest,'true')]'; % random points in a box

%%
    
%define which parmeters to test (if not testing, keep at 1 value)
%%%%%%%
%how to test several; or should test one at a time?
%%%%%%%
% alphaVals = [0,2,5,6,7,8,9];
alphaVals = [0,2,4,6, 8];
% alphaVals = 9;
% cVals = [nTrials/(nTrials*100), nTrials/(nTrials*100*2), nTrials/(nTrials*100*4), 999]; % larger c = less stochastic over trials (becomes det quite early on); smaller c = more stochastic over trials (still a bit stochastic by the end)
% cVals = cVals*10000;
% testVals = cVals;

testVals = alphaVals;

nClus = 40; colors = distinguishable_colors(nClus);
nTrials = 30000;
nIter = size(muAll,4);
%warpBox=0;

clusMu=nan(nClus,2,nIter,length(testVals));
sse=nan(1,nClus);
tsse=nan(nIter,length(testVals));
devAvgSSE=nan(nClus,nIter,length(testVals));
for iTestVal = 1:length(testVals)
    
    %set params, load in data
    epsMuOrig1000 = 75; %75, 100

    alpha10 = alphaVals(iTestVal);
%     alpha10 = alphaVals;
    
    fname=[saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_alpha%d_%diters',nClus,nTrials,epsMuOrig1000,alpha10,nIter)];
%     c10k = cVals(iDataSet);
%     fname=[saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_alpha%d_cVal%d_%diters',nClus,nTrials,epsMuOrig1000,alpha10,c10k,nIter)];
    load(fname);
    
    iterI = 1:nIter; %which of top3/bottom3 to use. all?
    for iterI = 1:nIter
        
        nTrlsToAvg = 10000; % also check 5k, 15k, 20k - could make a loop to check all later
        
        clusMu(:,:,iterI,iTestVal)=[mean(muAll(:,1,nTrials-nTrlsToAvg:end,iterI(iterI)),3),mean(muAll(:,2,nTrials-nTrlsToAvg:end,iterI(iterI)),3)];
        
        % compute sse with respect to all test data points
        distTrl=[];
        for iClus = 1:size(clusMu,1)
            distTrl(:,iClus)=sum([clusMu(iClus,1,iterI,iTestVal)-dataPtsTest(:,1), clusMu(iClus,2,iterI,iTestVal)-dataPtsTest(:,2)].^2,2);
        end
        [indValsTrl, indTmp]=min(distTrl,[],2); % find which clusters are points closest to
        for iClus = 1:size(clusMu,1)
            sse(iClus)=sum(sum([clusMu(iClus,1,iterI,iTestVal)-dataPtsTest(indTmp==iClus,1), clusMu(iClus,2,iterI,iTestVal)-dataPtsTest(indTmp==iClus,2)].^2,2)); %distance from each cluster from training set to datapoints closest to that cluster
        end
        tsse(iterI,iTestVal)=sum(sse);
        
        %compute 'spreaded-ness' - variance of SE across clusters is a measure
        %of this, assuing uniform data points
        devAvgSSE(:,iterI,iTestVal)   = sse-mean(sse);
        stdAcrossClus(iterI,iTestVal) = std(devAvgSSE(:,iterI,iTestVal));
        varAcrossClus(iterI,iTestVal) = var(devAvgSSE(:,iterI,iTestVal));
    end
end

% figure; plot(tsse');
% figure; plot(mean(tsse,1));
% figure; plot(mean(tsse(1:3,:),1));
% figure; plot(mean(tsse(4:6,:),1));
% 
% 
% figure; plot(stdAcrossClus');
% figure; plot(mean(stdAcrossClus,1));
% figure; plot(mean(stdAcrossClus(1:3,:),1));
% figure; plot(mean(stdAcrossClus(4:6,:),1));


% plot cluster centres, print their sse and stdAcrossClus
for iterI=1:nIter
    % figure('units','normalized','outerposition',[0 0 1 1]); hold on; iPlot=0;
    figure; hold on; iPlot=0;
    for iTestVal = 1:length(testVals)
        iPlot = iPlot+1;
        subplot(2,2,iPlot);
        scatter(clusMu(:,1,iterI,iTestVal),clusMu(:,2,iterI,iTestVal),20e+2,colors,'.')
        xlim(locRange); ylim(locRange);
        title(sprintf('%d alpha, %.1f tsse, %.2f. spreadedness', testVals(iTestVal), tsse(iterI,iTestVal), stdAcrossClus(iterI,iTestVal)))
%         title(sprintf('%d cVal, %.1f tsse, %.2f. spreadedness', testVals(iDataSet), tsse(iterI,iDataSet), stdAcrossClus(iterI,iDataSet)))
    end
end
%% make logical map of cluster positions, density plots

nTrlsToUse = 10000;
spacing = linspace(locRange(1),locRange(2),locRange(2)+1);
densityPlot = zeros(length(spacing),length(spacing),nIter);
for iterI=1:10
    clus = round(muAll(:,:,nTrials-nTrlsToUse+1:nTrials,iterI));
    for iTrl=1:nTrlsToUse
        for i=1:nClus
            densityPlot(clus(i,1,iTrl),clus(i,2,iTrl),iterI)=densityPlot(clus(i,1,iTrl),clus(i,2,iTrl),iterI)+1; % works, but better way / faster to vectorise?
        end
    end
    densityPlot(:,:,iterI) = imgaussfilt(densityPlot(:,:,iterI),12); %smooth
    figure;
    imagesc(densityPlot(:,:,iterI));
    % imagesc(densityPlot(:,:,iter),[100 800]);
end


%% Plot / check stochastic update params...

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
%     
%     
% end

%% Compute the cluster centres from the desnity map, then compute SSE and rank them

gaussSmooth=9; %smooth maps by x value

%define which parmeters to test (if not testing, keep at 1 value)
%%%%%%%
%how to test several; or should test one at a time?
%%%%%%%
alphaVals = [0,2,4,6,8]; 

stochasticType = 1; %1, 2, 3; % stochasticVals=[1,2,3] - later

cVals = round([3.25/nTrials, 3.5/nTrials, 4/nTrials].*1000000);
%%%%%%%
%set which one to test
%%%%%%
testVals = alphaVals; %alphaVals, cVals, stochasticVals
nTestVals = length(testVals);

nClus = 40; colors = distinguishable_colors(nClus);
nTrials = 30000;
nIter = size(muAll,4);
%warpBox=0;

clusMu=nan(nClus,2,nIter,length(testVals));
sse=nan(1,nClus);
tsse=nan(nIter,length(testVals));
devAvgSSE=nan(nClus,nIter,length(testVals));

for iTestVal = 1:nTestVals
    %set params, load in data
    epsMuOrig1000 = 75; %75, 100
    
    %alpha
    alpha10 = alphaVals(iTestVal);
%     alpha10 = alphaVals;
%     fname=[saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_alpha%d_%diters',nClus,nTrials,epsMuOrig1000,alpha10,nIter)];
    
    %stochastic parameters
    c =cVals(iTestVal);    
    fname = [saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_alpha%d_stype%d_cVal%d_%diters',nClus,nTrials,epsMuOrig1000,alpha10,stochasticType,c,nIter)];
    
    load(fname);
        
%         nTrlsToAvg = 10000; % also check 5k, 15k, 20k - could make a loop to check all later
    
    nTrlsToUse = 10000; % - for computing density map
    
    spacing = linspace(locRange(1),locRange(2),locRange(2)+1);
    densityPlot = zeros(length(spacing),length(spacing),nIter,nTestVals);
    for iterI = 1:nIter
        
        %new - here compute the smoothed density map - find the peaks,
        %treat those as the cluster centres, compute SSE
        
        %compute density map
        clus = round(muAll(:,:,nTrials-nTrlsToUse+1:nTrials,iterI));
        for iTrl=1:nTrlsToUse
            for i=1:nClus
                densityPlot(clus(i,1,iTrl),clus(i,2,iTrl),iterI,iTestVal)=densityPlot(clus(i,1,iTrl),clus(i,2,iTrl),iterI,iTestVal)+1; % works, but better way / faster to vectorise?
            end
        end
        densityPlot(:,:,iterI,iTestVal) = imgaussfilt(densityPlot(:,:,iterI,iTestVal),gaussSmooth); %smooth
        
        %find peaks - how to just get one centre value per peak? each
        %centre might have diff number of counts...
        
        
        
        
        
%         clusMu(:,:,iterI,iDataSet)=[mean(muAll(:,1,nTrials-nTrlsToAvg:end,iter(iterI)),3),mean(muAll(:,2,nTrials-nTrlsToAvg:end,iter(iterI)),3)];
        
        
        
        % compute sse with respect to all test data points
        distTrl=[];
        for iClus = 1:size(clusMu,1)
            distTrl(:,iClus)=sum([clusMu(iClus,1,iterI,iTestVal)-dataPtsTest(:,1), clusMu(iClus,2,iterI,iTestVal)-dataPtsTest(:,2)].^2,2);
        end
        [indValsTrl, indTmp]=min(distTrl,[],2); % find which clusters are points closest to
        for iClus = 1:size(clusMu,1)
            sse(iClus)=sum(sum([clusMu(iClus,1,iterI,iTestVal)-dataPtsTest(indTmp==iClus,1), clusMu(iClus,2,iterI,iTestVal)-dataPtsTest(indTmp==iClus,2)].^2,2)); %distance from each cluster from training set to datapoints closest to that cluster
        end
        tsse(iterI,iTestVal)=sum(sse);
        
        %compute 'spreaded-ness' - variance of SE across clusters is a measure
        %of this, assuing uniform data points
        devAvgSSE(:,iterI,iTestVal)   = sse-mean(sse);
        stdAcrossClus(iterI,iTestVal) = std(devAvgSSE(:,iterI,iTestVal));
        varAcrossClus(iterI,iTestVal) = var(devAvgSSE(:,iterI,iTestVal));
    end
end




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
