%% Plotting script 4 - quickPlot: run a quick sim (e.g. nIter=1) and plot results
% needs muAll and trials for some (see main script), and actOverTime for some 

spacing=linspace(locRange(1),locRange(2),locRange(2)+1)+1;
gaussSmooth = 1; %smoothing for density map
imageFilter = fspecial('gaussian',5,gaussSmooth); % param for smoothing (for nanconv - smoothing that deals with nans; this param is default for imgaussfilt)

%% plot centres at end of learning (requires muAll)

colors = distinguishable_colors(nClus); %function for making distinguishable colors for plotting
for iSet=1
    figure; hold on;
    scatter(muAll(:,1,end),muAll(:,2,end),20e+2,colors,'.');
    xlim(locRange); ylim(locRange);
end

%% plot autocorrelogram, gridness over time

for iSet=1:21
    for iterI = 1
        % densityPlotCentresSm = nanconv(densityPlot(:,:,iSet,iterI),imageFilter, 'nanout'); % just takes cluster positions and smooths
        densityPlotCentresSm = densityPlotActNorm(:,:,iSet,iterI);
        figure; hold on;
        subplot(1,2,1)
        imagesc(densityPlotCentresSm);
        aCorrMap=ndautoCORR(densityPlotCentresSm); %autocorrelogram
        subplot(1,2,2)
        [g,gdataA] = gridSCORE(aCorrMap,'allen',1);
    end
end

%% gridness - left vs right half

%trapz split into 2 equal areas - trapzKrupic only for now]
h = 50; hLeft  = 17;
hRight = 33;% - 33 = start from 18 from left % 27;
    
for iSet=1:21
    for iterI = 1
        %     densityPlotCentresSm = nanconv(densityPlot(:,:,iSet,iterI),imageFilter, 'nanout'); % just takes cluster positions and smooths
        densityPlotCentresSm = densityPlotActNorm(:,:,iSet,iterI);
        figure; hold on;
        subplot(1,3,1)
        imagesc(densityPlotCentresSm);
        %left half
        subplot(1,3,2)
        aCorrMap = ndautoCORR(densityPlotCentresSm(:,1:hLeft));
        [g,gdataA] = gridSCORE(aCorrMap,'allen',1);
        %right half
        subplot(1,3,3)
        aCorrMap = ndautoCORR(densityPlotCentresSm(:,h-hRight+1:end));
        [g,gdataA] = gridSCORE(aCorrMap,'allen',1);
    end
end
%% plot over time - single plot

if ~exist('dat')
    dat=dat2;
end

%if plotAgent/trials, need trials as an output argument
plotAgent = 0;

nClus = size(muAll,1);

colors = distinguishable_colors(nClus); %function for making distinguishable colors for plotting
colAgent = [.75, .75, .75];
colAgentCh = fliplr(linspace(0,.9,nTrials)); % colAgent gets darker over time 

%if plotting agent over batches - make cluster positions same over trials in a batch
% if plotAgent
muTrls = nan(nClus,2,nTrials/4); %nTrials/2 for trapzKfrmSq
for iBatch = 1:nBatches % nBatches/4 for trapzKfrmSq
    muTrls(:,1,batchSize*(iBatch-1)+1:batchSize*(iBatch-1)+batchSize)=repmat(muAll(:,1,iBatch),1,batchSize);
    muTrls(:,2,batchSize*(iBatch-1)+1:batchSize*(iBatch-1)+batchSize)=repmat(muAll(:,2,iBatch),1,batchSize);
end
% end

figure; hold on;
xlim([0, 50]);
ylim([0,50]); 
clear h1
for iTrl = 1000:nTrials
    if mod(iTrl,500)==0 %plot centers after x trials
        % plot agent
        if plotAgent
            %have to start from 1000 trials min
            plotTrls=iTrl-999:iTrl;
            if strcmp(dat(1:3),'cat')
                h1(iTrl)=plot(trials(plotTrls,1),trials(plotTrls,2),'.','Color',repmat(colAgentCh(iTrl),1,3)); hold on;
            else
                h1(iTrl)=plot(trials(plotTrls,1),trials(plotTrls,2),'Color',repmat(colAgentCh(iTrl),1,3)); hold on;
            end
            if mod(iTrl,5000)==0
                delete(h1);
            end
        end
        % plot clusters
        scatter(squeeze(muTrls(:,1,iTrl)),squeeze(muTrls(:,2,iTrl)),1000,colors,'.'); hold on;
        drawnow;
    end
end
%% plot over time - subplots (need muAll as output arg)

savePlots=0;

iterI = 1;

nClus = size(muAll,1);
colors = distinguishable_colors(nClus); %function for making distinguishable colors for plotting
figure('units','normalized','outerposition',[0 0 1 1]);
iPlot = 1; subplot(3,4,iPlot); hold on;%subplot
for iTrl = 1:nBatches
    if mod(iTrl,500)==0
        iPlot=iPlot+1;
        subplot(3,4,iPlot); hold on;
    end
    xlim(locRange); ylim(locRange);
    if mod(iTrl,500)==0 || iTrl == 1 %plot centers after x trials
        for i=1:nClus
            plot(squeeze(muAll(i,1,iTrl,iterI)),squeeze(muAll(i,2,iTrl,iterI)),'.','Color',colors(i,:),'MarkerSize',20); hold on;
        end
        drawnow;
    end
end

if savePlots
    fname=[saveDir, sprintf('/covering_map_%d_clus_locs_plotOverTime_top_%d_trls_eps%d_alpha%d_%diters',nClus,nTrials,epsMuOrig100,alpha10,nIter)];
    if warpBox
        fname = [fname sprintf('_warpBox_resetEps%d',resetEps)];
    end
    saveas(gcf,fname,'png');
end

%% plot navigation example - over 3 time points

% allows saving plot over 3 time points
% needs muAll and trials 
% needs to be done one (nIter=1) at a time - since muAll is not saved across iters

savePlots = 0;

figDir = [wd '/data_gridCell/figs'];
fntSiz = 15;
gaussSmooth=1;

iterI     = 1;
prtTrls = [.0005, .05, .99];
trls2plt = {1:prtTrls(1)*nTrials, prtTrls(2)*nTrials:(prtTrls(2)*nTrials)+2500, prtTrls(3)*nTrials+1:nTrials}; %1:500, 50000:52500 (2.5k trials; 5% to 5.25%), then 990000:nTrials (10k trials; 99%-100%)

%if plotAgent, need trials as an output argument
plotAgent = 1;
nClus = size(muAll,1);
colAgent = [.75, .75, .75];
colGreyClus = [.85, .85, .85];
colors(:,:,1) = repmat(colGreyClus,nClus,1);
colors(:,:,2) = distinguishable_colors(nClus)+(1-distinguishable_colors(nClus)).*.5; %lighter
colors(:,:,3) = distinguishable_colors(nClus);

%if plotting agent over batches - need make cluster positions same over trials in a batch
muTrls = nan(nClus,2,nTrials); %nTrials/2 - for trapzKfrmSq
for iBatch = 1:nBatches
    muTrls(:,1,batchSize*(iBatch-1)+1:batchSize*(iBatch-1)+batchSize)=repmat(muAll(:,1,iBatch),1,batchSize);
    muTrls(:,2,batchSize*(iBatch-1)+1:batchSize*(iBatch-1)+batchSize)=repmat(muAll(:,2,iBatch),1,batchSize);
end
cTime=datestr(now,'HHMMSS'); %so each iter has the same ending

for iPlot = 1:length(trls2plt)
    figure; hold on;
    plot(trials(trls2plt{iPlot},1),trials(trls2plt{iPlot},2),'Color',colAgent); hold on;
    scatter(squeeze(muTrls(:,1,trls2plt{iPlot}(end))),squeeze(muTrls(:,2,trls2plt{iPlot}(end))),1000,colors(:,:,iPlot),'.'); hold on;    
    xlim([locRange(1) locRange(2)+1]);
    ylim([locRange(1) locRange(2)+1]);
    xticks([]); xticklabels({''}); yticks([]); yticklabels({''});
    if iPlot==3
        densityPlotCentresSm = imgaussfilt(densityPlotActNorm(:,:,end,iterI),gaussSmooth);
        aCorrMap=ndautoCORR(densityPlotCentresSm); %autocorrelogram
        [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
        title(sprintf('Grid Score %.2f',g));
        set(gca,'FontSize',fntSiz,'fontname','Arial');
    end
    fname = [figDir, sprintf('/navExamples_%s_eps%d_batchSiz%d_nClus%d_iter%d_%s_t%d',dat,epsMuOrig1000,batchSizeVals(iBvals),clus2run(1),iterI,cTime,iPlot)];
    if savePlots
        set(gcf,'Renderer','painters');
        print(gcf,'-depsc2',fname)
        saveas(gcf,fname,'png');
        close all
    end
end