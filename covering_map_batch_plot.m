clear all;

% wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
cd(wd);

codeDir = [wd '/code_gridCell'];
saveDir = [wd '/data_gridCell'];
addpath(codeDir); addpath(saveDir);
addpath(genpath([codeDir '/gridSCORE_packed']));



% load
nSet        = 6;
gaussSmooth = 1; 
fixBatchSize = 1; %fixed batch size or depend on nClus (for fname)

% clus2run = [7,8,10,12]; 
% nTrials = 100000; 
% nIter=250;
% batchSizeVals = [1, 50, 100, 200, 500];
% epsMuVals=[.01, .05, .075, .1, .2, .3];% learning rate / starting learning rate 


% sims 1 - one single LARGE batchSize, large nTrials, all nClus conds
clus2run = 3:30;%[7,8,10,12]; 
nTrials=10000000;
batchSizeVals=1000;
nIter=200;
epsMuVals=.075;

% sims 2 - a smaller val of trials; testing batch sizes (works fine) - also
% have some sims with ntrials = 5000000 (less batchSizeVals)
% clus2run = 20;
% nTrials=2500000;
% batchSizeVals=[13, 25, 83, 125, 167, 250, 333, 500, 1000, 2000];

% new 3 - fixed batch sizes across clusters
nTrials=2500000;
% clus2run = [18:2:30]; 
% batchSizeVals=[333, 500, 1000];

% new 4 - batchSizes based on mean updates per clus per batch (avgBatch)
fixBatchSize = 0;
clus2run = [18:2:30]; 
% clus2run = [18 20 26]; 
batchSizeVals=[10, 25, 35, 50]; %avgBatchVals

%load loop
for iClus2run = 1:length(clus2run) 
    nClus = clus2run(iClus2run);
    for iEps = 1:length(epsMuVals) 
        epsMuOrig=epsMuVals(iEps);
        epsMuOrig1000=epsMuOrig*1000;
        for iBvals = 1:length(batchSizeVals)
            
            %fixed batch
            if fixBatchSize
                batchSize = batchSizeVals(iBvals);
                fprintf('Loading nClus=%d, epsMu=%d, batchSize=%d\n',nClus,epsMuOrig1000,batchSize)
                fname = [saveDir, sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters*',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter)];
            else %avgBatch
                avgBatch = batchSizeVals(iBvals);
                fprintf('Loading nClus=%d, epsMu=%d, avgBatchSize=%d\n',nClus,epsMuOrig1000,avgBatch)
                fname = [saveDir, sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_avgBatch%d_batchSiz*_%diters*',nClus,round(nTrials/1000),epsMuOrig1000,avgBatch,nIter)];
%                 fname = [saveDir, sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_avgBatch%d_batchSiz*_%diters_useSameTrls*',nClus,round(nTrials/1000),epsMuOrig1000,avgBatch,nIter)];
            end
            
            %edit if want to load more than one file per sim, merge
            f = dir(fname); filesToLoad = cell(1,length(f));
            for iF = 1%:length(f)
                filesToLoad{iF} = f(iF).name;
                load(f(iF).name);
            end
            
            for iterI = 1:nIter
                for iSet=1:nSet
                    densityPlotAll(:,:,iSet,iterI,iEps,iBvals,iClus2run) = imgaussfilt(densityPlot(:,:,iSet,iterI),gaussSmooth);
                end
            end
            %organise gridness values (allen vs willis method)
            gA_gAll(:,:,iEps,iBvals,iClus2run)   = gA(:,:,1);
            gA_oAll(:,:,iEps,iBvals,iClus2run)   = gA(:,:,2);
            gA_radAll(:,:,iEps,iBvals,iClus2run) = gA(:,:,3);
            gA_wavAll(:,:,iEps,iBvals,iClus2run) = gA(:,:,4);
            gW_gAll(:,:,iEps,iBvals,iClus2run) = gW(:,:,1);
            gW_oAll(:,:,iEps,iBvals,iClus2run) = gW(:,:,2);
            gW_radAll(:,:,iEps,iBvals,iClus2run) = gW(:,:,3);
            gW_wavAll(:,:,iEps,iBvals,iClus2run) = gW(:,:,4);
        end
    end
end
%% plot univar scatters

gridMsrType = 'a'; % 'a' or 'w' for allen or willis method

gridMeasure = 'grid';

switch gridMsrType
    case 'a'
        gridness    = gA_gAll;
        orientation = gA_oAll;
        rad         = gA_radAll;
        wav         = gA_wavAll;
    case 'w'
        gridness    = gW_gAll;
        orientation = gW_oAll;
        rad         = gW_radAll;
        wav         = gW_wavAll;
end
switch gridMeasure
    case 'grid'
        datTmp=gridness;
    case 'angle'
        datTmp=orientation;
    case 'rad'
        datTmp=rad;
    case 'wav'
        datTmp=wav;
end


iSet=6;

%plot univar scatters - over clusters (e.g. one batch size)
if size(gridness,4)==1 %only 1 batchSize
    figure; hold on;
    dat1     = squeeze(datTmp(iSet,:,:,:,:));
    barpos  = .25:.5:.5*size(dat1,2);
    colors  = distinguishable_colors(size(dat1,2));
    colgrey = [.5, .5, .5];
    mu      = mean(dat1,1);
    sm      = std(dat1)./sqrt(size(dat1,1));
    ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
    plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
    errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
    % scatter(barpos,mu,750,colors,'x');
    % scatter(barpos,mu,750,colors,'.');
    scatter(barpos,mu,100,colors,'d','filled');
    xlim([barpos(1)-.5, barpos(end)+.5]);
    ylim([-.5,1.25]);
    title(sprintf('%s - eps=%d',gridMeasure,epsMuVals(iEps)*1000))
    
    %plot hist
    figure; pltCount=1;
    iToPlot = 2:length(clus2run);
    for iClus2Run = iToPlot
        subplot(3,3,pltCount)
        hist(dat1(:,iClus2Run),40); %50
        xlim([-.5,1.25]);
        title(sprintf('%d',iClus2Run+2));
        if mod(pltCount,9)==0 && iClus2Run~=iToPlot(end)
            pltCount=1; figure;
        else
            pltCount = pltCount+1;
        end
    end
end


% %plot univar scatters - comparing batchSizeVals, with clusters in separate plots

if size(gridness,4)>1 % more than 1 batchSize
    
figure; hold on;
for iClus2Run = 1:length(clus2run)
    % comparing batch vals, with cluster nums in subplots; can edit for eps
%     figure; hold on;
    subplot(3,3,iClus2Run)
    for iEps = 1:length(epsMuVals)
        %     subplot(2,3,iEps);
        dat1     = squeeze(datTmp(iSet,:,iEps,:,iClus2Run));
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
        title(sprintf('%s - eps=%d, nClus=%d',gridMeasure,epsMuVals(iEps)*1000,clus2run(iClus2Run)))
    end
   
    
    % % comparing learning rate, with batch vals in subplots
    % figure; hold on;
    % for iBvals = 1:length(batchSizeVals)
    %     subplot(2,3,iBvals);
    %     dat1     = squeeze(datTmp(iSet,:,:,iBvals));
    %     barpos  = .25:.5:.5*size(dat1,2);
    %     colors  = distinguishable_colors(size(dat1,2));
    %     colgrey = [.5, .5, .5];
    %     mu      = mean(dat1,1);
    %     sm      = std(dat1)./sqrt(size(dat1,1));
    %     ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
    %     plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
    %     errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
    %     scatter(barpos,mu,750,colors,'x');
    %     xlim([barpos(1)-.5, barpos(end)+.5]);
    %     ylim([-.5,1.25]);
    %     title(sprintf('%s - batchVal=%d',gridMeasure,batchSizeVals(iBvals)))
    % end
end

end





%testing one sim only

% dat1=squeeze(gA(iSet,:,1));
% figure;
% barpos  = .25:.5:.5*size(dat1,1);
% colors  = distinguishable_colors(size(dat1,1));
% colgrey = [.5, .5, .5];
% mu      = mean(dat1);
% sm      = std(dat1)./sqrt(size(dat1,2));
% ci      = sm.*tinv(.025,size(dat1,2)-1); %compute conf intervals
% plotSpread(dat1','xValues',barpos);
% errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
% scatter(barpos,mu,750,'x');
% xlim([barpos(1)-.5, barpos(end)+.5]);
% ylim([-.5,1.25]);
%% density plots

iSet = 5;


for iClus2run = 1:length(clus2run) 
    
    for iEps = 3%1:length(epsMuVals) 
        
        for iBvals = 1:length(batchSizeVals)
            figure;
            subplot(1,3,1);
            imagesc(densityPlotAll(:,:,iSet,iterI,iEps,iBvals,iClus2run));
            aCorrMap=ndautoCORR(densityPlotAll(:,:,iSet,iterI,iEps,iBvals,iClus2run)); %autocorrelogram
            subplot(1,3,2);
            imagesc(aCorrMap,[-.45 .45]);
            subplot(1,3,3);
            [g,gdataA] = gridSCORE(aCorrMap,'allen',1);
            % [g,gdataW] = gridSCORE(aCorrMap,'wills',1);

        end
    end
end





