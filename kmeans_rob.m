%% k means

clear all;

wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
% wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
cd(wd);

codeDir = [wd '/code_gridCell'];
saveDir = [wd '/data_gridCell'];
addpath(codeDir); addpath(saveDir);
addpath(genpath([wd '/gridSCORE_packed']));

% kVals=[7, 9, 15, 20, 25, 30];
kVals = 3:30;%:10;
kVals = [5, 7, 12];
nKvals = length(kVals);

dat = 'rand'; %'rand' points in a box, or 'cat'

nKmeans = 1000;  % run k means n times %1000

nPoints = 10000;  %how many locations - atm not used for 'rand'

%add sth so can run with all unique locs vs random sample?



locRange  = [0, 49];%.9; from -locRange to locRange
spacing   = linspace(locRange(1),locRange(2),locRange(2)+1); 
gaussSmooth=1; % smoothing for density plot / autocorrelogram 


% does it matter how many points there are if all the same points? e.g.
% same if just have each trialsUnique twice/x10? - i think not
% dataPts = repmat(dataPts,50,1);

% load([saveDir '/randTrialsBox_40k']); %load in same data with same trial sequence so same for each sim
% dataPts = trials;

% % normal k means / 'cat learning'
nCats = 2;
sigmaG = [3 0; 0 3]; R = chol(sigmaG);    % isotropic

switch dat
    case 'rand'
        %all unique points in box
        load([saveDir '/randTrialsBox_trialsUnique']);
        dataPts = trialsUnique;
        
        %uniformly sample the box
%         dataPts = [randsample(linspace(locRange(1),locRange(2),50),nPoints,'true'); randsample(linspace(locRange(1),locRange(2),50),nPoints,'true')]';

    case 'cat'
        % draw points from 2 categories (gaussian) from a 2D feature space
        nPointsCat = floor(nPoints/nCats); % points to sample
        for iCat = 1:nCats
            mu(iCat,:)=randsample(locRange(1)+10:locRange(2)-10,2,'true'); % ±10 so category centres are not on the edge
            datPtsGauss(:,:,iCat) = round(repmat(mu(iCat,:),nPointsCat,1) + randn(nPointsCat,2)*R); % key - these are the coordinates of the points
        end
        dataPts = reshape(datPtsGauss,nPoints,2);
        dataPts = dataPts(randperm(length(dataPts)),:);
end

%% k means

nUpdSteps   =  30;    % update steps in the k means algorithm - 40 for random init, 25 for forgy; kmeans++ 20 fine, 22 safe

densityPlotCentres = zeros(length(spacing),length(spacing),nKmeans,nKvals);
for iKvals = 1:nKvals    
nK=kVals(iKvals);
muAll=nan(nK,2,nKmeans);
tsseAll = nan(1,nKmeans);
fprintf('nK = %d \n',nK);
tic
for kMeansIter=1:nKmeans
    if mod(kMeansIter,200)==0
        fprintf('Running k means iteration %d \n',kMeansIter);
    end
    mu   = nan(nK,2,nUpdSteps+1);
    densityPlotClus = zeros(length(spacing),length(spacing),nK);
    for i = 1:nK
        switch dat
            case 'rand'
                mu(:,:,1) = kmplusInit(dataPts,nK); %k means ++ initiatialization
            case 'cat'
                mu(i,:,1) = dataPts(randi(length(dataPts)),:);   %intiate each cluster with one data point - Forgy method
%                 mu(i,:,1) = round(locRange(1) + locRange(2).*rand(1,2));  %initiate clusters with a random point in the box
        end
    end
    
    %run kmeans
    [muEnd,tsse] = kmeans_rm(mu,dataPts,nK,nUpdSteps);
    muAll(:,:,kMeansIter)=muEnd;
    tsseAll(kMeansIter)=tsse;

    %compute gridness
    if strcmp(dat,'rand') %if cat learning, no need to compute gridness
        for iClus=1:nK
            clusTmp  = squeeze(round(muAll(iClus,:,kMeansIter)))';
            for iTrlUpd=1:size(clusTmp,2)
                densityPlotClus(clusTmp(1,iTrlUpd),clusTmp(2,iTrlUpd),iClus) = densityPlotClus(clusTmp(1,iTrlUpd),clusTmp(2,iTrlUpd),iClus)+1;
            end
        end
        %make combined (grid cell) plot, smooth
        densityPlotCentres(:,:,kMeansIter,iKvals) = sum(densityPlotClus,3); 
        densityPlotCentresSm = imgaussfilt(densityPlotCentres(:,:,kMeansIter,iKvals),gaussSmooth);
        aCorrMap=ndautoCORR(densityPlotCentresSm); %autocorrelogram
        [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
        [g,gdataW] = gridSCORE(aCorrMap,'wills',0);
        gA(kMeansIter,:,iKvals) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius];
        gW(kMeansIter,:,iKvals) = [gdataW.g_score, gdataW.orientation, gdataW.wavelength, gdataW.radius];
    end
end
toc
    muAllkVals{iKvals}=muAll; %need this since number of k increases; see if better way to code this
    tssekVals(iKvals,:)=tsseAll;
    
    %compute sse
    [indVal, indSSE] = sort(tsseAll);
    [y, indSSEtmp] = sort(indSSE);
    indSSE2(iKvals,:) = indSSEtmp;
end

% save('kmeans_nK_3-30_uniquePts','muAllkVals','tssekVals', 'gA','gW','densityPlotCentres')
% save('kmeans_nK_3-30_randomPts','muAllkVals','tssekVals', 'gA','gW','densityPlotCentres')

%need?
% [indVal, indSSE(iKvals)] = sort(tsse);
% [y, indSSE2(iKvals)] = sort(indSSE(iKvals));
%%

gridMsr = 'a'; % 'a' or 'w' for allen or willis method

switch gridMsr
    case 'a'
        g = gA;
    case 'w'
        g = gW;
end
%plot hist and density plots

%plot hist
figure; pltCount=1;
for iK = 1:nKvals
    k = kVals(iK);
    
    subplot(3,2,pltCount)
    hist(squeeze(g(:,1,iK)),50);
    xlim([-.5,1.25]);
    
    pltCount = pltCount+1;
    if mod(iK,6)==0
    pltCount=1; figure;
    end
end

%plot univar scatters
figure;
dat     = squeeze(g(:,1,:));
barpos  = .25:.5:.5*size(dat,2);
colors  = distinguishable_colors(size(dat,2));
colgrey = [.5, .5, .5];
mu      = mean(dat,1);
sm      = std(dat)./sqrt(size(dat,1));
ci      = sm.*tinv(.025,size(dat,1)-1); %compute conf intervals
plotSpread(dat,'xValues',barpos,'distributionColors',colors);
errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
scatter(barpos,mu,750,colors,'x');
xlim([barpos(1)-.5, barpos(end)+.5]);
ylim([-.5,1.25]);

% gridness corr with tsse
% figure;
% scatter(g(:,1,iK),tssekVals(iK,:))
% [r p] = corr(g(:,1),tssekVals(iK,:)')
% [r p] = corr(g(:,1),tssekVals(iK,:)','type','spearman')


%plot some grids




%%
% Crossvalidation on clusters from k means, assess error on hex / sq maps
if 1
    
nDataSets = 10;
for iDataSets = 1:nDataSets
    fprintf('xVal dataset %d \n',iDataSets);
    switch dat
        case 'rand' % draw random points in a box (uniform distribution)
            dataPtsTest = [randsample(linspace(locRange(1),locRange(2),50),nPoints,'true'); randsample(linspace(locRange(1),locRange(2),50),nPoints,'true')]'; % random points in a box            
        case 'gauss' % points from clusters of 2D gaussians
            % draw points from 2 categories (gaussian) from a 2D feature space
            nPointsCat = floor(nPoints/nCats); % points to sample
            for iCat = 1:nCats
                mu(iCat,:)=randsample(locRange(1)+10:locRange(2)-10,2,'true'); % ±10 so category centres are not on the edge
                datPtsGauss(:,:,iCat) = round(repmat(mu(iCat,:),nPointsCat,1) + randn(nPointsCat,2)*R); % key - these are the coordinates of the points
            end
            dataPtsTest = reshape(datPtsGauss,nPoints,2);
            dataPtsTest = dataPts(randperm(length(dataPts)),:);
    end
    
    for iKvals = 1:nKvals
    muCurr = muAllkVals{iKvals};
    nK=kVals(iKvals);

    % Crossvalidation start
%     iToXval=[1,2,3,nKmeans-2,nKmeans-1,nKmeans]; %top 3, bottom 3
    iToXval=1:nKmeans; % do all
    for clusMeans=1:length(iToXval)
        %compute distance of existing clusters with new datapoints
        for iClus = 1:nK
            distXval(:,iClus)=sum([muCurr(iClus,1,indSSE2(iKvals,:)==iToXval(clusMeans))-dataPtsTest(:,1), muCurr(iClus,2,indSSE2(iKvals,:)==iToXval(clusMeans))-dataPtsTest(:,2)].^2,2);
        end
        [indValsTest, indTest]=min(distXval,[],2); % find which clusters are points closest to
        
        for iClus = 1:nK
            sseXval(iClus)=sum(sum([(muCurr(iClus,1,indSSE2(iKvals,:)==iToXval(clusMeans)))-dataPtsTest(indTest==iClus,1), (muCurr(iClus,2,indSSE2(iKvals,:)==iToXval(clusMeans)))-dataPtsTest(indTest==iClus,2)].^2,2)); %distance from each cluster from training set to datapoints closest to that cluster
        end
        tsseXval(clusMeans,iDataSets,iKvals)=sum(sseXval);
    end
    %save top 3 bottom 3
    bestWorst3=[1,2,3,nKmeans-2,nKmeans-1,nKmeans];    
    muBest = nan(nK,2,length(bestWorst3),nKvals);
    for iterI=1:length(bestWorst3)
        muBest(:,:,iterI,iKvals) = muCurr(:,:,indSSE2(iKvals,:)==bestWorst3(iterI));
    end
    end
    
end

end



% 
% %save
% saveDat=0;
% if saveDat
% % fname=[savDir, sprintf('/kmeans_clus_dat_k_%d_datPts%dk',k,nPoints/1000)];
% % save(fname,'muAllBest','k','hexPts','nPoints','nKmeans','nDataSets','tsse','indSSE') %'sqPts','
% end



%% plots

saveplots=0;

colgrey = [.5, .5, .5];

% iToPlot=[1,2,3,nKmeans-2,nKmeans-1,nKmeans];

iK = 3;
nK=kVals(iK);
colors = distinguishable_colors(nK); %function for making distinguishable colors for plotting
muAll = muAllkVals{iK};

figure;
for i = 1:6
    subplot(2,3,i); hold on;

    voronoi(muAll(:,1,i),muAll(:,2,i),'k');

    for iClus = 1:nK
        plot(muAll(iClus,1,i),muAll(iClus,2,i),'.','Color',colors(iClus,:),'MarkerSize',25); hold on; %plot cluster final point
%         plot(muAll(iClus,1,(indSSE2==iToPlot(i))),muAll(iClus,2,(indSSE2==iToPlot(i))),'.','Color',colors(iClus,:),'MarkerSize',25); hold on; %plot cluster final point
    end
    xlim([locRange(1),locRange(2)]); ylim([locRange(1),locRange(2)]);
    hold on;
    if i==2
        title(sprintf('Lowest 3 and Highest 3 SSE cluster locations (k=%d)',nK));
    end
end
if saveplots
    fname=[wd, sprintf('/kmeans_clus_k_%d_locs_top_bottom_3_datPts%dk',nK,nPoints/1000)];
    saveas(gcf,fname,'png');
end

%%
% SSE on the test dataset, ranked from the training data set
% note: tsseXval is already ordered from low to high SSE on train set

% SSE sorted on training data
figure; plot((tsse(indSSE))); title(sprintf('SSE over diff initializations sorted (training data) (k=%d)',nK));
if saveplots
    fname=[wd, sprintf('/kmeans_clus_k_%d_sse_train_datPts%dk',nK,nPoints/1000)];
%     fname = [fname '_limDatPtsTo85'];
    saveas(gcf,fname,'png');
end

%on test data
%need to add: for iKvals = 1:nKvals..
figure;
plot(tsseXval(:,:,iKvals)','Color',[.75 .75 .75]); hold on;
% plot(mean(tsseXval(:,:,iKvals)',2),'k','LineWidth',5);
title(sprintf('SSE on Test data sorted by training data (k=%d) (nTest=%d)',nK,nDataSets));
if saveplots
    fname=[wd, sprintf('/kmeans_clus_k_%d_sse_test_sortedByTrain_datPts%dk',nK,nPoints/1000)];
    saveas(gcf,fname,'png');
end

% cross-validation - plot SSE over datasets

figure; hold on;
plot(tsseXval','Color',[.75 .75 .75]);
plot(tsseXval(min(mean(tsseXval,2))==mean(tsseXval,2),:),'Color',[.5 .5 .5],'LineWidth',2); %lowest mean SSE
plot(tsseXval(1,:),'Color',[.4 .4 .4],'LineWidth',2); %lowest mean SSE from training data
% plot([tsseSq;tsseHex]');

title(sprintf('SSE over all Test data sets (k=%d) (nTest=%d)(datPts=%dk)',nK,nDataSets,nPoints/1000));
if saveplots
    fname=[wd, sprintf('/kmeans_clus_k_%d_sse_test_overDatasets_datPts%dk',nK,nPoints/1000)];
%     fname = [fname '_limDatPtsTo85'];
    saveas(gcf,fname,'png');
end


% figure;hist([tsseHex;tsseXval(min(mean(tsseXval,2))==mean(tsseXval,2),:)]',20);



%%
gaussSmooth=1;

nSets=6; %top and bottom 3 SSE
spacing=linspace(locRange(1),locRange(2),locRange(2)+1); 
densityPlotClus      = zeros(length(spacing),length(spacing),nK,nSets);

for iSet=1:nSets
figure; hold on;
for iClus=1:nK
    clusTmp  = squeeze(round(muBest(iClus,:,iSet)))';
    for iTrlUpd=1:size(clusTmp,2)
        densityPlotClus(clusTmp(1,iTrlUpd),clusTmp(2,iTrlUpd),iClus,iSet) = densityPlotClus(clusTmp(1,iTrlUpd),clusTmp(2,iTrlUpd),iClus,iSet)+1;
    end
end

%make combined (grid cell) plot, smooth
densityPlotCentres(:,:,iSet) = sum(densityPlotClus(:,:,:,iSet),3);
densityPlotCentresSm = imgaussfilt(densityPlotCentres(:,:,iSet),gaussSmooth);

subplot(1,3,1);
imagesc(densityPlotCentresSm);


aCorrMap=ndautoCORR(densityPlotCentresSm); %autocorrelogram
subplot(1,3,2);
imagesc(aCorrMap,[-.45 .45]);

subplot(1,3,3);
[g,gdataA] = gridSCORE(aCorrMap,'allen',1);
% [g,gdataW] = gridSCORE(aCorrMap,'wills',1);
    
end


