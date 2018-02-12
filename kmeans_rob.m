%% k means - hexagonal maps?

clear all;

% wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
cd(wd);

codeDir = [wd '/code_gridCell'];
saveDir = [wd '/data_gridCell'];
addpath(codeDir); addpath(saveDir);
addpath(genpath([wd '/gridSCORE_packed']));

% kVals=[7, 9, 15, 20, 25, 30];
kVals = 3:30;%:10;
kVals = [10, 20, 30];

dat = 'rand'; %'rand' points in a box, or 'cat'

nKmeans = 1000;  % run k means n times %1000


nPoints = 10000;  %how many locations - atm not used for 'rand'

locRange = [0, 49];%.9; from -locRange to locRange

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

for iK = 1:length(kVals)    
k=kVals(iK);
muAll=nan(k,2,nKmeans);
tsseAll = nan(1,nKmeans);
fprintf('Running k means nK = %d \n',k);
tic
for kMeansIter=1:nKmeans
    if mod(kMeansIter,100)==0
        fprintf('iteration %d \n',kMeansIter);
    end
    mu   = nan(k,2,nUpdSteps+1);
    for i = 1:k
        switch dat
            case 'rand'
                %k means ++ initiatialization
                mu(:,:,1) = kmplusInit(dataPts,k);
            case 'cat'
                mu(i,:,1) = dataPts(randi(length(dataPts)),:);   %intiate each cluster with one data point - Forgy method
%                 mu(i,:,1) = round(locRange(1) + locRange(2).*rand(1,2));  %initiate clusters with a random point in the box
        end
    end
    
    %run kmeans
    [muEnd,tsse] = kmeans_rm(mu,dataPts,nK,nUpdSteps);
    muAll(:,:,kMeansIter)=muEnd;
    tsseAll(kMeansIter)=tsse;
end
toc
    muAllkVals{iK}=muAll; %need this since number of k increases; see if better way to code this
    tssekVals(iK,:)=tsseAll;
%need?
% [indVal, indSSE] = sort(tsse);
% [y, indSSE2] = sort(indSSE);

end




%%
gaussSmooth=1;

spacing         = linspace(locRange(1),locRange(2),locRange(2)+1); 
densityPlotClus = zeros(length(spacing),length(spacing),k,nKmeans);

for iK=1:length(kVals)
    k=kVals(iK);
    muAll=muAllkVals{iK};
    fprintf('Processing k means gridness: nK = %d \n',k);
    
    tic
    for kMeansIter=1:nKmeans
        if mod(kMeansIter,100)==0
            fprintf('iteration %d \n',kMeansIter);
        end
        for iClus=1:k
            clusTmp  = squeeze(round(muAll(iClus,:,kMeansIter)))';
            for iTrlUpd=1:size(clusTmp,2)
                densityPlotClus(clusTmp(1,iTrlUpd),clusTmp(2,iTrlUpd),iClus,kMeansIter) = densityPlotClus(clusTmp(1,iTrlUpd),clusTmp(2,iTrlUpd),iClus,kMeansIter)+1;
            end
        end
        
        %make combined (grid cell) plot, smooth
        densityPlotCentres(:,:,kMeansIter) = sum(densityPlotClus(:,:,:,kMeansIter),3);
        densityPlotCentresSm = imgaussfilt(densityPlotCentres(:,:,kMeansIter),gaussSmooth);
        
        aCorrMap=ndautoCORR(densityPlotCentresSm); %autocorrelogram
        
        [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
        [g,gdataW] = gridSCORE(aCorrMap,'wills',0);
        gA(kMeansIter,:,iK) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius];
        gW(kMeansIter,:,iK) = [gdataW.g_score, gdataW.orientation, gdataW.wavelength, gdataW.radius];
        
    end
    toc
end


%%
% Crossvalidation on clusters from k means, assess error on hex / sq maps
if 0
    
nDataSets = 100;
for iDataSets = 1:nDataSets
    switch dat
        case 'rand' % draw random points in a box (uniform distribution)
            dataPtsTest = [randsample(linspace(locRange(1),locRange(2),50),nPoints,'true'); randsample(linspace(locRange(1),locRange(2),50),nPoints,'true')]'; % random points in a box
%             dataPtsTest = [randsample(linspace(-locRange+.15,locRange-.15,101),nPoints,'true'); randsample(linspace(-locRange+.15,locRange-.15,101),nPoints,'true')]'; % random points in a box
            
        case 'gauss' % points from clusters of 2D gaussians
            dataPtsTest = ind2grid(pointsInGauss(randi(length(pointsInGauss),nPoints,1),:));
    end
    
    % Crossvalidation start
%     iToXval=[1,2,3,nKmeans-2,nKmeans-1,nKmeans]; %top 3, bottom 3
    iToXval=1:nKmeans; % do all
    for clusMeans=1:length(iToXval)
        %compute distance of existing clusters with new datapoints
        for iClus = 1:k
            distXval(:,iClus)=sum([muAll(iClus,1,indSSE2==iToXval(clusMeans))-dataPtsTest(:,1), muAll(iClus,2,indSSE2==iToXval(clusMeans))-dataPtsTest(:,2)].^2,2);
        end
        [indValsTest, indTest]=min(distXval,[],2); % find which clusters are points closest to
        
        for iClus = 1:k
            sseXval(iClus)=sum(sum([(muAll(iClus,1,indSSE2==iToXval(clusMeans)))-dataPtsTest(indTest==iClus,1), (muAll(iClus,2,indSSE2==iToXval(clusMeans)))-dataPtsTest(indTest==iClus,2)].^2,2)); %distance from each cluster from training set to datapoints closest to that cluster
        end
        tsseXval2(clusMeans,iDataSets)=sum(sseXval);
    end    
    
end

end


%save top 3 bottom 3
bestWorst3=[1,2,3,nKmeans-2,nKmeans-1,nKmeans];
muAllBest = nan(k,2,length(bestWorst3));
for iterI=1:length(bestWorst3)
    muAllBest(:,:,iterI) = muAll(:,:,indSSE2==bestWorst3(iterI));
end

%save
saveDat=0;
if saveDat
% fname=[savDir, sprintf('/kmeans_clus_dat_k_%d_datPts%dk',k,nPoints/1000)];
% save(fname,'muAllBest','k','hexPts','nPoints','nKmeans','nDataSets','tsse','indSSE') %'sqPts','
end



%% plots

saveplots=0;

colors = distinguishable_colors(k); %function for making distinguishable colors for plotting
colgrey = [.5, .5, .5];

% iToPlot=[1,2,3,nKmeans-2,nKmeans-1,nKmeans];


iK=20;
k = iK;
muAll = muAllkVals{k};



figure;
for i = 1:6
    subplot(2,3,i); hold on;

    voronoi(muAll(:,1,i),muAll(:,2,i),'k');

    for iClus = 1:k
        plot(muAll(iClus,1,i),muAll(iClus,2,i),'.','Color',colors(iClus,:),'MarkerSize',25); hold on; %plot cluster final point
%         plot(muAll(iClus,1,(indSSE2==iToPlot(i))),muAll(iClus,2,(indSSE2==iToPlot(i))),'.','Color',colors(iClus,:),'MarkerSize',25); hold on; %plot cluster final point
    end
    xlim([locRange(1),locRange(2)]); ylim([locRange(1),locRange(2)]);
    hold on;
    if i==2
        title(sprintf('Lowest 3 and Highest 3 SSE cluster locations (k=%d)',k));
    end
end
if saveplots
    fname=[wd, sprintf('/kmeans_clus_k_%d_locs_top_bottom_3_datPts%dk',k,nPoints/1000)];
%     fname = [fname '_limDatPtsTo85'];
    saveas(gcf,fname,'png');
end

%%
% SSE on the test dataset, ranked from the training data set
% note: tsseXval2 is already ordered from low to high SSE on train set

% SSE sorted on training data
figure; plot((tsse(indSSE))); title(sprintf('SSE over diff initializations sorted (training data) (k=%d)',k));
if saveplots
    fname=[wd, sprintf('/kmeans_clus_k_%d_sse_train_datPts%dk',k,nPoints/1000)];
%     fname = [fname '_limDatPtsTo85'];
    saveas(gcf,fname,'png');
end

%on test data
figure;
plot(tsseXval2,'Color',[.75 .75 .75]); hold on;
plot(mean(tsseXval2,2),'k','LineWidth',5);
title(sprintf('SSE on Test data sorted by training data (k=%d) (nTest=%d)',k,nDataSets));
if saveplots
    fname=[wd, sprintf('/kmeans_clus_k_%d_sse_test_sortedByTrain_datPts%dk',k,nPoints/1000)];
    saveas(gcf,fname,'png');
end

% cross-validation - plot SSE over datasets

figure; hold on;
plot(tsseXval2','Color',[.75 .75 .75]);
plot(tsseXval2(min(mean(tsseXval2,2))==mean(tsseXval2,2),:),'Color',[.5 .5 .5],'LineWidth',2); %lowest mean SSE
plot(tsseXval2(1,:),'Color',[.4 .4 .4],'LineWidth',2); %lowest mean SSE from training data
% plot([tsseSq;tsseHex]');

title(sprintf('SSE over all Test data sets (k=%d) (nTest=%d)(datPts=%dk)',k,nDataSets,nPoints/1000));
if saveplots
    fname=[wd, sprintf('/kmeans_clus_k_%d_sse_test_overDatasets_datPts%dk',k,nPoints/1000)];
%     fname = [fname '_limDatPtsTo85'];
    saveas(gcf,fname,'png');
end


% figure;hist([tsseHex;tsseXval2(min(mean(tsseXval2,2))==mean(tsseXval2,2),:)]',20);



%%
gaussSmooth=1;

nSets=6; %top and bottom 3 SSE
spacing=linspace(locRange(1),locRange(2),locRange(2)+1); 
densityPlotClus      = zeros(length(spacing),length(spacing),k,nSets);

for iSet=1:nSets
figure; hold on;
for iClus=1:k
    clusTmp  = squeeze(round(muAllBest(iClus,:,iSet)))';
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


