%% Hexagonal & Square maps
% do hexagonally organised clusters (or squares - or other controls) do the best?

clear all;

wd='/Users/robert.mok/Dropbox/Grid_cell_model';
% wd='/Users/robertmok/Dropbox/Grid_cell_model';
cd(wd);

locRange = 1;%.9; from -locRange to locRange

nPoints = 10000;  %how many locations

dat = 'rand'; %'rand' points in a box, or 'gauss' - n gaussians in a box
% dataPts = [randsample(linspace(-locRange,locRange,101),nPoints,'true'); randsample(linspace(-locRange,locRange,101),nPoints,'true')]'; % random points in a box
dataPts = [randsample(linspace(-locRange+.15,locRange-.15,101),nPoints,'true'); randsample(linspace(-locRange+.15,locRange-.15,101),nPoints,'true')]'; % random points in a box
        
spacingVal = 14;%16; %smaller spacing = more clusters     %for plotting  i did 18,14 10, 8

% %square
sqSpacing=linspace(-locRange,locRange,101/spacingVal);
sqSpacing(end)=[]; %remove and 1
h=diff(sqSpacing(1:2));
sqSpacing=sqSpacing+h/2; %start a bit after border, similar to hex

sqPts = [];
for i=1:length(sqSpacing),
    for j=1:length(sqSpacing),
        sqPts = [sqPts ;[sqSpacing(i),sqSpacing(j)]];
    end
end
sqPts = [];
for i=1:length(sqSpacing),
    for j=1:length(sqSpacing),
        sqPts = [sqPts ;[sqSpacing(i),sqSpacing(j)]];
    end
end

% %hexagonal 
%hex spacing - same as the square's diag (greater spacing than square
%spacing - to show hex is still better. could also do with the square's
% length)
% hexSpac1=sqrt(h^2+h^2);

%same as sq spacing (not diag)
% hexSpac=h;
% 
% hexSpac=hexSpac*.9; % if need fewer/more points for hex 
% hexSpac=hexSpac*1.3;

% % no shift - %old diamonds
% hexSpacing1=-locRange+hexSpac:hexSpac:locRange;
% hexSpacing2=-locRange+hexSpac/2:hexSpac:locRange;

%hexSpace1 will be length between y-axis points
%hexSpace2 is distance between x-axis points

%v1
% hexSpac1 = sqrt(h^2+h^2);
% hexSpac1 = h;
% % hexSpac1 = hexSpac1*1.05;% if need fewer/more points for hex 
% hexSpac2 = (sqrt(hexSpac1^2-(.5*hexSpac1)^2))*2;

%v2
hexSpac2 = sqrt(h^2+h^2);
hexSpac2 = h;
% hexSpac2 = hexSpac2*1.05;% if need fewer/more points for hex 
hexSpac1 = (sqrt(hexSpac2^2-(.5*hexSpac2)^2))*2;

hexSpacing1a=-locRange+hexSpac1:hexSpac1:locRange;
hexSpacing2a=-locRange+hexSpac2:hexSpac2:locRange;
hexSpacing1b=-locRange+hexSpac1/2:hexSpac1:locRange;
hexSpacing2b=-locRange+hexSpac2/2:hexSpac2:locRange;

hexPts = [];
for i=1:length(hexSpacing1a),
    for j=1:length(hexSpacing2a),
        hexPts = [hexPts ;[hexSpacing1a(i),hexSpacing2a(j)]];
    end
end
for i=1:length(hexSpacing1b),
    for j=1:length(hexSpacing2b),
        hexPts = [hexPts ;[hexSpacing1b(i),hexSpacing2b(j)]];
    end
end

figure; 
subplot(1,2,1); voronoi(sqPts(:,1),sqPts(:,2)); title(sprintf('Square, %d clusters',length(sqPts))); xlim([-1, 1]); ylim([-1, 1]);
subplot(1,2,2); voronoi(hexPts(:,1),hexPts(:,2)); title(sprintf('Hexagonal, %d clusters',length(hexPts))); xlim([-1, 1]); ylim([-1, 1]);
if 0
    fname=[wd, sprintf('/kmeans_clus_hex_sq_maps_%dclus',length(hexPts))];
    saveas(gcf,fname,'png');
end



%% Hex vs square vs k means

k=length(hexPts);

nKmeans = 100;%250;  % run k means n times %1000
nUpdSteps   =  22;    % update steps in the k means algorithm - 40 for random init, 25 for forgy; kmeans++ 20 fine, 22 safe

for kMeansIter=1:nKmeans,
    
    if mod(kMeansIter,10)==0,
        fprintf('batch iteration %d \n',kMeansIter);
    end
    
%     mu=nan(k,2,nUpdSteps+1);
    for i = 1:k
        switch dat
            case 'rand'
%                 mu(i,:,1) = -locRange + locRange.*2.*rand(1,2);  %initiate clusters with a random point in the box
%                 mu(i,:,1) = dataPts(randi(length(dataPts)),:);   %intiate each cluster with one data point - Forgy method
                
                %k means ++ initiatialization
                %1 - trial 1) random datapoint as 1st cluster
                %2) find distance of each cluster to all points; assign
                %each datapoint to the closest cluster and keep those
                %distances
                %3) assign probability of each data point as a new cluster
                %(shortest dist from all current clusters lower, further higher)                
                %4) select random cluster from dataset (exclude all current
                % clusters), weighted by the distance
                %5) loop from 2 until k 
                
                if i==1,% random datapoint as 1st cluster
                   mu(i,:,1) = dataPts(randi(length(dataPts)),:); 
                end
                
                if i~=k, % no need update k+1
                    clear distInit
                    for iClus = 1:i,% loop over clusters that exist now
                        distInit(:,iClus)=sum([mu(iClus,1,1)-dataPts(:,1),  mu(iClus,2,1)-dataPts(:,2)].^2,2); %squared euclid for k means
                    end
                    [indValsInit, indInit]=min(distInit,[],2); % find which clusters are points closest to
                    
                    distClus=[];
                    for iClus = 1:i,
                        indOrig(:,i) = indInit==iClus;
                        distClusTmp = sum([(mu(iClus,1,1)-dataPts(indOrig(:,i),1)), (mu(iClus,2,1)-dataPts(indOrig(:,i),2))].^2,2);
                        distClus = [distClus; [distClusTmp, repmat(iClus,length(distClusTmp),1)]];
                    end
                            
                    %need to keep track of the indices of the original dist variable - get the
                    %datapoints that were the farthest from all clusters, get that cluster and see which datapoint that was relative to that cluster (since i just save the distance)
                    
                    distClusNorm = distClus(:,1)./sum(distClus(:,1));
                    distClusPr   = cumsum(distClusNorm(:,1)); %get cumsum, then generate rand val from 0 to 1 and choose smallest val - the larger the dis, the more likely the rand value will lie between it and its previous value in a cumsum plot
                    ind=find(rand(1)<distClusPr,1);% %find smallest value that is larger than the random value (0 to 1 uniform distr)
                                        
%                     tmp(i)=distClus(ind,1); %testing if getting the right values; if plot, see that it should be lower pr for closer values, higher for larger. note if very few high distances, this hist will look normally distributed
                    
                    clusInd = distClus(ind,2); %find which is the closest cluster
                    indDat = find(distInit(:,clusInd)==distClus(ind,1)); %find where the datapoint is in the original vector
                    
                    if size(indDat,1)>1,
                        indDat=indDat(randi(size(indDat,1),1));
                    end
                    mu(i+1,:,1) = dataPts(indDat,:);                    
                end
            case 'gauss'
                mu(i,:,1) = ind2grid(pointsInGauss(randi(length(pointsInGauss),1,1),:)); %iniate each cluster with one data point - Forgy method (kind of - here it's a rand pt in the gauss, should be one of the data points in train set)
        end
    end
    
    % start of iteration loop
    for upd = 1:nUpdSteps,
        for iClus = 1:k,% loop over clusters not dataPts, because nPts can be huge (so vectorise over points
            dist(:,iClus)=sum([mu(iClus,1,upd)-dataPts(:,1),  mu(iClus,2,upd)-dataPts(:,2)].^2,2); %squared euclid for k means
        end
        [indVals, ind]=min(dist,[],2); % find which clusters are points closest to
        for iClus = 1:k,
            if any(isnan([mean(dataPts(ind==iClus,1)), mean(dataPts(ind==iClus,2))]))
                mu(iClus,:,upd+1)=mu(iClus,:,upd);
            else
                mu(iClus,:,upd+1)=[mean(dataPts(ind==iClus,1)), mean(dataPts(ind==iClus,2))];
            end
        end
    end
    
    muAll(:,:,kMeansIter)=mu(:,:,end);
    
    %compute sum of squared errors across clusters
    for iClus = 1:k,
        sse(iClus)=sum(sum([(muAll(iClus,1,kMeansIter))-dataPts(ind==iClus,1), (muAll(iClus,2,kMeansIter))-dataPts(ind==iClus,2)].^2,2));
    end
    tsse(kMeansIter)=sum(sse);
end


[indVal, indSSE] = sort(tsse);
[y, indSSE2] = sort(indSSE);

% Crossvalidation on clusters from k means, assess error on hex / sq maps
nDataSets = 100;
for iDataSets = 1:nDataSets,
    switch dat
        case 'rand' % draw random points in a box (uniform distribution)
            dataPtsTest = [randsample(linspace(-locRange,locRange,101),nPoints,'true'); randsample(linspace(-locRange,locRange,101),nPoints,'true')]'; % random points in a box
%             dataPtsTest = [randsample(linspace(-locRange+.15,locRange-.15,101),nPoints,'true'); randsample(linspace(-locRange+.15,locRange-.15,101),nPoints,'true')]'; % random points in a box
            
        case 'gauss' % points from clusters of 2D gaussians
            dataPtsTest = ind2grid(pointsInGauss(randi(length(pointsInGauss),nPoints,1),:));
    end
    
    % Crossvalidation start
%     iToXval=[1,2,3,nKmeans-2,nKmeans-1,nKmeans]; %top 3, bottom 3
    iToXval=1:nKmeans; % do all
    for clusMeans=1:length(iToXval),
        %compute distance of existing clusters with new datapoints
        for iClus = 1:k,
            distXval(:,iClus)=sum([muAll(iClus,1,indSSE2==iToXval(clusMeans))-dataPtsTest(:,1), muAll(iClus,2,indSSE2==iToXval(clusMeans))-dataPtsTest(:,2)].^2,2);
        end
        [indValsTest, indTest]=min(distXval,[],2); % find which clusters are points closest to
        
        for iClus = 1:k,
            sseXval(iClus)=sum(sum([(muAll(iClus,1,indSSE2==iToXval(clusMeans)))-dataPtsTest(indTest==iClus,1), (muAll(iClus,2,indSSE2==iToXval(clusMeans)))-dataPtsTest(indTest==iClus,2)].^2,2)); %distance from each cluster from training set to datapoints closest to that cluster
        end
        tsseXval2(clusMeans,iDataSets)=sum(sseXval);
    end    
        
    
    %Square
    %compute distance of existing clusters with new datapoints
    for iClus = 1:length(sqPts),
        distSq(:,iClus)=sum([sqPts(iClus,1)-dataPtsTest(:,1), sqPts(iClus,2)-dataPtsTest(:,2)].^2,2);
    end
    [indValsSq, indSq]=min(distSq,[],2); % find which clusters are points closest to
    
    for iClus = 1:length(sqPts), %double check mistakes
        sseSq(iClus)=sum(sum([sqPts(iClus,1)-dataPtsTest(indSq==iClus,1), sqPts(iClus,2)-dataPtsTest(indSq==iClus,2)].^2,2)); %distance from each cluster from training set to datapoints closest to that cluster
    end
    tsseSq(iDataSets)=sum(sseSq);


    %hex
    %compute distance of existing clusters with new datapoints
    for iClus = 1:length(hexPts), 
        distSq(:,iClus)=sum([hexPts(iClus,1)-dataPtsTest(:,1), hexPts(iClus,2)-dataPtsTest(:,2)].^2,2);
    end
    [indValsHex, indHex]=min(distSq,[],2); % find which clusters are points closest to
    
    for iClus = 1:length(hexPts), %double check mistakes
        sseHex(iClus)=sum(sum([hexPts(iClus,1)-dataPtsTest(indHex==iClus,1), hexPts(iClus,2)-dataPtsTest(indHex==iClus,2)].^2,2)); %distance from each cluster from training set to datapoints closest to that cluster
    end
    tsseHex(iDataSets)=sum(sseHex);
    
end

%save top 3 bottom 3
bestWorst3=[1,2,3,nKmeans-2,nKmeans-1,nKmeans];
muAllBest = nan(k,2,length(bestWorst3));
for iterI=1:length(bestWorst3),
    muAllBest(:,:,iterI) = muAll(:,:,indSSE2==bestWorst3(iterI));
end

%save
saveDat=0;
if saveDat,
fname=[wd, sprintf('/kmeans_clus_dat_k_%d_datPts%dk',length(hexPts),nPoints/1000)];
save(fname,'muAllBest','k','hexPts','nPoints','nKmeans','nDataSets','tsse','indSSE') %'sqPts','
end
%% plots

saveplots=0;

colors = distinguishable_colors(k); %function for making distinguishable colors for plotting
colgrey = [.5, .5, .5];

% iToPlot=[1,2,3,nKmeans-2,nKmeans-1,nKmeans];
figure;
for i = 1:6
    subplot(2,3,i); hold on;
    plot(dataPtsTest(:,1),dataPtsTest(:,2),'.','Color',colgrey,'MarkerSize',2); hold on; 
    voronoi(muAllBest(:,1,i),muAllBest(:,2,i),'k');

    for iClus = 1:k
        plot(muAllBest(iClus,1,i),muAllBest(iClus,2,i),'.','Color',colors(iClus,:),'MarkerSize',25); hold on; %plot cluster final point
%         plot(muAll(iClus,1,(indSSE2==iToPlot(i))),muAll(iClus,2,(indSSE2==iToPlot(i))),'.','Color',colors(iClus,:),'MarkerSize',25); hold on; %plot cluster final point
    end
    xlim([-1.1,1.1]); ylim([-1.1,1.1]);
    hold on;
    if i==2,
        title(sprintf('Lowest 3 and Highest 3 SSE cluster locations (k=%d)',length(hexPts)));
    end
end
if saveplots,
    fname=[wd, sprintf('/kmeans_clus_k_%d_locs_top_bottom_3_datPts%dk',length(hexPts),nPoints/1000)];
%     fname = [fname '_limDatPtsTo85'];
    saveas(gcf,fname,'png');
end

% SSE on the test dataset, ranked from the training data set
% note: tsseXval2 is already ordered from low to high SSE on train set

% SSE sorted on training data
figure; plot((tsse(indSSE))); title(sprintf('SSE over diff initializations sorted (training data) (k=%d)',length(hexPts)));
if saveplots,
    fname=[wd, sprintf('/kmeans_clus_k_%d_sse_train_datPts%dk',length(hexPts),nPoints/1000)];
%     fname = [fname '_limDatPtsTo85'];
    saveas(gcf,fname,'png');
end

%on test data
figure;
plot(tsseXval2,'Color',[.75 .75 .75]); hold on;
plot(mean(tsseXval2,2),'k','LineWidth',5);
title(sprintf('SSE on Test data sorted by training data (k=%d) (nTest=%d)',length(hexPts),nDataSets));
if saveplots,
    fname=[wd, sprintf('/kmeans_clus_k_%d_sse_test_sortedByTrain_datPts%dk',length(hexPts),nPoints/1000)];
    saveas(gcf,fname,'png');
end

% cross-validation - plot SSE over datasets

figure; hold on;
plot(tsseXval2','Color',[.75 .75 .75]);
plot(tsseXval2(min(mean(tsseXval2,2))==mean(tsseXval2,2),:),'Color',[.5 .5 .5],'LineWidth',2); %lowest mean SSE
plot(tsseXval2(1,:),'Color',[.4 .4 .4],'LineWidth',2); %lowest mean SSE from training data
plot(tsseHex,'Color',[0 .447, .8],'LineWidth',3)
plot(tsseSq,'Color',[1,0,0],'LineWidth',3)
% plot([tsseSq;tsseHex]');

title(sprintf('SSE over all Test data sets (k=%d) (nTest=%d)(datPts=%dk)',length(hexPts),nDataSets,nPoints/1000));
if saveplots,
    fname=[wd, sprintf('/kmeans_clus_k_%d_sse_test_overDatasets_datPts%dk',length(hexPts),nPoints/1000)];
%     fname = [fname '_limDatPtsTo85'];
    saveas(gcf,fname,'png');
end


% figure;hist([tsseHex;tsseXval2(min(mean(tsseXval2,2))==mean(tsseXval2,2),:)]',20);






