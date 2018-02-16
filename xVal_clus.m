function [tsseXVal] = xVal_clus(dat, nIter, nXvalDataSets, mu)


% Crossvalidation function - save SSE and SSE per cluster (spreadoutness
% measure)





% input args: dat, nIter, nXvalDataSets, mu (cluster centre locations

% output:tsseXval, spreadoutness

% put in saveDir here, or leave for outside?




%from kmeans

sseXval  = nan(1,nK); 
tsseXval = nan(nKmeans,nXvalDataSets,nKvals);
% muBest = nan(nK,2,length(bestWorst3),nKvals);
for iDataSets = 1:nXvalDataSets
    fprintf('xVal dataset %d \n',iDataSets);
    
    
    
    
    switch dat
        case 'rand' % draw random points in a box (uniform distribution)
            dataPtsTest = [randsample(linspace(locRange(1),locRange(2),50),nPoints,'true'); randsample(linspace(locRange(1),locRange(2),50),nPoints,'true')]'; % random points in a box
        case 'randUnique' % draw random points in a box (uniform distribution)
            dataPtsTest = [randsample(linspace(locRange(1),locRange(2),50),nPoints,'true'); randsample(linspace(locRange(1),locRange(2),50),nPoints,'true')]'; % random points in a box
        case 'cat' % points from clusters of 2D gaussians
            nPointsCat = floor(nPoints/nCats); % points to sample
            for iCat = 1:nCats
                mu(iCat,:)=randsample(locRange(1)+10:locRange(2)-10,2,'true'); % ±10 so category centres are not on the edge
                datPtsGauss(:,:,iCat) = round(repmat(mu(iCat,:),nPointsCat,1) + randn(nPointsCat,2)*R); % key - these are the coordinates of the points
            end
            dataPtsTest = reshape(datPtsGauss,nPoints,2);
            dataPtsTest = dataPts(randperm(length(dataPtsTest)),:);
    end
    
    
    
    
    
    for iKvals = 1:nKvals
        muCurr = muAllkVals{iKvals};
        nK=kVals(iKvals);
        % Crossvalidation start
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
    end
    
end

if saveDat
    switch dat
        case 'randUnique'
            fname = [saveDir, sprintf('/kmeans_nK_%d-%d_uniquePts_xVal_%ddatasets_%diters',kVals(1),kVals(end),nXvalDataSets,nKmeans)];
        case 'rand'
            fname = [saveDir, sprintf('/kmeans_nK_%d-%d_randPts_xVal_%ddatasets_%diters',kVals(1),kVals(end),nXvalDataSets,nKmeans)];
    end
    save(fname,'tsseXval','kVals')
end
end