function [tsseXVal,dataPtsTest] = xVal_clus(dat,nXvalDataSets, mu, dataPts)

% input args: dat, nIter, nXvalDataSets, mu (cluster centre locations

% output:tsseXval, spreadoutness

% put in saveDir here, or leave for outside?

%%%%%%%%
%inputs%
%%%%%%%%
% dat:            string - what is the data: 'rand', 'randUnique','cat'
% nXvalDataSets:  number of datasets to generate to do xVal over
% mu:             dims should be nClus x 2 (xy pos) x nIter

% dataPts:        optional (?) - if want to plot the SSE for data used for training should be nTrials x 2 (xy pos) 


nClus=size(mu,1);
nIter=size(mu,3);

% Crossvalidation function - save SSE and SSE per cluster (spreadoutness
% measure)





%compute SSE and spreadoutness from the original data

for iXval=1:nIter
    sse=nan(1,nClus);
%     for iClus = 1:nClus
%         distTrl(:,iClus)=sum([mu(iClus,1,iterI)-dataPts(:,1), mu(iClus,2,iterI)-dataPts(:,2)].^2,2);
%     end
    distTrl=(mu(:,1,iterI)'-dataPts(:,1)).^2+(mu(:,2,iterI)'-dataPts(:,2)).^2; % vectorised - check
    [indVals, indTmp]=min(distTrl,[],2); % find which clusters are points closest to
    for iClus = 1:nClus
        sse(iClus)=sum(sum([mu(iClus,1,iterI)-dataPts(indTmp==iClus,1), mu(iClus,2,iterI)-dataPts(indTmp==iClus,2)].^2,2)); %distance from each cluster from training set to datapoints closest to that cluster
    end
    tsse(iterI)=sum(sse);
    
    
    %SSE per cluster (spreadoutness measure)
    devAvgSSE            = sse-mean(sse);
    stdAcrossClus(iterI) = std(devAvgSSE); % may be better since normalises by nClus (?)
    varAcrossClus(iterI) = var(devAvgSSE);
    
end

%sort by SSE on training dataset; check
[indVal, indSSE] = sort(tsse);
[y, indSSE1(iterI,:)] = sort(indSSE);
indSSE2(iterI,:) = indSSE1(iterI,:);

       

%%%%%
%CHECK - iKvals, iToXval.... iterI; double check I've substituted these
%correctly!
%%%%%

%Crossvalidation

sseXval  = nan(1,nClus); 
tsseXval = nan(nIter,nXvalDataSets);
% muBest = nan(nClus,2,length(bestWorst3));
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
    
        % Crossvalidation start
        for iXval=1:nIter
            %compute distance of existing clusters with new datapoints
%             for iClus = 1:nClus
%                 distXval(:,iClus)=sum([mu(iClus,1,indSSE2(iKvals,:)==iToXval(clusMeans))-dataPtsTest(:,1), mu(iClus,2,indSSE2(iKvals,:)==iToXval(clusMeans))-dataPtsTest(:,2)].^2,2);
%             end
            distXval=(mu(:,1,indSSE2(iterI,:)==iXval)-dataPtsTest(:,1)).^2+mu(:,2,indSSE2(iterI,:)==iXval)-dataPtsTest(:,2).^2;% vectorised - check

%             distTrl= (mu(:,1,iterI)'-dataPts(:,1)).^2                     +(mu(:,2,iterI)'-dataPts(:,2)).^2; % vectorised - from above



            
            [indValsTest, indTest]=min(distXval,[],2); % find which clusters are points closest to
            
            for iClus = 1:nClus
                sseXval(iClus)=sum(sum([(mu(iClus,1,indSSE2(iKvals,:)==iXval))-dataPtsTest(indTest==iClus,1), (mu(iClus,2,indSSE2(iKvals,:)==iXval))-dataPtsTest(indTest==iClus,2)].^2,2)); %distance from each cluster from training set to datapoints closest to that cluster
            end
            tsseXval(iXval,iDataSets,iKvals)=sum(sseXval);
        end
    
    
end

if saveDat
    switch dat
        case 'randUnique'
            fname = [saveDir, sprintf('/kmeans_nClus_%d-%d_uniquePts_xVal_%ddatasets_%diters',kVals(1),kVals(end),nXvalDataSets,nIter)];
        case 'rand'
            fname = [saveDir, sprintf('/kmeans_nClus_%d-%d_randPts_xVal_%ddatasets_%diters',kVals(1),kVals(end),nXvalDataSets,nIter)];
    end
    save(fname,'tsseXval','kVals')
end
end