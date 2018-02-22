function [tsseXval,dataPtsTest] = xVal_clus(dat,nXvalDataSets, mu, dataPts,nIter)
% Crossvalidation function - save SSE and SSE per cluster (spreadoutness
% measure) 

% input args: dat, nIter, nXvalDataSets, mu (cluster centre locations
% output:tsseXval, spreadoutness

%%%%%%%%
%inputs%
%%%%%%%%
% dat:            string - what is the data: 'rand', 'randUnique','cat'
% nXvalDataSets:  number of datasets to generate to do xVal over
% mu:             dims should be nClus x 2 (xy pos) x nIter

% dataPts:        optional (?) - if want to plot the SSE for data used for
% training should be nTrials x 2 (xy pos) - maybe need this in order to
% compare...

%saveDat? - put in saveDir here, or leave for outside?

 
nClus=size(mu,1);
%% compute SSE and spreadoutness from the original data

%*****
% atm i don't save orig data or compute sse; i suppose saving orig data will
% be easiest..
%%%%%

for iterI=1:nIter
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
[y, indSSE1] = sort(indSSE);
% indSSE2 = indSSE1; %? no need?

%sort by spreadoutness (SSE per cluster)
[indValSpr, indSpr] = sort(stdAcrossClus);
% [indValSpr, indSpr] = sort(varAcrossClus);
[y, indSpr1] = sort(indSpr);

%% Crossvalidation

%%%%%%%
%to do:
%%%%%%%
% - save the indices for late - indSSE1, indSpr1
% - edit code below: don't order by SSE on training set: save the indices so can later sort
% by SSE or SSE per cluster


sseXval  = nan(1,nClus); 
tsseXval = nan(nIter,nXvalDataSets);
for iDataSets = 1:nXvalDataSets
    fprintf('xVal dataset %d \n',iDataSets);
    switch dat
        case 'rand' % draw random points in a box (uniform distribution)
            dataPtsTest = [randsample(linspace(locRange(1),locRange(2),50),nTrials,'true'); randsample(linspace(locRange(1),locRange(2),50),nTrials,'true')]'; % random points in a box
        case 'randUnique' % draw random points in a box (uniform distribution)
            dataPtsTest = [randsample(linspace(locRange(1),locRange(2),50),nTrials,'true'); randsample(linspace(locRange(1),locRange(2),50),nTrials,'true')]'; % random points in a box
        case 'cat' % points from clusters of 2D gaussians
            nTrialsCat = floor(nTrials/nCats); % points to sample
            for iCat = 1:nCats
                mu(iCat,:)=randsample(locRange(1)+10:locRange(2)-10,2,'true'); % ±10 so category centres are not on the edge
                datPtsGauss(:,:,iCat) = round(repmat(mu(iCat,:),nTrialsCat,1) + randn(nTrialsCat,2)*R); % key - these are the coordinates of the points
            end
            dataPtsTest = reshape(datPtsGauss,nTrials,2);
            dataPtsTest = dataPts(randperm(length(dataPtsTest)),:);
    end
    
    % Crossvalidation start
    for iterI=1:nIter
        %compute distance of existing clusters with new datapoints
        %             distTrl= (mu(:,1,iterI)'-dataPts(:,1)).^2+(mu(:,2,iterI)'-dataPts(:,2)).^2; % vectorised - from before
        distXval=(mu(:,1,indSSE1==iterI)-dataPtsTest(:,1)').^2+mu(:,2,indSSE1==iterI)-dataPtsTest(:,2)'.^2;% vectorised - check
        
        [indValsTest, indTest]=min(distXval,[],1); % find which clusters are points closest to
        
        for iClus = 1:nClus
            sseXval(iClus)=sum(sum([(mu(iClus,1,indSSE1==iterI))-dataPtsTest(indTest==iClus,1), (mu(iClus,2,indSSE1==iterI))-dataPtsTest(indTest==iClus,2)].^2,2)); %distance from each cluster from training set to datapoints closest to that cluster
        end
        tsseXval(iterI,iDataSets)=sum(sseXval);
        
        %SSE per cluster (spreadoutness measure)
        devAvgSSEXval            = sseXval-mean(sseXval);
        stdAcrossClusXval(iterI,iDataSets) = std(devAvgSSEXval); % may be better since normalises by nClus (?)
        varAcrossClusXval(iterI,iDataSets) = var(devAvgSSEXval);
    end
end


% if saveDat
%     switch dat
%         case 'randUnique'
%             fname = [saveDir, sprintf('/kmeans_nClus_%d-%d_uniquePts_xVal_%ddatasets_%diters',kVals(1),kVals(end),nXvalDataSets,nIter)];
%         case 'rand'
%             fname = [saveDir, sprintf('/kmeans_nClus_%d-%d_randPts_xVal_%ddatasets_%diters',kVals(1),kVals(end),nXvalDataSets,nIter)];
%     end
%     save(fname,'tsseXval','kVals')
% end
end