function xVal_results = xVal_clus(dat,nXvalDataSets, locRange, mu, nIter, dataPts)
% Crossvalidation function - save SSE and SSE per cluster (spreadoutness
% measure) 
% tsse, tsseXval, sseSprdSd, sseSprdVar, sseSprdSdXval, sseSprdVarXval, indSSE, indSpr

% input:  dat, nIter, nXvalDataSets, mu (cluster centre locs)
% output: tsseXval, spreadoutness, indices for top SSE/spreadoutness on training data

%%%%%%%%
%inputs%
%%%%%%%%
% dat:            string - what is the data: 'rand', 'randUnique','cat'
% nXvalDataSets:  number of datasets to generate to do xVal over
% locRange:       range where the data lie (box is [0, 49, but might need
% to edit for different shapes)
% mu:             dims should be nClus x 2 (xy pos) x nIter

% dataPts:        optional (?) - if want to plot the SSE for data used for
% training should be nTrials x 2 (xy pos) - maybe need this in order to
% compare...

%saveDat? - put in saveDir here, or leave for outside?


%if didn't save orig data, just generate some data? should be similar....
if nargin < 5
    switch dat
        case 'rand' % draw random points in a box (uniform distribution)
            dataPts = [randsample(linspace(locRange(1),locRange(2),50),nTrials,'true'); randsample(linspace(locRange(1),locRange(2),50),nTrials,'true')]'; % random points in a box
        case 'randUnique' % draw random points in a box (uniform distribution)
            dataPts = [randsample(linspace(locRange(1),locRange(2),50),nTrials,'true'); randsample(linspace(locRange(1),locRange(2),50),nTrials,'true')]'; % random points in a box
    end
end

nTrials = size(dataPts,1);
nClus   = size(mu,1);

% %testing the function (run a covering map, e.g. iter = 6)
% mu=muAll(:,:,end,:);
% nXvalDataSets = 3;
% dataPts=trials;


%% compute SSE and spreadoutness from the original data

%*****
% atm i don't save orig data or compute sse; i suppose saving orig data will
% be easiest..
%%%%%

tsse       = nan(1,nIter);
sseSprdSd  = nan(1,nIter);
sseSprdVar = nan(1,nIter);
for iterI = 1:nIter
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
    sseSprdSd(iterI) = std(devAvgSSE); % may be better since normalises by nClus (?)
    sseSprdVar(iterI) = var(devAvgSSE);
    
end

%sort by SSE on training dataset; check
[indVal, indSSE] = sort(tsse);
[y, indSSE1] = sort(indSSE);
% indSSE2 = indSSE1; %? no need?

%sort by spreadoutness (SSE per cluster)
[indValSpr, indSpr] = sort(sseSprdSd);
% [indValSpr, indSpr] = sort(sseSprdVar);
[y, indSpr1] = sort(indSpr);

%% Crossvalidation

sseXval        = nan(1,nClus); 
tsseXval       = nan(nIter,nXvalDataSets);
sseSprdSdXval  = nan(nIter,nXvalDataSets);
sseSprdVarXval = nan(nIter,nXvalDataSets);
dataPtsTest    = nan(nTrials,2,nXvalDataSets);

for iDataSets = 1:nXvalDataSets
    fprintf('xVal dataset %d \n',iDataSets);
    switch dat
        case 'rand' % draw random points in a box (uniform distribution)
            dataPtsTest(:,:,iDataSets) = [randsample(linspace(locRange(1),locRange(2),50),nTrials,'true'); randsample(linspace(locRange(1),locRange(2),50),nTrials,'true')]'; % random points in a box
        case 'randUnique' % draw random points in a box (uniform distribution)
            dataPtsTest(:,:,iDataSets) = [randsample(linspace(locRange(1),locRange(2),50),nTrials,'true'); randsample(linspace(locRange(1),locRange(2),50),nTrials,'true')]'; % random points in a box
        case 'cat' % points from clusters of 2D gaussians
            nTrialsCat = floor(nTrials/nCats); % points to sample
            for iCat = 1:nCats
                mu(iCat,:)=randsample(locRange(1)+10:locRange(2)-10,2,'true'); % ±10 so category centres are not on the edge
                datPtsGauss(:,:,iCat) = round(repmat(mu(iCat,:),nTrialsCat,1) + randn(nTrialsCat,2)*R); % key - these are the coordinates of the points
            end
            dataPtsTestTmp = reshape(datPtsGauss,nTrials,2);
            dataPtsTest(:,2,iDataSets) = dataPts(randperm(length(dataPtsTestTmp)),:);
    end
    
    % Crossvalidation start
    for iterI=1:nIter
        %compute distance of existing clusters with new datapoints
        %distTrl= (mu(:,1,iterI)'-dataPts(:,1)).^2+(mu(:,2,iterI)'-dataPts(:,2)).^2; % vectorised - from before
        distXval=(mu(:,1,indSSE1==iterI)-dataPtsTest(:,1,iDataSets)').^2+mu(:,2,indSSE1==iterI)-dataPtsTest(:,2,iDataSets)'.^2;% vectorised - check
        
        [indValsTest, indTest]=min(distXval,[],1); % find which clusters are points closest to
        for iClus = 1:nClus
            sseXval(iClus)=sum(sum([mu(iClus,1,iterI)-dataPtsTest(indTest==iClus,1,iDataSets), mu(iClus,2,iterI)-dataPtsTest(indTest==iClus,2,iDataSets)].^2,2)); %distance from each cluster from training set to datapoints closest to that cluster
%             sseXval(iClus)=sum(sum([(mu(iClus,1,indSSE1==iterI))-dataPtsTest(indTest==iClus,1), (mu(iClus,2,indSSE1==iterI))-dataPtsTest(indTest==iClus,2)].^2,2)); %             % old - ordered by SSE on training set
        end
        tsseXval(iterI,iDataSets)          = sum(sseXval);
        devAvgSSEXval                      = sseXval-mean(sseXval);% compute SSE per cluster (spreadoutness measure)
        sseSprdSdXval(iterI,iDataSets) = std(devAvgSSEXval); % may be better since normalises by nClus (?)
        sseSprdVarXval(iterI,iDataSets) = var(devAvgSSEXval);
    end
end

%variables to save - xVal_results
% Q: do i need indSSE1 and indSpr1? i think not?

xVal_results.tsse           = tsse;
xVal_results.tsseXval       = tsseXval;
xVal_results.sseSprdSd      = sseSprdSd;
xVal_results.sseSprdVar     = sseSprdVar;
xVal_results.sseSprdSdXval  = sseSprdSdXval;
xVal_results.sseSprdVarXval = sseSprdVarXval;
xVal_results.indSSE         = indSSE;
xVal_results.indSpr         = indSpr;
xVal_results.dataPtsTest    = dataPtsTest;

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