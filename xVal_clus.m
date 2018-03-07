function xVal_results = xVal_clus(mu, dat,nXvalDataSets, nDataPtsTest, locRange, nIter, dataPts)
% Crossvalidation function - save SSE and SSE per cluster (spreadoutness
% measure) 

%%%%%%%%
%inputs%
%%%%%%%%
% mu:             cluster centre locs - dims should be nClus x 2 (xy pos) x nIter
% dat:            string - what is the data: 'rand', 'randUnique','cat'
% nXvalDataSets:  number of datasets to generate to do xVal over
% nDataPtsTest:   number of points to generate in test set (same as
% nTrials?)
% locRange:       range where the data lie (box is [0, 49, but might need
% to edit for different shapes)

% dataPts:        OPTIONAL - if want to plot the SSE for data used for
% training should be nTrials x 2 (xy pos) - maybe need this in order to
% compare...

%%%%%%%%
%output%
%%%%%%%%
% tsseXval, spreadoutness, indices for top SSE/spreadoutness on training data
% ++

%saveDat? - put in saveDir here, or leave for outside?

nClus   = size(mu,1);

%% compute SSE and spreadoutness from the original data (if exist)
tsse       = nan(1,nIter);
sseSprdSd  = nan(1,nIter);
sseSprdVar = nan(1,nIter);
indSSE     = nan(1,nIter);
indSprSd   = nan(1,nIter);
indSprVar  = nan(1,nIter);

if nargin > 6 % if no orig data, don't compute
for iterI = 1:nIter
    sse=nan(1,nClus);
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
% [y, indSSE1] = sort(indSSE);

%sort by spreadoutness (SSE per cluster)
[indValSpr, indSprSd] = sort(sseSprdSd);
[indValSpr, indSprVar] = sort(sseSprdVar);
% [y, indSpr1] = sort(indSpr);
end

%% Crossvalidation

sseXval        = nan(1,nClus); 
tsseXval       = nan(nIter,nXvalDataSets);
sseSprdSdXval  = nan(nIter,nXvalDataSets);
sseSprdVarXval = nan(nIter,nXvalDataSets);
dataPtsTest    = nan(nDataPtsTest,2,nXvalDataSets);

for iDataSets = 1:nXvalDataSets
    fprintf('xVal dataset %d \n',iDataSets);
    switch dat
        case 'square' % draw random points in a box (uniform distribution)
            dataPtsTest(:,:,iDataSets) = [randsample(linspace(locRange(1),locRange(2),50),nDataPtsTest,'true'); randsample(linspace(locRange(1),locRange(2),50),nDataPtsTest,'true')]'; % random points in a box
        case 'randUnique' % draw random points in a box (uniform distribution)
            dataPtsTest(:,:,iDataSets) = [randsample(linspace(locRange(1),locRange(2),50),nDataPtsTest,'true'); randsample(linspace(locRange(1),locRange(2),50),nDataPtsTest,'true')]'; % random points in a box
        case 'cat' % points from clusters of 2D gaussians
            nTrialsCat = floor(nDataPtsTest/nCats); % points to sample
            for iCat = 1:nCats
                mu(iCat,:)=randsample(locRange(1)+10:locRange(2)-10,2,'true'); % ±10 so category centres are not on the edge
                datPtsGauss(:,:,iCat) = round(repmat(mu(iCat,:),nTrialsCat,1) + randn(nTrialsCat,2)*R); % key - these are the coordinates of the points
            end
            dataPtsTestTmp = reshape(datPtsGauss,nDataPtsTest,2);
            dataPtsTest(:,2,iDataSets) = dataPts(randperm(length(dataPtsTestTmp)),:);
        case 'circ'
            % Create logical image of a circle
            imageSizeX = nSteps;
            [columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeX);
            centerX = nSteps/2; centerY = nSteps/2;
            radius = nSteps/2-1;
            circIm = (rowsInImage - centerY).^2 ...
                + (columnsInImage - centerX).^2 <= radius.^2;
            circPts=[]; % find circle points in XY coords
            for iX=1:length(circIm)
                yVals = find(circIm(iX,:));
                circPts = [circPts; ones(length(yVals),1)*iX, yVals'];
            end
            trialInd=randi(length(circPts),nTrials,1);
            dataPtsTest=circPts(trialInd,:);
            trialIndTest = randi(length(circPts),nTrials,1);
            dataPtsTest  = circPts(trialIndTest,:);
    end
    
    % Crossvalidation start
    for iterI=1:nIter
        %compute distance of existing clusters with new datapoints
        distXval=(mu(:,1,iterI)-dataPtsTest(:,1,iDataSets)').^2+mu(:,2,iterI)-dataPtsTest(:,2,iDataSets)'.^2;% vectorised - not sorted by SSE
        [indValsTest, indTest]=min(distXval,[],1); % find which clusters are points closest to
        for iClus = 1:nClus
            sseXval(iClus)=sum(sum([mu(iClus,1,iterI)-dataPtsTest(indTest==iClus,1,iDataSets), mu(iClus,2,iterI)-dataPtsTest(indTest==iClus,2,iDataSets)].^2,2)); %distance from each cluster from training set to datapoints closest to that cluster
        end
        tsseXval(iterI,iDataSets)          = sum(sseXval);
        devAvgSSEXval                      = sseXval-mean(sseXval); % compute SSE per cluster (spreadoutness measure)
        sseSprdSdXval(iterI,iDataSets)     = std(devAvgSSEXval); % may be better since normalises by nClus (?)
        sseSprdVarXval(iterI,iDataSets)    = var(devAvgSSEXval);
    end
end

%variables to save - xVal_results
xVal_results.tsse           = tsse;
xVal_results.tsseXval       = tsseXval;
xVal_results.sseSprdSd      = sseSprdSd;
xVal_results.sseSprdVar     = sseSprdVar;
xVal_results.sseSprdSdXval  = sseSprdSdXval;
xVal_results.sseSprdVarXval = sseSprdVarXval;
xVal_results.indSSE         = indSSE;
xVal_results.indSprSd       = indSprSd;
xVal_results.indSprVar      = indSprVar;
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