function xVal_results = xVal_clusConds(mu, dat,nXvalDataSets, nDataPtsTest, locRange, nIter, dataPts)
% Crossvalidation function - save SSE and SSE per cluster (spreadoutness
% measure) 


%NOTE: mu and dataPts have to be a cell array: mu{nClusConds}(nClus,2,nIter) and dataPts{nClusConds}


%%%%%%%%
%inputs%
%%%%%%%%
% mu:             cluster centre locs - new: should be a cell array mu{nClusConds}(nClus,2,nIter)
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

nClusConds = length(mu);

%%%%%
% NEED TO SEE IF THIS WORKS WITH DATAPTS AS A CELL ARRAY
%%%%
    
%% compute SSE and spreadoutness from the original data (if exist)
tsse       = nan(1,nIter,nClusConds);
sseSprdSd  = nan(1,nIter,nClusConds);
sseSprdVar = nan(1,nIter,nClusConds);
indSSE     = nan(1,nIter,nClusConds);
indSprSd   = nan(1,nIter,nClusConds);
indSprVar  = nan(1,nIter,nClusConds);

if nargin > 6 % if no orig data, don't compute
    for iClusCond = 1:nClusConds
        nClus   = size(mu{iClusCond},1);
        for iterI = 1:nIter
            sse=nan(1,nClus);
            distTrl=(mu{iClusCond}(:,1,iterI)'-dataPts{iClusCond}(:,1)).^2+(mu{iClusCond}(:,2,iterI)'-dataPts{iClusCond}(:,2)).^2; % vectorised - check
            [indVals, indTmp]=min(distTrl,[],2); % find which clusters are points closest to
            for iClus = 1:nClus
                sse(iClus)=sum(sum([mu{iClusCond}(iClus,1,iterI)-dataPts{iClusCond}(indTmp==iClus,1), mu{iClusCond}(iClus,2,iterI)-dataPts{iClusCond}(indTmp==iClus,2)].^2,2)); %distance from each cluster from training set to datapoints closest to that cluster
            end
            tsse(iterI,iClusCond)=sum(sse);
            %SSE per cluster (spreadoutness measure)
            devAvgSSE            = sse-mean(sse);
            sseSprdSd(iterI,iClusCond) = std(devAvgSSE); % may be better since normalises by nClus (?)
            sseSprdVar(iterI,iClusCond) = var(devAvgSSE);
        end
        
        %sort by SSE on training dataset; check
        [indVal, indSSE(:,iClusCond)] = sort(tsse(:,iClusCond));
        % [y, indSSE1] = sort(indSSE);
        
        %sort by spreadoutness (SSE per cluster)
        [indValSpr, indSprSd(:,iClusCond)] = sort(sseSprdSd(:,iClusCond));
        [indValSpr, indSprVar(:,iClusCond)] = sort(sseSprdVar(:,iClusCond));
        % [y, indSpr1] = sort(indSpr);
    end
end

%% Crossvalidation 
tsseXval       = nan(nIter,nXvalDataSets,nClusConds);
sseSprdSdXval  = nan(nIter,nXvalDataSets,nClusConds);
sseSprdVarXval = nan(nIter,nXvalDataSets,nClusConds);
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
                datPtsTmp(iCat,:)=randsample(locRange(1)+10:locRange(2)-10,2,'true'); % ±10 so category centres are not on the edge
                datPtsGauss(:,:,iCat) = round(repmat(datPtsTmp(iCat,:),nTrialsCat,1) + randn(nTrialsCat,2)*R); % key - these are the coordinates of the points
            end
            dataPtsTestTmp = reshape(datPtsGauss,nDataPtsTest,2);
            dataPtsTest(:,2,iDataSets) = dataPts(randperm(length(dataPtsTestTmp)),:);
        case 'circ'
            % Create logical image of a circle
            nSteps = length(linspace(locRange(1),locRange(2),locRange(2)+1));
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
            trialIndTest = randi(length(circPts),nDataPtsTest,1);
            dataPtsTest(:,:,iDataSets)  = circPts(trialIndTest,:);
    end
    
    % Crossvalidation start
    for iClusCond = 1:nClusConds
        nClus     = size(mu{iClusCond},1);
        sseXval   = nan(1,nClus);
        for iterI=1:nIter
            %compute distance of existing clusters with new datapoints
            distXval=(mu{iClusCond}(:,1,iterI)-dataPtsTest(:,1,iDataSets)').^2+mu{iClusCond}(:,2,iterI)-dataPtsTest(:,2,iDataSets)'.^2;% vectorised - not sorted by SSE
            [indValsTest, indTest]=min(distXval,[],1); % find which clusters are points closest to
            for iClus = 1:nClus
                sseXval(iClus)=sum(sum([mu{iClusCond}(iClus,1,iterI)-dataPtsTest(indTest==iClus,1,iDataSets), mu{iClusCond}(iClus,2,iterI)-dataPtsTest(indTest==iClus,2,iDataSets)].^2,2)); %distance from each cluster from training set to datapoints closest to that cluster
            end
            tsseXval(iterI,iDataSets,iClusCond)          = sum(sseXval);
            devAvgSSEXval                                = sseXval-mean(sseXval); % compute SSE per cluster (spreadoutness measure)
            sseSprdSdXval(iterI,iDataSets,iClusCond)     = std(devAvgSSEXval); % may be better since normalises by nClus (?)
            sseSprdVarXval(iterI,iDataSets,iClusCond)    = var(devAvgSSEXval);
        end
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