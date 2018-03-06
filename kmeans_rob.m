%% k means

clear all;

% wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
% wd='/home/robmok/Documents/Grid_cell_model'; %on love01

cd(wd);

codeDir = [wd '/code_gridCell'];
saveDir = [wd '/data_gridCell'];
addpath(codeDir); addpath(saveDir);
addpath(genpath([codeDir '/gridSCORE_packed']));

kVals = 3:11; 
kVals = 12:17; 
% kVals = 18:25; 
% kVals = 26:30; 
% kVals = 15;
nKvals = length(kVals);

saveDat = 1;

doXval        = 1;  %do (and save) crossvalidation
nXvalDataSets = 20; %if do xVal, specify how many datasets to generate

dat = 'circ'; %'rand' points in a box, randUnique (all unique points in box), or 'cat'
nKmeans = 1000;  % run k means n times %1000
nPoints = 3000;  %how many locations/datapoints - not used for 'randUnique', though used for xVal


%%%%
locRange  = [0, 49];%.9; from -locRange to locRange
spacing   = linspace(locRange(1),locRange(2),locRange(2)+1);
nSteps = length(spacing);
gaussSmooth=1; % smoothing for density plot / autocorrelogram 

% % normal k means / 'cat learning'
nCats = 2;
sigmaG = [3 0; 0 3]; R = chol(sigmaG);    % isotropic

switch dat
    case 'randUnique'
        %all unique points in box
        load([saveDir '/randTrialsBox_trialsUnique']);
        dataPts = trialsUnique;
        % does it matter how many points there are if all the same points? e.g.
        % same if just have each trialsUnique twice/x10? - i think not
        % dataPts = repmat(dataPts,50,1);
    case 'square'
        %uniformly sample the box
        dataPts = [randsample(linspace(locRange(1),locRange(2),50),nPoints,'true'); randsample(linspace(locRange(1),locRange(2),50),nPoints,'true')]';
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
        trialInd=randi(length(circPts),nPoints,1);
        dataPts=circPts(trialInd,:);
        trialIndTest = randi(length(circPts),nPoints,1);
        dataPts  = circPts(trialIndTest,:);
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

%% run k means

nUpdSteps   =  30;    % update steps in the k means algorithm - 40 for random init, 25 for forgy; kmeans++ 20 fine, 22 safe

densityPlotCentres = zeros(length(spacing),length(spacing),nKmeans,nKvals);

% check
tssekVals = nan(nKvals,nKmeans);
indSSE1 = nan(nKvals,nKmeans);
indSSE2 = nan(nKvals,nKmeans);
muAllkVals = cell(1,nKvals); % check - this should work
gA = nan(nKmeans,4,nKvals);
gW = nan(nKmeans,4,nKvals);

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
%         switch dat
%             case 'square'
%                 mu(:,:,1) = kmplusInit(dataPts,nK); %k means ++ initiatialization
%             case 'randUnique'
%                 mu(:,:,1) = kmplusInit(dataPts,nK); %k means ++ initiatialization
%             case 'cat'
%                 mu(i,:,1) = dataPts(randi(length(dataPts)),:);   %intiate each cluster with one data point - Forgy method
% %                 mu(i,:,1) = round(locRange(1) + locRange(2).*rand(1,2));  %initiate clusters with a random point in the box
%         end
        mu(:,:,1) = kmplusInit(dataPts,nK); %k means ++ initiatialization
    end
    
    %run kmeans
    [muEnd,tsse] = kmeans_rm(mu,dataPts,nK,nUpdSteps);
    muAll(:,:,kMeansIter) = muEnd;
    tsseAll(kMeansIter)   = tsse;

    %compute gridness
    if ~strcmp(dat,'cat') %if cat learning, no need to compute gridness
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
    [y, indSSE1(iKvals,:)] = sort(indSSE);
    indSSE2(iKvals,:) = indSSE1(iKvals,:);
end

if saveDat
%     switch dat
%         case 'randUnique'
%             fname = [saveDir, sprintf('/kmeans_nK_%d-%d_uniquePts_%diters',kVals(1),kVals(end),nKmeans)];
%         case 'rand'
%             fname = [saveDir, sprintf('/kmeans_nK_%d-%d_randPts_%diters',kVals(1),kVals(end),nKmeans)];
%     end
    fname = [saveDir, sprintf('/kmeans_nK_%d-%d_%s_nPoints%d_%diters',kVals(1),kVals(end),dat,nPoints,nKmeans)];
    save(fname,'muAllkVals','tssekVals', 'gA','gW','densityPlotCentres','indSSE1','indSSE2','kVals')
end

%% Crossvalidation on clusters from k means, assess and save SSE
if doXval
    
sseXval  = nan(1,nK); 
tsseXval = nan(nKmeans,nXvalDataSets,nKvals);
% muBest = nan(nK,2,length(bestWorst3),nKvals);
for iDataSets = 1:nXvalDataSets
    fprintf('xVal dataset %d \n',iDataSets);
    switch dat
        case 'square' % draw random points in a box (uniform distribution)
            dataPtsTest = [randsample(linspace(locRange(1),locRange(2),50),nPoints,'true'); randsample(linspace(locRange(1),locRange(2),50),nPoints,'true')]'; % random points in a box
        case 'randUnique' % draw random points in a box (uniform distribution)
            dataPtsTest = [randsample(linspace(locRange(1),locRange(2),50),nPoints,'true'); randsample(linspace(locRange(1),locRange(2),50),nPoints,'true')]'; % random points in a box
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
            trialInd=randi(length(circPts),nPoints,1);
            trialIndTest     = randi(length(circPts),nPoints,1);
            dataPtsTest      = circPts(trialIndTest,:);
            trialIndTest = randi(length(circPts),nPoints,1);
            dataPtsTest  = circPts(trialIndTest,:);
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
%         %save top 3 bottom 3
%         bestWorst3=[1,2,3,nKmeans-2,nKmeans-1,nKmeans];
%         for iterI=1:length(bestWorst3)
%             muBest(:,:,iterI,iKvals) = muCurr(:,:,indSSE2(iKvals,:)==bestWorst3(iterI));
%         end
    end
    
end

if saveDat
%     switch dat
%         case 'randUnique'
%             fname = [saveDir, sprintf('/kmeans_nK_%d-%d_uniquePts_xVal_%ddatasets_%diters',kVals(1),kVals(end),nXvalDataSets,nKmeans)];
%         case 'rand'
%             fname = [saveDir, sprintf('/kmeans_nK_%d-%d_randPts_xVal_%ddatasets_%diters',kVals(1),kVals(end),nXvalDataSets,nKmeans)];
%     end
    fname = [saveDir, sprintf('/kmeans_nK_%d-%d_%s_nPoints%d_xVal_%ddatasets_%diters',kVals(1),kVals(end),dat,nPoints,nXvalDataSets,nKmeans)];

    save(fname,'tsseXval','kVals')
end
end