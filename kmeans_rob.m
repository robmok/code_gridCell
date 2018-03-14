%% k means

clear all;

wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
% wd='/home/robmok/Documents/Grid_cell_model'; %on love01

cd(wd);

codeDir = [wd '/code_gridCell'];
saveDir = [wd '/data_gridCell'];
addpath(codeDir); addpath(saveDir);
addpath(genpath([codeDir '/gridSCORE_packed']));

runScriptTimes=2; %rerun whole script x times
saveDat = 1;

kVals = 3:22; 
 kVals = 23:30; 
% kVals = 10;
nKvals = length(kVals);

dat = 'circ'; %'rand' points in a box, randUnique (all unique points in box), or 'cat'
nKmeans = 200;  % run k means n times 
% nPoints = 5000;  %how many locations/datapoints - not used for 'randUnique'
nPointsVals = [5000, 10000]; % also did 3000 before

%%%%
locRange  = [0, 49];
spacing   = linspace(locRange(1),locRange(2),locRange(2)+1);
nSteps = length(spacing);
gaussSmooth=1; % smoothing for density plot / autocorrelogram 

% normal k means / 'cat learning'
nCats = 2;
sigmaG = [3 0; 0 3]; % isotropic
R = chol(sigmaG); 

%% run k means

nUpdSteps   =  30;    % update steps in the k means algorithm - 40 for random init, 25 for forgy; kmeans++ 20 fine, 22 safe

for iRun=1:runScriptTimes
for iNpts = 1:length(nPointsVals)
    nPoints = nPointsVals(iNpts);
    %start
    densityPlotCentres = zeros(length(spacing),length(spacing),nKmeans,nKvals);
    tssekVals = nan(nKvals,nKmeans);
    indSSE1 = nan(nKvals,nKmeans);
    indSSE2 = nan(nKvals,nKmeans);
    muAllkVals = cell(1,nKvals);
    gA = nan(nKmeans,4,nKvals);
    gW = nan(nKmeans,4,nKvals);
    for iKvals = 1:nKvals
        nK=kVals(iKvals);
        muAll=nan(nK,2,nKmeans);
        tsseAll = nan(1,nKmeans);
        fprintf('nK = %d \n',nK);
        
        %generate/sample data
        switch dat
            case 'randUnique'  %all unique points in square box
                load([saveDir '/randTrialsBox_trialsUnique']);
                dataPts = trialsUnique;
                % does it matter how many points there are if all the same points?
                % dataPts = repmat(dataPts,50,1);
            case 'square' %uniformly sample the box
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
        
        
        tic
        for kMeansIter=1:nKmeans
            if mod(kMeansIter,200)==0
                fprintf('Running k means iteration %d \n',kMeansIter);
            end
            mu   = nan(nK,2,nUpdSteps+1);
            densityPlotClus = zeros(length(spacing),length(spacing),nK);
            for i = 1:nK
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
                %combine/make density plot
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
        fname = [saveDir, sprintf('/kmeans_nK_%d-%d_%s_nPoints%d_%diters',kVals(1),kVals(end),dat,nPoints,nKmeans)];
        cTime=datestr(now,'HHMMSS'); fname = sprintf([fname '_%s'],cTime);
        save(fname,'muAllkVals','tssekVals', 'gA','gW','densityPlotCentres','indSSE1','indSSE2','kVals')
    end
end
end
