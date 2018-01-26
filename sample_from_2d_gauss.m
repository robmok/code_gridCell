nTrials=40000;

spacing=linspace(locRange(1),locRange(2),locRange(2)+1);
densityPlot = zeros(length(spacing),length(spacing));

nCats   = 2; %2 categories
nPoints = floor(nTrials/nCats); % points to sample
sigmaG = [3 0; 0 3]; R = chol(sigmaG);    % isotropic
% sigmaG = [1 .5; .5 2]; R = chol(sigmaG);  % non-isotropic


for iCat = 1:nCats
%     mu(iCat,:)=randsample(locRange(1)+10:locRange(2)-10,2,'true'); % ±15 so category centres are not on the edge
    datPtsGauss(:,:,iCat) = round(repmat(mu(iCat,:),nPoints,1) + randn(nPoints,2)*R); % key - these are the coordinates of the points
    for iPts=1:nPoints
        densityPlot(datPtsGauss(iPts,1,iCat),datPtsGauss(iPts,2,iCat)) = densityPlot(datPtsGauss(iPts,1,iCat),datPtsGauss(iPts,2,iCat))+1;
    end
end

figure;
imagesc(densityPlot); %% visualise the points on a density plot

trials = reshape(datPtsGauss,nTrials,2);
trials=trials(randperm(length(trials)),:);







%% for main script, no need plot density map - just copy this over

% nTrials=40000;

% draw points from 2 categories (gaussian) from a 2D feature space
nCats   = 2; %2 categories
sigmaG = [3 0; 0 3]; R = chol(sigmaG);    % isotropic
% sigmaG = [1 .5; .5 2]; R = chol(sigmaG);  % non-isotropic

nPoints = floor(nTrials/nCats); % points to sample
for iCat = 1:nCats
    mu(iCat,:)=randsample(locRange(1)+10:locRange(2)-10,2,'true'); % ±10 so category centres are not on the edge
    datPtsGauss(:,:,iCat) = round(repmat(mu(iCat,:),nPoints,1) + randn(nPoints,2)*R); % key - these are the coordinates of the points
end

trials = reshape(datPtsGauss,nTrials,2);
trials = trials(randperm(length(trials)),:);