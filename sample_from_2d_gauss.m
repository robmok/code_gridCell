nTrials=40000;

spacing=linspace(locRange(1),locRange(2),locRange(2)+1);
densityPlot = zeros(length(spacing),length(spacing));

nCats   = 2; %2 categories
nPoints = floor(nTrials/nCats); % points to sample
for iCat = 1:nCats
    mu(iCat,:)=randsample(locRange(1)+10:locRange(2)-10,2,'true'); % ±15 so category centres are not on the edge
end

sigmaG = [3 0; 0 3]; R = chol(sigmaG);    % isotropic
% sigmaG = [1 .5; .5 2]; R = chol(sigmaG);  % non-isotropic

for iCat = 1:nCats
    clus = round(repmat(mu(iCat,:),nPoints,1) + randn(nPoints,2)*R);
    for iPts=1:nPoints
        densityPlot(clus(iPts,1),clus(iPts,2)) = densityPlot(clus(iPts,1),clus(iPts,2))+1;
    end
end

imagesc(densityPlot)

