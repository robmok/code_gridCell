%% make a box of gaussians

%current problems:
 
clear all;
% close all;

nGaussians = 50;
nTrials    = 7500;  %how many locations in the box / trials
borderClusters = 0; %add border clusters / 'cells' or not
twoLearningRates = 0;

updWinnerOnly = 1; %if 0, update both winner and others

if ~borderClusters
    nBorderClusters = 0;
end 
colgrey = [.5, .5, .5];

%box
sigmaGauss = 1.5;%.25; %std of the gaussian %.25
nSteps = 50; %to define spacing beween each loc in box
locRange = .9;%.9; from -locRange to locRange


[xxx, yyy] = meshgrid(linspace(-1,1,nSteps));
w = length(xxx); % length of the box
nPix = w^2; % total number of pixels
pixCoords = reshape(cat(3,xxx,yyy),[nPix 2]);

mixVec = zeros(nPix,1);
for gaussianI = 1:nGaussians
    % draw a center location
    c(:,:,gaussianI) = -locRange + locRange.*2.*rand(1,2);  %uniform distr from -1 to 1 (see help rand)
        
    % draw a covariance matrix
    S = eye(2)*sigmaGauss; % exactly isotropic, fixed-size Gaussian
    
    % compute the Gaussian
    xmc = pixCoords - repmat(c(:,:,gaussianI),[nPix 1]);
    g = exp(sum((-0.5.*(xmc*inv(S).*xmc))/S.*(2.*pi).^2,2));
%     g= mvnpdf(pixCoords, repmat(c(:,:,gaussianI),[nPix 1]), S);
    
    % add it to the mixture
    mixVec = mixVec + g;
    
    %store each separately
    mixVecSep(:,gaussianI) = g;
end

%add some border clusters
if borderClusters
    clustersPerSide=5;
    sides=4;
    borderLocs=linspace(-1,1,clustersPerSide+2);
    borderLocs=borderLocs(2:end-1)';
    borderClusters  = [-ones(clustersPerSide,1), borderLocs; ones(clustersPerSide,1), borderLocs; borderLocs, -ones(clustersPerSide,1); borderLocs, ones(clustersPerSide,1)];
    nBorderClusters = size(borderClusters,1);
    for gaussianI = 1:nBorderClusters,
        % draw a center location
        c(:,:,nGaussians+gaussianI) = borderClusters(gaussianI,:);
        
        % draw a covariance matrix
        S = eye(2)*sigmaGauss; % exactly isotropic, fixed-size Gaussian
        
        % compute the Gaussian
        xmc = pixCoords - repmat(c(:,:,nGaussians+gaussianI),[nPix 1]);
        g = exp(sum((-0.5.*(xmc*inv(S).*xmc))/S.*(2.*pi).^2,2));
        %     g= mvnpdf(pixCoords, repmat(c(:,:,gaussianI),[nPix 1]), S);
        
        % add it to the mixture
        mixVec = mixVec + g;
        
        %store each separately
        mixVecSep(:,nGaussians+gaussianI) = g;
    end
end

mixVol = reshape(mixVec,[w w]); % reshape vector box

% figure;
% imagesc(mixVol); colormap bone
% title('Gaussians in a box')

%find where the gaussians are
for iCentre = 1:nGaussians+nBorderClusters,
    ind=find(mixVecSep(:,iCentre)>.01);
%     ind=find(mixVecSep(:,iCentre)>1e-10);
    boxVec=zeros(nPix,1);
    boxVec(ind)=1;
    locInds(:,:,iCentre)=reshape(boxVec,[w,w]);
end
% figure;
% imagesc(sum(locInds,3));
% title('Gaussians in a box - indices')

clear xyInds
for iCentre = 1:nGaussians+nBorderClusters,
    xyInds{iCentre}=[]; %SPECIFY how large this will be? - why are they different number of pixels?
    for i = 1:size(locInds,1)
        if any(locInds(:,i,iCentre)),
            ind=find(locInds(:,i,iCentre));
            %find x and y coords for the gauss
            xyInds{iCentre} = [xyInds{iCentre}; repmat(i,length(ind),1), ind]; % [x; y]            
        end
    end
end

%%%%% 
%need to make the xyInds into the space the box is in? no need to change
%xyInds, but have another index that tells us what value the ind is
ind2grid = xxx(1,:);

%if border cells - add them to nGaussians
nGaussians = nGaussians+nBorderClusters;
colors = distinguishable_colors(nGaussians); %function for making distinguishable colors for plotting
%% 2d gauss model

% current state:
% - with random locations on each trial, epsMu=9.8, alpha = .05 - what std gauss?
% - with continuous locations, epsMu=9.8, alpha=.05, std gauss = 0.15

% close all;

% parameters to fit changes in mean
epsMu=15;%.02; 
epsMu1=.01;
epsMu2=.0001;
% alpha=1;%.001; %  target - how much to shift (weighted by diff from cluster mean) if feedback says 'correct category'

sigma = sigmaGauss;
fitDeltaMuX=@(epsMu,x,y,fdbck,act)(epsMu.*(fdbck-act).*exp(-(x.^2+y.^2)./2.*sigma.^2).*(x./(2.*pi.*sigma.^4)));
%alt? % (-x.(exp(-(x.^2+y.^2)./2.*sigma.^2))/(2.*pi.^4)

fitDeltaMuY=@(epsMu,x,y,fdbck,act)(epsMu.*(fdbck-act).*exp(-(x.^2+y.^2)./2.*sigma.^2).*(y./(2.*pi.*sigma.^4)));

stepSize=diff(linspace(-1,1,nSteps)); stepSize=stepSize(1); %smallest diff between locs

% random points
trials=[randsample(linspace(-1,1,101),nTrials,'true'); randsample(linspace(-1,1,101),nTrials,'true')]';

% % points on a path, like a mouse in a box
% %in a box
% % trials = nan(nTrials,2);
% trials(1,:)=randsample(linspace(-1,1,nSteps),2,'true');
% for i=1:nTrials-1
%     
%     moveDir=randsample([-stepSize*2,-stepSize,0,stepSize,stepSize*2],2,'true'); %move in a random direction or stay
%         
%     %if edge, 1 or 101 (check if these are the edges) stay or move
%     if trials(i,1)<=-1 
%         moveDir(1) = randsample([0,stepSize,stepSize*2],1,'true');
%     end
%     if trials(i,2)<=-1,
%         moveDir(2) = randsample([0,stepSize,stepSize*2],1,'true');
%     end
%     if trials(i,1)>=1,
%         moveDir(1) = randsample([0,-stepSize,-stepSize*2],1,'true');
%     end
%     if trials(i,2)>=1,
%         moveDir(2) = randsample([0,-stepSize,-stepSize*2],1,'true');
%     end
%     
%     trials(i+1,:)=trials(i,:)+moveDir;% add 1 or -1 
% end

%initial cluster at each center (later initialise at first trial at that
%cluster's center) - need also to remove previous trials? in a way that
%cluster has not started/been activated yet?

mu=nan(nGaussians,2);
for iClus = 1:nGaussians,
%     mu(iClus,:)=c(:,:,iClus);
    mu(iClus,:)=ind2grid(xyInds{iClus}(randi(length(xyInds{iClus}),1),:)); %random point at their own cluster - converting from xyInds to linspace(-1,1,101) location
end


clusInd = 1:nGaussians;
alpha = zeros(nGaussians,1);
epsMuVec = zeros(nGaussians,1);
for iTrl=1:length(trials)
    for iClus = 1:nGaussians,
        
        act(:,iTrl)=mvnpdf(trials(iTrl,:),mu(:,:,iTrl),S);
        
        %compute feedback v4 - closest cluster target (feedback) = 1, others target 0
%         dist2Clus = sqrt((c(:,:,iClus)-trials(iTrl,:))*(c(:,:,iClus)-trials(iTrl,:))');
%         dist2Clus = sqrt(sum([squeeze(c(:,1,:))-trials(iTrl,1), squeeze(c(:,2,:))-trials(iTrl,2)].^2,2)); % vectorising euclid dist - sqrt(sum((a-b).^2)), since can't use xval method
        dist2Clus = sqrt(sum([mu(:,1,iTrl)-trials(iTrl,1),  mu(:,2,iTrl)-trials(iTrl,2)].^2,2)); % vectorising euclid dist - sqrt(sum((a-b).^2)), since can't use xval method

        closestC=find(min(dist2Clus)==dist2Clus);
        if numel(closestC)>1, %if more than 1, randomly choose one
            closestC = randsample(closestC,1); 
        end
        alpha = zeros(nGaussians,1);
        alpha(closestC) = 1;
        if twoLearningRates,
            epsMuVec(closestC) = epsMu1;
            epsMuVec(clusInd~=closestC) = epsMu2;
        else
            epsMuVec = repmat(epsMu,numel(epsMuVec),1);
        end
        
        if updWinnerOnly
            epsMuVec(clusInd~=closestC)=0;
        end
        
        deltaMu(:,1,iTrl)=fitDeltaMuX(epsMuVec,trials(iTrl,1)-mu(:,1,iTrl),trials(iTrl,2)-mu(:,2,iTrl),alpha,act(:,iTrl));
        deltaMu(:,2,iTrl)=fitDeltaMuY(epsMuVec,trials(iTrl,1)-mu(:,1,iTrl),trials(iTrl,2)-mu(:,2,iTrl),alpha,act(:,iTrl));

        % update mean estimates
        if iTrl~=length(trials) %no need to update for last trial +1)
            mu(:,1,iTrl+1) = mu(:,1,iTrl) + deltaMu(:,1,iTrl);
            mu(:,2,iTrl+1) = mu(:,2,iTrl) + deltaMu(:,2,iTrl);
        end
        
        
    end
    
end
% 
% % plot
% figure; 
% plot(trials(:,1),trials(:,2),'Color',colgrey); hold on;
% for iClus = 1:nGaussians, plot(c(:,1,iClus),c(:,2,iClus),'Color',colors(iClus,:),'LineWidth',10); end
% xlim([-1.1,1.1]); ylim([-1.1,1.1]);
% 
% figure;
% for iClus = 1:nGaussians
% %     plot(squeeze(mu(iClus,1,:)),squeeze(mu(iClus,2,:)),'Color',colors(iClus,:)); hold on;
%     plot(squeeze(mu(iClus,1,round(nTrials/4):end)),squeeze(mu(iClus,2,round(nTrials/4):end)),'Color',colors(iClus,:)); hold on;
% %     plot(squeeze(mu(iClus,1,round(end-30:end))),squeeze(mu(iClus,2,end-30:end)),'Color',colors(iClus,:)); hold on;
% 
%     plot(c(:,1,iClus),c(:,2,iClus),'.','Color',colors(iClus,:),'MarkerSize',40); plot(c(:,1,iClus),c(:,2,iClus),'kx','LineWidth',2,'MarkerSize',12); 
% end
% xlim([-1.1,1.1]); ylim([-1.1,1.1]);
% 
%plot without the starting centres
figure;
for iClus = 1:nGaussians

    plot(squeeze(mu(iClus,1,round(nTrials/4):end)),squeeze(mu(iClus,2,round(nTrials/4):end)),'Color',colors(iClus,:)); hold on;

%         plot(squeeze(mu(iClus,1,:)),squeeze(mu(iClus,2,:)),'Color',colors(iClus,:)); hold on;
%     plot(squeeze(mu(iClus,1,round(nTrials/4):end)),squeeze(mu(iClus,2,round(nTrials/4):end)),'Color',colors(iClus,:)); hold on;
%         plot(c(:,1,iClus),c(:,2,iClus),'.','Color',colors(iClus,:),'MarkerSize',40); plot(c(:,1,iClus),c(:,2,iClus),'kx','LineWidth',2,'MarkerSize',12); 
        plot(squeeze(mean(mu(iClus,1,end))),squeeze(mean(mu(iClus,2,end))),'.','Color',colors(iClus,:),'MarkerSize',40); hold on;
%         plot(squeeze(mean(mu(iClus,1,round(nTrials/8):end))),squeeze(mean(mu(iClus,2,round(nTrials/8):end))),'.','Color',colors(iClus,:),'MarkerSize',40); hold on;
   
end
xlim([-1.1,1.1]); ylim([-1.1,1.1]);


%plot over time
figure;
for iTrl = 1:nTrials
%     if mod(iTrl,50)==0, fprintf('Trial %d \n',iTrl); end
%     % plot all points
%     plot(squeeze(mu(1,1,iTrl)),squeeze(mu(1,2,iTrl)),'.','Color',colors(1,:),'MarkerSize',10); hold on;
%     plot(squeeze(mu(2,1,iTrl)),squeeze(mu(2,2,iTrl)),'.','Color',colors(2,:),'MarkerSize',10); hold on;
%     plot(squeeze(mu(3,1,iTrl)),squeeze(mu(3,2,iTrl)),'.','Color',colors(3,:),'MarkerSize',10); hold on;
%     plot(squeeze(mu(4,1,iTrl)),squeeze(mu(4,2,iTrl)),'.','Color',colors(4,:),'MarkerSize',10); hold on;
%     plot(squeeze(mu(5,1,iTrl)),squeeze(mu(5,2,iTrl)),'.','Color',colors(5,:),'MarkerSize',10); hold on;
%     drawnow;
    if mod(iTrl,50)==0, %plot centers after x trials
        for i=1:nGaussians %10
            plot(squeeze(mu(i,1,iTrl)),squeeze(mu(i,2,iTrl)),'.','Color',colors(i,:),'MarkerSize',20); hold on;
        end
        drawnow;
    end
end
xlim([-1.1,1.1]); ylim([-1.1,1.1]);

%%
% %euclid dist from each loc to next
% tmpDist=diff(trials);
% euclidXy=[0; sqrt(sum(tmpDist.^2,2))]; %first trial no dist

%euclid dist to mu -double check if this is what i want? i think so
muTmp=permute(mu,[3,2,1]);

for iClus = 1:3,%nGaussians
    tmpDist=squeeze(muTmp(:,:,iClus))-trials;
    euclidXy(:,iClus)=sqrt(sum(tmpDist.^2,2));
    
    
    figure; hold on;
    subplot(4,2,1); plot(trials(:,1),'Color',colgrey); title('x Location');
    subplot(4,2,2); plot(trials(:,2),'Color',colgrey); title('y Location');
    
    subplot(4,2,3); plot(euclidXy(:,iClus),'Color',colgrey); title('Euclidean dist from previous mu');
    
    subplot(4,2,4); plot(act(iClus,:),'Color',colors(iClus,:)); title(sprintf('Activation - Cluster %d',iClus));
    subplot(4,2,5); plot(squeeze(mu(iClus,1,:)),'Color',colors(iClus,:)); title(sprintf('mu: x-axis - Cluster %d',iClus));
    subplot(4,2,6); plot(squeeze(mu(iClus,2,:)),'Color',colors(iClus,:)); title(sprintf('mu: y-axis - Cluster %d',iClus));
    subplot(4,2,7); plot(squeeze(deltaMu(iClus,1,:)),'Color',colors(iClus,:)); title(sprintf('deltaMu: x-axis- Cluster %d',iClus));
    subplot(4,2,8); plot(squeeze(deltaMu(iClus,2,:)),'Color',colors(iClus,:)); title(sprintf('deltaMu: y-axis- Cluster %d',iClus));
end

% % %together
% figure; hold on;
% subplot(4,2,1); plot(trials(:,1),'k'); title('x Location');
% subplot(4,2,2); plot(trials(:,2),'k'); title('y Location');
% subplot(4,2,3); hold on; plot(act(:,1)); plot(act(:,2),'r'); plot(act(:,3),'g'); title('Activation'); 
% subplot(4,2,5); hold on; plot(squeeze(mu(1,1,:))); plot(squeeze(mu(2,1,:)),'r'); plot(squeeze(mu(3,1,:)),'g'); title('mu - x-axis');
% subplot(4,2,6); hold on; plot(squeeze(mu(1,2,:))); plot(squeeze(mu(2,2,:)),'r'); plot(squeeze(mu(3,2,:)),'g'); title('mu - y-axis');
% subplot(4,2,7); hold on; plot(squeeze(deltaMu(:,1,1))); plot(squeeze(deltaMu(:,2,1)),'r'); plot(squeeze(deltaMu(:,3,1)),'g'); title('deltaMu - x-axis');
% subplot(4,2,8); hold on; plot(squeeze(deltaMu(:,1,2))); plot(squeeze(deltaMu(:,2,2)),'r'); plot(squeeze(deltaMu(:,3,2)),'g'); title('deltaMu - y-axis');

%%
for iClus = 1:nGaussians
figure;  hold on;
subplot(1,3,1);
plot(squeeze(mu(iClus,1,:)),squeeze(mu(iClus,2,:))); hold on;
plot(c(:,1,iClus),c(:,2,iClus),'ro'); 
xlim([-1.75,1.75]); ylim([-1.75,1.75]);
subplot(1,3,2);
plot(squeeze(mu(iClus,1,:))); ylim([-1,1]);
subplot(1,3,3);
plot(squeeze(mu(iClus,2,:))); ylim([-1,1]);

end

%% other plots

%activation over time; high initially then goes down (away from own
%cluster?)
figure; hold on;
subplot(1,3,1); plot(act(:,1))
subplot(1,3,2); plot(act(:,2));
subplot(1,3,3); plot(act(:,3));

%deltaMu
figure; hold on;
subplot(1,3,1); plot(squeeze(deltaMu(:,1,:)));
subplot(1,3,2); plot(squeeze(deltaMu(:,2,:)))
subplot(1,3,3); plot(squeeze(deltaMu(:,3,:)));

figure; plot(act(:,1),squeeze(deltaMu(:,1,1)),'b.');
figure; plot(act(:,1),squeeze(deltaMu(:,1,2)),'g.');


figure; plot(trials(:,1),squeeze(deltaMu(:,1,1)),'b.'); hold on;
plot(trials(:,2),squeeze(deltaMu(:,1,2)),'.g')


figure;  plot(trials(:,1)-squeeze(mu(1,1,:)),squeeze(deltaMu(:,1,1)),'b.');
hold on; plot(trials(:,2)-squeeze(mu(1,2,:)),squeeze(deltaMu(:,1,2)),'g.');


