function [muEnd, muAll, tsseTrls,sseTrl,cParams] = covering_map_sim(nClus,locRange,box,warpType,epsMuOrig,nTrials,nIter,warpBox,alpha,trials,stochasticType,c)

% if dont save all muAll and muEnd, function output is muEndBest and muAll
% Best, and uncomment the bit at the end of the script

spacing=linspace(locRange(1),locRange(2),locRange(2)+1); 
stepSize=diff(spacing(1:2));

%for crossvalidation - 
nTrialsTest = nTrials;
% nTrialsTest = 5000; %keep it constant if want to compare with nTrials above

muAll    = nan(nClus,2,nTrials,nIter);
tsseTrls = nan(nTrials,nIter);
muEnd    = nan(nClus,2,nIter);
sseTrl = zeros(nClus,nTrials,nIter);


for iterI = 1:nIter
    
    fprintf('iter %d \n',iterI);
    epsMu = epsMuOrig; %revert learning rate back to original if reset
    
    %%%
    %atm putting datapoints - 'trials' - outside to keep it the same 
    %%%
    switch box
        case 'square'
%             trials      = [randsample(spacing,nTrials,'true'); randsample(spacing,nTrials,'true')]';
            dataPtsTest = [randsample(linspace(locRange(1),locRange(2),locRange(2)+1),nTrialsTest,'true'); randsample(linspace(locRange(1),locRange(2),locRange(2)+1),nTrialsTest,'true')]'; % random points in a box
        case 'rect'
%             trials      = [randsample(-1:diff(spacing(1:2)):2,nTrials,'true'); randsample(spacing,nTrials,'true')]';
            dataPtsTest = [randsample(locRange(1):diff(spacing(1:2)):locRange(2)*2,nTrialsTest,'true'); randsample(spacing,nTrialsTest,'true')]';
        case 'trapz'
            trapY=locRange(2).*trapmf(spacing,[spacing(1), spacing(round(length(spacing)*.25)), spacing(round(length(spacing)*.75)),spacing(end)]);
            trapX=spacing;
            trapPts=[];
            for i=1:length(trapY)
               trapPts = [trapPts, [repmat(trapX(i),1,length(0:stepSize:trapY(i))); 0:stepSize:trapY(i)]];
            end
%             trapPts(2,:)=trapPts(2,:).*2-1; %put it back into -1 to 1            
            % use this to select from the PAIR in trapPts
            trialInd     = randi(length(trapPts),nTrials,1);
            trials       = trapPts(:,trialInd)';
            trialIndTest = randi(length(trapPts),nTrials,1);
            dataPtsTest  = trapPts(:,trialIndTest)';
            
        case 'trapzSq' %probably need a more narrow trapezium!
            trapY=locRange(2)/2.*trapmf(spacing,[spacing(1), spacing(round(length(spacing)*.25)), spacing(round(length(spacing)*.75)),spacing(end)]);
            trapY=trapY+floor(length(trapY)./2);
            trapX=spacing;
            trapPts=[];
            for i=1:length(trapY)
               trapPts = [trapPts, [repmat(trapX(i),1,length(0:stepSize:trapY(i))); 0:stepSize:trapY(i)]];
            end
            %make square box attached to it
            sqX=spacing;
            sqY=spacing(1:floor(length(spacing)/2));
            for i=1:length(sqX)
                tic
                for j=1:length(sqY)
                    trapPts = [trapPts, [sqX(i); sqY(j)]];
                end
                toc
            end
            % use this to select from the PAIR in trapPts
            trialInd=randi(length(trapPts),nTrials,1);
            trials=trapPts(:,trialInd)';
            trialIndTest = randi(length(trapPts),nTrials,1);
            dataPtsTest  = trapPts(:,trialIndTest)';
    end
    
    % if expand box
    switch warpType
        case 'sq2rect'
            trialsExpand = [randsample(locRange(1):diff(spacing(1:2)):locRange(2)*2,nTrials*.75,'true'); randsample(spacing,nTrials*.75,'true')]';
        case 'rect2sq'
            trialsExpand = [randsample(spacing,nTrials*.75,'true'); randsample(spacing,nTrials*.75,'true')]'; %this doesn't work - actually, maybe this isn't a thing?
    end

    
    %initialise each cluster location
    mu = nan(nClus,2,nTrials);
    for clusterI = 1:nClus
        %     mu(clusterI,:,1) = -locRange + locRange.*2.*rand(1,2);  %random location in box - uniform distr from -1 to 1 (see help rand)
        %     mu(clusterI,:,1) = trials(randi(length(trials)),:);   %intiate each cluster with one data point - Forgy method
        
        % k means++
        % initializing using dataPtsTest - not trials. prob doesn't matter
        % but this way it's just any random points in the box
        %%%%%%
        
        if clusterI==1% random datapoint as 1st cluster
            mu(clusterI,:,1) = dataPtsTest(randi(length(dataPtsTest)),:); 
        end
        if clusterI~=nClus % no need update k+1
            clear distInit
            for iClus = 1:clusterI% loop over clusters that exist now
                distInit(:,iClus)=sum([mu(iClus,1,1)-dataPtsTest(:,1),  mu(iClus,2,1)-dataPtsTest(:,2)].^2,2); %squared euclid for k means
            end
            [indValsInit, indInit]=min(distInit,[],2); % find which clusters are points closest to
            
            distClus=[];
            for iClus = 1:clusterI
                indOrig(:,clusterI) = indInit==iClus;
                distClusTmp = sum([(mu(iClus,1,1)-dataPtsTest(indOrig(:,clusterI),1)), (mu(iClus,2,1)-dataPtsTest(indOrig(:,clusterI),2))].^2,2);
                distClus = [distClus; [distClusTmp, repmat(iClus,length(distClusTmp),1)]];
            end
            
            %need to keep track of the indices of the original dist variable - get the
            %datapoints that were the farthest from all clusters, get that cluster and see which datapoint that was relative to that cluster (since i just save the distance)
            distClusNorm = distClus(:,1)./sum(distClus(:,1));
            distClusPr   = cumsum(distClusNorm(:,1)); %get cumsum, then generate rand val from 0 to 1 and choose smallest val - the larger the dis, the more likely the rand value will lie between it and its previous value in a cumsum plot
            ind=find(rand(1)<distClusPr,1);% %find smallest value that is larger than the random value (0 to 1 uniform distr)
            
            %                     tmp(i)=distClus(ind,1); %testing if getting the right values; if plot, see that it should be lower pr for closer values, higher for larger. note if very few high distances, this hist will look normally distributed
            
            clusInd = distClus(ind,2); %find which is the closest cluster
            indDat = find(distInit(:,clusInd)==distClus(ind,1)); %find where the datapoint is in the original vector
            
            if size(indDat,1)>1
                indDat=indDat(randi(size(indDat,1),1));
            end
            mu(clusterI+1,:,1) = dataPtsTest(indDat,:);
        end
    end
    %%

    updatedC = nan(nTrials,1);
    epsMuAll = nan(nTrials,2);
    deltaMu  = zeros(nClus,2,nTrials);
    distTrl  = nan(nTrials,nClus);
    
    clusUpdates = zeros(nClus,2); %acutally starting at 0 is OK, since there was no momentum from last trial
    
%     clusUpdAll{nClus}=[];
    for iTrl=1:nTrials
            
            %if change size of box half way
            if iTrl == nTrials*.75 && warpBox
                trials(nTrials*.75+1:end,:) = trialsExpand;
            end
            
            %compute distances
            dist2Clus = sqrt(sum([mu(:,1,iTrl)'-trials(iTrl,1); mu(:,2,iTrl)'-trials(iTrl,2)].^2)); % vectorising euclid dist - sqrt(sum((a-b).^2)), since can't use xval method

            %stochastic update - sel 1 of the closest clusters w random element - stochastic parameter c - large is deterministic, 0 - random            

            if stochasticType %stochastic update
                if stochasticType==1
%                     beta=c*(iTrl-1);         % so this gets bigger, and more deterministic with more trials
                    beta=c*(iTrl+nTrials/50);  % maybe no need to be so stochatic at start
                elseif stochasticType==2
%                     beta=c*(iTrl-1);
                    beta=c*(iTrl+nTrials/50);  % maybe no need to be so stochatic at start
                    if beta >= c*500 %it might be worth checking if this depends on nClusters - distances will change
                        beta = c*500;
                    end
                elseif stochasticType==3
                    beta = c*500; % as a function of c %prev:.185; 
                end
                dist2Clus2 = exp(-beta.*dist2Clus)./ sum(exp(-beta.*dist2Clus));
                distClusPr = cumsum(dist2Clus2);
                closestC=find(rand(1)<distClusPr,1);
                if isempty(closestC)  %if beta is too big, just do deterministic update (already will be near deterministic anyway)
                    closestC=find(min(dist2Clus)==dist2Clus);
                end
                cParams.betaAll(iTrl)       = beta;
            else %deterministic update
                closestC=find(min(dist2Clus)==dist2Clus);
            end
            if numel(closestC)>1 %if more than 1, randomly choose one
                closestC = randsample(closestC,1);
            end
            
            %if stochastic, closestC might not be if more than 1 actual
            %closest, have to check if stoch update chose one of them
            actualClosestC = find(min(dist2Clus)==dist2Clus); 
            closestMatch = zeros(1,length(actualClosestC));
            for iC = 1:length(actualClosestC)
                closestMatch(iC) = actualClosestC(iC) == closestC;
            end
            cParams.closestChosen(iTrl) = any(closestMatch); %was (one of) the closest cluster selected?
            cParams.closestDist(iTrl)   = dist2Clus(closestC);

            %log which cluster has been updated
            updatedC(iTrl) = closestC;

            epsMu = epsMuOrig;
            epsMuAll(iTrl,:) = [epsMu,closestC]; 

            %update (with momemtum-like parameter)
            deltaMu(closestC,1,iTrl) = ((1-alpha)*(epsMu*(trials(iTrl,1)-mu(closestC,1,iTrl))))+(alpha*clusUpdates(closestC,1));
            deltaMu(closestC,2,iTrl) = ((1-alpha)*(epsMu*(trials(iTrl,2)-mu(closestC,2,iTrl))))+(alpha*clusUpdates(closestC,2));
            
            clusUpdates(closestC,1)=deltaMu(closestC,1,iTrl);
            clusUpdates(closestC,2)=deltaMu(closestC,2,iTrl);

            deltaMuVec = zeros(nClus,2);
            deltaMuVec(closestC,:) = deltaMu(closestC,:,iTrl); % only update winner
            
            %weighted neighbour update
            % - need to change deltaMu to also compute changes from last
            % trial, not just closest cluster
            % - amount neighbors update should be proportion to distance;
            % use a beta value to make the neighbours update less; can try
            % a few of these; cumsum like above?
            % - note that this will interact with the momentum thing; if no
            % momentum, then fine, but if not, the neighbors will move less
            % esp if they didn't move much in the last few trials. this
            % interaction might be problematic since previously only takes
            % into account previous actual proper sized updates, not
            % neighbourhood ones... should the neighbourhood updates ignore
            % the momentum param..??? or only run without momentum?
            
            % update mean estimates
            if iTrl~=nTrials %no need to update for last trial +1)
                mu(:,1,iTrl+1) = mu(:,1,iTrl) + deltaMuVec(:,1);
                mu(:,2,iTrl+1) = mu(:,2,iTrl) + deltaMuVec(:,2);
            end        
        
        
        % compute sse on each trial with respect to 'all trials' - independent data - looks similar if you change 'dataPtsTest' to 'trials'; also useful if want to test less vs more trials on training so that you equate SSE by number of test points
        for iClus = 1:size(mu,1)
            distTrl(:,iClus)=sum([mu(iClus,1,iTrl)-dataPtsTest(:,1), mu(iClus,2,iTrl)-dataPtsTest(:,2)].^2,2);
        end
        [indValsTrl, indTrl]=min(distTrl,[],2); % find which clusters are points closest to
        
        for iClus = 1:size(mu,1)
            sseTrl(iClus,iTrl,iterI)=sum(sum([mu(iClus,1,iTrl)-dataPtsTest(indTrl==iClus,1), mu(iClus,2,iTrl)-dataPtsTest(indTrl==iClus,2)].^2,2)); %distance from each cluster from training set to datapoints closest to that cluster
        end
        tsseTrls(iTrl,iterI)=sum(sseTrl(:,iTrl,iterI));

        
    end
    muEnd(:,:,iterI)=mu(:,:,end);
    muAll(:,:,:,iterI) = mu;
end

%just save the best and worst 3 (since when nTrials are high, files get big)
% tsseIter = tsseTrls(end,:); %only use last trial

% xTrials = round(nTrials/3); %500;
% tsseIter = mean(tsseTrls(end-xTrials:end,:)); %average over last x trials 
% [indVal, indSSE1] = sort(tsseIter);
% [y, indSSE2] = sort(indSSE1); %use indSSE2 to sort for plotting
% 
% if size(muEnd,3) >= 6,
%     bestWorst3=[1,2,3,nIter-2,nIter-1,nIter];
%     muEndBest = nan(nClus,2,length(bestWorst3)); % save best.worse 3 clusters at end of learning
%     muAllBest = nan(nClus,2,nTrials,length(bestWorst3)); % save trials to plot over time
%     for iterI = 1:length(bestWorst3),
%         muEndBest(:,:,iterI) = muEnd(:,:,indSSE2==bestWorst3(iterI));
%         muAllBest(:,:,:,iterI) = muAll(:,:,:,indSSE2==bestWorst3(iterI));
%     end
% else
%     muEndBest = muEnd;
%     muAllBest = muAll;
% end
% 
% % also get top3 bottom3 for the spread outness measure?
% 
% %sort according to spread now
% tsseVarIter = mean(stdAcrossClus(end-xTrials:end,:)); %average over last x trials 
% [indVal, indSSEvar1] = sort(tsseVarIter);
% [y, indSSEvar2] = sort(indSSEvar1); %use indSSE2 to sort for plotting


%... another muEndBest and muAllBest variable?


%few things to note when updating this bit
% - if not saving all clusters maps, no need to save all sseTrl, etc.?
% - better: save all the cluster maps (final cluster positions, averaged
% over 1k, 5k, 10k, 15k, 20k, 25k. Also could save a 'downsampled' version;  with spacing of 100/500 trials. 
% - in any case, still save all for top and bottom 3


% to save:
% - all iters, all cluster positions over time

% - sort SSE according without averaging, last 1k, 5k, 10k, 15k, 20k, 25k 
% AND their indices: indSSE1  and indSSE2

end