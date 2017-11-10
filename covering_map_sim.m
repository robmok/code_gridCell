function [muEndBest, muAllBest, tsseTrls,sseTrl,epsMuAll,deltaMu] = covering_map_sim(nClus,locRange,box,warpType,epsMuOrig,nTrials,nIter,warpBox,alpha,c)

spacing=linspace(-1,1,101); 
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
    
    switch box
        case 'square'
            trials      = [randsample(spacing,nTrials,'true'); randsample(spacing,nTrials,'true')]';
            dataPtsTest = [randsample(linspace(-locRange,locRange,101),nTrialsTest,'true'); randsample(linspace(-locRange,locRange,101),nTrialsTest,'true')]'; % random points in a box
        case 'rect'
            trials      = [randsample(-1:diff(spacing(1:2)):2,nTrials,'true'); randsample(spacing,nTrials,'true')]';
            dataPtsTest = [randsample(-1:diff(spacing(1:2)):2,nTrialsTest,'true'); randsample(spacing,nTrialsTest,'true')]';
        case 'trapz'
            trapY=trapmf(spacing,[spacing(1), spacing(round(length(spacing)*.25)), spacing(round(length(spacing)*.75)),spacing(end)]);
            trapX=spacing;
            trapPts=[];
            for i=1:length(trapY),
               trapPts = [trapPts, [repmat(trapX(i),1,length(0:stepSize:trapY(i))); 0:stepSize:trapY(i)]];
            end
            trapPts(2,:)=trapPts(2,:).*2-1; %put it back into -1 to 1            
            % use this to select from the PAIR in trapPts
            trialInd     = randi(length(trapPts),nTrials,1);
            trials       = trapPts(:,trialInd)';
            trialIndTest = randi(length(trapPts),nTrials,1);
            dataPtsTest  = trapPts(:,trialIndTest)';
            
        case 'trapzSq' %probably need a more narrow trapezium!
            trapY=trapmf(spacing,[spacing(1), spacing(round(length(spacing)*.25)), spacing(round(length(spacing)*.75)),spacing(end)]);
            trapX=spacing;
            trapPts=[];
            for i=1:length(trapY),
               trapPts = [trapPts, [repmat(trapX(i),1,length(0:stepSize:trapY(i))); 0:stepSize:trapY(i)]];
            end
            %make square box attached to it
            sqX=spacing;
            sqY=spacing(1:floor(length(spacing)/2));
            for i=1:length(sqX),
                for j=1:length(sqY),
                trapPts = [trapPts, [sqX(i); sqY(j)]];
                end
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
            trialsExpand = [randsample(-1:diff(spacing(1:2)):2,nTrials*.75,'true'); randsample(spacing,nTrials*.75,'true')]';
        case 'rect2sq'
            trialsExpand = [randsample(spacing,nTrials/2,'true'); randsample(spacing,nTrials/2,'true')]'; %this doesn't work - actually, maybe this isn't a thing?
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
        
        if clusterI==1,% random datapoint as 1st cluster
            mu(clusterI,:,1) = dataPtsTest(randi(length(dataPtsTest)),:); 
        end
        if clusterI~=nClus, % no need update k+1
            clear distInit
            for iClus = 1:clusterI,% loop over clusters that exist now
                distInit(:,iClus)=sum([mu(iClus,1,1)-dataPtsTest(:,1),  mu(iClus,2,1)-dataPtsTest(:,2)].^2,2); %squared euclid for k means
            end
            [indValsInit, indInit]=min(distInit,[],2); % find which clusters are points closest to
            
            distClus=[];
            for iClus = 1:clusterI,
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
            
            if size(indDat,1)>1,
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
    
    clusUpdates = zeros(nClus,2)+.01; %acutally starting at 0 is OK, since there was no momentum from last trial
    
    clusUpdAll{nClus}=[];
    
    
    for iTrl=1:length(trials)
            
            %if change size of box half way
            if iTrl == nTrials*.75 && warpBox,
                trials(nTrials*.75+1:end,:) = trialsExpand;
            end
            
            %compute distances
            dist2Clus = sqrt(sum([mu(:,1,iTrl)'-trials(iTrl,1); mu(:,2,iTrl)'-trials(iTrl,2)].^2)); % vectorising euclid dist - sqrt(sum((a-b).^2)), since can't use xval method

%             %deterministic update
%             closestC=find(min(dist2Clus)==dist2Clus);
%             if numel(closestC)>1, %if more than 1, randomly choose one
%                 closestC = randsample(closestC,1);
%             end


            %stochastic update - sel 1 of the closest clusters w random element - stochastic parameter c - large is deterministic, 0 - random            
            
%             c = .0025; %constant - the larger, the more deterministic            
            % because this is trial by trial (not like k means), constant
            % must be much smaller than their recommendation (2). need it
            % to be much slower
            
%             beta=c*(iTrl-1); % so this gets bigger, and more deterministic with more trials
            beta=c*(iTrl+nTrials/10); %new - maybe no need to be so stochatic at start -BUT should this AND c be determined by nTrials?
            
            %or keep beta constant so always a bit stochastic
            if c == 999, %hack just to test this
                c = 15; %2 k is deterministic
                beta = c;
            end
            
            dist2Clus2 = exp(-beta.*dist2Clus)./ sum(exp(-beta.*dist2Clus));
            distClusPr = cumsum(dist2Clus2);
            closestC=find(rand(1)<distClusPr,1);    
            
% %             testing and plotting to see what's best - try to define c
% %             according to nTrials
%             for iTrl=1:10000,
% %                 beta=c*(iTrl-1); % so this gets bigger, and more deterministic with more trials
% %                 beta=c*(iTrl+nTrials/10); %new - no need to be so stochatic at start
% 
%                 %constantly stochastic at some num
% %                 c = 5; %2k = determinstic
%                 beta = 15;
%                 
%                 dist2Clus2 = exp(-beta.*dist2Clus)./ sum(exp(-beta.*dist2Clus));
%                 distClusPr = cumsum(dist2Clus2);
%                 closestC=find(rand(1)<distClusPr,1);
%                 
%                 xx(iTrl) = dist2Clus2(closestC); % this allows us to plot if it selects the most probable ones or not, so slowly goes toward 1 = always choose the closest one - should be random initially then going to 1
%                 yy(iTrl) = dist2Clus(closestC); 
%             end
%             figure; plot(xx);
%             figure; plot(yy);             


            %log which cluster has been updated
            updatedC(iTrl) = closestC;

            epsMu = epsMuOrig;
            epsMuAll(iTrl,:) = [epsMu,closestC];

                
%             new - adaptive learning rate (momentum-like)
%             need to save each clusters'previous update and position
%             
%             momentumUpdate:
%             alpha = 0.2% ? - if alpha is 0, then same as now, and move a good amount, not weightnig previous trials. if alpha is .99, then
%             weighting the previous trial a lot and not moving much
%             
%             currUpd = epsMu*(currClusLoc-dataPts); % on current trial
%             deltaMuCurrTrial = (1-alpha)*currUpd + alpha*prevUpd;
            
            
            % - set alpha (e.g. 0.2)
            % - update closest cluster
            % - save each cluster's update (save for each update/trial, or
            % just save the last update)
            
            
            deltaMu(closestC,1,iTrl) = ((1-alpha)*(epsMu*(trials(iTrl,1)-mu(closestC,1,iTrl))))+(alpha*clusUpdates(closestC,1));
            deltaMu(closestC,2,iTrl) = ((1-alpha)*(epsMu*(trials(iTrl,2)-mu(closestC,2,iTrl))))+(alpha*clusUpdates(closestC,2));
            
            clusUpdates(closestC,1)=deltaMu(closestC,1,iTrl);
            clusUpdates(closestC,2)=deltaMu(closestC,2,iTrl);
            
%             clusUpdAll{closestC}(length(clusUpdAll{closestC})+1)=deltaMu(closestC,1,iTrl);
%             %outputted this to check how each cluster updated
            
            
            deltaMuVec = zeros(nClus,2);
            deltaMuVec(closestC,:) = deltaMu(closestC,:,iTrl); % only update winner
            
            % update mean estimates
            if iTrl~=length(trials) %no need to update for last trial +1)
                mu(:,1,iTrl+1) = mu(:,1,iTrl) + deltaMuVec(:,1);
                mu(:,2,iTrl+1) = mu(:,2,iTrl) + deltaMuVec(:,2);
            end        
        
        
        % compute sse on each trial with respect to 'all trials' - independent data - looks similar if you change 'dataPtsTest' to 'trials'; also useful if want to test less vs more trials on training so that you equate SSE by number of test points
        for iClus = 1:size(mu,1),
            distTrl(:,iClus)=sum([mu(iClus,1,iTrl)-dataPtsTest(:,1), mu(iClus,2,iTrl)-dataPtsTest(:,2)].^2,2);
        end
        [indValsTrl, indTrl]=min(distTrl,[],2); % find which clusters are points closest to
        
        for iClus = 1:size(mu,1),
            sseTrl(iClus,iTrl,iterI)=sum(sum([mu(iClus,1,iTrl)-dataPtsTest(indTrl==iClus,1), mu(iClus,2,iTrl)-dataPtsTest(indTrl==iClus,2)].^2,2)); %distance from each cluster from training set to datapoints closest to that cluster
        end
        tsseTrls(iTrl,iterI)=sum(sseTrl(:,iTrl,iterI));

        
    end
    muEnd(:,:,iterI)=mu(:,:,end);
    muAll(:,:,:,iterI) = mu;
end

%just save the best and worst 3 (since when nTrials are high, files get big)
tsseIter = tsseTrls(end,:);
[indVal, indSSE1] = sort(tsseIter);
[y, indSSE2] = sort(indSSE1); %use indSSE2 to sort for plotting

if size(muEnd,3) >= 6,
    bestWorst3=[1,2,3,nIter-2,nIter-1,nIter];
    muEndBest = nan(nClus,2,length(bestWorst3)); % save best.worse 3 clusters at end of learning
    muAllBest = nan(nClus,2,nTrials,length(bestWorst3)); % save trials to plot over time
    for iterI = 1:length(bestWorst3),
        muEndBest(:,:,iterI) = muEnd(:,:,indSSE2==bestWorst3(iterI));
        muAllBest(:,:,:,iterI) = muAll(:,:,:,indSSE2==bestWorst3(iterI));
    end
else
    muEndBest = muEnd;
    muAllBest = muAll;
end

end