function [muEndBest, muAllBest, tsseTrls, epsMuAll] = covering_map_sim(nClus,locRange,box,warpType,epsMuOrig,deltaEpsMu,nTrials,nIter,warpBox,resetEps)

spacing=linspace(-1,1,101); 

%for crossvalidation
nTrialsTest = nTrials;
% nTrialsTest = 5000; %keep it constant if want to compare with nTrials above
dataPtsTest = [randsample(linspace(-locRange,locRange,101),nTrialsTest,'true'); randsample(linspace(-locRange,locRange,101),nTrialsTest,'true')]'; % random points in a box

muAll    = nan(nClus,2,nTrials,nIter);
tsseTrls = nan(nTrials,nIter);
muEnd    = nan(nClus,2,nIter);



for iterI = 1:nIter
    
    fprintf('iter %d \n',iterI);
    epsMu = epsMuOrig; %revert learning rate back to original if reset
    
    %reset trials after iteration (problem esp when box expands)
    %%%%%%
    %Q - should i want the datapoints to be the same each time, so all that
    %changes are the starting points?
    %%%%%%
    switch box
        case 'square'
            trials=[randsample(spacing,nTrials,'true'); randsample(spacing,nTrials,'true')]';
        case 'rect'
            trials=[randsample(-1:diff(spacing(1:2)):2,nTrials,'true'); randsample(spacing,nTrials,'true')]';
    end
    
    % if expand box
    switch warpType
        case 'sq2rect'
            trialsExpand = [randsample(-1:diff(spacing(1:2)):2,nTrials/2,'true'); randsample(spacing,nTrials/2,'true')]';
        case 'rect2sq'
            trialsExpand = [randsample(spacing,nTrials/2,'true'); randsample(spacing,nTrials/2,'true')]'; %this doesn't work - actually, maybe this isn't a thing?
    end

    
    %initialise each cluster location
    mu = nan(nClus,2,nTrials);
    for clusterI = 1:nClus
        %     mu(clusterI,:,1) = -locRange + locRange.*2.*rand(1,2);  %random location in box - uniform distr from -1 to 1 (see help rand)
        mu(clusterI,:,1) = trials(randi(length(trials)),:);   %intiate each cluster with one data point - Forgy method
        
        % k means++?
        
    end
    
    %%
    clusInd  = 1:nClus;
    epsMuVec = zeros(nClus,1);
%     updatedC = nan(nTrials,1);
    magChg   = nan(nTrials,1);
    epsMuAll = nan(nTrials,2);
    deltaMu  = nan(nClus,2,nTrials);
    distTrl  = nan(nTrials,nClus);
    
    updatedC = ones(nClus,1);
    
    for iTrl=1:length(trials)
        for iClus = 1:nClus,
            
            %if change size of box half way
            if iTrl == nTrials/2 && warpBox,
                trials(nTrials/2+1:end,:) = trialsExpand;
            end
            
            dist2Clus = sqrt(sum([mu(:,1,iTrl)'-trials(iTrl,1); mu(:,2,iTrl)'-trials(iTrl,2)].^2)); % vectorising euclid dist - sqrt(sum((a-b).^2)), since can't use xval method
            closestC=find(min(dist2Clus)==dist2Clus);
            if numel(closestC)>1, %if more than 1, randomly choose one
                closestC = randsample(closestC,1);
            end
            
            
            %moved learning rate change to up here - now update epsMu(iTrl), not epsMu(iTrl-1)
            
            %need to think how to update per cluster - since before just
            %defined x proprtion of trials. i suppose this is ok too - since
            %resetting them all a bit is better than just resetting one (prob
            %doesnt do anything).
            % the Q is whether should reset them all to a similar level - now
            % the ones that have been updated the most will have reset more and
            % move less...
            %         if iTrl == nTrials/2,
            %             updatedC = [];
            %         end
            
            %log which cluster has been updated
%             updatedC(iTrl) = closestC;
            
            updatedC(closestC)=updatedC(closestC)+1;
            
            
            if iTrl~=1,
            
                
                %reset learning rate here now - reset to a set value
                if resetEps == 1,
                    if iTrl > nTrials*.5
                        updatedC(:,1) = 100;
                    end
                elseif resetEps == 2,
                    if iTrl > nTrials*.25 && iTrl <= nTrials*.5

                    
                    elseif iTrl > nTrials*.5 && iTrl <= nTrials*.75

                    elseif iTrl > nTrials*.75 && iTrl <= nTrials

                    end
                end
                
                %slow down learning rate over time
%                 magChg(iTrl)=nnz(updatedC==closestC); %magnitude of change proportional to how many times the cluster has moved
                magChg(iTrl)=updatedC(closestC); 
                %             epsMu(iTrl+1) = epsMu(1); %no slowing down
                %             epsMu(iTrl) = epsMu(iTrl-1)*deltaEpsMu;
                
                
                %resetting learning rate
                %still need to think if the amount i'm resetting to is
                %calculated correctly, and whether these are good numbers (prob
                %doesn't matter too much, but sth like high, med low sounds OK)
%                 if resetEps == 1,
%                     if iTrl > nTrials*.5
%                         magChg(iTrl)=nnz(updatedC==closestC)*.5;
%                     end
%                 elseif resetEps == 2,
%                     if iTrl > nTrials*.25 && iTrl <= nTrials*.5
%                         magChg(iTrl)=nnz(updatedC==closestC)*.75;
%                     elseif iTrl > nTrials*.5 && iTrl <= nTrials*.75
%                         magChg(iTrl)=nnz(updatedC==closestC)*.25;
%                     elseif iTrl > nTrials*.75 && iTrl <= nTrials
%                         magChg(iTrl)=nnz(updatedC==closestC)*.75;
%                     end
%                 end
                
                
                epsMu = epsMuOrig*deltaEpsMu^magChg(iTrl); %squared because it's deltaMu (e.g. .99) to the power of nTimes it was updated, making it slightly smaller each time
            else
                epsMu = epsMuOrig;
            end
            epsMuAll(iTrl,:) = [epsMu,closestC];
            
            epsMuVec = zeros(nClus,1);
            epsMuVec(closestC) = epsMu;
            deltaMu(:,1,iTrl) = epsMuVec.*(trials(iTrl,1)-mu(:,1,iTrl));
            deltaMu(:,2,iTrl) = epsMuVec.*(trials(iTrl,2)-mu(:,2,iTrl));
            % update mean estimates
            if iTrl~=length(trials) %no need to update for last trial +1)
                mu(:,1,iTrl+1) = mu(:,1,iTrl) + deltaMu(:,1,iTrl);
                mu(:,2,iTrl+1) = mu(:,2,iTrl) + deltaMu(:,2,iTrl);
            end
        end
        
        % compute sse on each trial with respect to 'all trials' - independent data - looks similar if you change 'dataPtsTest' to 'trials'; also useful if want to test less vs more trials on training so that you equate SSE by number of test points
        for iClus = 1:size(mu,1),
            distTrl(:,iClus)=sum([mu(iClus,1,iTrl)-dataPtsTest(:,1), mu(iClus,2,iTrl)-dataPtsTest(:,2)].^2,2);
        end
        [indValsTrl, indTrl]=min(distTrl,[],2); % find which clusters are points closest to
        
        sseTrl = zeros(iClus,1);
        for iClus = 1:size(mu,1),
            sseTrl(iClus)=sum(sum([mu(iClus,1,iTrl)-dataPtsTest(indTrl==iClus,1), mu(iClus,2,iTrl)-dataPtsTest(indTrl==iClus,2)].^2,2)); %distance from each cluster from training set to datapoints closest to that cluster
        end
        tsseTrls(iTrl,iterI)=sum(sseTrl);
        
        
        % idea: resetting the learning rate at some point to 'adapt' to new
        % environment. useful if morph the box shape in real time.
        
        %     % might also be interesting to know if this actually improves SSE, since
        %     % allows the first go/first few gos where things get stuck to improve
        %     if iTrl == round(length(trials)*.25) || iTrl == round(length(trials)*.5) || iTrl == round(length(trials)*.75)
        %     if iTrl == round(length(trials)*.5)
        %         epsMu(iTrl+1)=epsMu(1); %reset learning rate at some point, see what happens
        %     end
        %
        
        % %
        %     %     %could also reset to half the previous original learning rate each time...
        %     if iTrl == round(length(trials)*.25)
        %         epsMu(iTrl+1)=epsMu(1)*.75;
        %     elseif iTrl == round(length(trials)*.5)
        %         epsMu(iTrl+1)=epsMu(1)*.5;
        %     elseif iTrl == round(length(trials)*.75)
        %         epsMu(iTrl+1)=epsMu(1)*.25;
        %     end
        
        %     if iTrl == round(length(trials)*.25)
        %         epsMu(iTrl+1)=epsMu(1)*.75;
        %     elseif iTrl == round(length(trials)*.33)
        %         epsMu(iTrl+1)=epsMu(1)*.66;
        %     elseif iTrl == round(length(trials)*.5)
        %         epsMu(iTrl+1)=epsMu(1)*.5;
        %     elseif iTrl == round(length(trials)*.66)
        %         epsMu(iTrl+1)=epsMu(1)*.33;
        %     elseif iTrl == round(length(trials)*.75)
        %         epsMu(iTrl+1)=epsMu(1)*.25;
        %     end
        
        %     %put a limit of lowest learning rate
        %     if epsMu(iTrl+1)<=.1,
        %         epsMu(iTrl+1) = .1;
        %     end
        
    end
    muEnd(:,:,iterI)=mu(:,:,end);
    muAll(:,:,:,iterI) = mu;
end

%just save the best and worst 3 (since when nTrials are high, files get big)
tsseIter = tsseTrls(end,:);
[indVal, indSSE1] = sort(tsseIter);
[y, indSSE2] = sort(indSSE1); %use indSSE2 to sort for plotting

bestWorst3=[1,2,3,nIter-2,nIter-1,nIter];
muEndBest = nan(nClus,2,length(bestWorst3)); % save best.worse 3 clusters at end of learning
muAllBest = nan(nClus,2,nTrials,length(bestWorst3)); % save trials to plot over time
for iterI = 1:length(bestWorst3),
    muEndBest(:,:,iterI) = muEnd(:,:,indSSE2==bestWorst3(iterI));
    muAllBest(:,:,:,iterI) = muAll(:,:,:,indSSE2==bestWorst3(iterI));
end

end