function [densityPlot,clusMu,muAvg,gA,gW,muAll]  = covering_map_sim_neigh(nClus,locRange,box,warpType,epsMuOrig,nTrials,nIter,warpBox,trials,beta)

%neighbour-weighted update 
% - no momentum, no stochastic update

spacing=linspace(locRange(1),locRange(2),locRange(2)+1); 
stepSize=diff(spacing(1:2));

%for crossvalidation - 
nTrialsTest = nTrials;
% nTrialsTest = 5000; %keep it constant if want to compare with nTrials above

gaussSmooth=1; %smoothing for density map


% muAll    = nan(nClus,2,nTrials,nIter);
% averaging over trials to compute density map
% 5k - 20k:25k, 25k:30k, 30k:35k, 35k:40k
% 10k - 20k:30k, 25:35k, 30k:40k,
% 15k - 20k:35k; 25k:40k
% 20k - 20k:40k
% fromTrlI = [2.0e+4, 2.5e+4, 3.0e+4, 3.5e+4, 2.0e+4, 2.5e+4, 3.0e+4, 2.0e+4, 2.5e+4, 3.0e+4];
% toTrlN   = [2.5e+4, 3.0e+4, 3.5e+4, 4.0e+4, 3.0e+4, 3.5e+4, 4.0e+4, 3.5e+4, 4.0e+4, 4.0e+4];

%new
% 5k  - 10k:15k, 20k:25k 35k:40k
% 10k - 10k:20k, 15:25k, 30k:40k,
% 15k - 10k:25k, 15k:30k; 25k:40k
% 20k - 10k:30k, 20k:40k
fromTrlI = [1.0e+4, 2.0e+4, 3.5e+4, 1.0e+4, 1.5e+4, 3.0e+4, 1.0e+4,  1.5e+4, 2.5e+4, 1.0e+4, 2.0e+4];
toTrlN   = [1.5e+4, 2.5e+4, 4.0e+4, 2.0e+4, 2.5e+4, 4.0e+4, 2.50e+4, 3.0e+4, 4.0e+4, 3.0e+4, 4.0e+4];
nSets = length(fromTrlI);
densityPlotClus  = zeros(length(spacing),length(spacing),nClus,nSets,nIter);
densityPlot      = zeros(length(spacing),length(spacing),nSets,nIter);
clusMu           = nan(nClus,2,nSets,nIter);
muAvg            = nan(nClus,2,nSets,nIter);
for iterI = 1:nIter
    
    fprintf('iter %d \n',iterI);
    epsMu = epsMuOrig; %revert learning rate back to original if reset
    
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
                        
            clusInd = distClus(ind,2); %find which is the closest cluster
            indDat = find(distInit(:,clusInd)==distClus(ind,1)); %find where the datapoint is in the original vector
            
            if size(indDat,1)>1
                indDat=indDat(randi(size(indDat,1),1));
            end
            mu(clusterI+1,:,1) = dataPtsTest(indDat,:);
        end
    end
    muUpd = mu;  %variable that only logs updated clusters mean value % - get initialised locs

    %%

    updatedC = nan(nTrials,1);
    epsMuAll = nan(nTrials,2);
    deltaMu  = zeros(nClus,2,nTrials);
    
    clusUpdates = zeros(nClus,2); %acutally starting at 0 is OK, since there was no momentum from last trial
    
    for iTrl=1:nTrials
            
            %if change size of box half way
            if iTrl == nTrials*.75 && warpBox
                trials(nTrials*.75+1:end,:) = trialsExpand;
            end
            
            %compute distances
            dist2Clus = sqrt(sum([mu(:,1,iTrl)'-trials(iTrl,1); mu(:,2,iTrl)'-trials(iTrl,2)].^2)); % vectorising euclid dist - sqrt(sum((a-b).^2)), since can't use xval method

            %deterministic update
            closestC=find(min(dist2Clus)==dist2Clus);
            if numel(closestC)>1 %if more than 1, randomly choose one
                closestC = randsample(closestC,1);
            end

            %log which cluster has been updated
            updatedC(iTrl) = closestC;

            epsMu = epsMuOrig;
%             epsMuAll(iTrl,:) = [epsMu,closestC]; 

            %softmax
            updW = exp(-beta.*dist2Clus)./ sum(exp(-beta.*dist2Clus));

            %manual check
%             updW = 0./dist2Clus; %last few do look a bit like this - need
%             to compare gridness

% %             updW = .00001./dist2Clus;
% %             updW = .00005./dist2Clus; % possibly starting to look like no weight?
% %             updW = .0001./dist2Clus; % possibly starting to look like no weight?
% %             updW = .0005./dist2Clus; %works well
%             updW = .001./dist2Clus; %works well
% %             updW = .005./dist2Clus; %works well - maybe a bit for a lower learning rate from .075?
% %             updW = .01./dist2Clus; %works well - tiny bit toward centr
% 
% %              updW = .015./dist2Clus; %tiny bit of gravity toward centre
% %             updW = .025./dist2Clus; %some gravity towards the centre
% %             updW = .05./dist2Clus; %starts going into the centre
%             
% %             updW = 1-max(distW)+distW; % so that closest clus weight = 1
%             updW(closestC) = 1;
            
            %update (no momemtum-like parameter)
            deltaMu(:,1,iTrl) = updW'.*(epsMu.*(trials(iTrl,1)-mu(:,1,iTrl)));
            deltaMu(:,2,iTrl) = updW'.*(epsMu.*(trials(iTrl,2)-mu(:,2,iTrl)));
            
%             % %checking
%             delta_w=updW'.*(epsMu.*(trials(iTrl,1)-mu(:,1,iTrl)));
%             delta_orig=(epsMu.*(trials(iTrl,1)-mu(:,1,iTrl)));
%             figure; scatter(delta_orig,delta_w);

%             deltaMuVec = zeros(nClus,2);
            deltaMuVec = deltaMu(:,:,iTrl);
            
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
                %log updated value only
                muUpd(closestC,1,iTrl+1)=mu(closestC,1,iTrl+1);
                muUpd(closestC,2,iTrl+1)=mu(closestC,2,iTrl+1);
            end        
    end
    muAll(:,:,:,iterI) = mu;
    
    for iSet = 1:nSets
        %compute density map
%         clus = round(mu(:,:,fromTrlI(iSet):toTrlN(iSet)));
%         ind=clus<=0; clus(ind)=1; %indices <= 0 make to 1
%         for iClus=1:nClus
%             for iTrl=1:size(clus,3)
%                 densityPlotClus(clus(iClus,1,iTrl),clus(iClus,2,iTrl),iClus,iSet,iterI) = densityPlotClus(clus(iClus,1,iTrl),clus(iClus,2,iTrl),iClus,iSet, iterI)+1;
%             end
%         end
        clus = round(muUpd(:,:,fromTrlI(iSet):toTrlN(iSet)));
        ind=clus<=0; clus(ind)=1; %indices <= 0 make to 1
        for iClus=1:nClus
            ntNanInd = squeeze(~isnan(clus(iClus,1,:)));
            clusTmp  = squeeze(clus(iClus,:,ntNanInd));
            nTrlsUpd(iClus,iSet,iterI)=nnz(ntNanInd);
            for iTrlUpd=1:size(clusTmp,2)
                densityPlotClus(clusTmp(1,iTrlUpd),clusTmp(2,iTrlUpd),iClus,iSet,iterI) = densityPlotClus(clusTmp(1,iTrlUpd),clusTmp(2,iTrlUpd),iClus,iSet, iterI)+1;
            end
        end     
        %now also compute clus means
        densityPlotClusSmth = zeros(length(spacing),length(spacing),nClus);
        for iClus=1:nClus
            %find peaks
            densityPlotClusSmth(:,:,iClus)=imgaussfilt(densityPlotClus(:,:,iClus,iSet,iterI),gaussSmooth);
            [peakX, peakY] = find(densityPlotClusSmth(:,:,iClus)==max(max((densityPlotClusSmth(:,:,iClus)))));
            if length(peakX)>1 || length(peakY)>1 %if more than one peak rand sel one; normally next to each other
                randInd=randi(length(peakX));
                peakX=peakX(randInd);
                peakY=peakY(randInd);
            end
            clusMu(iClus,:,iSet,iterI) = [peakX, peakY];
            
            %make combined (grid cell) plot, smooth
            densityPlot(:,:,iSet,iterI) = sum(densityPlotClus(:,:,:,iSet,iterI),3); %save this
            densityPlotSm = imgaussfilt(densityPlot(:,:,iSet,iterI),gaussSmooth);
        end
        %compute autocorrmap, no need to save
        aCorrMap = ndautoCORR(densityPlotSm);
        %compute gridness
        [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
        [g,gdataW] = gridSCORE(aCorrMap,'wills',0);

        gA(iSet,iterI,:) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius];
        gW(iSet,iterI,:) = [gdataW.g_score, gdataW.orientation, gdataW.wavelength, gdataW.radius];
        
        
        %save average cluster positions (to compare with above)
        muAvg(:,:,iSet,iterI) = mean(mu(:,:,fromTrlI(iSet):toTrlN(iSet)),3);
    end    
end
end