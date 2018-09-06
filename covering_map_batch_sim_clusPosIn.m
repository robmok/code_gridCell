function [densityPlot,densityPlotActNorm,gA,gW,gA_actNorm,gW_actNorm,rSeed,muAll,trials] = covering_map_batch_sim_clusPosIn(clusPos,nClus,locRange,epsMuOrig,epsMuTrapz,nTrialsOrig,nTrials,batchSize,nIter,dat,annEps,jointTrls)

spacing=linspace(locRange(1),locRange(2),locRange(2)+1); 

gaussSmooth=1; %smoothing for density map
imageFilter = fspecial('gaussian',5,gaussSmooth); %this is default for imgaussfilt

nBatch=round(nTrials/batchSize); 
batchSize = floor(batchSize); % when have decimal points, above needed
nBatchOrig = nTrialsOrig/batchSize;

%if decrease learning rate over time: 1/(1+decay+timestep); Decay - set a param
if annEps
%     eps stays high till 1/annEpsDecay batches - 
    annC = (1/nBatchOrig)/nBatchOrig; % 1/annC*nBatch = nBatch: constant to calc 1/annEpsDecay
    if epsMuOrig == 0.25
        annEpsDecay = annC*(nBatchOrig*100); %  epsMuOrig =.25;ends with .0025
    end
else %if not annEps, use fixed slower learning rate
    epsMuOrig=epsMuTrapz;
end

%compute gridness over time %20 timepoints; - 21 sets now - last one is last quarter
fromTrlI  = round(nTrials.*.75);
toTrlN    = nTrials;

if nargout > 7
    muAll            = nan(nClus,2,nBatch+1,nIter);
end
nSets              = length(fromTrlI);

gA = nan(nSets,nIter,9); %if saving the 5 r values, 9. if not, 4.
gW = nan(nSets,nIter,9);
gA_actNorm = nan(nSets,nIter,9);
gW_actNorm = nan(nSets,nIter,9);

%if trapz - compute gridness of left/right half of boxes too
if strcmp(dat(1:4),'trap')
    gA = nan(nSets,nIter,9,3);
    gW = nan(nSets,nIter,9,3);
    gA_actNorm = nan(nSets,nIter,9,3);
    gW_actNorm = nan(nSets,nIter,9,3);
end

%set densityPlot array size
b = length(spacing);
h = length(spacing);
%new correct trapzKrupic scale
hLeft  = 17;% 12;
hRight = 33;% - 33 = start from 18 from left % 27;

% now setting/refreshing some of these over sets below
densityPlot        = zeros(b,h,nSets,nIter);
densityPlotActNorm = zeros(b,h,nSets,nIter);

for iterI = 1:nIter
    
    fprintf('iter %d \n',iterI);
    catsInfo =[];
    [trials,~, rSeed(iterI),ntInSq] = createTrls(dat,nTrials,locRange,0,jointTrls,catsInfo);

    %initialise each cluster location - at learned positions
    mu = nan(nClus,2,nBatch+1);
    mu(:,:,1) = clusPos(:,:,iterI); % load in cluster positions from prev 
    
    deltaMu  = zeros(nClus,2,nBatch);        
    for iBatch=1:nBatch
        batchInd=batchSize*(iBatch-1)+1:batchSize*(iBatch-1)+batchSize; %trials to average over
        trls2Upd = trials(batchInd,:); %trials to use this batch

            %compute distances - vectorise both clusters and trials (in batch)
            dist2Clus = sqrt(sum(reshape([mu(:,1,iBatch)'-trls2Upd(:,1), mu(:,2,iBatch)'-trls2Upd(:,2)].^2,batchSize,nClus,2),3));% reshapes it into batchSize x nClus x 2 (x and y locs)
            
            %update cluster positions
            closestC = nan(1,batchSize);
            for iTrlBatch = 1:batchSize
                closestTmp = find(min(dist2Clus(iTrlBatch,:))==dist2Clus(iTrlBatch,:));
                
                if numel(closestTmp)>1
                    closestC(iTrlBatch) = randsample(closestTmp,1);
                else
                    closestC(iTrlBatch) = closestTmp;
                end
            end
            
            %learning rate
            if annEps %if use annealed learning rate
                epsMu = epsMuOrig/(1+(annEpsDecay*(iBatch+nBatchOrig))); %new
            else
                epsMu = epsMuOrig;
            end
            
            %batch update
            for iClus = 1:nClus
                updInd = closestC==iClus;
                if any(nnz(updInd)) %if not, no need to update
                deltaMu(iClus,1,iBatch) = nanmean(epsMu*(trls2Upd(updInd,1)-mu(iClus,1,iBatch)));
                deltaMu(iClus,2,iBatch) = nanmean(epsMu*(trls2Upd(updInd,2)-mu(iClus,2,iBatch)));
                end
            end

            % update mean estimates
            mu(:,1,iBatch+1) = mu(:,1,iBatch) + deltaMu(:,1,iBatch);
            mu(:,2,iBatch+1) = mu(:,2,iBatch) + deltaMu(:,2,iBatch);
    end
    if nargout > 7
        muAll(:,:,:,iterI)      = mu;
    end

    % densityplot over time (more samples)
    for iSet = 1:nSets
        densityPlotClus      = zeros(b,h,nClus);
        clus = round(mu(:,:,round(toTrlN(iSet)./batchSize))); %mu is in batches     
        for iClus=1:nClus
            ntNanInd = squeeze(~isnan(clus(iClus,1,:)));
            clusTmp = []; %clear else dimensions change over clus/sets
            clusTmp(1,:) = squeeze(clus(iClus,1,ntNanInd)); %split into two to keep array dim constant - when only 1 location, the array flips.
            clusTmp(2,:) = squeeze(clus(iClus,2,ntNanInd));
            for iTrlUpd=1:size(clusTmp,2)
                densityPlotClus(clusTmp(1,iTrlUpd),clusTmp(2,iTrlUpd),iClus) = densityPlotClus(clusTmp(1,iTrlUpd),clusTmp(2,iTrlUpd),iClus)+1;
            end
        end
        densityPlotTmp = nansum(densityPlotClus,3);

        %turn 0s outside of the shape into nans
        if ~strcmp(dat(1:2),'sq') % if circ/trapz
            for i=1:length(ntInSq)
                densityPlotTmp(ntInSq(i,1)+1,ntInSq(i,2)+1) = nan;
            end
        end
        densityPlot(:,:,iSet,iterI) = densityPlotTmp;
        densityPlotSm = nanconv(densityPlotTmp,imageFilter, 'nanout'); %new smooths ignoring nans
        
        if ~strcmp(dat(1:3),'cat') %if finding cats, won't be gridlike
            %compute autocorrmap
            aCorrMap = ndautoCORR(densityPlotSm);
            %compute gridness
            [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
            gA(iSet,iterI,:,1) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
            [g,gdataA] = gridSCORE(aCorrMap,'wills',0);
            gW(iSet,iterI,:,1) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
            
            %trapz left right of box
            if  strcmp(dat(1:4),'trap')
                %left half of box
                aCorrMap = ndautoCORR(densityPlotSm(:,1:hLeft));
                [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
                gA(iSet,iterI,:,2) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
                [g,gdataA] = gridSCORE(aCorrMap,'wills',0);
                gW(iSet,iterI,:,2) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
                
                %right half of box
                aCorrMap = ndautoCORR(densityPlotSm(:,h-hRight+1:end));
                [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
                gA(iSet,iterI,:,3) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
                [g,gdataA] = gridSCORE(aCorrMap,'wills',0);
                gW(iSet,iterI,:,3) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
            end
        end
    end
end
end