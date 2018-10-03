function [densityPlot,densityPlotActNorm,gA,gW,gA_actNorm,gW_actNorm,muInit,rSeed,clusDistB,muAll, trials] = covering_map_batch_sim(nClus,locRange,catsInfo,epsMuOrig,nTrials,batchSize,nIter,trials,useSameTrls,dat,annEps,jointTrls,actOverTime)
% run clustering algorithm with a batch update

spacing=linspace(locRange(1),locRange(2),locRange(2)+1); 
stepSize=diff(spacing(1:2)); 
gaussSmooth = 1; %smoothing for density map
imageFilter = fspecial('gaussian',5,gaussSmooth); % param for smoothing (for nanconv - smoothing that deals with nans; this param is default for imgaussfilt)
sigmaGauss = stepSize; %computing cluster activation
nBatch=round(nTrials/batchSize);
batchSize = floor(batchSize); % when have decimal points, above needed

% if decrease learning rate over time: 1/(1+decay+timestep); set decay
% param. eps stays high till 1/annEpsDecay batches
if annEps
    annC = (1/nBatch)/nBatch; % 1/annC*nBatch = nBatch: constant to calc 1/annEpsDecay
    if epsMuOrig == 0.25
        annEpsDecay = annC*(nBatch*100); %  epsMuOrig =.25;ends with .0025
    end
end

%compute gridness over time 20 timepoints; - 21 sets now - last one is last quarter
if actOverTime
    fromTrlI  = round([1,               nTrials.*.05+1,  nTrials.*.1+1,  nTrials.*.15+1,  nTrials.*.2+1,  nTrials.*.25+1,   nTrials.*.3+1,  nTrials.*.35+1,  nTrials.*.4+1,  nTrials.*.45+1,  nTrials.*.5+1,  nTrials.*.55+1, nTrials.*.6+1,  nTrials.*.65+1, nTrials.*.7+1, nTrials.*.75+1, nTrials.*.8+1,  nTrials.*.85+1,  nTrials.*.9+1,  nTrials.*.95+1, nTrials.*.75]);
    toTrlN    = round([nTrials.*.05,    nTrials.*.1,     nTrials.*.15,   nTrials.*.2,     nTrials.*.25,   nTrials.*.3,      nTrials.*.35,   nTrials.*.4,     nTrials.*.45,   nTrials.*.5,     nTrials.*.55,   nTrials.*.6,    nTrials.*.65,   nTrials.*.7,    nTrials.*.75,  nTrials.*.8,    nTrials.*.85,   nTrials.*.9,     nTrials.*.95,   nTrials,        nTrials]);
else %just compute last one
    fromTrlI  = round(nTrials.*.75);
    toTrlN    = nTrials;
end

if nargout > 9
    muAll            = nan(nClus,2,nBatch+1,nIter);
end
nSets              = length(fromTrlI);
muInit               = nan(nClus,2,nIter);

gA = nan(nSets,nIter,9); %if saving the extra 5 r values, 9. if not, 4.
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

%compute distance of each cluster to itself over time (batches)
clusDistB = nan(nSets-1,nIter);

if strcmp(dat(1:4),'trap') && length(dat)>10
    if strcmp(dat(1:11),'trapzKrupic')
        spacingTrapz = spacing(14:37);
        if length(dat)==11 % 'trapzKrupic'
            a = length(spacingTrapz);
            b = length(spacing);
            h = length(spacing);
        end
        % split trapz
        hLeft=17;% 1:17
        hRight=33;% - 33 means start from 18 from left, 18:50
    end 
else
    b=length(spacing);
    h=length(spacing);
end

densityPlot        = zeros(b,h,nSets,nIter);
densityPlotActNorm = zeros(b,h,nSets,nIter);
updatedC      = nan(1,nTrials);       
%start
for iterI = 1:nIter
    fprintf('iter %d \n',iterI);

    [trials,dataPtsTest, rSeed(iterI),ntInSq] = createTrls(dat,nTrials,locRange,useSameTrls,jointTrls,catsInfo);
    
    %initialise each cluster location  
    mu = nan(nClus,2,nBatch+1); 
    if ~strcmp(dat(1:3),'cat')
        mu(:,:,1) = dataPtsTest(randi(nTrials,nClus,1),:); %forgy method
    else
        mu(:,:,1) = reshape(randsample(locRange(1):diff(spacing(1:2)):locRange(2),nClus*2,'true'),nClus,2); %random
    end
    % mu(:,:,1) = kmplusInit(dataPtsTest,nClus); %kmeans++ initialisation
    
    muInit(:,:,iterI) = mu(:,:,1);
    actTrl = zeros(nClus,batchSize);
    %%
    
    deltaMu   = zeros(nClus,2,nBatch);
    actTrlAll = nan(nClus,batchSize,nBatch);

    for iBatch = 1:nBatch
        batchInd=batchSize*(iBatch-1)+1:batchSize*(iBatch-1)+batchSize; %trials to average over
        trls2Upd = trials(batchInd,:); %trials to use this batch

            %compute distances - vectorise both clusters and trials (in batch)
            dist2Clus = sqrt(sum(reshape([mu(:,1,iBatch)'-trls2Upd(:,1), mu(:,2,iBatch)'-trls2Upd(:,2)].^2,batchSize,nClus,2),3));% reshapes it into batchSize x nClus x 2 (x and y locs)
            
            %update cluster positions
            closestC = nan(1,batchSize);
            for iTrlBatch = 1:batchSize
                closestTmp = find(min(dist2Clus(iTrlBatch,:))==dist2Clus(iTrlBatch,:));
                
                if numel(closestTmp)>1 % if more than 1, randomly sel 1
                    closestC(iTrlBatch) = randsample(closestTmp,1);
                else
                    closestC(iTrlBatch) = closestTmp;
                end
                %compute activation
                actTrl(closestC(iTrlBatch),iTrlBatch)=mvnpdf(trls2Upd(iTrlBatch,:),mu(closestC(iTrlBatch),:,iBatch),eye(2)*sigmaGauss); % save only the winner
                
                %log which cluster has been updated
                updatedC((iBatch-1)*batchSize+iTrlBatch) = closestC(iTrlBatch);
            end
%             if ~strcmp(dat(1:4),'trap') %not computing for trapz
                actTrlAll(:,:,iBatch) = actTrl;
%             end
            
            %learning rate
            if annEps %if use annealed learning rate
                epsMu = epsMuOrig/(1+(annEpsDecay*iBatch)); 
            else
                epsMu = epsMuOrig;
            end
            
            %batch update - save all updates for each cluster for X trials, update - goes through each cluster. Compute distances, average, then update
            for iClus = 1:nClus
                updInd = closestC==iClus;
                if any(nnz(updInd)) %if not, no need to update
                    deltaMu(iClus,1,iBatch) = nanmean(epsMu*(trls2Upd(updInd,1)-mu(iClus,1,iBatch)));
                    deltaMu(iClus,2,iBatch) = nanmean(epsMu*(trls2Upd(updInd,2)-mu(iClus,2,iBatch)));
                end
            end
            % update
            mu(:,1,iBatch+1) = mu(:,1,iBatch) + deltaMu(:,1,iBatch);
            mu(:,2,iBatch+1) = mu(:,2,iBatch) + deltaMu(:,2,iBatch);
    end
    if nargout > 9 %save and output all positions over trials (variable is huge with high nIter)
        muAll(:,:,:,iterI)      = mu;
    end
    actAll  = reshape(actTrlAll,nClus,nTrials); %save trial-by-trial act over blocks, here unrolling it

    % compute densityPlotClus - density plot with each cluster in third dimension - more like a place cell map
    % densityplot over time (nSets)
    for iSet = 1:nSets
        densityPlotClus      = zeros(b,h,nClus);
        densityPlotAct       = zeros(b,h);
        densityPlotActUpd    = zeros(b,h); %start with ones, so won't divide by 0
        clus = round(mu(:,:,round(toTrlN(iSet)./batchSize))); %mu is in batches     
        for iClus=1:nClus
            ntNanInd = squeeze(~isnan(clus(iClus,1,:)));
            clusTmp = []; %clear else dimensions change over clus/sets
            clusTmp(1,:) = squeeze(clus(iClus,1,ntNanInd)); %split into two to keep array dim constant - when only 1 location, the array flips.
            clusTmp(2,:) = squeeze(clus(iClus,2,ntNanInd));
            if any(clusTmp>49) || any(clusTmp<1) %if out of box, don't add to densityPlot
                clusTmp=[];
            end
            for iTrlUpd=1:size(clusTmp,2)
                densityPlotClus(clusTmp(1,iTrlUpd)+1,clusTmp(2,iTrlUpd)+1,iClus) = densityPlotClus(clusTmp(1,iTrlUpd)+1,clusTmp(2,iTrlUpd)+1,iClus)+1;
            end
        end
        densityPlotTmp = nansum(densityPlotClus,3); % sum over third dimension to make it like a grid cell map

        %densityPlotActNorm
        for iTrl = fromTrlI(iSet):toTrlN(iSet)
            densityPlotAct(trials(iTrl,1)+1, trials(iTrl,2)+1)    = densityPlotAct(trials(iTrl,1)+1, trials(iTrl,2)+1)+ nansum(actAll(:,iTrl));
            densityPlotActUpd(trials(iTrl,1)+1, trials(iTrl,2)+1) = densityPlotActUpd(trials(iTrl,1)+1, trials(iTrl,2)+1)+1; %log nTimes loc was visited
        end
        
        if ~strcmp(dat(1:2),'sq')
            for i=1:length(ntInSq)% turn 0s outside of the shape (square) into nans
                densityPlotTmp(ntInSq(i,1)+1,ntInSq(i,2)+1) = nan;
                densityPlotAct(ntInSq(i,1)+1,ntInSq(i,2)+1) = nan;
            end
        end
        densityPlot(:,:,iSet,iterI) = densityPlotTmp;
        % smooth, compute normalized activation maps
        densityPlotSm         = nanconv(densityPlotTmp,imageFilter, 'nanout'); % smoothing (that deals with nans properly 8/8/18)
        densityPlotActNormTmp = densityPlotAct./densityPlotActUpd; %divide by number of times that location was visited
        densityPlotActNormTmp = nanconv(densityPlotActNormTmp,imageFilter, 'nanout');
        densityPlotActNorm(:,:,iSet,iterI) = densityPlotActNormTmp;
        
        %compute the sum of the distances between each cluster and itself over batches 
        if iSet>1 
           clusDistB(iSet-1,iterI)=sum(sqrt(sum([(mu(:,1,round(toTrlN(iSet)./batchSize)))-(mu(:,1,round(toTrlN(iSet-1)./batchSize))),(mu(:,2,round(toTrlN(iSet)./batchSize)))-(mu(:,2,round(toTrlN(iSet-1)./batchSize)))].^2,2)));
        end
        
        if ~strcmp(dat(1:3),'cat') %if finding cats, won't be gridlike
            %compute autocorrmap
            aCorrMap = ndautoCORR(densityPlotSm);
            %compute gridness
            [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
            gA(iSet,iterI,:,1) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
            [g,gdataA] = gridSCORE(aCorrMap,'wills',0);
            gW(iSet,iterI,:,1) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];

            %compute gridness over normalised activation map - normalised by times loc visited
                aCorrMap = ndautoCORR(densityPlotActNormTmp);
                [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
                gA_actNorm(iSet,iterI,:,1) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
                
                [g,gdataA] = gridSCORE(aCorrMap,'wills',0);
                gW_actNorm(iSet,iterI,:,1) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
            
            %trapz left vs right - split in half then compute gridness for each half
            if  strcmp(dat(1:4),'trap')
                %left half of box
                aCorrMap = ndautoCORR(densityPlotSm(:,1:hLeft));
                [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
                gA(iSet,iterI,:,2) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
                [g,gdataA] = gridSCORE(aCorrMap,'wills',0);
                gW(iSet,iterI,:,2) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
                
                %right half of box
                aCorrMap = ndautoCORR(densityPlotSm(:,h-hRight:end));
                [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
                gA(iSet,iterI,:,3) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
                [g,gdataA] = gridSCORE(aCorrMap,'wills',0);
                gW(iSet,iterI,:,3) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
            end
        end
    end
    
end
end