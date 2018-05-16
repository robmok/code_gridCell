function [densityPlot,densityPlotActNorm,gA,gW,gA_actNorm,gW_actNorm,muInit,rSeed,clusDistB,muAll, trials] = covering_map_batch_sim(nClus,locRange,catsInfo,warpType,epsMuOrig,nTrials,batchSize,nIter,warpBox,trials,useSameTrls,stoch,c,dat,boxSize,annEps,jointTrls,actOverTime)

spacing=linspace(locRange(1),locRange(2),locRange(2)+1); 
stepSize=diff(spacing(1:2)); 
% nSteps = length(spacing);
% nTrialsTest = nTrials;

gaussSmooth=1; %smoothing for density map

nBatch=round(nTrials/batchSize); %nBatch 7500 this better?
batchSize = floor(batchSize); % when have decimal points, above needed

%if decrease learning rate over time: 1/(1+decay+timestep); Decay - set a param
if annEps
    nBatches = nTrials./batchSize; %2500, 5000, 250000, 50000
%     epsMuOrig = nBatches/100;
% %     epsMuOrig = nBatches/150; %200
%     annEpsDecay = nBatches/40;
    
    epsMuOrig = nBatches*.6;
    annEpsDecay = nBatches/10;
end

%compute gridness over time %20 timepoints; - 21 sets now - last one is last quarter
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
% actAll               = nan(nClus,nTrials); %keep this trial by trial %no need to declare if not saving over iters
muInit               = nan(nClus,2,nIter);

gA = nan(nSets,nIter,9); %if saving the 5 r values, 9. if not, 4.
gW = nan(nSets,nIter,9);
% gA_act = nan(nSets,nIter,9);
% gW_act = nan(nSets,nIter,9);
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

%set densityPlot array size
trapzSpacing{1} = spacing(10:41);
trapzSpacing{2} = spacing(7:44); 
trapzSpacing{3} = spacing(4:47);
if strcmp(dat(1:4),'trap') && length(dat)>10
    if strcmp(dat(1:11),'trapzScaled')
        spacingTrapz = trapzSpacing{str2double(dat(12))}; %trapzScaled1,2,3
        a=length(spacingTrapz); %trapz length1
        b=locRange(2)+1; %50 - trapz length2 - +1 to let density plot go from 1:50 (rather than 0:49)
        h=round(((locRange(2)+1)^2)/((a+b)/2))+1; %trapz height (area = 50^2; like square)
    elseif strcmp(dat(1:11),'trapzKrupic')
        spacingTrapz = spacing(14:37);
        if boxSize==1.5
            spacingTrapz = spacing(14:37+length(14:37)*.5);
        elseif boxSize==2
            spacingTrapz = spacing(14:37+length(14:37));
        end
        if length(dat)==11 % 'trapzKrupic'
            a = length(spacingTrapz);
            b = length(spacing);
            h = length(spacing);
        end
    end
    
    %split left and right half the trapz with equal areas
    % c is the length of the line when split the trapz in half
%     halfArea = (((a+b)/2)*h)/2;
%     c = sqrt(((a^2)+(b^2))/2);
%     %(b+c)/2)*h=c
%     hLeft  = floor(halfArea/((b+c)/2)); %bigger side
%     hRight = ceil(halfArea/((a+c)/2))+1; %smaller side

    % new to equalize area - Krupic
    %equal number of points in trapz (301 each; nPoints in trap 877/2=438.5)
    hLeft = 14;
    hRight = 20;    
else
    b=length(spacing);
    h=length(spacing);
    spacingTrapz = spacing;
end

% %now setting/refreshing some of these over sets below
% densityPlotClus    = zeros(b,h,nClus,nSets,nIter); 
densityPlot        = zeros(b,h,nSets,nIter);
% densityPlotAct     = zeros(b,h,nSets,nIter);
densityPlotActNorm = zeros(b,h,nSets,nIter);
closestChosen = nan(1,nTrials);
updatedC      = nan(1,nTrials);
            
for iterI = 1:nIter
%     epsMu = epsMuOrig; %revert learning rate back to original if reset
    
    fprintf('iter %d \n',iterI);

    [trials,dataPtsTest, rSeed(iterI),ntInSq] = createTrls(dat,nTrials,locRange,useSameTrls,jointTrls,boxSize,catsInfo);
        
    % if expand box
    switch warpType
        case 'sq2rect'
            trialsExpand = [randsample(locRange(1):diff(spacing(1:2)):locRange(2)*2,nTrials*.75,'true'); randsample(spacing,nTrials*.75,'true')]';
        case 'rect2sq'
            trialsExpand = [randsample(spacing,nTrials*.75,'true'); randsample(spacing,nTrials*.75,'true')]'; %this doesn't work - actually, maybe this isn't a thing?
    end
    
    %initialise each cluster location  
    mu = nan(nClus,2,nBatch+1); 
    
    if ~strcmp(dat(1:3),'cat')
        mu(:,:,1) = dataPtsTest(randi(nTrials,nClus,1),:); %forgy method
    else
        mu(:,:,1) = reshape(randsample(locRange(1):diff(spacing(1:2)):locRange(2),nClus*2,'true'),nClus,2); %random
    end
    % mu(:,:,1) = kmplusInit(dataPtsTest,nClus); %kmeans++ initialisation
    
    %add one where all clusters are randomly placed within the box (need to
    %figure out how to do this for irregular shapes)
    %++
   
    muInit(:,:,iterI) = mu(:,:,1);
    actTrl = zeros(nClus,batchSize);
    %%
    
%     updatedC = nan(nTrials,1);
    deltaMu  = zeros(nClus,2,nBatch);
%     clusUpdates = zeros(nClus,2); %acutally starting at 0 is OK, since there was no momentum from last trial
    
    actTrlAll = nan(nClus,batchSize,nBatch); %check

    
    for iBatch=1:nBatch
        batchInd=batchSize*(iBatch-1)+1:batchSize*(iBatch-1)+batchSize; %trials to average over
        trls2Upd = trials(batchInd,:); %trials to use this batch

%             %if change size of box half way
%             if iTrl == nTrials*.75 && warpBox
%                 trials(nTrials*.75+1:end,:) = trialsExpand;
%             end

            %compute distances - vectorise both clusters and trials (in batch)
            dist2Clus = sqrt(sum(reshape([mu(:,1,iBatch)'-trls2Upd(:,1), mu(:,2,iBatch)'-trls2Upd(:,2)].^2,batchSize,nClus,2),3));% reshapes it into batchSize x nClus x 2 (x and y locs)
            
            %update cluster positions
            closestC = nan(1,batchSize);
            for iTrlBatch = 1:batchSize
                if stoch %stochastic update
%                     beta=c*(iBatch-1);  % by batch  % so this gets bigger, and more deterministic with more trials
                    beta=c*((iBatch-1)*nTrials/1000);  % by batch, no need to be so stochatic at start
                    dist2Clus2 = exp(-beta.*dist2Clus(iTrlBatch,:))./ sum(exp(-beta.*dist2Clus(iTrlBatch,:)));
                    distClusPr = cumsum(dist2Clus2);
                    closestTmp=find(rand(1)<distClusPr,1);
                    if isempty(closestTmp)  %if beta is too big, just do deterministic update (already will be near deterministic anyway)
                        closestTmp = find(min(dist2Clus(iTrlBatch,:))==dist2Clus(iTrlBatch,:));
                    end
%                     betaAll((iBatch-1)*batchSize+iTrlBatch)     = beta; %to check
                else %deterministic update
                    closestTmp = find(min(dist2Clus(iTrlBatch,:))==dist2Clus(iTrlBatch,:));
                end
                
                if numel(closestTmp)>1 % if more than 1, randomly sel 1
                    closestC(iTrlBatch) = randsample(closestTmp,1);
                else
                    closestC(iTrlBatch) = closestTmp;
                end
                %compute activation
%                 if ~strcmp(dat(1:4),'trap') %not computing for trapz
                    sigmaGauss = stepSize;%/3.5; %move up later - 1 seems to be fine
                    actTrl(closestC(iTrlBatch),iTrlBatch)=mvnpdf(trls2Upd(iTrlBatch,:),mu(closestC(iTrlBatch),:,iBatch),eye(2)*sigmaGauss); % save only the winner
%                 end
                
%               %if stochastic, closestC might not be if more than 1 actual
%               %closest, have to check if stoch update chose one of them
                actualClosestC = find(min(dist2Clus(iTrlBatch,:))==dist2Clus(iTrlBatch,:));
                closestMatch = zeros(1,length(actualClosestC));
                for iC = 1:length(actualClosestC)
                    closestMatch(iC) = actualClosestC(iC) == closestC(iTrlBatch);
                end
                closestChosen((iBatch-1)*batchSize+iTrlBatch) = any(closestMatch); %was (one of) the closest cluster selected? plot to check
                %log which cluster has been updated
                updatedC((iBatch-1)*batchSize+iTrlBatch) = closestC(iTrlBatch);
            end
%             if ~strcmp(dat(1:4),'trap') %not computing for trapz
                actTrlAll(:,:,iBatch) = actTrl;
%             end
            
            %learning rate
            if annEps %if use annealed learning rate
                epsMu = epsMuOrig./(1+annEpsDecay+iBatch*250); %need to check if it's actually epsMuOrig*(1/(1+annEpsDecay+iBatch);
%                 clear epsAll
% %                 for iBatch=1:nBatches, epsAll(iBatch)=epsMuOrig/(1+annEpsDecay+iBatch); end
%                 for iBatch=1:nBatches, epsAll(iBatch)=epsMuOrig/(1+annEpsDecay+iBatch*250); end
%                 figure; plot(epsAll);
%                 epsAll(nBatches.*[.05, .25, .5, .75, .95]) % eps at 25%, 50%, 75% of trials: 0.0396    0.0199    0.0133
            else
                epsMu = epsMuOrig;
            end
            
            %batch update - save all updates for each cluster for X trials, update it, then again - goes through each cluster,
            %compute distances, average, then update
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
    if nargout > 9
        muAll(:,:,:,iterI)      = mu;
    end
%     if ~strcmp(dat(1:4),'trap') %not computing for trapz
        actAll  = reshape(actTrlAll,nClus,nTrials); %save trial-by-trial act over blocks, here unrolling it
%     end 

    % densityPlotClus - density plot with each cluster in dim 3 - more like
    % a place cell map - also then use to make it a gridcell map: densityPlot=sum(densityPlotClus,3); - to compute autocorrelogram
    
    %densityplot over time
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
            if any(clusTmp>49) || any(clusTmp<1) %if out of box, don't add to densityplot
                clusTmp=[];
            end
            for iTrlUpd=1:size(clusTmp,2)
                densityPlotClus(clusTmp(1,iTrlUpd)+1,clusTmp(2,iTrlUpd)+1,iClus) = densityPlotClus(clusTmp(1,iTrlUpd)+1,clusTmp(2,iTrlUpd)+1,iClus)+1;
            end
        end
        densityPlotTmp = nansum(densityPlotClus,3);

        %densityPlotActNorm
%         if ~strcmp(dat(1:4),'trap') %not computing for trapz
        for iTrl = fromTrlI(iSet):toTrlN(iSet)
            densityPlotAct(trials(iTrl,1)+1, trials(iTrl,2)+1)    = densityPlotAct(trials(iTrl,1)+1, trials(iTrl,2)+1)+ nansum(actAll(:,iTrl));
            densityPlotActUpd(trials(iTrl,1)+1, trials(iTrl,2)+1) = densityPlotActUpd(trials(iTrl,1)+1, trials(iTrl,2)+1)+1; %log nTimes loc was visited
        end
        %turn 0s outside of the shape into nans
        if ~strcmp(dat(1:2),'sq') % if circ/trapz
            for i=1:length(ntInSq)
                densityPlotTmp(ntInSq(i,1)+1,ntInSq(i,2)+1) = nan;
                densityPlotAct(ntInSq(i,1)+1,ntInSq(i,2)+1) = nan;
            end
        end
        densityPlot(:,:,iSet,iterI) = densityPlotTmp;
%         densityPlotTmp(densityPlotTmp==0) = nan; %for circ, and prob trapz, to ignore points not in the shape - do this above now; for actNorm already does it by dividing by 0
        densityPlotSm               = imgaussfilt(densityPlotTmp,gaussSmooth); %smooth
        densityPlotActNormTmp = densityPlotAct./densityPlotActUpd; %divide by number of times that location was visited

%%%%%%%%%%%%%%%
         % IF COMPUTING GRIDNESS FOR ACT ITSELF
%         densityPlotAct(densityPlotAct==0) =  nan; %for circ, and prob trapz, to ignore points not in the shape - if computing this at all (atm not)

          densityPlotActNorm(:,:,iSet,iterI) = densityPlotActNormTmp;
%         end
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
            %compute gridness over activation map
            %             aCorrMap = ndautoCORR(densityPlotActTSm);
            %             [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
            %             gA_act_t(iSet,iterI,:,1) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
            
            %compute gridness over normalised activation map - normalised by times loc visited
%             if ~strcmp(dat(1:4),'trap') %not computing for trapz
                aCorrMap = ndautoCORR(densityPlotActNormTmp);
                [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
                gA_actNorm(iSet,iterI,:,1) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
                
                [g,gdataA] = gridSCORE(aCorrMap,'wills',0);
                gW_actNorm(iSet,iterI,:,1) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
%             end
            
            %trapz left right of box
            %split in half then compute gridness for each half
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
                %                 %actNorm                 
%                 aCorrMap = ndautoCORR(densityPlotActNormSm(:,1:hLeft));
%                 [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
%                 gA_actNorm(iSet,iterI,:,2) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
%                 [g,gdataA] = gridSCORE(aCorrMap,'wills',0);
%                 gW_actNorm(iSet,iterI,:,2) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
%                                 
%                 %right half of box
%                 aCorrMap = ndautoCORR(densityPlotActNormSm(:,h-hRight:end));
%                 [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
%                 gA_actNorm(iSet,iterI,:,3) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
%                 [g,gdataA] = gridSCORE(aCorrMap,'wills',0);
%                 gW_actNorm(iSet,iterI,:,3) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
            end
        end
    end
    
end
end