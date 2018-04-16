% function [densityPlot,densityPlotAct,densityPlotActNorm,clusMu,gA,gW,gA_act,gW_act,gA_actNorm,gW_actNorm,rSeed,muAll] = covering_map_batch_sim(nClus,locRange,warpType,epsMuOrig,nTrials,batchSize,nIter,warpBox,alpha,trials,useSameTrls,trialsUnique,stochasticType,c,dat,weightEpsSSE)
function [densityPlot,densityPlotActNorm,gA,gA_actNorm,muInit,rSeed,clusDistB,permPrc,muAll, trials] = covering_map_batch_sim(nClus,locRange,warpType,epsMuOrig,nTrials,batchSize,nIter,warpBox,trials,useSameTrls,stochasticType,c,dat,boxSize,annEps,jointTrls,doPerm)

%if end up not using desityPlotAct and gA_act - edit out below, or no need
%to save the iters, etc.


spacing=linspace(locRange(1),locRange(2),locRange(2)+1); 
stepSize=diff(spacing(1:2)); 
% nSteps = length(spacing);
% nTrialsTest = nTrials;

gaussSmooth=1; %smoothing for density map

% nBatch=floor(nTrials/batchSize);
nBatch=round(nTrials/batchSize); %nBatch 7500 this better?
batchSize = floor(batchSize); % when have decimal points, above needed

%if decrease learning rate over time: 1/(1+decay+timestep); Decay - set a param
if annEps
    nBatches = nTrials./batchSize; %2500, 5000, 250000, 50000
    epsMuOrig = nBatches/100;
    annEpsDecay = nBatches/40;
end
% trlSel = ceil([nBatch*.25, nBatch*.5, nBatch*.67, nBatch*.75, nBatch*.9, nBatch+1]);

%also get a bunch of trials to plot activations from current trial (gauss
%func of clus loc) - make it 6 like sets above
%note that this is averaging over trials, not just batches though - might
%be to show activations as clusters are stationary as well as move?

%compute gridness over time
% fromTrlI  = round([1,                  nTrials.*.1+1,    nTrials.*.25+1,   nTrials.*.30+1,  nTrials.*.40+1,  nTrials.*.50+1, nTrials.*.75+1, nTrials.*.90+1]);
% toTrlN    = round([nTrials.*.01+1,  nTrials.*.11+1,  nTrials.*.26+1,  nTrials.*.31+1, nTrials.*.41+1, nTrials.*.51+1, nTrials.*.76+1, nTrials.*.91+1]);%
% toTrlN    = round([nTrials.*.025+1,  nTrials.*.125+1,  nTrials.*.275+1,  nTrials.*.325+1, nTrials.*.425+1, nTrials.*.525+1, nTrials.*.775+1, nTrials.*.925+1]);%
% toTrlN    = round([nTrials.*.045+1,    nTrials.*.145+1,  nTrials.*.295+1,  nTrials.*.345+1, nTrials.*.445+1, nTrials.*.545+1, nTrials.*.795+1, nTrials.*.945+1]);%
% toTrlN    = round([nTrials.*.1+1,    nTrials.*.2+1,  nTrials.*.35+1,  nTrials.*.4+1, nTrials.*.5+1, nTrials.*.6+1, nTrials.*.85+1, nTrials]);%

%20 timepoints; last 2: first quarter, last quarter
%half, last quarter
fromTrlI  = round([1,               nTrials.*.05+1,  nTrials.*.1+1,  nTrials.*.15+1,  nTrials.*.2+1,  nTrials.*.25+1,   nTrials.*.3+1,  nTrials.*.35+1,  nTrials.*.4+1,  nTrials.*.45+1,  nTrials.*.5+1,  nTrials.*.55+1, nTrials.*.6+1,  nTrials.*.65+1, nTrials.*.7+1, nTrials.*.75+1, nTrials.*.8+1,  nTrials.*.85+1,  nTrials.*.9+1,  nTrials.*.95+1, nTrials.*.5, nTrials.*.75]);
toTrlN    = round([nTrials.*.05,    nTrials.*.1,     nTrials.*.15,   nTrials.*.2,     nTrials.*.25,   nTrials.*.3,      nTrials.*.35,   nTrials.*.4,     nTrials.*.45,   nTrials.*.5,     nTrials.*.55,   nTrials.*.6,    nTrials.*.65,   nTrials.*.7,    nTrials.*.75,  nTrials.*.8,    nTrials.*.85,   nTrials.*.9,     nTrials.*.95,   nTrials,        nTrials,      nTrials]);

if nargout > 8
    muAll            = nan(nClus,2,nBatch+1,nIter);
end
nSets              = length(fromTrlI);
% actAll               = nan(nClus,nTrials); %keep this trial by trial %no need to declare if not saving over iters
muInit               = nan(nClus,2,nIter);

% tsseAll              = nan(nBatch+1,nIter); %not sure if this is +1 or not; not tested
% rSeed = struct(1,nIter); %how to initialise this struct?
gA = nan(nSets,nIter,9); %if saving the 5 r values, 9. if not, 4.
% gW = nan(nSets,nIter,9);
% gA_act = nan(nSets,nIter,9);
% gW_act = nan(nSets,nIter,9);
gA_actNorm = nan(nSets,nIter,9);
% gW_actNorm = nan(nSets,nIter,9);

%if trapz - compute gridness of left/right half of boxes too
if strcmp(dat(1:4),'trap')
    gA = nan(nSets,nIter,9,3);
    gA_actNorm = nan(nSets,nIter,9,3);
end


%perm testing
nPerm = 1000;
permPrc = nan(nIter,4);
% gA_actNormPerm = nan(nPerm,nIter,9,1);

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
        if boxSize==2
            spacingTrapz = spacing(14:37+length(14:37));
        end
        if length(dat)==11 % 'trapzKrupic', no number following
            a = length(spacingTrapz);
            b = length(spacing);
            h = length(spacing);
        else
            if strcmp(dat(12),'2') % 'trapzKrupic2' - if do larger, add here
                b=length(spacing)*1.5; %make smaller; since datapoints dont reach out there
                h=length(spacing)*2;
            end
        end
    end
    
    %split left and right half the trapz with equal areas
    % c is the length of the line when split the trapz in half

    %krupic should work, check if valid for otehrs - esp rounding below
    
    halfArea = (((a+b)/2)*h)/2;
    c = sqrt(((a^2)+(b^2))/2);
    %(b+c)/2)*h=c
%     hLeft  = round(halfArea/((b+c)/2)); %bigger side

% new to equalize area - Krupic
    hLeft  = floor(halfArea/((b+c)/2)); %bigger side
    hRight = ceil(halfArea/((a+c)/2))+1; %smaller side

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

for iterI = 1:nIter
%     epsMu = epsMuOrig; %revert learning rate back to original if reset
    
    fprintf('iter %d \n',iterI);

    [trials,dataPtsTest, rSeed(iterI)] = createTrls(dat,nTrials,locRange,useSameTrls,jointTrls,boxSize,h);
        
    % if expand box
    switch warpType
        case 'sq2rect'
            trialsExpand = [randsample(locRange(1):diff(spacing(1:2)):locRange(2)*2,nTrials*.75,'true'); randsample(spacing,nTrials*.75,'true')]';
        case 'rect2sq'
            trialsExpand = [randsample(spacing,nTrials*.75,'true'); randsample(spacing,nTrials*.75,'true')]'; %this doesn't work - actually, maybe this isn't a thing?
    end
    
    %initialise each cluster location  
    mu = nan(nClus,2,nBatch+1);
%     mu(:,:,1) = kmplusInit(dataPtsTest,nClus); %kmeans++ initialisation
    mu(:,:,1) = dataPtsTest(randi(nTrials,nClus,1),:); %forgy method
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

            %compute distances
%             dist2Clus = sqrt(sum([mu(:,1,iTrl)'-trials(iTrl,1); mu(:,2,iTrl)'-trials(iTrl,2)].^2)); % vectorising euclid dist - sqrt(sum((a-b).^2)), since can't use xval method
            
            %compute distances - vectorise both clusters and trials (in batch)
            dist2Clus = sqrt(sum(reshape([mu(:,1,iBatch)'-trls2Upd(:,1), mu(:,2,iBatch)'-trls2Upd(:,2)].^2,batchSize,nClus,2),3));% reshapes it into batchSize x nClus x 2 (x and y locs)
%             dist2Clus = sqrt(sum(dist2Clus,3));
            
            %stochastic update - select 1 of the closest clusters w random element - stochastic parameter c - large is deterministic, 0 - random            
%             if stochasticType %stochastic update
%                 if stochasticType==1
% %                     beta=c*(iTrl-1);         % so this gets bigger, and more deterministic with more trials
%                     beta=c*(iTrl+nTrials/50);  % maybe no need to be so stochatic at start
%                 elseif stochasticType==2
%                     beta=c*(iTrl+nTrials/50);  % maybe no need to be so stochatic at start
%                     if beta >= c*500 %it might be worth checking if this depends on nClusters - distances will change
%                         beta = c*500;
%                     end
%                 elseif stochasticType==3
%                     beta = c*500; % as a function of c %prev:.185; 
%                 end
%                 dist2Clus2 = exp(-beta.*dist2Clus)./ sum(exp(-beta.*dist2Clus));
%                 distClusPr = cumsum(dist2Clus2);
%                 closestC=find(rand(1)<distClusPr,1);
%                 if isempty(closestC)  %if beta is too big, just do deterministic update (already will be near deterministic anyway)
%                     closestC=find(min(dist2Clus)==dist2Clus);
%                 end
%                 cParams.betaAll(iTrl)       = beta;
%             else %deterministic update
%                 closestC=find(min(dist2Clus)==dist2Clus);
%             end
%             if numel(closestC)>1 %if more than 1, randomly choose one
%                 closestC = randsample(closestC,1);
%             end
            
            %deterministic update
%             closestC=find(min(dist2Clus)==dist2Clus);

            %deterministic update - batch; save closestC for each trial
            %%%% -need to find min dist cluster for each trial; better way?
            %%%% / vectorize?
            closestC = nan(1,batchSize);
            for iTrlBatch = 1:batchSize
                closestTmp = find(min(dist2Clus(iTrlBatch,:))==dist2Clus(iTrlBatch,:));
                if numel(closestTmp)>1
                    closestC(iTrlBatch) = randsample(closestTmp,1);
                elseif numel(closestTmp)<1
                    a=1;
                else
                    closestC(iTrlBatch) = closestTmp;
                    
                end
                
                %compute activation
                sigmaGauss = stepSize;%/3.5; %move up later - 1 seems to be fine
                actTrl(closestC(iTrlBatch),iTrlBatch)=mvnpdf(trls2Upd(iTrlBatch,:),mu(closestC(iTrlBatch),:,iBatch),eye(2)*sigmaGauss); % save only the winner
                
                
                
            end
            actTrlAll(:,:,iBatch) = actTrl;

%             %if stochastic, closestC might not be if more than 1 actual
%             %closest, have to check if stoch update chose one of them
%             actualClosestC = find(min(dist2Clus)==dist2Clus); 
%             closestMatch = zeros(1,length(actualClosestC));
%             for iC = 1:length(actualClosestC)
%                 closestMatch(iC) = actualClosestC(iC) == closestC;
%             end
%             cParams.closestChosen(iTrl) = any(closestMatch); %was (one of) the closest cluster selected?
%             cParams.closestDist(iTrl)   = dist2Clus(closestC);

            %log which cluster has been updated
%             updatedC(iTrl) = closestC;
            
            %learning rate
            if annEps %if use annealed learning rate
                epsMu = epsMuOrig./(1+annEpsDecay+iBatch); %need to check if it's actually epsMuOrig*(1/(1+annEpsDecay+iBatch);
                %                     clear epsAll
                %                     for iBatch=1:nBatches, epsAll(iBatch)=epsMuOrig/(1+annEpsDecay+iBatch); end
                %                     figure; plot(epsAll);
                %                     epsAll(nBatches.*[.25, .5, .75]) % eps at 25%, 50%, 75% of trials: 0.0396    0.0199    0.0133
            else
                epsMu = epsMuOrig;
            end
            
            %batch update - save all updates for each cluster for X
            %trials, update it, then again
            % - goes through each cluster, compute distances, average,
            %then update
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
    if nargout > 8
        muAll(:,:,:,iterI)      = mu;
    end
    actAll  = reshape(actTrlAll,nClus,nTrials); %save trial-by-trial act over blocks, here unrolling it
    
%     tsseAll(:,iterI)        = tsse;
%     sseSpreadSd(:,iterI)    = stdAcrossClus;
%     sseSpreadVar(:,iterI)   = varAcrossClus;

    % densityPlotClus - density plot with each cluster in dim 3 - more like
    % a place cell map - leave it here so can use diff smoothing values outside . also then use to make
    % it a gridcell map: densityPlot=sum(densityPlotClus,3); - to compute autocorrelogram
    
    %densityplot over time (more samples)
    for iSet = 1:nSets
        densityPlotClus      = zeros(b,h,nClus);
        densityPlotAct       = zeros(b,h);
        densityPlotActUpd    = zeros(b,h); %start with ones, so won't divide by 0
      
        clus = round(mu(:,:,round(toTrlN(iSet)./batchSize))); %mu is in batches
        ind=clus<=0; clus(ind)=1; %indices <= 0 make to 1
        for iClus=1:nClus
            ntNanInd = squeeze(~isnan(clus(iClus,1,:)));
            clusTmp = []; %clear else dimensions change over clus/sets
            clusTmp(1,:) = squeeze(clus(iClus,1,ntNanInd)); %split into two to keep array dim constant - when only 1 location, the array flips.
            clusTmp(2,:) = squeeze(clus(iClus,2,ntNanInd));
            for iTrlUpd=1:size(clusTmp,2)
                densityPlotClus(clusTmp(1,iTrlUpd),clusTmp(2,iTrlUpd),iClus) = densityPlotClus(clusTmp(1,iTrlUpd),clusTmp(2,iTrlUpd),iClus)+1;
            end
        end
        densityPlot(:,:,iSet,iterI) = nansum(densityPlotClus,3);
        
        %densityPlotActNorm
        for iTrl = fromTrlI(iSet):toTrlN(iSet)
            densityPlotAct(trials(iTrl,1)+1, trials(iTrl,2)+1)    = densityPlotAct(trials(iTrl,1)+1, trials(iTrl,2)+1)+ nansum(actAll(:,iTrl));
            densityPlotActUpd(trials(iTrl,1)+1, trials(iTrl,2)+1) = densityPlotActUpd(trials(iTrl,1)+1, trials(iTrl,2)+1)+1; %log nTimes loc was visited
        end
        densityPlotActNorm(:,:,iSet,iterI) = densityPlotAct./densityPlotActUpd; %divide by number of times that location was visited
        densityPlotSm                      = imgaussfilt(densityPlot(:,:,iSet,iterI),gaussSmooth); %smooth
        densityPlotActNormSm               = imgaussfilt(densityPlotActNorm(:,:,iSet,iterI),gaussSmooth);   
        
        %compute the sum of the distances between each cluster and itself over batches 
        if iSet>1 
           clusDistB(iSet-1,iterI)=sum(sqrt(sum([(mu(:,1,round(toTrlN(iSet)./batchSize)))-(mu(:,1,round(toTrlN(iSet-1)./batchSize))),(mu(:,2,round(toTrlN(iSet)./batchSize)))-(mu(:,2,round(toTrlN(iSet-1)./batchSize)))].^2,2)));
        end
        
        if ~strcmp(dat,'cat') %if finding cats, won't be gridlike
            %compute autocorrmap
            aCorrMap = ndautoCORR(densityPlotSm);
            %compute gridness
            [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
            gA(iSet,iterI,:,1) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
            %compute gridness over activation map
            %             aCorrMap = ndautoCORR(densityPlotActTSm);
            %             [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
            %             gA_act_t(iSet,iterI,:,1) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
            
            %compute gridness over normalised activation map - normalised by times loc visited
            aCorrMap = ndautoCORR(densityPlotActNormSm);
            [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
            gA_actNorm(iSet,iterI,:,1) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
            
            %trapz left right of box
            %split in half then compute gridness for each half
            if  strcmp(dat(1:4),'trap')
                
                %left half of box
%                 aCorrMap = ndautoCORR(densityPlotSm(:,spacing(1):spacing(ceil(length(spacing)/2))));
                aCorrMap = ndautoCORR(densityPlotSm(:,1:hLeft));
                [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
                gA(iSet,iterI,:,2) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
                
                %right half of box
%                 aCorrMap = ndautoCORR(densityPlotSm(:,spacing(ceil(length(spacing)/2))+1:spacing(end)));
                aCorrMap = ndautoCORR(densityPlotSm(:,hRight:end));
                [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
                gA(iSet,iterI,:,3) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
                
                %act
                %left half of box
%                 aCorrMap = ndautoCORR(densityPlotActNormSm(:,spacing(1):spacing(ceil(length(spacing)/2))));
%                 aCorrMap = ndautoCORR(densityPlotActNormSm(:,1:hLeft));
%                 [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
%                 gA_act(iSet,iterI,:,2) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
                
%                 aCorrMap = ndautoCORR(densityPlotActNormSm(:,spacing(ceil(length(spacing)/2))+1:spacing(end)));
                aCorrMap = ndautoCORR(densityPlotActNormSm(:,1:hLeft));
                [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
                gA_actNorm(iSet,iterI,:,2) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
                
                %right half of box
%                 aCorrMap = ndautoCORR(densityPlotActNormSm(:,ceil(length(spacing)/2)+1:end));
                aCorrMap = ndautoCORR(densityPlotActNormSm(:,hRight:end));
                [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
                gA_actNorm(iSet,iterI,:,3) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
                
            end
            
        end
        
    end
    
    
    % for perm test, need actAll (nClus x nTrials), and mu (nClus x 2 x
    % nBatch). 
    % - note that since size of mu is nBatch, not nTrials, need to reformat
    % this so the it is nTrials size (so each mu cluster value is repeated
    % batchSize times).
    % - need to permulate the activations, then map them onto mu
    % (positions) or vice versal; make 1k/5k/10k
    % - then make the firing maps (with fromTrlI toTrlN) - assess gridness;
    % 95% percentile; i suppose only from positive side?
    
    %for now: shuffle is fine since not temporal dependencies between the
    %trials. 
    % LATER, with connected trials: need to cut "20s" periods, maybe longer than 20 timesteps
    % here (even 100/200?) and permute those
    
    if doPerm
    % Reformat mu (batchSize long) so it is nTrials long (so
    % each mu cluster value is repeated batchSize times) 
    muTrls = nan(nClus,2,nTrials);
        for iBatch = 1:nBatch
            muTrls(:,1,batchSize*(iBatch-1)+1:batchSize*(iBatch-1)+batchSize)=repmat(mu(:,1,iBatch),1,batchSize);
            muTrls(:,2,batchSize*(iBatch-1)+1:batchSize*(iBatch-1)+batchSize)=repmat(mu(:,2,iBatch),1,batchSize);
        end
        %perm testing
        gA_actNormPerm = nan(1,nPerm);
        
        for iPerm = 1:nPerm
            fprintf('Perm %d\n',iPerm);
            randInd=randperm(length(muTrls));
            %         muTrlsPerm = muTrls(:,:,randInd); % check
            actAllPerm = actAll(:,randInd); % check
            
            iSet = 20; %just final set enough?
            
            densityPlotActPerm       = zeros(b,h);
            densityPlotActUpdPerm    = zeros(b,h);
            % just permuted activations make more sense
            for iTrl = fromTrlI(iSet):toTrlN(iSet)
                densityPlotActPerm(trials(iTrl,1)+1, trials(iTrl,2)+1) = densityPlotActPerm(trials(iTrl,1)+1, trials(iTrl,2)+1)+ sum(actAllPerm(:,iTrl)); %.^2 to make it look better?
                densityPlotActUpdPerm(trials(iTrl,1)+1, trials(iTrl,2)+1) = densityPlotActUpdPerm(trials(iTrl,1)+1, trials(iTrl,2)+1)+1; %log nTimes loc was visited
            end
            densityPlotActNormPerm = densityPlotActPerm./densityPlotActUpdPerm; %divide by number of times that location was visited
            % smooth
            %         densityPlotSmPerm = imgaussfilt(densityPlotPerm(:,:,iSet,iterI),gaussSmooth);
            %         densityPlotActTSm     = imgaussfilt(densityPlotAct(:,:,iSet,iterI),gaussSmooth);
            densityPlotActNormSmPerm = imgaussfilt(densityPlotActNormPerm,gaussSmooth);
            
            %compute gridness
            aCorrMap = ndautoCORR(densityPlotActNormSmPerm);
            [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
            gA_actNormPerm(iPerm) = gdataA.g_score;
        end
        permPrc(iterI,:) = prctile(gA_actNormPerm,[2.5, 5, 95, 97.5]);
    end
end
end