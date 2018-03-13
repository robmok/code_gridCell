% function [actAll,densityPlot,densityPlotAct,clusMu,muAvg,nTrlsUpd,gA,gW,gA_act,gW_act,gA_actNorm,gW_actNorm,muAll] = gauss_2d_batch_sim(nClus,locRange,box,warpType,epsMuOrig,sigmaGauss,nTrials,batchSize,nIter,warpBox,alpha,trials,stochasticType,c,plotGrids,dat,weightEpsSSE)
% function [densityPlot,clusMu,gA,gW,rSeed,muAll] = gauss_2d_batch_sim(nClus,locRange,box,warpType,epsMuOrig,sigmaGauss,nTrials,batchSize,nIter,warpBox,alpha,trials,stochasticType,c,plotGrids,dat,weightEpsSSE)
function [actAll, densityPlot,clusMu,gA,gW,rSeed,muAll] = gauss_2d_batch_sim(nClus,locRange,box,warpType,epsMuOrig,sigmaGauss,nTrials,batchSize,nIter,warpBox,alpha,trials,stochasticType,c,plotGrids,dat,weightEpsSSE)

spacing=linspace(locRange(1),locRange(2),locRange(2)+1); 
stepSize=diff(spacing(1:2)); nSteps = length(spacing);

% epsMu = epsMuOrig;

%for crossvalidation - 
nTrialsTest = nTrials;% nTrialsTest = 5000; %keep it constant if want to compare with nTrials above

gaussSmooth=1; %smoothing for density map

% 2d gaussian gradient descent based update
% sigmaGauss = stepSize/3;
% sigmaGauss = stepSize/3.5;
% sigmaGauss = stepSize/4; % becomes square map

fitDeltaMuX=@(epsMuVec,x,y,fdbck,act)(epsMuVec.*(fdbck-act).*exp(-(x.^2+y.^2)./2.*sigmaGauss.^2).*(x./(2.*pi.*sigmaGauss.^4)));
fitDeltaMuY=@(epsMuVec,x,y,fdbck,act)(epsMuVec.*(fdbck-act).*exp(-(x.^2+y.^2)./2.*sigmaGauss.^2).*(y./(2.*pi.*sigmaGauss.^4)));

nBatch=floor(nTrials/batchSize);
batchSize = floor(batchSize); % when have decimal points, above needed

% selecting cluster postions from a certain phase to test gridness
% note: with batch, may be less purposeful to average over trials. it would
% also mean averaging over batches (which would be even more stable). only
% thing if just looking at some point in time, small batches = less stable.
trlSel = ceil([nBatch*.25, nBatch*.5, nBatch*.67, nBatch*.75, nBatch*.9, nBatch+1]);

if nargout > 4
    muAll            = nan(nClus,2,nBatch+1,nIter);
end
nSets                = length(trlSel);

clusMu               = nan(nClus,2,nSets,nIter);
muAvg                = nan(nClus,2,nSets,nIter);
actAll               = nan(nClus,nTrials,nIter); %keep this trial by trial
nTrlsUpd             = nan(nClus,nSets,nIter);
gA = nan(nSets,nIter,4);
gW = nan(nSets,nIter,4);
% gA_act = nan(nSets,nIter,4);
% gW_act = nan(nSets,nIter,4);
% gA_actNorm = nan(nSets,nIter,4);
% gW_actNorm = nan(nSets,nIter,4);

%if trapz - compute gridness of left/right half of boxes too
if strcmp(dat(1:4),'trap')
    gA = nan(nSets,nIter,4,3);
    gW = nan(nSets,nIter,4,3);
end

%set densityPlot array size
trapzSpacing{1} = spacing(10:41);
trapzSpacing{2} = spacing(7:44);
trapzSpacing{3} = spacing(4:47);
if strcmp(dat(1:4),'trap') && length(dat)>10
    if strcmp(dat(1:11),'trapzScaled')
        spacingTrapz = trapzSpacing{str2double(dat(12))}; %trapzScaled1,2,3
        a=length(spacingTrapz); %trapz length1
        b=locRange(2)+1; %50 - trapz length2
        h=round(((locRange(2)+1)^2)/((a+b)/2)); %trapz height (area = 50^2; like square)
    elseif strcmp(dat,'trapzKrupic2')
        b=length(spacing)*1.5; %make smaller; since datapoints dont reach out there
        h=length(spacing)*2;
    end
else
    b=length(spacing);
    h=length(spacing);
end
densityPlotClus      = zeros(b,h,nClus,nSets,nIter);
% densityPlotClusAct   = zeros(b,h,nClus,nSets,nIter);
densityPlot          = zeros(b,h,nSets,nIter);
% densityPlotAct       = zeros(b,h,nSets,nIter);

%     rSeed = 

for iterI = 1:nIter
    
    fprintf('iter %d \n',iterI);    
    %make seed so can regenerate the trials later without saving - load up
    %the seed using: rng(s); then run trials = ...
    rSeed(iterI)=rng;
    
    
    switch box
        case 'square'
            trials      = [randsample(linspace(locRange(1),locRange(2),locRange(2)+1),nTrialsTest,'true'); randsample(linspace(locRange(1),locRange(2),locRange(2)+1),nTrialsTest,'true')]'; % random points in a box
            dataPtsTest = [randsample(linspace(locRange(1),locRange(2),locRange(2)+1),nTrialsTest,'true'); randsample(linspace(locRange(1),locRange(2),locRange(2)+1),nTrialsTest,'true')]';
        case 'rect'
            trials      = [randsample(locRange(1):diff(spacing(1:2)):locRange(2)*2,nTrialsTest,'true'); randsample(spacing,nTrialsTest,'true')]';
            dataPtsTest = [randsample(locRange(1):diff(spacing(1:2)):locRange(2)*2,nTrialsTest,'true'); randsample(spacing,nTrialsTest,'true')]';
        case 'trapz' % not scale to Krupic: (this was trapzNorm)
            spacingTrapz = spacing;
            trapY=locRange(2).*trapmf(spacingTrapz,[spacingTrapz(1), spacingTrapz(round(length(spacingTrapz)*.25)), spacingTrapz(round(length(spacingTrapz)*.75)),spacingTrapz(end)]);
            trapX=spacingTrapz;
            trapPts=[];
            for i=1:length(trapY)
                trapPts = [trapPts, [repmat(trapX(i),1,length(0:stepSize:trapY(i))); 0:stepSize:trapY(i)]];
            end
            % trapPts(2,:)=trapPts(2,:).*2-1; %put it back into -1 to 1
            % use this to select from the PAIR in trapPts
            trialInd     = randi(length(trapPts),nTrials,1);
            trials       = trapPts(:,trialInd)';
            %dataPtsTest
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
            trialInd     = randi(length(trapPts),nTrials,1);
            trials       = trapPts(:,trialInd)';
            %dataPtsTest
            trialIndTest = randi(length(trapPts),nTrials,1);
            dataPtsTest  = trapPts(:,trialIndTest)';
        case 'circ'
            % Create logical image of a circle
            imageSizeX = nSteps;
            [columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeX);
            centerX = nSteps/2; centerY = nSteps/2;
            radius = nSteps/2-1;
            circIm = (rowsInImage - centerY).^2 ...
                + (columnsInImage - centerX).^2 <= radius.^2;
            circPts=[]; % find circle points in XY coords
            for iX=1:length(circIm)
                yVals = find(circIm(iX,:));
                circPts = [circPts; ones(length(yVals),1)*iX, yVals'];
            end
            trialInd=randi(length(circPts),nTrials,1);
            trials=circPts(trialInd,:);
            %dataPtsTest
            trialIndTest = randi(length(circPts),nTrials,1);
            dataPtsTest  = circPts(trialIndTest,:);
        case 'cat'
            % draw points from 2 categories (gaussian) from a 2D feature space
            nTrials = floor(nTrials/nCats); % points to sample
            for iCat = 1:nCats
                mu(iCat,:)=randsample(locRange(1)+10:locRange(2)-10,2,'true'); % ±10 so category centres are not on the edge
                datPtsGauss(:,:,iCat) = round(repmat(mu(iCat,:),nTrials,1) + randn(nTrials,2)*R); % key - these are the coordinates of the points
            end
            trials = reshape(datPtsGauss,nTrials,2);
            trials = trials(randperm(length(trials)),:);
            trialsUnique=[];
            
            for iCat = 1:nCats
                mu(iCat,:)=randsample(locRange(1)+10:locRange(2)-10,2,'true'); % ±10 so category centres are not on the edge
                datPtsGauss(:,:,iCat) = round(repmat(mu(iCat,:),nTrials,1) + randn(nTrials,2)*R); % key - these are the coordinates of the points
            end
            dataPtsTest = reshape(datPtsGauss,nTrials,2);
            dataPtsTest = dataPtsTest(randperm(length(dataPtsTest)),:);
            
    end
    
    % if expand box
    switch warpType
        case 'sq2rect'
            trialsExpand = [randsample(locRange(1):diff(spacing(1:2)):locRange(2)*2,nTrials*.75,'true'); randsample(spacing,nTrials*.75,'true')]';
        case 'rect2sq'
            trialsExpand = [randsample(spacing,nTrials*.75,'true'); randsample(spacing,nTrials*.75,'true')]'; %this doesn't work - actually, maybe this isn't a thing?
    end

    %initialise each cluster location  
    mu = nan(nClus,2,nBatch+1);
    mu(:,:,1) = kmplusInit(dataPtsTest,nClus); %kmeans++ initialisation
    
    actUpd = zeros(nClus,batchSize);
    muUpd = mu;  %variable that only logs updated clusters mean value % - get initialised locs

    %%

    updatedC     =   nan(nTrials,1);
    deltaMu      =   zeros(nClus,2,nTrials);  
    deltaMuBatch =   zeros(nClus,2,batchSize);  %new batch thing
    clusUpdates  =   zeros(nClus,2); %acutally starting at 0 is OK, since there was no momentum from last trial
    act          =   nan(nClus,batchSize);
    
    tsse            = nan(nTrials,1);
%     stdAcrossClus   = nan(nTrials,1);
%     varAcrossClus   = nan(nTrials,1);
    sseW            = ones(nTrials,1);
%     spreadW         = ones(nTrials,1);
%     spreadVarW      = ones(nTrials,1);
    
    for iBatch=1:nBatch
        batchInd=batchSize*(iBatch-1)+1:batchSize*(iBatch-1)+batchSize; %trials to average over

        trls2Upd = trials(batchInd,:); %trials to use this batch

        epsMuVec = zeros(nClus,1);

%             %if change size of box half way
%             if iTrl == nTrials*.75 && warpBox
%                 trials(nTrials*.75+1:end,:) = trialsExpand;
%             end
%             
            %compute distances
%             dist2Clus = sqrt(sum([mu(:,1,iTrl)'-trials(iTrl,1); mu(:,2,iTrl)'-trials(iTrl,2)].^2)); % vectorising euclid dist - sqrt(sum((a-b).^2)), since can't use xval method
            dist2Clus = sqrt(sum(reshape([mu(:,1,iBatch)'-trls2Upd(:,1), mu(:,2,iBatch)'-trls2Upd(:,2)].^2,batchSize,nClus,2),3));% reshapes it into batchSize x nClus x 2 (x and y locs)
            
%             %stochastic update - sel 1 of the closest clusters w random element - stochastic parameter c - large is deterministic, 0 - random            
% 
%             if stochasticType %stochastic update
%                 if stochasticType==1
% %                     beta=c*(iTrl-1);         % so this gets bigger, and more deterministic with more trials
%                     beta=c*(iTrl+nTrials/50);  % maybe no need to be so stochatic at start
%                 elseif stochasticType==2
% %                     beta=c*(iTrl-1);
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

            %deterministic update - batch; save closestC for each trial
            closestC = nan(1,batchSize);
            for iTrlBatch = 1:batchSize
                closestTmp = find(min(dist2Clus(iTrlBatch,:))==dist2Clus(iTrlBatch,:));
                if numel(closestTmp)>1
                    closestC(iTrlBatch) = randsample(closestTmp,1);
                else
                    closestC(iTrlBatch) = closestTmp;
                end
                %save the activation for each trial to update
                actUpd(closestC(iTrlBatch),iTrlBatch)=mvnpdf(trls2Upd(iTrlBatch,:),mu(closestC(iTrlBatch),:,iBatch),eye(2)*sigmaGauss); % save only the winner
                
                % hmm, better way?
%                 fbk=zeros(nClus,1); %put above later
%                 fbk(closestC(iTrlBatch))=1;
                epsMuVec(closestC(iTrlBatch))=epsMuOrig;
                
                % so this does it for each trial, set epsMuOrig, fbk at 1
                % - if want to use vector, need to note which clusters
                % might not have an update (unlikely but useful; esp for
                % cat learning) and put epsMuVec to 0 on those batches,
                % fdbck to 0 [though if have epsmuvec, no need bother?
                deltaMuBatch(closestC(iTrlBatch),1,iTrlBatch)=fitDeltaMuX(epsMuOrig,trls2Upd(iTrlBatch,1)-mu(closestC(iTrlBatch),1,iBatch),trls2Upd(iTrlBatch,2)-mu(closestC(iTrlBatch),2,iBatch),1,actUpd(closestC(iTrlBatch),iTrlBatch));
                deltaMuBatch(closestC(iTrlBatch),2,iTrlBatch)=fitDeltaMuY(epsMuOrig,trls2Upd(iTrlBatch,1)-mu(closestC(iTrlBatch),1,iBatch),trls2Upd(iTrlBatch,2)-mu(closestC(iTrlBatch),2,iBatch),1,actUpd(closestC(iTrlBatch),iTrlBatch));
                

            end
            deltaMuBatchAll(:,:,:,iBatch) = deltaMuBatch;
            actUpdAll(:,:,iBatch) = actUpd;
            
%             closestCbatch(:,iBatch)=closestC;



            
            %but need to also get deltaMuX and deltaMuY for each of these
            %trials; put this above? if so will be computing trial by trial
            %is there a way to vectorize the activations? (complicated bit
            %is that on each trial, it's a different one. 
            
            
            %one idea, save all the activations for each cluster when it's
            %closest, then run the deltaMuX and Y thing all in one go per
            %cluster? maybe faster
            
            %alt (but slow?) is to get the activation trial by trial as
            %above, then compute deltaX and Y as well
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            %%%
            %no need?
            %act - need save all activations over trials, but update only
            %over batch 
%             for iClus=1:nClus %vectorize over trials, since usually more trials in batch; faster than looping it above iTrlBatch loop -BUT no need to update every clus...
%                 act(iClus,:)=mvnpdf(trls2Upd,mu(iClus,:,iBatch),eye(2)*sigmaGauss);
%             end
            
%             actUpd(closestC,iTrl) =act(closestC,iTrl); % save only the winner

         
%             fbk=zeros(nClus,1); %put above later
%             fbk(closestC)=1;
%             % compute deltaMu first, then add the momentum value
%             deltaMuX=fitDeltaMuX(epsMuVec,trials(iTrl,1)-mu(:,1,iTrl),trials(iTrl,2)-mu(:,2,iTrl),fbk,act(:,iTrl));
%             deltaMuY=fitDeltaMuY(epsMuVec,trials(iTrl,1)-mu(:,1,iTrl),trials(iTrl,2)-mu(:,2,iTrl),fbk,act(:,iTrl));
%             
                

            
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

% %             epsMuVec(closestC)=epsMuOrig;
            
            %weight learning rate by sse
%             epsMuVec(closestC) = epsMuOrig*sseW(iTrl); %weight learning rate by prop SSE reduce from start
%             epsMuVec(closestC) = epsMuOrig*spreadW(iTrl); %weight learning rate by prop spreadout-ness reduced from start
%             epsMuVec(closestC) = epsMuOrig*spreadVarW(iTrl); %weight learning rate by prop spreadout-ness reduced from start
%             epsMuVec(closestC) = epsMuOrig*sseW(iTrl)*spreadW(iTrl); %weight learning rate by both the above
%             epsMuVec(closestC) = epsMuOrig*mean([sseW(iTrl),spreadW(iTrl)]); %weight learning rate by average of both the above
%             epsMuAll(iTrl)=epsMuVec(closestC);


            % update mean estimates
            mu(:,1,iBatch+1) = mu(:,1,iBatch) + mean(squeeze(deltaMuBatch(:,1,:)),2);
            mu(:,2,iBatch+1) = mu(:,2,iBatch) + mean(squeeze(deltaMuBatch(:,2,:)),2);
                                  
            %%%%
            
            
            
            
            
%             %update (with momemtum-like parameter)
%             deltaMu(:,1,iTrl) = ((1-alpha)*deltaMuX)+(alpha*clusUpdates(closestC,1));
%             deltaMu(:,2,iTrl) = ((1-alpha)*deltaMuY)+(alpha*clusUpdates(closestC,2));
%             
%             clusUpdates(closestC,1)=deltaMu(closestC,1,iTrl);
%             clusUpdates(closestC,2)=deltaMu(closestC,2,iTrl);
% 
%             deltaMuVec = zeros(nClus,2);
%             deltaMuVec(closestC,:) = deltaMu(closestC,:,iTrl); % only update winner            
%             
%             % update mean estimates
%             if iTrl~=nTrials %no need to update for last trial +1)
%                 mu(:,1,iTrl+1) = mu(:,1,iTrl) + deltaMuVec(:,1);
%                 mu(:,2,iTrl+1) = mu(:,2,iTrl) + deltaMuVec(:,2);
%                 %log updated value only
%                 % don't want to save the cluster positions when they don't move
%                 % but need to keep track of the last cluster position
%                 % - have two variables, one with all, one without
%                 % - the one without; this is with nans just remove
%                 % all the nans then put into densityPlotClus).
%                 muUpd(closestC,1,iTrl+1)=mu(closestC,1,iTrl+1);
%                 muUpd(closestC,2,iTrl+1)=mu(closestC,2,iTrl+1);
%             end
%             
            
            % compute sse on each trial with respect to 'all trials'
            % trials - since values are all points in the box, no need to use a
            % trialsTest, juse use all unique locations (unique pairs of xy) from trials
            
            if weightEpsSSE %need to edit if useing other measures (spreadW)
                
                sse=nan(1,nClus);
                distTrl=(mu(:,1,iTrl)'-trialsUnique(:,1)).^2+(mu(:,2,iTrl)'-trialsUnique(:,2)).^2; %vectorised
                
                [indValsTrl, indTmp]=min(distTrl,[],2); % find which clusters are points closest to
                for iClus = 1:size(clusMu,1)
                    sse(iClus)=sum(sum([mu(iClus,1,iTrl)-trialsUnique(indTmp==iClus,1), mu(iClus,2,iTrl)-trialsUnique(indTmp==iClus,2)].^2,2)); %distance from each cluster from training set to datapoints closest to that cluster
                    %                         sse(iClus)=sum(sum([clusMu(iClus,1,iSet,iterI)-dataPtsTest(indTmp==iClus,1), clusMu(iClus,2,iSet,iterI)-dataPtsTest(indTmp==iClus,2)].^2,2)); %distance from each cluster from training set to datapoints closest to that cluster
                end
                tsse(iTrl)=sum(sse);
                
                %             %compute 'spreaded-ness' - variance of SE across clusters is a measure
                %             %of this, assuming uniform data points
                %             devAvgSSE           = sse-mean(sse);
                %             stdAcrossClus(iTrl) = std(devAvgSSE); % may be better since normalises by nClus?
                %             varAcrossClus(iTrl) = var(devAvgSSE);
                sseW(iTrl+1) = tsse(iTrl)./tsse(1);% weight next learning rate by prop of sse from the start
                %                 if iTrl>100
                %                     sseW(iTrl+1) = tsse(iTrl)./mean(tsse(1:100));
                %                 end
            end
%             spreadW(iTrl+1) = stdAcrossClus(iTrl)./stdAcrossClus(1);
%             spreadW(iTrl+1) = (stdAcrossClus(iTrl)./stdAcrossClus(1)).^0.5; %0.75 %make it smaller
%             spreadVarW(iTrl+1) = varAcrossClus(iTrl)./varAcrossClus(1);
            
    end
    muAll(:,:,:,iterI) = mu;
    actAll(:,:,iterI)  = reshape(actUpdAll,nClus,nTrials); %save trial-by-trial act over blocks, here unrolling it    
    
%     closestClusAll(:,iterI) = reshape(closestCbatch,1,nTrials);
    
    for iSet = 1:nSets
        %compute density map
%         clus = round(muUpd(:,:,fromTrlI(iSet):toTrlN(iSet)));
        clus = round(mu(:,:,trlSel(iSet)));
%         actClus = actAll(:,fromTrlI(iSet)-1:toTrlN(iSet)-1,iterI); % NB: -1 trial - this is to match up with the 'updated' location on the next trial indexed by clus - if just use act, might want to match simply to mu itself
%         actClus = actAll(:,trlSel(iSet),iterI); 
        
        ind=clus<=0; clus(ind)=1;  % indices <= 0 make to 1
        ind=clus>50; clus(ind)=50; % if overshoot out of the box, make it 50
        for iClus=1:nClus
            ntNanInd = squeeze(~isnan(clus(iClus,1,:)));
            clusTmp = []; %clear else dimensions change over clus/sets
            clusTmp(1,:) = squeeze(clus(iClus,1,ntNanInd)); %split into two to keep array dim constant - when only 1 location, the array flips.
            clusTmp(2,:) = squeeze(clus(iClus,2,ntNanInd));
            nTrlsUpd(iClus,iSet,iterI) = nnz(ntNanInd);
%             actTmp   = squeeze(actClus(iClus,ntNanInd));
            for iTrlUpd=1:size(clusTmp,2)
                densityPlotClus(clusTmp(1,iTrlUpd),clusTmp(2,iTrlUpd),iClus,iSet,iterI)     = densityPlotClus(clusTmp(1,iTrlUpd),clusTmp(2,iTrlUpd),iClus,iSet, iterI)   + 1;
%                 densityPlotClusAct(clusTmp(1,iTrlUpd),clusTmp(2,iTrlUpd),iClus,iSet,iterI)  = densityPlotClusAct(clusTmp(1,iTrlUpd),clusTmp(2,iTrlUpd),iClus,iSet,iterI) + actTmp(iTrlUpd);
            end
        end
        
        %compute density map
        clus = round(mu(:,:,trlSel(iSet)));
        ind=clus<=0; clus(ind)=1; %indices <= 0 make to 1
        for iClus=1:nClus
            ntNanInd = squeeze(~isnan(clus(iClus,1,:)));
            clusTmp = []; %clear else dimensions change over clus/sets
            clusTmp(1,:) = squeeze(clus(iClus,1,ntNanInd)); %split into two to keep array dim constant - when only 1 location, the array flips.
            clusTmp(2,:) = squeeze(clus(iClus,2,ntNanInd));
            %             clusTmp = squeeze(clus(iClus,:,ntNanInd));
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
        end
        
        %need to compute clus mean of the activation map? might get sth
        %diff
        
        
        
        
        
        
        %make combined (grid cell) plot, smooth
        densityPlot(:,:,iSet,iterI) = sum(densityPlotClus(:,:,:,iSet,iterI),3); %save this
        densityPlotSm = imgaussfilt(densityPlot(:,:,iSet,iterI),gaussSmooth);
        
%         %make activation plot
%         densityPlotAct(:,:,iSet,iterI) = sum(densityPlotClusAct(:,:,:,iSet,iterI),3); % or average activation? normalise by number of updates? (indexed in densityPlot)
%         densityPlotActSm = imgaussfilt(densityPlotAct(:,:,iSet,iterI),gaussSmooth);
%         
%         %make normalised activation plot (normalised by number of updates
%         %at that location)
%         densityPlotActNormTmp=densityPlotAct(:,:,iSet,iterI)./densityPlot(:,:,iSet,iterI); %this creates nans at where desnityPlot is zero, so replace nans with zeros
%         ind=isnan(densityPlotActNormTmp); densityPlotActNormTmp(ind)=0;
%         densityPlotActNormSm = imgaussfilt(densityPlotActNormTmp,gaussSmooth);
%        
        if ~strcmp(dat,'cat') %if finding cats, won't be gridlike
            if plotGrids && nIter < 8
                %compute autocorrmap, no need to save
                aCorrMap = ndautoCORR(densityPlotSm);
                %compute gridness
                figure; hold on;
                subplot(3,3,1);
                imagesc(densityPlotSm);
                subplot(3,3,2);
                imagesc(aCorrMap);
                subplot(3,3,3);
                [g,gdataA] = gridSCORE(aCorrMap,'allen',1);
                [g,gdataW] = gridSCORE(aCorrMap,'wills',0);
                gA(iSet,iterI,:) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius];
                gW(iSet,iterI,:) = [gdataW.g_score, gdataW.orientation, gdataW.wavelength, gdataW.radius];
                
%                 %compute gridness for act
%                 aCorrMap = ndautoCORR(densityPlotActSm);
%                 subplot(3,3,4);
%                 imagesc(densityPlotActSm);
%                 subplot(3,3,5);
%                 imagesc(aCorrMap);
%                 subplot(3,3,6);
%                 [g,gdataA] = gridSCORE(aCorrMap,'allen',1);
%                 [g,gdataW] = gridSCORE(aCorrMap,'wills',0);
%                 gA_act(iSet,iterI,:) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius];
%                 gW_act(iSet,iterI,:) = [gdataW.g_score, gdataW.orientation, gdataW.wavelength, gdataW.radius];
%                 
%                 %compute gridness for normalised act
%                 aCorrMap = ndautoCORR(densityPlotActNormSm);
%                 subplot(3,3,7);
%                 imagesc(densityPlotActNormSm);
%                 subplot(3,3,8);
%                 imagesc(aCorrMap);
%                 subplot(3,3,9);
%                 [g,gdataA] = gridSCORE(aCorrMap,'allen',1);
%                 [g,gdataW] = gridSCORE(aCorrMap,'wills',0);
%                 gA_actNorm(iSet,iterI,:) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius];
%                 gW_actNorm(iSet,iterI,:) = [gdataW.g_score, gdataW.orientation, gdataW.wavelength, gdataW.radius];
            else
                aCorrMap = ndautoCORR(densityPlotSm);
                %compute gridness
                [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
                [g,gdataW] = gridSCORE(aCorrMap,'wills',0);
                gA(iSet,iterI,:) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius];
                gW(iSet,iterI,:) = [gdataW.g_score, gdataW.orientation, gdataW.wavelength, gdataW.radius];
                
%                 %compute gridness for act
%                 aCorrMap = ndautoCORR(densityPlotActSm);
%                 [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
%                 [g,gdataW] = gridSCORE(aCorrMap,'wills',0);
%                 gA_act(iSet,iterI,:) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius];
%                 gW_act(iSet,iterI,:) = [gdataW.g_score, gdataW.orientation, gdataW.wavelength, gdataW.radius];
%                 
%                 %compute gridness for normalised act
%                 aCorrMap = ndautoCORR(densityPlotActNormSm);
%                 [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
%                 [g,gdataW] = gridSCORE(aCorrMap,'wills',0);
%                 gA_actNorm(iSet,iterI,:) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius];
%                 gW_actNorm(iSet,iterI,:) = [gdataW.g_score, gdataW.orientation, gdataW.wavelength, gdataW.radius];
            end
        end
        
    end    
end
end

