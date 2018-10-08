function [permPrc_gA_act, permPrc_gW_act,permPrc_gA_actNorm, permPrc_gW_actNorm, gA_act,gA_actNorm,gW_act,gW_actNorm,gA_actNormPerm, gW_actNormPerm, densityPlotAct,densityPlotActNorm, rSeedTest] = gridnessTestData_Perm(densityPlot,dat,locRange,nClus,nTrialsTest,nPerm,nIters2run,doPerm)
% computes activation maps after learning, also can perform shuffling for
% permtation tests for stats for grid scores

spacing=linspace(locRange(1),locRange(2),locRange(2)+1);
useSameTrls = 0;
jointTrls   = 1;
stepSize = 1;
sigmaGauss = stepSize;

gaussSmooth = 1;
imageFilter=fspecial('gaussian',5,gaussSmooth); %this is default for imgaussfilt

b = length(spacing);
h = length(spacing);%for trapz
%trapzKrupic
hLeft=17;% 12;
hRight=33;% - 33 = start from 18 from left % 27;

%perm testing
permPrc_gA_act = nan(nIters2run,4);
permPrc_gW_act = nan(nIters2run,4);
permPrc_gA_actNorm = nan(nIters2run,4);
permPrc_gW_actNorm = nan(nIters2run,4);

%% compute actNorm maps after 'training' (on a new test set of locations)
[trialsTest,~, rSeedTest,ntInSq] = createTrls(dat,nTrialsTest,locRange,useSameTrls,jointTrls,h);

densityPlotAct     = zeros(b,h,nIters2run);
densityPlotActNorm = zeros(b,h,nIters2run);
gA_act = nan(nIters2run,9);
gW_act = nan(nIters2run,9);
gA_actNorm = nan(nIters2run,9);
gW_actNorm = nan(nIters2run,9);
if strcmp(dat(1:4),'trap') %if trapz - compute gridness of left/right half of boxes too
    gA_act = nan(nIters2run,9,3);
    gW_act = nan(nIters2run,9,3);
    gA_actNorm = nan(nIters2run,9,3);
    gW_actNorm = nan(nIters2run,9,3);
end
muTrain = nan(nClus,2,nIters2run);

gA_actPerm = nan(nPerm,nIters2run);
gW_actPerm = nan(nPerm,nIters2run);
gA_actNormPerm = nan(nPerm,nIters2run);
gW_actNormPerm = nan(nPerm,nIters2run);
for iterI=1:nIters2run
    fprintf('Running iter %d\n',iterI);
        actTrl = zeros(nClus,nTrialsTest);
        densityPlotActTmp = zeros(b,h);
        densityPlotActUpd = zeros(b,h);
        indTmp=find(densityPlot(:,:,end,iterI)>0); %sometimes clus not in shape, not counted, so fewer clusters than expected
        [muTrain(1:length(indTmp),1,iterI), muTrain(1:length(indTmp),2,iterI)] = find(densityPlot(:,:,end,iterI)>0); %find cluster positions - added >0 since counts nans as well
        
        dist2Clus = sqrt(sum(reshape([muTrain(:,1,iterI)'-trialsTest(:,1), muTrain(:,2,iterI)'-trialsTest(:,2)].^2,nTrialsTest,nClus,2),3)); %compute distances
        closestC = nan(1,nTrialsTest);
        for iTrl = 1:nTrialsTest
            closestTmp = find(min(dist2Clus(iTrl,:))==dist2Clus(iTrl,:));
            if numel(closestTmp)>1
                closestC(iTrl) = randsample(closestTmp,1);
            else
                closestC(iTrl) = closestTmp;
            end
            %compute activation
            actTrl(closestC(iTrl),iTrl)=mvnpdf(trialsTest(iTrl,:),muTrain(closestC(iTrl),:,iterI),eye(2)*sigmaGauss); % save only the winner
        end
        %densityPlotActNorm
        for iTrl = 1:nTrialsTest
            densityPlotActTmp(trialsTest(iTrl,1)+1, trialsTest(iTrl,2)+1)    = densityPlotActTmp(trialsTest(iTrl,1)+1, trialsTest(iTrl,2)+1)+ nansum(actTrl(:,iTrl));
            densityPlotActUpd(trialsTest(iTrl,1)+1, trialsTest(iTrl,2)+1)       = densityPlotActUpd(trialsTest(iTrl,1)+1, trialsTest(iTrl,2)+1)+1; %log nTimes loc was visited
        end

        %turn 0s outside of the shape into nans
        if ~strcmp(dat(1:2),'sq') % if circ/trapz
            for i=1:length(ntInSq)
                densityPlotActTmp(ntInSq(i,1)+1,ntInSq(i,2)+1) = nan;
            end
        end
        densityPlotActNormTmp = densityPlotActTmp./densityPlotActUpd; %divide by number of times that location was visited

        %smooth - not necessary but do anyway to compare with perm (which needs smoothing)
        densityPlotAct(:,:,iterI)               =  nanconv(densityPlotActTmp,imageFilter, 'nanout');
        densityPlotActNorm(:,:,iterI)           =  nanconv(densityPlotActNormTmp,imageFilter, 'nanout');        
        
        aCorrMap = ndautoCORR(densityPlotActTmp);
        [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
        gA_act(iterI,:,1) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
        
        [g,gdataA] = gridSCORE(aCorrMap,'wills',0);
        gW_act(iterI,:,1) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
        
        %compute gridness over normalised activation map - normalised by times loc visited
        aCorrMap = ndautoCORR(densityPlotActNormTmp);
        [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
        gA_actNorm(iterI,:,1) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
        [g,gdataA] = gridSCORE(aCorrMap,'wills',0);
        gW_actNorm(iterI,:,1) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
        
        if  strcmp(dat(1:4),'trap')
            %left
            aCorrMap = ndautoCORR(densityPlotActTmp(:,1:hLeft));
            [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
            gA_act(iterI,:,2) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
            [g,gdataA] = gridSCORE(aCorrMap,'wills',0);
            gW_act(iterI,:,2) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
            %right
            aCorrMap = ndautoCORR(densityPlotActTmp(:,h-hRight+1:end));
            [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
            gA_act(iterI,:,3) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
            [g,gdataA] = gridSCORE(aCorrMap,'wills',0);
            gW_act(iterI,:,3) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
            %left
            aCorrMap = ndautoCORR(densityPlotActNormTmp(:,1:hLeft));
            [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
            gA_actNorm(iterI,:,2) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
            [g,gdataA] = gridSCORE(aCorrMap,'wills',0);
            gW_actNorm(iterI,:,2) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
            %right
            aCorrMap = ndautoCORR(densityPlotActNormTmp(:,h-hRight+1:end));
            [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
            gA_actNorm(iterI,:,3) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
            [g,gdataA] = gridSCORE(aCorrMap,'wills',0);
            gW_actNorm(iterI,:,3) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
        end

        %permutation testing
        minTime = 20; % shuffle setting - each activation at least 20 time points away from orig
        if doPerm
            for iPerm = 1:nPerm                
                %shuffle trials, and log trials < 20 timepoints to orig
                ind2short = zeros(1,nTrialsTest);
                start=0; %get the while loop running
                while ~start || length(find(ind2short))==1 %if only 1, can't shuffle
                    randVec = randperm(nTrialsTest); %shuffle
                    for iTrl = 1:nTrialsTest
                        if (abs(randVec(iTrl)-iTrl)<minTime)
                            ind2short(iTrl) = iTrl; %log shuffled trials that are <20 timepoints to orig
                        end
                    end
                    start=1;
                end
                % shuffle  timepoints logged <20s, apply shuffled timepoints, double check all <20s
                ind2shuf = find(ind2short);
                start2=0;
                while ~start2 || any(ind2short2) %if any still < min time, perm again
                    permInd  = randperm(length(ind2shuf));
                    while any(permInd==1:length(permInd)) %if any stay in the same place, shuffle again
                        permInd  = randperm(length(ind2shuf));
                    end
                    indShuf = ind2shuf(permInd);
                    randVec(ind2shuf) = randVec(indShuf);
                    %double check if any shuffled timings too short
                    ind2short2 = zeros(1,length(ind2shuf));
                    for i = ind2shuf
                        if (abs(randVec(ind2shuf)-ind2shuf)<minTime)
                            ind2short2(i) = i;
                            warning('perm shuffling stuck on < min time?');
                        end
                    end
                    start2=1;
                end
                actAllPerm = actTrl(:,randVec);
                
                densityPlotActPerm       = zeros(b,h);
                densityPlotActUpdPerm    = zeros(b,h);
                for iTrl = 1:nTrialsTest
                    densityPlotActPerm(trialsTest(iTrl,1)+1, trialsTest(iTrl,2)+1)    = densityPlotActPerm(trialsTest(iTrl,1)+1, trialsTest(iTrl,2)+1)+ sum(actAllPerm(:,iTrl));
                    densityPlotActUpdPerm(trialsTest(iTrl,1)+1, trialsTest(iTrl,2)+1) = densityPlotActUpdPerm(trialsTest(iTrl,1)+1, trialsTest(iTrl,2)+1)+1; %log nTimes loc was visited
                end
                densityPlotActNormPerm   = densityPlotActPerm./densityPlotActUpdPerm; %divide by number of times that location was visited
                %smooth
                densityPlotActPerm     = nanconv(densityPlotActPerm,imageFilter, 'nanout'); 
                densityPlotActNormPerm = nanconv(densityPlotActNormPerm,imageFilter, 'nanout');

                %compute gridness
                aCorrMap = ndautoCORR(densityPlotActPerm);
                [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
                gA_actPerm(iPerm,iterI) = gdataA.g_score;
                [g,gdataW] = gridSCORE(aCorrMap,'wills',0);
                gW_actPerm(iPerm,iterI) = gdataW.g_score;
                aCorrMap = ndautoCORR(densityPlotActNormPerm);
                [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
                gA_actNormPerm(iPerm,iterI) = gdataA.g_score;
                [g,gdataW] = gridSCORE(aCorrMap,'wills',0);
                gW_actNormPerm(iPerm,iterI) = gdataW.g_score;
            end
            permPrc_gA_act(iterI,:) = prctile(gA_actPerm(:,iterI),[2.5, 5, 95, 97.5]);
            permPrc_gW_act(iterI,:) = prctile(gW_actPerm(:,iterI),[2.5, 5, 95, 97.5]);
            permPrc_gA_actNorm(iterI,:) = prctile(gA_actNormPerm(:,iterI),[2.5, 5, 95, 97.5]);
            permPrc_gW_actNorm(iterI,:) = prctile(gW_actNormPerm(:,iterI),[2.5, 5, 95, 97.5]);
        end
end
end