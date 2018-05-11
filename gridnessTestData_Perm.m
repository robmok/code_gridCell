function [permPrc_gA_act, permPrc_gW_act,permPrc_gA_actNorm, permPrc_gW_actNorm, gA_act,gA_actNorm,gW_act,gW_actNorm,gA_actNormPerm, gW_actNormPerm, densityPlotAct,densityPlotActNorm, rSeedTest] = gridnessTestData_Perm(densityPlot,dat,locRange,nClus,nTrialsTest,nPerm,nIters2run,doPerm)


%input - 
% densityPlot(unsmoothed)
% dat - shape
% locRange
% nClus
% nTrialsTest - number of trials to test here

% output
% densityPlotAct,densityPlotActNorm,gA_act,gA_actNorm,gW_act,gW_actNorm, rSeedTest

%tmp - inputs
% locRange = [0 49];
% nTrialsTest = 100000; %?
% dat = 'trapzKrupic';
% nClus=20;

spacing=linspace(locRange(1),locRange(2),locRange(2)+1);
boxSize     = 1;
useSameTrls = 0;
jointTrls   = 1;
stepSize = 1;
sigmaGauss = stepSize;

gaussSmooth = 1;
% nSets = size(densityPlot,3);
% nIters2run = size(densityPlot,4);

% nSetsTest = 1; %last 1/2 only for now

b = 50;
h = 50; %for trapz - height

if strcmp(dat(1:4),'trap') && length(dat)>10
    if strcmp(dat(1:6),'trapzS')
        spacingTrapz = trapzSpacing{str2double(dat(12))}; %trapzScaled1,2,3
        a=length(spacingTrapz); %trapz length1
        b=locRange(2)+1; %50 - trapz length2 - +1 to let density plot go from 1:50 (rather than 0:49)
        h=round(((locRange(2)+1)^2)/((a+b)/2))+1; %trapz height (area = 50^2; like square)
    elseif strcmp(dat(1:6),'trapzK')
        spacingTrapz = spacing(14:37);
        if boxSize==1.5
            spacingTrapz = spacing(14:37+length(14:37)*.5);
        elseif boxSize==2
            spacingTrapz = spacing(14:37+length(14:37));
        end
        a = length(spacingTrapz);
        b = length(spacing);
        h = length(spacing);
    end
    halfArea = (((a+b)/2)*h)/2;
    c = sqrt(((a^2)+(b^2))/2);
    % new to equalize area - Krupic
    hLeft  = floor(halfArea/((b+c)/2)); %bigger side
    hRight = ceil(halfArea/((a+c)/2))+1; %smaller side
    
    %equal number of points in trapz (301 each; nPoints in trap 877/2=438.5)
    hLeft = 14;
    hRight = 20;
else
    b=length(spacing);
    h=length(spacing);
    spacingTrapz = spacing;
end


%perm testing
% nPerm = 1000;
permPrc_gA_act = nan(nIters2run,4);
permPrc_gW_act = nan(nIters2run,4);
permPrc_gA_actNorm = nan(nIters2run,4);
permPrc_gW_actNorm = nan(nIters2run,4);

%% compute actNorm maps after 'training' (on a new test set of locations)
[trialsTest,~, rSeedTest,ntInSq] = createTrls(dat,nTrialsTest,locRange,useSameTrls,jointTrls,boxSize,h);

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

gA_actNormPerm = nan(nPerm,nIters2run);
gW_actNormPerm = nan(nPerm,nIters2run);
for iterI=1:nIters2run
    fprintf('Running iter %d\n',iterI);
        actTrl = zeros(nClus,nTrialsTest);
        densityPlotActTmp = zeros(b,h);
        densityPlotActUpd = zeros(b,h);
        
        if length(dat)==12 %sometimes clus not in shape, not counted, so fewer clusters than expected
            indTmp=find(densityPlot(:,:,end,iterI)>0);
            [muTrain(1:length(indTmp),1,iterI), muTrain(1:length(indTmp),2,iterI)] = find(densityPlot(:,:,end,iterI)>0); %find cluster positions - added >0 since counts nans as well
        else
        [muTrain(:,1,iterI), muTrain(:,2,iterI)] = find(densityPlot(:,:,end,iterI)>0); %find cluster positions - added >0 since counts nans as well
        end
        
        
        dist2Clus = sqrt(sum(reshape([muTrain(:,1,iterI)'-trialsTest(:,1), muTrain(:,2,iterI)'-trialsTest(:,2)].^2,nTrialsTest,nClus,2),3));
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
        end  %find closestC
        %densityPlotActNorm
        for iTrl = 1:nTrialsTest
            densityPlotActTmp(trialsTest(iTrl,1)+1, trialsTest(iTrl,2)+1)    = densityPlotActTmp(trialsTest(iTrl,1)+1, trialsTest(iTrl,2)+1)+ nansum(actTrl(:,iTrl));
            densityPlotActUpd(trialsTest(iTrl,1)+1, trialsTest(iTrl,2)+1)       = densityPlotActUpd(trialsTest(iTrl,1)+1, trialsTest(iTrl,2)+1)+1; %log nTimes loc was visited
        end
%         % need the nans for places that are not in shape -
%         for act:
%         densityPlotActTmp(densityPlotActTmp==0) =  nan;

        %turn 0s outside of the shape into nans
        if ~strcmp(dat(1:2),'sq') % if circ/trapz
            for i=1:length(ntInSq)
                densityPlotActTmp(ntInSq(i,1)+1,ntInSq(i,2)+1) = nan;
            end
        end
        densityPlotActNormTmp = densityPlotActTmp./densityPlotActUpd; %divide by number of times that location was visited

%         densityPlotAct(:,:,iterI)               =  densityPlotActTmp;
%         densityPlotActNorm(:,:,iterI)           =  densityPlotActNormTmp;

        %smooth - not nec but to compare with perm (which need smoothing)
        densityPlotAct(:,:,iterI)               =  imgaussfilt(densityPlotActTmp,gaussSmooth);
        densityPlotActNorm(:,:,iterI)           =  imgaussfilt(densityPlotActNormTmp,gaussSmooth);
        
        
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
            
            %right half of box
            aCorrMap = ndautoCORR(densityPlotActTmp(:,h-hRight:end));
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
            
            %right half of box
            aCorrMap = ndautoCORR(densityPlotActNormTmp(:,h-hRight:end));
            [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
            gA_actNorm(iterI,:,3) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
            [g,gdataA] = gridSCORE(aCorrMap,'wills',0);
            gW_actNorm(iterI,:,3) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
        end

        
        %permutation testing
        if doPerm
            for iPerm = 1:nPerm
                %             fprintf('Perm %d\n',iPerm);
                %             randInd=randperm(nTrialsTest); % randomise each point
                %20s perm - better way to code this?
                randInd20 = randperm(nTrialsTest/20); %atm has to be divisible by 20
                randInd = [];
                for i=1:length(randInd20)
                    randInd = [randInd randInd20(i):randInd20(i)+20];
                end
                
                actAllPerm = actTrl(:,randInd);
                
                densityPlotActPerm       = zeros(b,h);
                densityPlotActUpdPerm    = zeros(b,h);
                for iTrl = 1:nTrialsTest
                    densityPlotActPerm(trialsTest(iTrl,1)+1, trialsTest(iTrl,2)+1) = densityPlotActPerm(trialsTest(iTrl,1)+1, trialsTest(iTrl,2)+1)+ sum(actAllPerm(:,iTrl));
                    densityPlotActUpdPerm(trialsTest(iTrl,1)+1, trialsTest(iTrl,2)+1) = densityPlotActUpdPerm(trialsTest(iTrl,1)+1, trialsTest(iTrl,2)+1)+1; %log nTimes loc was visited
                end
                densityPlotActNormPerm   = densityPlotActPerm./densityPlotActUpdPerm; %divide by number of times that location was visited
                
                densityPlotActPerm = imgaussfilt(densityPlotActPerm,gaussSmooth);
                densityPlotActNormPerm = imgaussfilt(densityPlotActNormPerm,gaussSmooth);
                
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