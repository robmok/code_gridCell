function [densityPlotAct,densityPlotActNorm,gA_act,gA_actNorm,gW_act,gW_actNorm, rSeedTest] = gridnessTestData_Perm(densityPlot,dat,locRange,nClus,nTrialsTest)


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
nSets = size(densityPlot,3);
nIter = size(densityPlot,4);

nSetsTest = 2; %last 2 only for now


b = 50;
h = 50; %for trapz - height

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
        a = length(spacingTrapz);
        b = length(spacing);
        h = length(spacing);
    end
    halfArea = (((a+b)/2)*h)/2;
    c = sqrt(((a^2)+(b^2))/2);
    % new to equalize area - Krupic
    hLeft  = floor(halfArea/((b+c)/2)); %bigger side
    hRight = ceil(halfArea/((a+c)/2))+1; %smaller side
else
    b=length(spacing);
    h=length(spacing);
    spacingTrapz = spacing;
end

%% compute actNorm maps after 'training' (on a new test set of locations)
[trialsTest,~, rSeedTest] = createTrls(dat,nTrialsTest,locRange,useSameTrls,jointTrls,boxSize,h);

densityPlotAct     = zeros(b,h,nSetsTest,nIter);
densityPlotActNorm = zeros(b,h,nSetsTest,nIter);
gA_act = nan(nSetsTest,nIter,9);
gW_act = nan(nSetsTest,nIter,9);
gA_actNorm = nan(nSetsTest,nIter,9);
gW_actNorm = nan(nSetsTest,nIter,9);
if strcmp(dat(1:4),'trap') %if trapz - compute gridness of left/right half of boxes too
    gA_act = nan(nSetsTest,nIter,9,3);
    gW_act = nan(nSetsTest,nIter,9,3);
    gA_actNorm = nan(nSetsTest,nIter,9,3);
    gW_actNorm = nan(nSetsTest,nIter,9,3);
end
muTrain = nan(nClus,2,nSetsTest,nIter);

for iterI=1:nIter
    fprintf('Running iter %d\n',iterI);
    for iSet=1:nSetsTest%1:nSets
%         fprintf('Running set %d\n',iSet);

        actTrl = zeros(nClus,nTrialsTest);
        densityPlotActUpd = zeros(b,h);
        
        [muTrain(:,1,iSet,iterI), muTrain(:,2,iSet,iterI)] = find(densityPlot(:,:,nSets-2+iSet,iterI)); %find cluster positions
        
        dist2Clus = sqrt(sum(reshape([muTrain(:,1,iSet,iterI)'-trialsTest(:,1), muTrain(:,2,iSet,iterI)'-trialsTest(:,2)].^2,nTrialsTest,nClus,2),3));
        closestC = nan(1,nTrialsTest);
        for iTrl = 1:nTrialsTest
            closestTmp = find(min(dist2Clus(iTrl,:))==dist2Clus(iTrl,:));
            if numel(closestTmp)>1
                closestC(iTrl) = randsample(closestTmp,1);
            else
                closestC(iTrl) = closestTmp;
            end
            %compute activation
            actTrl(closestC(iTrl),iTrl)=mvnpdf(trialsTest(iTrl,:),muTrain(closestC(iTrl),:,iSet),eye(2)*sigmaGauss); % save only the winner
        end  %find closestC
        %densityPlotActNorm
        for iTrl = 1:nTrialsTest
            densityPlotAct(trialsTest(iTrl,1)+1, trialsTest(iTrl,2)+1,iSet,iterI)    = densityPlotAct(trialsTest(iTrl,1)+1, trialsTest(iTrl,2)+1,iSet,iterI)+ nansum(actTrl(:,iTrl));
            densityPlotActUpd(trialsTest(iTrl,1)+1, trialsTest(iTrl,2)+1) = densityPlotActUpd(trialsTest(iTrl,1)+1, trialsTest(iTrl,2)+1)+1; %log nTimes loc was visited
        end
        densityPlotActNorm(:,:,iSet,iterI) = densityPlotAct(:,:,iSet,iterI)./densityPlotActUpd; %divide by number of times that location was visited
        densityPlotActSm                   = imgaussfilt(densityPlotAct(:,:,iSet,iterI),gaussSmooth);
        densityPlotActNormSm               = imgaussfilt(densityPlotActNorm(:,:,iSet,iterI),gaussSmooth);
        
        aCorrMap = ndautoCORR(densityPlotActSm);
        [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
        gA_act(iSet,iterI,:,1) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
        
        [g,gdataA] = gridSCORE(aCorrMap,'wills',0);
        gW_act(iSet,iterI,:,1) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
        
        %compute gridness over normalised activation map - normalised by times loc visited
        aCorrMap = ndautoCORR(densityPlotActNormSm);
        [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
        gA_actNorm(iSet,iterI,:,1) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
        
        [g,gdataA] = gridSCORE(aCorrMap,'wills',0);
        gW_actNorm(iSet,iterI,:,1) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
        
        if  strcmp(dat(1:4),'trap')
            %left
            aCorrMap = ndautoCORR(densityPlotActSm(:,1:hLeft));
            [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
            gA_act(iSet,iterI,:,2) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
            [g,gdataA] = gridSCORE(aCorrMap,'wills',0);
            gW_act(iSet,iterI,:,2) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
            
            %right half of box
            aCorrMap = ndautoCORR(densityPlotActSm(:,h-hRight:end));
            [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
            gA_act(iSet,iterI,:,3) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
            [g,gdataA] = gridSCORE(aCorrMap,'wills',0);
            gW_act(iSet,iterI,:,3) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
            
            %left
            aCorrMap = ndautoCORR(densityPlotActNormSm(:,1:hLeft));
            [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
            gA_actNorm(iSet,iterI,:,2) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
            [g,gdataA] = gridSCORE(aCorrMap,'wills',0);
            gW_actNorm(iSet,iterI,:,2) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
            
            %right half of box
            aCorrMap = ndautoCORR(densityPlotActNormSm(:,h-hRight:end));
            [g,gdataA] = gridSCORE(aCorrMap,'allen',0);
            gA_actNorm(iSet,iterI,:,3) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
            [g,gdataA] = gridSCORE(aCorrMap,'wills',0);
            gW_actNorm(iSet,iterI,:,3) = [gdataA.g_score, gdataA.orientation, gdataA.wavelength, gdataA.radius, gdataA.r'];
        end
    end
end
end