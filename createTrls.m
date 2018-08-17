function [trials,dataPtsTest, rSeed,ntInSq] = createTrls(dat,nTrials,locRange,useSameTrls,jointTrls,catsInfo,rSeed)

spacing =linspace(locRange(1),locRange(2),locRange(2)+1); 
spacingOrig=spacing;
stepSize=diff(spacing(1:2)); nSteps = length(spacing);


if exist('rSeed','var')
    rng(rSeed); %to recreate trials
else
    rSeed=rng; % seed so can regenerate the trials later; do: rng(rSeed(iterI)); then run trials = ...
end
% first get points in the shape
switch dat
    case 'square' %square
        trials = nan(nTrials,2);
        sq=linspace(locRange(1),locRange(2),nSteps);
        sqPts=[];
        for i=1:length(sq)
            for j=1:length(sq)
                sqPts = [sqPts; [sq(i), sq(j)]];
            end
        end
        shapePts = sqPts;
    case 'catLearn'
        % draw points from 2 categories (gaussian) from a 2D feature space
        nTrialsCat = floor(nTrials/catsInfo.nCats); % points to sample
        catCentres = nan(catsInfo.nCats,2);
        datPtsGauss = nan(nTrialsCat,2,catsInfo.nCats);
        for iCat = 1:catsInfo.nCats
            if catsInfo.msExample %if generating example for ms - 2 gaussians at opposite sides of the square env
                catCentres = [15, 35; 35, 15];
            else
                catCentres(iCat,:)=randsample(locRange(1)+10:locRange(2)-10,2,'true'); % ±10 so category centres are not on the edge
            end
            datPtsGauss(:,:,iCat) = round(repmat(catCentres(iCat,:),nTrialsCat,1) + randn(nTrialsCat,2)*catsInfo.R); % key - these are the coordinates of the points
        end
        
        trials = reshape(datPtsGauss,numel(datPtsGauss)/2,2);
        %if odd number of categories, create the remaining number of trials 
        trlsRem=nTrials-nTrialsCat*catsInfo.nCats;
        if trlsRem>0 % if remainder
            extraTrlsCat=nan(1,trlsRem);
            extraTrlsCat=randi(catsInfo.nCats,1,trlsRem); %select cats randomly
            extraTrls   = nan(trlsRem,2);
            for iExtra=1:length(extraTrlsCat)
                extraTrls(iExtra,:) = round(catCentres(extraTrlsCat(iExtra),:) + randn(1,2)*catsInfo.R);
            end
            trials = [trials; extraTrls];
        end
        trials = trials(randperm(length(trials)),:);
%         for iCat = 1:catsInfo.nCats
%             catCentres(iCat,:)=randsample(locRange(1)+10:locRange(2)-10,2,'true'); % ±10 so category centres are not on the edge
%             datPtsGauss(:,:,iCat) = round(repmat(catCentres(iCat,:),nTrialsCat,1) + randn(nTrialsCat,2)*catsInfo.R); % key - these are the coordinates of the points
%         end
%         dataPtsTest = reshape(datPtsGauss,nTrials,2);
%         dataPtsTest = dataPtsTest(randperm(length(dataPtsTest)),:);
        dataPtsTest=[];
    case 'rect' 
        %EDIT
%         trials      = [randsample(locRange(1):diff(spacing(1:2)):locRange(2)/1.5,nTrials,'true'); randsample(spacing,nTrials,'true')]';
%         dataPtsTest = [randsample(locRange(1):diff(spacing(1:2)):locRange(2)/1.5,nTrials,'true'); randsample(spacing,nTrials,'true')]';
    case 'trapzKrupic'
        %if scale to Krupic:
        spacing = spacing(14:37);
%         trapY=locRange(2).*trapmf(spacing,[spacing(1), spacing(round(length(spacing)*.25)), spacing(round(length(spacing)*.75)),spacing(end)]); %old, a=18
        trapY=locRange(2).*trapmf(spacing,[spacing(1), spacingOrig(round(length(spacingOrig)*.46)), spacingOrig(round(length(spacingOrig)*.54)),spacing(end)]);% new, a=5; 23:27
        trapX=spacing;
        trapPts=[];
        for i=1:length(trapY)
            trapPts = [trapPts, [repmat(trapX(i),1,length(0:stepSize:trapY(i))); 0:stepSize:trapY(i)]];
        end
        shapePts = trapPts';
    case 'trapzKfrmSq1'
        spacing = spacing(14:37);
%         trapY=locRange(2).*trapmf(spacing,[spacing(1), spacing(round(length(spacing)*.25)), spacing(round(length(spacing)*.75)),spacing(end)]);
trapY=locRange(2).*trapmf(spacing,[spacing(1), spacingOrig(round(length(spacingOrig)*.46)), spacingOrig(round(length(spacingOrig)*.54)),spacing(end)]);% new, a=5; 23:27

        trapX=spacing;
        trapPts=[];
        for i=1:length(trapY)
            trapPts = [trapPts, [repmat(trapX(i),1,length(0:stepSize:trapY(i))); 0:stepSize:trapY(i)]];
        end
        shapePts = trapPts';
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
        shapePts = circPts;
end

if ~strcmp(dat(1:3),'cat')
%select points from specified shape - randomly sample the space, or joint trials
if ~useSameTrls && ~jointTrls
    trialInd=randi(length(shapePts),nTrials,1);
    trials=shapePts(trialInd,:);
elseif ~useSameTrls && jointTrls  %joined trials
    %move step sizes - possible step sizes; if larger, e.g. double, then need bigger steps?
    
    selStepSiz = [-stepSize*4,-stepSize*2,-stepSize,-stepSize,0,stepSize,stepSize,stepSize*2,stepSize*4]; %stepSiz

    % select points trial by trial, if not in shape go back into the shape
    trials       =  nan(nTrials,2);
    trials(1,:)  =  shapePts(randi(length(shapePts),1,1),:);
    if ~strcmp(dat(1:4),'trap') % regular shapes (circ, square, rect)
        for i=2:nTrials
            locRangeX = locRange;
            locRangeY = locRange;
            moveDir=randsample(selStepSiz,2); %move in a random direction or stay
            %if edge, stay or move
            while ~any(trials(i-1,1)+moveDir(1)==shapePts(:,1) & trials(i-1,2)+moveDir(2)==shapePts(:,2)) % if not in the shape
                if trials(i-1,1)+moveDir(1) < median(locRangeX(1):locRangeX(2)) %if less than half the box, +
                    moveDir(1) = randsample(selStepSiz(end-4:end),1);
                end
                if trials(i-1,2)+moveDir(2) < median(locRangeY(1):locRangeY(2))
                    moveDir(2) = randsample(selStepSiz(end-4:end),1);
                end
                if trials(i-1,1)+moveDir(1) > median(locRangeX(1):locRangeX(2))
                    moveDir(1) = randsample(selStepSiz(1:end-4),1);
                end
                if trials(i-1,2)+moveDir(2) > median(locRangeY(1):locRangeY(2))
                    moveDir(2) = randsample(selStepSiz(1:end-4),1);
                end
%                 if trials(i-1,1)+moveDir(1) == median(locRangeX(1):locRangeX(2)) %if middle, move a little or stay
%                     moveDir(1) = randsample([-1,0,1],1);
%                 end
%                 if trials(i-1,2)+moveDir(2) == median(locRangeY(1):locRangeY(2)) %if middle, move a little or stay
%                     moveDir(2) = randsample([-1,0,1],1);
%                 end
            end
            trials(i,:)=trials(i-1,:)+moveDir;% add 1 or -1
        end
        
    else % if trapz, y-axis is a bit different (always + if not in area, unless <0)
        % turns out too often going 'down' when out of trapz, and
        % rarely going up. reduce this by having more 0s/'stays', and
        % one more selStepsize in the + direction for up
        
        %trapz - left right don't make too big steps - since height is more than
        %double the width
        if strcmp(dat(1:6),'trapzK') %trapzKfrmSq as well
            selStepSizLR = [-stepSize,-stepSize,0,0,0,stepSize,stepSize]; %stepSizLR - combined with stepSiz
        end
        
        locRangeX = [min(trapX), max(trapX)];
        for i=2:nTrials
            moveDir=randsample(selStepSiz,2); %move in a random direction or stay
%             moveDir(1)=randsample(selStepSiz,1); %move in a random direction or stay (x-axis: left/right)
%             moveDir(2)=randsample(selStepSiz,1); %move in a random direction or stay
%             moveDir(2)=randsample([selStepSiz selStepSiz(end-4:end-2)],1); %move in a random direction or stay; added -  more likely to go up (y-axis) - trapzKrupic 
            
            %if edge, stay or move
            while ~any(trials(i-1,1)+moveDir(1)==shapePts(:,1) & trials(i-1,2)+moveDir(2)==shapePts(:,2)) % if not in the shape
                if trials(i-1,1)+moveDir(1) < median(locRangeX(1):locRangeX(2)) %if less than half the box, +
                    moveDir(1) = randsample(selStepSizLR(end-3:end),1);
                end
                if trials(i-1,2)+moveDir(2) < 0 % only if goes 'below' the shape
                    moveDir(2) = randsample(selStepSiz(end-4:end),1);
                end
                if trials(i-1,1)+moveDir(1) > median(locRangeX(1):locRangeX(2))
                    moveDir(1) = randsample(selStepSizLR(1:end-3),1);
                end
                if trials(i-1,2)+moveDir(2) > 0 % if above 0 and out
                    moveDir(2) = randsample(selStepSiz,1); %this works fine
%                     moveDir(2) = 0;% bias it not to go down
                end
                if trials(i-1,1)+moveDir(1) == median(locRangeX(1):locRangeX(2)) %if middle, move a little or stay
                    moveDir(1) = randsample([-1,0,1],1);
                end
            end
            trials(i,:)=trials(i-1,:)+moveDir;% add 1 or -1
        end
    end
end

%dataPtsTest - for initiating cluster centres, better to not select from the
%joint trials
trialIndTest = randi(length(shapePts),nTrials,1);
dataPtsTest  = shapePts(trialIndTest,:);
end
if nargout > 3 && ~strcmp(dat(1:3),'cat')
    sq=linspace(locRange(1),locRange(2),nSteps);
    sqPts=[];
    for i=1:length(sq)
        for j=1:length(sq)
            sqPts = [sqPts; [sq(i), sq(j)]];
        end
    end
    ntInSq = setdiff(sqPts,shapePts,'rows');
else
    ntInSq =[];
end
end
% %checking agent's visited locations
% densityPlot=zeros(locRange(2)+1,locRange(2)+1);
% % densityPlot=zeros(50,62); %trapzScaled1
% for iTrl = 1:nTrials
%     densityPlot(trials(iTrl,1)+1,trials(iTrl,2)+1) = densityPlot(trials(iTrl,1)+1,trials(iTrl,2)+1)+1;
% end
% figure;imagesc(densityPlot)

%need add this if want to check
% figure; hist(moveAll(:,2),50)
% nnz(moveAll(:,2)<0)
% nnz(moveAll(:,2)>0)