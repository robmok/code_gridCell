clear all;
% close all;

% wd='/Users/robert.mok/Dropbox/Grid_cell_model';
wd='/Users/robertmok/Dropbox/Grid_cell_model';
cd(wd);

nClus   = 60;
nTrials = 2500; %how many locations in the box / trials - 2.5k ; 5k if reset

colgrey = [.5, .5, .5];
colors = distinguishable_colors(nClus); %function for making distinguishable colors for plotting

%box
nSteps = 50; %to define spacing beween each loc in box
locRange = 1;%1; from -locRange to locRange
stepSize=diff(linspace(-1,1,nSteps)); stepSize=stepSize(1); %smallest diff between locs


% parameters
epsMuOrig=.6;% %learning rate / starting learning rate %.6
deltaEpsMu = .96;% %change in learning rate over time (slow down with 'learning')

% deltaEpsMu = .99; % slower decrease in learning rate for expanding (if no
% reset)

% for saving simulations - multiple by values to save the files with params
epsMuOrig10=epsMuOrig*10;
deltaEpsMu100 = deltaEpsMu*100;

% testing with 60 clusters, 5k trials
% starting learning rate = .6
% deltaEpsMu:
    %.99   - SSE ~ 60 - fluctuating more at end, not stabilising
    %.98   - SSE ~.58
    %.96/7 - SSE ~ 56-7
    %.95   - SSE ~56-7 (decrease in epsmu curve looks gd, exponential, half near zero; 
    %.94   - SSE ~ 56 (curve gd, steeper, most near zero at end)
    %.935  - SSE ~ 59-60 (curve gd, steeper, maybe too steep? most near zero at end)

    
%diff number of clusters
% 30 clusters
% deltaEpsMu:
    %.99  - SSE ~120-124 - curve not near 0
    %.985 - SSE 120
    %.98  - SSE 116/117 - curve now steep but not all 0 at end
    %.97  - SSE ~ 117 - curve pretty steep
    %.95  - SSE ~ 114 - curve too steep, at 0 
    
    % although curves are different when less clusters, clusters look hexagonal anyway - less
    % clusters easier?
%40 cluster, .96; SSE~90; .97; SSE~88-91; .98; SSE~90
    

%different starting learning rate
%30 clusters
%.99, SSE = 124 (too shallow,), .98; SSE~118-123, .96; SSE~115 (steeper), .95, SSE~ 120 (v steep)


%reset learning rate?
resetEps=0; %0 - no, 1 - once halfway, 2 - twice (quarter way and half way)  % ATM resetting to halfway

% resetMag=2; %1 - back to orig, 2 - 75% of the orig (if resetEps=2, 75%
% the orig then 50% the orig) %- not using this yet


%define box / environement - random points in a box
box = 'square';

spacing=linspace(-1,1,101); 
switch box
    case 'square'
        trials=[randsample(spacing,nTrials,'true'); randsample(spacing,nTrials,'true')]';
    case 'rect'
        trials=[randsample(-1:diff(spacing(1:2)):2,nTrials,'true'); randsample(spacing,nTrials,'true')]';
end

% change box shape during learning
%rectangle
warpBox = 0; %1 or 0
warpType = 'sq2rect';
switch warpType
    case 'sq2rect'
        trialsExpand = [randsample(-1:diff(spacing(1:2)):2,nTrials/2,'true'); randsample(spacing,nTrials/2,'true')]';
    case 'rect2sq'
        trialsExpand = [randsample(spacing,nTrials/2,'true'); randsample(spacing,nTrials/2,'true')]'; %this doesn't work - actually, maybe this isn't a thing?
end

%for crossvalidation
nTrialsTest = nTrials;
% nTrialsTest = 5000; %keep it constant if want to compare with nTrials above
dataPtsTest = [randsample(linspace(-locRange,locRange,101),nTrialsTest,'true'); randsample(linspace(-locRange,locRange,101),nTrialsTest,'true')]'; % random points in a box

%%

%save simulations - cluster centres and tsse
saveDat=0;

nIter=1; %how many iterations

muAll    = nan(nClus,2,nTrials,nIter);
tsseTrls = nan(nTrials,nIter);
muEnd    = nan(nClus,2,nIter);

tic
for iterI = 1:nIter
    
    fprintf('iter %d \n',iterI);
    epsMu = epsMuOrig; %revert learning rate back to original if reset

    %initialise each cluster location
    mu = nan(nClus,2,nTrials);
    for clusterI = 1:nClus
        %     mu(clusterI,:,1) = -locRange + locRange.*2.*rand(1,2);  %random location in box - uniform distr from -1 to 1 (see help rand)
        mu(clusterI,:,1) = trials(randi(length(trials)),:);   %intiate each cluster with one data point - Forgy method
        
        % k means++?
        
    end

    %%
    clusInd  = 1:nClus;
    epsMuVec = zeros(nClus,1);
    updatedC = nan(nTrials,1);
    magChg   = nan(nTrials,1);
    epsMuAll = nan(nTrials,2);
    deltaMu  = nan(nClus,2,nTrials);
    distTrl  = nan(nTrials,nClus);

    for iTrl=1:length(trials)
        for iClus = 1:nClus,
            
            %if change size of box half way
            if iTrl == nTrials/2 && warpBox,
               trials(nTrials/2+1:end,:) = trialsExpand;
            end
            
            
            
            
            dist2Clus = sqrt(sum([mu(:,1,iTrl)'-trials(iTrl,1); mu(:,2,iTrl)'-trials(iTrl,2)].^2)); % vectorising euclid dist - sqrt(sum((a-b).^2)), since can't use xval method
            closestC=find(min(dist2Clus)==dist2Clus);
            if numel(closestC)>1, %if more than 1, randomly choose one
                closestC = randsample(closestC,1);
            end
            
            
            %moved learning rate change to up here - now update epsMu(iTrl), not epsMu(iTrl-1)
            
            %need to think how to update per cluster - since before just
            %defined x proprtion of trials. i suppose this is ok too - since
            %resetting them all a bit is better than just resetting one (prob
            %doesnt do anything).
            % the Q is whether should reset them all to a similar level - now
            % the ones that have been updated the most will have reset more and
            % move less...
            %         if iTrl == nTrials/2,
            %             updatedC = [];
            %         end
            
            %log which cluster has been updated
            updatedC(iTrl) = closestC;
            
            if iTrl~=1,
                %slow down learning rate over time
                magChg(iTrl)=nnz(updatedC==closestC); %magnitude of change proportional to how many times the cluster has moved
                %             epsMu(iTrl+1) = epsMu(1); %no slowing down
                %             epsMu(iTrl) = epsMu(iTrl-1)*deltaEpsMu;
                
                
                %resetting learning rate
                %still need to think if the amount i'm resetting to is
                %calculated correctly, and whether these are good numbers (prob
                %doesn't matter too much, but sth like high, med low sounds OK)
                if resetEps == 1,
                    if iTrl > nTrials*.5
                        magChg(iTrl)=nnz(updatedC==closestC)*.75;
                    end
                elseif resetEps == 2,
                    if iTrl > nTrials*.25 && iTrl <= nTrials*.5
                        magChg(iTrl)=nnz(updatedC==closestC)*.75;
                    elseif iTrl > nTrials*.5 && iTrl <= nTrials*.75
                        magChg(iTrl)=nnz(updatedC==closestC)*.25;
                    elseif iTrl > nTrials*.75 && iTrl <= nTrials
                        magChg(iTrl)=nnz(updatedC==closestC)*.25;
                    end
                    
                end
                                
                epsMu = epsMuOrig*deltaEpsMu^magChg(iTrl); %squared because it's deltaMu (e.g. .99) to the power of nTimes it was updated, making it slightly smaller each time
            else
                epsMu = epsMuOrig;
            end
            epsMuAll(iTrl,:) = [epsMu,closestC];
            
            epsMuVec = zeros(nClus,1);
            epsMuVec(closestC) = epsMu;
            deltaMu(:,1,iTrl) = epsMuVec.*(trials(iTrl,1)-mu(:,1,iTrl));
            deltaMu(:,2,iTrl) = epsMuVec.*(trials(iTrl,2)-mu(:,2,iTrl));
            % update mean estimates
            if iTrl~=length(trials) %no need to update for last trial +1)
                mu(:,1,iTrl+1) = mu(:,1,iTrl) + deltaMu(:,1,iTrl);
                mu(:,2,iTrl+1) = mu(:,2,iTrl) + deltaMu(:,2,iTrl);
            end
        end
        
        % compute sse on each trial with respect to 'all trials' - independent data - looks similar if you change 'dataPtsTest' to 'trials'; also useful if want to test less vs more trials on training so that you equate SSE by number of test points
        for iClus = 1:size(mu,1),
            distTrl(:,iClus)=sum([mu(iClus,1,iTrl)-dataPtsTest(:,1), mu(iClus,2,iTrl)-dataPtsTest(:,2)].^2,2);
        end
        [indValsTrl, indTrl]=min(distTrl,[],2); % find which clusters are points closest to
        
        sseTrl = zeros(iClus,1);
        for iClus = 1:size(mu,1),
            sseTrl(iClus)=sum(sum([mu(iClus,1,iTrl)-dataPtsTest(indTrl==iClus,1), mu(iClus,2,iTrl)-dataPtsTest(indTrl==iClus,2)].^2,2)); %distance from each cluster from training set to datapoints closest to that cluster
        end
        tsseTrls(iTrl,iterI)=sum(sseTrl);
        
        
        % idea: resetting the learning rate at some point to 'adapt' to new
        % environment. useful if morph the box shape in real time.
        
        %     % might also be interesting to know if this actually improves SSE, since
        %     % allows the first go/first few gos where things get stuck to improve
        %     if iTrl == round(length(trials)*.25) || iTrl == round(length(trials)*.5) || iTrl == round(length(trials)*.75)
        %     if iTrl == round(length(trials)*.5)
        %         epsMu(iTrl+1)=epsMu(1); %reset learning rate at some point, see what happens
        %     end
        %
        
        % %
        %     %     %could also reset to half the previous original learning rate each time...
        %     if iTrl == round(length(trials)*.25)
        %         epsMu(iTrl+1)=epsMu(1)*.75;
        %     elseif iTrl == round(length(trials)*.5)
        %         epsMu(iTrl+1)=epsMu(1)*.5;
        %     elseif iTrl == round(length(trials)*.75)
        %         epsMu(iTrl+1)=epsMu(1)*.25;
        %     end
        
        %     if iTrl == round(length(trials)*.25)
        %         epsMu(iTrl+1)=epsMu(1)*.75;
        %     elseif iTrl == round(length(trials)*.33)
        %         epsMu(iTrl+1)=epsMu(1)*.66;
        %     elseif iTrl == round(length(trials)*.5)
        %         epsMu(iTrl+1)=epsMu(1)*.5;
        %     elseif iTrl == round(length(trials)*.66)
        %         epsMu(iTrl+1)=epsMu(1)*.33;
        %     elseif iTrl == round(length(trials)*.75)
        %         epsMu(iTrl+1)=epsMu(1)*.25;
        %     end
        
        %     %put a limit of lowest learning rate
        %     if epsMu(iTrl+1)<=.1,
        %         epsMu(iTrl+1) = .1;
        %     end
        
    end
    muEnd(:,:,iterI)=mu(:,:,end);
    muAll(:,:,:,iterI) = mu;
end
toc

if saveDat,
    % save simulations
    % might not need all tsseTrls if not viewing them or epsMuAll
%     save(sprintf('covering_map_dat_%dclus_%dktrls_eps%d_deltaMu%d',nClus,nTrials/1000,epsMuOrig10,deltaEpsMu100),'muEnd', 'tsseTrls', 'epsMuAll','nIter');
    fname = sprintf('covering_map_dat_%dclus_%dktrls_eps%d_deltaMu%d',nClus,nTrials,epsMuOrig10,deltaEpsMu100);
    if resetEps,
       fname = [fname 'resetHalfwayHalfEps']; 
    end
    save(fname,'muEnd','muAll', 'tsseTrls', 'epsMuAll','nIter');
end

%% Load in data, plot

% save plots?
saveplots=0;


%set to load in data (note divide by value in variable)
if 1
    epsMuOrig10 = 6;
    deltaEpsMu100 = 96;
    nClus = 30;
    colors = distinguishable_colors(nClus); % if changing the number of clusters
    nTrials = 2500;
    load(sprintf('covering_map_dat_%dclus_%dktrls_eps%d_deltaMu%d',nClus,nTrials,epsMuOrig10,deltaEpsMu100));
end


figure; plot(tsseTrls);
title(sprintf('SSE for each iteration (nClus=%d)',nClus));
if saveplots,
    fname=[wd, sprintf('/covering_map_%d_clus_%d_trls_tsseTrls',nClus,nTrials)];
    saveas(gcf,fname,'png');
end

%plot change in learning rate for each cluster:
% figure;
% plot(epsMuAll(:,1),'.');


% plot lowest and highest SSE
tsseIter = tsseTrls(end,:);

[indVal, indSSE1] = sort(tsseIter);
[y, indSSE2] = sort(indSSE1); %use indSSE2 to sort for plotting

iToPlot=[1,2,3,nIter-2,nIter-1,nIter];
figure;
for i = 1:6
    subplot(2,3,i); hold on;
%     plot(dataPtsTest(:,1),dataPtsTest(:,2),'.','Color',colgrey,'MarkerSize',2); hold on;
    voronoi(squeeze(muEnd(:,1,indSSE2==iToPlot(i))),squeeze(muEnd(:,2,indSSE2==iToPlot(i))),'k')
    for iClus = 1:nClus
        plot(squeeze(mean(muEnd(iClus,1,(indSSE2==iToPlot(i))))),squeeze(mean(muEnd(iClus,2,(indSSE2==iToPlot(i))))),'.','Color',colors(iClus,:),'MarkerSize',25); hold on; %plot cluster final point
    end
%     xlim([min(trials(:,1))-.1,max(trials(:,1))+.1]); ylim([min(trials(:,2))-.1,max(trials(:,2))+.1]);
 
    hold on;
    if i==2,
        title(sprintf('Lowest 3 and Highest 3 SSE cluster locations (nClus=%d)',nClus));
    end
end
if saveplots,
    fname=[wd, sprintf('/covering_map_%d_clus_locs_top_bottom_3_%d_trls',nClus,nTrials)];
    saveas(gcf,fname,'png');
end

%% plot cluster centers over time

toPlot = find(indSSE2==1); %lowest sse

% figure;   
figure('units','normalized','outerposition',[0 0 1 1]);
iPlot = 1; subplot(3,4,iPlot); hold on;%subplot
for iTrl = 1:nTrials
    %     if mod(iTrl,50)==0, fprintf('Trial %d \n',iTrl); end
    %     % plot all points
    %     plot(squeeze(mu(1,1,iTrl)),squeeze(mu(1,2,iTrl)),'.','Color',colors(1,:),'MarkerSize',10); hold on;
    %     plot(squeeze(mu(2,1,iTrl)),squeeze(mu(2,2,iTrl)),'.','Color',colors(2,:),'MarkerSize',10); hold on;
    %     plot(squeeze(mu(3,1,iTrl)),squeeze(mu(3,2,iTrl)),'.','Color',colors(3,:),'MarkerSize',10); hold on;
    %     plot(squeeze(mu(4,1,iTrl)),squeeze(mu(4,2,iTrl)),'.','Color',colors(4,:),'MarkerSize',10); hold on;
    %     plot(squeeze(mu(5,1,iTrl)),squeeze(mu(5,2,iTrl)),'.','Color',colors(5,:),'MarkerSize',10); hold on;
    %     drawnow;
    
    if mod(iTrl,500)==0
        iPlot=iPlot+1;
        subplot(3,4,iPlot); hold on;
        voronoi(squeeze(muAll(:,1,iTrl,toPlot)),squeeze(muAll(:,2,iTrl,toPlot)),'k')
    end
    if mod(iTrl,50)==0, %plot centers after x trials
        for i=1:nClus
            plot(squeeze(muAll(i,1,iTrl,toPlot)),squeeze(muAll(i,2,iTrl,toPlot)),'.','Color',colors(i,:),'MarkerSize',20); hold on;
        end
        drawnow;
    end
end
% xlim([min(trials(:,1))-.1,max(trials(:,1))+.1]); ylim([min(trials(:,2))-.1,max(trials(:,2))+.1]);

%% plot final centres 
for i = 1:nIter
    figure;
    for iClus = 1:nClus
        plot(squeeze(muEnd(iClus,1,i)),squeeze(muEnd(iClus,2,i)),'.','Color',colors(iClus,:),'MarkerSize',25); hold on; %plot cluster final point
    end
    xlim([min(trials(:,1))-.1,max(trials(:,1))+.1]); ylim([min(trials(:,2))-.1,max(trials(:,2))+.1]);
end
voronoi(squeeze(muEnd(:,1,i)),squeeze(muEnd(:,2,i)))

