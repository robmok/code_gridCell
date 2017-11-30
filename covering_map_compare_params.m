clear all;

wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
% wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
cd(wd);

codeDir = [wd '/code_gridCell'];
saveDir = [wd '/data_gridCell'];
addpath(codeDir); addpath(saveDir);
addpath(genpath([wd '/gridSCORE_packed']));

% comparing parameters

%make test data - keep same for comparison
nSteps = 50;
locRange = [0, nSteps-1]; 
nTrials = 40000;
nTrialsTest = nTrials;
% nTrialsTest = 50000;
dataPtsTest = [randsample(linspace(locRange(1),locRange(2),nSteps),nTrialsTest,'true'); randsample(linspace(locRange(1),locRange(2),nSteps),nTrialsTest,'true')]'; % random points in a box

%% Compute the cluster centres from the desnity map, then compute SSE and rank them

doPlot = 0; %plot - note if loading in many value, this will make too many plots, can crash

gaussSmooth=1; %smooth maps by x value
% dsFactor=5; %downsample the map so faster computing autocorr map

%for computing density map
nTrlsToUse = 10000; %
toTrlN     = nTrials; %-10000; %nTrials if from nTrlstoUse to end,
fromTrlI   = toTrlN-nTrlsToUse+1;


%define which parmeters to test (if not testing, keep at 1 value)

%how to test several; or should test one at a time?
%%%%%%%
alphaVals = 2;
%[0,2,4,6,8]; 

stochasticType = 1; %1, 2, 3; % stochasticVals=[1,2,3] - later

cVals = ([.1/nTrials, .25/nTrials, .5/nTrials, 2/nTrials, 3/nTrials, 5/nTrials, 10/nTrials, 20/nTrials]); %.1 and stoch = towards middle
% cVals = [2/nTrials, 3/nTrials, 5/nTrials, 10/nTrials];
% cVals = [3/nTrials];
cVals = round(cVals.*1000000);

%%%%%%%
%set which one to test
%%%%%%
testVals = cVals; %alphaVals, cVals, stochasticVals
nTestVals = length(testVals);

nClus = 20; %20, 30
colors = distinguishable_colors(nClus);
nTrials = 40000;
nIter = 10;

%warpBox=0;

clusMu=nan(nClus,2,nIter,length(testVals));
sse=nan(1,nClus);
tsse=nan(nIter,length(testVals));
devAvgSSE=nan(nClus,nIter,length(testVals));

spacing = linspace(locRange(1),locRange(2),locRange(2)+1);
densityPlot = zeros(length(spacing),length(spacing),nIter,nTestVals);
aCorrMap = nan(length(spacing)*2-1,length(spacing)*2-1,nIter,nTestVals);


for iTestVal = 1:nTestVals
%     iPlot=0; %for subplot below
    
    %set params, load in data
    epsMuOrig1000 = 75; %75, 100
    
    %alpha
    %     alpha10 = alphaVals(iTestVal);
    alpha10 = alphaVals;

    %stochastic parameters
    c =cVals(iTestVal);
    fname = [saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_alpha%d_stype%d_cVal%d_%diters*',nClus,nTrials,epsMuOrig1000,alpha10,stochasticType,c,nIter)];
    
    %if ran more than 1 set of iters, hv multiple files w diff end bit;load
    f = dir(fname); filesToLoad = cell(1,length(f));
    muAllTmp={};
    nIterCount=0;
    for iF = 1:length(f)
%         filesToLoad{iF} = f(iF).name;
%         fnames = [fnames filesToLoad{iF}];
        load(f(iF).name);
        muAllTmp{iF}=muAll;
        nIterCount = [nIterCount, nIter+nIterCount(end)]; %count number of iters to index below when merging
    end
    %merge
    clear muAll
    for iF = 1:length(f)
        muAll(:,:,:,nIterCount(iF)+1:nIterCount(iF+1)) = muAllTmp{iF};
    end
   
    %calculate
    for iterI = 1:size(muAll,4) %nIter
        %compute density map
        clus = round(muAll(:,:,fromTrlI:toTrlN,iterI));
        densityPlotClus     = zeros(length(spacing),length(spacing),nClus);
        densityPlotClusSmth = zeros(length(spacing),length(spacing),nClus);
        for iClus=1:nClus
            for iTrl=1:size(clus,3)
                %                 densityPlot(clus(i,1,iTrl),clus(i,2,iTrl),iterI,iTestVal)=densityPlot(clus(i,1,iTrl),clus(i,2,iTrl),iterI,iTestVal)+1; % works, but better way / faster to vectorise?
                densityPlotClus(clus(iClus,1,iTrl),clus(iClus,2,iTrl),iClus) = densityPlotClus(clus(iClus,1,iTrl),clus(iClus,2,iTrl),iClus)+1;
            end
            %find peaks
            densityPlotClusSmth(:,:,iClus)=imgaussfilt(densityPlotClus(:,:,iClus),gaussSmooth);
            [peakX, peakY] = find(densityPlotClusSmth(:,:,iClus)==max(max((densityPlotClusSmth(:,:,iClus)))));
            clusMu(iClus,:,iterI) = [peakX, peakY];
            
            %make combined (grid cell) plot
            densityPlot(:,:,iterI,iTestVal) = sum(densityPlotClus,3);
            densityPlot(:,:,iterI,iTestVal) = imgaussfilt(densityPlot(:,:,iterI,iTestVal),gaussSmooth); %smooth
        end
        
        
        
        % compute sse with respect to all test data points
        distTrl=[];
        for iClus = 1:size(clusMu,1)
            distTrl(:,iClus)=sum([clusMu(iClus,1,iterI)-dataPtsTest(:,1), clusMu(iClus,2,iterI)-dataPtsTest(:,2)].^2,2);
        end
        [indValsTrl, indTmp]=min(distTrl,[],2); % find which clusters are points closest to
        for iClus = 1:size(clusMu,1)
            sse(iClus)=sum(sum([clusMu(iClus,1,iterI)-dataPtsTest(indTmp==iClus,1), clusMu(iClus,2,iterI)-dataPtsTest(indTmp==iClus,2)].^2,2)); %distance from each cluster from training set to datapoints closest to that cluster
        end
        tsse(iterI,iTestVal)=sum(sse);
        
        %compute 'spreaded-ness' - variance of SE across clusters is a measure
        %of this, assuing uniform data points
        devAvgSSE            = sse-mean(sse);
        stdAcrossClus(iterI,iTestVal) = std(devAvgSSE);
        varAcrossClus(iterI,iTestVal) = var(devAvgSSE);
        
        %compute autocorrelation map
        
        %downsample
%         dPlotTmp=densityPlot(:,:,iterI,iTestVal);
%         ind=1:dsFactor:size(dPlotTmp,1);
%         clear dPlotTmpDS
%         for iRow = 1:size(dPlotTmp,1)/dsFactor
%             for iCol = 1:size(dPlotTmp,1)/dsFactor
%                 dPlotTmpDS(iCol,iRow) = dPlotTmp(ind(iRow),ind(iCol));
%             end
%         end
%         densityPlotDS(:,:,iterI,iTestVal) =imgaussfilt(dPlotTmpDS,gaussSmooth);

        aCorrMap(:,:,iterI,iTestVal) = ndautoCORR(densityPlot(:,:,iterI,iTestVal));
        
        %plot
        if doPlot
            figure;
            subplot(2,2,1);
            imagesc(densityPlot(:,:,iterI,iTestVal));
            subplot(2,2,2);
        end
        [g,gdataA] = gridSCORE(aCorrMap(:,:,iterI,iTestVal),'allen',doPlot);
        if doPlot, subplot(2,2,3); end
        [g,gdataW] = gridSCORE(aCorrMap(:,:,iterI,iTestVal),'wills',doPlot);
        
        
        gA_g(iterI,iTestVal)   = gdataA.g_score;
        gA_o(iterI,iTestVal)   = gdataA.orientation;
        gA_wav(iterI,iTestVal) = gdataA.wavelength;
        gA_rad(iterI,iTestVal) = gdataA.radius;
        
    end
    gA{iTestVal}=gdataA;
    gW{iTestVal}=gdataW;
end


% % plot
% for iTestVal = 1:nTestVals
%     figure;
%     for iterI = 1:nIter
%         subplot(2,2,1);
%         imagesc(densityPlot(:,:,iterI,iTestVal));
%         subplot(2,2,2);
%         [g,gdata] = gridSCORE(aCorrMap(:,:,iterI,iTestVal),'allen'); %allen or wills
%         subplot(2,2,3);
%         [g,gdata] = gridSCORE(aCorrMap(:,:,iterI,iTestVal),'wills'); %allen or wills
%     end   
% end
%% plot to compare
% the way the files are set up & loaded in loses info about which params.
% however can recover it by using the same loops when loading in the param
% names? e.g. make a variable with the param names, then plot them

barpos = .25:.5:.5*nTestVals;
barw   = .25;

%compare with bar plots, etc.

dat=gA_g;

figure; hold on;
mu=mean(dat,1);
sm=std(dat)./sqrt(size(dat,1));
ci=sm.*tinv(.025,size(dat,1)-1); %compute conf intervals
for i=1:length(mu)
    h(i) = bar(barpos(i),mu(i),barw);
end
% he = errorbar(barpos,mu,sm,-sm,'k.');
he = errorbar(barpos,mu,ci,'k.');
xlim([barpos(1)-.5, barpos(end)+.5]);

% [h p c s] = ttest(gA_g(:,2),gA_g(:,4))




dat=gA_o;

figure; hold on;
mu=mean(dat,1);
sm=std(dat)./sqrt(size(dat,1));
ci=sm.*tinv(.025,size(dat,1)-1);
for i=1:length(mu)
    h(i) = bar(barpos(i),mu(i),barw);
end
% he = errorbar(barpos,mu,sm,-sm,'k.');
he = errorbar(barpos,mu,ci,'k.');
xlim([barpos(1)-.5, barpos(end)+.5]);



%univariate scatter?





% plot actual data

% for iF=1:length(gA)
%    for iterI = 1:length(gA{iF})
%         
%        %get data from gA and gW to plot - are they plotting the same
%        thing? looks diff sometimes
% 
%    end
% end


% % plot (from gridSCORE.m)
% imc = imagesc(im);
% set(imc,'alphadata',msk);
% hold on
% plot(gdata.near_peaks(:,1),gdata.near_peaks(:,2),'kx','MarkerSize',10);
% title(sprintf('g = %.2f, s = %.2f, r = %.2f, o = %.2f',g,gdata.wavelength,gdata.radius,gdata.orientation));
% drawLine(L);
% caxis([0 nanmax(imcent(:))])
% daspect([1 1 1]);
% axis off


%%
% % plot cluster centres, print their sse and stdAcrossClus
% for iterI=1:nIter
%     % figure('units','normalized','outerposition',[0 0 1 1]); hold on; iPlot=0;
%     figure; hold on; iPlot=0;
%     for iTestVal = 1:length(testVals)
%         iPlot = iPlot+1;
%         subplot(2,2,iPlot);
%         scatter(clusMu(:,1,iterI,iTestVal),clusMu(:,2,iterI,iTestVal),20e+2,colors,'.')
%         xlim(locRange); ylim(locRange);
%         title(sprintf('%d alpha, %.1f tsse, %.2f. spreadedness', testVals(iTestVal), tsse(iterI,iTestVal), stdAcrossClus(iterI,iTestVal)))
% %         title(sprintf('%d cVal, %.1f tsse, %.2f. spreadedness', testVals(iDataSet), tsse(iterI,iDataSet), stdAcrossClus(iterI,iDataSet)))
%     end
% end


