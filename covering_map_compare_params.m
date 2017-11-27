clear all;

% wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
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

doPlot = 1; %plot - note if loading in many value, this will make too many plots, can crash

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

cVals = ([.1/nTrials, .25/nTrials, .5/nTrials, 2/nTrials, 3/nTrials, 5/nTrials, 10/nTrials, 20/nTrials]);
% cVals = [2/nTrials, 3/nTrials, 5/nTrials, 10/nTrials];
cVals = [3/nTrials];
cVals = round(cVals.*1000000);

% cVals = cVals(end-1);

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


%get file names
%%%%%%%
% can probably load in all variables wanted like this
% still to add - calc how many files to load then set up size of fnames
%%%%%%

fnames={};
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
    for iF = 1:length(f)
        filesToLoad{iF} = f(iF).name;
        fnames = [fnames filesToLoad{iF}];
    end
end


for iF = 1:length(fnames)
    load(fnames{iF});

    for iterI = 1:nIter
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
            densityPlot(:,:,iterI) = sum(densityPlotClus,3);
            densityPlot(:,:,iterI) = imgaussfilt(densityPlot(:,:,iterI),gaussSmooth); %smooth
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
        tsse(iterI)=sum(sse);
        
        %compute 'spreaded-ness' - variance of SE across clusters is a measure
        %of this, assuing uniform data points
        devAvgSSE(:,iterI)   = sse-mean(sse);
        stdAcrossClus(iterI) = std(devAvgSSE(:,iterI));
        varAcrossClus(iterI) = var(devAvgSSE(:,iterI));
        
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

        aCorrMap(:,:,iterI) = ndautoCORR(densityPlot(:,:,iterI));
        
        %plot
        if doPlot
            figure;
            subplot(2,2,1);
            imagesc(densityPlot(:,:,iterI));
            subplot(2,2,2);
        end
        [g,gdataA] = gridSCORE(aCorrMap(:,:,iterI),'allen',doPlot);
        if doPlot, subplot(2,2,3); end
        [g,gdataW] = gridSCORE(aCorrMap(:,:,iterI),'wills',doPlot);
        
        
        gA_g(iF,iterI)   = gdataA.g_score;
        gA_o(iF,iterI)   = gdataA.orientation;
        gA_wav(iF,iterI) = gdataA.wavelength;
        gA_rad(iF,iterI) = gdataA.radius;
        
    end
    gA{iF}=gdataA;
    gW{iF}=gdataW;
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



%compare with bar plots, etc.


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


