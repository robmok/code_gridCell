clear all;

% wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
cd(wd);

codeDir = [wd '/code_gridCell'];
saveDir = [wd '/data_gridCell'];
addpath(codeDir); addpath(saveDir);
addpath(genpath([wd '/gridSCORE_packed']));

%make test data - keep same for comparison
nSteps = 50;
locRange = [0, nSteps-1]; 
nTrials = 40000;
nTrialsTest = nTrials; % nTrialsTest = 50000;
dataPtsTest = [randsample(linspace(locRange(1),locRange(2),nSteps),nTrialsTest,'true'); randsample(linspace(locRange(1),locRange(2),nSteps),nTrialsTest,'true')]'; % random points in a box

%% Compute the cluster centres from the desnity map, then compute SSE and rank them

dat = 'cov'; % cov or gauss - covering map or gauss data

doPlot = 0; %plot - note if loading in many value, this will make too many plots, can crash

gaussSmooth=1; %smooth maps by x value

nClus = 30; %10, 20 (30)


%params to compare
epsMuVals = [25, 50 75 100];

switch dat
    case 'cov'
        epsMuVals = [.015, .025, .05, .075, .1, .2].*1000;
    case 'gauss'
        epsMuVals = [0.0025, 0.005, .01, .02]*10000;
end

alphaVals = 0;%[0, 2, 5, 8];

stochasticType = 0; %1, 2, 3; % 0, 1 - 2 and 3 are too stoch
% cVals = ([2/nTrials, 5/nTrials, 10/nTrials, 20/nTrials]); %.1 and stoch = towards middle % cVals = [20/nTrials];
cVals = 20/nTrials;
cVals = round(cVals.*1000000);
if ~stochasticType, cVals = 0; end

%%%%%%
sigmaGauss100 = round(1/3.5*100); %need to add the other one if decide to test this - see gauss_compare_params
%%%%%% 

sseWvals= [0,1];
% weightEpsSSE


nEps   = length(epsMuVals);
nAlpha = length(alphaVals);
nC     = length(cVals);
nsseW  = length(sseWvals);

nTrials = 40000;
nIter = 100;
%warpBox=0;

nSet=4;

spacing = linspace(locRange(1),locRange(2),locRange(2)+1);
tsse=nan(nIter,nC,nsseW);
devAvgSSE=nan(nClus,nIter,nEps,nAlpha,nC,nsseW);
densityPlotAll = zeros(length(spacing),length(spacing),nSet,nIter,nEps,nAlpha,nC,nsseW);
aCorrMap = nan(length(spacing)*2-1,length(spacing)*2-1,nSet,nIter,nEps,nAlpha,nC,nsseW);

for iEps = 1:nEps
    for iAlpha = 1:nAlpha
        fprintf('iEps = %d, iAlpha = %d\n',iEps, iAlpha);
        for iC = 1:nC
            fprintf('iC = %d\n',iC);
            for isseW = 1:2
                fprintf('isseW = %d\n',isseW);
                
                %set params, load in data
                %learning rate
                epsMuOrig1000 = epsMuVals(iEps);
                %alpha - momentum parameter
                alpha10 = alphaVals(iAlpha);
                %stochastic parameters
                c =cVals(iC);
                
                %weight eps by sse
                weightEpsSSE = sseWvals(isseW);
                
                %load
                switch dat
                    case 'cov'
                        fname = [saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_alpha%d_stype%d_cVal%d_sseW%d_%diters*',nClus,nTrials,epsMuOrig1000,alpha10,stochasticType,c,weightEpsSSE,nIter)];
                    case 'gauss'
                        fname = [saveDir, sprintf('/covering_map_dat_gauss_%dclus_%dsigma_%dtrls_eps%d_alpha%d_stype%d_cVal%d_sseW%d_%diters*',nClus,sigmaGauss100,nTrials,epsMuOrig1000,alpha10,stochasticType,c,weightEpsSSE,nIter)];
                end
                
                %if ran more than 1 set of iters, hv multiple files w diff end bit
                %for now, only load 1st file
                f = dir(fname); filesToLoad = cell(1,length(f));
                %             muAllTmp={}; nIterCount=0;
                for iF = 1%:length(f)
                    filesToLoad{iF} = f(iF).name;
                    % fnames = [fnames filesToLoad{iF}];
                    load(f(iF).name);
                    % muAllTmp{iF}=muAll;
                    % nIterCount = [nIterCount, nIter+nIterCount(end)]; %count number of iters to index below when merging
                end
                %           %merge... % if want to do more iters and merge
                
                for iterI = 1:nIter
                    for iSet=1:nSet
                        densityPlotAll(:,:,iSet,iterI,iEps,iAlpha,iC) = imgaussfilt(densityPlot(:,:,iSet,iterI),gaussSmooth);
                    end
                end
                %organise gridness values (allen vs willis method)
                gA_gAll(:,:,iEps,iAlpha,iC,isseW)   = gA(:,:,1);
                gA_oAll(:,:,iEps,iAlpha,iC,isseW)   = gA(:,:,2);
                gA_radAll(:,:,iEps,iAlpha,iC,isseW) = gA(:,:,3);
                gA_wavAll(:,:,iEps,iAlpha,iC,isseW) = gA(:,:,4);
                
                gW_gAll(:,:,iEps,iAlpha,iC,isseW) = gW(:,:,1);
                gW_oAll(:,:,iEps,iAlpha,iC,isseW) = gW(:,:,2);
                gW_radAll(:,:,iEps,iAlpha,iC,isseW) = gW(:,:,3);
                gW_wavAll(:,:,iEps,iAlpha,iC,isseW) = gW(:,:,4);
            end
        end
    end
end
%% plot to compare
% the way the files are set up & loaded in loses info about which params.
% however can recover it by using the same loops when loading in the param
% names? e.g. make a variable with the param names, then plot them

barw   = .25;
colgrey = [.25 .25 .25];

gridMeasure = 'grid'; % grid, angle (orientation), rad (radius), wav (wavelength)
gridMsrType = 'a'; % a/w - allen or willis


switch gridMsrType
    case 'a'
        gridness    = gA_gAll;
        orientation = gA_oAll;
        rad         = gA_radAll;
        wav         = gA_wavAll;
    case 'w'
        gridness    = gW_gAll;
        orientation = gW_oAll;
        rad         = gW_radAll;
        wav         = gW_wavAll;
end

switch gridMeasure
    case 'grid'
        datTmp=gridness;
    case 'angle'
        datTmp=orientation;
    case 'rad'
        datTmp=rad;
    case 'wav'
        datTmp=wav;
end


for iSet = 1:4
    
    figure; hold on;
    for isseW = 1:nsseW
        for iAlpha=1:nAlpha
            for iC=1:nC
                subplot(1,2,isseW);
                dat=squeeze(datTmp(iSet,:,:,iAlpha,iC,isseW));
                barpos = .25:.5:.5*size(dat,2);
                colors = distinguishable_colors(size(dat,2));
                mu=mean(dat,1);
                sm=std(dat)./sqrt(size(dat,1));
                ci=sm.*tinv(.025,size(dat,1)-1); %compute conf intervals
                plotSpread(dat,'xValues',barpos,'distributionColors',colors);
                errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
                scatter(barpos,mu,750,colors,'.');
                xlim([barpos(1)-.5, barpos(end)+.5]);
                if strcmp(gridMeasure,'grid'), ylim([-.5, 1]); end
                title(sprintf('%s - alpha=%d, c=%d, sseW=%d',gridMeasure,alphaVals(iAlpha),cVals(iC),sseWvals(isseW)))
            end
        end
    end
    
    
    
    % comparing learning rate, with alpha vals in subplots
%     figure; hold on;
%     for iAlpha = 1:nAlpha
%         for iC=1:nC
%             subplot(1,4,iAlpha);
%             dat=squeeze(datTmp(iSet,:,:,iAlpha,iC));
%             barpos = .25:.5:.5*size(dat,2);
%             colors = distinguishable_colors(size(dat,2));
%             mu=mean(dat,1);
%             sm=std(dat)./sqrt(size(dat,1));
%             ci=sm.*tinv(.025,size(dat,1)-1); %compute conf intervals
%             plotSpread(dat,'xValues',barpos,'distributionColors',colors);
%             errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
%             scatter(barpos,mu,750,colors,'.');
%             xlim([barpos(1)-.5, barpos(end)+.5]);
%             if strcmp(gridMeasure,'grid'), ylim([-.5, 1]); end
%             %         title(sprintf('%s - eps=%d, alpha=%d, c=%d',gridMeasure,epsMuVals(iEps),alphaVals(iAlpha),cVals(iC)))
%             title(sprintf('%s - alpha=%d, c=%d',gridMeasure,alphaVals(iAlpha),cVals(iC)))
%         end
%     end
    
%     % comparing alpha vals, with learning rate in subplots
%     figure; hold on;
%     for iAlpha = 1:nAlpha
%         for iEps=1:nEps
%             subplot(1,4,iEps);
%             dat=squeeze(datTmp(iSet,:,iEps,:,iC));
%             barpos = .25:.5:.5*size(dat,2);
%             colors = distinguishable_colors(size(dat,2));
%             mu=mean(dat,1);
%             sm=std(dat)./sqrt(size(dat,1));
%             ci=sm.*tinv(.025,size(dat,1)-1); %compute conf intervals
%             % figure; hold on;
%             plotSpread(dat,'xValues',barpos,'distributionColors',colors);
%             errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
%             scatter(barpos,mu,750,colors,'.');
%             xlim([barpos(1)-.5, barpos(end)+.5]);
%             if strcmp(gridMeasure,'grid'), ylim([-.5, 1]); end
%             title(sprintf('%s - eps=%d',gridMeasure,epsMuVals(iEps)))
%         end
%     end
end



%%

%%%%%%%%%%%%
% ++ more ways of plotting / comparing? maybe add the mean gridness score
% to the plot?
%%%%%%%%%%%%













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