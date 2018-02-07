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

% load in gridness value computed on: simple, activation of gauss, or normalised version of activation of gauss
densityPlotMsr = 'actNorm';  % 'simple','act','actNorm'

doPlot = 0; %plot - note if loading in many value, this will make too many plots, can crash

gaussSmooth=1; %smooth maps by x value

%params to compare
epsMuVals = [0.0005, 0.0008, 0.001, 0.0015, 0.002, .003, .0035, .004].*10000; 

stepSize = 1;
% sigmaGaussVals = round([stepSize/3.5, stepSize/4].*100);
sigmaGaussVals = round([stepSize/3, stepSize/3.5, stepSize/4].*100);

alphaVals = [0, 2, 5, 8];

stochasticType = 0; %1, 2, 3; % 0, 1 - 2 and 3 are too stoch
cVals = ([1/nTrials, 5/nTrials, 10/nTrials, 20/nTrials]); %.1 and stoch = towards middle % cVals = [20/nTrials];
% cVals = 20/nTrials;

cVals = round(cVals.*1000000);
if ~stochasticType, cVals = 0; end

nEps   = length(epsMuVals);
nAlpha = length(alphaVals);
nC     = length(cVals);
nSigma = length(sigmaGaussVals);

nClus = 10; %10, 20

nTrials = 40000;
nIter = 200;
%warpBox=0;

nSet=8;

spacing = linspace(locRange(1),locRange(2),locRange(2)+1);
% tsse=nan(nIter,nC);
% devAvgSSE=nan(nClus,nIter,nEps,nAlpha,nC,nSigma);
densityPlotAll = zeros(length(spacing),length(spacing),nSet,nIter,nEps,nAlpha,nC,nSigma);
aCorrMap = nan(length(spacing)*2-1,length(spacing)*2-1,nSet,nIter,nEps,nAlpha,nC,nSigma);

for iEps = 1:nEps
    for iAlpha = 1:nAlpha
        fprintf('iEps = %d, iAlpha = %d\n',iEps, iAlpha);
        for iC = 1:nC
            fprintf('iC = %d\n',iC);
            for iSigma = 1:nSigma
                fprintf('sigma = %d\n',iSigma);

                %set params, load in data
                %learning rate
                epsMuOrig10000 = epsMuVals(iEps);
                %alpha - momentum parameter
                alpha10 = alphaVals(iAlpha);
                %stochastic parameters
                c =cVals(iC);
                
                %std of gaussian function
                sigmaGauss100  = sigmaGaussVals(iSigma);
                
                %load
                fname = [saveDir, sprintf('/covering_map_dat_gauss_%dclus_%dsigma_%dtrls_eps%d_alpha%d_stype%d_cVal%d_%diters*',nClus,sigmaGauss100,nTrials,epsMuOrig10000,alpha10,stochasticType,c,nIter)];
                
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
                        densityPlotAll(:,:,iSet,iterI,iEps,iAlpha,iC,iSigma) = imgaussfilt(densityPlot(:,:,iSet,iterI),gaussSmooth);
                    end
                end
                
                % load in gridness value computed on: simple, activation of
                % gauss, or normalised version of activation of gauss
                switch densityPlotMsr
                    case 'simple'
                        gAtmp = gA;
                        gWtmp = gW;
                    case 'act'
                        gAtmp = gA;
                        gWtmp = gW;
                    case 'actNorm'
                        gAtmp = gA;
                        gWtmp = gW;
                end
                %organise gridness values (allen vs willis method)
                gA_gAll(:,:,iEps,iAlpha,iC,iSigma)   = gAtmp(:,:,1);
                gA_oAll(:,:,iEps,iAlpha,iC,iSigma)   = gAtmp(:,:,2);
                gA_radAll(:,:,iEps,iAlpha,iC,iSigma) = gAtmp(:,:,3);
                gA_wavAll(:,:,iEps,iAlpha,iC) = gAtmp(:,:,4);
                
                gW_gAll(:,:,iEps,iAlpha,iC,iSigma) = gWtmp(:,:,1);
                gW_oAll(:,:,iEps,iAlpha,iC,iSigma) = gWtmp(:,:,2);
                gW_radAll(:,:,iEps,iAlpha,iC,iSigma) = gWtmp(:,:,3);
                gW_wavAll(:,:,iEps,iAlpha,iC,iSigma) = gWtmp(:,:,4);

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
        

for iSet = 3%:8

% comparing learning rate, with alpha in subplots; with sigmaGauss in
% separate figures
for iSigma = 1:nSigma
    figure; hold on;
    for iAlpha = 1:nAlpha
        for iC=1:nC
            subplot(1,nAlpha,iAlpha);
            dat=squeeze(datTmp(iSet,:,:,iAlpha,iC,iSigma));
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
            title(sprintf('%s - alpha=%d, c=%d, sig=%d',gridMeasure,alphaVals(iAlpha),cVals(iC),sigmaGaussVals(iSigma)))
        end
    end
end


% comparing alpha, with learning rate in subplots, with sigmaGauss in
% separate figures
for iSigma = 1:nSigma
    figure; hold on;
    for iEps=1:nEps
        subplot(1,nEps,iEps);
        dat=squeeze(datTmp(iSet,:,iEps,:,:,iSigma));
        barpos = .25:.5:.5*size(dat,2);
        colors = distinguishable_colors(size(dat,2));
        mu=mean(dat,1);
        sm=std(dat)./sqrt(size(dat,1));
        ci=sm.*tinv(.025,size(dat,1)-1); %compute conf intervals
        % figure; hold on;
        plotSpread(dat,'xValues',barpos,'distributionColors',colors);
        errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
        scatter(barpos,mu,750,colors,'.');
        xlim([barpos(1)-.5, barpos(end)+.5]);
        if strcmp(gridMeasure,'grid'), ylim([-.5, 1]); end
        title(sprintf('%s - eps=%d, c=%d, sig=%d',gridMeasure,epsMuVals(iEps),cVals(iC),sigmaGaussVals(iSigma)))
    end
end

end



%%
for iSet = 1%:8

% %remove nans - this is a bit hacky ...
% if any(isnan(gA_g(1,:)))
%     ind=isnan(dat(1,:)); gA_g(:,ind)=[]; gA_o(:,ind)=[];gA_rad(:,ind)=[];gA_wav(:,ind)=[]; barpos(ind)=[]; colors(ind,:)=[];
% end
% if length(barpos)==length(ind); barpos(ind)=[]; colors(ind,:)=[]; end %this is reset at top..

% comparing learning rate, with stochastic cVals in subplots, with alpha
% vals per fig; with sigmaGauss doubling the number of figures
for iSigma = 1:nSigma
    for iAlpha = 1:nAlpha
        figure; hold on;
        for iC=1:nC
            subplot(1,nC,iC);
            dat=squeeze(datTmp(iSet,:,:,iAlpha,iC,iSigma));
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
            title(sprintf('%s - alpha=%d, c=%d, sig=%d',gridMeasure,alphaVals(iAlpha),cVals(iC),sigmaGaussVals(iSigma)))
        end
    end
end


% comparing stochastic cVals, with learning rate in subplots, with alpha vals per fig; with sigmaGauss doubling the number of figures
for iSigma = 1:nSigma
    for iAlpha = 1:nAlpha
        figure; hold on;
        for iEps=1:nEps
            subplot(1,nEps,iEps);
            dat=squeeze(datTmp(iSet,:,iEps,iAlpha,:,iSigma));
            barpos = .25:.5:.5*size(dat,2);
            colors = distinguishable_colors(size(dat,2));
            mu=mean(dat,1);
            sm=std(dat)./sqrt(size(dat,1));
            ci=sm.*tinv(.025,size(dat,1)-1); %compute conf intervals
            % figure; hold on;
            plotSpread(dat,'xValues',barpos,'distributionColors',colors);
            errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
            scatter(barpos,mu,750,colors,'.');
            xlim([barpos(1)-.5, barpos(end)+.5]);
            if strcmp(gridMeasure,'grid'), ylim([-.5, 1]); end
            title(sprintf('%s - eps=%d, alpha=%d, sig=%d',gridMeasure,epsMuVals(iEps),alphaVals(iAlpha),sigmaGaussVals(iSigma)))
        end
    end
end

% % comparing stochastic cVals, with alpha in subplots, with learning rate vals per fig
% for iEps=1:nEps
%     figure; hold on;
%     for iAlpha = 1:nAlpha
%         subplot(1,4,iAlpha);
%         
%         switch gridMeasure
%             case 'grid'
%                 datTmp=gridness;
%             case 'angle'
%                 datTmp=orientation;
%             case 'rad'
%                 datTmp=rad;
%             case 'wav'
%                 datTmp=wav;
%         end
%         dat=squeeze(datTmp(iSet,:,:,iAlpha,iC));
%         barpos = .25:.5:.5*size(dat,2);
%         colors = distinguishable_colors(size(dat,2));
%         
%         mu=mean(dat,1);
%         sm=std(dat)./sqrt(size(dat,1));
%         ci=sm.*tinv(.025,size(dat,1)-1); %compute conf intervals
%         % figure; hold on;
%         plotSpread(dat,'xValues',barpos,'distributionColors',colors);
%         errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
%         scatter(barpos,mu,750,colors,'.');
%         xlim([barpos(1)-.5, barpos(end)+.5]);
%         if strcmp(gridMeasure,'grid'), ylim([-.5, 1]); end
%         title(sprintf('%s - eps=%d, alpha=%d',gridMeasure,epsMuVals(iEps),alphaVals(iAlpha)))
%     end
% end





end

















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

%% compute sse with respect to all test data points


