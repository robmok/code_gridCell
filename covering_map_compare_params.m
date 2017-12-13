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

doPlot = 0; %plot - note if loading in many value, this will make too many plots, can crash

gaussSmooth=1; %smooth maps by x value

% %for computing density map
% nTrlsToUse = 10000; %
% toTrlN     = nTrials; %-10000; %nTrials if from nTrlstoUse to end,
% fromTrlI   = toTrlN-nTrlsToUse+1;

%params to compare
epsMuVals = [50, 75, 100]; %atm, lowest learning rate looks best
% epsMuOrig1000 = 75; %50, 75, 100 - atm, lowest learning rate looks best

alphaVals = [2, 5, 8];

stochasticType = 0; %1, 2, 3; % stochasticVals=[1,2,3] - later
cVals = ([.25/nTrials, .5/nTrials, 2/nTrials, 3/nTrials, 5/nTrials, 10/nTrials, 20/nTrials]); %.1 and stoch = towards middle % cVals = [20/nTrials];
cVals = round(cVals.*1000000);
if ~stochasticType, cVals = 0; end

nEps   = length(epsMuVals);
nAlpha = length(alphaVals);
nC     = length(cVals);

nClus = 20; %20, 30

nTrials = 40000;
nIter = 500;
%warpBox=0;

nSet=10;

spacing = linspace(locRange(1),locRange(2),locRange(2)+1);
tsse=nan(nIter,nC);
devAvgSSE=nan(nClus,nIter,nEps,nAlpha,nC);
densityPlot = zeros(length(spacing),length(spacing),nSet,nIter,nEps,nAlpha,nC);
aCorrMap = nan(length(spacing)*2-1,length(spacing)*2-1,nSet,nIter,nEps,nAlpha,nC);

%- atm takes about 18s per condition (500 iters- compute density plot, cluster cetres, sse, acorr+gridness)
%- need to merge files, or merge the iterations from separate files

% better to compute de

for iEps = 1:nEps
    for iAlpha = 1:nAlpha
        fprintf('iEps = %d, iAlpha = %d\n',iEps, iAlpha);
        for iC = 1:nC
%             fprintf('iC = %d',iC);

            %     iPlot=0; %for subplot below
            
            %set params, load in data
            
            %learning rate
            epsMuOrig1000 = epsMuVals(iEps);
            
            %alpha - momentum parameter
            alpha10 = alphaVals(iAlpha);
            
            %stochastic parameters
            c =cVals(iC);
            
            fname = [saveDir, sprintf('/covering_map_dat_%dclus_%dtrls_eps%d_alpha%d_stype%d_cVal%d_%diters*',nClus,nTrials,epsMuOrig1000,alpha10,stochasticType,c,nIter)];
            
%             %if ran more than 1 set of iters, hv multiple files w diff end bit;load
            f = dir(fname); filesToLoad = cell(1,length(f));
%             muAllTmp={};
%             nIterCount=0;
            for iF = 1%:length(f)
                        filesToLoad{iF} = f(iF).name;
%                         fnames = [fnames filesToLoad{iF}];
                load(f(iF).name);
%                 muAllTmp{iF}=muAll;
%                 nIterCount = [nIterCount, nIter+nIterCount(end)]; %count number of iters to index below when merging
            end
%             %merge
%             clear muAll
%             for iF = 1:length(f)
%                 muAll(:,:,:,nIterCount(iF)+1:nIterCount(iF+1)) = muAllTmp{iF};
%             end
            
            %calculate            
            %load in densityPlotClus 
%             load(fname);
            
            clusMu=nan(nClus,2,nIter);
            for iterI = 1:nIter
                for iSet=1:10
                    
%                     fprintf('iter %d, set %d\n',iterI,iSet);
                    %                 %compute density map
                    %                 clus = round(muAll(:,:,fromTrlI:toTrlN,iterI));
                    %                 densityPlotClus     = zeros(length(spacing),length(spacing),nClus);
                    
%                     densityPlotClusSmth = zeros(length(spacing),length(spacing),nClus);
%                     for iClus=1:nClus
%                         %                     for iTrl=1:size(clus,3)
%                         %                         densityPlotClus(clus(iClus,1,iTrl),clus(iClus,2,iTrl),iClus) = densityPlotClus(clus(iClus,1,iTrl),clus(iClus,2,iTrl),iClus)+1;
%                         %                     end
%                         %find peaks
%                         densityPlotClusSmth(:,:,iClus)=imgaussfilt(densityPlotClus(:,:,iClus,iSet,iterI),gaussSmooth);
%                         [peakX, peakY] = find(densityPlotClusSmth(:,:,iClus)==max(max((densityPlotClusSmth(:,:,iClus)))));
%                         clusMu(iClus,:,iSet,iterI) = [peakX, peakY];
%                         
%                         %make combined (grid cell) plot, smooth
%                         densityPlot(:,:,iSet,iterI,iEps,iAlpha,iC) = imgaussfilt(sum(densityPlotClus(:,:,:,iSet,iterI),3),gaussSmooth);
%                     end

                    % compute sse with respect to all test data points
                    distTrl=[];
                    sse=nan(1,nClus);
                    for iClus = 1:size(clusMu,1)
                        distTrl(:,iClus)=sum([clusMu(iClus,1,iSet,iterI)-dataPtsTest(:,1), clusMu(iClus,2,iSet,iterI)-dataPtsTest(:,2)].^2,2);
                    end
                    [indValsTrl, indTmp]=min(distTrl,[],2); % find which clusters are points closest to
                    for iClus = 1:size(clusMu,1)
                        sse(iClus)=sum(sum([clusMu(iClus,1,iSet,iterI)-dataPtsTest(indTmp==iClus,1), clusMu(iClus,2,iSet,iterI)-dataPtsTest(indTmp==iClus,2)].^2,2)); %distance from each cluster from training set to datapoints closest to that cluster
                    end
                    tsse(iSet,iterI,iEps,iAlpha,iC)=sum(sse);
                    
                    %compute 'spreaded-ness' - variance of SE across clusters is a measure
                    %of this, assuing uniform data points
                    devAvgSSE            = sse-mean(sse);
                    stdAcrossClus(iSet,iterI,iEps,iAlpha,iC) = std(devAvgSSE);
                    varAcrossClus(iSet,iterI,iEps,iAlpha,iC) = var(devAvgSSE);
                    
                    %compute autocorrelation map - compute it at sim for
                    %gridness, but here compute again to visualise (not
                    %save for size)
%                     aCorrMap(:,:,iSet,iterI,iEps,iAlpha,iC) = ndautoCORR(densityPlot(:,:,iSet,iterI,iEps,iAlpha,iC));
                    
%                     [g,gdataA] = gridSCORE(aCorrMap(:,:,iSet,iterI,iEps,iAlpha,iC),'allen',0);
%                     [g,gdataW] = gridSCORE(aCorrMap(:,:,iSet,iterI,iEps,iAlpha,iC),'wills',0);
                    
                    gA_gA(iSet,iterI,iEps,iAlpha,iC)   = gdataA.g_score;
                    gA_oA(iSet,iterI,iEps,iAlpha,iC)   = gdataA.orientation;
                    gA_wavA(iSet,iterI,iEps,iAlpha,iC) = gdataA.wavelength;
                    gA_radA(iSet,iterI,iEps,iAlpha,iC) = gdataA.radius;
                    
                    gA_gW(iSet,iterI,iEps,iAlpha,iC)   = gdataW.g_score;
                    gA_oW(iSet,iterI,iEps,iAlpha,iC)   = gdataW.orientation;
                    gA_wavW(iSet,iterI,iEps,iAlpha,iC) = gdataW.wavelength;
                    gA_radW(iSet,iterI,iEps,iAlpha,iC) = gdataW.radius;
                end
                gA{iEps,iAlpha,iC}=gdataA;
                gW{iEps,iAlpha,iC}=gdataW;
            end
        end
    end
end

% % plot
% for iTestVal = 1:nCvals
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

barw   = .25;
colgrey = [.25 .25 .25];

gridMeasure = 'grid'; % grid, angle (orientation), rad (radius), wav (wavelength)
gridMsrType = 'a'; % a/w - allen or willis

%sets 
% 5k - 20k:25k, 25k:30k, 30k:35k, 35k:40k
% 10k - 20k:30k, 25:35k, 30k:40k,
% 15k - 20k:35k; 25k:40k
% 20k - 20k:40k

%in general, low learning rate best, 0.1 is bad; could do even slower

%- for some reason, the first ones are best..? averages are same, but much
%more at the top of the distribution around1... could be that theres more
%movement in early trials so average is more spread out?

for iSet = 1
    
switch gridMsrType
    case 'a'
        %temp - inspect how average over trials one at a time; or loop?
%         iSet=1;
        
        gA_g = squeeze(gA_gA(iSet,:,:,:,:));
        gA_o = squeeze(gA_oA(iSet,:,:,:,:));
        gA_rad = squeeze(gA_radA(iSet,:,:,:,:));
        gA_wav = squeeze(gA_wavA(iSet,:,:,:,:));
%         gA_g = gA_gA;
%         gA_o = gA_oA;
%         gA_rad = gA_radA;
%         gA_wav = gA_wavA;
    case 'w'
        gA_g = gA_gW;
        gA_o = gA_oW;
        gA_rad = gA_radW;
        gA_wav = gA_wavW;        
end


% %remove nans - this is a bit hacky now...
% if any(isnan(gA_g(1,:)))
%     ind=isnan(dat(1,:)); gA_g(:,ind)=[]; gA_o(:,ind)=[];gA_rad(:,ind)=[];gA_wav(:,ind)=[]; barpos(ind)=[]; colors(ind,:)=[];
% end
% if length(barpos)==length(ind); barpos(ind)=[]; colors(ind,:)=[]; end %this is reset at top..

% compare alpha and c for each learning rate
for iAlpha = 1:nAlpha %comparing cVals
    figure; hold on;    
    for iC=1:nC
%         subplot(1,7,iC);
% for iC=1:nC %comparing alpha vals
%     figure; hold on;
%     for iAlpha = 1:nAlpha
%         subplot(1,3,iAlpha);
%         
        switch gridMeasure
            case 'grid'
                datTmp=gA_g;
            case 'angle'
                datTmp=gA_o;
            case 'rad'
                datTmp=gA_rad;
            case 'wav'
                datTmp=gA_wav;
        end
        dat=datTmp(:,:,iAlpha,iC);
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
        %         title(sprintf('%s - eps=%d, alpha=%d, c=%d',gridMeasure,epsMuVals(iEps),alphaVals(iAlpha),cVals(iC)))
        title(sprintf('%s - alpha=%d, c=%d',gridMeasure,alphaVals(iAlpha),cVals(iC)))
    end
end

end


% for iEps=1:nEps
%     figure; hold on;
%     for iAlpha = 1:nAlpha
%         subplot(1,3,iAlpha);
%         
%         switch gridMeasure
%             case 'grid'
%                 datTmp=gA_g;
%             case 'angle'
%                 datTmp=gA_o;
%             case 'rad'
%                 datTmp=gA_rad;
%             case 'wav'
%                 datTmp=gA_wav;
%         end
%         dat=squeeze(datTmp(:,iEps,iAlpha,:));
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


%%
% % plot cluster centres, print their sse and stdAcrossClus
% for iterI=1:nIter
%     % figure('units','normalized','outerposition',[0 0 1 1]); hold on; iPlot=0;
%     figure; hold on; iPlot=0;
%     for iTestVal = 1:nC
%         iPlot = iPlot+1;
%         subplot(2,2,iPlot);
%         scatter(clusMu(:,1,iterI,iTestVal),clusMu(:,2,iterI,iTestVal),20e+2,colors,'.')
%         xlim(locRange); ylim(locRange);
%         title(sprintf('%d alpha, %.1f tsse, %.2f. spreadedness', testVals(iTestVal), tsse(iterI,iTestVal), stdAcrossClus(iterI,iTestVal)))
% %         title(sprintf('%d cVal, %.1f tsse, %.2f. spreadedness', testVals(iDataSet), tsse(iterI,iDataSet), stdAcrossClus(iterI,iDataSet)))
%     end
% end


