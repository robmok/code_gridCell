clear all;

% wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
cd(wd);

codeDir = [wd '/code_gridCell'];
saveDir = [wd '/data_gridCell'];
addpath(codeDir); addpath(saveDir);
addpath(genpath([codeDir '/gridSCORE_packed']));

% load
nSet        = 22;
gaussSmooth = 1; 
fixBatchSize = 1; %fixed batch size or depend on nClus (for fname)

dat='circ';
dat='square';
annEps=0;
boxSize=1;

% joined trials
jointTrls=1;
% clus2run = [8, 12, 16, 20,24, 28]; 
% clus2run = [8:2:28]; 
epsMuVals=.025;
nTrials=1000000;
% batchSizeVals = [1000, 400, 100]; 
% batchSizeVals=400;
annEps=0;
nIter=200;

% clus2run = [3:16, 18, 20:26];  %missed 17, 19?


dat='trapzKrupic';

clus2run = [3:26]; 
% clus2run=3:26;
% clus2run=26;
clus2run=4:2:26;

batchSizeVals = 400; %100, 125, 200,400, 1000

%new - slower learning rate
% epsMuVals=.015;
% batchSizeVals = 100; %100, 125, 200, 400
% clus2run = [12, 16, 24, 28]; %batchSize200 missed 20?

rHex=0; %if choose raw 60deg corr values, not gridness


%perm
nPerm=500;
nIters2run=nIter;

%load loop
for iClus2run = 1:length(clus2run) 
    nClus = clus2run(iClus2run);
    for iEps = 1:length(epsMuVals) 
        epsMuOrig=epsMuVals(iEps);
        epsMuOrig1000=epsMuOrig*1000;
        for iBvals = 1:length(batchSizeVals)
            batchSize = batchSizeVals(iBvals);
            fprintf('Loading %s, nClus=%d, epsMu=%d, batchSize=%d\n',dat,nClus,epsMuOrig1000,batchSize)
            if ~strcmp(dat(1:4),'trap')
                fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wAct_jointTrls_stepSiz_perm_%dpermsOn%diters',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat,nPerm,nIters2run)];
            else
                fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wAct_jointTrls_stepSiz_noPerm',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat)];
%                 fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wAct_jointTrls_stepSiz_noPerm_noJointTrlsTest',nClus,round(nTrials/1000),epsMuOrig1000,batchSize,nIter,dat)];
            end
                        
%             if boxSize>1
%                 fname = [fname sprintf('_boxSizex%d',boxSize)];
%             end
%             if annEps %epsMu is different here
%                 fname = [sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps*_batchSiz%d_%diters_%s_wAct_jointTrls_stepSiz_annEps',nClus,round(nTrials/1000),batchSize,nIter,dat)];
%             end            
%             finish with directory and * for date/time

            fname = [saveDir, fname '*'];
            
            %edit if want to load more than one file per sim, merge
            f = dir(fname); filesToLoad = cell(1,length(f));
            for iF = 1%:length(f)
                filesToLoad{iF} = f(iF).name;
                load(f(iF).name);
            end
            
            %organise gridness values (allen vs willis method)
            gA_gAll_act(:,iEps,iBvals,iClus2run,:)   = gA_act(:,1,:);
            gA_oAll_act(:,iEps,iBvals,iClus2run,:)   = gA_act(:,2,:);
            gA_radAll_act(:,iEps,iBvals,iClus2run,:) = gA_act(:,3,:);
            gA_wavAll_act(:,iEps,iBvals,iClus2run,:) = gA_act(:,4,:);
            gW_gAll_act(:,iEps,iBvals,iClus2run,:) = gW_act(:,1,:);
            gW_oAll_act(:,iEps,iBvals,iClus2run,:) = gW_act(:,2,:);
            gW_radAll_act(:,iEps,iBvals,iClus2run,:) = gW_act(:,3,:);
            gW_wavAll_act(:,iEps,iBvals,iClus2run,:) = gW_act(:,4,:);
            %
            gA_gAll_actNorm(:,iEps,iBvals,iClus2run,:)   = gA_actNorm(:,1,:);
            gA_oAll_actNorm(:,iEps,iBvals,iClus2run,:)   = gA_actNorm(:,2,:);
            gA_radAll_actNorm(:,iEps,iBvals,iClus2run,:) = gA_actNorm(:,3,:);
            gA_wavAll_actNorm(:,iEps,iBvals,iClus2run,:) = gA_actNorm(:,4,:);
            gW_gAll_actNorm(:,iEps,iBvals,iClus2run,:) = gW_actNorm(:,1,:);
            gW_oAll_actNorm(:,iEps,iBvals,iClus2run,:) = gW_actNorm(:,2,:);
            gW_radAll_actNorm(:,iEps,iBvals,iClus2run,:) = gW_actNorm(:,3,:);
            gW_wavAll_actNorm(:,iEps,iBvals,iClus2run,:) = gW_actNorm(:,4,:);  
            
            %note this perm is only with actNorm - forgot to save for act
            if ~strcmp(dat(1:4),'trap')
            gA_act_permPrc(:,iEps,iBvals,iClus2run,:)   = permPrc_gA_act(:,3,:);
%             gW_act_permPrc(:,iEps,iBvals,iClus2run,:)   = permPrc_gW_act(:,3,:);
            
            gA_actNorm_permPrc(:,iEps,iBvals,iClus2run,:)   = permPrc_gA_actNorm(:,1,:);
%             gW_actNorm_permPrc(:,iEps,iBvals,iClus2run,:)   = permPrc_gW_actNorm(:,1,:);
            end
            
        end 
    end
end
%% Making figs - univar scatters 1 - test set

figsDir = [wd '/grid_figs'];

savePlots=0;

clusPosAct = 'act'; %'act' or 'actNorm'

gridMsrType = 'a'; % 'a' or 'w' for allen or willis method - a preferred

gridMeasure = 'grid';

switch clusPosAct
% case 'clus'
%     switch gridMsrType
%         case 'a'
%             gridness    = gA_gAll;
%             orientation = gA_oAll;
%             rad         = gA_radAll;
%             wav         = gA_wavAll;
%         case 'w'
%             gridness    = gW_gAll;
%             orientation = gW_oAll;
%             rad         = gW_radAll;
%             wav         = gW_wavAll;
%     end
    case 'act'
        switch gridMsrType
            case 'a'
                gridness    = gA_gAll_act;
                orientation = gA_oAll_act;
                rad         = gA_radAll_act;
                wav         = gA_wavAll_act;
            case 'w'
                gridness    = gW_gAll_act;
                orientation = gW_oAll_act;
                rad         = gW_radAll_act;
                wav         = gW_wavAll_act;
        end
    case 'actNorm'
        switch gridMsrType
            case 'a'
                gridness    = gA_gAll_actNorm;
                orientation = gA_oAll_actNorm;
                rad         = gA_radAll_actNorm;
                wav         = gA_wavAll_actNorm;
            case 'w'
                gridness    = gW_gAll_actNorm;
                orientation = gW_oAll_actNorm;
                rad         = gW_radAll_actNorm;
                wav         = gW_wavAll_actNorm;
        end
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


iEps=1;

% clus2plot=(3:26)-2;
clus2plot=(6:26)-2;
% clus2plot=([6:24])-2;

clus2plot=1:length(clus2run);

%trapz
% clus2plot = 1:length(clus2run);

iBatchVals=1;

%fig specs
xTickLabs = num2cell(clus2run(clus2plot));
fontSiz=15;

figure; hold on;
    for iEps = 1:length(epsMuVals)
        %     subplot(2,3,iEps);
        dat1     = squeeze(datTmp(:,iEps,iBatchVals,clus2plot,1));
        barpos  = .25:.5:.5*size(dat1,2);
        colors  = distinguishable_colors(size(dat1,2));
        colgrey = [.6, .6, .6];
        mu      = mean(dat1,1);
        sm      = std(dat1)./sqrt(size(dat1,1));
        ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
        plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
%         errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
%         scatter(barpos,mu,50,colors,'filled','d');
        scatter(barpos,mu,50,colgrey,'filled','d');
        xticklabels(xTickLabs);
        xlim([barpos(1)-.5, barpos(end)+.5]);
        %         ylim([0,1]);
        if strcmp(gridMsrType,'a')
            ylim([-.45,1.4]);
        elseif strcmp(gridMsrType,'w')
            ylim([-1.25,1.4]);
        end
        xlabel('Number of Clusters');
        ylabel('Grid Score');
        
%         title(sprintf('%s, %s - eps=%d, batchSize=%d',dat, gridMeasure,epsMuVals(iEps)*1000,batchSizeVals(iBatchVals)))
        if strcmp(dat(1:2),'ci')
        title('Circular box')
        elseif strcmp(dat(1:2),'sq')
        title('Square box')
        end        
        
        set(gca,'FontSize',fontSiz,'fontname','Arial')
    end

    fname = [figsDir sprintf('/gridness_%s_univarScatters_testSet_nClus%d-%d_eps%d_batchSiz%d_%s_%s',dat,clus2run(clus2plot(1)),clus2run(clus2plot(end)),epsMuVals(iEps)*1000,batchSizeVals(iBatchVals),clusPosAct,gridMsrType)];

if savePlots
   set(gcf,'Renderer','painters');
   print(gcf,'-depsc2',fname)
   saveas(gcf,fname,'png');
end


% %prop grid cells, averaged across clusters
if ~strcmp(dat(1:4),'trap')
    maxThresh=squeeze(max(max(gA_act_permPrc(:,1,1,:))));
propGrid=nan(1,size(dat1,2));
for i=1:size(dat1,2) %length(clus2run)
propGrid(i)=nnz(dat1(:,i)>maxThresh)/200;
end
mean(propGrid)
std(propGrid)/sqrt(length(propGrid))
end
% 
% 
% %plot 3:5
% clus2plot=(3:5)-2;
% %fig specs
% xTickLabs = num2cell(clus2run(clus2plot));
% fontSiz=15;
% figure; hold on;
%     for iEps = 1:length(epsMuVals)
%         %     subplot(2,3,iEps);
%         dat1     = squeeze(datTmp(:,iEps,iBatchVals,clus2plot,1));
%         barpos  = .25:.5:.5*size(dat1,2);
%         colors  = distinguishable_colors(size(dat1,2));
%         colgrey = [.5, .5, .5];
%         mu      = mean(dat1,1);
%         sm      = std(dat1)./sqrt(size(dat1,1));
%         ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
%         plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
%         errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
%         scatter(barpos,mu,50,colors,'filled','d');
%         xticklabels(xTickLabs);
%         xlim([barpos(1)-.5, barpos(end)+.5]);
%         %         ylim([0,1]);
%         
%         if strcmp(gridMsrType,'a')
%             ylim([-.45,1.4]);
%         elseif strcmp(gridMsrType,'w')
%             ylim([-1.25,1.4]);
%         end
%         xlabel('Number of Clusters');
%         ylabel('Grid Score');
%         
%         title(sprintf('%s, %s - eps=%d, batchSize=%d',dat, gridMeasure,epsMuVals(iEps)*1000,batchSizeVals(iBatchVals)))
%         set(gca,'FontSize',fontSiz,'fontname','Arial')
%     end
% 
%     fname = [figsDir sprintf('/gridness_%s_univarScatters_nClus%d-%d_eps%d_batchSiz%d_%s',dat,clus2run(clus2plot(1)),clus2run(clus2plot(end)),epsMuVals(iEps)*1000,batchSizeVals(iBatchVals),gridMsrType)];
% if savePlots
%    set(gcf,'Renderer','painters');
%    print(gcf,'-depsc2',fname)
%    saveas(gcf,fname,'png');
% end

%%

if strcmp(dat(1:4),'trap')
    figure; hold on;
    dat1    = squeeze(datTmp(:,:,:,:,2));
    mu      = nanmean(dat1,1);
    sm      = nanstd(dat1)./sqrt(size(dat1,1));
    ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
    plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
    errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
    scatter(barpos,mu,100,colors,'d','filled');
    xlim([barpos(1)-.5, barpos(end)+.5]);
    ylim([-1.5,1.5]);
    title(sprintf('Left half of box %s - eps=%d',gridMeasure,epsMuVals(iEps)*1000))
    
    figure; hold on;
    dat1    = squeeze(datTmp(:,:,:,:,3));
    mu      = nanmean(dat1,1);
    sm      = nanstd(dat1)./sqrt(size(dat1,1));
    ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
    plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
    errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
    scatter(barpos,mu,100,colors,'d','filled');
    xlim([barpos(1)-.5, barpos(end)+.5]);
    ylim([-1.5,1.5]);
    title(sprintf('Right half of box %s - eps=%d',gridMeasure,epsMuVals(iEps)*1000))
    
    figure; hold on;
    dat1    = squeeze(datTmp(:,:,:,:,2))-squeeze(datTmp(:,:,:,:,3));
    mu      = nanmean(dat1,1);
    sm      = nanstd(dat1)./sqrt(size(dat1,1));
    ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
    plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
    errorbar(barpos,mu,ci,'Color',colgrey,'LineStyle','None','LineWidth',1);
    scatter(barpos,mu,100,colors,'d','filled');
    xlim([barpos(1)-.5, barpos(end)+.5]);
    ylim([-1.5,2]);
    title(sprintf('Left-right half of box %s - eps=%d',gridMeasure,epsMuVals(iEps)*1000))
    
    
end
%% Thresholds from perm stats
figsDir = [wd '/grid_figs'];

savePlots=0;

clusPosAct = 'actNorm'; %'act' or 'actNorm'

gridMsrType = 'a'; % 'a' or 'w' for allen or willis method - a preferred

gridMeasure = 'grid';

switch clusPosAct
    case 'act'
        switch gridMsrType
            case 'a'
                permThresh = gA_act_permPrc;
            case 'w'
                permThresh = gW_act_permPrc;
        end
    case 'actNorm'
        switch gridMsrType
            case 'a'
                %ATM WRONG NAMES
                
                permThresh = gA_act_permPrc;
%                 permThresh = gA_actNorm_permPrc;
            case 'w'
                permThresh = gW_act_permPrc;
%                 permThresh = gW_actNorm_permPrc;
        end
end

% clus2plot=([6:26])-2;
clus2plot=([6:24])-2;

iBatchVals=1; %'medium' one

%fig specs
xTickLabs = num2cell(clus2run(clus2plot));
fontSiz=15;

% gA_act_permPrc(:,1,1,:)

dSize=400;

figure; hold on;
    for iEps = 1:length(epsMuVals)
        %     subplot(2,3,iEps);
        dat1     = squeeze(permThresh(:,iEps,iBatchVals,clus2plot,1));
        barpos  = .25:.5:.5*size(dat1,2);
        colors  = distinguishable_colors(size(dat1,2));
        colgrey = [.5, .5, .5];
        mu      = mean(dat1,1);
        sm      = std(dat1)./sqrt(size(dat1,1));
        ci      = sm.*tinv(.025,size(dat1,1)-1); %compute conf intervals
        plotSpread(dat1,'xValues',barpos,'distributionColors',colors);
        scatter(barpos,min(dat1),dSize,colors,'d');  %min
        scatter(barpos,mean(dat1),dSize,colors,'d'); %mean
        scatter(barpos,max(dat1),dSize,colors,'d');  %max
        xticklabels(xTickLabs);
        xlim([barpos(1)-.5, barpos(end)+.5]);
%         if strcmp(gridMsrType,'a')
%             ylim([-.45,1.4]);
%         elseif strcmp(gridMsrType,'w')
%             ylim([-1.25,1.4]);
%         end
        xlabel('Number of Clusters');
        ylabel('Grid Score Threshold (95% of permuted distribution)');
        
        
%         title(sprintf('%s, %s - eps=%d, batchSize=%d',dat, gridMeasure,epsMuVals(iEps)*1000,batchSizeVals(iBatchVals)))
        
        set(gca,'FontSize',fontSiz,'fontname','Arial')
    end

    fname = [figsDir sprintf('/gridness_%s_univarScatters_permThresh_nClus%d-%d_eps%d_batchSiz%d_%s',dat,clus2run(clus2plot(1)),clus2run(clus2plot(end)),epsMuVals(iEps)*1000,batchSizeVals(iBatchVals),gridMsrType)];
if savePlots
   set(gcf,'Renderer','painters');
   print(gcf,'-depsc2',fname)
   saveas(gcf,fname,'png');
end

   


% squeeze(min(min(permThresh(:,1,1,:))))
% squeeze(mean(mean(permThresh(:,1,1,:))))
% squeeze(max(max(permThresh(:,1,1,:))))
    
    