clear all;

% wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
cd(wd);

codeDir = [wd '/code_gridCell'];
savDir = [wd '/data_gridCell'];
figDir = [wd '/data_gridCell/figs'];
addpath(codeDir); addpath(savDir); addpath(figDir);
addpath(genpath([codeDir '/gridSCORE_packed']));

dat='catLearn';
% annEps=0;
boxSize=1;
nIter=50; %20, 50
locRange = [0, 49];

% clus2run = 2:26; 
% clus2run = [2,3,4,5,8]; 
% clus2run = [10, 15,20,25,30]; 

% clus2run = [2,3,4]; 
% clus2run = [5,8,10]; 
% clus2run = [15, 25, 30]; 

clus2run = 18;

jointTrls=0;
epsMuVals=.025;
nTrials=50000;
% batchSizeVals= [5, 10, 25]; %5, 10, 25?
batchSizeVals= 10; 

nCats=2; %2,3,4
stoch=0;
% cVals = [2, 4, 10, 40];
cVals = 0;

catsInfo.nCats=2; %2 categories
% sigmaG = [5 0; 0 5];   % isotropic % sigmaG = [1 .5; .5 2]; R = chol(sigmaG);  % non-isotropic
sigmaG = [3 0; 0 3];
catsInfo.R=chol(sigmaG);

catsInfo.msExample = 1;


%load loop
muAllClus = cell(1,length(clus2run));
rSeedAll  = cell(length(clus2run),length(batchSizeVals),length(cVals));

 for iClus = 1:length(clus2run)
    nClus = clus2run(iClus);
    for iEps = 1:length(epsMuVals)
        epsMuOrig=epsMuVals(iEps);
        epsMuOrig1000=epsMuOrig*1000;
        for iBvals = 1:length(batchSizeVals)
            for iC = 1:length(cVals)
                fprintf('Loading %s, nClus=%d, epsMu=%d, c=%d, batchSize=%d\n',dat,nClus,epsMuOrig1000,cVals(iC),batchSizeVals(iBvals))

%                 fname = [savDir, sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wAct_%dcats_stoch%d_c%d_*',nClus,round(nTrials/1000),epsMuOrig1000,batchSizeVals(iBvals),nIter,dat,nCats,stoch,cVals(iC))];
                fname = [savDir, sprintf('/covering_map_batch_dat_%dclus_%dktrls_eps%d_batchSiz%d_%diters_%s_wActNorm_%dcats_stoch%d_c%d_*',nClus,round(nTrials/1000),epsMuOrig1000,batchSizeVals(iBvals),nIter,dat,nCats,stoch,cVals(iC))];

                f=dir(fname);
                load(f.name);
                
                %just load in those to plot
                nBatches = size(muAll,3)-1;
%                 trls2Plt = [1, nBatches*.25, nBatches*.5, nBatches*.75, nBatches+1];
                trls2plt = [1, nBatches*.5, nBatches+1];
                muAllClus{iClus}(:,:,:,:,iBvals,iC)=muAll(:,:,trls2plt,:);
                rSeedAll{iClus,iBvals,iC} = rSeed;

%                 muAllClus{iClus}(:,:,:,:,iBvals,iC)=muAll;
%                 muAllClus{iClus}(:,:,:,:,iEps,iBvals,iC)=muAll;                
                
                
            end
        end
    end
end
    
%% 

%in the end used iter 13 and 14
savePlots = 0;

fontSiz=15;
datPtSiz=130;

% iterI=1;

trls2plt = {1:25, 1:100, 1:2000};

catCentres = [15, 35; 35, 15];

colGreyClus = [.85, .85, .85];
colGreyDat  = [.45 .45 .45];

catAcol = [0, 0, 1];
catBcol = [1, 0, 0];
catAcol = [catAcol+((1-catAcol).*.5); catAcol];
catBcol = [catBcol+((1-catBcol).*.5); catBcol];

% colors(:,:,2) = distinguishable_colors(nClus)+(1-distinguishable_colors(nClus)).*.5; %lighter


for iterI=1:20
    
    for iBvals= 1:length(batchSizeVals)
    for iC = 1:length(cVals)
    
    ctr=0;
%     figure; hold on;
    for iClus = 1:length(clus2run)       
%         colors    = distinguishable_colors(clus2run(iClus));
        %     figure; hold on;
        for iPlot = 2%1:length(trls2plt)
            ctr=ctr+1;
%             subplot(length(clus2run),length(trls2plt),ctr); hold on;
%             subplot(1,length(trls2Plt),iPlot);
            figure;
            trials = createTrls(dat,nTrials,locRange,1,jointTrls,boxSize,catsInfo,rSeedAll{iClus,iBvals,iC}(iterI)); hold on;
            
%             scatter(trials(:,1),trials(:,2),5,[.5 .5 .5],'.');
            scatter(trials(trls2plt{iPlot},1),trials(trls2plt{iPlot},2),datPtSiz,colGreyDat,'.');
                        
            colors = repmat(colGreyClus,clus2run(iClus),1);
            if iPlot==2 % if moved and closer to one category than the other, assign intermediate category colour (blue-grey, red-grey)
                moved = sqrt(sum([muAllClus{iClus}(:,1,iPlot,iterI,iBvals,iC)-muAllClus{iClus}(:,1,iPlot-1,iterI,iBvals,iC), muAllClus{iClus}(:,2,iPlot,iterI,iBvals,iC)-muAllClus{iClus}(:,2,iPlot-1,iterI,iBvals,iC)].^2,2));
                dist2catA = sqrt(sum([muAllClus{iClus}(:,1,iPlot,iterI,iBvals,iC)-catCentres(1,1), muAllClus{iClus}(:,2,iPlot,iterI,iBvals,iC)-catCentres(1,2)].^2,2));
                dist2catB = sqrt(sum([muAllClus{iClus}(:,1,iPlot,iterI,iBvals,iC)-catCentres(2,1), muAllClus{iClus}(:,2,iPlot,iterI,iBvals,iC)-catCentres(2,2)].^2,2));
                indA = moved & (dist2catA < dist2catB);
                indB = moved & (dist2catB < dist2catA);
                colors(indA,:) = repmat(catAcol(iPlot-1,:),size(colors(indA,:),1),1);
                colors(indB,:) = repmat(catBcol(iPlot-1,:),size(colors(indB,:),1),1);
                
            elseif iPlot==3 % if within the vicinity of data, assign category color (blue, red)
                dist2catA = sqrt(sum([muAllClus{iClus}(:,1,iPlot,iterI,iBvals,iC)-catCentres(1,1), muAllClus{iClus}(:,2,iPlot,iterI,iBvals,iC)-catCentres(1,2)].^2,2));
                dist2catB = sqrt(sum([muAllClus{iClus}(:,1,iPlot,iterI,iBvals,iC)-catCentres(2,1), muAllClus{iClus}(:,2,iPlot,iterI,iBvals,iC)-catCentres(2,2)].^2,2));
                indA = dist2catA < 7;
                indB = dist2catB < 7;
                colors(indA,:) = repmat(catAcol(iPlot-1,:),size(colors(indA,:),1),1);
                colors(indB,:) = repmat(catBcol(iPlot-1,:),size(colors(indB,:),1),1);
                
            end
            scatter(muAllClus{iClus}(:,1,iPlot,iterI,iBvals,iC),muAllClus{iClus}(:,2,iPlot,iterI,iBvals,iC),1000,colors,'.');hold on;
            
            xlim([locRange(1) locRange(2)+1]);
            ylim([locRange(1) locRange(2)+1]);
            %         xticks([0, 50]); xticklabels({'0', '50'}); yticks(50); yticklabels({'50'});
            xticks([]); xticklabels({''}); yticks([]); yticklabels({''});
            
%             if iClus==1 && iPlot==ceil(length(trls2plt)/2)
% %                 title(sprintf('Category Learning, %d categories; batchSize=%d',nCats, batchSizeVals(iBvals)))
%                 title(sprintf('Category Learning, %d categories',nCats))
%             end
%             if iPlot==1
%                 ylabel(sprintf('%d clus',clus2run(iClus)));
%             end
%             if iClus==3 && iPlot==2
%                 xlabel('Timesteps (start, middle, end)')
%             end
            if iPlot==2
                xlabel('Timesteps (start, middle, end)')
            end
            set(gca,'FontSize',fontSiz,'fontname','Arial')
        
            %save individual figs
            fname = [figDir, sprintf('/catLearn_eps%d_batchSiz%d_%dcats_stoch%d_c%d_nClus%d_iter%d_t%d',epsMuOrig1000,batchSizeVals(iBvals),nCats,stoch,cVals(iC),clus2run(1),iterI,ctr)];
            %     fname = [figDir, sprintf('/catLearn_eps%d_batchSiz%d_%dcats_stoch%d_c%d_nClus%d-%d-%d_iter%d',epsMuOrig1000,batchSizeVals(iBvals),nCats,stoch,cVals(iC),clus2run(1),clus2run(2),clus2run(3),iterI)];
            if savePlots
                set(gcf,'Renderer','painters');
                print(gcf,'-depsc2',fname)
                saveas(gcf,fname,'png');
                close all
            end
        end
    end
    
%     fname = [figDir, sprintf('/catLearn_eps%d_batchSiz%d_%dcats_stoch%d_c%d_nClus%d_iter%d',epsMuOrig1000,batchSizeVals(iBvals),nCats,stoch,cVals(iC),clus2run(1),iterI)];
% %     fname = [figDir, sprintf('/catLearn_eps%d_batchSiz%d_%dcats_stoch%d_c%d_nClus%d-%d-%d_iter%d',epsMuOrig1000,batchSizeVals(iBvals),nCats,stoch,cVals(iC),clus2run(1),clus2run(2),clus2run(3),iterI)];
%     if savePlots
%         set(gcf,'Renderer','painters');
%         print(gcf,'-depsc2',fname)
%         saveas(gcf,fname,'png');
%         close all
%     end
    end
    end
    
end
                