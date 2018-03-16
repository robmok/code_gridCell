%kmeans_merge

clear all;

wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
cd(wd);
codeDir = [wd '/code_gridCell'];
saveDir = [wd '/data_gridCell'];

saveDat = 0;

dat='square'; %square, circ
nPoints = 10000; %3k, 5k, 10k
nIters  = 200; %nKmeans

mergeAcrossClus = 0; % if more than 2 files to merge, can edit below to use this as the value

%% merge

if ~mergeAcrossClus
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % merge across iters, but not clus%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    fname=[saveDir sprintf('/kmeans_nK_3-30_%s_nPoints%d_%diters*',dat,nPoints,nIters)];
    f = dir(fname);
    filesToLoad = cell(1,length(f));
    for iF = 1:length(f)
        filesToLoad{iF} = f(iF).name;
        d(iF)=load(f(iF).name);
    end
    densityPlotCentres = cat(3,d.densityPlotCentres);
    gA                 = cat(1,d.gA);
    gW                 = cat(1,d.gW);
    indSSE1            = cat(2,d.indSSE1);
    indSSE2            = cat(2,d.indSSE2);
    kVals              = d(1).kVals;
    tssekVals          = cat(2,d.tssekVals);
    muAllkVals         = cell(1,length(kVals)); %merge across iters, same nClus conds
    for iK = 1:length(kVals)
        for iBlk=1:length(f) %nFiles
            muAllkVals{iK} = cat(3,muAllkVals{iK}, d(iBlk).muAllkVals{iK});
        end
    end
    
else 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % merge multiple same nclus conds then across multiple nclus conds%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %first set of nClus conds
    fname=[saveDir sprintf('/kmeans_nK_3-22_%s_nPoints%d_%diters*',dat,nPoints,nIters)];
    f = dir(fname);
    filesToLoad = cell(1,length(f));
    for iF = 1:length(f)
        filesToLoad{iF} = f(iF).name;
        d(iF)=load(f(iF).name);
    end
    densityPlotCentres1 = cat(3,d.densityPlotCentres);
    gA1                 = cat(1,d.gA);
    gW1                 = cat(1,d.gW);
    indSSE11            = cat(2,d.indSSE1);
    indSSE21            = cat(2,d.indSSE2);
    kVals1              = d(1).kVals;
    tssekVals1          = cat(2,d.tssekVals);
    muAllkVals1         = cell(1,length(kVals1)); %merge across iters, same nClus conds
    for iK = 1:length(kVals1)
        for iBlk=1:length(f) %nFiles
            muAllkVals1{iK} = cat(3,muAllkVals1{iK}, d(iBlk).muAllkVals{iK});
        end
    end
    
    %second set of nClus conds
    fname=[saveDir sprintf('/kmeans_nK_23-30_%s_nPoints%d_%diters*',dat,nPoints,nIters)];
    f = dir(fname);
    filesToLoad = cell(1,length(f));
    for iF = 1:length(f)
        filesToLoad{iF} = f(iF).name;
        d(iF)=load(f(iF).name);
    end
    densityPlotCentres2 = cat(3,d.densityPlotCentres);
    gA2                 = cat(1,d.gA);
    gW2                 = cat(1,d.gW);
    indSSE12            = cat(2,d.indSSE1);
    indSSE22            = cat(2,d.indSSE2);
    kVals2              = d(1).kVals;
    tssekVals2          = cat(2,d.tssekVals);
    muAllkVals2         = cell(1,length(kVals2)); %merge across iters, same nClus conds
    for iK = 1:length(kVals2)
        for iBlk=1:length(f) %nFiles
            muAllkVals2{iK} = cat(3,muAllkVals2{iK}, d(iBlk).muAllkVals{iK});
        end
    end
    
    %then merge across clus
    densityPlotCentres = cat(4,densityPlotCentres1,densityPlotCentres2);
    gA                 = cat(3,gA1,gA2);
    gW                 = cat(3,gW1,gW2);
    indSSE1            = cat(1,indSSE11,indSSE12);
    indSSE2            = cat(1,indSSE21,indSSE22);
    kVals              = cat(2,kVals1,kVals2);
    tssekVals          = cat(1,tssekVals1,tssekVals2);
    muAllkVals         = cat(2,muAllkVals1,muAllkVals2);
    
end

if saveDat
    fname=[saveDir sprintf('/kmeans_nK_3-30_%s_nPoints%d_%diters_%dsims_merged_2',dat,nPoints,nIters,length(f))]; %note assuming length(f) is the num of sims run (and same for set 1 and set 2)
    save(fname,'densityPlotCentres','gA','gW','indSSE1','indSSE2','kVals','tssekVals','muAllkVals');
end

%%

% d1=load('/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model/data_gridCell/kmeans_nK_3-30_square_nPoints5000_200iters_5sims_merged');
% d2=load('/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model/data_gridCell/kmeans_nK_3-30_square_nPoints5000_200iters_5sims_merged_2');
% % d1=load('/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model/data_gridCell/kmeans_nK_3-30_square_nPoints10000_200iters_5sims_merged');
% % d2=load('/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model/data_gridCell/kmeans_nK_3-30_square_nPoints10000_200iters_5sims_merged_2');
% 
% densityPlotCentres = cat(3,d1.densityPlotCentres,d2.densityPlotCentres);
% gA                 = cat(1,d1.gA,d2.gA);
% gW                 = cat(1,d1.gW,d2.gW);
% indSSE1            = cat(2,d1.indSSE1,d2.indSSE1);
% indSSE2            = cat(2,d1.indSSE2,d2.indSSE2);
% kVals              = d1.kVals;
% tssekVals          = cat(2,d1.tssekVals,d2.tssekVals);
% muAllkVals         = cell(1,length(kVals)); %merge across iters, same nClus conds
% for iK = 1:length(kVals)
%     muAllkVals{iK} = cat(3,d1.muAllkVals{iK}, d2.muAllkVals{iK});
% end
% 
%     fname=[saveDir '/kmeans_nK_3-30_square_nPoints5000_200iters_10sims_merged'];
% % fname=[saveDir '/kmeans_nK_3-30_square_nPoints10000_200iters_10sims_merged'];
% save(fname,'densityPlotCentres','gA','gW','indSSE1','indSSE2','kVals','tssekVals','muAllkVals');