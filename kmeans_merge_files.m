%kmeans_merge
wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';
cd(wd);
codeDir = [wd '/code_gridCell'];
saveDir = [wd '/data_gridCell'];

saveDat = 1;

dat='circ'; %square, circ
nPoints=3000; %3k, 5k, 10k

fname=[saveDir sprintf('/kmeans_nK_3-17_%s_nPoints%d_1000iters',dat,nPoints)];
d1=load(fname);
fname=[saveDir sprintf('/kmeans_nK_18-25_%s_nPoints%d_1000iters',dat,nPoints)];
d2=load(fname);
fname=[saveDir sprintf('/kmeans_nK_26-30_%s_nPoints%d_1000iters',dat,nPoints)];
d3=load(fname);

% % new - edit
% fname=[saveDir sprintf('/kmeans_nK_3-17_%s_nPoints%d_1000iters_*',dat,nPoints)];
% f=dir(fname);
% d1=load([saveDir '/' f(1).name]); %need edit since can load in more than 1
% fname=[saveDir sprintf('/kmeans_nK_18-25_%s_nPoints%d_1000iters_*',dat,nPoints)];
% f=dir(fname);
% d2=load([saveDir '/' f(1).name]);
% fname=[saveDir sprintf('/kmeans_nK_26-30_%s_nPoints%d_1000iters_*',dat,nPoints)];
% f=dir(fname);
% d3=load([saveDir '/' f(1).name]);
densityPlotCentres  = cat(4,d1.densityPlotCentres,d2.densityPlotCentres,d3.densityPlotCentres);
gA                  = cat(3,d1.gA,d2.gA,d3.gA);
gW                  = cat(3,d1.gW,d2.gW,d3.gW);
indSSE1             = cat(1,d1.indSSE1,d2.indSSE1,d3.indSSE1);
indSSE2             = cat(1,d1.indSSE2,d2.indSSE2,d3.indSSE2);
kVals               = cat(2,d1.kVals,d2.kVals,d3.kVals);
muAllkVals          = cat(2,d1.muAllkVals,d2.muAllkVals,d3.muAllkVals);
tssekVals           = cat(1,d1.tssekVals,d2.tssekVals,d3.tssekVals);
if saveDat
    fname=[saveDir sprintf('/kmeans_nK_3-30_%s_nPoints%d_1000iters',dat,nPoints)];
%     fname=[saveDir sprintf('/kmeans_nK_3-30_%s_nPoints%d_1000iters_2',dat,nPoints)];
    save(fname,'muAllkVals','tssekVals', 'gA','gW','densityPlotCentres','indSSE1','indSSE2','kVals')
end



%% new - merge multiple same nclus conds then across multiple nclus conds


%first set of nClus conds

fname=[saveDir sprintf('/kmeans_nK_3-22_%s_nPoints%d_1000iters*',dat,nPoints)];
f = dir(fname); 
filesToLoad = cell(1,length(f));
% d1=load(fname);
%             muAllTmp={}; nIterCount=0;
for iF = 1%:length(f)
    filesToLoad{iF} = f(iF).name;
    dTmp=load(f(iF).name);
    d(iF)=load(f(iF).name); %would this work? or below
%     d(iF) = dTmp;
    



%     need?
%     nIterCount = [nIterCount, nIter+nIterCount(end)]; %count number of iters to index below when merging
end

% merge here e.g. 
% cat(2,d(:).gA); %consider which dim to cat over

% note muAllkVals will be struct array








%second set of nClus conds
fname=[saveDir sprintf('/kmeans_nK_23-30_%s_nPoints%d_1000iters*',dat,nPoints)];
f = dir(fname); 



%then merge across clus (can use original code in above cell)



