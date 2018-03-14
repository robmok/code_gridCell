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
d3=load(fname

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