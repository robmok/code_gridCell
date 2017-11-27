wd='/Users/robertmok/Documents/Postdoc_ucl/Grid_cell_model';
% wd='/Users/robert.mok/Documents/Postdoc_ucl/Grid_cell_model';

cd(wd);
addpath(genpath([wd '/gridSCORE_packed']));

load('ex_amap.mat')
im=amap;
[g,gdata] = gridSCORE(im,'allen'); %allen or wills



%% testing it out

nTrials=size(muAll,3);
nIter = size(muAll,4);

nTrlsToUse = 10000; 
spacing = linspace(locRange(1),locRange(2),locRange(2)+1);
densityPlot = zeros(length(spacing),length(spacing),nIter);
for iterI=1:nIter
    clus = round(muAll(:,:,nTrials-nTrlsToUse+1:nTrials,iterI));
    for iTrl=1:nTrlsToUse
        for i=1:nClus
            densityPlot(clus(i,1,iTrl),clus(i,2,iTrl),iterI)=densityPlot(clus(i,1,iTrl),clus(i,2,iTrl),iterI)+1; % works, but better way / faster to vectorise?
        end
    end
    densityPlot(:,:,iterI) = imgaussfilt(densityPlot(:,:,iterI),1); %smooth
%     aCorrMap(:,:,iterI) = ndautoCORR(densityPlot(:,:,iterI));

    figure;
    imagesc(densityPlot(:,:,iterI));
%     imagesc(densityPlot(:,:,iterI),[1000 4000]);
end


%%

figure;
iPlot=0;
for iter=1:nIter
    iPlot=iPlot+1;
%     subplot(2,2,iPlot);
figure;
    im=aCorrMap(:,:,iter);
    [g,gdata] = gridSCORE(im,'allen'); %allen or wills
%     [g,gdata] = gridSCORE(im,'wills'); %allen or wills
end