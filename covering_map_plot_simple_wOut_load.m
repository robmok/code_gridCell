%% plot centres 

nClus = size(muAvg,1);
colors = distinguishable_colors(nClus); %function for making distinguishable colors for plotting

iterI=1; 

%%%%%
% currently looking at 1 iteration, comparing how i'm averaging over
% trials;. may want to compare iters as well, or plot differently to show
% differences better
%%%%


for iSet=1:size(muAvg,3) %plot - diff averaging over nTrials
    figure; hold on;
    scatter(muAvg(:,1,iSet,iterI),muAvg(:,2,iSet,iterI),20e+2,colors,'.');
    xlim(locRange); ylim(locRange);
    voronoi(muAvg(:,1,iSet,iterI),muAvg(:,2,iSet,iterI),'k');
end


%% density map, autocorrelogram

gaussSmooth=1;

for iSet=1:size(muAvg,3) %plot - diff averaging over nTrials

    densityPlot = sum(densityPlotClus(:,:,:,iSet,iterI),3);
    densityPlot = imgaussfilt(densityPlot,gaussSmooth);
    aCorrMap=ndautoCORR(densityPlot); %autocorrelogram

    figure; hold on;
    subplot(1,2,1); imagesc(densityPlot);
    subplot(1,2,2); imagesc(aCorrMap);

end
