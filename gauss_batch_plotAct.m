% plot gauss activations (batch)

densityPlotAct        = zeros(50,50);
densityPlotAct1       = zeros(50,50); %one cluster
figure; hold on;
xlim([0,50]); ylim([0,50]);

zlims=[15 40];
for iTrl = nTrials*.5:nTrials
    
    densityPlotAct(trials(iTrl,1)+1, trials(iTrl,2)+1) = densityPlotAct(trials(iTrl,1)+1, trials(iTrl,2)+1)+ sum(actAll(:,iTrl)).^2; %squaring might make it look better..
    if mod(iTrl,2000)==0 %plot centers after x trials
%     imagesc(densityPlotAct,zlims); drawnow;
    imagesc(imgaussfilt(densityPlotAct,1)); drawnow;
    end
    
%     densityPlotAct1(trials(iTrl,1)+1, trials(iTrl,2)+1) = densityPlotAct1(trials(iTrl,1)+1, trials(iTrl,2)+1)+ (actAll(1,iTrl)).^2;
%     if mod(iTrl,2000)==0 %plot centers after x trials
%         imagesc(imgaussfilt(densityPlotAct1,1)); drawnow;
%     end
    
    
end

%%
zlims=[10 20];

figure;
imagesc(densityPlotAct,zlims);
figure;
imagesc(imgaussfilt(densityPlotAct,1),zlims)
% 
% zlims1=[.5 1.5];
% figure;
% imagesc(densityPlotAct1,zlims1);
% figure;
% imagesc(imgaussfilt(densityPlotAct1,1),zlims1)
%     
%%

aCorrMap=ndautoCORR(imgaussfilt(densityPlotAct,1)); %autocorrelogram
[g,gdataA] = gridSCORE(aCorrMap,'allen',1);