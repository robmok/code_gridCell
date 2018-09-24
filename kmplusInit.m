function muInit = kmplusInit(dataPts,nK)
%kmeans++ initiatialization

for i=1:nK    
    if i==1% random datapoint as 1st cluster
        muInit(i,:,1) = dataPts(randi(length(dataPts)),:);
    end
    if i~=nK % no need update k+1
        clear distInit
        for iClus = 1:i% loop over clusters that exist now
            distInit(:,iClus)=sum([muInit(iClus,1,1)-dataPts(:,1),  muInit(iClus,2,1)-dataPts(:,2)].^2,2); %squared euclid for k means
        end
        [indValsInit, indInit]=min(distInit,[],2); % find which clusters are points closest to
        
        distClus=[];
        for iClus = 1:i
            indOrig(:,i) = indInit==iClus;
            distClusTmp = sum([(muInit(iClus,1,1)-dataPts(indOrig(:,i),1)), (muInit(iClus,2,1)-dataPts(indOrig(:,i),2))].^2,2);
            distClus = [distClus; [distClusTmp, repmat(iClus,length(distClusTmp),1)]];
        end
        
        % keep track of the indices of the original dist variable - get datapoints
        % that were the furthest from all clusters, get that cluster and see which
        % datapoint that was relative to that cluster (since i just save the distance)
        distClusNorm = distClus(:,1)./sum(distClus(:,1));
        distClusPr   = cumsum(distClusNorm(:,1)); %get cumsum, then generate rand val from 0 to 1 and choose smallest val - the larger the dis, the more likely the rand value will lie between it and its previous value in a cumsum plot
        ind=find(rand(1)<distClusPr,1);% %find smallest value that is larger than the random value (0 to 1 uniform distr)
        %tmp(i)=distClus(ind,1); %testing if getting the right values; if plot, see that it should be lower pr for closer values, higher for larger. note if very few high distances, this hist will look normally distributed
        
        clusInd = distClus(ind,2); %find which is the closest cluster
        indDat = find(distInit(:,clusInd)==distClus(ind,1)); %find where the datapoint is in the original vector
        
        if size(indDat,1)>1
            indDat=indDat(randi(size(indDat,1),1));
        end
        muInit(i+1,:,1) = dataPts(indDat,:);
    end
end