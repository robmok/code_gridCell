function [muEnd,tsse] = kmeans_rm(mu,dataPts,nK,nUpdSteps)
%kmeans function 
%first trial of mu has to be specified (e.g. by kmplusInit.m)

for upd = 1:nUpdSteps
%     dist = nan(length(dataPts),nK);
    %vectorized
    dist = sum(reshape([mu(:,1,upd)'-dataPts(:,1),  mu(:,2,upd)'-dataPts(:,2)].^2,length(dataPts),nK,2),3); %squared euclid for k means

    [indVals, ind]=min(dist,[],2); % find which clusters are points closest to
    for iK = 1:nK
        if any(isnan([mean(dataPts(ind==iK,1)), mean(dataPts(ind==iK,2))]))
            mu(iK,:,upd+1)=mu(iK,:,upd);
        else
            mu(iK,:,upd+1)=[mean(dataPts(ind==iK,1)), mean(dataPts(ind==iK,2))];
        end
    end
end

muEnd=mu(:,:,end);

%compute sum of squared errors across clusters
sse  = nan(1,nK);
for iK = 1:nK
    sse(iK)=nansum(nansum([(muEnd(iK,1))-dataPts(ind==iK,1), (muEnd(iK,2))-dataPts(ind==iK,2)].^2,2));
end
tsse = nansum(sse);


end
