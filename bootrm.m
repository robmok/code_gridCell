%function for getting bootstrap intervals on any value
% dat is a vector
function [bootstrapVals] = bootrm(dat,nBoot)

if vargin==1
    nBoot = numel(dat); 
end

%bootstrapping

    
%add an option for what function you want to do, means, percent, std;
%maybe use stuf like @mean for the function
    
    
means=nan(1,nBoot);
for iBoot = 1:nBoot
    s = randsample(dat,length(dat),true); %sample with replacement
    means(iBoot) = nanmean(s);
end

bootstrapVals = prctile(means,[2.5, 97.5]);
end