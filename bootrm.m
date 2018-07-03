function bootstrapVals = bootrm(dat,nBoot,type,threshVal)
%function for getting bootstrap intervals on any value
% dat is a vector nSamples x 1
% type is a string - mean, percentGr, percentLs
% threshVal is a single number to find percent greater/less than threshVal

if nargin==1
    nBoot = numel(dat); 
end
if nargin<3 % if compute mean only, nanmean 2nd arg is dimension
    threshVal = 1; 
end   

%what to bootstrap: means, percentGr, percentLs, std
switch type
    case 'mean'
        func=@nanmean;
    case 'percentGr'
        func=@(n,thresh)nanmean(n>thresh);
    case 'percentLs'
        func=@(n,thresh)nanmean(n<thresh);
end

%do bootstrap
bootVal=nan(1,nBoot);
for iBoot = 1:nBoot
    s = randsample(dat,length(dat),true); %sample with replacement
    bootVal(iBoot) = func(s,threshVal);    
end
bootstrapVals = prctile(bootVal,[2.5, 97.5]);
end