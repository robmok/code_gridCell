%% Set up model

% Sakamoto et al., 2008 Gaussian model
% 10 blocks of training trials; each block had each training example once
% (12)in random order

clear all

trainSetA = 140:2:150;
trainSetB = 230:20:330;
trainSet  = [trainSetA, trainSetB];
nCondsTrain = numel(trainSet);

allConds=[120 130 trainSetA 160 170 180 190 200 210 220 trainSetB 340 350];

%functions
% fitGauss=@(P,X)(normpdf(X,P(1),P(2)));
fitGauss=@(x,mu,sigma)exp(-0.5 .* ((x - mu)./sigma).^2)/(sqrt(2.*pi) .* sigma); % (-0.5 * ((x - mu)./sigma).^2) is equavalent to simply writing it out: -((x-mu).^2)./(2.*sigma.^2)

% parameters to fit changes in mean and sigma
epsMu=98000;
epsSigma=70000;
alpha=.05; % how much to shift (weighted by diff from cluster mean) if feedback says 'correct category'

%fit change in mean
% fitDeltaMu=@(epsMu,fdbck,clusAct,x,mu,sigma)epsMu.*(fdbck-clusAct).*(((x-mu)./((sqrt(2.*pi)).*sigma.^3)).*exp((-((x-mu).^2))./(2.*sigma.^2)));
%shorter (ish), more similar to above % exp((-((x-mu).^2))./(2.*sigma.^2)) is equavalent to: exp(-0.5*((x-mu)./sigma).^2)
fitDeltaMu=@(x,mu,sigma,fdbck,act)exp(-0.5*((x-mu)/sigma)^2)*(x-mu)/((sqrt(2*pi))*sigma^3)*(epsMu*(fdbck-act));

% fit change in variance:
% fitDeltaSigma=@(epsSigma,fdbck,clusAct,x,mu,sigma)(epsSigma.*(fdbck-clusAct)).*(((x-mu).^2)-sigma.^2)./((sqrt(2.*pi)).*sigma.^4).*(exp(-((x-mu).^2)/(2.*sigma.^2)))
fitDeltaSigma=@(x,mu,sigma,fdbck,act)exp(-0.5*((x-mu)/sigma)^2)*(((x-mu)^2)-sigma^2)/((sqrt(2*pi))*sigma^4)*(epsSigma*(fdbck-act));

% epsilon is learning rate, fdbck is equal to clusAct if stim x is in category, and 0 if not (so if correct then no adjustment)
% episilon is for each category m. fdbck, clusAct, mu, sd are per trial

%% run

for iter=1:1000
    fprintf('Running iteration %d \n',iter);
    trials = [];
    for iBlock=1:10,
        trials=[trials, trainSet(randperm(length(trainSet)))];
    end
    
    %initialise cluster at first instance
    firstA = nan(1,7); firstB=nan(1,7);
    for iTrl=1:7, %will be at least one of each cat at 7th
        firstA(iTrl)=any(trials(iTrl)==trainSetA);
        firstB(iTrl)=any(trials(iTrl)==trainSetB);
    end
    firstAind=find(firstA);
    firstBind=find(firstB);
    
    % define size of variables
    actA = nan(1,length(trials)); actB = nan(1,length(trials));
    fdbckA = nan(1,length(trials)); fdbckB = nan(1,length(trials));
    deltaMuA = nan(1,length(trials)); deltaMuB = nan(1,length(trials));
    deltaSigmaA = nan(1,length(trials)); deltaSigmaB = nan(1,length(trials));
    muA = nan(1,length(trials)); muB = nan(1,length(trials));
    sigmaA = nan(1,length(trials)); sigmaB = nan(1,length(trials));
    
    %parameter starting points
    muA(1)    =  trials(firstAind(1));
    sigmaA(1) =  20;
    muB(1)    =  trials(firstBind(1));
    sigmaB(1) =  20;
    
    for iTrl=1:length(trials)
        
        %compute activation
        actA(iTrl)=fitGauss(trials(iTrl),muA(iTrl),sigmaA(iTrl));
        actB(iTrl)=fitGauss(trials(iTrl),muB(iTrl),sigmaB(iTrl));
        %compute feedback
        if any(trials(iTrl)==trainSetA)
            fdbckA(iTrl)=alpha;
            fdbckB(iTrl)=0;
        elseif any(trials(iTrl)==trainSetB)
            fdbckA(iTrl)=0;
            fdbckB(iTrl)=alpha;
        end
        
        %compute change in Mu and Sigma
        deltaMuA(iTrl)=fitDeltaMu(trials(iTrl),muA(iTrl),sigmaA(iTrl),fdbckA(iTrl),actA(iTrl));
        deltaMuB(iTrl)=fitDeltaMu(trials(iTrl),muB(iTrl),sigmaB(iTrl),fdbckB(iTrl),actB(iTrl));
        deltaSigmaA(iTrl)=fitDeltaSigma(trials(iTrl),muA(iTrl),sigmaA(iTrl),fdbckA(iTrl),actA(iTrl));
        deltaSigmaB(iTrl)=fitDeltaSigma(trials(iTrl),muB(iTrl),sigmaB(iTrl),fdbckB(iTrl),actB(iTrl));
        
        % update mean and sigma estimates
        if iTrl~=length(trials) %no need to update for last trial +1)
            muA(iTrl+1) = muA(iTrl) + deltaMuA(iTrl);
            muB(iTrl+1) = muB(iTrl) + deltaMuB(iTrl);
            sigmaA(iTrl+1) = sigmaA(iTrl) + deltaSigmaA(iTrl);
            sigmaB(iTrl+1) = sigmaB(iTrl) + deltaSigmaB(iTrl);
        end
        
    end
    
    muA_blk(iter)=muA(end);
    muB_blk(iter)=muB(end);
    sigmaA_blk(iter)=sigmaA(end);
    sigmaB_blk(iter)=sigmaB(end);
    
    
%     % plot trajectory of mu and sigma over one block
%     if iter == 1 | iter == 100 | iter == 500,
%         figure;
%         subplot(1,2,1); hold on;
%         plot(muA);
%         subplot(1,2,2);
%         plot(muB);
%         
%         figure;
%         subplot(1,2,1); hold on;
%         plot(sigmaA);
%         subplot(1,2,2);
%         plot(sigmaB);
%     end
end



% figure;
% subplot(1,2,1); hold on;
% plot(deltaMuA);
% subplot(1,2,2);
% plot(deltaMuB);
% 
% figure;
% subplot(1,2,1); hold on;
% plot(deltaSigmaA);
% subplot(1,2,2);
% plot(deltaSigmaB);

%% Test & plot

trials = allConds;
% trials = 120:1:350; % plot nicer curves

%compute activation on new stimuli
actA_test=fitGauss(trials,mean(muA_blk),mean(sigmaA_blk));
actB_test=fitGauss(trials,mean(muB_blk),mean(sigmaB_blk));


figure;
plot(trials,[actA_test;actB_test]')

% border item greater act for cat B
actA_test(12)<actB_test(12)
