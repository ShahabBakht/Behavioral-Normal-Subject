function SuppIdx = simulateSuppIdx(slope,integR,numSamples)

sizes=sqrt(0)+[3,10];%(.2:2:10.2); % the avg eccentricity of the stimuli is sqrt(5) deg
popSize=zeros(1,2);
errSize=zeros(1,2); % reviewer #2 asked us to calculate the error of our estimates
% fun=@(x1) (1.14*x1.^-0.76); % Albright & Desimone 87 Exp Brain Res (mm/deg)
fun=@(x1) (6*x1.^-0.9); % Erickson (mm/deg)
for i=1:length(sizes) % calculate the # of  neurons for each size
%     popSize(i)=integral(fun,0.5,sizes(i))+1; % integral from 0.5 to edge, between 0 to 0.5 is roughly 1 mm (2mm/deg)
    popSize(i)=integral(fun,12 - sizes(i),12 + sizes(i)); % integral from 1 to edge, between 0 to 1 is roughly (9mm/deg) erickson et al. EBR
%     popSize(i)=integral(fun,.5,sizes(i)) + 9;

    errSize(i)=integral(fun,0.5,sizes(i)*0.16+0.51)+1; % integral of the error
end
% popSize=round(1*popSize.^2); % multiple by a constant of 10, but the constant did not affect the trend
popSize=round(20*popSize); % new constant


%%

tauD = .25;
tauS = .5;
rmax = .25;

% set decoding parameters -  read-out
param(1) = integR;
param(2) = slope;
param(3) = .5;

for trial = 1:numSamples
    fprintf('.')
    Target.MotionSpeed = 15 * ones(1000,1);
    Target.MotionDirection = 0 * ones(1000,1);
    
    k = 0;
    for i = popSize%max(popSize)./popSize%[1,1.5,2,3,4]%log2(2:(62/5):64-62/5)
        %  i = 1;
        k = k + 1;
        if sizes(k) > 4000
            thisSize = 4;
        else
            thisSize = sizes(k);
        end
        Target.Size = thisSize * ones(1000,1);
        
        nps(k) = round(sqrt(i));%floor(46/(1*i));
        npd(k) = round(sqrt(i));%floor(90/(1*i));
        
        
        [mtpopulation, R, COV] = MTpopulation(Target,nps(k),npd(k),rmax,tauD,tauS);
        [mtpopulation_denom, R_denom, COV_denom] = MTpopulation(Target,nps(k),npd(k),rmax,tauD,tauS);
        [TargetEstimate] = DecodeMTpopulation(mtpopulation, R, mtpopulation_denom, R_denom,'weighted_uncorr_norm',param);
        
        SPD(k,trial,:) = TargetEstimate.SPDest;
        DIR(k,trial,:) = TargetEstimate.DIRest;
        COVest = cov(real(R)');
        FiringRate_sum(k,trial,:) = sum(real(R),1);
        FiringRate_mean(k,trial,:) = nanmean(real(R),1);
        SI = cellfun(@(x)(x.SuppressionIndex),mtpopulation);
        
        FR_sum_ss(k,trial,:) = sum(real(R(SI>quantile(SI,.4),:)),1);
        FR_sum_nss(k,trial,:) = sum(real(R(SI<=quantile(SI,.4),:)),1);
        FR_mean_ss(k,trial,:) = nanmean(real(R(SI>quantile(SI,.4),:)),1);
        FR_mean_nss(k,trial,:) = nanmean(real(R(SI<=quantile(SI,.4),:)),1);
        
        
        COVss(k,trial) = mean2(COVest(SI>quantile(SI,.4),SI>quantile(SI,.4)));
        COVnss(k,trial) = mean2(COVest(SI<=quantile(SI,.4),SI<=quantile(SI,.4)));
        
        SNRnss(k,trial) = nanmean(sum(real(R(SI<=quantile(SI,.4),:)),1))./sqrt(sum(sum(triu(COVest(SI<=quantile(SI,.4),SI<=quantile(SI,.4))))));
        SNRss(k,trial) = nanmean(sum(real(R(SI>quantile(SI,.4),:)),1))./sqrt(sum(sum(triu(COVest(SI>quantile(SI,.4),SI>quantile(SI,.4))))));
        SNR(k,trial) = nanmean(sum(real(R),1))./sqrt(sum(sum(triu(COVest))));
        
        
        
        Rmean{k} = nanmean(R,2);
        
        clear mtpopulation R COV TargetEstimate
        
    end
    
    
    
end

SPD(abs(SPD)>100) = nan;
SPDvar = squeeze(nanstd(SPD(:,:,:,:),[],3));
SuppIdx = SPDvar(1,:)./SPDvar(2,:);


fprintf('| \n')
end