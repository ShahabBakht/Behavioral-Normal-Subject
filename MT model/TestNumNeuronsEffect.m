%% Number of neurons as a function of stimulus size

sizes=sqrt(0)+[1,3,10];%(.2:2:10.2); % the avg eccentricity of the stimuli is sqrt(5) deg
popSize=zeros(1,3);
errSize=zeros(1,3); % reviewer #2 asked us to calculate the error of our estimates
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
errSize=round(10*errSize); % multiple the error a constant of 10



% figure % plot for the manuscript
% plot(sizes-sqrt(0),popSize,'k')
% set(gca,'XTick',[0 5 10 15]);
% set(gca,'YTick',[400 600 800 1000]);
% ylim([400 1100])
% set(gca,'TickDir','out')
% xlabel(['Size (' char(176) ')']); ylabel('Estimated pool size')
% 
% box off

%% Single Speed for all the trials
% % tauD = 1.5*.125;%0.1125;%0.1125;%
% % tauS = 1.5*.15;%0.135;%0.1350;%
% tauD = .25;
% tauS = .5;
% % Target.MotionSpeed        =   [15 * ones(500,1);15 * ones(500,1);...
% %     15 * ones(500,1); 15 * ones(500,1); ...
% %     15 * ones(500,1)];
% % Target.MotionDirection    =   [40 * oxnes(500,1);40 * ones(500,1);...
% %     40 * ones(500,1); 40 * ones(500,1); ...
% %     40 * ones(500,1)];
% % 
% % Target.Size = [sizes(end-4) * ones(500,1);sizes(end-3) * ones(500,1);...
% %     sizes(end-2) * ones(500,1); sizes(end-1) * ones(500,1); ...
% %     sizes(end) * ones(500,1)];
% 
% for trial = 1:1000
% Target.MotionSpeed = 15 * ones(1000,1);
% Target.MotionDirection = 0 * ones(1000,1);
% 
% j = 0;
% % DIRstd = nan(length([1,1.5,2,3,4]),length([0,0.09,0.18,0.27,0.36,0.45]));
% % DIRm = nan(length([1,1.5,2,3,4]),length([0,0.09,0.18,0.27,0.36,0.45]));
% % SPDstd = nan(length([1,1.5,2,3,4]),length([0,0.09,0.18,0.27,0.36,0.45]));
% % SPDm = nan(length([1,1.5,2,3,4]),length([0,0.09,0.18,0.27,0.36,0.45]));
% for rmax = [.25]
%     
%     
%     j = j + 1;
%     k = 0;
% for i = popSize%max(popSize)./popSize%[1,1.5,2,3,4]%log2(2:(62/5):64-62/5)
% %  i = 1;    
%     k = k + 1;
%     if sizes(k) > 4000
%        thisSize = 4;
%     else
%         thisSize = sizes(k);
%     end
%     Target.Size = thisSize * ones(1000,1);
%     
%     nps(k) = round(sqrt(i));%floor(46/(1*i));
%     npd(k) = round(sqrt(i));%floor(90/(1*i));
%     
%     % set decoding parameters -  read-out
%     param(1) = 3;
%     param(2) = 2;
%     param(3) = .5;
% %     
% %     nps(k) = min(round(sqrt(i)),12);
% %     npd(k) = min(round(sqrt(i)),12);
% %     
%     [mtpopulation, R, COV] = MTpopulation(Target,nps(k),npd(k),rmax,tauD,tauS);
%     [mtpopulation_denom, R_denom, COV_denom] = MTpopulation(Target,nps(k),npd(k),rmax,tauD,tauS);
%     [TargetEstimate] = DecodeMTpopulation(mtpopulation, R, mtpopulation_denom, R_denom,'weighted_uncorr_norm',param);
% %     [TargetEstimate] = DecodeMTpopulation(mtpopulation, real(R));
% %     DIRstd(k,j,trial) = nanstd(TargetEstimate.DIRest);
% %     SPDstd(k,j,trial) = nanstd(TargetEstimate.SPDest);
% %     DIRm(k,j,trial) = nanmedian(TargetEstimate.DIRest);
% %     SPDm(k,j,trial) = nanmedian(TargetEstimate.SPDest);
%     SPD(k,j,trial,:) = TargetEstimate.SPDest;
%     DIR(k,j,trial,:) = TargetEstimate.DIRest;
%     COVest = cov(real(R)');
%     FiringRate_sum(k,j,trial,:) = sum(real(R),1);
%     FiringRate_mean(k,j,trial,:) = nanmean(real(R),1);
%     SI = cellfun(@(x)(x.SuppressionIndex),mtpopulation);
%    
%     FR_sum_ss(k,j,trial,:) = sum(real(R(SI>quantile(SI,.4),:)),1);
%     FR_sum_nss(k,j,trial,:) = sum(real(R(SI<=quantile(SI,.4),:)),1);
%     FR_mean_ss(k,j,trial,:) = nanmean(real(R(SI>quantile(SI,.4),:)),1);
%     FR_mean_nss(k,j,trial,:) = nanmean(real(R(SI<=quantile(SI,.4),:)),1);
%     
%     
%     COVss(k,j,trial) = mean2(COVest(SI>quantile(SI,.4),SI>quantile(SI,.4)));
%     COVnss(k,j,trial) = mean2(COVest(SI<=quantile(SI,.4),SI<=quantile(SI,.4)));
%     
%     SNRnss(k,j,trial) = nanmean(sum(real(R(SI<=quantile(SI,.4),:)),1))./sqrt(sum(sum(triu(COVest(SI<=quantile(SI,.4),SI<=quantile(SI,.4))))));
%     SNRss(k,j,trial) = nanmean(sum(real(R(SI>quantile(SI,.4),:)),1))./sqrt(sum(sum(triu(COVest(SI>quantile(SI,.4),SI>quantile(SI,.4))))));
%     SNR(k,j,trial) = nanmean(sum(real(R),1))./sqrt(sum(sum(triu(COVest))));
%   
%     
%     
%     Rmean{k} = nanmean(R,2);
%     
%     clear mtpopulation R COV TargetEstimate
%     
% end
% 
% end
% end
%% Plots - Speed and Direction variance as a function of #neurons
% 
% 
% SPD(abs(SPD)>100) = nan;
% SPDvar = squeeze(nanstd(SPD(:,:,:,:),[],4));
% SPDvar_std = nanstd(SPDvar,[],2);
% SPDvar_mean = nanmedian(SPDvar,2);
% figure(1);hold on;hh = ploterr(sizes, SPDvar_mean, [], SPDvar_std./sqrt(1000),'r','abshhxy', .3);
% hold on;plot(sizes, SPDvar_mean,'.r','MarkerSize',30);
% xlabel('size (degree)');ylabel('standard deviation (degree)');%legend(num2str(0),num2str(0.09),num2str(0.18),num2str(0.27),num2str(0.36),num2str(0.45))

% SPD(abs(SPD)>100) = nan;
% SPDmean = squeeze(nanmean(SPD(:,:,:,:),4));
% SPDmean_std = nanstd(SPDmean,[],2);
% SPDmean_mean = nanmedian(SPDmean,2);
% figure(2);hold on;hh = ploterr(sizes, SPDmean_mean, [], SPDmean_std./sqrt(1000),'g','abshhxy', .3);
% hold on;plot(sizes, SPDmean_mean,'.g','MarkerSize',30);
% xlabel('size (degree)');ylabel('speed (degree/s)');%legend(num2str(0),num2str(0.09),num2str(0.18),num2str(0.27),num2str(0.36),num2str(0.45))


% % 
% % SNR(abs(SNR)>100) = nan;
% SNRmean = squeeze(nanmean(SNR,3));
% % SNRmean = nanmean(SNRmean1,2);
% SNRstd = squeeze(nanstd(SNR,[],3));
% 
% figure;hh = ploterr(sizes, SNRmean, [], SNRstd./sqrt(1000),'k','abshhxy', .3);
% hold on;plot(sizes, SNRmean,'.k','MarkerSize',30);
% xlabel('size (degree)');ylabel('S/N (degree)');%legend(num2str(0),num2str(0.09),num2str(0.18),num2str(0.27),num2str(0.36),num2str(0.45))

% FiringRate(abs(FiringRate)>100) = nan;
% FiringRatemean = squeeze(nanmean(FiringRate,3));
% % SNRmean = nanmean(SNRmean1,2);
% FiringRatestd = squeeze(nanstd(FiringRate,[],3));

% figure;hh = ploterr(sizes, FiringRatemean, [], FiringRatestd./sqrt(1000),'k','abshhxy', .3);
% hold on;plot(sizes, FiringRatemean,'.k','MarkerSize',30);
% xlabel('size (degree)');ylabel('FR');

% SPDstd_mean = median(SPDstd,3);
% SPDstd_mean(SPDstd_mean>100) = nan;
% figure;plot(sizes,SPDstd_mean,'-o','LineWidth',2);
% xlabel('size');ylabel('speed variance')
% % legend(num2str(0),num2str(0.09),num2str(0.18),num2str(0.27))
% figure;plot(sizes,median(DIRstd,3),'-o','LineWidth',2);
% xlabel('size');ylabel('direction variance')
% 
% figure;plot(nps.*npd,median(SPDm,3),'-o','LineWidth',2);
% xlabel('number of neurons');ylabel('speed mean')
% figure;plot(nps.*npd,median(DIRm,3),'-o','LineWidth',2);
% xlabel('number of neurons');ylabel('direction mean')


%% Calculate spem suppression
% % tauD = 1.5*.125;%0.1125;%0.1125;%
% % tauS = 1.5*.15;%0.135;%0.1350;%
% tauD = 1.0*.125;
% tauS = 1.0*.15;
% % Target.MotionSpeed        =   [15 * ones(500,1);15 * ones(500,1);...
% %     15 * ones(500,1); 15 * ones(500,1); ...
% %     15 * ones(500,1)];
% % Target.MotionDirection    =   [40 * ones(500,1);40 * ones(500,1);...
% %     40 * ones(500,1); 40 * ones(500,1); ...
% %     40 * ones(500,1)];
% % 
% % Target.Size = [sizes(end-4) * ones(500,1);sizes(end-3) * ones(500,1);...
% %     sizes(end-2) * ones(500,1); sizes(end-1) * ones(500,1); ...
% %     sizes(end) * ones(500,1)];
% 
% for iter = 1:10
% for trial = 1:100
% Target.MotionSpeed = 15 * ones(1000,1);
% Target.MotionDirection = 0 * ones(1000,1);
% 
% j = 0;
% % DIRstd = nan(length([1,1.5,2,3,4]),length([0,0.09,0.18,0.27,0.36,0.45]));
% % DIRm = nan(length([1,1.5,2,3,4]),length([0,0.09,0.18,0.27,0.36,0.45]));
% % SPDstd = nan(length([1,1.5,2,3,4]),length([0,0.09,0.18,0.27,0.36,0.45]));
% % SPDm = nan(length([1,1.5,2,3,4]),length([0,0.09,0.18,0.27,0.36,0.45]));
% for rmax = [85]
%     
%     
%     j = j + 1;
%     k = 0;
% for i = popSize%max(popSize)./popSize%[1,1.5,2,3,4]%log2(2:(62/5):64-62/5)
% %  i = 1;    
%     k = k + 1;
%     if sizes(k) > 4000
%        thisSize = 4;
%     else
%         thisSize = sizes(k);
%     end
%     Target.Size = thisSize * ones(1000,1);
%     nps(k) = round(sqrt(i));%floor(46/(1*i));
%     npd(k) = round(sqrt(i));%floor(90/(1*i));
% %     
% %     nps(k) = min(round(sqrt(i)),9);
% %     npd(k) = min(round(sqrt(i)),9);
%     
%     [mtpopulation, R, COV] = MTpopulation(Target,nps(k),npd(k),rmax,tauD,tauS);
%     [mtpopulation_denom, R_denom, COV_denom] = MTpopulation(Target,nps(k),npd(k),rmax,tauD,tauS);
%     [TargetEstimate] = DecodeMTpopulation(mtpopulation, R, mtpopulation_denom, R_denom,'weighted_uncorr_norm');
% %     [TargetEstimate] = DecodeMTpopulation(mtpopulation, real(R));
% %     DIRstd(k,j,trial) = nanstd(TargetEstimate.DIRest);
% %     SPDstd(k,j,trial) = nanstd(TargetEstimate.SPDest);
% %     DIRm(k,j,trial) = nanmedian(TargetEstimate.DIRest);
% %     SPDm(k,j,trial) = nanmedian(TargetEstimate.SPDest);
%     SPD(k,j,trial,:) = TargetEstimate.SPDest;
%     DIR(k,j,trial,:) = TargetEstimate.DIRest;
%     COVest = cov(real(R)');
%     FiringRate_sum(k,j,trial,:) = sum(real(R),1);
%     FiringRate_mean(k,j,trial,:) = nanmean(real(R),1);
%     SI = cellfun(@(x)(x.SuppressionIndex),mtpopulation);
%    
%     FR_sum_ss(k,j,trial,:) = sum(real(R(SI>quantile(SI,.4),:)),1);
%     FR_sum_nss(k,j,trial,:) = sum(real(R(SI<=quantile(SI,.4),:)),1);
%     FR_mean_ss(k,j,trial,:) = nanmean(real(R(SI>quantile(SI,.4),:)),1);
%     FR_mean_nss(k,j,trial,:) = nanmean(real(R(SI<=quantile(SI,.4),:)),1);
%     
%     
%     COVss(k,j,trial) = mean2(COVest(SI>quantile(SI,.4),SI>quantile(SI,.4)));
%     COVnss(k,j,trial) = mean2(COVest(SI<=quantile(SI,.4),SI<=quantile(SI,.4)));
%     
%     SNRnss(k,j,trial) = nanmean(sum(real(R(SI<=quantile(SI,.4),:)),1))./sqrt(sum(sum(triu(COVest(SI<=quantile(SI,.4),SI<=quantile(SI,.4))))));
%     SNRss(k,j,trial) = nanmean(sum(real(R(SI>quantile(SI,.4),:)),1))./sqrt(sum(sum(triu(COVest(SI>quantile(SI,.4),SI>quantile(SI,.4))))));
%     SNR(k,j,trial) = nanmean(sum(real(R),1))./sqrt(sum(sum(triu(COVest))));
%   
%     
%     
%     Rmean{k} = nanmean(R,2);
%     
%     clear mtpopulation R COV TargetEstimate
%     
% end
% 
% end
% end
% 
% SPD(abs(SPD)>100) = nan;
% SPDvar = squeeze(nanstd(SPD(:,:,:,:),[],4));
% SPDvar_std = nanstd(SPDvar,[],2);
% SPDvar_mean = nanmedian(SPDvar,2);
% 
% spemSI(iter) = min(SPDvar_mean)./SPDvar_mean(end);
% 
% end

%% Random Target Speed for any trial
Target.MotionSpeed        =   10 * rand(10000,1) + 10;
Target.MotionDirection    =   40 * ones(10000,1) + 10;

for trial = 1:100



% CorrSPD = nan(length(max(popSize)./popSize),length([0,0.09,0.18,0.27,0.36,0.45]));
% CorrDIR= nan(length(max(popSize)./popSize),length([0,0.09,0.18,0.27,0.36,0.45]));
% rmax = [0,0.09,0.18,0.27,0.36,0.45];
rmax = .25;%[0,0.09,0.18,0.27];
tauD = .25;
tauS = .5;

for j = 1:length(rmax);
    r = rmax(j);
    
    
    %populationRatio = max(popSize)./popSize;
    populationRatio = popSize;
for k = 1:length(populationRatio)
    i = populationRatio(k);
    
    Target.Size = sizes(k) * ones(10000,1);
    param(1) = 3;
    param(2) = 2;
    param(3) = .5;
    
    nps(k) = round(sqrt(i));%floor(46/(1*i));
    npd(k) = round(sqrt(i));%floor(90/(1*i));
    [mtpopulation, R, COV] = MTpopulation(Target,nps(k),npd(k),rmax,tauD,tauS);
    [mtpopulation_denom, R_denom, COV] = MTpopulation(Target,nps(k),npd(k),rmax,tauD,tauS);
    [TargetEstimate] = DecodeMTpopulation(mtpopulation, R, mtpopulation_denom, R_denom,'weighted_uncorr_norm',param);%     Cs = corrcoef(TargetEstimate.SPDest,Target.MotionSpeed');
%     Cd = corrcoef(TargetEstimate.DIRest,Target.MotionDirection');
    
%     mdl = LinearModel.fit(TargetEstimate.DIRest,Target.MotionDirection');
%     biasDIR(k,j,trial) = mdl.Coefficients.Estimate(1);
%     sensitivityDIR(k,j,trial) = mdl.Coefficients.Estimate(2);
%     RDIR(k,j,trial) = mdl.Rsquared.Ordinary;
    TargetEstimate.SPDest(abs(TargetEstimate.SPDest)>100) = nan;
    mdl = LinearModel.fit(Target.MotionSpeed',TargetEstimate.SPDest);
%     biasSPD(k,j,trial) = mdl.Coefficients.Estimate(1);
%     sensitivitySPD(k,j,trial) = mdl.Coefficients.Estimate(2);
%     RSPD(k,j,trial) = mdl.Rsquared.Adjusted;
    RMSE(k,j,trial) = mdl.RMSE;
    
    
%     estimateSPD(k,j,trial,:) = TargetEstimate.SPDest;
%     estimateDIR(k,j,trial,:) = TargetEstimate.DIRest;
%     CorrSPD(k,j,trial) = Cs(1,2);
%     CorrDIR(k,j,trial) = Cd(1,2);
    clear mtpopulation R COV TargetEstimate
    
end

end
end
%%

figure;plot(2*sizes,median(RMSE,3),'.k','MarkerSize',25);
hold on;ploterr(2*sizes, median(RMSE,3), [], std(RMSE,[],3)./sqrt(100),'-k','abshhxy', .3);
xlabel('size (degree)');ylabel('RMSE (degree)');legend(num2str(0),num2str(0.09),num2str(0.18),num2str(0.27),num2str(0.36),num2str(0.45))
% 
% figure;plot(2*sizes,median(sensitivitySPD,3),'.k','MarkerSize',25);
% hold on;ploterr(2*sizes, median(sensitivitySPD,3), [], std(RMSE,[],3)./sqrt(1000),'-k','abshhxy', .3);
% xlabel('size (degree)');ylabel('sensitivity');legend(num2str(0),num2str(0.09),num2str(0.18),num2str(0.27),num2str(0.36),num2str(0.45))


% figure;plot(nps.*npd,median(biasSPD,3),'-o','LineWidth',2);
% xlabel('number of neurons');ylabel('Speed Bias');legend(num2str(0),num2str(0.09),num2str(0.18),num2str(0.27),num2str(0.36),num2str(0.45))
% 
% figure;plot(nps.*npd,median(RSPD,3),'-o','LineWidth',2);
% xlabel('number of neurons');ylabel('Speed R^2');legend(num2str(0),num2str(0.09),num2str(0.18),num2str(0.27),num2str(0.36),num2str(0.45))

%% estimate Vector Averaging parameters

for trial = 1:1000
    param(1) = 3;
    param(2) = 100*rand;
    sampleP(1,trial) = param(2);
    param(3) = 3*rand;
    sampleP(2,trial) = param(3);


% CorrSPD = nan(length(max(popSize)./popSize),length([0,0.09,0.18,0.27,0.36,0.45]));
% CorrDIR= nan(length(max(popSize)./popSize),length([0,0.09,0.18,0.27,0.36,0.45]));
% rmax = [0,0.09,0.18,0.27,0.36,0.45];
rmax = .25;%[0,0.09,0.18,0.27];
tauD = .25;
tauS = .5;

for j = 1:length(rmax);
    r = rmax(j);
    
    
    %populationRatio = max(popSize)./popSize;
    populationRatio = popSize;
for k = 1:length(populationRatio)
    i = populationRatio(k);
    
    Target.Size = sizes(k) * ones(10000,1);
    
    
    nps(k) = round(sqrt(i));%floor(46/(1*i));
    npd(k) = round(sqrt(i));%floor(90/(1*i));
    [mtpopulation, R, COV] = MTpopulation(Target,nps(k),npd(k),rmax,tauD,tauS);
    [mtpopulation_denom, R_denom, COV] = MTpopulation(Target,nps(k),npd(k),rmax,tauD,tauS);
    [TargetEstimate] = DecodeMTpopulation(mtpopulation, R, mtpopulation_denom, R_denom,'weighted_uncorr_norm',param);%     Cs = corrcoef(TargetEstimate.SPDest,Target.MotionSpeed');
%     Cd = corrcoef(TargetEstimate.DIRest,Target.MotionDirection');
    
%     mdl = LinearModel.fit(TargetEstimate.DIRest,Target.MotionDirection');
%     biasDIR(k,j,trial) = mdl.Coefficients.Estimate(1);
%     sensitivityDIR(k,j,trial) = mdl.Coefficients.Estimate(2);
%     RDIR(k,j,trial) = mdl.Rsquared.Ordinary;
    TargetEstimate.SPDest(abs(TargetEstimate.SPDest)>100) = nan;
    mdl = LinearModel.fit(Target.MotionSpeed',TargetEstimate.SPDest);
%     biasSPD(k,j,trial) = mdl.Coefficients.Estimate(1);
%     sensitivitySPD(k,j,trial) = mdl.Coefficients.Estimate(2);
%     RSPD(k,j,trial) = mdl.Rsquared.Adjusted;
    RMSEsample(k,j,trial) = mdl.RMSE;
    
    
%     estimateSPD(k,j,trial,:) = TargetEstimate.SPDest;
%     estimateDIR(k,j,trial,:) = TargetEstimate.DIRest;
%     CorrSPD(k,j,trial) = Cs(1,2);
%     CorrDIR(k,j,trial) = Cd(1,2);
    clear mtpopulation R COV TargetEstimate
    
end

end
end

