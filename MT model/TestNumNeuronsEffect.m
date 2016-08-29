%% Number of neurons as a function of stimulus size

sizes=sqrt(0)+(1:2:21); % the avg eccentricity of the stimuli is sqrt(5) deg
popSize=zeros(1,5);
errSize=zeros(1,5); % reviewer #2 asked us to calculate the error of our estimates
% fun=@(x1) (1.14*x1.^-0.76); % Albright & Desimone 87 Exp Brain Res (mm/deg)
fun=@(x1) (6*x1.^-0.9); % Erickson (mm/deg)
for i=1:length(sizes) % calculate the # of  neurons for each size
%     popSize(i)=integral(fun,0.5,sizes(i))+1; % integral from 0.5 to edge, between 0 to 0.5 is roughly 1 mm (2mm/deg)
    popSize(i)=integral(fun,22 - sizes(i),22 + sizes(i)); % integral from 1 to edge, between 0 to 1 is roughly (9mm/deg) erickson et al. EBR
%     popSize(i)=integral(fun,.5,sizes(i)) + 9;

    errSize(i)=integral(fun,0.5,sizes(i)*0.16+0.51)+1; % integral of the error
end
% popSize=round(1*popSize.^2); % multiple by a constant of 10, but the constant did not affect the trend
popSize=round(20*popSize); % new constant`
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
% tauD = 1.5*.125;%0.1125;%0.1125;%
% tauS = 1.5*.15;%0.135;%0.1350;%
tauD = 1.0*.125;
tauS = 1.0*.15;
% Target.MotionSpeed        =   [15 * ones(500,1);15 * ones(500,1);...
%     15 * ones(500,1); 15 * ones(500,1); ...
%     15 * ones(500,1)];
% Target.MotionDirection    =   [40 * ones(500,1);40 * ones(500,1);...
%     40 * ones(500,1); 40 * ones(500,1); ...
%     40 * ones(500,1)];
% 
% Target.Size = [sizes(end-4) * ones(500,1);sizes(end-3) * ones(500,1);...
%     sizes(end-2) * ones(500,1); sizes(end-1) * ones(500,1); ...
%     sizes(end) * ones(500,1)];

for trial = 1:1000
Target.MotionSpeed = 15 * ones(1000,1);
Target.MotionDirection = 40 * ones(1000,1);

j = 0;
% DIRstd = nan(length([1,1.5,2,3,4]),length([0,0.09,0.18,0.27,0.36,0.45]));
% DIRm = nan(length([1,1.5,2,3,4]),length([0,0.09,0.18,0.27,0.36,0.45]));
% SPDstd = nan(length([1,1.5,2,3,4]),length([0,0.09,0.18,0.27,0.36,0.45]));
% SPDm = nan(length([1,1.5,2,3,4]),length([0,0.09,0.18,0.27,0.36,0.45]));
for rmax = [200]
    
    
    j = j + 1;
    k = 0;
for i = popSize%max(popSize)./popSize%[1,1.5,2,3,4]%log2(2:(62/5):64-62/5)
%  i = 1;    
    k = k + 1;
    Target.Size = sizes(k) * ones(1000,1);
    nps(k) = round(sqrt(i));%floor(46/(1*i));
    npd(k) = round(sqrt(i));%floor(90/(1*i));
    
    [mtpopulation, R, COV] = MTpopulation(Target,nps(k),npd(k),rmax,tauD,tauS);
    [mtpopulation_denom, R_denom, COV] = MTpopulation(Target,nps(k),npd(k),rmax,tauD,tauS);
    [TargetEstimate] = DecodeMTpopulation(mtpopulation, R, R_denom,'uncorr_norm');
%     [TargetEstimate] = DecodeMTpopulation(mtpopulation, real(R));
%     DIRstd(k,j,trial) = nanstd(TargetEstimate.DIRest);
%     SPDstd(k,j,trial) = nanstd(TargetEstimate.SPDest);
%     DIRm(k,j,trial) = nanmedian(TargetEstimate.DIRest);
%     SPDm(k,j,trial) = nanmedian(TargetEstimate.SPDest);
    SPD(k,j,trial,:) = TargetEstimate.SPDest;
    
    
    
%     DIRstd(j,trial,:) = [nanstd(TargetEstimate.DIRest(1:500)),...
%         nanstd(TargetEstimate.DIRest(501:1000)),...
%         nanstd(TargetEstimate.DIRest(1001:1500)),...
%         nanstd(TargetEstimate.DIRest(1501:2000)),...
%         nanstd(TargetEstimate.DIRest(2001:2500))];
%     SPDstd(j,trial,:) = [nanstd(TargetEstimate.SPDest(1:500)),...
%         nanstd(TargetEstimate.SPDest(501:1000)),...
%         nanstd(TargetEstimate.SPDest(1001:1500)),...
%         nanstd(TargetEstimate.SPDest(1501:2000)),...
%         nanstd(TargetEstimate.SPDest(2001:2500))];
    
    
    
    Rmean{k} = nanmean(R,2);
    
    clear mtpopulation R COV TargetEstimate
    
end

end
end
%% Plots - Speed and Direction variance as a function of #neurons
% 

SPD(abs(SPD)>100) = nan;
SPDvar = squeeze(nanstd(SPD(:,:,:,:),[],4));
SPDvar_std = std(SPDvar,[],2);
SPDvar_mean = nanmedian(SPDvar,2);
hh = ploterr(sizes, SPDvar_mean, [], SPDvar_std./sqrt(1000));
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

%% Random Target Speed for any trial
% for trial = 1:100
% Target.MotionSpeed        =   10 * rand(10000,1) + 10;
% Target.MotionDirection    =   40 * ones(10000,1) + 10;
% 
% 
% % CorrSPD = nan(length(max(popSize)./popSize),length([0,0.09,0.18,0.27,0.36,0.45]));
% % CorrDIR= nan(length(max(popSize)./popSize),length([0,0.09,0.18,0.27,0.36,0.45]));
% % rmax = [0,0.09,0.18,0.27,0.36,0.45];
% rmax = 250;%[0,0.09,0.18,0.27];
% tauD = 0.125;%0.2;%
% tauS = 0.15;%0.24;%
% for j = 1:length(rmax);
%     r = rmax(j);
%     
%     
%     %populationRatio = max(popSize)./popSize;
%     populationRatio = popSize;
% for k = 1:length(populationRatio)
%     i = populationRatio(k);
%     
%     Target.Size = sizes(k) * ones(10000,1);
%     
%     
%     nps(k) = round(sqrt(i));%floor(46/(1*i));
%     npd(k) = round(sqrt(i));%floor(90/(1*i));
%     [mtpopulation, R, COV] = MTpopulation(Target,nps(k),npd(k),rmax,tauD,tauS);
%     [mtpopulation_denom, R_denom, COV] = MTpopulation(Target,nps(k),npd(k),rmax,tauD,tauS);
%     [TargetEstimate] = DecodeMTpopulation(mtpopulation, R, R_denom,'num&denom_opp');%     Cs = corrcoef(TargetEstimate.SPDest,Target.MotionSpeed');
% %     Cd = corrcoef(TargetEstimate.DIRest,Target.MotionDirection');
%     
% %     mdl = LinearModel.fit(TargetEstimate.DIRest,Target.MotionDirection');
% %     biasDIR(k,j,trial) = mdl.Coefficients.Estimate(1);
% %     sensitivityDIR(k,j,trial) = mdl.Coefficients.Estimate(2);
% %     RDIR(k,j,trial) = mdl.Rsquared.Ordinary;
%     TargetEstimate.SPDest(abs(TargetEstimate.SPDest)>100) = nan;
%     mdl = LinearModel.fit(Target.MotionSpeed',TargetEstimate.SPDest);
% %     biasSPD(k,j,trial) = mdl.Coefficients.Estimate(1);
% %     sensitivitySPD(k,j,trial) = mdl.Coefficients.Estimate(2);
%     RSPD(k,j,trial) = mdl.Rsquared.Adjusted;
%     
%     
%     estimateSPD(k,j,trial,:) = TargetEstimate.SPDest;
%     estimateDIR(k,j,trial,:) = TargetEstimate.DIRest;
% %     CorrSPD(k,j,trial) = Cs(1,2);
% %     CorrDIR(k,j,trial) = Cd(1,2);
%     clear mtpopulation R COV TargetEstimate
%     
% end
% 
% end
% end
%%

% figure;plot(nps.*npd,median(sensitivitySPD,3),'-o','LineWidth',2);
% xlabel('number of neurons');ylabel('Speed Sensitivity');legend(num2str(0),num2str(0.09),num2str(0.18),num2str(0.27),num2str(0.36),num2str(0.45))
% 
% figure;plot(nps.*npd,median(biasSPD,3),'-o','LineWidth',2);
% xlabel('number of neurons');ylabel('Speed Bias');legend(num2str(0),num2str(0.09),num2str(0.18),num2str(0.27),num2str(0.36),num2str(0.45))
% 
% figure;plot(nps.*npd,median(RSPD,3),'-o','LineWidth',2);
% xlabel('number of neurons');ylabel('Speed R^2');legend(num2str(0),num2str(0.09),num2str(0.18),num2str(0.27),num2str(0.36),num2str(0.45))

%% Effect of tauD and tauS on the variance 
% 
% for trial = 1:10
% Target.MotionSpeed        =   15 * ones(500,1);
% Target.MotionDirection    =   0 * ones(500,1);
% 
% j = 0;
% % DIRstd = nan(length([1,1.5,2,3,4]),length([0,0.09,0.18,0.27,0.36,0.45]));
% % DIRm = nan(length([1,1.5,2,3,4]),length([0,0.09,0.18,0.27,0.36,0.45]));
% % SPDstd = nan(length([1,1.5,2,3,4]),length([0,0.09,0.18,0.27,0.36,0.45]));
% % SPDm = nan(length([1,1.5,2,3,4]),length([0,0.09,0.18,0.27,0.36,0.45]));
% TauD        =   0.01:0.05:1;
% TauS        =   0.01:0.05:1;
% Tau         =   [TauD;TauS];
% rmax        =   0.18;
% 
% for taucount = 1:size(Tau,2)
%     tauD = Tau(1,taucount);
%     tauS = Tau(2,taucount);
%     j = j + 1;
%     k = 0;
% for i = max(popSize)./popSize(end)%[1,1.5,2,3,4]%log2(2:(62/5):64-62/5)
%     k = k + 1;
%     Target.Size = sizes(k) * ones(500,1);
%     nps(k) = floor(46/(1*i));
%     npd(k) = floor(90/(1*i));
%     
%     [mtpopulation, R, COV] = MTpopulation(Target,nps(k),npd(k),rmax,tauD,tauS);
%     [TargetEstimate] = DecodeMTpopulation(mtpopulation, real(R));
%     DIRstd(k,j,trial) = std(TargetEstimate.DIRest);
%     SPDstd(k,j,trial) = std(TargetEstimate.SPDest);
%     DIRm(k,j,trial) = median(TargetEstimate.DIRest);
%     SPDm(k,j,trial) = median(TargetEstimate.SPDest);
%     
%     clear mtpopulation R COV TargetEstimate
%     
% end
% 
% end
% end

%% plot the results of TauD and TauS changes
% figure;plot(nps.*npd,median(SPDstd,3),'-o','LineWidth',2);
% xlabel('number of neurons');ylabel('speed variance')
% figure;plot(nps.*npd,median(DIRstd,3),'-o','LineWidth',2);
% xlabel('number of neurons');ylabel('direction variance')
% 
% figure;plot(nps.*npd,median(SPDm,3),'-o','LineWidth',2);
% xlabel('number of neurons');ylabel('speed mean')
% figure;plot(nps.*npd,median(DIRm,3),'-o','LineWidth',2);
% xlabel('number of neurons');ylabel('direction mean')