%% Number of neurons as a function of stimulus size

sizes=sqrt(0)+(.2:2:10.2); % the avg eccentricity of the stimuli is sqrt(5) deg
popSize=zeros(1,5);
errSize=zeros(1,5); % reviewer #2 asked us to calculate the error of our estimates
% fun=@(x1) (1.14*x1.^-0.76); % Albright & Desimone 87 Exp Brain Res (mm/deg)
fun=@(x1) (6*x1.^-0.9); % Erickson (mm/deg)
for i=1:length(sizes) % calculate the # of neurons for each size
    %     popSize(i)=integral(fun ,0.5,sizes(i))+1; % integral from 0.5 to edge, between 0 to 0.5 is roughly 1 mm (2mm/deg)
    popSize(i)=integral(fun,12 - sizes(i),12 + sizes(i)); % integral from 1 to edge, between 0 to 1 is roughly (9mm/deg) erickson et al. EBR
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


%% Simulate the population activity
clc;
TauD = (.1:.2:5)*.125;%.5;% 
TauS = (.1:.2:5)*.15;%.5;%
% TauD = .2;
% TauS = .5;
% Rmax = .1:.01:.5;
[meshTauD,meshTauS] = meshgrid(TauD,TauS);

for iter = 1:size(meshTauD,1)*size(meshTauD,2)%1%1::length(Rmax)%

tauD = meshTauD(iter);
tauS = meshTauS(iter);
% tauD = TauD;
% tauS = TauS;
% r0 = .0;
% b = 13.4;

Target.MotionSpeed        =   15 * ones(1000,1);
Target.MotionDirection    =   40 * ones(1000,1);

rmax = .25;%Rmax(iter);%[0,0.09,0.18,0.27]

i = popSize(end);

Target.Size = sizes(end) * ones(1000,1);
nps = round(sqrt(i));%floor(46/(1*i));
npd = round(sqrt(i));%floor(90/(1*i));

[mtpopulation, R, COV] = MTpopulation(Target,nps,npd,rmax,tauD,tauS);
[mtpopulation_denom, R_denom, COV] = MTpopulation(Target,nps,npd,rmax,tauD,tauS);

[TargetEstimate] = DecodeMTpopulation(mtpopulation, (R), mtpopulation_denom, R_denom, 'weighted_uncorr_norm');
DIRstd = nanstd(TargetEstimate.DIRest);
SPDstd = nanstd(TargetEstimate.SPDest);
DIRm = nanmedian(TargetEstimate.DIRest);
SPDm = nanmedian(TargetEstimate.SPDest);

SI = cellfun(@(x)(x.SuppressionIndex),mtpopulation)';
surrSuppx=false(1,length(SI)); surrSuppy=surrSuppx;
surrSuppx(SI>median(SI))= true; surrSuppy(SI>median(SI))=1; % surround suppressed or not
[x,y] = meshgrid(surrSuppx,surrSuppy);
surrSuppMat=x+y; surrSuppMat(logical(eye(size(surrSuppMat))))=3; % categorize to SS-SS, SS-NS, and NS-NS
surrSuppRep = repmat(surrSuppx',1,size(R,2));

RmeanSS = mean(R(surrSuppRep));RstdSS = std(R(surrSuppRep));
RmeanNSS = mean(R(~surrSuppRep));RstdNSS = std(R(~surrSuppRep));
% figure;errorbar([RmeanSS,RmeanNSS],[RstdSS/2,RstdNSS/2],'--x');title('Firing Rate SS and nSS neurons')

PD = cellfun(@(x)(x.PreferredDirection),mtpopulation)';
PS = cellfun(@(x)(x.PreferredSpeed),mtpopulation)';
PreferredProp = [PD,PS];

Theta = 0:359;%Theta = Theta * pi / 180;
S = .5:256;
DTW     =       40 * pi/180;
STW     =       1.5;
% SignalCorr = zeros(length(mtpopulation),length(mtpopulation));
for i = 1:length(mtpopulation)

    
    PD_i = PD(i);
    PS_i = PS(i);

%     DDiff1 = ((Theta - PD_i) > pi) .* (-pi + ((Theta - PD_i) - pi));
%     DDiff2 = ((Theta - PD_i) < -pi) .* pi - (-pi - (Theta - PD_i));
%     DDiff3 = ((Theta - PD_i) <= pi & (Theta - PD_i) >= -pi) .* (Theta -
%     PD_i);
%     
%     DDiff = DDiff1 + DDiff2 + DDiff3;
    DDiff = AngDiff(Theta,PD_i) * pi/180;
    
    DirectionTune(i,:) = exp(-( DDiff.^2 )/( 2 * DTW^2 ));
    SpeedTune(i,:) = exp(-( (log2(S) - (PS_i)).^2 )./( 2 * STW^2 ));

end


NoiseCorr = corrcoef(R');
SignalCorrSPD = corrcoef(SpeedTune');
SignalCorrDIR = corrcoef(DirectionTune');
SignalCorr = SignalCorrSPD.*SignalCorrDIR;

simCorr = nan((length(mtpopulation))*(length(mtpopulation) - 1)./2,4);
counter = 0;
for i = 1:length(mtpopulation)
    for j = i:length(mtpopulation)
        counter = counter + 1;
        simCorr(counter,1) = SignalCorrDIR(i,j); % signal correlation
        simCorr(counter,2) = NoiseCorr(i,j); % noise correlation
        simCorr(counter,3) = SI(i); % suppression index first neuron
        simCorr(counter,4) = SI(j); % suppresion index secodn neuron
    end
end


sigcorrSS = SignalCorrDIR(surrSuppMat == 2);
noisecorrSS = NoiseCorr(surrSuppMat == 2);
sigcorrNSS = SignalCorrDIR(surrSuppMat == 0);
noisecorrNSS = NoiseCorr(surrSuppMat == 0);
sigcorrSSnSS = SignalCorrDIR(surrSuppMat == 1);
noisecorrSSnSS = NoiseCorr(surrSuppMat == 1);

CorrelationMatrix_ss = corrcoef(noisecorrSS,sigcorrSS);
CorrelationMatrix_nss = corrcoef(noisecorrNSS,sigcorrNSS);
CorrelationMatrix_ssnss = corrcoef(noisecorrSSnSS,sigcorrSSnSS);


mdlSS = LinearModel.fit(sigcorrSS,noisecorrSS);
noisPredSS = predict(mdlSS,[min(sigcorrSS):.01:max(sigcorrSS)]');
mdlNSS = LinearModel.fit(sigcorrNSS,noisecorrNSS);
noisPredNSS = predict(mdlNSS,[min(sigcorrNSS):.01:max(sigcorrNSS)]');
mdlSSnSS = LinearModel.fit(sigcorrSSnSS,noisecorrSSnSS);
noisPredSSnSS = predict(mdlSSnSS,[min(sigcorrSSnSS):.01:max(sigcorrSSnSS)]');

coeffSS(iter) = mdlSS.Coefficients.Estimate(2);
intrcptSS = mdlSS.Coefficients.Estimate(1);
coeffNSS(iter) = mdlNSS.Coefficients.Estimate(2);
intrcptNSS = mdlNSS.Coefficients.Estimate(1);
coeffSSnSS(iter) = mdlSSnSS.Coefficients.Estimate(2);
intrcptSSnSS = mdlSSnSS.Coefficients.Estimate(1);


% fprintf('coeff_{SS} = %1.6f and coeff_{nSS} = %1.6f and coeff_{SSnSS} = %1.6f \n intercept_{SS} = %1.6f and intercept_{nSS} = %1.6f and intercept_{SSnSS} = %1.6f \n',coeffSS,coeffNSS,coeffSSnSS,intrcptSS,intrcptNSS,intrcptSSnSS); 
% figure(2);
% plot(SignalCorrDIR(surrSuppMat == 2),NoiseCorr(surrSuppMat == 2),'.','Color',[0,0,1]);
% % grid on;
% hold on;
% plot(SignalCorrDIR(surrSuppMat == 0),NoiseCorr(surrSuppMat == 0),'.','Color',[1,0,0]);
% plot(SignalCorrDIR(surrSuppMat == 1),NoiseCorr(surrSuppMat == 1),'.','Color',[0,0,0]);

% plot([min(sigcorrSS):.01:max(sigcorrSS)],noisPredSS,'k','LineWidth',3)
% plot([min(sigcorrNSS):.01:max(sigcorrNSS)],noisPredNSS,'r','LineWidth',3)
% legend('with SS','no SS', 'SS and nSS')

[IDX,Distance] = knnsearch(simCorr,CorrMatrix,'K',5);
meanDistance(iter) = 1./mean2(max(Distance,[],2));


end
% figure;surf(meshTauD,meshTauS,reshape(meanDistance,sqrt(length(meanDistance)),sqrt(length(meanDistance))));

%% Effect of SS on Firirng Rate

% tauD = 1*.125;%
% tauS = 1*.15;%
% 
% 
% Target.MotionSpeed        =   [15 * ones(500,1);15 * ones(500,1);...
%     15 * ones(500,1); 15 * ones(500,1); ...
%     15 * ones(500,1)];
% Target.MotionDirection    =   [40 * ones(500,1);40 * ones(500,1);...
%     40 * ones(500,1); 40 * ones(500,1); ...
%     40 * ones(500,1)];
% 
% rmax = .27;%[0,0.09,0.18,0.27]
% 
% i = 1;
% 
% Target.Size = [sizes(end) * ones(500,1);sizes(end-1) * ones(500,1);...
%     sizes(end-2) * ones(500,1); sizes(end-3) * ones(500,1); ...
%     sizes(end-4) * ones(500,1)];
% nps = floor(46/(1*i));
% npd = floor(90/(1*i));
% 
% [mtpopulation, R, COV] = MTpopulation(Target,nps,npd,rmax,tauD,tauS);
% 
% SI = cellfun(@(x)(x.SuppressionIndex),mtpopulation)';
% surrSuppx=false(1,length(SI)); surrSuppy=surrSuppx;
% surrSuppx(SI>median(SI))= true; surrSuppy(SI>median(SI))=1; % surround suppressed or not
% 
% 
% figure;plot(1,nanmean(nanmean(R(surrSuppx,1:500),2)),'ok')
% hold on;plot(2,nanmean(nanmean(R(surrSuppx,501:1000),2)),'-ok')
% plot(3,nanmean(nanmean(R(surrSuppx,1001:1500),2)),'-ok')
% plot(4,nanmean(nanmean(R(surrSuppx,1501:2000),2)),'-ok')
% plot(5,nanmean(nanmean(R(surrSuppx,2001:2500),2)),'-ok')
% plot(1,nanmean(nanmean(R(~surrSuppx,1:500),2)),'-or')
% hold on;plot(2,nanmean(nanmean(R(~surrSuppx,501:1000),2)),'-or')
% plot(3,nanmean(nanmean(R(~surrSuppx,1001:1500),2)),'-or')
% plot(4,nanmean(nanmean(R(~surrSuppx,1501:2000),2)),'-or')
% plot(5,nanmean(nanmean(R(~surrSuppx,2001:2500),2)),'-or')
% 
% 
% %% Try find the best tauD and tauS
% 
% clc;
% TauD = .1:.05:.5;%.125
% TauS = .1:.05:.5;%.15
% 
% 
% Target.MotionSpeed        =   15 * ones(1000,1);
% Target.MotionDirection    =   40 * ones(1000,1);
% 
% rmax = .45;%[0,0.09,0.18,0.27]
% 
% i = 1;
% 
% Target.Size = sizes(end) * ones(1000,1);
% nps = floor(46/(1*i));
% npd = floor(90/(1*i));
% 
% for k = 1:length(TauD)
%     for j = 1:length(TauS)
%         tauD = TauD(k);
%         tauS = TauS(j);
%         [mtpopulation, R, COV] = MTpopulation(Target,nps,npd,rmax,tauD,tauS);
%         % [mtpopulation_denom, R_denom, COV] = MTpopulation(Target,nps,npd,rmax,tauD,tauS);
%         
%         % [TargetEstimate] = DecodeMTpopulation(mtpopulation, (R), R_denom, 'uncorr_norm');
%         % DIRstd = nanstd(TargetEstimate.DIRest);
%         % SPDstd = nanstd(TargetEstimate.SPDest);
%         % DIRm = nanmedian(TargetEstimate.DIRest);
%         % SPDm = nanmedian(TargetEstimate.SPDest);
%         
%         SI = cellfun(@(x)(x.SuppressionIndex),mtpopulation)';
%         surrSuppx=false(1,length(SI)); surrSuppy=surrSuppx;
%         surrSuppx(SI>median(SI))= true; surrSuppy(SI>median(SI))=1; % surround suppressed or not
%         [x,y] = meshgrid(surrSuppx,surrSuppy);
%         surrSuppMat=x+y; surrSuppMat(logical(eye(size(surrSuppMat))))=3; % categorize to SS-SS, SS-NS, and NS-NS
%         surrSuppRep = repmat(surrSuppx',1,size(R,2));
%         
%         % RmeanSS = mean(R(surrSuppRep));RstdSS = std(R(surrSuppRep));
%         % RmeanNSS = mean(R(~surrSuppRep));RstdNSS = std(R(~surrSuppRep));
%         % figure;errorbar([RmeanSS,RmeanNSS],[RstdSS/2,RstdNSS/2],'--x');title('Firing Rate SS and nSS neurons')
%         
%         PD = cellfun(@(x)(x.PreferredDirection),mtpopulation)';
%         PS = cellfun(@(x)(x.PreferredSpeed),mtpopulation)';
%         PreferredProp = [PD,PS];
%         
%         Theta = 0:359;%Theta = Theta * pi / 180;
%         S = .5:256;
%         DTW     =       40 * pi/180;
%         STW     =       1.5;
%         % SignalCorr = zeros(length(mtpopulation),length(mtpopulation));
%         for i = 1:length(mtpopulation)
%             
%             
%             PD_i = PD(i);
%             PS_i = PS(i);
%             
%             %     DDiff1 = ((Theta - PD_i) > pi) .* (-pi + ((Theta - PD_i) - pi));
%             %     DDiff2 = ((Theta - PD_i) < -pi) .* pi - (-pi - (Theta - PD_i));
%             %     DDiff3 = ((Theta - PD_i) <= pi & (Theta - PD_i) >= -pi) .* (Theta - PD_i);
%             %
%             %     DDiff = DDiff1 + DDiff2 + DDiff3;
%             DDiff = AngDiff(Theta,PD_i) * pi/180;
%             
%             DirectionTune(i,:) = exp(-( DDiff.^2 )/( 2 * DTW^2 ));
%             SpeedTune(i,:) = exp(-( (log2(S) - (PS_i)).^2 )./( 2 * STW^2 ));
%             
%         end
%         
%         
%         NoiseCorr = corrcoef(R');
%         SignalCorrSPD = corrcoef(SpeedTune');
%         SignalCorrDIR = corrcoef(DirectionTune');
%         SignalCorr = SignalCorrSPD.*SignalCorrDIR;
%         
%         
%         
%         sigcorrSS = SignalCorrSPD(surrSuppMat == 2);
%         noisecorrSS = NoiseCorr(surrSuppMat == 2);
%         sigcorrNSS = SignalCorrSPD(surrSuppMat == 0);
%         noisecorrNSS = NoiseCorr(surrSuppMat == 0);
%         
%         CorrelationMatrix_ss = corrcoef(noisecorrSS,sigcorrSS);
%         CorrelationMatrix_nss = corrcoef(noisecorrNSS,sigcorrNSS);
%         
%         
%         mdlSS = LinearModel.fit(sigcorrSS,noisecorrSS);
%         noisPredSS = predict(mdlSS,[min(sigcorrSS):.01:max(sigcorrSS)]');
%         mdlNSS = LinearModel.fit(sigcorrNSS,noisecorrNSS);
%         noisPredNSS = predict(mdlNSS,[min(sigcorrNSS):.01:max(sigcorrNSS)]');
%         
%         coeffSS(k,j) = mdlSS.Coefficients.Estimate(2);
%         intrcptSS(k,j) = mdlSS.Coefficients.Estimate(1);
%         coeffNSS(k,j) = mdlNSS.Coefficients.Estimate(2);
%         intrcptNSS(k,j) = mdlNSS.Coefficients.Estimate(1);
%         
%     end
% end
% 
