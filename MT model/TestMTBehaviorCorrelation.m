%% Number of neurons as a function of stimulus size

sizes=sqrt(0)+(1:2:11); % the avg eccentricity of the stimuli is sqrt(5) deg
popSize=zeros(1,5);
errSize=zeros(1,5); % reviewer #2 asked us to calculate the error of our estimates
% fun=@(x1) (1.14*x1.^-0.76); % Albright & Desimone 87 Exp Brain Res (mm/deg)
fun=@(x1) (6*x1.^-0.9); % Erickson (mm/deg)
for i=1:length(sizes) % calculate the # of neurons for each size
    %     popSize(i)=integral(fun,0.5,sizes(i))+1; % integral from 0.5 to edge, between 0 to 0.5 is roughly 1 mm (2mm/deg)
    popSize(i)=integral(fun,12 - sizes(i),12 + sizes(i)); % integral from 1 to edge, between 0 to 1 is roughly (9mm/deg) erickson et al. EBR
    %     popSize(i)=integral(fun,.5,sizes(i)) + 9;
    
    errSize(i)=integral(fun,0.5,sizes(i)*0.16+0.51)+1; % integral of the error
end
% popSize=round(1*popSize.^2); % multiple by a constant of 10, but the constant did not affect the trend
popSize=round(1*popSize); % new constant`
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

%%
tauD = 1.2*.125;%
tauS = 1.2*.15;%


Target.MotionSpeed        =   15 * ones(1000,1);
Target.MotionDirection    =   40 * ones(1000,1);

rmax = 250;%[0,0.09,0.18,0.27]

i = 1.6154;

Target.Size = sizes(end-1) * ones(1000,1);
nps = floor(46/(1*i));
npd = floor(90/(1*i));

[mtpopulation, FR, COV] = MTpopulation(Target,nps,npd,rmax,tauD,tauS);
[~, FR_denom, COV] = MTpopulation(Target,nps,npd,rmax,tauD,tauS);
[TargetEstimate_uncorrnorm] = DecodeMTpopulation(mtpopulation, FR, FR_denom,'uncorr_norm');
[TargetEstimate_numdenomopp] = DecodeMTpopulation(mtpopulation, FR, FR_denom,'num&denom_opp');
[TargetEstimate_numopp] = DecodeMTpopulation(mtpopulation, FR, FR_denom,'num_opp');



DIRest_uncorrnorm = TargetEstimate_uncorrnorm.DIRest;
SPDest_uncorrnorm = TargetEstimate_uncorrnorm.SPDest;

DIRest_numdenomopp = TargetEstimate_numdenomopp.DIRest;
SPDest_numdenomopp = TargetEstimate_numdenomopp.SPDest;

DIRest_numopp = TargetEstimate_numopp.DIRest;
SPDest_numopp = TargetEstimate_numopp.SPDest;




for cellcount = 1:size(FR,1)
    
    thisFRnormalized = (FR(cellcount,:) - nanmean(FR(cellcount,:)))./(nanstd(FR(cellcount,:)));
    
    DIRestnormalized = (DIRest_uncorrnorm - nanmean(DIRest_uncorrnorm))./(nanstd(DIRest_uncorrnorm));
    [c,p] = corrcoef(thisFRnormalized,DIRest_uncorrnorm);
    MTBehaviorCorrDIR_uncorrnorm(cellcount) = c(2);
    MTBehaviorCorrDIR_p_uncorrnorm(cellcount) = p(2);
    
    DIRestnormalized = (DIRest_numdenomopp - nanmean(DIRest_numdenomopp))./(nanstd(DIRest_numdenomopp));
    [c,p] = corrcoef(thisFRnormalized,DIRest_numdenomopp);
    MTBehaviorCorrDIR_numdenomopp(cellcount) = c(2);
    MTBehaviorCorrDIR_p_numdenomopp(cellcount) = p(2);
    
    DIRestnormalized = (DIRest_numopp - nanmean(DIRest_numopp))./(nanstd(DIRest_numopp));
    [c,p] = corrcoef(thisFRnormalized,DIRest_numopp);
    MTBehaviorCorrDIR_numopp(cellcount) = c(2);
    MTBehaviorCorrDIR_p_numopp(cellcount) = p(2);
    
    SPDestnormalized = (SPDest_uncorrnorm - nanmean(SPDest_uncorrnorm))./(nanstd(SPDest_uncorrnorm));
    [c,p] = corrcoef(thisFRnormalized,SPDestnormalized);
    MTBehaviorCorrSPD_uncorrnorm(cellcount) = c(2);
    MTBehaviorCorrSPD_p_uncorrnorm(cellcount) = p(2);
    
    SPDestnormalized = (SPDest_numopp - nanmean(SPDest_numopp))./(nanstd(SPDest_numopp));
    [c,p] = corrcoef(thisFRnormalized,SPDestnormalized);
    MTBehaviorCorrSPD_numopp(cellcount) = c(2);
    MTBehaviorCorrSPD_p_numopp(cellcount) = p(2);
    
    SPDestnormalized = (SPDest_numdenomopp - nanmean(SPDest_numdenomopp))./(nanstd(SPDest_numdenomopp));
    [c,p] = corrcoef(thisFRnormalized,SPDestnormalized);
    MTBehaviorCorrSPD_numdenomopp(cellcount) = c(2);
    MTBehaviorCorrSPD_p_numdenomopp(cellcount) = p(2);
    
end




SI = cellfun(@(x)(x.SuppressionIndex),mtpopulation)';
surrSuppx=false(1,length(SI)); surrSuppy=surrSuppx;
surrSuppx(SI>median(SI))=true; surrSuppy(SI>median(SI))=1; % surround suppressed or not
[x,y] = meshgrid(surrSuppx,surrSuppy);
surrSuppMat=x+y; surrSuppMat(logical(eye(size(surrSuppMat))))=3; % categorize to SS-SS, SS-NS, and NS-NS





% MTBehaviorCorrDIRSSmeanstd_uncorrnorm = [nanmean(MTBehaviorCorrDIR_uncorrnorm(surrSuppx)),nanstd(MTBehaviorCorrDIR_uncorrnorm(surrSuppx))];
MTBehaviorCorrSPDSSmeanstd_uncorrnorm = [nanmean(MTBehaviorCorrSPD_uncorrnorm(surrSuppx)),nanstd(MTBehaviorCorrSPD_uncorrnorm(surrSuppx))];
% MTBehaviorCorrDIRNSSmeanstd_uncorrnorm = [nanmean(MTBehaviorCorrDIR_uncorrnorm(~surrSuppx)),nanstd(MTBehaviorCorrDIR_uncorrnorm(~surrSuppx))];
MTBehaviorCorrSPDNSSmeanstd_uncorrnorm = [nanmean(MTBehaviorCorrSPD_uncorrnorm(~surrSuppx)),nanstd(MTBehaviorCorrSPD_uncorrnorm(~surrSuppx))];

% MTBehaviorCorrDIRSSmeanstd_numdenomopp = [nanmean(MTBehaviorCorrDIR_numdenomopp(surrSuppx)),nanstd(MTBehaviorCorrDIR_numdenomopp(surrSuppx))];
MTBehaviorCorrSPDSSmeanstd_numdenomopp = [nanmean(MTBehaviorCorrSPD_numdenomopp(surrSuppx)),nanstd(MTBehaviorCorrSPD_numdenomopp(surrSuppx))];
% MTBehaviorCorrDIRNSSmeanstd_numdenomopp = [nanmean(MTBehaviorCorrDIR_numdenomopp(~surrSuppx)),nanstd(MTBehaviorCorrDIR_numdenomopp(~surrSuppx))];
MTBehaviorCorrSPDNSSmeanstd_numdenomopp = [nanmean(MTBehaviorCorrSPD_numdenomopp(~surrSuppx)),nanstd(MTBehaviorCorrSPD_numdenomopp(~surrSuppx))];

% MTBehaviorCorrDIRSSmeanstd_numopp = [nanmean(MTBehaviorCorrDIR_numopp(surrSuppx)),nanstd(MTBehaviorCorrDIR_numopp(surrSuppx))];
MTBehaviorCorrSPDSSmeanstd_numopp = [nanmean(MTBehaviorCorrSPD_numopp(surrSuppx)),nanstd(MTBehaviorCorrSPD_numopp(surrSuppx))];
% MTBehaviorCorrDIRNSSmeanstd_numopp = [nanmean(MTBehaviorCorrDIR_numopp(~surrSuppx)),nanstd(MTBehaviorCorrDIR_numopp(~surrSuppx))];
MTBehaviorCorrSPDNSSmeanstd_numopp = [nanmean(MTBehaviorCorrSPD_numopp(~surrSuppx)),nanstd(MTBehaviorCorrSPD_numopp(~surrSuppx))];

% fprintf('MT-behavior correlation for direction for SS neurons: %1.4f +- %1.4f \n',MTBehaviorCorrDIRSSmeanstd(1),MTBehaviorCorrDIRSSmeanstd(2));
% fprintf('MT-behavior correlation for direction for nSS neurons: %1.4f +- %1.4f \n',MTBehaviorCorrDIRNSSmeanstd(1),MTBehaviorCorrDIRNSSmeanstd(2));
% 
% fprintf('MT-behavior correlation for speed for SS neurons: %1.4f +- %1.4f \n',MTBehaviorCorrSPDSSmeanstd(1),MTBehaviorCorrSPDSSmeanstd(2));
% fprintf('MT-behavior correlation for speed for nSS neurons: %1.4f +- %1.4f \n',MTBehaviorCorrSPDNSSmeanstd(1),MTBehaviorCorrSPDNSSmeanstd(2));


figure;subplot(3,1,1);
histogram(MTBehaviorCorrSPD_uncorrnorm(surrSuppx),20);hold on;histogram(MTBehaviorCorrSPD_uncorrnorm(~surrSuppx),20);
legend('with SS','no SS');title('uncorrelated normalization')
subplot(3,1,2);
histogram(MTBehaviorCorrSPD_numdenomopp(surrSuppx),20);hold on;histogram(MTBehaviorCorrSPD_numdenomopp(~surrSuppx),20);
legend('with SS','no SS');title('numerator and denominator opponency')
subplot(3,1,3);
histogram(MTBehaviorCorrSPD_numopp(surrSuppx),20);hold on;histogram(MTBehaviorCorrSPD_numopp(~surrSuppx),20);
legend('with SS','no SS');title('numerator opponency')


figure;subplot(1,3,1);
errorbar([MTBehaviorCorrSPDSSmeanstd_uncorrnorm(1),MTBehaviorCorrSPDNSSmeanstd_uncorrnorm(1)], ...
    [MTBehaviorCorrSPDSSmeanstd_uncorrnorm(2)/2,MTBehaviorCorrSPDNSSmeanstd_uncorrnorm(2)/2],'x');
title('uncorrelated normalization')
subplot(1,3,2);
errorbar([MTBehaviorCorrSPDSSmeanstd_numdenomopp(1),MTBehaviorCorrSPDNSSmeanstd_numdenomopp(1)], ...
    [MTBehaviorCorrSPDSSmeanstd_numdenomopp(2)/2,MTBehaviorCorrSPDNSSmeanstd_numdenomopp(2)/2],'x');
title('numerator and denominator opponency')
subplot(1,3,3);
errorbar([MTBehaviorCorrSPDSSmeanstd_numopp(1),MTBehaviorCorrSPDNSSmeanstd_numopp(1)], ...
    [MTBehaviorCorrSPDSSmeanstd_numopp(2)/2,MTBehaviorCorrSPDNSSmeanstd_numopp(2)/2],'x');
title('numerator opponency')

PD = cellfun(@(x)(x.PreferredDirection),mtpopulation)';
PS = cellfun(@(x)(x.PreferredSpeed),mtpopulation)';

PSdiff = abs(PS' - log2(mean(Target.MotionSpeed)));
figure;
subplot(1,3,1);plot(PSdiff(MTBehaviorCorrSPD_p_uncorrnorm<0.01),MTBehaviorCorrSPD_uncorrnorm(MTBehaviorCorrSPD_p_uncorrnorm<0.01),'.','MarkerSize',15);title('uncorrelated normalization');
subplot(1,3,2);plot(PSdiff(MTBehaviorCorrSPD_p_numdenomopp<0.01),MTBehaviorCorrSPD_numdenomopp(MTBehaviorCorrSPD_p_numdenomopp<0.01),'.','MarkerSize',15);title('numerator and denominator opponency');
subplot(1,3,3);plot(PSdiff(MTBehaviorCorrSPD_p_numopp<0.01),MTBehaviorCorrSPD_numopp(MTBehaviorCorrSPD_p_numopp<0.01),'.','MarkerSize',15);title('numerator opponency');

PDdiff = abs(AngDiff(PD,mean(Target.MotionDirection)));
figure;
subplot(1,3,1);plot(PDdiff(MTBehaviorCorrSPD_p_uncorrnorm<inf),MTBehaviorCorrSPD_uncorrnorm(MTBehaviorCorrSPD_p_uncorrnorm<inf),'.','MarkerSize',15);title('uncorrelated normalization');
subplot(1,3,2);plot(PDdiff(MTBehaviorCorrSPD_p_numdenomopp<inf),MTBehaviorCorrSPD_numdenomopp(MTBehaviorCorrSPD_p_numdenomopp<inf),'.','MarkerSize',15);title('numerator and denominator opponency');
subplot(1,3,3);plot(PDdiff(MTBehaviorCorrSPD_p_numopp<inf),MTBehaviorCorrSPD_numopp(MTBehaviorCorrSPD_p_numopp<inf),'.','MarkerSize',15);title('numerator opponency');
