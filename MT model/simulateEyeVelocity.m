%% Number of neurons as a function of stimulus size

sizes=sqrt(0)+(1:2:9); % the avg eccentricity of the stimuli is sqrt(5) deg
popSize=zeros(1,5);
errSize=zeros(1,5); % reviewer #2 asked us to calculate the error of our estimates
% fun=@(x1) (1.14*x1.^-0.76); % Albright & Desimone 87 Exp Brain Res (mm/deg)
fun=@(x1) (6*x1.^-0.9); % Erickson (mm/deg)
for i=1:length(sizes) % calculate the # of neurons for each size
%     popSize(i)=integral(fun,0.5,sizes(i))+1; % integral from 0.5 to edge, between 0 to 0.5 is roughly 1 mm (2mm/deg)
    popSize(i)=integral(fun,10 - sizes(i),10 + sizes(i)); % integral from 1 to edge, between 0 to 1 is roughly (9mm/deg) erickson et al. EBR
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

gT = 0.5;
gD = 1-gT;
tauD = 0.125;%0.1125;%
tauS = 0.15;%0.1350;%

for trial = 1:40
Target.MotionSpeed              =   10 * rand(1000,1) + 10; % Target is the Dots here
Target.MotionDirection          =   30 * ones(1000,1) + 10;
TargetTarget.MotionSpeed        =   10 * rand(1000,1) + 10; % TargetTarget is the red main target
TargetTarget.MotionDirection    =   30 * ones(1000,1) + 10;
TargetTargetSize                =   sizes(1); % the size of the main target

j = 0;k= 0;
% CorrSPD = nan(length(max(popSize)./popSize),length([0,0.09,0.18,0.27,0.36,0.45]));
% CorrDIR= nan(length(max(popSize)./popSize),length([0,0.09,0.18,0.27,0.36,0.45]));
for rmax = 0.18%[0,0.09,0.18,0.27]
    j = j + 1;
    k = 0;
    
    % -----  simulating the response to main target
    TargetTarget.Size = TargetTargetSize * ones(1000,1);
    npsTarget = floor(46/(1*max(popSize)./popSize(1)));
    npdTarget = floor(90/(1*max(popSize)./popSize(1)));
    [mtpopulationTarget, RTarget, COVTarget] = MTpopulation(TargetTarget,npsTarget,npdTarget,rmax,tauD,tauS);
    [TargetEstimateTarget] = DecodeMTpopulation(mtpopulationTarget, real(RTarget));
    
    estimateSPDTarget(j,trial,:) = TargetEstimateTarget.SPDest;
    estimateDIRTarget(j,trial,:) = TargetEstimateTarget.DIRest;
    
    % -----
    
    gT = max(popSize)./popSize; % pursuit gain for the main target
    gT = (gT - min(gT))./(max(gT)-min(gT));
%     gT = zeros(size(gT));
    
for i = max(popSize)./popSize%log2(2:(62/5):64-62/5)
    k = k + 1;
    Target.Size = sizes(k) * ones(1000,1);
    
    
    npsDots(k) = floor(46/(1*i));
    npdDots(k) = floor(90/(1*i));
    
    
    [mtpopulationDots, RDots, COVDots] = MTpopulation(Target,npsDots(k),npdDots(k),rmax,tauD,tauS);
    [TargetEstimateDots] = DecodeMTpopulation(mtpopulationDots, real(RDots));
    
    estimateSPDDots(k,j,trial,:) = TargetEstimateDots.SPDest; %#ok
    estimateDIRDots(k,j,trial,:) = TargetEstimateDots.DIRest; %#ok
    meanCOV(k,j,trial) = mean2(COVDots); %#ok
    
    clear mtpopulationDots RDots COVDots TargetEstimateDots mtpopulationDots RTarget COVTarget TargetEstimateTarget
    
    eyeVelocity(k,j,trial,:) = gT(k) * squeeze(estimateSPDTarget(j,trial,:)) + (1 - gT(k)) * squeeze(estimateSPDDots(k,j,trial,:));
    eyeDirection(k,j,trial,:) = gT(k) * squeeze(estimateDIRTarget(j,trial,:)) + (1 - gT(k)) * squeeze(estimateDIRDots(k,j,trial,:));
    
    mdl = LinearModel.fit(squeeze(eyeVelocity(k,j,trial,:)),Target.MotionSpeed');
    bias(k,j,trial) = mdl.Coefficients.Estimate(1);
    sensitivity(k,j,trial) = mdl.Coefficients.Estimate(2);
    R(k,j,trial) = mdl.Rsquared.Adjusted;

end

end
end

%% Figures

figure;plot(npsDots.*npdDots,median(sensitivity,3),'-o','LineWidth',2);
xlabel('number of neurons');ylabel('Speed Sensitivity');legend(num2str(0),num2str(0.09),num2str(0.18),num2str(0.27),num2str(0.36),num2str(0.45))

figure;plot(npsDots.*npdDots,median(bias,3),'-o','LineWidth',2);
xlabel('number of neurons');ylabel('Speed Bias');legend(num2str(0),num2str(0.09),num2str(0.18),num2str(0.27),num2str(0.36),num2str(0.45))

figure;plot(npsDots.*npdDots,median(R,3),'-o','LineWidth',2);
xlabel('number of neurons');ylabel('Speed R^2');legend(num2str(0),num2str(0.09),num2str(0.18),num2str(0.27),num2str(0.36),num2str(0.45))

