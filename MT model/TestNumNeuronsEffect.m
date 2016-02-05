%% Number of neurons as a function of stimulus size

sizes=sqrt(5)+(1:2:19); % the avg eccentricity of the stimuli is sqrt(5) deg
popSize=zeros(1,8);
errSize=zeros(1,8); % reviewer #2 asked us to calculate the error of our estimates
% fun=@(x1) (1.14*x1.^-0.76); % Albright & Desimone 87 Exp Brain Res (mm/deg)
fun=@(x1) (6*x1.^-0.9); % Erickson (mm/deg)
for i=1:length(sizes) % calculate the # of neurons for each size
%     popSize(i)=integral(fun,0.5,sizes(i))+1; % integral from 0.5 to edge, between 0 to 0.5 is roughly 1 mm (2mm/deg)
    popSize(i)=integral(fun,0.5,sizes(i))+9; % integral from 1 to edge, between 0 to 1 is roughly (9mm/deg) erickson et al. EBR
    errSize(i)=integral(fun,0.5,sizes(i)*0.16+0.51)+1; % integral of the error
end
% popSize=round(1*popSize.^2); % multiple by a constant of 10, but the constant did not affect the trend
popSize=round(20.5*popSize); % new constant
errSize=round(10*errSize); % multiple the error a constant of 10

figure % plot for the manuscript
plot(sizes-sqrt(5),popSize,'k')
set(gca,'XTick',[0 5 10 15]);
set(gca,'YTick',[400 600 800 1000]);
ylim([400 1100])
set(gca,'TickDir','out')
xlabel(['Size (' char(176) ')']); ylabel('Estimated pool size')

box off

%% Single Speed for all the trials
% nps = 12;
% npd = 24;
for trial = 1:10
Target.MotionSpeed        =   16 * ones(500,1);
Target.MotionDirection    =   0 * ones(500,1);

j = 0;
% DIRstd = nan(length([1,1.5,2,3,4]),length([0,0.09,0.18,0.27,0.36,0.45]));
% DIRm = nan(length([1,1.5,2,3,4]),length([0,0.09,0.18,0.27,0.36,0.45]));
% SPDstd = nan(length([1,1.5,2,3,4]),length([0,0.09,0.18,0.27,0.36,0.45]));
% SPDm = nan(length([1,1.5,2,3,4]),length([0,0.09,0.18,0.27,0.36,0.45]));
for rmax = [0,0.09,0.18,0.27,0.36,0.45]
    j = j + 1;
    k = 0;
for i = max(popSize)./popSize%[1,1.5,2,3,4]%log2(2:(62/5):64-62/5)
    k = k + 1;
    Target.Size = sizes(k) * ones(500,1);
    nps(k) = floor(46/(2*i));
    npd(k) = floor(90/(2*i));
    
    [mtpopulation, R, COV] = MTpopulation(Target,500,nps(k),npd(k),rmax);
    [TargetEstimate] = DecodeMTpopulation(mtpopulation, real(R));
    DIRstd(k,j,trial) = std(TargetEstimate.DIRest);
    SPDstd(k,j,trial) = std(TargetEstimate.SPDest);
%     DIRm(k,j) = median(TargetEstimate.DIRest);
%     SPDm(k,j) = median(TargetEstimate.SPDest);
    
    clear mtpopulation R COV TargetEstimate
    
end

end
end
%% Plots - Speed and Direction variance as a function of #neurons

figure;plot(nps.*npd,median(SPDstd,3),'-o','LineWidth',2);
ax = gca;
set(ax,'XTick',[200,1200,2200,3200,4200]);
xlabel('number of neurons');ylabel('speed variance')
figure;plot(nps.*npd,median(DIRstd,3),'-o','LineWidth',2);
ax = gca;
set(ax,'XTick',[200,1200,2200,3200,4200]);
xlabel('number of neurons');ylabel('direction variance')


%% Random Target Speed for any trial
% Target.MotionSpeed        =   10 * rand(500,1) + 10;
% Target.MotionDirection    =   30 * ones(500,1) + 10;
% 
% j = 0;m = 0;
% CorrSPD = nan(length([1,1.5,2,3,4]),length([0,0.09,0.18,0.27,0.36,0.45]));
% CorrDIR= nan(length([1,1.5,2,3,4]),length([0,0.09,0.18,0.27,0.36,0.45]));
% for rmax = [0,0.09,0.18,0.27,0.36,0.45]
%     j = j + 1;
%     k = 0;
% for i = [1,1.5,2,3,4]%log2(2:(62/5):64-62/5)
%     k = k + 1;
%     m = m + 1;
%     nps(k) = floor(46/i);
%     npd(k) = floor(90/i);
%     [mtpopulation, R, COV] = MTpopulation(Target,1,nps(k),npd(k),rmax);
%     [TargetEstimate] = DecodeMTpopulation(mtpopulation, real(R));
%     Cs = corrcoef(TargetEstimate.SPDest,Target.MotionSpeed);
%     Cd = corrcoef(TargetEstimate.DIRest,Target.MotionDirection);
%     CorrSPD(k,j) = Cs(1,2);
%     CorrDIR(k,j) = Cd(1,2);
%     clear mtpopulation R COV TargetEstimate
%     
% end
% 
% end

%%

% figure;plot(nps.*npd,CorrSPD(:,:),'-o','LineWidth',2);
% ax = gca;
% set(ax,'XTick',[200,1200,2200,3200,4200]);
% xlabel('number of neurons');ylabel('speed correlation');legend(num2str(0),num2str(0.09),num2str(0.18),num2str(0.27),num2str(0.36),num2str(0.45))
% figure;plot(nps.*npd,CorrDIR(:,:),'-o','LineWidth',2);
% ax = gca;
% set(ax,'XTick',[200,1200,2200,3200,4200]);
% xlabel('number of neurons');ylabel('direction correlation');legend(num2str(0),num2str(0.09),num2str(0.18),num2str(0.27),num2str(0.36),num2str(0.45))