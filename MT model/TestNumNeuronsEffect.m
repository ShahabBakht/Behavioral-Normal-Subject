%% Single Speed for all the trials

Target.MotionSpeed        =   16 * ones(500,1);
Target.MotionDirection    =   135 * ones(500,1);

j = 0;
DIRstd = nan(length([1,1.5,2,3,4]),length([0,0.09,0.18,0.27,0.36,0.45]));
SPDstd = nan(length([1,1.5,2,3,4]),length([0,0.09,0.18,0.27,0.36,0.45]));
for rmax = [0,0.09,0.18,0.27,0.36,0.45]
    j = j + 1;
    k = 0;
for i = [1,1.5,2,3,4]%log2(2:(62/5):64-62/5)
    k = k + 1;
    
    nps(k) = floor(46/i);
    npd(k) = floor(90/i);
    [mtpopulation, R, COV] = MTpopulation(Target,500,nps(k),npd(k),rmax);
    [TargetEstimate] = DecodeMTpopulation(mtpopulation, real(R));
    DIRstd(k,j) = std(TargetEstimate.DIRest);
    SPDstd(k,j) = std(TargetEstimate.SPDest);
    
    clear mtpopulation R COV TargetEstimate
    
end

end

%% Plots - Speed and Direction variance as a function of #neurons

figure;plot(nps.*npd,SPDstd(:,:),'-o','LineWidth',2);
ax = gca;
set(ax,'XTick',[200,1200,2200,3200,4200]);
xlabel('number of neurons');ylabel('speed variance')
figure;plot(nps.*npd,DIRstd(:,:),'-o','LineWidth',2);
ax = gca;
set(ax,'XTick',[200,1200,2200,3200,4200]);
xlabel('number of neurons');ylabel('direction variance')


%% Random Target Speed for any trial
Target.MotionSpeed        =   10 * rand(500,1) + 10;
Target.MotionDirection    =   30 * rand(500,1) + 10;

j = 0;m = 0;
CorrSPD = nan(length([1,1.5,2,3,4]),length([0,0.09,0.18,0.27,0.36,0.45]));
CorrDIR= nan(length([1,1.5,2,3,4]),length([0,0.09,0.18,0.27,0.36,0.45]));
for rmax = [0,0.09,0.18,0.27,0.36,0.45]
    j = j + 1;
    k = 0;
for i = [1,1.5,2,3,4]%log2(2:(62/5):64-62/5)
    k = k + 1;
    m = m + 1;
    nps(k) = floor(46/i);
    npd(k) = floor(90/i);
    [mtpopulation, R, COV] = MTpopulation(Target,1,nps(k),npd(k),rmax);
    [TargetEstimate] = DecodeMTpopulation(mtpopulation, real(R));
    Cs = corrcoef(TargetEstimate.SPDest,Target.MotionSpeed);
    Cd = corrcoef(TargetEstimate.DIRest,Target.MotionDirection);
    CorrSPD(k,j) = Cs(1,2);
    CorrDIR(k,j) = Cd(1,2);
    clear mtpopulation R COV TargetEstimate
    
end

end

%%

figure;plot(nps.*npd,CorrSPD(:,:),'-o','LineWidth',2);
ax = gca;
set(ax,'XTick',[200,1200,2200,3200,4200]);
xlabel('number of neurons');ylabel('speed correlation');legend(num2str(0),num2str(0.09),num2str(0.18),num2str(0.27),num2str(0.36),num2str(0.45))
figure;plot(nps.*npd,CorrDIR(:,:),'-o','LineWidth',2);
ax = gca;
set(ax,'XTick',[200,1200,2200,3200,4200]);
xlabel('number of neurons');ylabel('direction correlation');legend(num2str(0),num2str(0.09),num2str(0.18),num2str(0.27),num2str(0.36),num2str(0.45))