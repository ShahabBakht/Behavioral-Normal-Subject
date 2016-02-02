Target.MotionSpeed        =   16;
Target.MotionDirection    =   135;

j = 0;
DIRstd = nan(length([1,1.5,2,4,6]),length([0,0.09,0.18,0.27,0.36,0.45]));
SPDstd = nan(length([1,1.5,2,4,6]),length([0,0.09,0.18,0.27,0.36,0.45]));
for rmax = [0,0.09,0.18,0.27,0.36,0.45]
    j = j + 1;
    k = 0;
for i = [1,1.5,2,4,6]%log2(2:(62/5):64-62/5)
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

figure;h=plot(nps.*npd,SPDstd(:,:),'-o','LineWidth',2);
ax = gca;
set(ax,'XTick',[200,1200,2200,3200,4200]);
xlabel('number of neurons');ylabel('speed variance')
figure;plot(nps.*npd,DIRstd(:,:),'-o','LineWidth',2);
ax = gca;
set(ax,'XTick',[200,1200,2200,3200,4200]);
xlabel('number of neurons');ylabel('direction variance')
    