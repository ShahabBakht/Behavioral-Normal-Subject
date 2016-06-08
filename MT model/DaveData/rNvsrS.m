function ri=rNvsrS(a,inputData)

SI=inputData(1,:)';
tuningCurve=inputData(2:end,:);

xvec = SI'; yvec = SI'; [x y] = meshgrid(xvec,yvec);
siMat=x+y; % SI matrix is sum of the 2 SIs
siMat=max(max(siMat))-siMat; % high SI is surround suppressed, flip the scale
siMat(logical(eye(size(siMat))))=1; % set diagonal to 1

surrSuppx=zeros(1,length(SI)); surrSuppy=surrSuppx;
surrSuppx(SI>median(SI))=1; surrSuppy(SI>median(SI))=1; % surround suppressed or not
[x y] = meshgrid(surrSuppx,surrSuppy);
surrSuppMat=x+y; surrSuppMat(logical(eye(size(surrSuppMat))))=3; % categorize to SS-SS, SS-NS, and NS-NS

signalCorr=corrcoef(tuningCurve); signalCorr(logical(eye(size(signalCorr))))=0;
noiseCorr=a(1)*(siMat-a(2))*0.0753.*signalCorr+0.711; % SI term
noiseCorr(logical(eye(size(noiseCorr))))=1; % set diagonal to 1

p_nsSim=polyfit(signalCorr(surrSuppMat==0),noiseCorr(surrSuppMat==0),1);
p_ssSim=polyfit(signalCorr(surrSuppMat==2),noiseCorr(surrSuppMat==2),1);
p_ssnsSim=polyfit(signalCorr(surrSuppMat==1),noiseCorr(surrSuppMat==1),1);
ri=[p_nsSim(1) p_ssnsSim(1) p_ssSim(1)];
