% %% Get the size tuning curves, f, for preferred and null directions
tuningCurve=squeeze([tuningcurve(1,:,:) tuningcurve(3,:,:)]); % size tuning curves, f, is firing rate for preferred and null directions
tuningCurve(1,psychParam(:,2)==0)=tuningcurve(1,2,psychParam(:,2)==0); % some with 0 as the first size,
tuningCurve(9,psychParam(:,2)==0)=tuningcurve(3,2,psychParam(:,2)==0); % these are 0:2:14
variance=squeeze([tuningcurve(2,:,:) tuningcurve(4,:,:)]); % get the variance
variance(1,psychParam(:,2)==0)=tuningcurve(2,2,psychParam(:,2)==0);% need to replace the size 0 responses with size 2 responses
variance(9,psychParam(:,2)==0)=tuningcurve(4,2,psychParam(:,2)==0);% fit the response function would solve the problem

%% fit parameters to the responses
paramDist=zeros(10,size(tuningCurve,2)); % the parameter distribution of the variables
for i=1:size(tuningCurve,2)
    if tuningCurve(1,i)==tuningCurve(2,i)
        [a,si,fitType]=mtSizeFit([0 2:2:14],[baselinePop(i) tuningCurve(2:8,i)'],1); % fit DoE function
        paramDist(1:5,i)=a;
        tuningCurve(1:8,i)=diffGauss(a,1:2:15);
        [a,si,fitType]=mtSizeFit([0 2:2:14],[baselinePop(i) tuningCurve(10:end,i)'],1); % fit DoE function
        paramDist(6:10,i)=a;
        tuningCurve(9:end,i)=diffGauss(a,1:2:15);
    elseif tuningCurve(1,i)~=tuningCurve(2,i)
        [a,si,fitType]=mtSizeFit([0 1:2:15],[baselinePop(i) tuningCurve(1:8,i)'],1); % fit DoE function
        paramDist(1:5,i)=a;
        tuningCurve(1:8,i)=diffGauss(a,1:2:15);
        [a,si,fitType]=mtSizeFit([0 1:2:15],[baselinePop(i) tuningCurve(9:end,i)'],1); % fit DoE function
        paramDist(6:10,i)=a;
        tuningCurve(9:end,i)=diffGauss(a,1:2:15);
    end
end
%% variance relative to mean
fanoFactor=variance./tuningCurve; % The ratio of variance to mean is around 1.2
mFanoFactor=mean(mean(fanoFactor(~isnan(fanoFactor)))); % high due to rate
%% add variance as a parameter (independent)
% paramDist=[paramDist(1:5,:); mean(variance./tuningCurve,1)/10; paramDist(6:end,:)];
%% add SI as a parameter (doesn't work)
% paramDist=[paramDist(1:5,:); SI'; paramDist(6:end,:)];
%% add eccentricity as a parameter (no need)
% paramDist=[paramDist(1:5,:); relEcc1'; paramDist(6:end,:)];
%% Plot parameters of the tuning curve
figure
paramLabel=['ExciAmpl';'ExciRadi';'InhiAmpl';'InhiRadi';'Baseline';'ExciAmpN';'ExciRadN';'InhiAmpN';'InhiRadN';'BaselinN'];
numP=5;
for i=1:numP
    subplot(numP,numP,i)
    hist(paramDist(i,:))% histogram
    xlim([0 max(paramDist(i,:))]);
    box off
    axis square
    for j=i+1:numP
        subplot(numP,numP,i*numP+j)
        plot(paramDist(j,:),paramDist(i,:),'.') % plot the parameters
        xlim([0 max(paramDist(j,:))]); ylim([0 max(paramDist(i,:))]);
        if j==i+1
            ylabel(paramLabel(i,:)); xlabel(paramLabel(j,:));
        end
        box off
        r=corrcoef(paramDist(j,:),paramDist(i,:));
        title(num2str(round(r(1,2)*100)/100))       
        axis square
    end
end
%% remove variance as a parameter
% paramDist=[paramDist(1:5,:); paramDist(7:end,:)];
%% Fit copula to data
paramDistT=zeros(size(paramDist)); % transform the data to the copula scale (unit square)
for i=1:size(paramDist,1)
    paramDistT(i,:)=ksdensity(paramDist(i,:),paramDist(i,:),'function','cdf');
end
rng default % For reproducibility
rho = copulafit('Gaussian',paramDistT');
u=copularnd('Gaussian',rho,164)';
for i=1:size(paramDist,1)
    u(i,:)=ksdensity(paramDist(i,:),u(i,:),'function','icdf');
end % transform the random sample back
u(u<0)=0; % set to 0

figure
numP=9;
for i=1:numP
    subplot(numP,numP,i)
    hist(u(i,:))% histogram
    xlim([0 max(u(i,:))]);
    box off
    axis square
    for j=i+1:numP
        subplot(numP,numP,i*numP+j)
        plot(u(j,:),u(i,:),'.') % plot the parameters
        xlim([0 max(u(j,:))]); ylim([0 max(u(i,:))]);
        if j==i+1
            ylabel(paramLabel(i,:)); xlabel(paramLabel(j,:));
        end
        box off
        r=corrcoef(u(j,:),u(i,:));
        title(num2str(round(r(1,2)*100)/100))       
        axis square
    end
end
%% find the offset and slope to build covariance matrix
inputData=[SI'; tuningCurve]; % the xdata is the SI and tuning curves
[a1,resnorm]=lsqcurvefit(@rNvsrS,[1.4 2.1],inputData,[p_ns(1) p_ssns(1) p_ss(1)],[1 1],[2 3]);
%% noise correlation (Build the covariance matrix)
xvec = SI'; yvec = SI'; [x y] = meshgrid(xvec,yvec);
siMat=x+y; % SI matrix is sum of the 2 SIs
siMat=max(max(siMat))-siMat; % high SI is surround suppressed, flip the scale
siMat(logical(eye(size(siMat))))=1; % set diagonal to 1

surrSuppx=zeros(1,length(SI)); surrSuppy=surrSuppx;
surrSuppx(SI>median(SI))=1; surrSuppy(SI>median(SI))=1; % surround suppressed or not
[x y] = meshgrid(surrSuppx,surrSuppy);
surrSuppMat=x+y; surrSuppMat(logical(eye(size(surrSuppMat))))=3; % categorize to SS-SS, SS-NS, and NS-NS

signalCorr=corrcoef(tuningCurve); signalCorr(logical(eye(size(signalCorr))))=0;
%% SI and eccentricity
plot(relEcc1,SI,'.')
p=polyfit(relEcc1,SI,1)
hold on
plot(0:15,p(1)*(0:15)+p(2))
%% Simulation
sizes=sqrt(5)+(1:2:15); % the avg eccentricity of the stimuli is sqrt(5) deg
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