function [mtpopulation, R, COV] = MTpopulation(Target,nps,npd,rmax,tauD,tauS)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Population Parameters
global Nps Npd MinSpeed MaxSpeed MinDirection MaxDirection NumNeurons;
Nps = nps;
Npd = npd;
MinSpeed = 0.5;
MaxSpeed = 256;
MinDirection = 0;
MaxDirection = 359;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% paramDist = LoadMTdata();
CurrentFolder = pwd;
cd('D:\Project Codes\Behavioral-Normal-Subject\MT model\DaveData');
data = load('paramDist.mat');
cd(CurrentFolder);
paramDist = data.paramDist;

% randomly choosing the preferred speed and direction
% PSs = rand(1,Nps)*(log2(MaxSpeed) - log2(MinSpeed)) + log2(MinSpeed);
% PDs = rand(1,Npd)*(MaxDirection - MinDirection) + MinDirection;

% parcellating the range of speeds and directions for preferred speed and
% direction
% PSs = log2(MinSpeed):(log2(MaxSpeed) - log2(MinSpeed))/Nps:(log2(MaxSpeed) - (log2(MaxSpeed) - log2(MinSpeed))/Nps);
% PDs = MinDirection:(MaxDirection - MinDirection)/Npd:(MaxDirection - (MaxDirection - MinDirection)/Npd);

% [PSmesh, PDmesh] = meshgrid(PSs,PDs);
PSmesh = rand(1,Nps * Npd)*(log2(MaxSpeed) - log2(MinSpeed)) + log2(MinSpeed);
PDmesh = rand(1,Npd * Nps)*(MaxDirection - MinDirection) + MinDirection;
% PSmesh = randi((log2(MaxSpeed) - log2(MinSpeed))+1,1,Nps * Npd) - 2;
% PDmesh = randi(((MaxDirection - MinDirection)+1)/8,1,Npd * Nps)*8 - 1 + MinDirection;
% NumNeurons = size(PSs,2)*size(PDs,2);
NumNeurons = length(PSmesh);
RFLocations = setRFLocations(Target);
tic;
% simulate population firing rate and variance
mtpopulation = cell(1,NumNeurons);
for neuroncount = 1:NumNeurons
    if neuroncount == NumNeurons
        fprintf([num2str(neuroncount), ' neurons simulated \n'])
    end
    S.PS      =       PSmesh(neuroncount);      % Preferred Speed
    S.PD      =       PDmesh(neuroncount);      % Preferred Direction
    S.STW     =       (1.5);                    % Width of Speed Tuning Curve 1.45
    S.DTW     =       40;                       % Width of Direction Tuning Curve 98*pi/180
    
    SampledMTneuron = Sample(paramDist,1,2);
    
    S.EA      =       SampledMTneuron(1);       % Excitation amplitude for preferred direction;
    S.PR      =       SampledMTneuron(2);       % Preferred Size for preferred direction
    S.IA      =       SampledMTneuron(3);       % Inhibition amplitude for preferred direction
    S.nEA     =       SampledMTneuron(6);       % Excitation amplitude for null direction;
    S.nPR     =       SampledMTneuron(7);       % Preferred Size for null direction
    S.nIA     =       SampledMTneuron(8);       % Inhibition amplitude for null direction
    S.RTW     =       SampledMTneuron(4);       % Width of Size Tuning Curve for preferred direction
    S.nRTW    =       SampledMTneuron(9);       % Width of Size Tuning Curve for null direction
    S.B       =       SampledMTneuron(5);       % Baseline firing rate for preferred direction;
    S.nB      =       SampledMTneuron(10);      % Baseline firing rate for null direction;
    S.SI      =       SampledMTneuron(11);      % Suppression Index
    
    % The following lines should be uncommented for no surround suppression
%     S.EA      =       [];
%     S.PR      =       [];       % Preferred Size for preferred direction
%     S.IA      =       [];
%     S.nEA     =       [];
%     S.nPR     =       [];       % Preferred Size for null direction
%     S.nIA     =       [];
%     S.RTW     =       [];       % Width of Size Tuning Curve for preferred direction
%     S.nRTW    =       [];       % Width of Size Tuning Curve for null direction
%     S.B       =       [];
%     S.nB      =       [];
%     S.SI      =       [];
    S.G       =       1;%8;%.1;%100;4             % Gain
    S.B0      =       0;%2;                       % Base line activity
    S.rfloc = RFLocations(neuroncount);            % Location of the receptive field (degree)
    
    MT = MTneuron(S);
    MT.SimulateMeanFiringRate(Target);
    MT.SimulateVariance;
    
    mtpopulation{neuroncount} = MT;
    clear MT S
end
fprintf('------------------------------------------ \n')
COV = ConstructCovariance(mtpopulation,rmax,tauD,tauS);

% R = nan(Nps*Npd,NumTrials);

R = MCsimulate(mtpopulation,COV)';
% R = MCsimulate(mtpopulation,COV,NumTrials)';


toc;
end

function RFLocations = setRFLocations(Target)
global NumNeurons
Sizes = 0:0.1:(Target.Size(1) - 0.1);
% popSize=zeros(1,5);
% errSize=zeros(1,5); % reviewer #2 asked us to calculate the error of our estimates

fun=@(x1) (6*x1.^-0.9); % Erickson (mm/deg)
for l=1:length(Sizes) % calculate the # of  neurons for each size
    %     popSize(i)=integral(fun,0.5,sizes(i))+1; % integral from 0.5 to edge, between 0 to 0.5 is roughly 1 mm (2mm/deg)
    totPopSize(l)=integral(fun,12 - Sizes(l),12 + Sizes(l)); % integral from 1 to edge, between 0 to 1 is roughly (9mm/deg) erickson et al. EBR
    %     popSize(i)=integral(fun,.5,sizes(i)) + 9;
    
    %     errSize(i)=integral(fun,0.5,sizes(i)*0.16+0.51)+1; % integral of the error
end

totPopSize=round(20*totPopSize); % new constant
for l = 2:length(Sizes)
    popSizeIncr(l-1) = totPopSize(l) - totPopSize(l-1);
end
popSizeIncr = [totPopSize(1),popSizeIncr];

RFLocations = [];
for l = 1:length(popSizeIncr)
    RFLocations = [RFLocations,Sizes(l)*ones(1,popSizeIncr(l))];
end

numExtraNeurons = (NumNeurons - length(RFLocations));
RFLocations = [RFLocations,Target.Size(1)*rand(1,numExtraNeurons)];

end

function COV = ConstructCovariance(mtpopulation,rmax,tauD,tauS)

global Nps Npd MinSpeed MaxSpeed NumNeurons;

fprintf(['Calculating Covariance Matrix ... '])
PD = cellfun(@(x)(x.PreferredDirection),mtpopulation)';
PS = cellfun(@(x)(x.PreferredSpeed),mtpopulation)';
% rOnDiag = cellfun(@(x)(x.Variance),mtpopulation).^0.5;

% uncomment this for considering surround suppression
SI = cellfun(@(x)(x.SuppressionIndex),mtpopulation)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% correlations parameters
% rmax    =      0.3;
% if no surround suppression
% td      =      tauD; 
% ts      =      tauS;

% SI-dependent correlation matrix (proportional with SIs)
% td = nan(392,392);td(SI<0.5,SI<0.5) = 0.1;td(SI>0.5,SI>0.5) = 2;td(SI>0.5,SI<0.5) = 0.5;td(SI<0.5,SI>0.5) = 0.5;
% ts = nan(392,392);ts(SI<0.5,SI<0.5) = 0.1;ts(SI>0.5,SI>0.5) = 2;ts(SI>0.5,SI<0.5) = 0.5;ts(SI<0.5,SI>0.5) = 0.5;
[SI1, SI2] = meshgrid(SI,SI);

% SI-dependent correlation matrix (proportional with SIs)
% td = tauD - (tauD-0.001)*(SI1 + SI2)./max(max(SI1 + SI2)); 
% ts = tauS - (tauS-0.001)*(SI1 + SI2)./max(max(SI1 + SI2));
% td(td < 0) = 1e-3;
% ts(ts < 0) = 1e-3;

% coeff_ss proportional to SIs
coeff_ss = (SI1 + SI2)./max(max(SI1 + SI2)); 

% coeff_ss soft threshold
% coeff_ss(SI1>median(SI) & SI2>median(SI)) = .4*randn(size(coeff_ss(SI1>median(SI) & SI2>median(SI))))+.6;
% coeff_ss(xor(SI1<median(SI),SI2<median(SI))) = .4*randn(size(coeff_ss(xor(SI1<median(SI),SI2<median(SI))))) + .3;
% coeff_ss(SI1<=median(SI) & SI2<=median(SI)) = .4*rand(size(coeff_ss(SI1<=median(SI) & SI2<=median(SI))));
% coeff_ss(SI1>median(SI) & SI2>median(SI)) = .2*randn(size(coeff_ss(SI1>median(SI) & SI2>median(SI))))+.8;
% coeff_ss(xor(SI1<median(SI),SI2<median(SI))) = .5*randn(size(coeff_ss(xor(SI1<median(SI),SI2<median(SI))))) + .3;
% coeff_ss(SI1<=median(SI) & SI2<=median(SI)) = .4*rand(size(coeff_ss(SI1<=median(SI) & SI2<=median(SI))));

% coeff_ss hard threshold
% coeff_ss(SI1>median(SI) & SI2>median(SI)) = 1;
% coeff_ss(xor(SI1<median(SI),SI2<median(SI))) = .5;
% coeff_ss(SI1<=median(SI) & SI2<=median(SI)) = 0;

% coeff_ss zero
coeff_ss = zeros(size(coeff_ss));
% coeff_ss = 1.3 * (max(maxx(SI1 + SI2)) - SI1 - SI2 - 2);

% SI-independent correlation matrix
td1 = 1*(tauD*ones(size(SI1)).*ones(size(SI2))); 
ts1 = 1*(tauS*ones(size(SI1)).*ones(size(SI2))); 

% these constants are being used for modeling no relationship between noise
% correlation and signal correlation
% td2 = 1*(2*ones(size(SI1)).*ones(size(SI2))); 
% ts2 = 1*(2*ones(size(SI1)).*ones(size(SI2))); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


PDrep = repmat(PD,1,length(PD));
PSrep = repmat(PS,1,length(PS));

PDdiff = abs(AngDiff(PDrep,PDrep'));
PSdiff = abs(PSrep - PSrep');
randomPDdiff = (max(max(PDdiff)) - min(min(PDdiff)))*rand(size(PDdiff));
randomPSdiff = (max(max(PSdiff)) - min(min(PSdiff)))*rand(size(PSdiff));
% rOffDiag = rmax .* exp((-(PDrep.^2)-(PDrep'.^2)+2.*PD*PD')./((360 * td).^2)) .* ...
%                 exp((-(PSrep.^2)-(PSrep'.^2)+2.*PS*PS')./(((log2(MaxSpeed) - log2(MinSpeed)) * ts).^2)) ;

% rOffDiag = rmax .* exp(-(PDdiff.^2)./((180 * td).^2)) .* ...
%                 exp(-(PSdiff.^2)./(((log2(MaxSpeed) - log2(MinSpeed)) * ts).^2)) ;
% rOffDiag1 = exp(-(PDdiff.^2)./((180 * td1).^2)) .* ...
%                 exp(-(PSdiff.^2)./(((log2(MaxSpeed) - log2(MinSpeed)) * ts1).^2)) ;

% Gaussian Kernel
rOffDiag1 = exp(-(((1-coeff_ss).*PDdiff + coeff_ss.*randomPDdiff).^2)./((180 * td1).^2)) .* ...
                exp(-(((1-coeff_ss).*PSdiff + coeff_ss.*randomPSdiff).^2)./(((log2(MaxSpeed) - log2(MinSpeed)) * ts1).^2)) ;

            
% Logistic Kernel
% rOffDiag1_1 = 1./((exp((((1-coeff_ss).*PDdiff + coeff_ss.*randomPDdiff).^2)./((180 * td1).^2))) + 2 + ...
%     (exp(-(((1-coeff_ss).*PDdiff + coeff_ss.*randomPDdiff).^2)./((180 * td1).^2))));
% rOffDiag2_1 = 1./(( exp((((1-coeff_ss).*PSdiff + coeff_ss.*randomPSdiff).^2)./(((log2(MaxSpeed) - log2(MinSpeed)) * ts1).^2))) + 2 + ...
%     ( exp(-(((1-coeff_ss).*PSdiff + coeff_ss.*randomPSdiff).^2)./(((log2(MaxSpeed) - log2(MinSpeed)) * ts1).^2))));
% rOffDiag1 = rOffDiag1_1 .* rOffDiag2_1;




% rOffDiag2 = exp(-(randomPDdiff.^2)./((180 * td1).^2)) .* ...
%                 exp(-(randomPSdiff.^2)./(((log2(MaxSpeed) - log2(MinSpeed)) * ts1).^2)) ;
% rOffDiag1 = rOffDiag2;
% rOffDiag2 = std2(rOffDiag1)*rand(size(rOffDiag1));
% rOffDiag = rmax * (coeff_ss .* rOffDiag2 + (1 - coeff_ss) .*  rOffDiag1);
% rOffDiag = (coeff_ss) .* 1 .* rOffDiag2 + (1 - coeff_ss) .* rmax .* rOffDiag1;
% rOffDiag = (exp(-15.5*coeff_ss)) .* rmax .* rOffDiag1; % exponential relationship
% rOffDiag = (-2./(1+exp(-15*coeff_ss)) + 2) .* rmax .* rOffDiag1; % sigmoid relationship
% rOffDiag = (-erf(b*coeff_ss)+1) .* rmax .* rOffDiag1 + r0; % erf relationship

% rOffDiag = nan(size(rOffDiag2));
% rOffDiag(coeff_ss>quantile(reshape(coeff_ss,size(coeff_ss,1)*size(coeff_ss,2),1),.8)) = rOffDiag2(coeff_ss>quantile(reshape(coeff_ss,size(coeff_ss,1)*size(coeff_ss,2),1),.8)) ;
% rOffDiag(coeff_ss<=quantile(reshape(coeff_ss,size(coeff_ss,1)*size(coeff_ss,2),1),.8)) = rOffDiag1(coeff_ss<=quantile(reshape(coeff_ss,size(coeff_ss,1)*size(coeff_ss,2),1),.8));

% rOffDiag = (exp(-12*(1-coeff_ss))) .* rmax .* rOffDiag2 + (exp(-12*coeff_ss)) .* rmax .* rOffDiag1;
% rOffDiag = (exp(-2*(coeff_ss))) .* rmax .* rOffDiag1 + (1-exp(-2*coeff_ss)) .* rmax .* rOffDiag2;

rmax = rmax./ (sum(sum(rOffDiag1 - diag(diag(rOffDiag1))))./(size(rOffDiag1,1)*(size(rOffDiag1,1)-1)));
rOffDiag = rmax .*rOffDiag1;
% rOffDiag = rOffDiag1 - (sum(sum(rOffDiag1 - diag(diag(rOffDiag1))))./(size(rOffDiag1,1)*(size(rOffDiag1,1)-1))) + rmax;



rOnDiag = sum(sum(rOffDiag - diag(diag(rOffDiag))))./(size(rOffDiag,1)*(size(rOffDiag,1)-1)) * ones(1,NumNeurons);
QonDiag = CovOnDiagonal(rOnDiag,NumNeurons);
QoffDiag = CovOffDiagoanl(rOffDiag,NumNeurons);
Q = QoffDiag - diag(diag(QoffDiag)) + diag(QonDiag);
COV = ((Q) * (Q)');
% COV = COV - diag(diag(COV)) + diag(rOnDiag);

end

function R = MCsimulate(mtpopulation,COV)
global NumNeurons

for ncount = 1:(NumNeurons)
    MT = mtpopulation{ncount};
    mu(:,ncount) = MT.FiringRate;
    sigma(:,ncount) = MT.Variance;
    
    clear MT
end

try
R = mvnrnd(zeros(size(mu)),COV);
catch
    COV
    pause
end
% R = mvnrnd(mu,COV,NumTrials);
R = R.*sqrt(sigma) + mu;
end

function v = CovOffDiagoanl(r,N)
v = (1/N^0.5) .* abs((2/N) + r - (2*r)/N - (2/N).*((1 - r).*(1 - r + r*N)).^0.5).^0.5;
end

function u = CovOnDiagonal(r,N)
u = (1./(r*N^0.5)) .* (2/N + r - 2.*r/N - (2/N).*((1 - r).*(1 - r + r*N)).^0.5).^0.5 .* ...
    (1 + ((1 - r).*(1 - r + r*N)).^0.5);
end

function s = Sample(v,n,dim) % this function randomly draws n samples from the matrix v along dimension dim
if dim==1
    L = size(v,1);
    s = nan(n,size(v,2));
else
    L = size(v,2);
    s = nan(n,size(v,1));
end

for samplenumber = 1:n
    i = randi(L);
    if dim == 1
        s(samplenumber,:) = v(i,:);
    else
        s(samplenumber,:) = v(:,i);
    end
end
end

function paramDist = LoadMTdata() % this function runs the codes for reading Dave's MT data - adopted from Dave's getData script
% uiopen;
load('D:\Data\Ephys\Dave Data\Data\dataRalfO.mat');
CurrentFolder = pwd;
DaveCodesFolder = 'D:\Data\Ephys\Dave Data\Data';
cd(DaveCodesFolder);

% form the size tuning curves
tuningCurve=squeeze([tuningcurve(1,:,:) tuningcurve(3,:,:)]); % size tuning curves, f, is firing rate for preferred and null directions
tuningCurve(1,psychParam(:,2)==0)=tuningcurve(1,2,psychParam(:,2)==0); % some with 0 as the first size,
tuningCurve(9,psychParam(:,2)==0)=tuningcurve(3,2,psychParam(:,2)==0); % these are 0:2:14

% fit parameters to the responses
disp('Fitting Size Tuning Curves ... ')
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

cd(CurrentFolder);

end
