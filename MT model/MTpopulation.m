function [mtpopulation, R, COV] = MTpopulation(Target,nps,npd,rmax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Population Parameters
global Nps Npd MinSpeed MaxSpeed MinDirection MaxDirection NumNeurons;
Nps = nps;
Npd = npd;
MinSpeed = 0.5;
MaxSpeed = 256;
MinDirection = 0;
MaxDirection = 360;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% paramDist = LoadMTdata();
CurrentFolder = pwd;
cd('D:\Data\Ephys\Dave Data\Data');
data = load('paramDist.mat');
cd(CurrentFolder);
paramDist = data.paramDist;

% randomly choosing the preferred speed and direction
% PSs = rand(1,Nps)*(log2(MaxSpeed) - log2(MinSpeed)) + log2(MinSpeed);
% PDs = rand(1,Npd)*(MaxDirection - MinDirection) + MinDirection;

% parcellating the range of speeds and directions for preferred speed and
% direction
PSs = log2(MinSpeed):(log2(MaxSpeed) - log2(MinSpeed))/Nps:(log2(MaxSpeed) - (log2(MaxSpeed) - log2(MinSpeed))/Nps);
PDs = MinDirection:(MaxDirection - MinDirection)/Npd:(MaxDirection - (MaxDirection - MinDirection)/Npd);

[PSmesh, PDmesh] = meshgrid(PSs,PDs);
NumNeurons = size(PSs,2)*size(PDs,2);
tic;
% simulate population firing rate and variance
mtpopulation = cell(1,NumNeurons);
for neuroncount = 1:NumNeurons
    if neuroncount == NumNeurons
        fprintf([num2str(neuroncount), ' neurons simulated \n'])
    end
    S.PS      =       PSmesh(neuroncount);      % Preferred Speed
    S.PD      =       PDmesh(neuroncount);      % Preferred Direction
    S.STW     =       (1.45);                   % Width of Speed Tuning Curve
    S.DTW     =       98*pi/180;                % Width of Direction Tuning Curve
    
    SampledMTneuron = Sample(paramDist,1,2);
    
    S.PR      =       SampledMTneuron(2);       % Preferred Size for preferred direction
    S.nPR     =       SampledMTneuron(7);       % Preferred Size for null direction
    S.RTW     =       SampledMTneuron(4);       % Width of Size Tuning Curve for preferred direction
    S.nRTW    =       SampledMTneuron(9);       % Width of Size Tuning Curve for null direction
    S.SI      =       SampledMTneuron(11);      % Suppression Index
    
    % The following lines should be uncommented for no surround suppression
%     S.PR      =       [];       % Preferred Size for preferred direction
%     S.nPR     =       [];       % Preferred Size for null direction
%     S.RTW     =       [];       % Width of Size Tuning Curve for preferred direction
%     S.nRTW    =       [];       % Width of Size Tuning Curve for null direction
%     S.SI      =       [];
    S.G       =       100;                     % Gain
    S.B0      =       2;                       % Base line activity
    
    MT = MTneuron(S);
    MT.SimulateMeanFiringRate(Target);
    MT.SimulateVariance;
    
    mtpopulation{neuroncount} = MT;
    clear MT S
end
fprintf('------------------------------------------ \n')
COV = ConstructCovariance(mtpopulation,rmax);

% R = nan(Nps*Npd,NumTrials);

R = MCsimulate(mtpopulation,COV)';
% R = MCsimulate(mtpopulation,COV,NumTrials)';


toc;
end

function COV = ConstructCovariance(mtpopulation,rmax)

global Nps Npd MinSpeed MaxSpeed NumNeurons;

fprintf(['Calculating Covariance Matrix ... '])
PD = cellfun(@(x)(x.PreferredDirection),mtpopulation)';
PS = cellfun(@(x)(x.PreferredSpeed),mtpopulation)';
rOnDiag = cellfun(@(x)(x.Variance),mtpopulation).^0.5;

% uncomment this for considering surround suppression
SI = cellfun(@(x)(x.SuppressionIndex),mtpopulation)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% correlations parameters
% rmax    =      0.3;
% if no surround suppression
% td      =      .45; 
% ts      =      .45;

% SI-dependent correlation matrix (proportional with SIs)
% td = nan(392,392);td(SI<0.5,SI<0.5) = 0.1;td(SI>0.5,SI>0.5) = 2;td(SI>0.5,SI<0.5) = 0.5;td(SI<0.5,SI>0.5) = 0.5;
% ts = nan(392,392);ts(SI<0.5,SI<0.5) = 0.1;ts(SI>0.5,SI>0.5) = 2;ts(SI>0.5,SI<0.5) = 0.5;ts(SI<0.5,SI>0.5) = 0.5;
[SI1, SI2] = meshgrid(SI,SI);

% SI-dependent correlation matrix (proportional with SIs)
td = 1*(SI1.*SI2 + .05); 
ts = 1*(SI1.*SI2 + .05);

% SI-independent correlation matrix
% td = 1*(.4*ones(size(SI1)).*ones(size(SI2)) + .05);
% ts = 1*(.4*ones(size(SI1)).*ones(size(SI2)) + .05);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


PDrep = repmat(PD,1,length(PD));
PSrep = repmat(PS,1,length(PS));

rOffDiag = rmax .* exp((-(PDrep.^2)-(PDrep'.^2)+2.*PD*PD')./((180 * td).^2)) .* ...
                exp((-(PSrep.^2)-(PSrep'.^2)+2.*PS*PS')./(((log2(MaxSpeed) - log2(MinSpeed)) * ts).^2)) ;

QonDiag = CovOnDiagonal(rOnDiag,NumNeurons);
QoffDiag = CovOffDiagoanl(rOffDiag,NumNeurons);
Q = QoffDiag - diag(diag(QoffDiag)) + diag(QonDiag);
COV = ((Q) * (Q)');

end

function R = MCsimulate(mtpopulation,COV)
global NumNeurons

for ncount = 1:(NumNeurons)
    MT = mtpopulation{ncount};
    mu(:,ncount) = MT.FiringRate;
    clear MT
end

try
R = mvnrnd(mu,COV);
catch
    COV
    pause
end
% R = mvnrnd(mu,COV,NumTrials);

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