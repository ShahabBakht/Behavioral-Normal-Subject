function [mtpopulation, R, COV] = MTpopulation(Target,NumTrials,nps,npd,rmax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Population Parameters
global Nps Npd MinSpeed MaxSpeed MinDirection MaxDirection;
Nps = nps;
Npd = npd;
MinSpeed = 0.5;
MaxSpeed = 256;
MinDirection = 0;
MaxDirection = 360;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PSs = log2(MinSpeed):(log2(MaxSpeed) - log2(MinSpeed))/Nps:(log2(MaxSpeed) - (log2(MaxSpeed) - log2(MinSpeed))/Nps);
PDs = MinDirection:(MaxDirection - MinDirection)/Npd:(MaxDirection - (MaxDirection - MinDirection)/Npd);
[PSmesh, PDmesh] = meshgrid(PSs,PDs);

tic;
% simulate population firing rate and variance
mtpopulation = cell(1,Nps*Npd);
for neuroncount = 1:Nps*Npd
    fprintf([num2str(neuroncount), ' neurons simulated \n'])
    S.PS      =       PSmesh(neuroncount);     % Preferred Speed
    S.PD      =       PDmesh(neuroncount);     % Preferred Direction
    S.STW     =       (1.45);                    % Width of Speed Tuning Curve
    S.DTW     =       98*pi/180;               % Width of Direction Tuning Curve
    S.G       =       11.5;                    % Gain
    
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

global Nps Npd MinSpeed MaxSpeed;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% correlations parameters
% rmax    =      0.3;
td      =      1;
ts      =      .5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(['Calculating Covariance Matrix ... '])
PD = cellfun(@(x)(x.PreferredDirection),mtpopulation)';
PS = cellfun(@(x)(x.PreferredSpeed),mtpopulation)';
rOnDiag = cellfun(@(x)(x.Variance),mtpopulation).^0.5;

PDrep = repmat(PD,1,length(PD));
PSrep = repmat(PS,1,length(PS));

rOffDiag = rmax .* exp((-(PDrep.^2)-(PDrep'.^2)+2.*PD*PD')./((180 * td)^2)) .* ...
                exp((-(PSrep.^2)-(PSrep'.^2)+2.*PS*PS')/(((log2(MaxSpeed) - log2(MinSpeed)) * ts)^2)) ;

QonDiag = CovOnDiagonal(rOnDiag,Nps*Npd);
QoffDiag = CovOffDiagoanl(rOffDiag,Nps*Npd);
Q = QoffDiag - diag(diag(QoffDiag)) + diag(QonDiag);
COV = Q * Q';

end

function R = MCsimulate(mtpopulation,COV)
global Nps Npd

for ncount = 1:(Nps*Npd)
    MT = mtpopulation{ncount};
    mu(:,ncount) = MT.FiringRate;
    clear MT
end

R = mvnrnd(mu,COV);
% R = mvnrnd(mu,COV,NumTrials);

end

function v = CovOffDiagoanl(r,N)
v = (1/N^0.5) .* abs((2/N) + r - (2*r)/N - (2/N).*((1 - r).*(1 - r + r*N)).^0.5).^0.5;
end

function u = CovOnDiagonal(r,N)
u = (1./(r*N^0.5)) .* (2/N + r - 2.*r/N - (2/N).*((1 - r).*(1 - r + r*N)).^0.5).^0.5 .* ...
    (1 + ((1 - r).*(1 - r + r*N)).^0.5);
end
