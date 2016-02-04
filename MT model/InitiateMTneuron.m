%% MT neuron tuning curve parameters

S.PS      =       10;     % Preferred Speed 
S.PD      =       0;     % Preferred Direction
S.STW     =       2;      % Width of Speed Tuning Curve
S.DTW     =       2;      % Width of Direction Tuning Curve
S.G       =       1;      % Gain

%% Target Motion Parameters;

Target.MotionSpeed        =   10 * rand(500,1) + 10;
Target.MotionDirection    =   150 * ones(500,1);

%%

MT = MTneuron(S);
MT.SimulateMeanFiringRate(Target);