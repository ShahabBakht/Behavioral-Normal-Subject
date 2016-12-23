%%  Simulate the model 
simSlope = 10;
integR = 2;
numSamples = 20;
SuppIdx_data = simulateSuppIdx(simSlope,integR,numSamples);

%% simulate Likelihood for grid of slopes
rcount = 0;
slopecount = 0;
r = [2,4,8,16];
slope = 0:2:20;
parpool(6);
parfor slopecount = 1:length(slope)
    for rcount = 1:length(r)
        
        FF(slopecount,rcount,:) = simulateSuppIdx(slope,2,5000);
    end
end
