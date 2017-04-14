
%% simulate Likelihood for grid of slopes
rcount = 0;
slopecount = 0;
r = [2,4,8,16];
slope = (0:2:20);
[rMesh,slopeMesh] = meshgrid(r,slope);
parpool(6);
parfor k = 1:length(r)*length(slope)
    thisslope = slopeMesh(k);
    thisr = rMesh(k);
    FF(k,:) = simulateSuppIdx(thisslope,thisr,5000);
end
% parfor slopecount = 1:length(slope)
%     for rcount = 1:length(r)
%         
%         FF(slopecount,rcount,:) = simulateSuppIdx(slope(slopecount),r(rcount),5000);
%     end
% end

