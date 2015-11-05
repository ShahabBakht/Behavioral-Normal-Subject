function V = OpenLoopVelocity(TestName,PatientName)

InitiationTime = 1500;
SampleRate = 0.001;
S = SaccadeDetection(TestName,PatientName);
SaccadeTimes = squeeze(S(:,:,2));
EndTime = min(min(SaccadeTimes));
I = Eye(TestName,PatientName);
I.LoadEyeFlag = true;
I.LoadPreProcessedEye;
X = I.PreProcessedEye.EyePreProcessed.Xtrunc;
NumConditions = size(X,1);
NumTrials = size(X,2);
for c = 1:NumConditions
    for tr = 1:NumTrials
        if SaccadeTimes(c,tr) ~= 0
            xnow = squeeze(X(c,tr,InitiationTime:SaccadeTimes(c,tr)));
            vnow = MeasureVelocity(xnow,SampleRate);
            V(c,tr,:) = vnow(InitiationTime:EndTime);
        end
    end
end
        


end

function v = MeasureVelocity(x,sr)
[b,a] = butter(6,20*2*sr);
xfit = filtfilt(b,a,x);
v = gradient(xfit,sr);

end