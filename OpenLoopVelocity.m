function V = OpenLoopVelocity(TestName,PatientName,S)

InitiationTime = 1500;
SampleRate = 0.001;
% S = SaccadeDetection(TestName,PatientName);
SaccadeTimes = squeeze(S(:,:,2));
SaccadeTimes(SaccadeTimes < 140) = nan;
EndTime = min(SaccadeTimes,250*ones(size(SaccadeTimes)));
I = Eye(TestName,PatientName);
I.LoadEyeFlag = true;
I.LoadPreProcessedEye;
X = I.PreProcessedEye.EyePreProcessed.Xtrunc;
NumConditions = size(X,1);
NumTrials = size(X,2);
for c = 1:NumConditions
    for tr = 1:NumTrials
        if ~isnan(SaccadeTimes(c,tr))
            xnow = squeeze(X(c,tr,InitiationTime:floor(SaccadeTimes(c,tr))+1600));
%             plot(xnow);pause;
        else
            xnow = squeeze(X(c,tr,InitiationTime:2000));
            
        end
%         
%         vnow = MeasureVelocity(xnow,SampleRate);
%         plot(vnow);pause
%         V(c,tr,:) = vnow(100:floor(EndTime)+100);
        if EndTime(c,tr) < 240
            V(c,tr,:) = 0;
        else
            V(c,tr,:) = (nanmean(xnow(240:260)) - xnow(floor(EndTime(c,tr))+100))./((EndTime(c,tr)/1000)-(0.15));
        end
    end
end
        


end

function v = MeasureVelocity(x,sr)
[b,a] = butter(6,50*2*sr);
xfit = filtfilt(b,a,x);
v = gradient(xfit,sr);

end