classdef MTneuron < handle
    
    properties
        PreferredSpeed
        PreferredDirection
        PreferredSizePD
        PreferredSizeND
        SpeedTuningWidth
        DirectionTuningWidth
        SizeTuningWidthPD
        SizeTuningWidthND
        InhibitionAmplitudePD
        InhibitionAmplitudeND
        ExcitationAmplitudePD
        ExcitationAmplitudeND
        BaseLinePD
        BaseLineND
        SuppressionIndex
        Gain
        BaseLine
        FiringRate
        Variance
        RFLocation
    end
    
    methods
        function MT = MTneuron(S)
            MT.PreferredSpeed = S.PS;
            MT.PreferredDirection = S.PD;
            MT.PreferredSizePD = S.PR;
            MT.PreferredSizeND = S.nPR;
            MT.SpeedTuningWidth = S.STW;
            MT.DirectionTuningWidth = S.DTW;
            MT.SizeTuningWidthPD = S.RTW;
            MT.SizeTuningWidthND = S.nRTW;
            MT.InhibitionAmplitudePD = S.IA;
            MT.InhibitionAmplitudeND = S.nIA;
            MT.ExcitationAmplitudePD = S.EA;
            MT.ExcitationAmplitudeND = S.nEA;
            MT.BaseLinePD = S.B;
            MT.BaseLineND = S.nB;
            MT.SuppressionIndex = S.SI;
            MT.Gain = S.G;
            MT.BaseLine = S.B0;
            MT.RFLocation = S.rfloc;
        end
        
        function SimulateMeanFiringRate(MT,Target)
            % neurons tuning characteristics
            G = MT.Gain;
            B0 = MT.BaseLine;
            PD = MT.PreferredDirection;
            PS = MT.PreferredSpeed;
            PR = MT.PreferredSizePD;
            nPR = MT.PreferredSizeND;
            STW = MT.SpeedTuningWidth;
            DTW = MT.DirectionTuningWidth * pi/180;
            RTW = MT.SizeTuningWidthPD;
            nRTW = MT.SizeTuningWidthND;
            IA = MT.InhibitionAmplitudePD;
            nIA = MT.InhibitionAmplitudeND;
            EA = MT.ExcitationAmplitudePD;
            nEA = MT.ExcitationAmplitudeND;
            B = MT.BaseLinePD;
            nB = MT.BaseLineND;
            rfLocation = MT.RFLocation;
            
            % target motion parameters
            Theta = Target.MotionDirection;
            S = Target.MotionSpeed;
            R = Target.Size;
            
            
            
%             Theta = Theta * pi / 180;
%             PD = PD * pi / 180;
%             if (Theta - PD) > pi
%                 DDiff = -pi + ((Theta - PD) - pi);
%             elseif (Theta - PD) < -pi
%                 DDiff = pi - (-pi - (Theta - PD));
%             else
%                 DDiff = (Theta - PD);
%             end
            DDiff = abs(AngDiff(Theta,PD) * pi/180);
            
            DirectionTune = exp(-( DDiff.^2 )/( 2 * DTW^2 ));
            SpeedTune = exp(-( (log2(S) - (PS)).^2 )./( 2 * STW^2 ));
            possibleSizes = 0:0.1:10;
            if ~isempty(PR) % if size tuning is defined for the neurons
                if (DDiff <= pi/2) 
                    SizeTune = EA*erf((R-rfLocation)./PR) - IA*erf((R-rfLocation)./(PR + RTW)) + B;
                else
                    SizeTune = nEA*erf((R-rfLocation)./nPR) - nIA*erf((R-rfLocation)./(nPR + nRTW)) + nB; 
                end
            else
                SizeTune = 1;
            end
%             fit = EA*erf((possibleSizes)./PR) - IA*erf((possibleSizes)./(PR + RTW)) + B; 
%             MT.SuppressionIndex = (max(fit)-fit(end))/(max(fit)-B);
            
            FR = G * DirectionTune .* SpeedTune .* SizeTune + B0; 
            
            MT.FiringRate = FR;
            
            clear fit;
        end
        
        function SimulateVariance(MT)
            
            c = 1.5;
            MT.Variance = c * (MT.FiringRate);
            
        end
    end
    
    
end