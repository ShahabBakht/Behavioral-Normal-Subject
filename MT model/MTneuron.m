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
        SuppressionIndex
        Gain
        BaseLine
        FiringRate
        Variance
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
            MT.SuppressionIndex = S.SI;
            MT.Gain = S.G;
            MT.BaseLine = S.B0;
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
            DTW = MT.DirectionTuningWidth;
            RTW = MT.SizeTuningWidthPD;
            nRTW = MT.SizeTuningWidthND;
            
            % target motion parameters
            Theta = Target.MotionDirection;
            S = Target.MotionSpeed;
            R = Target.Size;
            
            
            Theta = Theta * pi / 180;
            PD = PD * pi / 180;
            if (Theta - PD) > pi
                DDiff = -pi + ((Theta - PD) - pi);
            elseif (Theta - PD) < -pi
                DDiff = pi - (-pi - (Theta - PD));
            else
                DDiff = (Theta - PD);
            end
            
            DirectionTune = exp(-( DDiff.^2 )/( 2 * DTW^2 ));
            SpeedTune = exp(-( (log2(S) - (PS)).^2 )./( 2 * STW^2 ));
            
            if ~isempty(PR) % if size tuning is defined for the neurons
                if DDiff <= pi/2
                    SizeTune = erf(R./PR) - .5*erf(R./(PR + RTW));
                else
                    SizeTune = erf(R./nPR) - .5*erf(R./(nPR + nRTW));
                end
            else
                SizeTune = 1;
            end
            
            FR = G * DirectionTune .* SpeedTune .* SizeTune + B0; 
            
            MT.FiringRate = FR;
            
        end
        
        function SimulateVariance(MT)
            
            c = 1;
            MT.Variance = c * mean(MT.FiringRate);
            
        end
    end
    
    
end