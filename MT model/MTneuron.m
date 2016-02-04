classdef MTneuron < handle
    
    properties
        PreferredSpeed
        PreferredDirection
        SpeedTuningWidth
        DirectionTuningWidth
        Gain
        FiringRate
        Variance
    end
    
    methods
        function MT = MTneuron(S)
            MT.PreferredSpeed = S.PS;
            MT.PreferredDirection = S.PD;
            MT.SpeedTuningWidth = S.STW;
            MT.DirectionTuningWidth = S.DTW;
            MT.Gain = S.G;
        end
        
        function SimulateMeanFiringRate(MT,Target)
            % neurons tuning characteristics
            G = MT.Gain;
            PD = MT.PreferredDirection;
            PS = MT.PreferredSpeed;
            STW = (MT.SpeedTuningWidth);
            DTW = MT.DirectionTuningWidth;
            
            % target motion parameters
            Theta = Target.MotionDirection;
            S = Target.MotionSpeed;
            
            Theta = Theta * pi / 180;
            PD = PD * pi / 180;
            if (Theta - PD) > pi
                DDiff = -pi + ((Theta - PD) - pi);
            elseif (Theta - PD) < -pi
                DDiff = pi - (-pi - (Theta - PD));
            else
                DDiff = (Theta - PD);
            end
            
            FR = G * exp(-( DDiff.^2 )/( 2 * DTW^2 )) .* exp(-( (log2(S) - (PS)).^2 )./( 2 * STW^2 )); 
            
            MT.FiringRate = FR;
            
        end
        
        function SimulateVariance(MT)
            
            c = 1;
            MT.Variance = c * mean(MT.FiringRate);
            
        end
    end
    
    
end