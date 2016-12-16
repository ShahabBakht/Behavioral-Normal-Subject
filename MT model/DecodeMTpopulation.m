function [TargetEstimate] = DecodeMTpopulation(mtpopulation, R_numerator, mtpopulation_denom, R_denominator,method, param)

%%%%%%%%%%%%%%%%
epcilon = .1;
%%%%%%%%%%%%%%%%

% for i = 1:length(mtpopulation)
%     PS(i) = mtpopulation{i}.PreferredSpeed;
%     PD(i) = mtpopulation{i}.PreferredDirection;
% end

switch method
    case 'num&denom_opp'
        PD = cellfun(@(x)(x.PreferredDirection),mtpopulation);
        PS = cellfun(@(x)(x.PreferredSpeed),mtpopulation);
        
        Rx = sum(repmat(cos(PD.*pi/180)',1,size(R_numerator,2)).*R_numerator,1);
        Ry = sum(repmat(sin(PD.*pi/180)',1,size(R_numerator,2)).*R_numerator,1);
        
        Ex = ( sum(repmat(cos(PD.*pi/180)',1,size(R_numerator,2)).*R_numerator.*repmat((PS)',1,size(R_numerator,2)),1) )./( epcilon + sqrt(Rx.^2 + Ry.^2));
        Ey = ( sum(repmat(sin(PD.*pi/180)',1,size(R_numerator,2)).*R_numerator.*repmat((PS)',1,size(R_numerator,2)),1) )./( epcilon + sqrt(Rx.^2 + Ry.^2));
        
        
        TargetEstimate.SPDest = 2.^( sqrt(Ex.^2 + Ey.^2) );
        TargetEstimate.DIRest = atan( Ey./Ex ).* 180/pi;
        
    case 'num_opp'
        PD = cellfun(@(x)(x.PreferredDirection),mtpopulation);
        PS = cellfun(@(x)(x.PreferredSpeed),mtpopulation);
        
        Ex = ( sum(repmat(cos(PD.*pi/180)',1,size(R_numerator,2)).*R_numerator.*repmat((PS)',1,size(R_numerator,2)),1) )./( epcilon + sum(R_numerator,1));
        Ey = ( sum(repmat(sin(PD.*pi/180)',1,size(R_numerator,2)).*R_numerator.*repmat((PS)',1,size(R_numerator,2)),1) )./( epcilon + sum(R_numerator,1));
        
        
        TargetEstimate.SPDest = 2.^( sqrt(Ex.^2 + Ey.^2) );
        TargetEstimate.DIRest = atan( Ey./Ex ).* 180/pi;
        
        
    case 'uncorr_norm'
        PD = cellfun(@(x)(x.PreferredDirection),mtpopulation);
        PS = cellfun(@(x)(x.PreferredSpeed),mtpopulation);
        
        Ex = ( sum(repmat(cos(PD.*pi/180)',1,size(R_numerator,2)).*R_numerator.*repmat((PS)',1,size(R_numerator,2)),1) )./( epcilon + abs(sum(R_denominator,1)));
        Ey = ( sum(repmat(sin(PD.*pi/180)',1,size(R_numerator,2)).*R_numerator.*repmat((PS)',1,size(R_numerator,2)),1) )./( epcilon + abs(sum(R_denominator,1)));
        
        
        TargetEstimate.SPDest = 2.^( sqrt(Ex.^2 + Ey.^2) );
        TargetEstimate.DIRest = atan( Ey./Ex ).* 180/pi;
        
    case 'weighted_uncorr_norm'
        PD = cellfun(@(x)(x.PreferredDirection),mtpopulation);
        PS = cellfun(@(x)(x.PreferredSpeed),mtpopulation);
%         RTW = cellfun(@(x)(x.SizeTuningWidthPD),mtpopulation);
        SI = cellfun(@(x)(x.SuppressionIndex),mtpopulation);
        RFLocation = cellfun(@(x)(x.RFLocation),mtpopulation);
        
        SIdenom = cellfun(@(x)(x.SuppressionIndex),mtpopulation_denom);
        RFLocationdenom = cellfun(@(x)(x.RFLocation),mtpopulation);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% 
        %%% window for Suppresson Index
%         Weight = nan(size(SI));
%         Weight(SI>quantile(SI,.5)) = 1;
%         Weight(SI<=quantile(SI,.5)) = 1;
%         
%         Weight_denom = nan(size(SIdenom));
%         Weight_denom(SIdenom>quantile(SIdenom,.5)) = 1;
%         Weight_denom(SIdenom<=quantile(SIdenom,.5)) = 1;
        
        %%% window for RF location
%         rf = param(1);
%         Weight = nan(size(RFLocation));
%         Weight(RFLocation<=rf) = 1;
%         Weight(RFLocation>rf) = 0;
%         
%         Weight_denom = nan(size(RFLocationdenom));
%         Weight_denom(RFLocationdenom<=rf) = 1;
%         Weight_denom(RFLocationdenom>rf) = 0;
        
        %%% window for RF location and Suppresson Index
%         rf = param(1);  
%         Weight = nan(size(RFLocation));
%         Weight(RFLocation<=rf & SI>quantile(SI,.5)) = 1;
%         Weight(RFLocation>rf | SI<=quantile(SI,.5)) = 0;
%         
%         Weight_denom = nan(size(RFLocationdenom));
%         Weight_denom(RFLocationdenom<=rf & SIdenom>quantile(SIdenom,.5)) = 1;
%         Weight_denom(RFLocationdenom>rf | SIdenom<=quantile(SIdenom,.5)) = 0;

        %%% window for RF location (hard) and Suppresson Index (soft parametric)
        
        rf = param(1);
        a = param(2);
        c = param(3);
        
        Weight_RF = nan(size(RFLocation));
        Weight_RF(RFLocation<=rf) = 1;
        Weight_RF(RFLocation>rf) = 0;
        Weight_SI = sigmf(SI,[a quantile(SIdenom,.5)]);
        Weight = Weight_RF.*Weight_SI;
        
        Weight_RF = nan(size(RFLocationdenom));
        Weight_RF(RFLocationdenom<=rf) = 1;
        Weight_RF(RFLocationdenom>rf) = 0;
        Weight_SI = sigmf(SIdenom,[a quantile(SIdenom,.5)]);
        Weight_denom = Weight_RF.*Weight_SI;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
      
        
%         Rx = sum(repmat(cos(PD.*pi/180)',1,size(R_numerator,2)).*R_numerator.*repmat((Weight)',1,size(R_numerator,2)),1);
%         Ry = sum(repmat(sin(PD.*pi/180)',1,size(R_numerator,2)).*R_numerator.*repmat((Weight)',1,size(R_numerator,2)),1);
        
%         Ex = ( sum(repmat(cos(PD.*pi/180)',1,size(R_numerator,2)).*repmat((Weight)',1,size(R_numerator,2)).*R_numerator.*repmat((PS)',1,size(R_numerator,2)),1) )./( epcilon + abs(sum(R_numerator.*repmat((Weight)',1,size(R_numerator,2)),1)) );
%         Ey = ( sum(repmat(sin(PD.*pi/180)',1,size(R_numerator,2)).*repmat((Weight)',1,size(R_numerator,2)).*R_numerator.*repmat((PS)',1,size(R_numerator,2)),1) )./( epcilon + abs(sum(R_numerator.*repmat((Weight)',1,size(R_numerator,2)),1)) );
        Ex = ( sum(repmat(cos(PD.*pi/180)',1,size(R_numerator,2)).*repmat((Weight)',1,size(R_numerator,2)).*R_numerator.*repmat((PS)',1,size(R_numerator,2)),1) )./...
            ( epcilon + abs(sum(repmat((Weight_denom)',1,size(R_denominator,2)).*R_denominator,1)));
        Ey = ( sum(repmat(sin(PD.*pi/180)',1,size(R_numerator,2)).*repmat((Weight)',1,size(R_numerator,2)).*R_numerator.*repmat((PS)',1,size(R_numerator,2)),1) )./...
            ( epcilon + abs(sum(repmat((Weight_denom)',1,size(R_denominator,2)).*R_denominator,1)));
        
       
        
        TargetEstimate.SPDest = 2.^( sqrt(Ex.^2 + Ey.^2) );
        TargetEstimate.DIRest = atan( Ey./Ex ).* 180/pi;
end


end