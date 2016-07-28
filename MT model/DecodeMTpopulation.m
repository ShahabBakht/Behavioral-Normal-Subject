function [TargetEstimate] = DecodeMTpopulation(mtpopulation, R_numerator, R_denominator,method)

%%%%%%%%%%%%%%%%
epcilon = 1;
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
end


end