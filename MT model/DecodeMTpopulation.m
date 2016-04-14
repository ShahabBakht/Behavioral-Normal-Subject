function [TargetEstimate] = DecodeMTpopulation(mtpopulation, R)

%%%%%%%%%%%%%%%%
epcilon = .01;
%%%%%%%%%%%%%%%%

for i = 1:length(mtpopulation)
    PS(i) = mtpopulation{i}.PreferredSpeed;
    PD(i) = mtpopulation{i}.PreferredDirection;
end

Rx = sum(repmat(cos(PD.*pi/180)',1,size(R,2)).*R,1);
Ry = sum(repmat(sin(PD.*pi/180)',1,size(R,2)).*R,1);

Ex = ( sum(repmat(cos(PD.*pi/180)',1,size(R,2)).*R.*repmat((PS)',1,size(R,2)),1) )./( epcilon + sqrt(Rx.^2 + Ry.^2));
Ey = ( sum(repmat(sin(PD.*pi/180)',1,size(R,2)).*R.*repmat((PS)',1,size(R,2)),1) )./( epcilon + sqrt(Rx.^2 + Ry.^2));

TargetEstimate.SPDest = 2.^( sqrt(Ex.^2 + Ey.^2) );
TargetEstimate.DIRest = atan( Ey./Ex ).* 180/pi;


end