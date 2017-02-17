function [a,si,fitType,resnorm]=mtSizeFit(x,y,force)
% a(1) = excitation amplitude
% a(2) = excitation size
% a(3) = inhibition amplitude
% a(4) = inhibition - excitation size
% a(5) = baseline response
warning off
lastx=x(end);
indmax=find(y==max(y));
if indmax~=1
    maxx=x(indmax(end)-1);
else
    maxx=x(indmax(end));
end
if maxx==0
    maxx=0.1;   % hack to make the initial start non-zero
end
% ao1 = [max(y) maxx min(y)];
fitType=0;
options=optimset('Display','off');
% [a1,resnorm]=lsqcurvefit(@erfSize,ao1,x,y,[0 0 0],[1.5*(max(y)-min(y)) max(x) y(1)+0.01],options);

ao2 = [1.5*(max(y)-min(y)) maxx 1.5*(max(y)-min(y)) lastx y(1)];

[a2,resnorm]=lsqcurvefit(@diffGauss,ao2,x,y,[max(y)-min(y) 0 0 0 0],[1.5*(max(y)-min(y)) max(x) 1.5*(max(y)-min(y)) max(x) y(1)+0.01],options);

fit2=diffGauss(a2,x); % SI calculation from fit
% si=(max(fit2)-fit2(end))/(max(fit2)-0.5);
% si=(max(fit2)-fit2(end))/(max(fit2)); % for d'
si=(max(fit2)-fit2(end))/(max(fit2)-a2(5));

% si=(diffGauss(a2,1.163*a2(2))-diffGauss(a2,max(x)))/(diffGauss(a2,1.163*a2(2))-a2(5));
% if diffGauss(a2,max(x))>diffGauss(a2,1.163*a2(2))
%     si=0;
% end

% p=fRatioTest(y,erfSize(a1,x),length(a1),diffGauss(a2,x),length(a2));

% if p<0.05
%     a=a2;
% %     disp('DoE is a better fit')
%     fitType=1;
% else
%     a=a1;
% %     disp('DoE is NOT a better fit')
% end

if force==1
    a=a2;
%     disp('Forced to fit DoE')
end