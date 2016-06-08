function p=fRatioTest(data,model1,coeff1,model2,coeff2)
err1=sum((model1-data).^2);
err2=sum((model2-data).^2);
n1=coeff2-coeff1;
n2=length(data)-coeff2;
fRatio=n2*(err1-err2)/(n1*err2);

% if fRatio<1
%     fRatio=1/fRatio;
% end

p=1-fcdf(fRatio,n1,n2);
