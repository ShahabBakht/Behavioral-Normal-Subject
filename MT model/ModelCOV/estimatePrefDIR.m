[FileName,PathName] = uigetfile('*.mat','Select the noise correlation data');
addpath('D:\Project Codes\Behavioral-Normal-Subject\MT model');
load([PathName,FileName]);
delta = .0001;
DIRlr = 0:pi/4:(2*pi);
DIRhr = 0:delta:(2*pi);
numPairs = size(noiseCorrDataDirectionTuning,1);
DIRpref = nan(numPairs,3);
for paircount = 1:numPairs
    FRlr1 = [noiseCorrDataDirectionTuning(paircount,4:11),noiseCorrDataDirectionTuning(paircount,4)];
    FRlr2 = [noiseCorrDataDirectionTuning(paircount,12:end),noiseCorrDataDirectionTuning(paircount,12)];
    FRhr1 = spline(DIRlr,FRlr1,DIRhr);
    FRhr2 = spline(DIRlr,FRlr2,DIRhr);
    c = corrcoef(noiseCorrDataDirectionTuning(paircount,4:11),noiseCorrDataDirectionTuning(paircount,12:end));
    
    
    [~,DIRprefIdx1] = max(FRhr1);
    [~,DIRprefIdx2] = max(FRhr2);
    
    if DIRhr(DIRprefIdx1) ~= 2*pi
        DIRpref(paircount,1) = DIRhr(DIRprefIdx1) * 180/pi;
    else
        DIRpref(paircount,1) = 0 * 180/pi;
    end
    
    if DIRhr(DIRprefIdx2) ~= 2*pi
        DIRpref(paircount,2) = DIRhr(DIRprefIdx2) * 180/pi;
    else
        DIRpref(paircount,2) = 0 * 180/pi;
    end
    
    DIRpref(paircount,3) = abs(AngDiff(DIRpref(paircount,1),DIRpref(paircount,2)));
    DIRpref(paircount,4) = c(2);
        
    
end

save('DIRpref.mat','DIRpref');


SS1 = noiseCorrDataDirectionTuning(:,2)>median(noiseCorrDataDirectionTuning(:,2));
SS2 = noiseCorrDataDirectionTuning(:,3)>median(noiseCorrDataDirectionTuning(:,3));

Data_bothSS = [DIRpref(SS1&SS2,4),noiseCorrDataDirectionTuning(SS1&SS2,1)];
Data_bothnSS = [DIRpref(~SS1&~SS2,4),noiseCorrDataDirectionTuning(~SS1&~SS2,1)];
Data_SSnSS = [DIRpref(xor(SS1,SS2),4),noiseCorrDataDirectionTuning(xor(SS1,SS2),1)];

save('Data_bothSS.mat','Data_bothSS');
save('Data_bothnSS.mat','Data_bothnSS');
save('Data_SSnSS.mat','Data_SSnSS');

CorrMatrix(:,1) = DIRpref(:,4);
CorrMatrix(:,2) = noiseCorrDataDirectionTuning(:,1);
CorrMatrix(:,3) = noiseCorrDataDirectionTuning(:,2);
CorrMatrix(:,4) = noiseCorrDataDirectionTuning(:,3);

%%

figure;plot(DIRpref(SS1&SS2,3),noiseCorrDataDirectionTuning(SS1&SS2,1),'.b','MarkerSize',15);hold on;
hold on;plot(DIRpref(~SS1&~SS2,3),noiseCorrDataDirectionTuning(~SS1&~SS2,1),'.r','MarkerSize',15);hold on;
hold on;plot(DIRpref(xor(SS1,SS2),3),noiseCorrDataDirectionTuning(xor(SS1,SS2),1),'.k','MarkerSize',15);hold on;
title('noise correlation vs. tuning difference')

figure;plot(DIRpref(SS1&SS2,4),noiseCorrDataDirectionTuning(SS1&SS2,1),'.b','MarkerSize',15);hold on;
hold on;plot(DIRpref(~SS1&~SS2,4),noiseCorrDataDirectionTuning(~SS1&~SS2,1),'.r','MarkerSize',15);hold on;
hold on;plot(DIRpref(xor(SS1,SS2),4),noiseCorrDataDirectionTuning(xor(SS1,SS2),1),'.k','MarkerSize',15);hold on;
title('noise correlation vs. signal correlation')

mdl = LinearModel.fit(DIRpref(SS1&SS2,3),noiseCorrDataDirectionTuning(SS1&SS2,1));
Coeffs_bothSS = mdl.Coefficients.Estimate;
rnoise_predict = predict(mdl,(-180:180)');
hold on;mdl.plot;

mdl = LinearModel.fit(DIRpref(~SS1&~SS2,3),noiseCorrDataDirectionTuning(~SS1&~SS2,1));
Coeffs_bothnSS = mdl.Coefficients.Estimate;
hold on;mdl.plot;

mdl = LinearModel.fit(DIRpref(xor(SS1,SS2),3),noiseCorrDataDirectionTuning(xor(SS1,SS2),1));
Coeffs_SSnSS = mdl.Coefficients.Estimate;
hold on;mdl.plot;

