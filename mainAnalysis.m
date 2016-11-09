%% Load the Subjects data
display('Choose the Raw Eye mat file:');uiopen;
display('Choose the Saccade mat file:');uiopen;
listOFvariables = who;
X = [];
S = [];
TV = [];
for i = 1:length(listOFvariables)
    switch listOFvariables{i}(1)
        case 'X'
            Command = ['X = cat(2,X,',listOFvariables{i},'(:,:,1:1804));'];
            eval(Command);
        case 'S'
            Command = ['S = cat(2,S,',listOFvariables{i},');'];
            eval(Command);
        case 'T'
            Command = ['TV = cat(2,TV,',listOFvariables{i},');'];
            eval(Command);
    end
end

% Victoire data
% Uthreshold_L = ?; Lthreshold_L = ?;
% load('D:\Analysis\Behavioral-Normal-Subject\Raw Eye/VictoireData.mat')
% load('D:\Analysis\Behavioral-Normal-Subject\Saccades/VictoireSaccades.mat')
% X = cat(2, ...
%     X16042704(:,:,1:1805), ...
%     X16042705(:,:,1:1805), ...
%     X16050503(:,:,1:1805), ...
%     X16050504(:,:,1:1805), ...
%     X16050606(:,:,1:1805), ...
%     X16050607(:,:,1:1805), ...
%     X16051005(:,:,1:1805), ...
%     X16051006(:,:,1:1805) ...
%     );
% S = cat(2, ...
%     S16042704, ...
%     S16042705, ...
%     S16050503, ...
%     S16050504, ...
%     S16050606, ...
%     S16050607, ...
%     S16051005, ...
%     S16051006 ...
%     );
% TV = cat(2, ...
%     TV16042704, ...
%     TV16042705, ...
%     TV16050503, ...
%     TV16050504, ...
%     TV16050606, ...
%     TV16050607, ...
%     TV16051005, ...
%     TV16051006 ...
%     );
% 
% Travis data
% Uthreshold_L = ?; Lthreshold_L = ?;
% load('D:\Analysis\Behavioral-Normal-Subject\Raw Eye/TravisData.mat');
% load('D:\Analysis\Behavioral-Normal-Subject\Saccades/TravisSaccades.mat')
% X = cat(2, ...
%     X16042602(:,:,1:1805), ...
%     X16050501(:,:,1:1805), ...
%     X16050502(:,:,1:1805), ...
%     X16050604(:,:,1:1805), ...
%     X16050605(:,:,1:1805), ...
%     X16051003(:,:,1:1805), ...
%     X16051004(:,:,1:1805) ...
%     );
% S = cat(2, ...
%     S16042602, ...
%     S16050501, ...
%     S16050502, ...
%     S16050604, ...
%     S16050605, ...
%     S16051003, ...
%     S16051004 ...
%     );
% TV = cat(2, ...
%     TV16042602, ...
%     TV16050501, ...
%     TV16050502, ...
%     TV16050604, ...
%     TV16050605, ...
%     TV16051003, ...
%     TV16051004 ...
%     );
% % 
% Claire data
% Uthreshold_L = ?; Lthreshold_L = ?;
% load('/Users/shahab/MNI/Analysis/Normal-Subject-Behavior/RawEye/ClaireData.mat');
% load('/Users/shahab/MNI/Analysis/Normal-Subject-Behavior/Saccades/ClaireSaccades.mat');
% X = cat(2, ...
%     X16042908(:,:,1:1804), ...
%     X16050201(:,:,1:1804), ...
%     X16050202(:,:,1:1804), ...
%     X16050301(:,:,1:1804), ...
%     X16050302(:,:,1:1804), ...
%     X16051001(:,:,1:1804), ...
%     X16051002(:,:,1:1804) ...
%     );
% S = cat(2, ...
%     S16042908, ...
%     S16050201, ...
%     S16050202, ...
%     S16050301, ...
%     S16050302, ...
%     S16051001, ...
%     S16051002 ...
%     );
% TV = cat(2, ...
%     TV16042908, ...
%     TV16050201, ...
%     TV16050202, ...
%     TV16050301, ...
%     TV16050302, ...
%     TV16051001, ...
%     TV16051002 ...
%     );
% 
% Lucy data
% Uthreshold_L = 0; Lthreshold_L = -30;
% load('/Users/shahab/MNI/Analysis/Normal-Subject-Behavior/RawEye/LucyData.mat');
% load('/Users/shahab/MNI/Analysis/Normal-Subject-Behavior/Saccades/LucySaccades.mat');
% X = cat(2, ...
%     X16031608(:,:,1:1805), ...
%     X16032101(:,:,1:1805), ...
%     X16032102(:,:,1:1805), ...
%     X16032301(:,:,1:1805), ...
%     X16032302(:,:,1:1805), ...
%     X16033001(:,:,1:1805), ...
%     X16033002(:,:,1:1805) ...
%     );
% S = cat(2, ...
%     S16031608, ...
%     S16032101, ...
%     S16032102, ...
%     S16032301, ...
%     S16032302, ...
%     S16033001, ...
%     S16033002 ...
%     );
% TV = cat(2, ...
%     TV16031608, ...
%     TV16032101, ...
%     TV16032102, ...
%     TV16032301, ...
%     TV16032302, ...
%     TV16033001, ...
%     TV16033002 ...
%     );

% Helga data
% load('/Users/shahab/MNI/Analysis/Normal-Subject-Behavior/RawEye/HelgaData.mat');
% load('/Users/shahab/MNI/Analysis/Normal-Subject-Behavior/Saccades/HelgaSaccades.mat');
% X = cat(2, ...
%     X16031104(:,:,1:1807), ...
%     X16031705(:,:,1:1807), ...
%     X16031706(:,:,1:1807), ...
%     X16032203(:,:,1:1807), ...
%     X16032204(:,:,1:1807), ...
%     X16033101(:,:,1:1807), ...
%     X16033102(:,:,1:1807) ...
%     );
% S = cat(2, ...
%     S16031104, ...
%     S16031705, ...
%     S16031706, ...
%     S16032203, ...
%     S16032204, ...
%     S16033101, ...
%     S16033102 ...
%     );
% TV = cat(2, ...
%     TV16031104, ...
%     TV16031705, ...
%     TV16031706, ...
%     TV16032203, ...
%     TV16032204, ...
%     TV16033101, ...
%     TV16033102 ...
%     );
% 
% Melissa data
% Uthreshold_L = 5; Lthreshold_L = -30;
% load('/Users/shahab/MNI/Analysis/Normal-Subject-Behavior/RawEye/MelissaData.mat');
% load('/Users/shahab/MNI/Analysis/Normal-Subject-Behavior/Saccades/MelissaSaccades.mat');
% X = cat(2, ...
%     X16031408(:,:,1:1808), ...
%     X16031707(:,:,1:1808), ...
%     X16031708(:,:,1:1808), ...
%     X16032205(:,:,1:1808), ...
%     X16032206(:,:,1:1808), ...
%     X16032903(:,:,1:1808), ...
%     X16032904(:,:,1:1808) ...
%     );
% S = cat(2, ...
%     S16031408, ...
%     S16031707, ...
%     S16031708, ...
%     S16032205, ...
%     S16032206, ...
%     S16032903, ...
%     S16032904 ...
%     );
% TV = cat(2, ...
%     TV16031408, ...
%     TV16031707, ...
%     TV16031708, ...
%     TV16032205, ...
%     TV16032206, ...
%     TV16032903, ...
%     TV16032904 ...
%     );
% 
% Arnaud data
% Uthreshold_L = 0; Lthreshold_L = -25;
% load('D:\Analysis\Behavioral-Normal-Subject\Raw Eye/ArnaudData.mat');
% load('D:\Analysis\Behavioral-Normal-Subject\Saccades/ArnaudSaccades.mat');
% X = cat(2, ...
%     X16031003(:,:,1:1806), ...
%     X16031404(:,:,1:1806), ...
%     X16031504(:,:,1:1806), ...
%     X16031505(:,:,1:1806), ...
%     X16031604(:,:,1:1806), ...
%     X16031605(:,:,1:1806) ...
%     );
% S = cat(2, ...
%     S16031003, ...
%     S16031404, ...
%     S16031504, ...
%     S16031505, ...
%     S16031604, ...
%     S16031605 ...
%     );
% TV = cat(2, ...
%     TV16031003, ...
%     TV16031404, ...
%     TV16031504, ...
%     TV16031505, ...
%     TV16031604, ...
%     TV16031605 ...
%     );
% 
% Moneek data



% Gabriel data 
% Threshold: Uthreshold_L = 5;Lthreshold_L = -30;
% load('/Users/shahab/MNI/Analysis/Normal-Subject-Behavior/RawEye/GabrielData.mat');
% X = cat(2, ...
%     X2501(:,:,1:1806), ...
%     X2502(:,:,1:1806), ...
%     X2601(:,:,1:1806), ...
%     X2602(:,:,1:1806), ...
%     X0201(:,:,1:1806), ...
%     X0202(:,:,1:1806) ...
%     );
% S = cat(2, ...
%     S2501, ...
%     S2502, ...
%     S2601, ...
%     S2602, ...
%     S0201, ...
%     S0202 ...
%     );
% TV = cat(2, ...
%     TV2501, ...
%     TV2502, ...
%     TV2601, ...
%     TV2602, ...
%     TV0201, ...
%     TV0202 ...
%     );

% Azar data
% load('/Users/shahab/MNI/Analysis/Normal-Subject-Behavior/RawEye/AzarData.mat');
% X = cat(2 ...
%     ,X1601(:,:,1:1806) ...
%     ,X1602(:,:,1:1806) ...
%     ,X1701(:,:,1:1806) ...
%     ,X1702(:,:,1:1806) ...
%     ,X1901(:,:,1:1806) ...
%     ,X1902(:,:,1:1806) ...
%     );
% S = cat(2 ...
%     ,S1601 ...
%     ,S1602 ...
%     ,S1701 ...
%     ,S1702 ...
%     ,S1901 ...
%     ,S1902 ...
%     );
% TV = cat(2 ...
%     ,TV1601 ...
%     ,TV1602 ...
%     ,TV1701 ...
%     ,TV1702 ...
%     ,TV1901 ...
%     ,TV1902 ...
%     );
% 
% X = X1202(:,:,1:1806);
% X(:,31:60,:) = X1203(:,:,1:1806);
% X(:,61:90,:) = X1301(:,:,1:1806);
% X(:,91:120,:) = X1302(:,:,1:1806);
% X = X1301(:,:,1:1806);
% X(:,31:60,:) = X1302(:,:,1:1806);

% Shahab data
% load('/Users/shahab/MNI/Analysis/Normal-Subject-Behavior/RawEye/ShahabData.mat');
% X = cat(2 ...
%     ,X1801(:,:,1:1805) ...
%     ,X1802(:,:,1:1805) ...
%     ,X1803(:,:,1:1805) ...
%     ,X1804(:,:,1:1805) ...
%     ,X1805(:,:,1:1805) ...
%     ,X2001(:,:,1:1805) ...
%     ,X2002(:,:,1:1805) ...
%     );
% S = cat(2 ...
%     ,S1801 ...
%     ,S1802 ...
%     ,S1803 ...
%     ,S1804 ...
%     ,S1805 ...
%     ,S2001 ...
%     ,S2002 ...
%     );
% TV = cat(2 ...
%     ,TV1801 ...
%     ,TV1802 ...
%     ,TV1803 ...
%     ,TV1804 ...
%     ,TV1805 ...
%     ,TV2001 ...
%     ,TV2002 ...
%     );
% 
% my old data
% X = X1901(:,:,1:1805);
% X(:,31:60,:) = X1902(:,:,1:1805);
% X(:,61:90,:) = X1903(:,:,1:1805);
% X(:,91:120,:) = X1904(:,:,1:1805);
% X(:,121:150,:) = X1905(:,:,1:1805);
% X(:,211:240,:) = X2003(:,:,1:1805);
% X(:,151:180,:) = X1406(:,:,1:1808);
% X = X1103(:,:,1:1806);
% X(:,31:60,:) = X1104(:,:,1:1806);
% X(:,61:90,:) = X1105(:,:,1:1806);
% X(:,91:120,:) = X1106(:,:,1:1806);
% X(:,121:150,:) = X1107(:,:,1:1806);
% X(:,151:180,:) = X1201(:,:,1:1806);
% S = S1202;
% S(:,31:60,:) = S1203;
% S(:,61:90,:) = S1301;
% S(:,91:120,:) = S1302;
% S = S1301;
% S(:,31:60,:) = S1302;
% 
% S = S1901;
% S(:,31:60,:) = S1902;
% S(:,61:90,:) = S1903;
% S(:,91:120,:) = S1904;
% S(:,121:150,:) = S1905;
% S(:,181:210,:) = S2002;
% S(:,211:240,:) = S2003;
% S(:,151:180,:) = S1806;
% S = S1103;
% S(:,31:60,:) = S1104;
% S(:,61:90,:) = S1105;
% S(:,91:120,:) = S1106;
% S(:,121:150,:) = S1107;
% S(:,151:180,:) = S1201;


% TV = TV1202;
% TV(:,31:60) = TV1203;
% TV(:,61:90) = TV1301;
% TV(:,91:120) = TV1302;
% TV = TV1301;
% TV(:,31:60) = TV1302;
%
% TV = TV1901;
% TV(:,31:60) = TV1902;
% TV(:,61:90) = TV1903;
% TV(:,91:120) = TV1904;
% TV(:,121:150) = TV1905;
% 
% TV(:,181:210) = TV2002;
% TV(:,211:240) = TV2003;

% TV(:,151:180) = TV1806;
% TV = TV1103;
% TV(:,31:60) = TV1104;
% TV(:,61:90) = TV1105;
% TV(:,91:120) = TV1106;
% TV(:,121:150) = TV1107;
% TV(:,151:180) = TV1201;



%% Calculate initial velocity

for c = 1:6
    for tr = 1:size(X,2)
        x = squeeze(X(c,tr,:));
        Send = floor(S(c,tr,3));
        Sstart = floor(S(c,tr,2));
        if ~isnan(Sstart) 
        if length(x) < Send + 120
            x_posS = [];
            v_posS = [];
            X_posS(c,tr,:) = nan;
            V_posS(c,tr,:) = nan;
        else
            x_posS = x((Send + 40):(Send + 110));
%             plot((Send + 20):(Send + 120),x_posS,'r');hold on;plot(x,'--k');pause;close
            v_posS = (x_posS(end) - x_posS(1))./0.07;
            X_posS(c,tr,:) = x_posS;
            V_posS(c,tr,:) = v_posS;
        end
        
        x_preS = x((Sstart - 100):(Sstart - 20));
        v_preS = (x_preS(end) - x_preS(1))./0.08;
        
%         y = squeeze(Y(c,tr,:));
%         y_posS = y((Send + 40):(Send + 100));
%         y_preS = y((Sstart - 100):(Sstart - 20));
%         
        
        X_preS(c,tr,:) = x_preS;
        V_preS(c,tr,:) = v_preS;
        
%         Y_posS(c,tr,:) = y_posS;
%         Y_preS(c,tr,:) = y_preS;
        
        clear x x_posS x_preS v_posS v_preS Send Sstart %y y_posS y_preS 
        else
            V_posS(c,tr,:) = nan;
        end
        
    end 
end




%% plot 
                        
% 
figure;
plot(TV(1,:),V_posS(1,:),'.r','MarkerSize',10);hold on;
plot(TV(3,:),V_posS(3,:),'.b','MarkerSize',10);
plot(TV(5,:),V_posS(5,:),'.g','MarkerSize',10);
title('Rightward Pursuit')


figure;
plot(TV(2,:),-V_posS(2,:),'.r','MarkerSize',10);hold on;
plot(TV(4,:),-V_posS(4,:),'.b','MarkerSize',10);
plot(TV(6,:),-V_posS(6,:),'.g','MarkerSize',10);
title('Leftward Pursuit')
% 

%% Throwing out outliers
% Leftward
Uthreshold_L = 5;
Lthreshold_L = -15;
%Rightward
Uthreshold_R = 20;
Lthreshold_R = -10;

UThresholdMatrix(1:2:5,:) = Uthreshold_R * ones(3,size(V_posS,2));
UThresholdMatrix(2:2:6,:) = Uthreshold_L * ones(3,size(V_posS,2));
LThresholdMatrix(1:2:5,:) = Lthreshold_R * ones(3,size(V_posS,2));
LThresholdMatrix(2:2:6,:) = Lthreshold_L * ones(3,size(V_posS,2));


V_posS(V_posS > UThresholdMatrix) = nan;
V_posS(V_posS < LThresholdMatrix) = nan;



% post threshold plot
figure;
plot(TV(1,:),V_posS(1,:),'.r','MarkerSize',10);hold on;
plot(TV(3,:),V_posS(3,:),'.b','MarkerSize',10);
plot(TV(5,:),V_posS(5,:),'.g','MarkerSize',10);
title('Rightward Pursuit')


figure;
plot(TV(2,:),-V_posS(2,:),'.r','MarkerSize',10);hold on;
plot(TV(4,:),-V_posS(4,:),'.b','MarkerSize',10);
plot(TV(6,:),-V_posS(6,:),'.g','MarkerSize',10);
title('Leftward Pursuit')


%% prediction
% 
% Ydata6 = [V_posS(1,:),-V_posS(2,:)]';
% Xdata6 = [TV(1,:),TV(2,:)]';
Ydata6 = [-V_posS(2,:)]';
Xdata6 = [TV(2,:)]';
% Ydata6 = [V_posS(1,:)]';
% Xdata6 = [TV(1,:)]';
Xdata6(isnan(Ydata6)) = [];
Ydata6(isnan(Ydata6)) = [];
csvwrite('D:\Project Codes\Behavioral-Normal-Subject\Y6.csv',Ydata6);
csvwrite('D:\Project Codes\Behavioral-Normal-Subject\X6.csv',Xdata6);

% Ydata10 = [V_posS(3,:),-V_posS(4,:)]';
% Xdata10 = [TV(3,:),TV(4,:)]';
Ydata10 = [-V_posS(4,:)]';
Xdata10 = [TV(4,:)]';
% Ydata10 = [V_posS(3,:)]';
% Xdata10 = [TV(3,:)]';
Xdata10(isnan(Ydata10)) = [];
Ydata10(isnan(Ydata10)) = [];
csvwrite('D:\Project Codes\Behavioral-Normal-Subject\Y10.csv',Ydata10);
csvwrite('D:\Project Codes\Behavioral-Normal-Subject\X10.csv',Xdata10);


% Ydata20 = [V_posS(5,:),-V_posS(6,:)]';
% Xdata20 = [TV(5,:),TV(6,:)]';
Ydata20 = [-V_posS(6,:)]';
Xdata20 = [TV(6,:)]';
% Ydata20 = [V_posS(5,:)]';
% Xdata20 = [TV(5,:)]';
Xdata20(isnan(Ydata20)) = [];
Ydata20(isnan(Ydata20)) = [];
csvwrite('D:\Project Codes\Behavioral-Normal-Subject\Y20.csv',Ydata20);
csvwrite('D:\Project Codes\Behavioral-Normal-Subject\X20.csv',Xdata20);


% Ydata20 = [V_posS(7,:),-V_posS(8,:)]';
% Xdata20 = [TV(7,:),TV(8,:)]';
% Ydata20 = [-V_posS(8,:)]';
% Xdata20 = [TV(8,:)]';
% Ydata20 = [V_posS(7,:)]';
% Xdata20 = [TV(7,:)]';

% csvwrite('D:\Project Codes\Behavioral-Normal-Subject\Y20.csv',Ydata20);
% csvwrite('D:\Project Codes\Behavioral-Normal-Subject\X20.csv',Xdata20);

%% PCA and explained variance
% [PCA2.coeff,PCA2.score,PCA2.latent,PCA2.tsquared,PCA2.explained,PCA2.mu] = pca([Xdata6,Ydata6]);
% [PCA6.coeff,PCA6.score,PCA6.latent,PCA6.tsquared,PCA6.explained,PCA6.mu] = pca([Xdata10,Ydata10]);
% [PCA20.coeff,PCA20.score,PCA20.latent,PCA20.tsquared,PCA20.explained,PCA20.mu] = pca([Xdata20,Ydata20]);
% 
% 
% 
% listPCA = ls('D:\Analysis\Behavioral-Normal-Subject\PCAs\');
% sofarPCAs = (size(listPCA,1) - 2)./3;
% 
% save(['D:\Analysis\Behavioral-Normal-Subject\PCAs\PCA2_',num2str(sofarPCAs+1),'.mat'],'PCA2');
% save(['D:\Analysis\Behavioral-Normal-Subject\PCAs\PCA6_',num2str(sofarPCAs+1),'.mat'],'PCA6');
% save(['D:\Analysis\Behavioral-Normal-Subject\PCAs\PCA20_',num2str(sofarPCAs+1),'.mat'],'PCA20');

%%
% Leftward
% V_posS_L = V_posS(2:2:6,:)';
% V_posS_L = reshape(V_posS_L,size(V_posS_L,1)*size(V_posS_L,2),1);
% TV_L = TV(2:2:6,:)';
% TV_L = reshape(TV_L,size(V_posS_L,1)*size(V_posS_L,2),1);
% Cond = cell(size(V_posS_L,1)*size(V_posS_L,2),1);
% for i = 1:length(Cond)
%     if i <=size(V_posS_L,1)*size(V_posS_L,2)/3
%         Cond{i} = '2 degrees';
%     elseif i>size(V_posS_L,1)*size(V_posS_L,2)/3 && i<=size(V_posS_L,1)*size(V_posS_L,2)*2/3
%         Cond{i} = '6 degrees';
%     else
%         Cond{i} = '20 degrees';
%     end
% end

% scatterhist(TV_L,-V_posS_L,'Group',Cond,'Location','SouthEast',...
% 'Direction','out','Color','gbr','LineStyle',{'-'},...
% 'LineWidth',[2,2,2],'Marker','do*','MarkerSize',[6,6,6],'Kernel','on');
% scatterhist(TV_L,-V_posS_L,'Group',Cond,'LineStyle',{'-'},'Marker','.','MarkerSize',20);

% Rightward
% V_posS_L = V_posS(1:2:5,:)';
% V_posS_L = reshape(V_posS_L,540,1);
% TV_L = TV(1:2:5,:)';
% TV_L = reshape(TV_L,540,1);
% Cond = cell(540,1);
% for i = 1:length(Cond)
%     if i <=180
%         Cond{i} = '2 degrees';
%     elseif i>180 && i<=360
%         Cond{i} = '6 degrees';
%     else
%         Cond{i} = '20 degrees';
%     end
% end
% 
% scatterhist(TV_L,V_posS_L,'Group',Cond,'LineStyle',{'-'});

%% Saccades analysis
% Samp = squeeze(S(:,:,1));
% Samp_n = (Samp-min(min(Samp)))./(max(max(Samp))-min(min(Samp)));
% % Samp_n = Samp;
% TV_n = (TV-min(min(TV)))./(max(max(TV))-min(min(TV))); 
% % TV_n = TV;
% 
% SampL_n = Samp_n(:,:)';
% SampL_n = reshape(SampL_n,size(SampL_n,1)*size(SampL_n,2),1);
% TV_L = TV_n(:,:)';
% TV_L = reshape(TV_L,size(SampL_n,1)*size(SampL_n,2),1);
% Cond = cell(size(SampL_n,1)*size(SampL_n,2),1);
% for i = 1:length(Cond)
%     if i <=size(SampL_n,1)*size(SampL_n,2)/3
%         Cond{i} = '2 degrees';
%     elseif i>size(SampL_n,1)*size(SampL_n,2)/3 && i<=size(SampL_n,1)*size(SampL_n,2)*2/3
%         Cond{i} = '6 degrees';
%     else
%         Cond{i} = '20 degrees';
%     end
% end
% scatterhist(TV_L,SampL_n,'Group',Cond,'LineStyle',{'-'},'Marker','.','MarkerSize',20);
% 
% 
% figure;scatter(TV(1,:),squeeze(S(1,:,2)),'or');
% hold on;
% scatter(TV(3,:),squeeze(S(3,:,2)),'*k');
% scatter(TV(5,:),squeeze(S(5,:,2)),'xb');
% scatter(TV(2,:),squeeze(S(2,:,2)),'or');
% scatter(TV(4,:),squeeze(S(4,:,2)),'*k');
% scatter(TV(6,:),squeeze(S(6,:,2)),'xb');
% xlabel('target V');ylabel('saccade delay')
% 
% figure;scatter(squeeze(S(1,:,1)),squeeze(S(1,:,2)),'or');
% hold on;
% scatter(squeeze(S(3,:,1)),squeeze(S(3,:,2)),'*k');
% scatter(squeeze(S(5,:,1)),squeeze(S(5,:,2)),'xb');
% scatter(squeeze(S(2,:,1)),squeeze(S(2,:,2)),'or');
% scatter(squeeze(S(4,:,1)),squeeze(S(4,:,2)),'*k');
% scatter(squeeze(S(6,:,1)),squeeze(S(6,:,2)),'xb');
% xlabel('saccade amp');ylabel('saccade delay')
% 
% figure;scatter(TV(1,:),squeeze(S(1,:,1))*1000./(squeeze(S(1,:,3))-squeeze(S(1,:,2))),'or');
% hold on;
% scatter(TV(3,:),squeeze(S(3,:,1))*1000./(squeeze(S(3,:,3))-squeeze(S(3,:,2))),'*k');
% scatter(TV(5,:),squeeze(S(5,:,1))*1000./(squeeze(S(5,:,3))-squeeze(S(5,:,2))),'xb');
% scatter(TV(2,:),squeeze(S(5,:,1))*1000./(squeeze(S(2,:,3))-squeeze(S(2,:,2))),'or');
% scatter(TV(4,:),squeeze(S(5,:,1))*1000./(squeeze(S(4,:,3))-squeeze(S(4,:,2))),'*k');
% scatter(TV(6,:),squeeze(S(5,:,1))*1000./(squeeze(S(6,:,3))-squeeze(S(6,:,2))),'xb');
% 
% figure;scatter(V_posS(1,:),squeeze(S(1,:,1)),'or');
% hold on;
% scatter(V_posS(3,:),squeeze(S(3,:,1)),'*k');
% scatter(V_posS(5,:),squeeze(S(5,:,1)),'xb');
% scatter(-V_posS(2,:),squeeze(S(2,:,1)),'or');
% scatter(-V_posS(4,:),squeeze(S(4,:,1)),'*k');
% scatter(-V_posS(6,:),squeeze(S(6,:,1)),'xb');
% xlabel('target V');ylabel('saccade amp')
% 
% SVC = (squeeze(S(:,:,1))-10)*1000./(squeeze(S(:,:,3))-1000);
% SVC_L = SVC(:,:)';
% SVC_L = reshape(SVC_L,size(SVC_L,1)*size(SVC_L,2),1);
% TV_L = TV(:,:)';
% TV_L = reshape(TV_L,size(SVC_L,1)*size(SVC_L,2),1);
% Cond = cell(size(SVC_L,1)*size(SVC_L,2),1);
% for i = 1:length(Cond)
%     if i <=size(SVC_L,1)*size(SVC_L,2)/3
%         Cond{i} = '2 degrees';
%     elseif i>size(SVC_L,1)*size(SVC_L,2)/3 && i<=size(SVC_L,1)*size(SVC_L,2)*2/3
%         Cond{i} = '6 degrees';
%     else
%         Cond{i} = '20 degrees';
%     end
% end
% scatterhist(TV_L,SVC_L,'Group',Cond,'LineStyle',{'-'},'Marker','.','MarkerSize',20,'LineWidth',2);
% figure;scatter(TV(1,:),SVC(1,:),'or');hold on
% scatter(TV(3,:),SVC(3,:),'*k');
% scatter(TV(5,:),SVC(5,:),'xb');
% 
% Ydata6 = [SVC(1,:),SVC(2,:)]';
% Xdata6 = [TV(1,:),TV(2,:)]';
% Xdata6(isnan(Ydata6)) = [];
% Ydata6(isnan(Ydata6)) = [];
% 
% Ydata10 = [SVC(3,:),SVC(4,:)]';
% Xdata10 = [TV(3,:),TV(4,:)]';
% Xdata10(isnan(Ydata10)) = [];
% Ydata10(isnan(Ydata10)) = [];
% 
% Ydata20 = [SVC(5,:),SVC(6,:)]';
% Xdata20 = [TV(5,:),TV(6,:)]';
% Xdata20(isnan(Ydata20)) = [];
% Ydata20(isnan(Ydata20)) = [];
% 
% [c6,p6,rlo6,rup6] = corrcoef(Ydata6,Xdata6);c6 = c6(2);p6 = p6(2);rlo6 = rlo6(2);rup6 = rup6(2);
% [c10,p10,rlo10,rup10] = corrcoef(Ydata10,Xdata10);c10 = c10(2);p10 = p10(2);rlo10 = rlo10(2);rup10 = rup10(2);
% [c20,p20,rlo20,rup20] = corrcoef(Ydata20,Xdata20);c20 = c20(2);p20 = p20(2);rlo20 = rlo20(2);rup20 = rup20(2);
% figure;plot([2,6,20],[c6,c10,c20],'-o');hold on;plot([2,6,20],[rlo6,rlo10,rlo20],'--k');plot([2,6,20],[rup6,rup10,rup20],'--k');
% figure;plot([2,6,20],[p6,p10,p20],'-o');

%% percentile 

% SVC
% SVC = (squeeze(S(:,:,1))-10)*1000./(squeeze(S(:,:,3))-1000);
% prctlSVC_2 = [[prctile([TV(1,:),TV(2,:)],20), ...
% prctile([TV(1,:),TV(2,:)],40), ...
% prctile([TV(1,:),TV(2,:)],60), ...
% prctile([TV(1,:),TV(2,:)],80), ...
% prctile([TV(1,:),TV(2,:)],100)];...
% [prctile([SVC(1,:),SVC(2,:)],20), ...
% prctile([SVC(1,:),SVC(2,:)],40), ...
% prctile([SVC(1,:),SVC(2,:)],60), ...
% prctile([SVC(1,:),SVC(2,:)],80), ...
% prctile([SVC(1,:),SVC(2,:)],100)]];
% 
% prctlSVC_6 = [[prctile([TV(3,:),TV(4,:)],20), ...
% prctile([TV(3,:),TV(4,:)],40), ...
% prctile([TV(3,:),TV(4,:)],60), ...
% prctile([TV(3,:),TV(4,:)],80), ...
% prctile([TV(3,:),TV(4,:)],100)];...
% [prctile([SVC(3,:),SVC(4,:)],20), ...
% prctile([SVC(3,:),SVC(4,:)],40), ...
% prctile([SVC(3,:),SVC(4,:)],60), ...
% prctile([SVC(3,:),SVC(4,:)],80), ...
% prctile([SVC(3,:),SVC(4,:)],100)]];
% 
% prctlSVC_20 = [[prctile([TV(5,:),TV(6,:)],20), ...
% prctile([TV(5,:),TV(6,:)],40), ...
% prctile([TV(5,:),TV(6,:)],60), ...
% prctile([TV(5,:),TV(6,:)],80), ...
% prctile([TV(5,:),TV(6,:)],100)];...
% [prctile([SVC(5,:),SVC(6,:)],20), ...
% prctile([SVC(5,:),SVC(6,:)],40), ...
% prctile([SVC(5,:),SVC(6,:)],60), ...
% prctile([SVC(5,:),SVC(6,:)],80), ...
% prctile([SVC(5,:),SVC(6,:)],100)]];


% SPEM velocity
% prctlSPEMv_2 = [[prctile([TV(2,:)],20), ...
%     prctile([TV(2,:)],40), ...
%     prctile([TV(2,:)],60), ...
%     prctile([TV(2,:)],80), ...
%     prctile([TV(2,:)],100)];...
%     [prctile([-V_posS(2,:)],20), ...
%     prctile([-V_posS(2,:)],40), ...
%     prctile([-V_posS(2,:)],60), ...
%     prctile([-V_posS(2,:)],80), ...
%     prctile([-V_posS(2,:)],100)]];
% 
% prctlSPEMv_6 = [[prctile([TV(4,:)],20), ...
%     prctile([TV(4,:)],40), ...
%     prctile([TV(4,:)],60), ...
%     prctile([TV(4,:)],80), ...
%     prctile([TV(4,:)],100)];...
%     [prctile([-V_posS(4,:)],20), ...
%     prctile([-V_posS(4,:)],40), ...
%     prctile([-V_posS(4,:)],60), ...
%     prctile([-V_posS(4,:)],80), ...
%     prctile([-V_posS(4,:)],100)]];
% 
% prctlSPEMv_20 = [[prctile([TV(6,:)],20), ...
%     prctile([TV(6,:)],40), ...
%     prctile([TV(6,:)],60), ...
%     prctile([TV(6,:)],80), ...
%     prctile([TV(6,:)],100)];...
%     [prctile([-V_posS(6,:)],20), ...
%     prctile([-V_posS(6,:)],40), ...
%     prctile([-V_posS(6,:)],60), ...
%     prctile([-V_posS(6,:)],80), ...
%     prctile([-V_posS(6,:)],100)]];
% 



%% plot percentiles

% subjectsList = {'ag','az','cs','gc','hr','ls','mp','sb','tc','vs'};
% figure;
% for i = 1:length(subjectsList)
%     subplot(4,3,i);
%     Command = ['plot(prctlSPEMv_',subjectsList{i},'2(1,:),prctlSPEMv_',subjectsList{i},'2(2,:),''-ob'',''LineWidth'',2);hold on;'];
%     eval(Command);
%     Command = ['plot(prctlSPEMv_',subjectsList{i},'6(1,:),prctlSPEMv_',subjectsList{i},'6(2,:),''-or'',''LineWidth'',2);hold on;'];
%     eval(Command);
%     Command = ['plot(prctlSPEMv_',subjectsList{i},'20(1,:),prctlSPEMv_',subjectsList{i},'20(2,:),''-ok'',''LineWidth'',2);hold on;'];
%     eval(Command);
%     title(['subject: ',subjectsList{i},]);
%     xlabel('target velocity');ylabel('spem velocity')
% end
% 
% figure;
% for i = 1:length(subjectsList)
%     subplot(4,3,i);
%     Command = ['plot(prctlSVC_',subjectsList{i},'2(1,:),prctlSVC_',subjectsList{i},'2(2,:),''-ob'',''LineWidth'',2);hold on;'];
%     eval(Command);
%     Command = ['plot(prctlSVC_',subjectsList{i},'6(1,:),prctlSVC_',subjectsList{i},'6(2,:),''-or'',''LineWidth'',2);hold on;'];
%     eval(Command);
%     Command = ['plot(prctlSVC_',subjectsList{i},'20(1,:),prctlSVC_',subjectsList{i},'20(2,:),''-ok'',''LineWidth'',2);hold on;'];
%     eval(Command);
%     title(['subject: ',subjectsList{i},]);
%     xlabel('target velocity');ylabel('svc')
% end
% 
%% Variance calculations

% thisSubject = 'gc';
% SVC = (squeeze(S(:,:,1))-10)*1000./(squeeze(S(:,:,3))-1000);
% load('D:\Analysis\Behavioral-Normal-Subject\percentiles_SVCandSPEMv.mat');
% condList = [2,6,20];
% idxList = [2 4 6];
% for j = [1:3]
%     thisCond = condList(j);
%     thisIdx = idxList(j);
% for i = 1:size(V_posS,2)
%     eval(['maxV = prctlSPEMv_',thisSubject,num2str(thisCond),'(1,end);']);
%     eval(['minV = prctlSPEMv_',thisSubject,num2str(thisCond),'(1,1);']);
%     eval(['thisTV = TV(',num2str(thisIdx),',i);']);
%     if thisTV > minV && thisTV < maxV
%         eval(['y_spem',thisSubject,'(i,j) = pwLinearModel(TV(',num2str(thisIdx),',i),prctlSPEMv_',thisSubject,num2str(thisCond),'(1,:),prctlSPEMv_',thisSubject,num2str(thisCond),'(2,:));']);
%         eval(['error_spem_',thisSubject,'(i,j) = abs(y_spem',thisSubject,'(i,j) - V_posS(',num2str(thisIdx),',i));']);
%         eval(['yplus_spem_',thisSubject,'(i,j) = y_spem',thisSubject,'(i,j) + error_spem_',thisSubject,'(i,j)/2;']);
%         eval(['yminus_spem_',thisSubject,'(i,j) = y_spem',thisSubject,'(i,j) - error_spem_',thisSubject,'(i,j)/2;']);
%     else
%         eval(['y_spem',thisSubject,'(i,j) = nan;']);
%         eval(['error_spem_',thisSubject,'(i,j) = nan;']);
%         eval(['yplus_spem_',thisSubject,'(i,j) = nan;'])
%         eval(['yminus_spem_',thisSubject,'(i,j) = nan;'])
%     end
% end
% end
% 
% for j = [1:3]
%     thisCond = condList(j);
%     thisIdx = idxList(j);
% for i = 1:size(V_posS,2)
%     eval(['maxV = prctlSVC_',thisSubject,num2str(thisCond),'(1,end);']);
%     eval(['minV = prctlSVC_',thisSubject,num2str(thisCond),'(1,1);']);
%     eval(['thisTV = TV(',num2str(thisIdx),',i);']);
%     if thisTV > minV && thisTV < maxV
%         eval(['y_svc',thisSubject,'(i,j) = pwLinearModel(TV(',num2str(thisIdx),',i),prctlSVC_',thisSubject,num2str(thisCond),'(1,:),prctlSVC_',thisSubject,num2str(thisCond),'(2,:));']);
%         eval(['error_svc_',thisSubject,'(i,j) = abs(y_svc',thisSubject,'(i,j) - SVC(',num2str(thisIdx),',i));']);
%         eval(['yplus_svc_',thisSubject,'(i,j) = y_svc',thisSubject,'(i,j) + error_svc_',thisSubject,'(i,j)/2;']);
%         eval(['yminus_svc_',thisSubject,'(i,j) = y_svc',thisSubject,'(i,j) - error_svc_',thisSubject,'(i,j)/2;']);
%     else
%         eval(['y_svc',thisSubject,'(i,j) = nan;']);
%         eval(['error_svc_',thisSubject,'(i,j) = nan;']);
%         eval(['yplus_svc_',thisSubject,'(i,j) = nan;'])
%         eval(['yminus_svc_',thisSubject,'(i,j) = nan;'])
%     end
% end
% end
% subjectsList = {'ag','az','cs','gc','hr','ls','mp','sb','tc','vs'};
% figure;
% for i = 1:10
%     subplot(3,4,i);
%     eval(['plot(nanmedian(error_spem_',subjectsList{i},',1));']);
% end
% figure;
% for i = 1:10
%     subplot(3,4,i);
%     eval(['plot(nanmedian(error_svc_',subjectsList{i},',1));']);
% end

%% percentile 2

SVC = (squeeze(S(:,:,1))-10)*1000./(squeeze(S(:,:,3))-1000);
[TVsorted,sortedIdx]=sort([TV(1,:),TV(2,:)]);
SVCtemp = [SVC(1,:),SVC(2,:)];
SVCtemp = SVCtemp(sortedIdx);
k = 0;
for i = 1:length(SVCtemp)/5:length(SVCtemp)
    k = k + 1;
    prctlSVC_2(1,k) = nanmedian(SVCtemp(i:(i+(length(SVCtemp)/5)-1)));
    prctlSVC_2(2,k) = nanmedian(TVsorted(i:(i+(length(SVCtemp)/5)-1)));
%     
%     varSVC_2(1,k) = nanvar(SVCtemp(i:(i+(length(SVCtemp)/5)-1)));
%     varSVC_2(2,k) = nanvar(TVsorted(i:(i+(length(SVCtemp)/5)-1)));
    
end

[TVsorted,sortedIdx]=sort([TV(3,:),TV(4,:)]);
SVCtemp = [SVC(3,:),SVC(4,:)];
SVCtemp = SVCtemp(sortedIdx);
k = 0;

for i = 1:length(SVCtemp)/5:length(SVCtemp)
    k = k + 1;
    prctlSVC_6(1,k) = nanmedian(SVCtemp(i:(i+(length(SVCtemp)/5)-1)));
    prctlSVC_6(2,k) = nanmedian(TVsorted(i:(i+(length(SVCtemp)/5)-1)));
%     
%     varSVC_6(1,k) = nanvar(SVCtemp(i:(i+(length(SVCtemp)/5)-1)));
%     varSVC_6(2,k) = nanvar(TVsorted(i:(i+(length(SVCtemp)/5)-1)));
end
    
[TVsorted,sortedIdx]=sort([TV(5,:),TV(6,:)]);
SVCtemp = [SVC(5,:),SVC(6,:)];
SVCtemp = SVCtemp(sortedIdx);
k = 0;

for i = 1:length(SVCtemp)/5:length(SVCtemp)
    k = k + 1;
    prctlSVC_20(1,k) = nanmedian(SVCtemp(i:(i+(length(SVCtemp)/5)-1)));
    prctlSVC_20(2,k) = nanmedian(TVsorted(i:(i+(length(SVCtemp)/5)-1)));
%     varSVC_20(1,k) = nanvar(SVCtemp(i:(i+(length(SVCtemp)/5)-1)));
%     varSVC_20(2,k) = nanvar(TVsorted(i:(i+(length(SVCtemp)/5)-1)));
end
% bias_20 = mean(prctlSVC_20(1,:));
mdl = LinearModel.fit(prctlSVC_20(2,:),prctlSVC_20(1,:));
bias_20 = mdl.Coefficients.Estimate(1);
sensitivity_sacc_20 = mdl.Coefficients.Estimate(2);
R_sacc_20 = mdl.Rsquared.Ordinary;
MSE_sacc_20 = mdl.MSE;
SSE_sacc_20 = mdl.SSE;
SSR_sacc_20 = mdl.SSR;
SST_sacc_20 = mdl.SST;
LL_sacc_20 = mdl.LogLikelihood;
AIC_sacc_20 = mdl.ModelCriterion.AIC;
RMSE_sacc_20 = mdl.RMSE;
% bias_20 = (prctlSVC_20(1,1));
% sensitivity_20 = ((prctlSVC_20(2,:) - mean(prctlSVC_20(2,:)))./std(prctlSVC_20(2,:)))*((prctlSVC_20(1,:) - mean(prctlSVC_20(1,:)))./std(prctlSVC_20(1,:)))';
% sensitivity_20 = sensitivity_20./((norm(prctlSVC_20(2,:) - mean(prctlSVC_20(2,:)))./std(prctlSVC_20(2,:)))*norm((prctlSVC_20(1,:) - mean(prctlSVC_20(1,:)))./std(prctlSVC_20(1,:))));

mdl = LinearModel.fit(prctlSVC_6(2,:),prctlSVC_6(1,:));
bias_6 = mdl.Coefficients.Estimate(1);
sensitivity_sacc_6 = mdl.Coefficients.Estimate(2);
R_sacc_6 = mdl.Rsquared.Ordinary;
MSE_sacc_6 = mdl.MSE;
SSE_sacc_6 = mdl.SSE;
SSR_sacc_6 = mdl.SSR;
SST_sacc_6 = mdl.SST;
LL_sacc_6 = mdl.LogLikelihood;
AIC_sacc_6 = mdl.ModelCriterion.AIC;
RMSE_sacc_6 = mdl.RMSE;
% bias_6 = (prctlSVC_6(1,1));
% sensitivity_6 = ((prctlSVC_6(2,:) - mean(prctlSVC_6(2,:)))./std(prctlSVC_6(2,:)))*((prctlSVC_6(1,:) - mean(prctlSVC_6(1,:)))./std(prctlSVC_6(1,:)))';
% sensitivity_6 = sensitivity_6./((norm(prctlSVC_6(2,:) - mean(prctlSVC_6(2,:)))./std(prctlSVC_6(2,:)))*norm((prctlSVC_6(1,:) - mean(prctlSVC_6(1,:)))./std(prctlSVC_6(1,:))));
% sensitivity_6 = corrcoef(((prctlSVC_6(2,:) - mean(prctlSVC_6(2,:)))./std(prctlSVC_6(2,:))),((prctlSVC_6(1,:) - mean(prctlSVC_6(1,:)))./std(prctlSVC_6(1,:))));
% sensitivity_6 = sensitivity_6(2);
% sensitivity_6 = tan(acos(sensitivity_6));
% bias_2 = mean(prctlSVC_2(1,:));
mdl = LinearModel.fit(prctlSVC_2(2,:),prctlSVC_2(1,:));
bias_2 = mdl.Coefficients.Estimate(1);
sensitivity_sacc_2 = mdl.Coefficients.Estimate(2);
R_sacc_2 = mdl.Rsquared.Ordinary;
MSE_sacc_2 = mdl.MSE;
SSE_sacc_2 = mdl.SSE;
SSR_sacc_2 = mdl.SSR;
SST_sacc_2 = mdl.SST;
LL_sacc_2 = mdl.LogLikelihood;
AIC_sacc_2 = mdl.ModelCriterion.AIC;
RMSE_sacc_2 = mdl.RMSE;
% bias_2 = (prctlSVC_2(1,1));
% sensitivity_2 = ((prctlSVC_2(2,:) - mean(prctlSVC_2(2,:)))./std(prctlSVC_2(2,:)))*((prctlSVC_2(1,:) - mean(prctlSVC_2(1,:)))./std(prctlSVC_2(1,:)))';
% sensitivity_2 = sensitivity_2./((norm(prctlSVC_2(2,:) - mean(prctlSVC_2(2,:)))./std(prctlSVC_2(2,:)))*norm((prctlSVC_2(1,:) - mean(prctlSVC_2(1,:)))./std(prctlSVC_2(1,:))));
% sensitivity_2 = corrcoef(((prctlSVC_2(2,:) - mean(prctlSVC_2(2,:)))./std(prctlSVC_2(2,:))),((prctlSVC_2(1,:) - mean(prctlSVC_2(1,:)))./std(prctlSVC_2(1,:))));
% sensitivity_2 = sensitivity_2(2);
% sensitivity_2 = tan(acos(sensitivity_2));
figure;plot(prctlSVC_20(2,:),prctlSVC_20(1,:),'xk');hold on;...
plot(prctlSVC_2(2,:),prctlSVC_2(1,:),'or');...
plot(prctlSVC_6(2,:),prctlSVC_6(1,:),'og');
% 
% figure;plot(prctlSVC_20(1,:)./sqrt(varSVC_20(1,:)),'-ob');
% hold on;
% plot(prctlSVC_2(1,:)./sqrt(varSVC_6(1,:)),'-or')
% plot(prctlSVC_6(1,:)./sqrt(varSVC_2(1,:)),'-og')
% figure;plot([sensitivity_2,sensitivity_6,sensitivity_20],'-o')


%%
% SPEM

[TVsorted,sortedIdx]=sort([TV(2,:)]);
V_posStemp = [-V_posS(2,:)];
V_posStemp = V_posStemp(sortedIdx);
k = 0;
for i = 1:length(V_posStemp)/length(V_posStemp):length(V_posStemp)
    k = k + 1;
    prctlSPEMv_2(1,k) = nanmean(V_posStemp(i:(i+(length(V_posStemp)/length(V_posStemp))-1)));
    prctlSPEMv_2(2,k) = nanmean(TVsorted(i:(i+(length(V_posStemp)/length(V_posStemp))-1)));
end

[TVsorted,sortedIdx]=sort([TV(4,:)]);
V_posStemp = [-V_posS(4,:)];
V_posStemp = V_posStemp(sortedIdx);
k = 0;
for i = 1:length(V_posStemp)/length(V_posStemp):length(V_posStemp)
    k = k + 1;
    prctlSPEMv_6(1,k) = nanmean(V_posStemp(i:(i+(length(V_posStemp)/length(V_posStemp))-1)));
    prctlSPEMv_6(2,k) = nanmean(TVsorted(i:(i+(length(V_posStemp)/length(V_posStemp))-1)));
end

[TVsorted,sortedIdx]=sort([TV(6,:)]);
V_posStemp = [-V_posS(6,:)];
V_posStemp = V_posStemp(sortedIdx);
k = 0;
for i = 1:length(V_posStemp)/length(V_posStemp):length(V_posStemp)
    k = k + 1;
    prctlSPEMv_20(1,k) = nanmean(V_posStemp(i:(i+(length(V_posStemp)/length(V_posStemp))-1)));
    prctlSPEMv_20(2,k) = nanmean(TVsorted(i:(i+(length(V_posStemp)/length(V_posStemp))-1)));
end   
figure;plot(prctlSPEMv_20(2,:),prctlSPEMv_20(1,:),'xk');hold on;...
plot(prctlSPEMv_2(2,:),prctlSPEMv_2(1,:),'or');...
plot(prctlSPEMv_6(2,:),prctlSPEMv_6(1,:),'og');
% 


mdl = LinearModel.fit(prctlSPEMv_20(2,:),prctlSPEMv_20(1,:));
bias_20 = mdl.Coefficients.Estimate(1);
sensitivity_spem_20 = mdl.Coefficients.Estimate(2);
R_20 = mdl.Rsquared.Ordinary;
MSE_spem_20 = mdl.MSE;
SSE_spem_20 = mdl.SSE;
SSR_spem_20 = mdl.SSR;
SST_spem_20 = mdl.SST;
LL_spem_20 = mdl.LogLikelihood;
AIC_spem_20 = mdl.ModelCriterion.AIC;
RMSE_spem_20 = mdl.RMSE;


mdl = LinearModel.fit(prctlSPEMv_6(2,:),prctlSPEMv_6(1,:));
bias_6 = mdl.Coefficients.Estimate(1);
sensitivity_spem_6 = mdl.Coefficients.Estimate(2);
R_6 = mdl.Rsquared.Ordinary;
MSE_spem_6 = mdl.MSE;
SSE_spem_6 = mdl.SSE;
SSR_spem_6 = mdl.SSR;
SST_spem_6 = mdl.SST;
LL_spem_6 = mdl.LogLikelihood;
AIC_spem_6 = mdl.ModelCriterion.AIC;
RMSE_spem_6 = mdl.RMSE;

mdl = LinearModel.fit(prctlSPEMv_2(2,:),prctlSPEMv_2(1,:));
bias_2 = mdl.Coefficients.Estimate(1);
sensitivity_spem_2 = mdl.Coefficients.Estimate(2);
R_2 = mdl.Rsquared.Ordinary;
MSE_spem_2 = mdl.MSE;
SSE_spem_2 = mdl.SSE;
SSR_spem_2 = mdl.SSR;
SST_spem_2 = mdl.SST;
LL_spem_2 = mdl.LogLikelihood;
AIC_spem_2 = mdl.ModelCriterion.AIC;
RMSE_spem_2 = mdl.RMSE;


%% learning rate
% WhichPartList = {1:10,11:40,41:70,71:100,101:130,131:160};
% WhichPartList = {1:30,31:60,61:90,91:120,121:150,151:180};
% WhichPartList = {1:10,11:40,41:70,71:100,101:130,131:160,161:190};
% WhichPartList = {1:20,21:50,51:80,81:110,111:140,141:170,171:200};
% WhichPartList = {1:30,31:60,61:90,91:120,121:150,151:180,181:210};
% for w = 1:length(WhichPartList)
%     WhichPart = WhichPartList{w};
% SVC = (squeeze(S(:,:,1))-10)*1000./(squeeze(S(:,:,3))-1000);
% [TVsorted,sortedIdx]=sort([TV(5,WhichPart),TV(6,WhichPart)]);
% SVCtemp = [SVC(5,WhichPart),SVC(6,WhichPart)];
% SVCtemp = SVCtemp(sortedIdx);
% k = 0;
% 
% for i = 1:length(SVCtemp)/5:length(SVCtemp)
%     k = k + 1;
%     prctlSVC_20(1,k) = nanmedian(SVCtemp(i:(i+(length(SVCtemp)/5)-1)));
%     prctlSVC_20(2,k) = nanmedian(TVsorted(i:(i+(length(SVCtemp)/5)-1)));
% %     varSVC_20(1,k) = nanvar(SVCtemp(i:(i+(length(SVCtemp)/5)-1)));
% %     varSVC_20(2,k) = nanvar(TVsorted(i:(i+(length(SVCtemp)/5)-1)));
% end
% mdl = LinearModel.fit(prctlSVC_20(2,:),prctlSVC_20(1,:));
% bias_20(w) = mdl.Coefficients.Estimate(1);
% sensitivity_20(w) = mdl.Coefficients.Estimate(2);
% R_20(w) = mdl.Rsquared.Ordinary;
% 
% % bias_20(w) = (prctlSVC_20(1,1));
% % sensitivity_20(w) = ((prctlSVC_20(2,:) - mean(prctlSVC_20(2,:)))./std(prctlSVC_20(2,:)))*((prctlSVC_20(1,:) - mean(prctlSVC_20(1,:)))./std(prctlSVC_20(1,:)))';
% % sensitivity_20(w) = sensitivity_20(w)./((norm(prctlSVC_20(2,:) - mean(prctlSVC_20(2,:)))./std(prctlSVC_20(2,:)))*norm((prctlSVC_20(1,:) - mean(prctlSVC_20(1,:)))./std(prctlSVC_20(1,:))));
% end
