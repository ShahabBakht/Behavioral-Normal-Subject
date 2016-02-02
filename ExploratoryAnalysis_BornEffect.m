%%

X = X1901(:,:,1:1805);
X(:,31:60,:) = X1902(:,:,1:1805);
X(:,61:90,:) = X1903(:,:,1:1805);
X(:,91:120,:) = X1904(:,:,1:1805);
X(:,121:150,:) = X1905(:,:,1:1805);
X(:,151:180,:) = X2001(:,:,1:1805);
X(:,181:210,:) = X2002(:,:,1:1805);
X(:,211:240,:) = X2003(:,:,1:1805);
% X(:,151:180,:) = X1406(:,:,1:1808);
% X = X1103(:,:,1:1806);
% X(:,31:60,:) = X1104(:,:,1:1806);
% X(:,61:90,:) = X1105(:,:,1:1806);
% X(:,91:120,:) = X1106(:,:,1:1806);
% X(:,121:150,:) = X1107(:,:,1:1806);
% X(:,151:180,:) = X1201(:,:,1:1806);


S = S1901;
S(:,31:60,:) = S1902;
S(:,61:90,:) = S1903;
S(:,91:120,:) = S1904;
S(:,121:150,:) = S1905;
S(:,151:180,:) = S2001;
S(:,181:210,:) = S2002;
S(:,211:240,:) = S2003;
% S(:,151:180,:) = S1806;
% S = S1103;
% S(:,31:60,:) = S1104;
% S(:,61:90,:) = S1105;
% S(:,91:120,:) = S1106;
% S(:,121:150,:) = S1107;
% S(:,151:180,:) = S1201;


TV = TV1901;
TV(:,31:60) = TV1902;
TV(:,61:90) = TV1903;
TV(:,91:120) = TV1904;
TV(:,121:150) = TV1905;
TV(:,151:180) = TV2001;
TV(:,181:210) = TV2002;
TV(:,211:240) = TV2003;

% TV(:,151:180) = TV1806;
% TV = TV1103;
% TV(:,31:60) = TV1104;
% TV(:,61:90) = TV1105;
% TV(:,91:120) = TV1106;
% TV(:,121:150) = TV1107;
% TV(:,151:180) = TV1201;



%%

for c = 1:6
    for tr = 1:240
        x = squeeze(X(c,tr,:));
        Send = floor(S(c,tr,3));
        Sstart = floor(S(c,tr,2));
        if length(x) < Send + 120
            x_posS = [];
            v_posS = [];
            X_posS(c,tr,:) = nan;
            V_posS(c,tr,:) = nan;
        else
            x_posS = x((Send + 60):(Send + 120));
            v_posS = (x_posS(end) - x_posS(1))./0.06;
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
    end 
end


%%
x1 = squeeze(X_posS(1,:,:))';
y1 = squeeze(Y_posS(1,:,:))';
x2 = squeeze(X_posS(2,:,:))';
y2 = squeeze(Y_posS(2,:,:))';
x3 = squeeze(X_posS(3,:,:))';
y3 = squeeze(Y_posS(3,:,:))';
x4 = squeeze(X_posS(4,:,:))';
y4 = squeeze(Y_posS(4,:,:))';
x5 = squeeze(X_posS(5,:,:))';
y5 = squeeze(Y_posS(5,:,:))';
x6 = squeeze(X_posS(6,:,:))';
y6 = squeeze(Y_posS(6,:,:))';
x7 = squeeze(X_posS(7,:,:))';
y7 = squeeze(Y_posS(7,:,:))';
x8 = squeeze(X_posS(8,:,:))';
y8 = squeeze(Y_posS(8,:,:))';


%% overlap intials

% X
x = [x1(1,:),x2(1,:),-x3(1,:),-x4(1,:),x5(1,:),-x6(1,:),x7(1,:),-x8(1,:)];
xinitmin = min(x);
x = x - min(x);
x1init = x(1:80);
x2init = x(81:160);
x3init = x(161:240);
x4init = x(241:320);
x5init = x(321:400);
x6init = x(401:480);
x7init = x(481:560);
x8init = x(561:640);


x1initall = repmat(x1init,61,1);
x2initall = repmat(x2init,61,1);
x3initall = repmat(x3init,61,1);
x4initall = repmat(x4init,61,1);
x5initall = repmat(x5init,61,1);
x6initall = repmat(x6init,61,1);
x7initall = repmat(x7init,61,1);
x8initall = repmat(x8init,61,1);

x1 = x1 - x1initall;
x2 = x2 - x2initall;
x3 = x3 - x3initall;
x4 = x4 - x4initall;
x5 = x5 - x5initall;
x6 = x6 - x6initall;
x7 = x7 - x7initall;
x8 = x8 - x8initall;

% Y
y = [y1(1,:),y2(1,:),-y3(1,:),-y4(1,:),y5(1,:),-y6(1,:),y7(1,:),-y8(1,:)];
yinitmin = min(y);
y = y - min(y);
y1init = y(1:60);
y2init = y(61:120);
y3init = y(121:180);
y4init = y(181:240);
y5init = y(241:300);
y6init = y(301:360);
y7init = y(361:420);
y8init = y(421:480);


y1initall = repmat(y1init,61,1);
y2initall = repmat(y2init,61,1);
y3initall = repmat(y3init,61,1);
y4initall = repmat(y4init,61,1);
y5initall = repmat(y5init,61,1);
y6initall = repmat(y6init,61,1);
y7initall = repmat(y7init,61,1);
y8initall = repmat(y8init,61,1);

y1 = y1 - y1initall;
y2 = y2 - y2initall;
y3 = y3 - y3initall;
y4 = y4 - y4initall;
y5 = y5 - y5initall;
y6 = y6 - y6initall;
y7 = y7 - y7initall;
y8 = y8 - y8initall;




%% plot

figure;plot(x1,'-','color',[0 0 0]);hold on;plot(x2,'-','color',[1 0 0]);
plot(x5,'-','color',[0 0 1]);plot(x7,'-','color',[0 1 1]);

%% plot 2
                        
V_posS(V_posS > 25) = nan;
V_posS(V_posS < -25) = nan;
% fitobject = fit([TV(2,:),TV(4,:)]',[V_posS(2,:),-V_posS(4,:)]','smoothingspline','SmoothingParam',.2);
mdl = LinearModel.fit([TV(1,:),TV(2,:)]',[V_posS(1,:),-V_posS(2,:)]');[vfit2, ~] = predict(mdl,(10:0.5:20)');
R2 = mdl.Rsquared.Adjusted;
% vfit2 = feval(fitobject,10:0.5:20);
% fitobject = fit([TV(1,:),TV(3,:)]',[V_posS(1,:),-V_posS(3,:)]','smoothingspline','SmoothingParam',.2);
% vfit1 = feval(fitobject,10:0.5:20);
mdl = LinearModel.fit([TV(3,:),TV(4,:)]',[V_posS(3,:),-V_posS(4,:)]');[vfit1, ~] = predict(mdl,(10:0.5:20)');
R1 = mdl.Rsquared.Adjusted;
% fitobject = fit([TV(5,:),TV(6,:)]',[V_posS(5,:),-V_posS(6,:)]','smoothingspline','SmoothingParam',.2);
% vfit5 = feval(fitobject,10:0.5:20);
mdl = LinearModel.fit([TV(5,:),TV(6,:)]',[V_posS(5,:),-V_posS(6,:)]');[vfit5, ~] = predict(mdl,(10:0.5:20)');
R5 = mdl.Rsquared.Adjusted;
% fitobject = fit([TV(7,:),TV(8,:)]',[V_posS(7,:),-V_posS(8,:)]','smoothingspline','SmoothingParam',.2);
% vfit7 = feval(fitobject,10:0.5:20);
% mdl = LinearModel.fit([TV(7,:),TV(8,:)]',[V_posS(7,:),-V_posS(8,:)]');[vfit7, ~] = predict(mdl,(10:0.5:20)');
% R7 = mdl.Rsquared.Adjusted;

figure;plot(TV(1,:),V_posS(1,:),'.r','MarkerSize',10);hold on;plot(TV(2,:),-V_posS(2,:),'.r','MarkerSize',10);plot(10:0.5:20,vfit2,'-r','LineWidth',3)

plot(TV(3,:),V_posS(3,:),'.b','MarkerSize',10);hold on;plot(TV(4,:),-V_posS(4,:),'.b','MarkerSize',10);plot(10:0.5:20,vfit1,'-b','LineWidth',3)
plot(TV(5,:),V_posS(5,:),'.g','MarkerSize',10);hold on;plot(TV(6,:),-V_posS(6,:),'.g','MarkerSize',10);plot(10:0.5:20,vfit5,'-c','LineWidth',3)
% plot(TV(7,:),V_posS(7,:),'.g','MarkerSize',10);hold on;plot(TV(8,:),-V_posS(8,:),'.g','MarkerSize',10);plot(10:0.5:20,vfit7,'-g','LineWidth',3)

% pre- vs post-

% figure;plot(V_preS(2,:),V_posS(2,:),'.r','MarkerSize',10);hold on;plot(-V_preS(4,:),-V_posS(4,:),'.r','MarkerSize',10)
% 
% plot(V_preS(1,:),V_posS(1,:),'.b','MarkerSize',10);hold on;plot(-V_preS(3,:),-V_posS(3,:),'.b','MarkerSize',10);
% plot(V_preS(5,:),V_posS(5,:),'.c','MarkerSize',10);hold on;plot(-V_preS(6,:),-V_posS(6,:),'.c','MarkerSize',10);
% plot(V_preS(7,:),V_posS(7,:),'.g','MarkerSize',10);hold on;plot(-V_preS(8,:),-V_posS(8,:),'.g','MarkerSize',10);

%% prediction

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
% 
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


%%
figure;plot(TV(2,:),V_posS(2,:),'.b','MarkerSize',10);hold on
plot(TV(1,:),V_posS(1,:),'.b','MarkerSize',10);
plot(TV(5,:),V_posS(5,:),'.b','MarkerSize',10);
plot(TV(7,:),V_posS(7,:),'.b','MarkerSize',10);
plot(TV(4,:),-V_posS(4,:),'.r','MarkerSize',10)
plot(TV(3,:),-V_posS(3,:),'.r','MarkerSize',10);
plot(TV(6,:),-V_posS(6,:),'.r','MarkerSize',10);
plot(TV(8,:),-V_posS(8,:),'.r','MarkerSize',10);

%% correlation
% [R,P,RLO,RUP]=corrcoef(Ydata1,Xdata1);C1 = R(2,1);C1_lo = RLO(2,1);C1_up = RUP(2,1);
[R,P,RLO,RUP]=corrcoef(Ydata6,Xdata6);C6 = R(2,1);C6_lo = RLO(2,1);C6_up = RUP(2,1);
[R,P,RLO,RUP]=corrcoef(Ydata10,Xdata10);C10 = R(2,1);C10_lo = RLO(2,1);C10_up = RUP(2,1);
[R,P,RLO,RUP]=corrcoef(Ydata20,Xdata20);C20 = R(2,1);C20_lo = RLO(2,1);C20_up = RUP(2,1);

figure;plot([4,12,20],[C6,C10,C20],'-r');hold on;plot([4,12,20],[C6_lo,C10_lo,C20_lo],'-b');plot([4,12,20],[C6_up,C10_up,C20_up],'-b');

%% SVD


[coeff6,score6,latent6,~,explained6,mu6] = pca([Ydata6,Xdata6]);Eidx6 = explained6(1)/100;
[coeff10,score10,latent10,~,explained10,mu10] = pca([Ydata10,Xdata10]);Eidx10 = explained10(1)/100;
[coeff20,score20,latent20,~,explained20,mu20] = pca([Ydata20,Xdata20]);Eidx20 = explained20(1)/100;

idx_ci6 = bootci(20000,{@(x)EllipticalIndex(x),[Ydata6,Xdata6]},'type','norm');
idx_ci10 = bootci(20000,{@(x)EllipticalIndex(x),[Ydata10,Xdata10]},'type','norm');
idx_ci20 = bootci(20000,{@(x)EllipticalIndex(x),[Ydata20,Xdata20]},'type','norm');

figure;plot([4,12,20],[Eidx6,Eidx10,Eidx20],'b')
figure;plot([4,12,20],[idx_ci6(1),idx_ci10(1),idx_ci20(1)],'r');hold on
plot([4,12,20],[idx_ci6(2),idx_ci10(2),idx_ci20(2)],'r');


