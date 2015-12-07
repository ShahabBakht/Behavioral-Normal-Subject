%%

X = X0701(:,:,1:2414);
X(:,21:40,:) = X0702(:,:,1:2414);
X(:,41:60,:) = X0703(:,:,1:2414);
X(:,61:80,:) = X0704(:,:,1:2414);

S = S0701;
S(:,21:40,:) = S0702;
S(:,41:60,:) = S0703;
S(:,61:80,:) = S0704;

TV = TV0701;
TV(:,21:40) = TV0702;
TV(:,41:60) = TV0703;
TV(:,61:80) = TV0704;


%%

for c = 1:8
    for tr = 1:80
        x = squeeze(X(c,tr,:));
        Send = floor(S(c,tr,3));
        Sstart = floor(S(c,tr,2));
        
        x_posS = x((Send + 40):(Send + 100));
        x_preS = x((Sstart - 100):(Sstart - 20));
        
        v_posS = (x_posS(end) - x_posS(1))./0.06;
        v_preS = (x_preS(end) - x_preS(1))./0.08;
        
%         y = squeeze(Y(c,tr,:));
%         y_posS = y((Send + 40):(Send + 100));
%         y_preS = y((Sstart - 100):(Sstart - 20));
%         
        X_posS(c,tr,:) = x_posS;
        X_preS(c,tr,:) = x_preS;
        V_posS(c,tr,:) = v_posS;
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
                        
V_posS(V_posS > 30) = 0;
V_posS(V_posS < -30) = 0;
% fitobject = fit([TV(2,:),TV(4,:)]',[V_posS(2,:),-V_posS(4,:)]','smoothingspline','SmoothingParam',.2);
mdl = LinearModel.fit([TV(2,:),TV(4,:)]',[V_posS(2,:),-V_posS(4,:)]');[vfit2, ~] = predict(mdl,(10:0.5:20)');
R2 = mdl.Rsquared.Adjusted;
% vfit2 = feval(fitobject,10:0.5:20);
% fitobject = fit([TV(1,:),TV(3,:)]',[V_posS(1,:),-V_posS(3,:)]','smoothingspline','SmoothingParam',.2);
% vfit1 = feval(fitobject,10:0.5:20);
mdl = LinearModel.fit([TV(1,:),TV(3,:)]',[V_posS(1,:),-V_posS(3,:)]');[vfit1, ~] = predict(mdl,(10:0.5:20)');
R1 = mdl.Rsquared.Adjusted;
% fitobject = fit([TV(5,:),TV(6,:)]',[V_posS(5,:),-V_posS(6,:)]','smoothingspline','SmoothingParam',.2);
% vfit5 = feval(fitobject,10:0.5:20);
mdl = LinearModel.fit([TV(5,:),TV(6,:)]',[V_posS(5,:),-V_posS(6,:)]');[vfit5, ~] = predict(mdl,(10:0.5:20)');
R5 = mdl.Rsquared.Adjusted;
% fitobject = fit([TV(7,:),TV(8,:)]',[V_posS(7,:),-V_posS(8,:)]','smoothingspline','SmoothingParam',.2);
% vfit7 = feval(fitobject,10:0.5:20);
mdl = LinearModel.fit([TV(7,:),TV(8,:)]',[V_posS(7,:),-V_posS(8,:)]');[vfit7, ~] = predict(mdl,(10:0.5:20)');
R7 = mdl.Rsquared.Adjusted;

figure;plot(TV(2,:),V_posS(2,:),'.r','MarkerSize',10);hold on;plot(TV(4,:),-V_posS(4,:),'.r','MarkerSize',10);plot(10:0.5:20,vfit2,'-r','LineWidth',3)

plot(TV(1,:),V_posS(1,:),'.b','MarkerSize',10);hold on;plot(TV(3,:),-V_posS(3,:),'.b','MarkerSize',10);plot(10:0.5:20,vfit1,'-b','LineWidth',3)
plot(TV(5,:),V_posS(5,:),'.c','MarkerSize',10);hold on;plot(TV(6,:),-V_posS(6,:),'.c','MarkerSize',10);plot(10:0.5:20,vfit5,'-c','LineWidth',3)
plot(TV(7,:),V_posS(7,:),'.g','MarkerSize',10);hold on;plot(TV(8,:),-V_posS(8,:),'.g','MarkerSize',10);plot(10:0.5:20,vfit7,'-g','LineWidth',3)

% pre- vs post-

% figure;plot(V_preS(2,:),V_posS(2,:),'.r','MarkerSize',10);hold on;plot(-V_preS(4,:),-V_posS(4,:),'.r','MarkerSize',10)
% 
% plot(V_preS(1,:),V_posS(1,:),'.b','MarkerSize',10);hold on;plot(-V_preS(3,:),-V_posS(3,:),'.b','MarkerSize',10);
% plot(V_preS(5,:),V_posS(5,:),'.c','MarkerSize',10);hold on;plot(-V_preS(6,:),-V_posS(6,:),'.c','MarkerSize',10);
% plot(V_preS(7,:),V_posS(7,:),'.g','MarkerSize',10);hold on;plot(-V_preS(8,:),-V_posS(8,:),'.g','MarkerSize',10);

%% prediction

Ydata1 = [V_posS(2,:),-V_posS(4,:)]';
Xdata1 = [TV(2,:),TV(4,:)]';
% Ydata1 = [-V_posS(4,:)]';
% Xdata1 = [TV(4,:)]';
csvwrite('D:\Project Codes\Behavioral-Normal-Subject\Y1.csv',Ydata1);
csvwrite('D:\Project Codes\Behavioral-Normal-Subject\X1.csv',Xdata1);

Ydata6 = [V_posS(1,:),-V_posS(3,:)]';
Xdata6 = [TV(1,:),TV(3,:)]';
% Ydata6 = [-V_posS(3,:)]';
% Xdata6 = [TV(3,:)]';
csvwrite('D:\Project Codes\Behavioral-Normal-Subject\Y6.csv',Ydata6);
csvwrite('D:\Project Codes\Behavioral-Normal-Subject\X6.csv',Xdata6);


Ydata10 = [V_posS(5,:),-V_posS(6,:)]';
Xdata10 = [TV(5,:),TV(6,:)]';
% Ydata10 = [-V_posS(6,:)]';
% Xdata10 = [TV(6,:)]';
csvwrite('D:\Project Codes\Behavioral-Normal-Subject\Y10.csv',Ydata10);
csvwrite('D:\Project Codes\Behavioral-Normal-Subject\X10.csv',Xdata10);


Ydata20 = [V_posS(7,:),-V_posS(8,:)]';
Xdata20 = [TV(7,:),TV(8,:)]';
% Ydata20 = [-V_posS(8,:)]';
% Xdata20 = [TV(8,:)]';
csvwrite('D:\Project Codes\Behavioral-Normal-Subject\Y20.csv',Ydata20);
csvwrite('D:\Project Codes\Behavioral-Normal-Subject\X20.csv',Xdata20);


%%
figure;plot(TV(2,:),V_posS(2,:),'.b','MarkerSize',10);hold on
plot(TV(1,:),V_posS(1,:),'.b','MarkerSize',10);
plot(TV(5,:),V_posS(5,:),'.b','MarkerSize',10);
plot(TV(7,:),V_posS(7,:),'.b','MarkerSize',10);
plot(TV(4,:),-V_posS(4,:),'.r','MarkerSize',10)
plot(TV(3,:),-V_posS(3,:),'.r','MarkerSize',10);
plot(TV(6,:),-V_posS(6,:),'.r','MarkerSize',10);
plot(TV(8,:),-V_posS(8,:),'.r','MarkerSize',10);