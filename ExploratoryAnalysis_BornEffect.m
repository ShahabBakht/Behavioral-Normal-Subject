%%

X = X2602;
Y = Y2602;
X(:,21:40,:) = X2603;
Y(:,21:40,:) = Y2603;
X(:,41:60,:) = X2604(:,:,1:2914);
Y(:,41:60,:) = Y2604(:,:,1:2914);

S = S2602;
S(:,21:40,:) = S2603;
S(:,41:60,:) = S2604;

%%

for c = 1:6
    for tr = 1:60
        x = squeeze(X(c,tr,:));
        Send = floor(S(c,tr,3));
        Sstart = floor(S(c,tr,2));
        
        x_posS = x((Send + 40):(Send + 100));
        x_preS = x((Sstart - 100):(Sstart - 20));
        
        y = squeeze(Y(c,tr,:));
        y_posS = y((Send + 40):(Send + 100));
        y_preS = y((Sstart - 100):(Sstart - 20));
        
        X_posS(c,tr,:) = x_posS;
        X_preS(c,tr,:) = x_preS;
        
        Y_posS(c,tr,:) = y_posS;
        Y_preS(c,tr,:) = y_preS;
        
        clear x x_posS x_preS y y_posS y_preS Send Sstart
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

%% overlap intials

x = [x3(1,:),x4(1,:)];
xinitmin = min(x);
x = x - min(x);
x3init = x(1:60);
x4init = x(61:120);
x3initall = repmat(x3init,61,1);
x4initall = repmat(x4init,61,1);
x3 = x3 - x3initall;
x4 = x4 - x4initall;


y = [y3(1,:),y4(1,:)];
yinitmin = min(y);
y = y - min(y);
y3init = y(1:60);
y4init = y(61:120);
y3initall = repmat(y3init,61,1);
y4initall = repmat(y4init,61,1);
y3 = y3 - y3initall;
y4 = y4 - y4initall;

%% plot

figure;plot(x1,'b');hold on;plot(x2,'r');
figure;plot(x3,y3,'b');hold on;plot(x4,y4,'r');
figure;plot(x1,'b');hold on;plot(x2,'r');

%%

figure( 'Position', [ 100, 380, 1380, 550 ] )
for i = 1:60
u(:,i) = TVRegDiff( x1(:,i),500, 1e-6, [], 'small', 1e-8, .001, 1, 1 );
end
