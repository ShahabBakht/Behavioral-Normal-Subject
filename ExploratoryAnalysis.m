%% Wrap all the tests
ListOfTests = {...
%                 '15080301','15080303','15080304' ...                 
%                 '15080305','15080306','15080307','15080308','15080309','15080310'...
%                 '15080501','15080502','15080503','15080504' ...
%                 '15080505','15080506','15080507','15080508' ...
%                 '15080601','15080602','15080603','15080604','15080605','15080606'...
%                 '15080607','15080608','15080609','15080610','15080611','15080612' ...
%                 '15081701','15081702','15081703','15081704' ...% GOOOD Shahab
%                 '15081801' ...
%                 '15082104','15082105' ... % Azar
%                 '15100901','15101601','15101701','15101702','15102001','15102002','15102003','15102004','15102005','15102101','15102201' ...
%                 '15102401','15102402','15102403','15103101','15103102','15103103' ... % Azar
%                 '15102601','15102602','15102603' ... % Bennett
%                 '15102801','15102802','15102803','15110401','15110402','15110403' ... % Shahab
                '15110501','15110502','15110503' ...
                };
InitWindow = 1930:1970;            
for testscount = 1:length(ListOfTests)
    TestName        =   ListOfTests{testscount};
    run InitiateEye;
    V = I.PreProcessedEye.EyePreProcessed.Vxtrunc;
    X = I.PreProcessedEye.EyePreProcessed.Xtrunc;
    NumConditions = I.StimulusObject.S.NumConditions;
    NumTrials = I.StimulusObject.S.NumTrials;
    for condcount = 1:NumConditions
        for trcount = 1:NumTrials
            A(condcount,trcount,:) = gradient(squeeze(V(condcount,trcount,:)),0.001);
        end
    end
%     FitError = nan(NumConditions,NumTrials);
    vend = nan(NumConditions,NumTrials);
%     for condcount = 1:NumConditions
%         for trcount = 1:NumTrials
%             
%             
%             % -------------------------------------------------------------------
%             % initiation time calculation
% %             T = 1.5:0.001:1.8 - 0.001;
% %             xsample = squeeze(X(condcount,trcount,1500:1799))';
% %             xsample = xsample(1:length(T));
% %             options = optimoptions(@lsqcurvefit,'Algorithm','levenberg-marquardt','TolFun',1e-50,'Display','off');
% %             for rep = 1:100
% %                 [p,resnorm,residual,exitflag,output] = lsqcurvefit(@SmoothHinge,[4*rand,1.6,10*rand,1.7],T,xsample);
% %                 error(rep) = resnorm;
% %                 param(rep,:) = p;
% %                 
% %             end
% %             [~,i] = min(error);
% %             pp(trcount,condcount,:) = param(i,:);
% %             h = SmoothHinge(squeeze(pp(trcount,condcount,:)),T);
% %             plot(T,h,'b');hold on;plot(T,xsample,'r');title(['error = ',num2str(min(error))]);pause(0.01);hold off
%             %-------------------------------------------------------------------
% %             InitWindow = floor((pp(trcount,condcount,2)*1000)):floor((pp(trcount,condcount,2)*1000 + 100));
% 
%             
%             FitError(condcount,trcount) = min(error);
%             vend(condcount,trcount) = param(i,3);
%             
%         end
%     end
%     
    for condcount = 1:NumConditions
        for trcount = 1:NumTrials
%             StartWindow = floor(nanmedian(pp(:,condcount,2),1)*1000);
%             EndWindow = StartWindow + 100;
%             InitWindow = StartWindow:EndWindow;
%             vend(condcount,trcount) = median(A(condcount,trcount,InitWindow),3);
%             vend(condcount,trcount) = mean(V(condcount,trcount,InitWindow),3);
             vend(condcount,trcount) = abs(nanmean(X(condcount,trcount,1930:1970),3) - nanmean(X(condcount,trcount,1700:1740),3))./0.23;
        end
    end
    % calculating the initiation time
%     for condcount = 1:NumConditions
%         for trcount = 1:NumTrials
%             T = 1.5:0.001:1.9 - 0.001;
%             xsample = squeeze(X(condcount,trcount,1500:1899))';
%             xsample = xsample(1:length(T));
%             options = optimoptions(@lsqcurvefit,'Algorithm','levenberg-marquardt','TolFun',1e-50,'Display','off');
%             for rep = 1:100
%                 [p,resnorm,residual,exitflag,output] = lsqcurvefit(@Hinge,[1.6,0,randn],T,xsample);
%                 error(rep) = resnorm;
%                 param(rep,:) = p;
%             end
%             [~,i] = min(error);
%             pp(trcount,condcount,:) = param(i,:);
%             h = Hinge(squeeze(pp(trcount,condcount,:)),T);
%         end
%     end
    
    Command1 = ['vend_',TestName,' = vend;'];
    Command2 = ['V_',TestName,' = V;'];
    Command3 = ['X_',TestName,' = X;'];
%     Command4 = ['Tinit_',TestName,' = pp(:,:,2);'];
%     Command4 = ['Tinit_',TestName,' = pp(:,:,1);'];
%     Command5 = ['FitError_',TestName,' = FitError;'];
    eval(Command1);
    eval(Command2);
    eval(Command3);
%     eval(Command4);
%     eval(Command5);
    clear I X A V vend pp FitError
end

vend = [];
Tinit = [];
for testscount = 1:length(ListOfTests)
    TestName        =   ListOfTests{testscount};
    eval(['vend = [vend,','vend_',TestName,'];']);
%     eval(['Tinit = [Tinit,','Tinit_',TestName,'];']);
end

%% Remove the saccadic pursuit trials
for testscount = 1:length(ListOfTests)
TestName        =   ListOfTests{testscount};
X = eval(['X_',TestName]);
RejectTrials = nan(NumConditions,NumTrials);
for condcount = 1:NumConditions
        for trcount = 1:NumTrials
            h = figure;
            set(h,'OuterPosition',[200   100   576   512]);
            plot(squeeze(X(condcount,trcount,1600:2000)));
            title([TestName,' condition #',num2str(condcount),' trial #',num2str(trcount)]);
            
            choice = questdlg('Is this a bad trial?', ...
                'YES or NO', ...
                'YES','NO','NO');
            if strcmp(choice,'YES')
                RejectTrials(condcount,trcount) = true;
            else
                RejectTrials(condcount,trcount) = false;
            end
            close;
        end
end


Command4 = ['RejectTrials_',TestName,' = RejectTrials;'];
eval(Command4);
clear RejectTrials X choice

end
vend = [];
Tinit = [];
for testscount = 1:length(ListOfTests)
    TestName        =   ListOfTests{testscount};
    Command5 = ['vend_',TestName,'(RejectTrials_',TestName,' == 1) = nan;'];
    eval(Command5);
    eval(['vend = [vend,','vend_',TestName,'];']);
    Command6 = ['Tinit_',TestName,'(RejectTrials_',TestName,' == 1) = nan;'];
    eval(Command6);
    eval(['Tinit = [Tinit,','Tinit_',TestName,'];']);
end


%%
% plot mean values
figure;k = 0;
vname = [5,10,25];
for i = 1:3
    k = k + 1;
%     vtot = [vend(i,:);vend(i+3,:);vend(i+6,:);vend(i+9,:);vend(i+12,:);];
%     Ttot = [Tinit(i,:);Tinit(i+3,:);Tinit(i+6,:);Tinit(i+9,:);Tinit(i+12,:)];
%      vtot = [vend(i+15,:);vend(i+18,:);vend(i+21,:);vend(i+24,:);vend(i+27,:)];
%     vtot = [vend(i,:);vend(i+3,:);vend(i+6,:);vend(i+9,:);vend(i+12,:);...
%         vend(i+15,:);vend(i+18,:);vend(i+21,:);vend(i+24,:);vend(i+27,:)];
%     vtot = [vend(i,:);vend(i+2,:);vend(i+4,:);vend(i+6,:);vend(i+8,:)];
%     vtot = [[vend(i,:);vend(i+3,:);vend(i+6,:);vend(i+9,:);vend(i+12,:)],[-vend(i+15,:);-vend(i+18,:);-vend(i+21,:);-vend(i+24,:);-vend(i+27,:)]]; % long version
%     vtot = [[vend(i,:);vend(i+3,:);vend(i+6,:)],[-vend(i+9,:);-vend(i+12,:);-vend(i+15,:)]]; % short version
    vtot = [[vend(i,:);vend(i+3,:);vend(i+6,:)]];
%     vtot = [[-vend(i+9,:);-vend(i+12,:);-vend(i+15,:)]]
    
% vtot(vtot < 0) = nan;
    vtot(vtot > 10) = nan;
% %     vtot = [-vend(i+15,:);-vend(i+18,:);-vend(i+21,:);-vend(i+24,:);-vend(i+27,:)];
%     vtot = [vend(i,:);vend(i+1,:);vend(i+2,:);vend(i+3,:);vend(i+4,:);vend(i+5,:)];
% vtot = [-vend(i+6,:);-vend(i+7,:);-vend(i+8,:);-vend(i+9,:);-vend(i+10,:);-vend(i+11,:)];
%     e = iqr(vtot,2);
subplot(2,2,k); plot(vtot,'.k');
hold on; plot(nanmean(vtot,2),'-or','LineWidth',2)  
% subplot(2,2,k); plot(Ttot,'.k');
% hold on; plot(nanmedian(Ttot,2),'-or','LineWidth',2)  
% axis([0 6 10 200]);
% set(gca, ...
%   'Box'         , 'off'     ,           ...
%   'TickDir'     , 'out'     ,           ...
%   'TickLength'  , [.02 .02] ,           ...
%   'XMinorTick'  , 'on'      ,           ...
%   'YMinorTick'  , 'on'      ,           ...
%   'YGrid'       , 'on'      ,           ...
%   'XColor'      , [.3 .3 .3],           ...
%   'XTick'       , [1,2,3],            ...
%   'YTick'       , [0:50:200],            ...
%   'XTickLabel'  , {'6','10','14','18','22'},  ...
%   'YColor'      , [.3 .3 .3],           ...
%   'LineWidth'   , 1         );
% xlabel('target diameter ^o');
% ylabel('initial acceleration ^o/s^2');
% hTitle = title(sprintf('\\it{velocity = %g (degree/s)}',vname(k)));
% set( hTitle                    , ...
%     'FontSize'   , 8           ...
%     );

end

% plot measures of dispersion
figure;k=0;
for i = 1:3
    k = k + 1;
%     vtot = [vend(i,:);vend(i+3,:);vend(i+6,:);vend(i+9,:);vend(i+12,:)];
%     Ttot = [Tinit(i,:);Tinit(i+3,:);Tinit(i+6,:);Tinit(i+9,:);Tinit(i+12,:)];
% vtot = [-vend(i+15,:);-vend(i+18,:);-vend(i+21,:);-vend(i+24,:);-vend(i+27,:)];
% vtot = [vend(i,:);vend(i+3,:);vend(i+6,:);vend(i+9,:);vend(i+12,:);...
%         -vend(i+15,:);-vend(i+18,:);-vend(i+21,:);-vend(i+24,:);-vend(i+27,:)];
%     vtot = [vend(i,:);vend(i+2,:);vend(i+4,:);vend(i+6,:);vend(i+8,:)];
vtot = [[vend(i,:);vend(i+3,:);vend(i+6,:);vend(i+9,:);vend(i+12,:)],[-vend(i+15,:);-vend(i+18,:);-vend(i+21,:);-vend(i+24,:);-vend(i+27,:)]];
%     vtot = [vend(i,:);vend(i+1,:);vend(i+2,:);vend(i+3,:);vend(i+4,:);vend(i+5,:)];
% vtot = [-vend(i+6,:);-vend(i+7,:);-vend(i+8,:);-vend(i+9,:);-vend(i+10,:);-vend(i+11,:)];

subplot(2,2,k); 
% plot(iqr(vtot,2),'-ok','LineWidth',2)
plot(nanstd(vtot,[],2),'-ok','LineWidth',2)
% plot(iqr(Ttot,2),'-ok','LineWidth',2)
%  axis([0 6 20 100]);
set(gca, ...
  'Box'         , 'off'     ,           ...
  'TickDir'     , 'out'     ,           ...
  'TickLength'  , [.02 .02] ,           ...
  'XMinorTick'  , 'on'      ,           ...
  'YMinorTick'  , 'on'      ,           ...
  'YGrid'       , 'on'      ,           ...
  'XColor'      , [.3 .3 .3],           ...
  'XTick'       , [1,2,3,4,5],            ...
  'YTick'       , [0:0.004:0.1],            ...
  'XTickLabel'  , {'6','10','14','18','22'},  ...
  'YColor'      , [.3 .3 .3],           ...
  'LineWidth'   , 1         );

xlabel('target diameter ^o');
ylabel('initial acceleration variance');
hTitle = title(sprintf('\\it{velocity = %g (degree/s)}',vname(k)));
set( hTitle                    , ...
    'FontSize'   , 8           ...
    );
end

%% Plot velocity traces
% 
trials = 1:110;

% only one velocity
v1 = [vend(3,trials),-vend(18,trials)];
v2 = [vend(6,trials),-vend(21,trials)];v2(v2 > 4) = nan;
v3 = [vend(9,trials),-vend(24,trials)];
v4 = [vend(12,trials),-vend(27,trials)];
v5 = [vend(15,trials),-vend(30,trials)];

% all the velocities pooled for standard deviation measurement
% v1 = [vend(1,trials),-vend(16,trials),vend(2,trials),-vend(17,trials),vend(3,trials),-vend(18,trials)]; %v1(v1 > 140) = nan;
% v2 = [vend(4,trials),-vend(19,trials),vend(5,trials),-vend(20,trials),vend(6,trials),-vend(21,trials)]; v2(v2 > 3) = nan;
% v3 = [vend(7,trials),-vend(22,trials),vend(8,trials),-vend(23,trials),vend(9,trials),-vend(24,trials)];
% v4 = [vend(10,trials),-vend(25,trials),vend(11,trials),-vend(26,trials),vend(12,trials),-vend(27,trials)];
% v5 = [vend(13,trials),-vend(28,trials),vend(14,trials),-vend(29,trials),vend(15,trials),-vend(30,trials)];


s_ci5 = bootci(2000,{@(x)nanstd(x),v5},'type','norm');
s_ci4 = bootci(2000,{@(x)nanstd(x),v4},'type','norm');
s_ci3 = bootci(2000,{@(x)nanstd(x),v3},'type','norm');
s_ci2 = bootci(2000,{@(x)nanstd(x),v2},'type','norm');
s_ci1 = bootci(2000,{@(x)nanstd(x),v1},'type','norm');
deltas1 = abs(s_ci1(1) - s_ci1(2));
deltas2 = abs(s_ci2(1) - s_ci2(2));
deltas3 = abs(s_ci3(1) - s_ci3(2));
deltas4 = abs(s_ci4(1) - s_ci4(2));
deltas5 = abs(s_ci5(1) - s_ci5(2));
hold on;plot(trials(end)/10,deltas1,'.b');hold on;plot(trials(end)/10,deltas2,'.k');hold on;plot(trials(end)/10,deltas3,'.r');hold on;plot(trials(end)/10,deltas4,'.g');hold on;plot(trials(end)/10,deltas5,'.c')
%%
m_ci5 = bootci(2000,{@(x)nanmean(x),v5},'type','norm');
m_ci4 = bootci(2000,{@(x)nanmean(x),v4},'type','norm');
m_ci3 = bootci(2000,{@(x)nanmean(x),v3},'type','norm');
m_ci2 = bootci(2000,{@(x)nanmean(x),v2},'type','norm');
m_ci1 = bootci(2000,{@(x)nanmean(x),v1},'type','norm');

figure;plot(1:5,[nanstd(v1),nanstd(v2),nanstd(v3),nanstd(v4),nanstd(v5)],'g');
hold on
plot(ones(1,length(s_ci1(1):0.01:s_ci1(2))),s_ci1(1):0.01:s_ci1(2),'--r')
plot(2*ones(1,length(s_ci2(1):0.01:s_ci2(2))),s_ci2(1):0.01:s_ci2(2),'--r')
plot(3*ones(1,length(s_ci3(1):0.01:s_ci3(2))),s_ci3(1):0.01:s_ci3(2),'--r')
plot(4*ones(1,length(s_ci4(1):0.01:s_ci4(2))),s_ci4(1):0.01:s_ci4(2),'--r')
plot(5*ones(1,length(s_ci5(1):0.01:s_ci5(2))),s_ci5(1):0.01:s_ci5(2),'--r')

figure;plot(1:5,[nanmean(v1),nanmean(v2),nanmean(v3),nanmean(v4),nanmean(v5)],'g');
hold on
plot(ones(1,length(m_ci1(1):0.01:m_ci1(2))),m_ci1(1):0.01:m_ci1(2),'--r')
plot(2*ones(1,length(m_ci2(1):0.01:m_ci2(2))),m_ci2(1):0.01:m_ci2(2),'--r')
plot(3*ones(1,length(m_ci3(1):0.01:m_ci3(2))),m_ci3(1):0.01:m_ci3(2),'--r')
plot(4*ones(1,length(m_ci4(1):0.01:m_ci4(2))),m_ci4(1):0.01:m_ci4(2),'--r')
plot(5*ones(1,length(m_ci5(1):0.01:m_ci5(2))),m_ci5(1):0.01:m_ci5(2),'--r')
