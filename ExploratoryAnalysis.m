%% Wrap all the tests
ListOfTests = {...
%                 '15080301','15080303','15080304' ...                 
%                 '15080305','15080306','15080307','15080308','15080309','15080310'...
%                 '15080501','15080502','15080503','15080504' ...
%                 '15080505','15080506','15080507','15080508' ...
%                 '15080601','15080602','15080603','15080604','15080605','15080606'...
%                 '15080607','15080608','15080609','15080610','15080611','15080612' ...
                '15081701','15081702','15081703','15081704' ...% Shahab
%                 '15081801' ...
%                 '15082104','15082105' ... % Azar
                };
InitWindow = 1600:1700;            
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
    
    for condcount = 1:NumConditions
        for trcount = 1:NumTrials
            
            
            % -------------------------------------------------------------------
            % initiation time calculation
            T = 1.5:0.001:1.9 - 0.001;
            xsample = squeeze(X(condcount,trcount,1500:1899))';
            xsample = xsample(1:length(T));
            options = optimoptions(@lsqcurvefit,'Algorithm','levenberg-marquardt','TolFun',1e-50,'Display','off');
            for rep = 1:100
                [p,resnorm,residual,exitflag,output] = lsqcurvefit(@Hinge,[1.6,0,randn],T,xsample);
                error(rep) = resnorm;
                param(rep,:) = p;
            end
            [~,i] = min(error);
            pp(trcount,condcount,:) = param(i,:);
            h = Hinge(squeeze(pp(trcount,condcount,:)),T);
            % -------------------------------------------------------------------
            InitWindow = floor((pp(trcount,condcount,1)*1000 - 100)):floor((pp(trcount,condcount,1)*1000));
            vend(condcount,trcount) = median(A(condcount,trcount,InitWindow),3);
            
        end
    end
    
%     % calculating the initiation time
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
    Command4 = ['Tinit_',TestName,' = pp(:,:,1);'];
    eval(Command1);
    eval(Command2);
    eval(Command3);
    eval(Command4);
    clear I X V A vend pp
end

vend = [];
Tinit = [];
for testscount = 1:length(ListOfTests)
    TestName        =   ListOfTests{testscount};
    eval(['vend = [vend,','vend_',TestName,'];']);
    eval(['Tinit = [Tinit,','Tinit_',TestName,'];']);
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
vname = [15,20,25];
for i = 1:3
    k = k + 1;
    vtot = [vend(i,:);vend(i+3,:);vend(i+6,:);vend(i+9,:);vend(i+12,:)];
    Ttot = [Tinit(i,:);Tinit(i+3,:);Tinit(i+6,:);Tinit(i+9,:);Tinit(i+12,:)];
%      vtot = [vend(i+15,:);vend(i+18,:);vend(i+21,:);vend(i+24,:);vend(i+27,:)];
%     vtot = [vend(i,:);vend(i+3,:);vend(i+6,:);vend(i+9,:);vend(i+12,:);...
%         vend(i+15,:);vend(i+18,:);vend(i+21,:);vend(i+24,:);vend(i+27,:)];
%     vtot = [vend(i,:);vend(i+2,:);vend(i+4,:);vend(i+6,:);vend(i+8,:)];
%     vtot = [[vend(i,:);vend(i+3,:);vend(i+6,:);vend(i+9,:);vend(i+12,:)],[-vend(i+15,:);-vend(i+18,:);-vend(i+21,:);-vend(i+24,:);-vend(i+27,:)]];

%     vtot = [-vend(i+15,:);-vend(i+18,:);-vend(i+21,:);-vend(i+24,:);-vend(i+27,:)];
%     vtot = [vend(i,:);vend(i+1,:);vend(i+2,:);vend(i+3,:);vend(i+4,:);vend(i+5,:)];
% vtot = [-vend(i+6,:);-vend(i+7,:);-vend(i+8,:);-vend(i+9,:);-vend(i+10,:);-vend(i+11,:)];
%     e = iqr(vtot,2);
% subplot(2,2,k); plot(vtot,'.k');
% hold on; plot(nanmedian(vtot,2),'-or','LineWidth',2)  
subplot(2,2,k); %plot(Ttot,'.k');
hold on; plot(nanmedian(Ttot,2),'-or','LineWidth',2)  
% axis([0 6 10 150]);
set(gca, ...
  'Box'         , 'off'     ,           ...
  'TickDir'     , 'out'     ,           ...
  'TickLength'  , [.02 .02] ,           ...
  'XMinorTick'  , 'on'      ,           ...
  'YMinorTick'  , 'on'      ,           ...
  'YGrid'       , 'on'      ,           ...
  'XColor'      , [.3 .3 .3],           ...
  'XTick'       , [1,2,3,4,5,6,7,8,9,10],            ...
  'YTick'       , [100:20:150],            ...
  'XTickLabel'  , {'6','10','14','18','22'},  ...
  'YColor'      , [.3 .3 .3],           ...
  'LineWidth'   , 1         );
xlabel('target diameter ^o');
ylabel('initial acceleration ^o/s^2');
hTitle = title(sprintf('\\it{velocity = %g (degree/s)}',vname(k)));
set( hTitle                    , ...
    'FontSize'   , 8           ...
    );

end

% plot measures of dispersion
figure;k=0;
for i = 1:3
    k = k + 1;
    vtot = [vend(i,:);vend(i+3,:);vend(i+6,:);vend(i+9,:);vend(i+12,:)];
    Ttot = [Tinit(i,:);Tinit(i+3,:);Tinit(i+6,:);Tinit(i+9,:);Tinit(i+12,:)];
% vtot = [vend(i+15,:);vend(i+18,:);vend(i+21,:);vend(i+24,:);vend(i+27,:)];
% vtot = [vend(i,:);vend(i+3,:);vend(i+6,:);vend(i+9,:);vend(i+12,:);...
%         vend(i+15,:);vend(i+18,:);vend(i+21,:);vend(i+24,:);vend(i+27,:)];
%     vtot = [vend(i,:);vend(i+2,:);vend(i+4,:);vend(i+6,:);vend(i+8,:)];
% vtot = [[vend(i,:);vend(i+3,:);vend(i+6,:);vend(i+9,:);vend(i+12,:)],[-vend(i+15,:);-vend(i+18,:);-vend(i+21,:);-vend(i+24,:);-vend(i+27,:)]];
%     vtot = [vend(i,:);vend(i+1,:);vend(i+2,:);vend(i+3,:);vend(i+4,:);vend(i+5,:)];
% vtot = [-vend(i+6,:);-vend(i+7,:);-vend(i+8,:);-vend(i+9,:);-vend(i+10,:);-vend(i+11,:)];

subplot(2,2,k); 
% plot(iqr(vtot,2),'-ok','LineWidth',2)
% plot(nanstd(vtot,[],2),'-ok','LineWidth',2)
plot(iqr(Ttot,2),'-ok','LineWidth',2)
%  axis([0 6 20 100]);
set(gca, ...
  'Box'         , 'off'     ,           ...
  'TickDir'     , 'out'     ,           ...
  'TickLength'  , [.02 .02] ,           ...
  'XMinorTick'  , 'on'      ,           ...
  'YMinorTick'  , 'on'      ,           ...
  'YGrid'       , 'on'      ,           ...
  'XColor'      , [.3 .3 .3],           ...
  'XTick'       , [1,2,3,4,5,6,7,8,9,10],            ...
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
for condcount = 1:15
    for trcount = 1:5
        if RejectTrials_15080306(condcount,trcount) == 0
            
            plot(squeeze(V_15080306(condcount,trcount,:))','r');hold on
        end
    end
end


plot(squeeze(X_15080306(15,:,:))','c');hold on
plot(squeeze(X_15080307(15,:,:))','c');hold on
plot(squeeze(X_15080308(15,:,:))','c');hold on
plot(squeeze(X_15080309(15,:,:))','c');hold on
plot(squeeze(X_15080310(15,:,:))','c');hold on
