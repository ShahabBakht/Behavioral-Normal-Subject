function Results = decodeDirMTdata()

% load Dave's MT data - it contains a cell array named 'shahab'
[FileName,PathName] = uigetfile('*.mat','Select Dave''s MT data');
load([PathName,'/',FileName]);
Results.data = shahab;
numDays = length(shahab);

%% calculate suppresssion index SI for each neuron
allSI_withBaseline = [];
allSI_withoutBaseline = [];

for daycount = 1:numDays
%     fprintf(['day ',num2str(daycount), ' --> \n']);
    todayData = shahab{daycount};
    numNeurons = (size(todayData,2) - 3)/2; % 3 last columns are direction and size of the stimulus
    stimulusSizes = unique(todayData(:,end-1));
    stimulusSizesIndex = unique(todayData(:,end));
    Datapref = todayData(todayData(:,end-2) == 1,:);
    Datanull = todayData(todayData(:,end-2) == 0,:);
    for sizecount = 1:length(stimulusSizesIndex)
        DatapthisSize = Datapref(Datapref(:,end) == stimulusSizesIndex(sizecount),:);
        DatanthisSize = Datanull(Datanull(:,end) == stimulusSizesIndex(sizecount),:);
        RpthisSize = DatapthisSize(:,1:2:(2*numNeurons-1));
        RnthisSize = DatanthisSize(:,1:2:(2*numNeurons-1));
        BasepthisSize = DatapthisSize(:,2:2:(2*numNeurons));
        BasenthisSize = DatanthisSize(:,2:2:(2*numNeurons));
        
        % Dprime with baseline 
        RpthisSizeMean = nanmedian(RpthisSize + BasepthisSize,1);
        RnthisSizeMean = nanmedian(RnthisSize + BasenthisSize,1);
        RpthisSizeVar = var(RpthisSize + BasepthisSize,[], 1);
        RnthisSizeVar = var(RnthisSize + BasenthisSize, [], 1);
        thisDprime = (RpthisSizeMean - RnthisSizeMean)./sqrt((RpthisSizeVar+RnthisSizeVar)/2);
        Dprime_withBaseline(sizecount,:) = thisDprime;
        
        %Dprime without baseline
%         RpthisSizeMean = nanmedian(RpthisSize,1);
%         RnthisSizeMean = nanmedian(RnthisSize,1);
%         RpthisSizeVar = var(RpthisSize,[], 1);
%         RnthisSizeVar = var(RnthisSize, [], 1);
%         thisDprime = (RpthisSizeMean - RnthisSizeMean)./sqrt((RpthisSizeVar+RnthisSizeVar)/2);
%         Dprime_withoutBaseline(sizecount,:) = thisDprime;
        
        clear thisDprime;
        
        
    end
%   
    % fit difference of Gaussians to the Dprime curves for each neuron
    si_wb = nan(1,numNeurons);
    si_woutb = nan(1,numNeurons);
    for neuroncount = 1:numNeurons
%         fprintf([' neuorn ',num2str(neuroncount),'\n']);
        try
            %flip the size tuning if null dir is preferred
%             figure(1);subplot(3,3,neuroncount);...
            if sum(Dprime_withBaseline(:,neuroncount)<0) > sum(Dprime_withBaseline(:,neuroncount)>=0)
                neuronSigns{daycount}(neuroncount) = -1;
                Dprime_withBaseline(:,neuroncount) = -Dprime_withBaseline(:,neuroncount);
                [a,si_wb(neuroncount),~,resnorm]=mtSizeFit(stimulusSizes,Dprime_withBaseline(:,neuroncount),1);
%                 plot(stimulusSizes,Dprime_withBaseline(:,neuroncount),'or');hold on;plot(stimulusSizes,diffGauss(a,stimulusSizes),'-k');...
            else
                neuronSigns{daycount}(neuroncount) = 1;
                [a,si_wb(neuroncount),~,resnorm]=mtSizeFit(stimulusSizes,Dprime_withBaseline(:,neuroncount),1);
%                 plot(stimulusSizes,Dprime_withBaseline(:,neuroncount),'or');hold on;plot(stimulusSizes,diffGauss(a,stimulusSizes),'-k');...
            end
        
%             thisError = resnorm/sum(Dprime_withBaseline(:,neuroncount).^2);
%             fprintf([' fit residual norm = ',num2str(thisError),'\n'])
            
            
            
            
            
            
%             if max(Dprime_withoutBaseline(:,neuroncount)) > 0
%             [~,si_woutb(neuroncount),~,~]=mtSizeFit(stimulusSizes,Dprime_withoutBaseline(:,neuroncount),1);
%             end
            
        catch
            fprintf(['error in curve fitting'])
        end
    
    end
%     answer = inputdlg('Enter bad neuorns:');close;
%     si_wb(str2num(answer{:})) = nan;
%     SI_withBaseline{daycount} = si_wb;
%     SI_withoutBaseline{daycount} = si_woutb;

    SI_withBaseline{daycount} = (max(Dprime_withBaseline,[],1) - Dprime_withBaseline(end,:))./max(Dprime_withBaseline,[],1);
    if size(Dprime_withBaseline,1) == 8
        Dprime_withBaseline = [Dprime_withBaseline;nan(1,size(Dprime_withBaseline,2))];
    end
    SizeTuning(:,daycount) = nanmean(Dprime_withBaseline,2);
%     SI_withoutBaseline{daycount} = (max(Dprime_withoutBaseline,[],1) - Dprime_withoutBaseline(end,:))./max(Dprime_withoutBaseline,[],1);
    
    clear Dprime_withBaseline %Dprime_withoutBaseline 
    
    allSI_withBaseline = [allSI_withBaseline,SI_withBaseline{daycount}];%(SI_withBaseline{daycount}>0)];
%     allSI_withoutBaseline = [allSI_withoutBaseline,SI_withoutBaseline{daycount}(SI_withoutBaseline{daycount}>0)];
    
end

% separate size tunings to SS and nSS cells

Results.SizeTuning.SizeTunings = SizeTuning;
Results.neuronSign = neuronSigns;
Results.SI.withBaseline = SI_withBaseline;
% Results.SI.withoutBaseline = SI_withoutBaseline;
Results.SI.allwithB = allSI_withBaseline;
% Results.SI.allwithoutB = allSI_withoutBaseline;

%% Group the recording sessions to days and stimulus sizes
% Rsizesorted = response of the neuorns to different stimulus sizes in each
% trial for each day
% Bsizesorted = baseline of the neuorns to different stimulus sizes in each
% trial for each day
% StimDirsizesorted = the direction of the stimulus (0 or 1) for different
% stimulus sizes in each trial for each day

% Labels = cell(size(Results.data));
for daycount = 1:numDays

    todayNeuronSigns = neuronSigns{daycount};
    todayData = shahab{daycount};
    stimulusSizesIndex = unique(todayData(:,end));
    numNeurons = (size(todayData,2) - 3)/2;

    for sizecount = 1:length(stimulusSizesIndex)
        Rsizesorted{daycount,sizecount} = todayData(todayData(:,end) == stimulusSizesIndex(sizecount),1:2:(2*numNeurons-1));
        Bsizesorted{daycount,sizecount} = todayData(todayData(:,end) == stimulusSizesIndex(sizecount),2:2:(2*numNeurons));
        StimDirsizesorted{daycount,sizecount} = todayData(todayData(:,end) == stimulusSizesIndex(sizecount),end-2);
    end
    
    
    
    
    
end
Results.Rsizesorted = Rsizesorted;
Results.Bsizesorted = Bsizesorted;
Results.StimDirsizesorted = StimDirsizesorted;

%% 10-fold x-validation classification accuracy

for daycount = 1:numDays
    for sizecount = 1:size(Rsizesorted,2)
        testX = Results.Rsizesorted{daycount,sizecount} + Results.Bsizesorted{daycount,sizecount};
        testY = Results.StimDirsizesorted{daycount,sizecount};
        if ~isempty(testX)
        lda = fitcdiscr(testX,testY);
        testY(testY == 0) = 2;
        
        cp = cvpartition(testY,'KFold',10);
        cvlda = crossval(lda,'CVPartition',cp);
%         ldaCVErr(daycount,sizecount) = kfoldLoss(cvlda);
        ldaCVErr(daycount,sizecount) = kfoldLoss(cvlda,'lossfun',@decoderVar);
%         ldaCVErr(daycount,sizecount) = kfoldLoss(cvlda,'lossfun','binodeviance');
        end
    end

end

Results.xvalidation.xvalidation = ldaCVErr;

%% Compare x-validation error with average SI across days

for daycount = 1:numDays
    xvalidationToday = ldaCVErr(daycount,:);
    
%     if xvalidationToday(end) ~= 0
%         xvalidationToday = xvalidationToday(2:end);
%     else
%         xvalidationToday = xvalidationToday(1:end-1);
%     end
    if xvalidationToday(end) == 0
        xvalidationToday = [xvalidationToday(1),xvalidationToday(1:end-1)];
    
    end
    
    Results.xvalidation.xvalidation2(daycount,:) = xvalidationToday;
    ErrorSizeEffect = (min(xvalidationToday) - xvalidationToday(end))./min(xvalidationToday);
    
    SItodaymean = nanmedian(Results.SI.withBaseline{daycount});
    SItoday = (Results.SI.withBaseline{daycount});
%     allSImedian = nanmedian(Results.SI.allwithB);
    allSImedian = prctile(Results.SI.allwithB,25);
    numSStoday = sum(SItoday>=allSImedian);
    numnSStoday = sum(SItoday<allSImedian);
    
    ErrorSI(daycount,:) = [SItodaymean,ErrorSizeEffect,numSStoday,numnSStoday];
    
end


weightedXvalidationAverage_ss = ((1./nansum(repmat((ErrorSI(:,3)./(ErrorSI(:,3)+ErrorSI(:,4))),1,size(Results.xvalidation.xvalidation2,2))))) .* ...
    nansum(repmat((ErrorSI(:,3)./(ErrorSI(:,3)+ErrorSI(:,4))),1,size(Results.xvalidation.xvalidation2,2)) .* Results.xvalidation.xvalidation2);
weightedXvalidationAverage_nss = ((1./nansum(repmat((ErrorSI(:,4)./(ErrorSI(:,3)+ErrorSI(:,4))),1,size(Results.xvalidation.xvalidation2,2))))) .* ...
    nansum(repmat((ErrorSI(:,4)./(ErrorSI(:,3)+ErrorSI(:,4))),1,size(Results.xvalidation.xvalidation2,2)) .* Results.xvalidation.xvalidation2);


% weightedXvalidationAverage_ss = nanmean(repmat((ErrorSI(:,3)./(ErrorSI(:,3)+ErrorSI(:,4))),1,size(Results.xvalidation.xvalidation2,2)) .* Results.xvalidation.xvalidation2);
% weightedXvalidationAverage_nss = nanmean(repmat((ErrorSI(:,4)./(ErrorSI(:,3)+ErrorSI(:,4))),1,size(Results.xvalidation.xvalidation2,2)) .* Results.xvalidation.xvalidation2);

% weightedSizeTuningAverage_ss = nansum(repmat((ErrorSI(:,3)./(ErrorSI(:,3)+ErrorSI(:,4))),1,size(Results.SizeTuning.SizeTunings,1)) .* Results.SizeTuning.SizeTunings');
% weightedSizeTuningAverage_nss = nansum(repmat((ErrorSI(:,4)./(ErrorSI(:,3)+ErrorSI(:,4))),1,size(Results.SizeTuning.SizeTunings,1)) .* Results.SizeTuning.SizeTunings');


[a_ss,~,~,~]=mtSizeFit(stimulusSizes',weightedXvalidationAverage_ss,1);
xvalidationfit_ss = diffGauss(a_ss,stimulusSizes);

[a_nss,~,~,~]=mtSizeFit(stimulusSizes',weightedXvalidationAverage_nss,1);
xvalidationfit_nss = diffGauss(a_nss,stimulusSizes);

Results.ErrorSI = ErrorSI;
Results.xvalidation.weightedxvalidation_ss = weightedXvalidationAverage_ss;
Results.xvalidation.weightedxvalidation_nss = weightedXvalidationAverage_nss;
% Results.SizeTuning.weightedSizeTuningAverage_ss = weightedSizeTuningAverage_ss;
% Results.SizeTuning.weightedSizeTuningAverage_nss = weightedSizeTuningAverage_nss;
Results.xvalidation.weightedxvalidation_nssfit = xvalidationfit_nss;
Results.xvalidation.weightedxvalidation_ssfit = xvalidationfit_ss;


end

function VAR = decoderVar(Y,Yfit,W,C)
Yfitbin = Yfit(:,1) > Yfit(:,2);
% VAR = sum((Yfit(:,1) - 0.5).^2);

Ytargetbin = Y(:,1);
[~,cm,~,~] = confusion(Y(:,1)',Yfit(:,1)');
Sens = cm(2,2)/(cm(2,2) + cm(2,1));
Spec = cm(1,1)/(cm(1,1) + cm(1,2));
VAR = Sens + Spec - 1;
% VAR = .5 * ((cm(2,2)/(cm(2,/2)+cm(2,1))) + (cm(2,2)/(cm(2,2)+cm(1,2))));
% YshouldbeOne = Yfitbin(Ytargetbin == 1);
% p = binofit(sum(YshouldbeOne),length(YshouldbeOne));
% p = sum(YshouldbeOne)/length(YshouldbeOne);
% VAR = length(YshouldbeOne) * p * (1 - p);
% VAR = var(Yfit(:,1));
% VAR = median(Yfit(Y(:,1) == 1,1)./(Yfit(Y(:,1) == 1,1) + Yfit(Y(:,1) == 1,2)));

end