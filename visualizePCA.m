cd('D:\Analysis\Behavioral-Normal-Subject\PCAs');
subjectsIDs = [1,2,3,4,5,6,7,8,9];
for subjects = 1:length(subjectsIDs)
    
    for cond = 1:3
        if cond == 1
            SizeDots = 2;
        elseif cond == 2
            SizeDots = 6;
        else
            SizeDots = 20;
        end
            load(['PCA',num2str(SizeDots),'_',num2str(subjectsIDs(subjects)),'.mat']);
            Command1 = ['thisExplained = PCA',num2str(SizeDots),'.explained;'];
            Command2 = ['thisLatent = PCA',num2str(SizeDots),'.latent;'];
            eval(Command1);
            eval(Command2);
            Explained(subjects,cond,:) = thisExplained;
            Latent(subjects,cond,:) = thisLatent;
            explainedRatio(subjects,cond) = thisExplained(2)./thisExplained(1);
            LatentRatio(subjects,cond) = thisLatent(2)./thisLatent(1);
            clear thisExplained thisLatent
    end
end

figure;
for subjects = 1:length(subjectsIDs)
    subplot(1,2,1);plot([2,6,20],explainedRatio(subjects,:),'-','Color',[0.7,0.7,0.7]);hold on
    subplot(1,2,2);plot([2,6,20],squeeze(Explained(subjects,:,1)),'-','Color',[0.7,0.7,0.7]);hold on
end
subplot(1,2,1);ploterr([2,6,20],mean(explainedRatio(:,:),1),[],std(explainedRatio(:,:),1)./3,'-k');hold on
subplot(1,2,1);plot([2,6,20],mean(explainedRatio(:,:),1),'.','MarkerSize',30,'Color',[0 0 0]);hold on;
ylabel('\sigma_{pc2}^2/\sigma_{pc1}^2');xlabel('sizes (degree)')
set(gca,'XTick',[2,6,20]);box off
subplot(1,2,2);ploterr([2,6,20],nanmean(squeeze(Explained(:,:,1)),1),[],nanstd(squeeze(Explained(:,:,1)),1)./3,'-k');hold on            
subplot(1,2,2);plot([2,6,20],nanmean(squeeze(Explained(:,:,1)),1),'.-','MarkerSize',30,'Color',[0 0 0]);hold on            
ylabel('\sigma_{pc1}^2');xlabel('sizes (degree)')
set(gca,'XTick',[2,6,20]);box off