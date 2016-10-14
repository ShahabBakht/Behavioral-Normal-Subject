% sub = 10;
% load('D:\Analysis\Behavioral-Normal-Subject\GoFs\Gof_spem.mat');
% load('D:\Analysis\Behavioral-Normal-Subject\GoFs\Gof_sacc.mat');
% R_spem_all = GoF_spem.R_spem_all;
% MSE_spem_all = GoF_spem.MSE_spem_all;
% SSE_spem_all = GoF_spem.SSE_spem_all;
% SSR_spem_all = GoF_spem.SSR_spem_all;
% SST_spem_all = GoF_spem.SST_spem_all;
% LL_spem_all = GoF_spem.LL_spem_all;
% AIC_spem_all = GoF_spem.AIC_spem_all;
% RMSE_spem_all = GoF_spem.RMSE_spem_all;
% R_sacc_all = GoF_sacc.R_sacc_all;
% MSE_sacc_all = GoF_sacc.MSE_sacc_all;
% SSE_sacc_all = GoF_sacc.SSE_sacc_all;
% SSR_sacc_all = GoF_sacc.SSR_sacc_all;
% SST_sacc_all = GoF_sacc.SST_sacc_all;
% LL_sacc_all = GoF_sacc.LL_sacc_all;
% AIC_sacc_all = GoF_sacc.AIC_sacc_all;
% RMSE_sacc_all = GoF_sacc.RMSE_sacc_all;
% 
% 
% R_spem_all(sub,:) = [R_2,R_6,R_20];
% MSE_spem_all(sub,:) = [MSE_spem_2,MSE_spem_6,MSE_spem_20];
% SSE_spem_all(sub,:) = [SSE_spem_2,SSE_spem_6,SSE_spem_20];
% SSR_spem_all(sub,:) = [SSR_spem_2,SSR_spem_6,SSR_spem_20];
% SST_spem_all(sub,:) = [SST_spem_2,SST_spem_6,SST_spem_20];
% LL_spem_all(sub,:) = [LL_spem_2,LL_spem_6,LL_spem_20];
% AIC_spem_all(sub,:) = [AIC_spem_2,AIC_spem_6,AIC_spem_20];
% RMSE_spem_all(sub,:) = [RMSE_spem_2,RMSE_spem_6,RMSE_spem_20];
% 
% R_sacc_all(sub,:) = [R_sacc_2,R_sacc_6,R_sacc_20];
% MSE_sacc_all(sub,:) = [MSE_sacc_2,MSE_sacc_6,MSE_sacc_20];
% SSE_sacc_all(sub,:) = [SSE_sacc_2,SSE_sacc_6,SSE_sacc_20];
% SSR_sacc_all(sub,:) = [SSR_sacc_2,SSR_sacc_6,SSR_sacc_20];
% SST_sacc_all(sub,:) = [SST_sacc_2,SST_sacc_6,SST_sacc_20];
% LL_sacc_all(sub,:) = [LL_sacc_2,LL_sacc_6,LL_sacc_20];
% AIC_sacc_all(sub,:) = [AIC_sacc_2,AIC_sacc_6,AIC_sacc_20];
% RMSE_sacc_all(sub,:) = [RMSE_sacc_2,RMSE_sacc_6,RMSE_sacc_20];
% 
% GoF_spem.R_spem_all = R_spem_all;
% GoF_spem.MSE_spem_all = MSE_spem_all;
% GoF_spem.SSE_spem_all = SSE_spem_all;
% GoF_spem.SSR_spem_all = SSR_spem_all;
% GoF_spem.SST_spem_all = SST_spem_all;
% GoF_spem.LL_spem_all = LL_spem_all;
% GoF_spem.AIC_spem_all = AIC_spem_all;
% GoF_spem.RMSE_spem_all = RMSE_spem_all;
% 
% GoF_sacc.R_sacc_all = R_sacc_all;
% GoF_sacc.MSE_sacc_all = MSE_sacc_all;
% GoF_sacc.SSE_sacc_all = SSE_sacc_all;
% GoF_sacc.SSR_sacc_all = SSR_sacc_all;
% GoF_sacc.SST_sacc_all = SST_sacc_all;
% GoF_sacc.LL_sacc_all = LL_sacc_all;
% GoF_sacc.AIC_sacc_all = AIC_sacc_all;
% GoF_sacc.RMSE_sacc_all = RMSE_sacc_all;
% 
% save('D:\Analysis\Behavioral-Normal-Subject\GoFs\GoF_spem.mat','GoF_spem');
% save('D:\Analysis\Behavioral-Normal-Subject\GoFs\GoF_sacc.mat','GoF_sacc');
% 
% 
% clear all
% close all

load('D:\Analysis\Behavioral-Normal-Subject\GoFs\Gof_spem.mat');
load('D:\Analysis\Behavioral-Normal-Subject\GoFs\Gof_sacc.mat');

R_spem_all = GoF_spem.R_spem_all;
MSE_spem_all = GoF_spem.MSE_spem_all;
SSE_spem_all = GoF_spem.SSE_spem_all;
SSR_spem_all = GoF_spem.SSR_spem_all;
SST_spem_all = GoF_spem.SST_spem_all;
LL_spem_all = GoF_spem.LL_spem_all;
AIC_spem_all = GoF_spem.AIC_spem_all;
RMSE_spem_all = GoF_spem.RMSE_spem_all;

R_sacc_all = GoF_sacc.R_sacc_all;
MSE_sacc_all = GoF_sacc.MSE_sacc_all;
SSE_sacc_all = GoF_sacc.SSE_sacc_all;
SSR_sacc_all = GoF_sacc.SSR_sacc_all;
SST_sacc_all = GoF_sacc.SST_sacc_all;
LL_sacc_all = GoF_sacc.LL_sacc_all;
AIC_sacc_all = GoF_sacc.AIC_sacc_all;
RMSE_sacc_all = GoF_sacc.RMSE_sacc_all;

figure;subplot(3,3,1);plot(R_spem_all','o-','Color',[0.7,0.7,0.7]);hold on;plot(nanmean(R_spem_all,1),'o-k');title('R^2')
subplot(3,3,2);plot(MSE_spem_all','o-','Color',[0.7,0.7,0.7]);hold on;plot(nanmean(MSE_spem_all,1),'o-k');title('MSE')
subplot(3,3,3);plot(SSE_spem_all','o-','Color',[0.7,0.7,0.7]);hold on;plot(nanmean(SSE_spem_all,1),'o-k');title('SSE')
subplot(3,3,4);plot(SSR_spem_all','o-','Color',[0.7,0.7,0.7]);hold on;plot(nanmean(SSR_spem_all,1),'o-k');title('SSR')
subplot(3,3,5);plot(SST_spem_all','o-','Color',[0.7,0.7,0.7]);hold on;plot(nanmean(SST_spem_all,1),'o-k');title('SST')
subplot(3,3,6);plot(LL_spem_all','o-','Color',[0.7,0.7,0.7]);hold on;plot(nanmean(LL_spem_all,1),'o-k');title('LL')
subplot(3,3,7);plot(AIC_spem_all','o-','Color',[0.7,0.7,0.7]);hold on;plot(nanmean(AIC_spem_all,1),'o-k');title('AIC')
subplot(3,3,8);plot(RMSE_spem_all','o-','Color',[0.7,0.7,0.7]);hold on;plot(nanmean(RMSE_spem_all,1),'o-k');title('RMSE')

figure;subplot(3,3,1);plot(R_sacc_all','o-','Color',[0.7,0.7,0.7]);hold on;plot(nanmean(R_sacc_all,1),'o-k');title('R^2')
subplot(3,3,2);plot(MSE_sacc_all','o-','Color',[0.7,0.7,0.7]);hold on;plot(nanmean(MSE_sacc_all,1),'o-k');title('MSE')
subplot(3,3,3);plot(SSE_sacc_all','o-','Color',[0.7,0.7,0.7]);hold on;plot(nanmean(SSE_sacc_all,1),'o-k');title('SSE')
subplot(3,3,4);plot(SSR_sacc_all','o-','Color',[0.7,0.7,0.7]);hold on;plot(nanmean(SSR_sacc_all,1),'o-k');title('SSR')
subplot(3,3,5);plot(SST_sacc_all','o-','Color',[0.7,0.7,0.7]);hold on;plot(nanmean(SST_sacc_all,1),'o-k');title('SST')
subplot(3,3,6);plot(LL_sacc_all','o-','Color',[0.7,0.7,0.7]);hold on;plot(nanmean(LL_sacc_all,1),'o-k');title('LL')
subplot(3,3,7);plot(AIC_sacc_all','o-','Color',[0.7,0.7,0.7]);hold on;plot(nanmean(AIC_sacc_all,1),'o-k');title('AIC')
subplot(3,3,8);plot(RMSE_sacc_all','o-','Color',[0.7,0.7,0.7]);hold on;plot(nanmean(RMSE_sacc_all,1),'o-k');title('RMSE')
