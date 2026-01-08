%子问题1:社会成本最小化问题
clc
clear
close all

%% ADMM迭代参数设置
%拉格朗日乘子初始化
lambda_e_12=zeros(1,24);%1和2之间的拉格朗日乘子
lambda_e_13=zeros(1,24);%1和3之间的拉格朗日乘子
lambda_e_21=zeros(1,24);%2和1之间的拉格朗日乘子
lambda_e_23=zeros(1,24);%2和3之间的拉格朗日乘子
lambda_e_31=zeros(1,24);%3和1之间的拉格朗日乘子
lambda_e_32=zeros(1,24);%3和2之间的拉格朗日乘子
%% 辅助变量，初始值无所谓
P12d=200*ones(1,24);
P13d=200*ones(1,24);
P21d=200*ones(1,24);
P23d=200*ones(1,24);
P31d=200*ones(1,24);
P32d=200*ones(1,24);
maxIter=100;  %最大迭代次数
p=0.06;
iter=1;%迭代次数初始化
Ben_Store=[];%历史目标函数
toler1=[]; %残差1,电能交易部分
%% 迭代
for i = 1:maxIter  %maxloopnumber is the lower bound of loop
%      [P_e_12,P_e_13,Obj_MG1(iter)]=Fun_MG1(P12d,P13d,lambda_e_12,lambda_e_13);
%         [P_e_21,P_e_23,Obj_MG2(iter)]=Fun_MG2(P21d,P23d,lambda_e_21,lambda_e_23);
%         [P_e_31,P_e_32,Obj_MG3(iter)]=Fun_MG3(P31d,P32d,lambda_e_31,lambda_e_32);
%         P12d=(lambda_e_21-lambda_e_12+(P_e_12-P_e_21))/2;
%         P13d=(lambda_e_31-lambda_e_13+(P_e_13-P_e_31))/2;
%         P21d=(lambda_e_12-lambda_e_21+(P_e_21-P_e_12))/2;
%         P23d=(lambda_e_32-lambda_e_23+(P_e_23-P_e_32))/2;
%         P31d=(lambda_e_13-lambda_e_31+(P_e_31-P_e_13))/2;
%         P32d=(lambda_e_23-lambda_e_32+(P_e_32-P_e_23))/2;
        [P_e_12,P_e_13,P_e_cut1,Obj_MG1(iter),P_loss1,total_loss1,CVaR1]=Fun_MG1(P12d,P13d,lambda_e_12,lambda_e_13,p);
        [P_e_21,P_e_23,P_e_cut2,Obj_MG2(iter),P_loss2,total_loss2,CVaR2]=Fun_MG2(P21d,P23d,lambda_e_21,lambda_e_23,p);
        [P_e_31,P_e_32,P_e_cut3,Obj_MG3(iter),P_loss3,total_loss3,CVaR3]=Fun_MG3(P31d,P32d,lambda_e_31,lambda_e_32,p);
        P12d=(lambda_e_21-lambda_e_12+p*(P_e_12-P_e_21))/(2*p);
        P13d=(lambda_e_31-lambda_e_13+p*(P_e_13-P_e_31))/(2*p);
        P21d=(lambda_e_12-lambda_e_21+p*(P_e_21-P_e_12))/(2*p);
        P23d=(lambda_e_32-lambda_e_23+p*(P_e_23-P_e_32))/(2*p);
        P31d=(lambda_e_13-lambda_e_31+p*(P_e_31-P_e_13))/(2*p);
        P32d=(lambda_e_23-lambda_e_32+p*(P_e_32-P_e_23))/(2*p);
        
    Ben_Store=[Ben_Store,[Obj_MG1(iter);Obj_MG2(iter);Obj_MG3(iter)]];
    %残差计算
    AP2P = [P_e_12;P_e_13;P_e_21;P_e_23;P_e_31;P_e_32];
    P2P = [P12d;P13d;P21d;P23d;P31d;P32d];
%     toler1=[toler1,norm(P_e_12-P_e_12(iter,:))+norm(P_e_13-P_e_13(iter,:))+norm(P_e_23-P_e_23(iter,:))];%保存残差1 
    %判断收敛条件
    if compareStrategies1(AP2P,P2P)==0
       display(['迭代收敛,在第 ', num2str(iter),' 次收敛']);
       break; 
    else
        lambda_e_12=lambda_e_12+p*(P12d-P_e_12);
        lambda_e_13=lambda_e_13+p*(P13d-P_e_13);
        lambda_e_21=lambda_e_21+p*(P21d-P_e_21);
        lambda_e_23=lambda_e_23+p*(P23d-P_e_23);
        lambda_e_31=lambda_e_31+p*(P31d-P_e_31);
        lambda_e_32=lambda_e_32+p*(P32d-P_e_32);
    end
    iter=iter+1;
end
%% 画图
figure
plot(Ben_Store(1,:),'b-o','LineWidth',2);
xlabel('迭代次数');
ylabel('成本/元');
title('商场的分布式迭代情况');
box off

figure
plot(Ben_Store(2,:),'r-o','LineWidth',2);
xlabel('迭代次数');
ylabel('成本/元');
title('酒店的分布式迭代情况');
box off

figure
plot(Ben_Store(3,:),'k-o','LineWidth',2);
xlabel('迭代次数');
ylabel('成本/元');
title('办公中心的分布式迭代情况');
box off
%
Party=Ben_Store(1,:)+Ben_Store(2,:)+Ben_Store(3,:);
figure
plot(Party,'k-o','LineWidth',1.5);
xlabel('迭代次数');
ylabel('成本/元');
title('商业园区用户协同交易系统总效益值');
box off

%交易部分
P_e_MG1=P_e_12+P_e_13;  
P_e_MG2=P_e_21+P_e_23;
P_e_MG3=P_e_31+P_e_32;

%电能交易结果
figure
plot(P_e_MG1,'c-*','LineWidth',2);
hold on
plot(P_e_MG2,'b-^','LineWidth',2);
hold on
plot(P_e_MG3,'g-v','LineWidth',2);
axis([1 24 -400 400]);
grid minor;
xlabel('\fontname{宋体}时间\fontname{Times new roman}/h');
ylabel('\fontname{宋体}用户间的电能交易量\fontname{Times new roman}/kWh');
% title('用户间电能交易结果');
hold on
legend('\fontname{宋体}商场\fontname{Times new roman}','\fontname{宋体}酒店\fontname{Times new roman}','\fontname{宋体}办公中心\fontname{Times new roman}','Location', 'northwest');
grid minor;
legend('boxoff');
grid off
box off

%% 结果展示1
% 第一个图：各时段功率偏差分布（所有场景）
figure; % 创建新的图形窗口
T = 24;
N_w = 10;
boxplot(squeeze(P_loss1(1, :, :))', 'Labels', 1:T);
xlabel('时段');
ylabel('功率偏差 (kW)');
title('1-各时段功率偏差分布（所有场景）');
grid on;

fprintf('1条件风险值(CVaR1): %.2f\n', CVaR1);
fprintf('1风险成本: %.2f\n', 0.1*CVaR1);
%% 结果展示2
% 第一个图：各时段功率偏差分布（所有场景）
figure; % 创建新的图形窗口
T = 24;
boxplot(squeeze(P_loss2(1, :, :))', 'Labels', 1:T);
xlabel('时段');
ylabel('功率偏差 (kW)');
title('2-各时段功率偏差分布（所有场景）');
grid on;

fprintf('2条件风险值(CVaR2): %.2f\n', CVaR2);
fprintf('2风险成本: %.2f\n', 0.1*CVaR2);
%% 结果展示3
% 第一个图：各时段功率偏差分布（所有场景）
figure; % 创建新的图形窗口
T = 24;
boxplot(squeeze(P_loss3(1, :, :))', 'Labels', 1:T);
xlabel('时段');
ylabel('功率偏差 (kW)');
title('3-各时段功率偏差分布（所有场景）');
grid on;

fprintf('3条件风险值(CVaR3): %.2f\n', CVaR3);
fprintf('3风险成本: %.2f\n', 0.1*CVaR3);
%%
set(gca,'Fontname', 'Times New Roman');%坐标轴字体设置为Times New Roman
ax2 = axes('Position',get(gca,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
set(ax2,'YTick', []);
set(ax2,'XTick', []);

%% 保存变量数据
% save('3.mat')      % 保存所有变量的数据