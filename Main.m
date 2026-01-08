%子问题2:支付效益最大化

clc
clear
close all
%% ADMM迭代参数设置
%拉格朗日乘子初始化
lambda_e_12=zeros(1,24);%拉格朗日乘子
lambda_e_13=zeros(1,24);
lambda_e_21=zeros(1,24);
lambda_e_23=zeros(1,24);
lambda_e_31=zeros(1,24);
lambda_e_32=zeros(1,24);
maxIter=200;  %最大迭代次数
tolerant=1e-3;%收敛精度
iter=1;%迭代次数初始化
Ben_Store=[];%历史目标函数
toler1=[]; %残差1,电能交易部分
%交易量记录矩阵初始化
pri_e_12=zeros(maxIter+1,24);pri_e_21=zeros(maxIter+1,24);
pri_e_13=zeros(maxIter+1,24);pri_e_31=zeros(maxIter+1,24);
pri_e_23=zeros(maxIter+1,24);pri_e_32=zeros(maxIter+1,24);
%% 迭代
while 1
    if iter==maxIter  %限制迭代次数
       disp('迭代不收敛,参数有误');
       break; 
    end 
    display(['迭代还未收敛,当前迭代第 ', num2str(iter),' 次']);
    if iter==1  %第一次求解比较特殊,其他主体初值为0
        [pri_e_12(2,:),pri_e_13(2,:),Obj_MG1(iter)]=Fun_MG1(pri_e_21(iter,:),pri_e_31(iter,:),lambda_e_12,lambda_e_13);
        [pri_e_21(2,:),pri_e_23(2,:),Obj_MG2(iter)]=Fun_MG2(pri_e_12(2,:),pri_e_32(iter,:),lambda_e_21,lambda_e_23);
        [pri_e_31(2,:),pri_e_32(2,:),Obj_MG3(iter)]=Fun_MG3(pri_e_13(2,:),pri_e_23(2,:),lambda_e_31,lambda_e_32);
        lambda_e_12=lambda_e_12+10*(pri_e_12(2,:)-pri_e_21(2,:));
        lambda_e_13=lambda_e_13+10*(pri_e_13(2,:)-pri_e_31(2,:));
        lambda_e_21=lambda_e_21+10*(pri_e_21(2,:)-pri_e_12(2,:));
        lambda_e_23=lambda_e_23+10*(pri_e_23(2,:)-pri_e_32(2,:));
        lambda_e_31=lambda_e_31+10*(pri_e_31(2,:)-pri_e_13(2,:));
        lambda_e_32=lambda_e_32+10*(pri_e_32(2,:)-pri_e_23(2,:));
    else
        [pri_e_12(iter+1,:),pri_e_13(iter+1,:),Obj_MG1(iter)]=Fun_MG1(pri_e_21(iter,:),pri_e_31(iter,:),lambda_e_12,lambda_e_13);
        [pri_e_21(iter+1,:),pri_e_23(iter+1,:),Obj_MG2(iter)]=Fun_MG2(pri_e_12(iter+1,:),pri_e_32(iter,:),lambda_e_21,lambda_e_23);
        [pri_e_31(iter+1,:),pri_e_32(iter+1,:),Obj_MG3(iter)]=Fun_MG3(pri_e_13(iter+1,:),pri_e_23(iter+1,:),lambda_e_31,lambda_e_32);
        lambda_e_12=lambda_e_12+10*(pri_e_12(iter+1,:)-pri_e_21(iter+1,:));
        lambda_e_13=lambda_e_13+10*(pri_e_13(iter+1,:)-pri_e_31(iter+1,:));
        lambda_e_21=lambda_e_21+10*(pri_e_21(iter+1,:)-pri_e_12(iter+1,:));
        lambda_e_23=lambda_e_23+10*(pri_e_23(iter+1,:)-pri_e_32(iter+1,:));
        lambda_e_31=lambda_e_31+10*(pri_e_31(iter+1,:)-pri_e_13(iter+1,:));
        lambda_e_32=lambda_e_32+10*(pri_e_32(iter+1,:)-pri_e_23(iter+1,:));
    end
    %保存历史数据
    Ben_Store=[Ben_Store,[Obj_MG1(iter);Obj_MG2(iter);Obj_MG3(iter)]];
    %残差计算
    toler1=[toler1,norm(pri_e_12(iter+1,:)-pri_e_21(iter+1,:))+norm(pri_e_13(iter+1,:)-pri_e_31(iter+1,:))+norm(pri_e_23(iter+1,:)-pri_e_32(iter+1,:))];%保存残差1 
    %判断收敛条件
    if toler1(iter)<=tolerant
       display(['迭代收敛,在第 ', num2str(iter),' 次收敛']);
       break; 
    end
    iter=iter+1;
end
%% 画图
figure
plot(Ben_Store(1,:),'b-o','LineWidth',1.5);
xlabel('迭代次数');
ylabel('商场议价/元');
title('商场的分布式迭代情况');
box off
figure
plot(Ben_Store(2,:),'r-o','LineWidth',1.5);
xlabel('迭代次数');
ylabel('酒店议价/元');
title('酒店的分布式迭代情况');
box off
figure
plot(Ben_Store(3,:),'k-o','LineWidth',1.5);
xlabel('迭代次数');
ylabel('办公中心议价/元');
title('办公中心的分布式迭代情况');
box off
%各个微网的交易部分
pri_e=[0.4711*ones(1,8),0.8759*ones(1,6),1.0947*ones(1,3),0.8759*ones(1,2),1.0947*ones(1,3),0.8759*ones(1,2)];
pri_e_MG1toMG2=pri_e_12(iter+1,:); 
pri_e_MG1toMG3=pri_e_13(iter+1,:);  
pri_e_MG2toMG3=pri_e_23(iter+1,:); 
%微网之间的交易电价
figure
plot(pri_e_MG1toMG2,'c-o','LineWidth',1.5);
hold on
plot(pri_e_MG1toMG3,'b-^','LineWidth',1.5);
hold on
plot(pri_e_MG2toMG3,'g-v','LineWidth',1.5);
hold on
plot(pri_e,'m--','LineWidth',1.5);
axis([1 24 0.1 1.3]);
xlabel('\fontname{宋体}时间\fontname{Times new roman}/h');
ylabel('\fontname{宋体}交易电价\fontname{Times new roman}/(\fontname{宋体}元\fontname{Times new roman}/kWh)');
% title('\fontname{宋体}交易电价');
hold on
legend('\fontname{宋体}\fontname{Times new roman}1-2\fontname{宋体}电价','\fontname{宋体}\fontname{Times new roman}1-3\fontname{宋体}电价','\fontname{宋体}\fontname{Times new roman}2-3\fontname{宋体}电价','\fontname{宋体}电网电价','Location', 'northwest');
legend('boxoff');
grid off
box off
set(gca,'Fontname', 'Times New Roman');%坐标轴字体设置为Times New Roman
ax2 = axes('Position',get(gca,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
set(ax2,'YTick', []);
set(ax2,'XTick', []);