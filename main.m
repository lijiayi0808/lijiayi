clc
clear
close all

%% ADMM迭代参数设置
lambda_e_12 = zeros(1,96);
lambda_e_13 = zeros(1,96);
lambda_e_21 = zeros(1,96);
lambda_e_23 = zeros(1,96);
lambda_e_31 = zeros(1,96);
lambda_e_32 = zeros(1,96);

Lambda_e_12 = zeros(1,96);
Lambda_e_13 = zeros(1,96);
Lambda_e_21 = zeros(1,96);
Lambda_e_23 = zeros(1,96);
Lambda_e_31 = zeros(1,96);
Lambda_e_32 = zeros(1,96);

%% 辅助变量
P12d = 200 * ones(1,96);
P13d = 200 * ones(1,96);
P21d = 200 * ones(1,96);
P23d = 200 * ones(1,96);
P31d = 200 * ones(1,96);
P32d = 200 * ones(1,96);

maxIter = 50;

% 价格初始化
pri_e_12=zeros(maxIter+1,96);pri_e_21=zeros(maxIter+1,96);
pri_e_13=zeros(maxIter+1,96);pri_e_31=zeros(maxIter+1,96);
pri_e_23=zeros(maxIter+1,96);pri_e_32=zeros(maxIter+1,96);

p = 0.01;  % 初始惩罚因子
iter = 1;
Ben_Store = [];
Ben_Store2 = [];
alpha_k = 1;
prev_residual = inf;  % 上一次迭代的残差初始化为无穷大

tolerant=1e-1;%收敛精度
toler1=[]; %残差1,电能交易部分

% 设置μ+、μ-和υ的值
mu_plus = 1.14;
mu_minus = 1.14;
v = 1.1;

%% 迭代


    for k = 1:maxIter
        % Step 1: 更新 alpha_k+1
        alpha_k_next = (1 + sqrt(1 + 4 * alpha_k^2)) / 2;
    
        % 预测更新
        P12d_prev = P12d;
        P13d_prev = P13d;
        P21d_prev = P21d;
        P23d_prev = P23d;
        P31d_prev = P31d;
        P32d_prev = P32d;
    
        % 预测更新
        P12d = P12d + (alpha_k - 1) / alpha_k_next * (P12d - P12d_prev);
        P13d = P13d + (alpha_k - 1) / alpha_k_next * (P13d - P13d_prev);
        P21d = P21d + (alpha_k - 1) / alpha_k_next * (P21d - P21d_prev);
        P23d = P23d + (alpha_k - 1) / alpha_k_next * (P23d - P23d_prev);
        P31d = P31d + (alpha_k - 1) / alpha_k_next * (P31d - P31d_prev);
        P32d = P32d + (alpha_k - 1) / alpha_k_next * (P32d - P32d_prev);
    
        % 更新 alpha_k 为 alpha_k+1
        alpha_k = alpha_k_next;
    
        % 计算 Fun_MG1、Fun_MG2、Fun_MG3 函数
        [P_e_12, P_e_13, ~, Obj_MG1(iter) , trcost1, CVaR1] = Fun_MG1(P12d, P13d, lambda_e_12, lambda_e_13, p);
        [P_e_21, P_e_23, ~, Obj_MG2(iter) , trcost2, CVaR2] = Fun_MG2(P21d, P23d, lambda_e_21, lambda_e_23, p);
        [P_e_31, P_e_32, ~, Obj_MG3(iter) , trcost3, CVaR3] = Fun_MG3(P31d, P32d, lambda_e_31, lambda_e_32, p);
    
        % 预测-修正更新
        P12d_new = (lambda_e_21 - lambda_e_12 + p * (P_e_12 - P_e_21)) / (2 * p);
        P13d_new = (lambda_e_31 - lambda_e_13 + p * (P_e_13 - P_e_31)) / (2 * p);
        P21d_new = (lambda_e_12 - lambda_e_21 + p * (P_e_21 - P_e_12)) / (2 * p);
        P23d_new = (lambda_e_32 - lambda_e_23 + p * (P_e_23 - P_e_32)) / (2 * p);
        P31d_new = (lambda_e_13 - lambda_e_31 + p * (P_e_31 - P_e_13)) / (2 * p);
        P32d_new = (lambda_e_23 - lambda_e_32 + p * (P_e_32 - P_e_23)) / (2 * p);
    
        % 动态调整步长
        step_size = 1.25; % 增加步长的因子
        P12d = P12d + step_size * (P12d_new - P12d);
        P13d = P13d + step_size * (P13d_new - P13d);
        P21d = P21d + step_size * (P21d_new - P21d);
        P23d = P23d + step_size * (P23d_new - P23d);
        P31d = P31d + step_size * (P31d_new - P31d);
        P32d = P32d + step_size * (P32d_new - P32d);
        
        if iter==1  %第一次求解比较特殊,其他主体初值为0
        [pri_e_12(2,:),pri_e_13(2,:),Obj_MG1p(iter)]=Fun_MG1_p(pri_e_21(iter,:),pri_e_31(iter,:),P_e_12, P_e_13,Lambda_e_12,Lambda_e_13);
        [pri_e_21(2,:),pri_e_23(2,:),Obj_MG2p(iter)]=Fun_MG2_p(pri_e_12(2,:),pri_e_32(iter,:),P_e_21, P_e_23,Lambda_e_21,Lambda_e_23);
        [pri_e_31(2,:),pri_e_32(2,:),Obj_MG3p(iter)]=Fun_MG3_p(pri_e_13(2,:),pri_e_23(2,:),P_e_31, P_e_32,Lambda_e_31,Lambda_e_32);
        Lambda_e_12=Lambda_e_12+10*(pri_e_12(2,:)-pri_e_21(2,:));
        Lambda_e_13=Lambda_e_13+10*(pri_e_13(2,:)-pri_e_31(2,:));
        Lambda_e_21=Lambda_e_21+10*(pri_e_21(2,:)-pri_e_12(2,:));
        Lambda_e_23=Lambda_e_23+10*(pri_e_23(2,:)-pri_e_32(2,:));
        Lambda_e_31=Lambda_e_31+10*(pri_e_31(2,:)-pri_e_13(2,:));
        Lambda_e_32=Lambda_e_32+10*(pri_e_32(2,:)-pri_e_23(2,:));
        else
        [pri_e_12(iter+1,:),pri_e_13(iter+1,:),Obj_MG1p(iter)]=Fun_MG1_p(pri_e_21(iter,:),pri_e_31(iter,:),P_e_12, P_e_13,Lambda_e_12,Lambda_e_13);
        [pri_e_21(iter+1,:),pri_e_23(iter+1,:),Obj_MG2p(iter)]=Fun_MG2_p(pri_e_12(iter+1,:),pri_e_32(iter,:),P_e_21, P_e_23,Lambda_e_21,Lambda_e_23);
        [pri_e_31(iter+1,:),pri_e_32(iter+1,:),Obj_MG3p(iter)]=Fun_MG3_p(pri_e_13(iter+1,:),pri_e_23(iter+1,:),P_e_31, P_e_32,Lambda_e_31,Lambda_e_32);
        Lambda_e_12=Lambda_e_12+10*(pri_e_12(iter+1,:)-pri_e_21(iter+1,:));
        Lambda_e_13=Lambda_e_13+10*(pri_e_13(iter+1,:)-pri_e_31(iter+1,:));
        Lambda_e_21=Lambda_e_21+10*(pri_e_21(iter+1,:)-pri_e_12(iter+1,:));
        Lambda_e_23=Lambda_e_23+10*(pri_e_23(iter+1,:)-pri_e_32(iter+1,:));
        Lambda_e_31=Lambda_e_31+10*(pri_e_31(iter+1,:)-pri_e_13(iter+1,:));
        Lambda_e_32=Lambda_e_32+10*(pri_e_32(iter+1,:)-pri_e_23(iter+1,:));
        end
        
        Ben_Store = [Ben_Store, [Obj_MG1(iter); Obj_MG2(iter); Obj_MG3(iter)]];
    
        % 残差计算
        AP2P = [P_e_12; P_e_13; P_e_21; P_e_23; P_e_31; P_e_32];
        P2P = [P12d; P13d; P21d; P23d; P31d; P32d];
        primal_residual = norm(AP2P - P2P);  % 计算原始残差
        dual_residual = norm([lambda_e_12; lambda_e_13; lambda_e_21; lambda_e_23; lambda_e_31; lambda_e_32]);  % 计算对偶残差
    
        %保存历史数据
        Ben_Store2=[Ben_Store2,[Obj_MG1p(iter);Obj_MG2p(iter);Obj_MG3p(iter)]];
         % 残差计算
        P1 = [pri_e_12; pri_e_13; pri_e_23];
        P2 = [pri_e_21; pri_e_31; pri_e_32];
        toler1=[toler1,norm(pri_e_12-pri_e_21)+norm(pri_e_13-pri_e_31)+norm(pri_e_23-pri_e_32)];%保存残差1 

        % 判断收敛条件
        if primal_residual/48<=tolerant  &&  toler1(iter)/96<=tolerant
            display(['迭代收敛,在第 ', num2str(iter), ' 次收敛']);
            break;
        else
            % 更新拉格朗日乘子
            lambda_e_12 = lambda_e_12 + p * (P12d - P_e_12);
            lambda_e_13 = lambda_e_13 + p * (P13d - P_e_13);
            lambda_e_21 = lambda_e_21 + p * (P21d - P_e_21);
            lambda_e_23 = lambda_e_23 + p * (P23d - P_e_23);
            lambda_e_31 = lambda_e_31 + p * (P31d - P_e_31);
            lambda_e_32 = lambda_e_32 + p * (P32d - P_e_32);
            display(['迭代', num2str(iter), ' 次，不收敛']);
        end
    
        % 调整惩罚因子p
        if primal_residual > v * dual_residual  % 如果原始残差大于对偶残差的υ倍
            p = p * mu_plus;  % 更新惩罚因子p为p*μ+
        elseif dual_residual > v * primal_residual  % 如果对偶残差大于原始残差的υ倍
            p = p * mu_minus;  % 更新惩罚因子p为p*μ-
        end
    
        iter = iter + 1;
    end


%% 画图  图1-4
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

Party=Ben_Store(1,:)+Ben_Store(2,:)+Ben_Store(3,:);
figure
plot(Party,'k-o','LineWidth',1.5);
xlabel('迭代次数');
ylabel('成本/元');
title('商业园区用户协同交易系统总效益值');
box off

%% 图5 - 堆积柱状图
% 交易部分
P_e_MG1 = P_e_12 + P_e_13;  
P_e_MG2 = P_e_21 + P_e_23;
P_e_MG3 = P_e_31 + P_e_32;

% 使用所有96个数据点
time_points = 1:96;
P_e_MG1_sampled = P_e_MG1;
P_e_MG2_sampled = P_e_MG2;
P_e_MG3_sampled = P_e_MG3;

% 准备堆积柱状图数据 - 正值部分向上堆积
P_e_MG1_sampled_pos = max(P_e_MG1_sampled, 0);
P_e_MG2_sampled_pos = max(P_e_MG2_sampled, 0);
P_e_MG3_sampled_pos = max(P_e_MG3_sampled, 0);

% 负值部分单独处理 - 从0向下堆积
P_e_MG1_sampled_neg = min(P_e_MG1_sampled, 0);
P_e_MG2_sampled_neg = min(P_e_MG2_sampled, 0);
P_e_MG3_sampled_neg = min(P_e_MG3_sampled, 0);

% 创建图形窗口
figure('Position', [100, 100, 1200, 500])

% 定义三种固定颜色
mg1_color = [0, 0.5, 0.8];   % 商场 - 蓝色
mg2_color = [0.8, 0.4, 0];   % 酒店 - 橙色
mg3_color = [0.2, 0.7, 0.2]; % 办公中心 - 绿色

% 创建堆积柱状图 - 正值部分
bar_data_pos = [P_e_MG1_sampled_pos; P_e_MG2_sampled_pos; P_e_MG3_sampled_pos]';
bar_pos = bar(time_points, bar_data_pos, 'stacked', 'BarWidth', 0.8);

% 设置正值部分颜色
bar_pos(1).FaceColor = mg1_color;  % 商场 - 蓝色
bar_pos(2).FaceColor = mg2_color;  % 酒店 - 橙色
bar_pos(3).FaceColor = mg3_color;  % 办公中心 - 绿色

% 如果有负值，添加负值堆积部分
if any([P_e_MG1_sampled_neg, P_e_MG2_sampled_neg, P_e_MG3_sampled_neg] < 0)
    hold on
    bar_data_neg = [P_e_MG1_sampled_neg; P_e_MG2_sampled_neg; P_e_MG3_sampled_neg]';
    bar_neg = bar(time_points, bar_data_neg, 'stacked', 'BarWidth', 0.8);
    
    % 负值部分使用完全相同的颜色（无透明度变化）
    bar_neg(1).FaceColor = mg1_color;  % 商场 - 蓝色
    bar_neg(2).FaceColor = mg2_color;  % 酒店 - 橙色
    bar_neg(3).FaceColor = mg3_color;  % 办公中心 - 绿色
end

% 设置坐标轴
xlabel('时间/h', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('用户间的电能交易量/kWh', 'FontSize', 12, 'FontWeight', 'bold');
title('用户间电能交易结果 - 堆积柱状图', 'FontSize', 14, 'FontWeight', 'bold');
xlim([0, 97]);  % 设置X轴范围，包含所有96个点

% 设置X轴刻度
xticks([1, 24, 48, 72, 96]);  % 显示关键时间点

% 添加图例（只需要三个图例项）
legend({'商场', '酒店', '办公中心'}, 'Location', 'northwest', 'FontSize', 10);
legend('boxoff');

% 添加零线和网格
yline(0, 'k-', 'LineWidth', 1.5);
grid on
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.3);

hold off
%% 图6
%各个微网的交易部分
pri_e=[0.4711*ones(1,32),0.8759*ones(1,24),1.0947*ones(1,12),0.8759*ones(1,8),1.0947*ones(1,12),0.8759*ones(1,8)];
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
axis([1 96 0.1 1.3]);
xlabel('\fontname{宋体}时间\fontname{Times new roman}/h');
ylabel('\fontname{宋体}交易电价\fontname{Times new roman}/(\fontname{宋体}元\fontname{Times new roman}/kWh)');
% title('\fontname{宋体}交易电价');
hold on
legend('\fontname{宋体}\fontname{Times new roman}1-2\fontname{宋体}电价','\fontname{宋体}\fontname{Times new roman}1-3\fontname{宋体}电价','\fontname{宋体}\fontname{Times new roman}2-3\fontname{宋体}电价','\fontname{宋体}电网电价','Location', 'northwest');
legend('boxoff');
grid off
box off

%% 图7
% 创建柱状图
figure;
bar_data = [CVaR1, CVaR2, CVaR3];
bar(bar_data, 'FaceColor', [0.2 0.4 0.6]);
xlabel('园区编号', 'FontSize', 12);
ylabel('CVaR值', 'FontSize', 12);
title('各园区条件风险价值(CVaR)对比', 'FontSize', 14);
grid on;
box off;

% 设置x轴刻度标签
set(gca, 'XTickLabel', {'商场', '酒店', '办公中心'});
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);

% 在柱状图上添加数值标签
for i = 1:length(bar_data)
    text(i, bar_data(i)+max(bar_data)*0.02, ...
        sprintf('%.2f', bar_data(i)), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 10, 'FontWeight', 'bold');
end

% 添加图例
legend('CVaR值', 'Location', 'best');
legend('boxoff');



%% 最后
set(gca,'Fontname', 'Times New Roman');%坐标轴字体设置为Times New Roman
ax2 = axes('Position',get(gca,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
set(ax2,'YTick', []);
set(ax2,'XTick', []);

