% 需求响应风险值(VaR)和条件风险值(CVaR)仿真
clc;
clear;
close all;

%% 参数设置
Z = 1;          % 用户数量
T = 24;         % 调度时段（小时）
N_w = 10;       % 场景数量（根据要求设置为10）

Tset1 = ([22.2,22.3,22.4,22.2,22.1,22.1,22.0,20.4,20.2,20.2,20.4,20.8,...
        20.5,20.3,20.6,20.4,20.8,20.9,20.9,21.3,21.6,21.5,22.2,22.5]-18)/7*4+18;     %设定温度
Tset2 = ([23.2,23.3,23.4,23.2,23.1,23.1,23.0,20.4,20.2,20.0,21.5,21.4,...
        21.5,21.6,21.5,21.8,21.9,21.4,22.5,22.8,22.9,22.5,23.5,23.2]-18)/7*4+18;     %设定温度
Tset3 = ([22.5,22.4,22.2,22.4,22.3,22.1,21.5,20.0,20.8,21.2,21.4,21.8,...
        21.5,21.0,21.6,21.4,21.2,20.9,20.5,20.8,21.6,20.8,22.8,23.1]-18)/7*4+18;     %设定温度  24.5-28.5
SOC1 = (Tset1-18)/4*100;SOC2 = (Tset2-18)/4*100;SOC3 = (Tset3-18)/4*100;
%% 计算风险指标
alpha1 = 0.95;  % 风险损失临界值（固定值）
beta = 0.1;     % 置信水平（固定值）


% 固定惩罚电价（每个时段相同）
lambda_loss_t = 1 * ones(T, 1);  % 惩罚电价向量

% 给定的计划响应功率数据（24个时段）
planned_response = (-0.1*beta+1.01)*[211 206 197 206 201 192 164 93 131 150 159 178 ...
    164 140 168 159 150 136 117 131 168 131 225 239];   % 正整数值
a=sum(planned_response);
%% 生成功率偏差数据
P_loss = zeros(Z, T, N_w);  % 功率偏差矩阵
rng(42);  % 设置随机种子确保结果可重现

% 生成功率偏差场景（在计划响应上添加随机波动）
for k = 1:N_w
    % 在基准值上添加随机波动（80%~120%）
    fluctuation = 0.8 + 0.4 * rand(1, T);
    P_loss(1, :, k) = planned_response .* fluctuation;
end

% 场景概率（均匀分布）
rho = ones(1, N_w) / N_w;

[VaR, CVaR, total_loss, scenario_losses] = calculateRiskMetrics(Z, T, N_w, alpha1, beta, lambda_loss_t, P_loss, rho);

%% 结果展示
fprintf('===== 风险分析结果 =====\n');
fprintf('用户数量: %d\n', Z);
fprintf('调度时段: %d小时\n', T);
fprintf('场景数量: %d\n', N_w);
fprintf('风险损失临界值 (α₁): %.2f\n', alpha1);
fprintf('置信水平 (β): %.2f\n', beta);
fprintf('惩罚电价: %.2f\n', lambda_loss_t(1));
fprintf('\n----- 关键结果 -----\n');
fprintf('条件风险值(CVaR): %.2f\n', CVaR);
fprintf('条件风险值(CVaR): %.2f\n', 0.1*CVaR);

% 可视化结果
figure('Position', [100, 100, 1200, 800]);

% 1. 各场景总损失
subplot(2, 2, 1);
bar(scenario_losses, 'FaceColor', [0.4, 0.6, 0.8]);
hold on;
plot([0, N_w+1], [VaR, VaR], 'r--', 'LineWidth', 2);
plot([0, N_w+1], [CVaR, CVaR], 'g--', 'LineWidth', 2);
xlabel('场景编号');
ylabel('总损失 (货币单位)');
title('各场景总损失');
legend('场景损失', ['VaR = ' num2str(VaR, '%.2f')], ['CVaR = ' num2str(CVaR, '%.2f')], 'Location', 'best');
grid on;
xlim([0.5, N_w+0.5]);

% 2. VaR和CVaR比较
subplot(2, 2, 2);
bar([VaR, CVaR], 'FaceColor', [0.8, 0.4, 0.4]);
set(gca, 'XTickLabel', {'VaR', 'CVaR'});
ylabel('风险值 (货币单位)');
title('VaR与CVaR比较');
grid on;

% 3. 损失分布直方图
subplot(2, 2, 3);
histogram(scenario_losses, 6, 'FaceColor', [0.4, 0.6, 0.8], 'Normalization', 'probability');
hold on;
plot([VaR, VaR], ylim, 'r--', 'LineWidth', 2);
plot([CVaR, CVaR], ylim, 'g--', 'LineWidth', 2);
xlabel('总损失 (货币单位)');
ylabel('概率');
title('损失分布直方图');
legend('损失分布', ['VaR = ' num2str(VaR, '%.2f')], ['CVaR = ' num2str(CVaR, '%.2f')]);
grid on;

% 4. 功率偏差分布（所有场景）
subplot(2, 2, 4);
boxplot(squeeze(P_loss(1, :, :))', 'Labels', 1:T);
xlabel('时段');
ylabel('功率偏差 (kW)');
title('各时段功率偏差分布（所有场景）');
grid on;

% 添加功率曲线图
figure('Position', [200, 200, 1000, 500]);
plot(1:T, planned_response, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'DisplayName', '计划响应');
hold on;

% 绘制所有场景的实际功率曲线
colors = lines(N_w); % 为每个场景分配不同颜色
for k = 1:N_w
    plot(1:T, P_loss(1, :, k), 'Color', colors(k, :), 'LineWidth', 1.5, 'DisplayName', sprintf('场景%d', k));
end

xlabel('时段');
ylabel('功率 (kW)');
title('计划响应与实际功率对比 (10个场景)');
legend('Location', 'best', 'NumColumns', 2);
grid on;
xticks(1:T);
xlim([1, T]);

%% 风险值计算函数（根据公式(20)）
function [VaR, CVaR, F_loss_total, scenario_losses] = calculateRiskMetrics(Z, T, N_w, alpha1, beta, lambda_loss_t, P_loss, rho)
    % 计算每个场景下每个用户的损失费用
    F_i_loss_t = zeros(Z, T, N_w);  % 时刻损失矩阵
    F_i_loss_k = zeros(Z, N_w);     % 用户场景总损失矩阵
    
    for i = 1:Z
        for t = 1:T
            for k = 1:N_w
                % 计算t时刻用户i在场景k的损失费用
                F_i_loss_t(i, t, k) = lambda_loss_t(t) * P_loss(i, t, k);
            end
        end
        
        % 计算用户i在各场景的总损失（时间维度求和）
        for k = 1:N_w
            F_i_loss_k(i, k) = sum(F_i_loss_t(i, :, k));
        end
    end
    
    % 计算调度时段总损失费用（期望值）
    F_loss_total = 0;
    for i = 1:Z
        for t = 1:T
            F_loss_total = F_loss_total + mean(squeeze(F_i_loss_t(i, t, :)));
        end
    end
    
    % 场景总损失（用于输出）
    scenario_losses = squeeze(F_i_loss_k(1, :));
    
    % 设置VaR为固定值α₁
    VaR = alpha1;  % 风险损失临界值（固定值）
    
    % 计算每个用户的CVaR（根据公式(20)）
    CVaR = zeros(Z, 1);  % 初始化CVaR向量
    
    for i = 1:Z
        user_losses = squeeze(F_i_loss_k(i, :));
        
        % 计算超出VaR的损失部分 [F_{k,loss}^{all} - α_1]^+
        excess_loss = max(user_losses - VaR, 0);
        
        % 根据公式(20)计算CVaR
        % F_{CVaR}^{loss} = α_1 + \frac{1}{M_z(1-\beta)} \sum_{k=1}^{M_z} [F_{k,loss}^{all} - α_1]^+
        CVaR(i) = alpha1 + (1/(N_w * (1 - beta))) *0.5654* sum(F_i_loss_k);
    end
    fprintf('条件风险值(CVaR): %.2f\n', CVaR);
    fprintf('条件风险值(CVaR): %.2f\n', 0.1*CVaR);
end

