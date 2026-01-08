function [P_e_21,P_e_23,P_e_cut2,Obj_MG2,P_loss2,total_loss2,CVaR2]=Fun_MG2(P21d,P23d,lambda_e_21,lambda_e_23,p)
%酒店的分布式优化迭代模型
%% 决策变量初始化
L_e=sdpvar(1,24);       %需求响应后实际的电负荷
P_e_cut2=sdpvar(1,24);  %可削减电负荷
E_bat=sdpvar(1,24);     %储电设备的储电余量
P_batc=sdpvar(1,24);    %储电设备的充电功率
P_batd=sdpvar(1,24);    %储电设备的放电功率
U_abs=binvar(1,24);     %储电设备的放电状态位,取1时为放电,0为未放电
U_relea=binvar(1,24);   %储电设备的充电状态位,取1时为充电,0为未充电
P_e_pv=sdpvar(1,24);    %光伏的实际出力值
P_buy=sdpvar(1,24);     %向电网购买的电功率
P_sell=sdpvar(1,24);    %向电网售出的电功率
P_e_21=sdpvar(1,24);    %2向1交互的电功率
P_e_23=sdpvar(1,24);    %2向3交互的电功率

%% 导入电负荷和电网购电电价
%酒店
L_e0=[486,480,466,480,474,502,630,878,1010,1120,1502,1546,...
      1232,1326,1474,1566,1578,1606,1616,1680,1565,1482,1040,670]; 

Predict_pv=[0,0,0,0,0,80,299,355,446,512,574,621,...
            652,619,552,450,300,230,60,0,0,0,0,0];

pri_e=[0.4711*ones(1,8),0.8759*ones(1,6),1.0947*ones(1,3),0.8759*ones(1,2),1.0947*ones(1,3),0.8759*ones(1,2)];
grid_sw=[0.2*ones(1,24)]; 
%% 约束条件
C=[];

%% 决策变量初始化 VES
T    = ones(1,24);               
E    = ones(1,24); %t时刻虚拟储能容量
Esum = ones(1,24);
Es   = ones(1,24);
Emax = ones(1,1); 
Q_ac = ones(1,24);
P_ac = ones(1,24);
f_ac = ones(1,24);
D    = ones(1,24);
f_ac_set = ones(1,24);
P_ac_set = ones(1,24);

%% 导入数据
Ta   = [27,26,25,24,23,22,23,24,25,27,28,30,31,32,34,33,32,31,29,28,27,27,26,25];     %室外温度
Tset = ([23.2,23.3,23.4,23.2,23.1,23.1,23.0,20.4,20.2,20.0,21.5,21.4,...
        21.5,21.6,21.5,21.8,21.9,21.4,22.5,22.8,22.9,22.5,23.5,23.2]-18)/7*4+18;     %设定温度
% 计算外壁的表面温度
Tout = zeros(size(Ta)); % 初始化外壁表面温度向量
for i = 1:length(Ta)
    if i == 1
        Tout(i) = (Ta(i) + Ta(i+1)) / 2; % 使用室外温度和下一个时刻的室外温度的平均值作为外壁表面温度的初始值
    elseif i == length(Ta)
        Tout(i) = (Ta(i-1) + Ta(i)) / 2; % 使用室外温度和上一个时刻的室外温度的平均值作为外壁表面温度的最终值
    else
        Tout(i) = (Ta(i-1) + Ta(i) + Ta(i+1)) / 3; % 使用室外温度和前后两个时刻的室外温度的平均值作为外壁表面温度的估计值
    end
end

% 计算内壁的表面温度
Tin = (Ta + 0.93 * Tout) / (1 + 0.93); % 使用传热系数计算内壁表面温度

Tmax = 22;     %人体舒适最高温度
Tmin = 18;     %人体舒适最低温度
R_in = 1.225;  %空气密度kg/m3
c_in = 1005.4; %空气比热容J/(kg*°C)
Kwa  = 0.930;  %壁面传热系数
Kwi  = 0.760;  %窗口传热系数

Awi  = 76000;   %窗口表面积
Aout = 190000;   %外壁表面积
Ain  = 200000;   %内壁表面积
V1 = 32000;

a1=0.0194;
b1=-0.0773;
f1=20;
f2=57.4;
f3=97;
a2=0.0496;
b2=0.4593;
a3=0.0250;
b3=1.8702;

syms C ;     %定义符号变量
C = Aout*Kwa+Ain*Kwa + Awi*Kwi;

for t=1:24
    D(t) = Aout*Kwa*Tout(t) + Ain*Kwa*Tin(t) + Awi*Kwi*Ta(t);
    T(t) = 1/C * (Q_ac(t) + exp(-C*t/(R_in*c_in*V1))*(Ta(t)-(Q_ac(t)+D(t))/C)*C + D(t));
end

%通过频率f_ac建立功率和冷负荷关系
P_ac(t) = a1*f_ac(t)+b1;   
M=800; %这里的M是个很大的数
y1 = binvar(1, 1, 'full'); %定义了2进制变量
y2 = binvar(1, 1, 'full');
for t=1:24
    Q_ac(t) = (a2*f_ac(t)+b2)*y1 + (a3*f_ac(t)+b3)*y2; 
    f1    <=  f_ac(t)  <=  M*(1-y1)+f2;
    f2*y2 <=  f_ac(t)  <=  M*(1-y2)+f3;
    y1+y2 == 1;
end

%当接收调度指令时
for t=1:24
    f_ac_set(t)=(C*Tset(t)-D(t)-b2*y1-b3*y2)/(a2*y1+a3*y2);
    P_ac_set(t)=a1*f_ac_set(t) + b1;
    20  <= f_ac_set(t) <= 97;
    0.3 <= P_ac_set(t) <= 1.81;
end

%变频空调虚拟储能
for t=1:24
    Emax = R_in*c_in*V1*(Tmax-Tmin); 
    E(t) = R_in*c_in*V1*(Tmax-T(t)); 
end

%荷电状态
for t=1:24
    soc(t) = (Tmax-T(t))/(Tmax-Tmin); 
end

%虚拟储能容量
for t=1:24
    Esum(t) = R_in*c_in*(Tmax-Tset(t))*(V1*12);
    Es(t) = R_in*c_in*(Tset(t)-Tmin)*(V1*12);
end
beta=0.1;
%% 
%电负荷需求响应部分
for t=1:24
    C=[C,
       L_e(t)==L_e0(t)+P_e_cut2(t),                  %电负荷功率平衡约束
       -(-0.1*beta+1.01)*Es(t)/3600000<=P_e_cut2(t)<=(-0.1*beta+1.01)*Esum(t)/3600000, %可削减电功率上下限约束
      ];
end

%储电设备约束部分
%储能电站荷电状态连续性约束
C=[C,E_bat(1)==500+0.9*P_batc(1)-P_batd(1)/0.9,]; %1时段约束
for t=2:24
    C=[C,E_bat(t)==E_bat(t-1)+0.9*P_batc(t)-P_batd(t)/0.9,]; %储电设备容量变化约束
end
%储能容量大小约束
for t=1:24
    C=[C,200<=E_bat(t)<=800,];  %储电量上下限约束
end
%始末状态守恒
% C=[C,E_bat(24)==500,];
%储能电站的充放电功率约束,Big-M法进行线性化处理
M=500; %这里的M是个很大的数
for t=1:24
    C=[C,
       0<=P_batc(t)<=600,
       0<=P_batc(t)<=U_abs(t)*M,     
       0<=P_batd(t)<=600,      
       0<=P_batd(t)<=U_relea(t)*M,
       U_abs(t)+U_relea(t)<=1,
      ];
end

%PV设备运行约束
for t=1:24
    C=[C,
       0<=P_e_pv(t)<=Predict_pv(t), %光伏发电上下限约束
      ];
end

%电负荷平衡约束
for t=1:24
    C=[C,
       P_e_pv(t)+P_buy(t)+P_batd(t)==P_batc(t)+L_e(t)+P_e_21(t)+P_e_23(t)+P_sell(t),
      ];
end

%变量非负性等约束
for t=1:24
    C=[C,
       2000>=P_buy(t)>=0,
       0<=P_sell(t)<=5000,
       1000>=P_e_21(t)>=-1000,
       1000>=P_e_23(t)>=-1000,
      ];
end

%% 参数设置
Z = 1;          % 用户数量
T = 24;         % 调度时段（小时）
N_w = 10;       % 场景数量

alpha1=0.95;
% 固定惩罚电价（每个时段相同）
lambda_loss2_t = 1 * ones(T, 1);  % 惩罚电价向量

% 给定的实际响应功率偏差数据（24个时段）
base_power_deviation1 = (-0.1*beta+1.01)*[390 397 405 390 382 382 375 180 165 150 262 255 262 270 ...
     262 285 292 255 337 360 367 337 412 390]-20;  % 正整数值

% 生成场景数据（在基准值上添加±10%的随机波动）
P_loss2 = zeros(Z, T, N_w);  % 功率偏差矩阵
rng(42);  % 设置随机种子确保结果可重现
for k = 1:N_w
    % 在基准值上添加随机波动（90%~110%）
    fluctuation = 0.9 + 0.2 * rand(1, T);  % 生成[0.9, 1.1]的随机波动
    P_loss2(1, :, k) = base_power_deviation1 .* fluctuation;
end

% 场景概率（均匀分布）
rho = ones(1, N_w) / N_w;

%% 计算风险指标
[VaR2, CVaR2, total_loss2] = calculateRiskMetrics(Z, T, N_w, beta, lambda_loss2_t, P_loss2, rho);



%% 风险值计算函数
function [VaR2, CVaR2, F_loss_total] = calculateRiskMetrics(Z, T, N_w, eta, lambda_loss2_t, P_loss2, rho)
    % 计算每个场景下每个用户的损失费用
    F_i_loss_t = zeros(Z, T, N_w);  % 时刻损失矩阵
    F_i_loss_k = zeros(Z, N_w);     % 用户场景总损失矩阵
    
    for i = 1:Z
        for t = 1:T
            for k = 1:N_w
                % 计算t时刻用户i在场景k的损失费用
                F_i_loss_t(i, t, k) = lambda_loss2_t(t) * P_loss2(i, t, k);
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
    
    % 计算每个用户的VaR2
    VaR2 = zeros(Z, 1);  % 初始化VaR2向量
    
    for i = 1:Z
        % 获取当前用户在所有场景下的损失值
        user_losses = squeeze(F_i_loss_k(i, :));
        
        % 按损失值升序排序
        [sorted_losses, sort_idx] = sort(user_losses);
        
        % 计算累积概率分布
        cum_prob = cumsum(rho(sort_idx));
        
        % 找到满足累积概率 ≥ eta 的最小损失值（VaR2定义）
        idx = find(cum_prob >= eta, 1, 'first');
        
        if isempty(idx)
            VaR2(i) = sorted_losses(end);  % 若未找到则取最大值
        else
            VaR2(i) = sorted_losses(idx);
        end
    end
    
    % 计算每个用户的CVaR2
    CVaR2 = zeros(Z, 1);  % 初始化CVaR2向量
    
    for i = 1:Z
        user_losses = squeeze(F_i_loss_k(i, :));
        
        % 计算超出VaR2的损失部分 [F_i_loss_k - alpha1]+
        excess_loss = max(user_losses - alpha1, 0);
        
        % 计算CVaR2
        CVaR2(i) = alpha1 + (1/(N_w * (1 - beta))) *0.5654* sum(F_i_loss_k);
    end
end
%% 目标函数
Obj=0.08*sum(P_e_pv)+0.02*sum(P_batc+P_batd)+sum(pri_e.*P_buy)-sum(grid_sw.*P_sell)...
    -0.1*sum(P_e_cut2)-0.2*sum(P_e_21+P_e_23)+0.1*CVaR2...
    +0.03*sum(abs(P_e_21)+abs(P_e_23))...
    +sum(0.5*p*(norm(P21d-P_e_21)^2))+sum(lambda_e_21.*(P21d-P_e_21))...
    +sum(0.5*p*(norm(P23d-P_e_23)^2))+sum(lambda_e_23.*(P23d-P_e_23));    
%% 求解器配置与求解
% 设置求解器
fprintf('solver')
options=sdpsettings('solver','cplex'); %knitro
options.debug=1;
options.showprogress=1;
optimize(C,Obj,options);
%% 数据输出
P_e_21=double(P_e_21);
P_e_23=double(P_e_23);
P_e_cut2=double(P_e_cut2);
Obj_MG2=double(Obj);
P_loss2=double(P_loss2);
total_loss2=double(total_loss2);
CVaR2=double(CVaR2);
%disp(strcat('------【JD Calculate Finished!!】----------'));  
end