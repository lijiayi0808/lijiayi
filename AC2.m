clc
clear
%
%% 决策变量初始化
a    = ones(1,24); 
b    = ones(1,24);
T    = ones(1,24);               
E    = ones(1,24); %t时刻虚拟储能容量
Esum = ones(1,24);
Es = ones(1,24);
Emax = ones(1,1); 
Q_ac = ones(1,24);
P_ac = ones(1,24);
f_ac = ones(1,24);
D    = ones(1,24);
f_ac_set = ones(1,24);
P_ac_set = ones(1,24);
%% 导入数据
Ta   = [27,26,25,24,23,22,23,24,25,27,28,30,31,32,34,33,32,31,29,28,27,27,26,25]; 
Tset = ([23.2,23.3,23.4,23.2,23.1,23.1,23.0,20.4,20.2,20.0,21.5,21.4,...
        21.5,21.6,21.5,21.8,21.9,21.4,22.5,22.8,22.9,22.5,23.5,23.2]-18)/7*4+18; 
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
%% 求解

%% 求解
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
    Emax = R_in*c_in*(V1*12)*(Tmax-Tmin)/3600000; 
    E(t) = R_in*c_in*V1*(Tmax-T(t)); 
end

%荷电状态
for t=1:24
    soc(t) = (Tmax-T(t))/(Tmax-Tmin); 
end

%虚拟储能容量
for t=1:24
    Esum(t) = R_in*c_in*(Tmax-Tset(t))*(V1*12); 
    a(t) = Esum(t)/3600000;
    Es(t) = R_in*c_in*(Tset(t)-Tmin)*(V1*12); 
    b(t) = Es(t)/3600000;
end



% 绘制虚拟储能量随时间变化的柱状图
figure;
plot(1:24, Esum/3600000, 'b');
xlabel('时间/h');
ylabel('虚拟储能量/kW');
title('虚拟储能量随时间变化');
grid on;

