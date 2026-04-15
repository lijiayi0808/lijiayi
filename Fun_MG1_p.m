function [pri_e_12,pri_e_13,Obj_MG1p]=Fun_MG1_p(pri_e_21,pri_e_31,P_e_12, P_e_13,Lambda_e_12,Lambda_e_13)
%微网1(MG1)的分布式优化迭代模型
%% 决策变量初始化
pri_e_12=sdpvar(1,96); %微网1向微网2交互的电价
pri_e_13=sdpvar(1,96); %微网1向微网3交互的电价

%% 导入电/热负荷和电网购电电价
pri_e = [0.4711*ones(1,32),0.8759*ones(1,24),1.0947*ones(1,12),0.8759*ones(1,8),1.0947*ones(1,12),0.8759*ones(1,8)];
%% 约束条件
C=[];
for t=1:96
    C=[C,
       pri_e(t)>=pri_e_12(t)>=0.2,
       pri_e(t)>=pri_e_13(t)>=0.2,

      ];
end

%% 目标函数
C_trade=sum(pri_e_12.*P_e_12+pri_e_13.*P_e_13);
Obj=-log(C_trade)...
    +0.5*1e1*(norm(pri_e_12-pri_e_21)^2)+sum(Lambda_e_12.*(pri_e_12-pri_e_21))...
    +0.5*1e1*(norm(pri_e_13-pri_e_31)^2)+sum(Lambda_e_13.*(pri_e_13-pri_e_31));    
%% 求解器配置与求解
ops=sdpsettings('solver','mosek','verbose',0,'usex0',0);
result=solvesdp(C,Obj,ops);
%% 数据输出
pri_e_12=double(pri_e_12);
pri_e_13=double(pri_e_13);
Obj_MG1p=double(C_trade);
end