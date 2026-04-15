function [pri_e_21,pri_e_23,Obj_MG2p]=Fun_MG2_p(pri_e_12,pri_e_32,P_e_21,P_e_23,Lambda_e_21,Lambda_e_23)
%微网2(MG2)的分布式优化迭代模型
%% 决策变量初始化
pri_e_21=sdpvar(1,96); %微网2向微网1交互的电价
pri_e_23=sdpvar(1,96); %微网2向微网3交互的电价
%% 导入电/热负荷和电网购电电价
pri_e = [0.4711*ones(1,32),0.8759*ones(1,24),1.0947*ones(1,12),0.8759*ones(1,8),1.0947*ones(1,12),0.8759*ones(1,8)];

%% 约束条件
C=[];
for t=1:96
    C=[C,
       pri_e(t)>=pri_e_21(t)>=0.2,
       pri_e(t)>=pri_e_23(t)>=0.2,

      ];
end
%% 目标函数
C_trade=sum(pri_e_21.*P_e_21+pri_e_23.*P_e_23);
Obj=-log(C_trade)...
    +0.5*1e1*(norm(pri_e_21-pri_e_12)^2)+sum(Lambda_e_21.*(pri_e_21-pri_e_12))...
    +0.5*1e1*(norm(pri_e_23-pri_e_32)^2)+sum(Lambda_e_23.*(pri_e_23-pri_e_32));    
%% 求解器配置与求解
ops=sdpsettings('solver','mosek','verbose',0,'usex0',0);
result=solvesdp(C,Obj,ops);
%% 数据输出
pri_e_21=double(pri_e_21);
pri_e_23=double(pri_e_23);
Obj_MG2p=double(C_trade);
end