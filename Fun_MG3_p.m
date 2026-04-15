function [pri_e_31,pri_e_32,Obj_MG3p]=Fun_MG3_p(pri_e_13,pri_e_23,P_e_31, P_e_32,Lambda_e_31,Lambda_e_32)
%微网3(MG3)的分布式优化迭代模型
%% 决策变量初始化
pri_e_31=sdpvar(1,96); %微网3向微网1交互的电价
pri_e_32=sdpvar(1,96); %微网3向微网2交互的电价

%% 导入电/热负荷和电网购电电价
pri_e = [0.4711*ones(1,32),0.8759*ones(1,24),1.0947*ones(1,12),0.8759*ones(1,8),1.0947*ones(1,12),0.8759*ones(1,8)];

%% 约束条件
C=[];
for t=1:96
    C=[C,
       pri_e(t)>=pri_e_31(t)>=0.2,
       pri_e(t)>=pri_e_32(t)>=0.2,

      ];
end
%% 目标函数
C_trade=sum(pri_e_31.*P_e_31+pri_e_32.*P_e_32);
Obj=-log(C_trade)...
    +0.5*1e1*(norm(pri_e_31-pri_e_13)^2)+sum(Lambda_e_31.*(pri_e_31-pri_e_13))...
    +0.5*1e1*(norm(pri_e_32-pri_e_23)^2)+sum(Lambda_e_32.*(pri_e_32-pri_e_23));    
%% 求解器配置与求解
ops=sdpsettings('solver','mosek','verbose',0,'usex0',0);
result=solvesdp(C,Obj,ops);
%% 数据输出
pri_e_31=double(pri_e_31);
pri_e_32=double(pri_e_32);
Obj_MG3p=double(C_trade);
end