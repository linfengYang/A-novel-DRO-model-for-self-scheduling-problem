%%  DR_CVaR   by yy 2020.11.10

% function [k,time,epsilon]=DRO_CVaR_Alg1(pathAndFilename)

% FileName = 'OPF_DAT/IEEE4.dat';
% FileName = 'SCUC_dat/SCUC6_4period.txt';
% FileName = 'SCUC_dat/SCUC118.txt';
% FileName = 'SCUC_dat/SCUC30.txt';
FileName = 'SCUC_dat/SCUC6.txt';
SCUC_data = ReadDataSCUC(FileName);

% % MonteCarlo_Price(FileName); %要更新lamda_q_NT时再打开，不然浪费时间
% load 'lamda_c_q_NT';
% lamda_q_NT = lamda_c_q_NT;
load 'lamda_q_6N24T';
% load 'lamda_q_30N24T';

T = SCUC_data.totalLoad.T;  % 时段数T
G = SCUC_data.units.N;      % 发电机数
N = SCUC_data.baseparameters.busN;  % 节点总数

all_branch.I = [ SCUC_data.branch.I; SCUC_data.branchTransformer.I ]; %所有支路起点 前是支路起点 后是变压器支路起点
all_branch.J = [ SCUC_data.branch.J; SCUC_data.branchTransformer.J ]; %所有支路终点
all_branch.P = [ SCUC_data.branch.P; SCUC_data.branchTransformer.P ]; %支路功率上限

beta_CVaR = 0.99;

% 形成直流潮流系数矩阵B
type_of_pf = 'DC';
Y = SCUC_nodeY(SCUC_data,type_of_pf);
B = -Y.B; %因为是直流方程 所以B忽略了电阻 只考虑电抗

% 定义一些方便用的常量
diag_E_T = sparse(1:T,1:T,1); %T*T的对角阵，对角线全1
low_trianle = diag_E_T(2:T,:) - diag_E_T(1:T-1,:); 

% 定义变量y                               
for i=1:N    % 逐 个 处理每个节点   
    if ismember(i, SCUC_data.units.bus_G)    % 如果是发电机节点 决策变量包括 功率，相角，功率约束z
        y{i}.PG = sdpvar(T,1); %sdpvar创建实数型决策变量
        y{i}.theta = sdpvar(T,1);
    else
        y{i}.theta = sdpvar(T,1);
    end
end


cons = [];  % 所有的约束集合

% 按照节点开始建立各类约束
% 直流潮流不等式 每个循环构造T行不等式，总的构造了N*T个不等式
for i = 1:N  
    %构造中间项
    cons_PF = sparse(T,1); %一个t一个约束 T行约束
    if ismember(i, SCUC_data.units.bus_G) %如果i是发电机节点，则中间项应该加上PG
        cons_PF = cons_PF +  y{i}.PG;
    end
    for j = 1:N
        cons_PF = cons_PF -B(i,j) .* y{j}.theta;
%         if ismember(j, SCUC_data.units.bus_G) && j == i   % 如果是发电机节点 所在的对角子矩阵位置
%             cons_PF = cons_PF +  y{j}.PG;
%         end
    end
    
    %构造右侧项
    %find函数返还的是下标号
    index_loadNode = find(SCUC_data.busLoad.bus_PDQR==i); % 节点i是否为负荷节点
    if index_loadNode>0
        b_tmp = SCUC_data.busLoad.node_P(:,index_loadNode); %按下标取出该负荷节点负载
    else
        b_tmp = sparse(T,1);
    end
    cons = [cons, ...
        0 <= cons_PF <= b_tmp     ];
        
end
% end of 直流潮流不等式


% 参考节点
%每个模型都只设定第一个节点为参考节点？
cons = [cons, ...
    y{1}.theta == sparse(T,1)      ];  
% end of 参考节点


% 处理发电机节点约束
for i = 1:G
    bus_index_Gen = SCUC_data.units.bus_G(i);
    % 机组出力上下界
    cons = [cons, ...
        SCUC_data.units.PG_low(i) * ones(T,1) <= y{bus_index_Gen}.PG <= SCUC_data.units.PG_up(i) * ones(T,1)     ];
    
    % 爬坡 Pup和Pdown相等 且不考虑P0
    cons = [cons, ...
        -SCUC_data.units.ramp(i) * ones(T-1,1) <= low_trianle * y{bus_index_Gen}.PG <= SCUC_data.units.ramp(i) * ones(T-1,1)  ];
    
%     % 发电费用线性化
%     for l = 0:(L-1)
%         p_i_l =SCUC_data.units.PG_low(i) +  ( SCUC_data.units.PG_up(i) - SCUC_data.units.PG_low(i) ) / L * l;
%         cons = [cons, ...
%             (2* p_i_l * SCUC_data.units.gamma(i) +SCUC_data.units.beta(i) ) * y{bus_index_Gen}.PG - y{bus_index_Gen}.z <= (p_i_l^2 * SCUC_data.units.gamma(i) - SCUC_data.units.alpha(i)) * ones(T,1)  ];
%     end
end
% end of 处理发电机节点约束


% 按照支路构造各类约束
%线路潮流约束
for i = 1:size(all_branch.I,1) %从第一条支路开始循环遍历所有支路
    left = all_branch.I(i); %支路起点和终点即可得到B电纳
    right = all_branch.J(i);
    abs_x4branch = abs(1/B(left,right));  % 当前支路的阻抗的绝对值 |x_ij|    
    cons = [ cons, ...
      -all_branch.P(i) * abs_x4branch * ones(T,1) <= y{left}.theta - y{right}.theta <= all_branch.P(i) * abs_x4branch * ones(T,1) ];  
end
% end of 按照支路构造各类约束

% 下面计算 miu_hat均值的估计值 sigema_hat协方差的估计值
miu_hat_G_T = zeros(G,T); %每个机组每个时段电价都不同，G*T个
for q = 1:q_line % q_line标记样本长度
    miu_hat_G_T = miu_hat_G_T + reshape(lamda_q_NT(q,:,:),G,T); %按样本数q依次累加 G*T矩阵型的电价
end
miu_hat_G_T = 1/q_line * miu_hat_G_T;
miu_hat = reshape(miu_hat_G_T',G*T,1); %按每个机组所有时段排列为列向量

%定义变量
t = sdpvar(1);
uk = sdpvar(q_line,T);

Constraints = [];
Constraints = [Constraints,cons]; 

% P_vector = []; %形如[p11 p12 ...p1t;p21...p2t...pnt]的列向量
% cost = 0;
% for g = 1:G
%     P_vector = [P_vector;  y{SCUC_data.units.bus_G(g)}.PG  ]; %列向量用;分隔
%     cost = cost + SCUC_data.units.alpha(g) + SCUC_data.units.beta(g)* y{SCUC_data.units.bus_G(g)}.PG + (SCUC_data.units.gamma(g) * y{SCUC_data.units.bus_G(g)}.PG * y{SCUC_data.units.bus_G(g)}.PG);
% end

% for q = 1:q_line
%     Constraints = [Constraints, uk(q) <= 0 ];
%     lamda_k = reshape(lamda_q_NT(q,:,1),G,1);
%     Constraints = [Constraints, uk(q) <= lamda_k' * P_vector - cost -t]; 
% end

for h = 1:T
    P_vector = []; 
    cost = 0;
    for g = 1:G
        P_vector = [P_vector;  y{SCUC_data.units.bus_G(g)}.PG(h)  ]; %列向量用;分隔
        cost = cost + SCUC_data.units.alpha(g) + SCUC_data.units.beta(g)* y{SCUC_data.units.bus_G(g)}.PG(h) + (SCUC_data.units.gamma(g) * y{SCUC_data.units.bus_G(g)}.PG(h) * y{SCUC_data.units.bus_G(g)}.PG(h));
    end
    
    for q = 1:q_line
        Constraints = [Constraints, uk(q,h) <= 0 ];
        lamda_k = reshape(lamda_q_NT(q,:,h),G,1);
        Constraints = [Constraints, uk(q,h) <= lamda_k' * P_vector - cost -t]; 
    end
end

% W = 0;
% for g = 1:G
%     W = W + (SCUC_data.units.gamma(g) * y{SCUC_data.units.bus_G(g)}.PG * y{SCUC_data.units.bus_G(g)}.PG);
% end
% Constraints = [Constraints, p >= W ]; 

Objective = t + sum(sum(uk)) / ((1 - beta_CVaR) * (q_line));

options = sdpsettings('verbose',1,'solver','mosek','debug',1,'savesolveroutput',1,'savesolverinput',1);
% options.cplex.exportmodel='xyz.lp';

sol = optimize(Constraints,-Objective,options);

% Analyze error flags
if sol.problem == 0
    % Extract and display value
    Obj = value(Objective);

    
%     %输出调度方案
%     for i=1:N  
%         if ismember(i, SCUC_data.units.bus_G) 
%             disp(i);
%             disp(value(y{i}.PG));
%             disp(value(y{i}.theta));
%         else
%             disp(i);
%             disp(value(y{i}.theta));
%         end
%     end

  %输出发电费用
    P_vector_re = []; 
    all_cost = 0;
    for g = 1:G
        P_vector_re = [P_vector_re;  y{SCUC_data.units.bus_G(g)}.PG  ]; %列向量用;分隔
        for h = 1:T
            all_cost = all_cost + SCUC_data.units.alpha(g) + SCUC_data.units.beta(g)* y{SCUC_data.units.bus_G(g)}.PG(h) + (SCUC_data.units.gamma(g) * y{SCUC_data.units.bus_G(g)}.PG(h) * y{SCUC_data.units.bus_G(g)}.PG(h));
        end
    end
    profit = miu_hat' * P_vector_re - all_cost;
    disp(value(profit));
    disp(Obj);
    disp(sol.solvertime);
else
    display('Oh shit!, something was wrong!');
    disp(value(t));
    sol.info
    yalmiperror(sol.problem)
end




