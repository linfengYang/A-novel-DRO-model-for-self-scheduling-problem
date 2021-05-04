%%  DRO_CVaR   by ylf 2020.7.2

% function [k,time,epsilon]=DRO_CVaR_Alg1(pathAndFilename)

% FileName = 'OPF_DAT/IEEE4.dat';
% FileName = 'SCUC_dat/SCUC6_4period.txt';
FileName = 'SCUC_dat/SCUC6.txt';
SCUC_data = ReadDataSCUC(FileName);
L = 4; %给定的参数，用于线性化发电费用

% MonteCarlo_Price(FileName); %要更新lamda_q_NT时再打开，不然浪费时间
% load 'lamda_c_q_NT';
% lamda_q_NT = lamda_c_q_NT;
load 'lamda_q_6N24T';

T = SCUC_data.totalLoad.T;  % 时段数T
G = SCUC_data.units.N;      % 发电机数
N = SCUC_data.baseparameters.busN;  % 节点总数

all_branch.I = [ SCUC_data.branch.I; SCUC_data.branchTransformer.I ]; %所有支路起点 前是支路起点 后是变压器支路起点
all_branch.J = [ SCUC_data.branch.J; SCUC_data.branchTransformer.J ]; %所有支路终点
all_branch.P = [ SCUC_data.branch.P; SCUC_data.branchTransformer.P ]; %支路功率上限

beta_CVaR = 0.90;

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
        y{i}.z = sdpvar(T,1);
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
        disp(value(b_tmp));
    else
        b_tmp = sparse(T,1);
    end
    cons = [cons, ...
        0 <= cons_PF <= b_tmp     ];
        
end
% end of 直流潮流不等式  

i = 2;
cons_PF = zeros(T,1); %一个t一个约束 T行约束
for j = 1:N
    cons_PF = cons_PF -B(i,j) .* y{j}.theta;

end

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
    
    % 发电费用线性化
    for l = 0:(L-1)
        p_i_l =SCUC_data.units.PG_low(i) +  ( SCUC_data.units.PG_up(i) - SCUC_data.units.PG_low(i) ) / L * l;
        cons = [cons, ...
            (2* p_i_l * SCUC_data.units.gamma(i) +SCUC_data.units.beta(i) ) * y{bus_index_Gen}.PG - y{bus_index_Gen}.z <= (p_i_l^2 * SCUC_data.units.gamma(i) - SCUC_data.units.alpha(i)) * ones(T,1)  ];
    end
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

% Sigma_hat = sparse(1:G*T,1:G*T,1);
% Sigma_hat = full(Sigma_hat);
Sigma_hat = zeros(G*T,G*T);
for q = 1:q_line
    tmp1 = reshape(lamda_q_NT(q,:,:),G,T) - miu_hat_G_T;
    tmp2 = reshape(tmp1', G*T,1);
    Sigma_hat = Sigma_hat + tmp2  * tmp2';
end
Sigma_hat = 1/q_line * Sigma_hat;
Sigma_hat_neg_half = Sigma_hat^(-1/2);
% end of 均值 协方差计算


% 构造模糊集S G*T维的列向量，按每个G所有时段排列
tmp = max(lamda_q_NT,[],1); %取各个时段，所有样本中G个机组的最大值
% tmp = max(lamda_q_NT(:,:,1),[],1); %取第一时段，所有样本G个机组的最大值
tmp = reshape(tmp,G,T)';
lamda_positive = reshape(tmp, G*T,1);  % 计算 lamda正
tmp = min(lamda_q_NT,[],1);
% tmp = min(lamda_q_NT(:,:,1),[],1);
tmp = reshape(tmp,G,T)';
lamda_negative = reshape(tmp, G*T,1);  % 计算 lamda负
% end of 构造模糊集S


%构造gama
% qline = 200000000; %30
qline = 10000000;
deta = 0.2; %δ在0，1之间任取
deta_bar = 1 - sqrt(1-deta);

%  R ?=max┬(λ∈S)?‖Σ ?^(-1?2) (λ-μ ? )‖
part1 = abs(Sigma_hat_neg_half * (lamda_positive - miu_hat ));
part2 = abs(Sigma_hat_neg_half * (lamda_negative - miu_hat ));
part1 = norm(part1);
part2 = norm(part2);
% R_hat = max(part2, part1)/10; & 30-24
R_hat = max(part2, part1)-10;

yy1 = (1 -  ( (R_hat^2 + 2) * (2+sqrt(2*log(4/deta_bar))) / sqrt(qline) )   )^(-1/2);
R_bar = R_hat * yy1 ;

a_hua_bar = (R_bar^2/sqrt(qline)) * (sqrt(1-G*T/R_bar^4) + sqrt(log(4/deta_bar)));
b_hua_bar = (R_bar^2/sqrt(qline)) * ( 2+sqrt(2*log(2/deta_bar)) )^2;
M_hat = max( (R_hat^2 + 2)^2 * (2 + sqrt(2* log(4/deta_bar)))^2  , (8+sqrt(32*log(4/deta_bar)))^2 / (sqrt(R_hat+4) - R_hat)^4 );

gamma_bar_1 = b_hua_bar/(1- a_hua_bar - b_hua_bar);                     
gamma_bar_2 = (1+b_hua_bar)/(1- a_hua_bar - b_hua_bar);
% end of 构造gama

A = [sparse(1:G*T,1:G*T,1); -sparse(1:G*T,1:G*T,1)];
B_set = [lamda_positive; -lamda_negative ];

% 下面开始形成第一个SDP模型，其中变量排列为：y(包括P,theta, z), alpha, Q,r,t,tao1,tao2
% index_P = [];
% for i = 1:N
%     index_P = [index_P, (i-1) * T +1 :i * T];
% end
% y = sdpvar(2*G*T,1);

r = sdpvar(1);
t = sdpvar(1);
alpha_CVaR = sdpvar(1);
Q = sdpvar(G*T, G*T);
q = sdpvar(G*T,1);
tao1 = sdpvar(2*G*T, 1);
tao2 = sdpvar(2*G*T, 1);
Constraints = [];
Constraints = [Constraints, tao1>=0, tao2>=0];  
Constraints = [Constraints,cons];  % y(- Y
Constraints = [Constraints, Q>=0];

% %.*表示元素对位相乘

Constraints = [Constraints, t>= sum(sum(   (gamma_bar_2 * Sigma_hat + miu_hat * miu_hat') .* Q   )) + miu_hat' * q + sqrt(gamma_bar_1) * norm(Sigma_hat^(1/2) * (q + 2 * Q * miu_hat))];

Constraints = [Constraints, [Q, 1/2 * q;  1/2 *q' , r-alpha_CVaR] >= -1/2 * [sparse(G*T,G*T), A' * tao1;  tao1' * A,  -2 * tao1' * B_set]];

% Constraints = [Constraints, [Q, 1/2 * ( q + A' * tao1 );  1/2 * ( q + A' * tao1 )' , r - alpha_CVaR -  tao1' * B] >= 0];


sum_z = 0;
P_vector = []; %形如[p11 p12 ...p1t;p21...p2t...pnt]的列向量
for g = 1:G
    sum_z = sum_z +sum( y{SCUC_data.units.bus_G(g)}.z );
    P_vector = [P_vector;  y{SCUC_data.units.bus_G(g)}.PG  ]; %列向量用;分隔
end

Constraints = [Constraints, [Q, 1/2 * q;  1/2 *q' , r] + 1 / (1-beta_CVaR) * [sparse(G*T,G*T), 1/2 * P_vector ;  1/2 * P_vector' ,beta_CVaR * alpha_CVaR - sum_z  ] >= -1/2 * [sparse(G*T,G*T), A'* tao2;  tao2' * A,  -2 * tao2' * B_set]];

% Z2 = sdpvar(1);
% Z2 = r +  1 / (1-beta_CVaR) * (beta_CVaR*alpha_CVaR - sum_z ) -  tao2' * B ;
% Constraints = [Constraints, Z2 == r +  1 / (1-beta_CVaR) * (beta_CVaR*alpha_CVaR - sum_z ) -  tao2' * B ;];
% Constraints = [Constraints, [Q, 1/2 * q +  1/2 * (1 / (1-beta_CVaR) * P_vector) + 1/2 * A' * tao2; 1/2 *q' + 1/2 * (1 / (1-beta_CVaR) * P_vector') + 1/2 * tao2' * A  ,  Z2] >= 0];

% Constraints = [Constraints, [Q, 1/2 * q +  1/2 * (1 / (1-beta_CVaR) * P_vector) + 1/2 * A' * tao2; 1/2 *q' + 1/2 * (1 / (1-beta_CVaR) * P_vector') + 1/2 * tao2' * A  ,  r +  1 / (1-beta_CVaR) * (beta_CVaR*alpha_CVaR - sum_z ) -  tao2' * B] >= 0];
% Constraints = [Constraints,r + (beta_CVaR*alpha_CVaR - sum_z ) * 1 / (1-beta_CVaR) -  tao2' * B + Z2 >= 0];


Objective =  r + t;
% Objective = sum_z;

% options = sdpsettings('verbose',1,'debug',1,'solver','cplex','savesolveroutput',1,'savesolverinput',1);
% options.cplex.exportmodel='abc.lp';%这句就是用来输出线性模型文件，保存在根目录下。
% options = sdpsettings('verbose',1,'solver','sdpt3','debug',1,'sdpt3.gaptol',1e-4,'savesolveroutput',1,'savesolverinput',1);
options = sdpsettings('verbose',2,'solver','mosek','debug',1,'savesolveroutput',1,'savesolverinput',1);

% options = sdpsettings('verbose',1,'solver','sdpt3');
sol = optimize(Constraints,Objective,options);


% Analyze error flags
if sol.problem == 0
    % Extract and display value
    Obj = value(Objective);
    
    %输出调度方案
    for i=1:N  
        if ismember(i, SCUC_data.units.bus_G) 
            disp(i);
            disp(value(y{i}.PG));
%             disp(value(y{i}.theta));
%             disp(value(y{i}.z));
%         else
%             disp(i);
%             disp(value(y{i}.theta));
        end
    end
    

    %输出发电费用
    cost = P_vector' * miu_hat - sum_z;
    disp(Obj);
    disp(value(cost));
    
    disp(sol.solvertime);
else
    disp('Oh shit!, something was wrong!');
    disp(value(r));
    disp(value(t));
    sol.info
    yalmiperror(sol.problem)
end

