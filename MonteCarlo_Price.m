function MonteCarlo_Price(FileName)
% FileName = 'SCUC_dat/SCUC6_4period.txt';
SCUC_data = ReadDataSCUC(FileName);
L = 4;

% G = size(SCUC_data.busLoad.bus_PDQR,1);
% T = size(SCUC_data.busLoad.node_P,1);
% N = SCUC_data.baseparameters.busN;

T = SCUC_data.totalLoad.T;  % 时段数T
G = SCUC_data.units.N;      % 发电机数
N = SCUC_data.baseparameters.busN;  % 节点总数
all_branch.I = [ SCUC_data.branch.I; SCUC_data.branchTransformer.I ]; %所有支路起点 前是支路起点 后是变压器支路起点
all_branch.J = [ SCUC_data.branch.J; SCUC_data.branchTransformer.J ]; %所有支路终点
all_branch.P = [ SCUC_data.branch.P; SCUC_data.branchTransformer.P ]; %支路功率上限

type_of_pf = 'DC';
Y = SCUC_nodeY(SCUC_data,type_of_pf);
B = -Y.B;

% 定义一些方便用的常量
diag_E_T = sparse(1:T,1:T,1); %T*T的对角阵，对角线全1
low_trianle = diag_E_T(2:T,:) - diag_E_T(1:T-1,:); 

cons = [];  % 所有的约束集合
for i=1:N    % 逐 个 处理每个节点   
    if ismember(i, SCUC_data.units.bus_G)    % 如果是发电机节点 决策变量包括 功率，相角，功率约束z
        y{i}.PG = sdpvar(T,1); %sdpvar创建实数型决策变量
        y{i}.theta = sdpvar(T,1);
        y{i}.z = sdpvar(T,1);
    else
        y{i}.theta = sdpvar(T,1);
    end
end

% 直流潮流不等式 每个循环构造T行不等式，总的构造了N*T个不等式
for i = 1:N  
    %构造中间项
    cons_PF = sparse(T,1); %一个t一个约束 T行约束
    if ismember(i, SCUC_data.units.bus_G) %如果i是发电机节点，则中间项应该加上PG
        cons_PF = cons_PF +  y{i}.PG;
    end
    for j = 1:N
        cons_PF = cons_PF -B(i,j) .* y{j}.theta;
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
      cons_PF - b_tmp == sparse(T,1)   ];
        
end

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

% 参考节点
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
end
% end of 处理发电机节点约束

Constraints = [];
Constraints = [Constraints,cons];

f_G = 0;
g=0;
for i = 1:N
    if ismember(i, SCUC_data.units.bus_G)
        g = g+1;
        for t = 1:T
            f_G = f_G + (SCUC_data.units.alpha(g) + SCUC_data.units.beta(g)*y{i}.PG(t) + SCUC_data.units.gamma(g)*y{i}.PG(t)^2 );
        end
    end
end

Objective =  f_G;

% options = sdpsettings('verbose',1,'solver','sdpt3','debug',1,'sdpt3.gaptol',1e-5);
options = sdpsettings('solver','cplex','debug',1,'savesolveroutput',1,'savesolverinput',1);
options.verbose = 1;%设置显示信息程度，1为适度显示，2为完全显示
options.cplex.exportmodel='abc.lp';

sol = optimize(Constraints,Objective,options);


% % Analyze error flags
% if sol.problem == 0
%     % Extract and display value
%     Obj = value(Objective);
%     disp(Obj);
% else
%     disp('Oh shit!, something was wrong!');
%     sol.info
%     yalmiperror(sol.problem)
% end


%开始扰动lamda
q_line = 168;

miu_T_0 = -1*sol.solveroutput.lambda.eqlin(1:N*T);
miu_T_0 = reshape(miu_T_0,N,T);
miu_T_1 = zeros(G,T);
j = 1;
for i = 1:N
    if ismember(i, SCUC_data.units.bus_G)
        miu_T_1(j,:) = miu_T_0(i,:);
        j = j+1;
    end
end
% miu_T_0 = 1:T;
% miu_T_1 = miu_T_1 / 10;
miu_T_1 = abs(miu_T_1);
disp(miu_T_1);
PD = zeros(1,T);
std_div_T_0 = zeros(1,T);
% std_div_T_0 = (PD/100-6);
% std_div_T_0(8) = std_div_T_0(8) +9;
% std_div_T_0(9) = std_div_T_0(9) +10;
% std_div_T_0(10) = std_div_T_0(10) +9;
% std_div_T_0(18) = std_div_T_0(18) +7;
% std_div_T_0(19) = std_div_T_0(19) +8;
% std_div_T_0(20) = std_div_T_0(20) +8;
miu_NT = zeros(G,T);
miu_NT_2 = zeros(G,T);
std_div_NT = zeros(G,T);


for i = 1:G
%     miu_NT(i,:) =  miu_T_1(i,:) - SCUC_data.busLoad.node_P(:,i)' * 2; %数据的均值
%     miu_NT_2(i,:) =  miu_T_1(i,:) - SCUC_data.busLoad.node_P(:,i)' * 5;
    miu_NT(i,:) =  miu_T_1(i,:) + 4.24; %数据的均值
    miu_NT_2(i,:) =  miu_T_1(i,:) - 0.235;
%     std_div_NT(i,:) =  std_div_T_0 * 0 +  1  +  randn*0; %数据的方差
end


lamda_q_NT = zeros(q_line,G,T);
for i = 1:G
    for t = 1:T
%           lamda_q_NT(:,i,t) = exprnd(miu_T_1(i,t),[q_line,1]); 
        lamda_q_NT(:,i,t) = unifrnd(miu_NT_2(i,t),miu_NT(i,t),[q_line,1]);  %均匀分布
%         lamda_q_NT(:,i,t) = random('norm',miu_NT(i,t),std_div_NT(i,t),[q_line,1]); %正态分布
%         lamda_q_NT(:,t) = price_miu(t) * ones(q_line,1);
    end
end
% for i =1:G
%     figure(i);
%     tmp = reshape(lamda_q_NT(:,i,:),q_line,T);
%     plot(tmp');
%     hold on
% end

save('lamda_q_118N24Ttest.mat','lamda_q_NT','q_line');