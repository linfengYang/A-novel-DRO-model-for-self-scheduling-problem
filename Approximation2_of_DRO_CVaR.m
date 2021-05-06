%%  DRO_CVaR_ADMM   by yy 2020.12.22

% function [k,time,epsilon]=DRO_CVaR_Alg1(pathAndFilename)

% FileName = 'OPF_DAT/IEEE4.dat';
FileName = 'SCUC_dat/SCUC6_4period.txt';
% FileName = 'SCUC_dat/SCUC30.txt';
SCUC_data = ReadDataSCUC(FileName);
L = 4;

% MonteCarlo_Price(SCUC_data)
load 'lamda_q_6N4Tterminal';

T = SCUC_data.totalLoad.T;  % 时段数T
G = SCUC_data.units.N;      % 发电机数
N = SCUC_data.baseparameters.busN;  % 节点总数

all_branch.I = [ SCUC_data.branch.I; SCUC_data.branchTransformer.I ]; %所有支路起点 前是支路起点 后是变压器支路起点
all_branch.J = [ SCUC_data.branch.J; SCUC_data.branchTransformer.J ]; %所有支路终点
all_branch.P = [ SCUC_data.branch.P; SCUC_data.branchTransformer.P ]; %支路功率上限

beta_CVaR = 0.95;

% 形成直流潮流系数矩阵B
type_of_pf = 'DC';
Y = SCUC_nodeY(SCUC_data,type_of_pf);
B = -Y.B; %因为是直流方程 所以B忽略了电阻 只考虑电抗

% 定义一些方便用的常量
diag_E_T = sparse(1:T,1:T,1); %T*T的对角阵，对角线全1
low_trianle = diag_E_T(2:T,:) - diag_E_T(1:T-1,:); 

% 定义变量 y 和 z                            
for i=1:N    % 逐 个 处理每个节点   
    if ismember(i, SCUC_data.units.bus_G)    % 如果是发电机节点 决策变量包括 功率，相角，功率约束z
        y{i}.PG = sdpvar(T,1); %sdpvar创建实数型决策变量
        y{i}.theta = sdpvar(T,1);
        y{i}.z = sdpvar(T,1);
    else
        y{i}.theta = sdpvar(T,1);
    end
end

for i=1:N    % 逐 个 处理每个节点   
    if ismember(i, SCUC_data.units.bus_G)    % 如果是发电机节点 决策变量包括 功率，相角，功率约束z
        z{i}.PG = sdpvar(T,1); %sdpvar创建实数型决策变量
        z{i}.theta = sdpvar(T,1);
        z{i}.z = sdpvar(T,1);
    else
        z{i}.theta = sdpvar(T,1);
    end
end
% end of 定义变量

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

%调用切分函数切分区域
[PI,PINumber,PIG] = portion(N);
D = length(PINumber); %确定划分块数

%区分各节点是内节点还是外节点
ext = [];
for d = 1:D
    for i = 1:size(all_branch.I,1) %从第一条支路开始循环遍历所有支路
        left = all_branch.I(i); %支路起点和终点即可得到B电纳
        right = all_branch.J(i);
        if ismember(left, PI{d}') && ~ismember(right, PI{d}') %判断该支路是否属于当前分区
            ext = [ext;left;right];
        end
    end
end
ext = unique(ext);
%定义其他变量
r = sdpvar(D,1);
t = sdpvar(D,1);
alpha_CVaR = sdpvar(D,1);
z_vector = [];
P_vector = [];
for g = 1:G
    z_vector = [z_vector;  y{SCUC_data.units.bus_G(g)}.z  ];
    P_vector = [P_vector;  y{SCUC_data.units.bus_G(g)}.PG ]; %列向量用;分隔
end

%设置ADMM参数
u_A_b = [];
rho = 5;
gap_pri = 1e-2;
gap_dual = 1e-2;

M = 100; %迭代次数
% core = D;
% p = parpool(core);
% p.IdleTimeout = 100;
for k = 1:M
    p_A = [];
    z_A = [];
    z_A_pri = [];
    u_A = [];
    p_A_k1 = [];
    u_r = 0; %标记分区变量长度
    %p-update
    for d = 1:D   
        PIi = PI{d};
        PIN = PINumber{d};
        d_g = PIG{d};
        Constraints = []; %每次将约束集置空
        miu_hat_d = [];
        lamda_positive_d = [];
        lamda_negative_d = [];
        y_d = [];
        zy_d = [];
        theta_d = [];
        ztheta_d = [];
        z_d = [];
        for i = PIi' %循环处理当前片区节点
            theta_d = [theta_d;y{i}.theta]; %取当前分区所有节点相角
            ztheta_d = [ztheta_d;z{i}.theta]; %取当前分区所有节点相角
            
            if ismember(i, SCUC_data.units.bus_G)
                y_d = [y_d;y{i}.PG]; %取当前片区所有电机发电量
                zy_d = [zy_d;z{i}.PG]; %取当前片区所有电机发电量
                z_d = [z_d;y{i}.z]; %取当前片区电机费用
                i_g = find(SCUC_data.units.bus_G == i); %求当前电机是第几个电机
                % 机组出力上下界
                Constraints = [Constraints, ...
                    SCUC_data.units.PG_low(i_g) * ones(T,1) <= y{i}.PG <= SCUC_data.units.PG_up(i_g) * ones(T,1)     ];
                
                % 爬坡 Pup和Pdown相等 且不考虑P0
                Constraints = [Constraints, ...
                    -SCUC_data.units.ramp(i_g) * ones(T-1,1) <= low_trianle * y{i}.PG <= SCUC_data.units.ramp(i_g) * ones(T-1,1)  ];
                
                % 发电费用线性化
                for l = 0:(L-1)
                    p_i_l =SCUC_data.units.PG_low(i_g) +  ( SCUC_data.units.PG_up(i_g) - SCUC_data.units.PG_low(i_g) ) / L * l;
                    Constraints = [Constraints, ...
                        (2* p_i_l * SCUC_data.units.gamma(i_g) +SCUC_data.units.beta(i_g) ) * y{i}.PG - y{i}.z <= (p_i_l^2 * SCUC_data.units.gamma(i_g) - SCUC_data.units.alpha(i_g)) * ones(T,1)  ];
                end  
                
                %计算分区电价均值
                miu_hat_m = zeros(1,T);
                for q = 1:q_line % q_line标记样本长度
                    miu_hat_m = miu_hat_m + reshape(lamda_q_NT(q,i_g,:),1,T); %按样本数q依次累加第i个机组的电价
                end
                miu_hat_m = 1/q_line * miu_hat_m;
                miu_hat_m = reshape(miu_hat_m',T,1); %按每个机组所有时段排列为列向量
                miu_hat_d = [miu_hat_d ; miu_hat_m];
                
                %取对应电机最大最小值
                lamda_positive_m = lamda_positive(1+T*(i_g-1):T*i_g,1);
                lamda_positive_d = [lamda_positive_d;lamda_positive_m];
                lamda_negative_m = lamda_negative(1+T*(i_g-1):T*i_g,1);
                lamda_negative_d = [lamda_negative_d;lamda_negative_m];
            end 
            
            if ~ismember(i, ext) %判断节点是否是内部节点
                %构造分区内部的直流潮流不等式
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
                Constraints = [Constraints, ...
                    0 <= cons_PF <= b_tmp     ];
                %end of 构造分区内部的支流潮流不等式         
            end
            
            
        end
        
        %分区内部线路功率约束
        for i = 1:size(all_branch.I,1) %从第一条支路开始循环遍历所有支路
            left = all_branch.I(i); %支路起点和终点即可得到B电纳
            right = all_branch.J(i);
            if ismember(left,PIi') && ismember(right,PIi') %判断该支路是否属于当前分区
                abs_x4branch = abs(1/B(left,right));  % 当前支路的阻抗的绝对值 |x_ij|
                Constraints = [ Constraints, ...
                    -all_branch.P(i) * abs_x4branch * ones(T,1) <= y{left}.theta - y{right}.theta <= all_branch.P(i) * abs_x4branch * ones(T,1) ];
            end 
        end
        % end of 分区内部线路功率约束
        
        %计算分区电价协方差
        Sigma_hat_d = zeros(d_g*T,d_g*T);
        for q = 1:q_line
            tmp2 = [];
            for i = PIi'
                if ismember(i, SCUC_data.units.bus_G)
                    i_g = find(SCUC_data.units.bus_G == i); 
                    tmp1 = reshape(lamda_q_NT(q,i_g,:),1*T,1);
                    tmp2 = [tmp2;tmp1];
                end
            end
            tmp2 = tmp2 - miu_hat_d;
            Sigma_hat_d = Sigma_hat_d + tmp2  * tmp2';
        end
        Sigma_hat_d = 1/q_line * Sigma_hat_d;
        Sigma_hat_neg_half = Sigma_hat_d^(-1/2);
     
        %取分块变量
        A_d = [sparse(1:T*d_g,1:T*d_g,1); -sparse(1:T*d_g,1:T*d_g,1)];
        B_d = [lamda_positive_d; -lamda_negative_d ];
        
        % 构造分区模糊集
        qline = 2000;
        deta = 0.3; %δ在0，1之间任取
        deta_bar = 1 - sqrt(1-deta);
        part1 = abs(Sigma_hat_neg_half * (lamda_positive_d - miu_hat_d ));
        part2 = abs(Sigma_hat_neg_half * (lamda_negative_d - miu_hat_d ));
        part1 = norm(part1);
        part2 = norm(part2);
        R_hat = max(part2, part1);
        
        yy1 = (1 -  ( (R_hat^2 + 2) * (2+sqrt(2*log(4/deta_bar))) / sqrt(qline) )   )^(-1/2);
        R_bar = R_hat * yy1 ;
        
        a_hua_bar = (R_bar^2/sqrt(qline)) * (sqrt(1-G*T/R_bar^4) + sqrt(log(4/deta_bar)));
        b_hua_bar = (R_bar^2/sqrt(qline)) * ( 2+sqrt(2*log(2/deta_bar)) )^2;
        M_hat = max( (R_hat^2 + 2)^2 * (2 + sqrt(2* log(4/deta_bar)))^2  , (8+sqrt(32*log(4/deta_bar)))^2 / (sqrt(R_hat+4) - R_hat)^4 );
        
        gamma_bar_d_1 = b_hua_bar/(1- a_hua_bar - b_hua_bar);
        gamma_bar_d_2 = (1+b_hua_bar)/(1- a_hua_bar - b_hua_bar);
        % end of 构造gama

        r_d = r(d);
        t_d = t(d);
        alpha_CVaR_d = alpha_CVaR(d);
        Q_d = sdpvar(d_g*T, d_g*T);
        q_d = sdpvar(d_g*T,1);
        tao1_d = sdpvar(2*d_g*T, 1);
        tao2_d = sdpvar(2*d_g*T, 1);
        p_A_d = [ y_d ; theta_d];
        p_A = [p_A;value(p_A_d)]; %按片区排列pk
%         z_A_d = [r_d ; r_d; zy_d ; ztheta_d];
%         z_A = [z_A;(z_A_d)];
%         u_A_d = u_A_0(1 + u_r:(2+d_g*T+PIN*T) + u_r);
%         u_r = u_r + (2+d_g*T+PIN*T);
%         u_A = [u_A;(u_A_d)];
        if k == 1
            z_A_d = zeros(d_g*T+PIN*T,1);
            z_A = [z_A;[ zy_d ; ztheta_d]];
            u_A_d = zeros(d_g*T+PIN*T,1);
            u_A = [u_A;(u_A_d)];
        else
            z_A_d = [ zy_d ; ztheta_d];
            z_A = [z_A;(z_A_d)];
            if d ==1
                u_A_d = u_A_b(1 : (d_g*T+PIN*T));
            elseif d == 2
                 u_A_d = u_A_b(1 +(PIG{1}*T+PINumber{1}*T): (d_g*T+PIN*T)+(PIG{1}*T+PINumber{1}*T));
            else
                u_A_d = u_A_b(1 +((PIG{1}+PIG{2})*T+(PINumber{1}+PINumber{2})*T): (d_g*T+PIN*T)+((PIG{1}+PIG{2})*T+(PINumber{1}+PINumber{2})*T));
            end
            u_A = [u_A;(u_A_d)];
        end
        z_A_pri = [z_A_pri;value(z_A_d)];
        
        % 参考节点
        Constraints = [Constraints, ...
            y{1}.theta == sparse(T,1)      ];
        % end of 参考节点
        
        Constraints = [Constraints, tao1_d>=0, tao2_d>=0];  
        Constraints = [Constraints, Q_d>=0];
        Constraints = [Constraints, t_d>= sum(sum(   (gamma_bar_d_2 * Sigma_hat_d + miu_hat_d * miu_hat_d') .* Q_d   )) + miu_hat_d' * q_d + sqrt(gamma_bar_d_1) * norm(Sigma_hat_d^(1/2) * (q_d + 2 * Q_d * miu_hat_d))];
        Constraints = [Constraints, [Q_d, 1/2 * q_d;  1/2 *q_d' , r_d-alpha_CVaR_d] >= -1/2 * [sparse(d_g*T,d_g*T), A_d' * tao1_d;  tao1_d' * A_d,  -2 * tao1_d' * B_d]];
        
        z_d = sum(z_d);
        P_vector_d = y_d;
        Constraints = [Constraints, [Q_d, 1/2 * q_d;  1/2 *q_d' , r_d] + 1 / (1-beta_CVaR) * [sparse(d_g*T,d_g*T), 1/2 * P_vector_d ;  1/2 * P_vector_d' ,beta_CVaR * alpha_CVaR_d - z_d  ] >= -1/2 * [sparse(d_g*T,d_g*T), A_d'* tao2_d;  tao2_d' * A_d,  -2 * tao2_d' * B_d]];
        
        Objective =  r_d + t_d +(rho / 2) * norm(p_A_d - value(z_A_d) + u_A_d)^2;
        options = sdpsettings('verbose',0,'solver','mosek','debug',1);
        sol = optimize(Constraints,Objective,options);
        p_A_k1 = [p_A_k1;p_A_d]; %按片区排列pk+1
        % Analyze error flags
        if sol.problem == 0
            % Extract and display value
%             Obj = value(Objective);
%             disp(Obj);
%             disp(sol.solvertime);
        else
            disp('Oh shit!, something was wrong!');
            sol.info
            yalmiperror(sol.problem)
        end
    end
    
    %z-update
    Constraints = [];
    
%     % 参考节点
%     Constraints = [Constraints, ...
%         z{1}.theta == sparse(T,1)      ];
%     % end of 参考节点

    % 分区之间约束条件
    for d = 1:D %按各个分区依次构造
        for i = PI{d}' % 分区之间的直流潮流约束
            if ismember(i,ext)
                cons_PF = sparse(T,1);
                if ismember(i, SCUC_data.units.bus_G) %如果i是发电机节点，则中间项应该加上PG
                    cons_PF = cons_PF +  z{i}.PG;
                end
                for j = 1:N
                    cons_PF = cons_PF -B(i,j) .* z{j}.theta;
                end
                
                %构造右侧项
                index_loadNode = find(SCUC_data.busLoad.bus_PDQR==i); % 节点i是否为负荷节点
                if index_loadNode>0
                    b_tmp = SCUC_data.busLoad.node_P(:,index_loadNode); %按下标取出该负荷节点负载
                else
                    b_tmp = sparse(T,1);
                end
                Constraints = [Constraints, ...
                    0 <= cons_PF <= b_tmp     ]; 
            end
        end
        % end of 分区之间的直流潮流约束
        
        %分区之间线路功率约束
        for s = 1:size(all_branch.I,1) %从第一条支路开始循环遍历所有支路
            left = all_branch.I(s); %支路起点和终点即可得到B电纳
            right = all_branch.J(s);
            if ismember(left, PI{d}') && ~ismember(right, PI{d}') %判断该支路是否属于分区间
                abs_x4branch = abs(1/B(left,right));  % 当前支路的阻抗的绝对值 |x_ij|
                Constraints = [ Constraints, ...
                    -all_branch.P(s) * abs_x4branch * ones(T,1) <= z{left}.theta - z{right}.theta <= all_branch.P(s) * abs_x4branch * ones(T,1) ];
            end
        end
        % end of 分区之间线路功率约束   
    
    end
    

    Objective = (rho / 2) * norm( z_A - value(p_A_k1) - u_A)^2;
    options = sdpsettings('verbose',0,'solver','mosek','debug',1);
    sol = optimize(Constraints,Objective,options);
    % Analyze error flags
    if sol.problem == 0
        % Extract and display value
%         Obj = value(Objective);
%         disp(Obj);
%         disp(sol.solvertime);
    else
        disp('Oh shit!, something was wrong!');
        sol.info
        yalmiperror(sol.problem)
    end
    
    %u-update
    u_A = u_A + (value(p_A_k1) - value(z_A));
    u_A_b = u_A;
    
    %判断精度
    pri = p_A - z_A_pri;
    pri = norm(pri);
    dual = rho * ( z_A_pri - z_A );
    dual = norm(value(dual));
    
    if ( pri <= gap_pri) && ( dual <= gap_dual)
        disp('ADMM over');
        break
    end
    msg = sprintf('第 %d 次循环. pri_gap : %d, dual_gap : %d',k,pri,dual);
    disp(msg);
    
end

%输出目标函数
obj = r + t;
obj = sum(obj);
disp(value(obj));

%输出发电费用
cost = P_vector' * miu_hat - sum(z_vector);
disp(value(cost));
