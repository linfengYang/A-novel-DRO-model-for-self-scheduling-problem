%% ***************************************************************
%             作者：ylf
%             原创日期：2017年12月7号
%             修改日期： 
%             函数说明：从SCUC数据文件读取各相关数据
%% ***************************************************************
%%%  读文件取数据  
function SCUC_data=ReadDataSCUC(FileName)
      %% 变量说明：******************************************************
% 输入参数说明：FileName：数据文件名
% 输出参数说明：SCUC_data.nodeNum:节点总数  SCUC_data.branchNum:支路数  standardP：基准容量  iterationMax：迭代最大数
%centerParameter：中心参数  IterationPrecision：迭代精度 objectFunc：目标函数类型 branchroute 线路号 BalanceI：平衡节点号  branchI：线路始节点i  branchJ：线路终节点j
% branchR：线路电阻  branchX：线路电抗  branchB：线路对地电纳  branchGroundI：接地支路节点i  branchGroundB：接地支路电纳 
% branchTransformerI:变压器节点i  branchTransformerJ：变压器节点j  branchTransformerR：变压器电阻  branchTransformerX：变压器电抗
% branchTransformerK：变压器变比  NodeI：运行参数节点i   NodePg:有功出力  NodeQg：无功出力  NodePl：有功负荷 
% NodeQl:无功负荷  PvI：PV节点i  PvV:PV节点电压  PvQmin：无功出力下限  PvQmax：无功出力上限  
% GenI:发电机节点号  GenC:耗量特性零次项系数  GenB:耗量特性一次项系数  GenA：耗量特性二次项系数   
% GenPmin：有功出力下限  GenPmax：有功出力上限  PvNum:Pv节点数  GenNum：发电机节点数
%% ***************************************************************

filedata=fopen(FileName); 

fgetl(filedata);     %忽略第一行注释
%%%  读取第2行数据
baseparameters = fscanf(filedata,'%f',[1,7]);
SCUC_data.baseparameters.busN          = baseparameters(1);              %%%  节点总数
SCUC_data.baseparameters.branchN        = baseparameters(2);              %%%  支路数
SCUC_data.baseparameters.balanceBus         = baseparameters(3);              %%%  平衡节点
SCUC_data.baseparameters.standardP        = baseparameters(4);              %%%  基准功率
SCUC_data.baseparameters.iterationMax     = baseparameters(5);              %%%  内点法最大迭代次数
SCUC_data.baseparameters.centerParameter  = baseparameters(6);              %%%  中心参数
SCUC_data.baseparameters.unitN          = baseparameters(7);              %%%  机组数
standardP = SCUC_data.baseparameters.standardP;

img1=fgetl(filedata);     %忽略注释
img1=fgetl(filedata);     %忽略注释

%%%  读取机组数据
checkZero = fscanf(filedata,'%d',1);
k=1;
while checkZero~=0
    paraG(k,:) = fscanf(filedata,'%f',[1,16]); %读取带排放费用的数据
%     paraG(k,:) = fscanf(filedata,'%f',[1,13]); %无排放数据
    k=k+1;
    checkZero  = fscanf(filedata,'%d',1);
end
SCUC_data.units.bus_G      = paraG(:,1);
SCUC_data.units.alpha      = paraG(:,2);
SCUC_data.units.beta       = paraG(:,3);
SCUC_data.units.gamma      = paraG(:,4);
SCUC_data.units.PG_up      = paraG(:,5)/standardP;
SCUC_data.units.PG_low     = paraG(:,6)/standardP;
SCUC_data.units.QG_up      = paraG(:,7)/standardP;
SCUC_data.units.QG_low     = paraG(:,8)/standardP;
SCUC_data.units.U_ini      = paraG(:,9);            % 初始连续运行时间？我不需要
SCUC_data.units.T_off      = paraG(:,10);           % 最小停机时段
SCUC_data.units.ramp       = paraG(:,12)/standardP; % 爬坡约束
SCUC_data.units.start_cost = paraG(:,13);
% SCUC_data.units.a      = paraG(:,16);
% SCUC_data.units.b       = paraG(:,15);
% SCUC_data.units.c      = paraG(:,14); %排放系数，没有时就注销
SCUC_data.units.N          = k - 1;                 % 机组数目


% %%%%%%%%%%%%%%%%%%%%%    目标函数系数回乘基准值     %%%%%%%%%%%%%%%%%%%%
% beta  = beta*standardP;                                              %
% gamma = gamma*(standardP*standardP);                                 %
SCUC_data.units.alpha = SCUC_data.units.alpha/SCUC_data.baseparameters.standardP;                                             %
SCUC_data.units.gamma = SCUC_data.units.gamma*SCUC_data.baseparameters.standardP;     
% SCUC_data.units.a = SCUC_data.units.a/SCUC_data.baseparameters.standardP; 
% SCUC_data.units.c = SCUC_data.units.c*SCUC_data.baseparameters.standardP; 
%
% %%%%%%%%%%%%%%%%%%%%%    目标函数系数回乘基准值     %%%%%%%%%%%%%%%%%%%%



%%%  读取节点电压上下界数据
fgetl(filedata);     %忽略注释
fgetl(filedata);     %忽略注释
checkZero = fscanf(filedata,'%d',1); 
k=1;
while checkZero~=0
    paraV(k,:) = fscanf(filedata,'%f',[1,2]);
    k=k+1;
    checkZero  = fscanf(filedata,'%d',1);
end
SCUC_data.busV.v_up  = paraV(:,1);     %%% 节点电压上界
SCUC_data.busV.v_low = paraV(:,2);     %%% 节点电压下界

%%%  读取支路数据 (不包括变压器支路 )
fgetl(filedata);     %忽略注释
fgetl(filedata);     %忽略注释
checkZero = fscanf(filedata,'%d',1);
k=1;
while checkZero~=0
    branch(k,:) = fscanf(filedata,'%f',[1,7]);
    k           = k+1;
    checkZero   = fscanf(filedata,'%d',1);
end
SCUC_data.branch.I = branch(:,1);     %%%  支路始点
SCUC_data.branch.J = branch(:,2);     %%%  支路终点
SCUC_data.branch.R = branch(:,4);     %%%  支路电阻
SCUC_data.branch.X = branch(:,5);     %%%  支路电抗
SCUC_data.branch.B = branch(:,6);     %%%  支路对地电纳  对应输电线pi型等值电路的对地电纳
SCUC_data.branch.P = branch(:,7)/standardP;   %%% 支路传输功率上界

%%%   读取变压器支路数据
fgetl(filedata);     %忽略注释
fgetl(filedata);     %忽略注释
checkZero = fscanf(filedata,'%d',1);
k=1;
while checkZero~=0
    branchTransformer(k,:) = fscanf(filedata,'%f',[1,8]);
    k                      = k+1;
    checkZero              = fscanf(filedata,'%d',1);
end
SCUC_data.branchTransformer.I = branchTransformer(:,1);
SCUC_data.branchTransformer.J = branchTransformer(:,2);
SCUC_data.branchTransformer.R = branchTransformer(:,4);
SCUC_data.branchTransformer.X = branchTransformer(:,5);
SCUC_data.branchTransformer.K = 1./branchTransformer(:,8);
SCUC_data.branchTransformer.P = branchTransformer(:,7)/standardP;


%%%  24时段负荷数据及旋转备用
fgetl(filedata);     %忽略注释
fgetl(filedata);     %忽略注释
checkZero = fscanf(filedata,'%d',1);
k=1;
while checkZero~=0
    load(k,:) = fscanf(filedata,'%f',[1,7]);
    k         = k+1;
    checkZero = fscanf(filedata,'%d',1);
end
SCUC_data.totalLoad.PD_T = load(:,1)/standardP;    %%%  表示 T 时段各个时段的总负荷
SCUC_data.totalLoad.QD_T = load(:,2)/standardP;
R_T  = load(:,5:7)/standardP;
SCUC_data.totalLoad.R_T  = sum(R_T,2);                  %% 关于负荷的参数还不太了解。。。。????   理解为三个加起来为发电机的旋转备用约束
SCUC_data.totalLoad.T = size(SCUC_data.totalLoad.PD_T,1);

%%%  节点负荷数据
fgetl(filedata);     %忽略注释
fgetl(filedata);     %忽略注释
checkZero=fscanf(filedata,'%d',1);
k=1;
while checkZero~=0
    nodeSequence(k) = checkZero;
    nodePQ(k,:)     = fscanf(filedata,'%f',[1,2]);
    k               = k+1;
    checkZero       = fscanf(filedata,'%d',1);
end
SCUC_data.busLoad.bus_PDQR = nodeSequence';
SCUC_data.busLoad.node_P   = nodePQ(:,1)/standardP;
SCUC_data.busLoad.node_Q   = nodePQ(:,2)/standardP;

SCUC_data.busLoad.bus_P_factor = SCUC_data.busLoad.node_P/sum(SCUC_data.busLoad.node_P);  %%%  节点负荷因子   各个节点占负荷的比例，按照这个比例从总负荷中分配负荷。
SCUC_data.busLoad.bus_Q_factor = SCUC_data.busLoad.node_Q/sum(SCUC_data.busLoad.node_Q);

SCUC_data.busLoad.node_P = SCUC_data.totalLoad.PD_T * SCUC_data.busLoad.bus_P_factor';
SCUC_data.busLoad.node_Q = SCUC_data.totalLoad.QD_T * SCUC_data.busLoad.bus_Q_factor';

% 读取排放数据 
fgetl(filedata);     %忽略注释
fgetl(filedata);     %忽略注释
% Emission(1,:)     = fscanf(filedata,'%f',[1,5]);
% SCUC_data.emission.E0 = Emission(1);
% SCUC_data.emission.pi_b = Emission(2);
% SCUC_data.emission.pi_s = Emission(3);
% SCUC_data.emission.deta_E_b_max = Emission(4);
% SCUC_data.emission.deta_E_s_max = min(Emission(5),SCUC_data.emission.E0);


