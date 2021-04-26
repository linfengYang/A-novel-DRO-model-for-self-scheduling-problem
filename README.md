A novel DRO model with CVaR for self-scheduling problem
================

What about this project/study?

      This study proposes a novel moment-based distributionally robust optimization (DRO) model with 
    conditional value-at-risk (CVaR) for self-scheduling problem under elecyticity price uncertainty.
    This model comprehensively considers electricity price fluctuations, unit parameters and load re-
    quirement, and it can be adjusted by adjusting the size of the ambiguity set, so it provides a 
    suitable and adjustable self-scheduling strategy for generation companies (GENCOs), decision makers
    can adjust the strategy according to the actual scenarios for the scheme to be accepted by independents
    system operators (ISOs) while maxmizing the generation profit. Such DRO models are usually translated 
    into semi-definite pro-gramming (SDP) for solution, however, solving large-scale SDP needs a lot of com-
    putational time and resources. For this short-coming, two effective approximate models are proposed: 
    one ap-proximate model based on vector splitting and another based on alternate direction multiplier me-
    thod (ADMM) algorithm, both can greatly reduce the calculation time and resources. In order to verify the
    correctness and effectiveness of our proposed model, we apply them on a various of numerical case studies,
    and compare with the model in [1].Simulations of three IEEE test systems (6-bus system, 30-bus system and 
    118-bus system) are conducted.
    
    [1]	R. A. Jabr, "Robust self-scheduling under price uncertainty using conditional value-at-risk," 
    IEEE Trans. Power Syst., vol. 20, no. 4, pp. 1852-1858, Nov. 2005.


User Guide
-----------

The description of folder 

    SCUC_dat : The numerical case studies data. 

    Picture : The figures and tables in this paper. 数据结果图表。




The description of implement code files  (函数文件说明)

    DCOPF_ADMM.m : The main funciton.  主函数。

    ReadDataSCUC :  Read the SCUC6.txt and SCUC1062.txt.  SCUC6.txt和SCUC1062.txt的读取函数。

    ReadDataDCDOPF : Read the DDOPF118.txt和RTS48.txt.  DDOPF118.txt和RTS48.txt的读取函数。

    SCUC_nodeY :  Construct network admittance matrix. 形成导纳矩阵的函数。

    partitionNode :  Set the partition of the system. 设置分区的函数。

    partitionDataPI :  The procedure of identifying "real" coupling boudnary branches and brandary buses. 识别耦合节点和耦合支路的函数。

    formMatrixA :  Corresponding to the constraint (17) in manuscript. 构造文章中约束(17)的系数矩阵。

    formMatrixM :  Corresponding to the constraint (18) in manuscript. 构造文章中约束(18)的系数矩阵。

    formQC4Emission : Corresponding to the constraint (19) in manuscript.  构造文章中约束(19)的系数矩阵。

    yanZheng :  Using Cplex to solve DC-DOPF-CET problem. 使用Cplex求解DC-DOPF-CET。

    formQCP_PI_x_i : Corresponding to the constraint (30) in manuscript. 构造文章中约束(30)的系数矩阵。





Prerequisite:
-----------

    Matlab R2018a
    Cplex 12.7.1
    Mosek 9.2




Publication:
-----------
    If you use our study in academic work then please consider citing our papers.




About Us 
-----------
    Authors：Lingfeng Yang (ylf@gxu.edu.cn),Ying Yang (907803678@qq.com),Guo Chen,Zhaoyang Dong
    Team：www.scholat.com/team/eidp
    Webpage: http://jians.gxu.edu.cn/default.do
