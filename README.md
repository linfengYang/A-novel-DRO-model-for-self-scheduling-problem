A novel DRO model for self-scheduling problem
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

The description of implement code files  (函数文件说明)

    DRO_CVaR_Alg1 : The DRO model.  
    
    RO_CVaR : The RO model.
    
    Approximation1_of_DRO_CVaR : The APP1 model.
    
    Approximation2_of_DRO_CVaR : The APP1 mode2.
    
    ReadDataSCUC :  Read the data.
    
    SCUC_nodeY :  Construct network admittance matrix.
    
    acquire_lamda : Calculate price of electricity.
    
    MonteCarlo_Price : Generating fluctuating electricity price.
    
    portion : Divide the area.
    
    SCUC_X_Y ：IEEE X-bus Y-periods test system.



Prerequisite:
-----------

    Matlab R2018a
    Cplex 12.7.1
    Mosek 9.2




Publication:
-----------
    If you find our paper/code is useful to you, please cite it:
    "L. Yang, Y. Yang, G. Chen and Z. Dong, 'Distributionally Robust Framework and its Approximations 
    Based on Vector and Region Split for Self-Scheduling of Generation Companies,' 
    in IEEE Transactions on Industrial Informatics, vol. 18, no. 8, pp. 5231-5241, Aug. 2022, 
    doi: 10.1109/TII.2021.3125964."




About Us 
-----------
    Authors：Lingfeng Yang (ylf@gxu.edu.cn),Ying Yang (yingyoung1997@gmail.com),Guo Chen,Zhaoyang Dong
    Team：www.scholat.com/team/eidp
    Webpage: http://jians.gxu.edu.cn/default.do
