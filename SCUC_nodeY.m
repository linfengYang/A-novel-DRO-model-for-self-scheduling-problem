%%%   ylf编写于2017.12.7
%%%   求节点导纳矩阵
function Y = SCUC_nodeY(SCUC_data, type) 

Y.type = type;  % type = 'DC'为直流潮流   type = 'AC'为交流潮流的Y
if strcmp(type, 'DC') == 1
    SCUC_data.branch.R = 0;
    SCUC_data.branchTransformer.R = 0;
end


n = SCUC_data.baseparameters.busN;  %%%  

%%%   形成支路导纳矩阵
Y1        = sparse(1./(SCUC_data.branch.R + 1i * SCUC_data.branch.X));
Y11       = sparse(SCUC_data.branch.I,SCUC_data.branch.J,Y1,n,n);
branchYij = sparse(-Y11-Y11.');                     %%%  支路导纳的非对角元素

Ya        = sparse(SCUC_data.branch.I,SCUC_data.branch.J,1i*SCUC_data.branch.B,n,n);
Yc        = sparse(Ya+Ya.');
branchYii = sparse(diag(sum(-branchYij)+sum(Yc)));  %%%  形成支路导纳矩阵对角元素
branchY   = sparse(branchYij+branchYii);            %%%  形成支路导纳


%%%   变压器支路导纳
Y2                 = sparse(1./(SCUC_data.branchTransformer.R+1i*SCUC_data.branchTransformer.X));
transformerYij     = sparse(SCUC_data.branchTransformer.I,SCUC_data.branchTransformer.J,-SCUC_data.branchTransformer.K.*Y2,n,n);
transformerYij     = transformerYij+transformerYij.';                                                %%%  变压器导纳非对角元素
transformerYii     = sparse(SCUC_data.branchTransformer.I,SCUC_data.branchTransformer.I,(SCUC_data.branchTransformer.K.^2).*Y2,n,n);  %%%  变压器导纳对应i对角元素
transformerYjj     = sparse(SCUC_data.branchTransformer.J,SCUC_data.branchTransformer.J,Y2,n,n);                           %%%  变压器导纳对应J对角元素
branchTransformerY = sparse(transformerYii+transformerYjj+transformerYij);                           %%%  变压器支路导纳

%%%   形成节点导纳矩阵
YY = sparse(branchY+branchTransformerY);  %%%  支路导纳＋变压器导纳
Y.G = sparse(real(YY));
Y.B = sparse(imag(YY));