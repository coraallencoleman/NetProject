import PhyloNetworks;
import Base.Test;
import StatsFuns;
import DataFrames;

@everywhere begin
    using PhyloNetworks;
    using Base.Test;
    using StatsFuns;
    using DataFrames;
    include("ticr_new.jl")
end

N = 400
net1 = readTopology("(s17:4.588553191489361,(((s3:3.660276595744681,(s4:2.99627659574468,s5:2.99627659574468)I1:0.664)I2:0.157,(((s6:0.7692765957446803,s7:0.7692765957446803)I3:1.342,(s8:1.6582765957446806,#H24:0.0::0.279)I4:0.45299999999999985)I5:1.213,((s9:2.7632765957446805,((s10:0.7922765957446805,s11:0.7922765957446805)I6:1.006,(s12:0.8922765957446805,s13:0.8922765957446805)I7:0.906)I8:0.965)I9:0.07,((s14:0.9442765957446806,(s15:0.35227659574468095,s16:0.35227659574468095)I10:0.592)I11:0.714)#H24:1.175::0.721)I12:0.491)I13:0.493)I14:0.419,(((s18:1.82027659574468,s19:1.82027659574468)I15:0.196,(s20:1.5752765957446804,(s21:0.8002765957446805,s22:0.8002765957446805)I16:0.775)I17:0.441)I18:0.894,(s23:2.85527659574468,(s1:1.5472765957446803,s2:1.5472765957446803)I19:1.308)I20:0.055)I21:1.326)I22:0.3522765957446809);");
pMat = SharedArray{Float64,2}(N,4)
pAry_1 = SharedArray{Float64,2}(8855,N)
pAry_2 = SharedArray{Float64,2}(8855,N)
pAry_3 = SharedArray{Float64,2}(8855,N)
pAry_4 = SharedArray{Float64,2}(8855,N)

@sync @parallel for i in 1:N
    sed = rand(Int64,1)[1]
    hlcommand = `./apps/hybrid-Lambda/src/hybrid-Lambda -spcu '(s17:4.588553191489361,(((s3:3.660276595744681,(s4:2.99627659574468,s5:2.99627659574468)I1:0.664)I2:0.157,(((s6:0.7692765957446803,s7:0.7692765957446803)I3:1.342,(s8:1.6582765957446806,h24#0.279:0.0)I4:0.45299999999999985)I5:1.213,((s9:2.7632765957446805,((s10:0.7922765957446805,s11:0.7922765957446805)I6:1.006,(s12:0.8922765957446805,s13:0.8922765957446805)I7:0.906)I8:0.965)I9:0.07,((s14:0.9442765957446806,(s15:0.35227659574468095,s16:0.35227659574468095)I10:0.592)I11:0.714)h24#0.279:1.175)I12:0.491)I13:0.493)I14:0.419,(((s18:1.82027659574468,s19:1.82027659574468)I15:0.196,(s20:1.5752765957446804,(s21:0.8002765957446805,s22:0.8002765957446805)I16:0.775)I17:0.441)I18:0.894,(s23:2.85527659574468,(s1:1.5472765957446803,s2:1.5472765957446803)I19:1.308)I20:0.055)I21:1.326)I22:0.3522765957446809)r;' -num 304 -seed $sed -o ./output/hybridLambda/gt-$N-$i.`;
    run(hlcommand);
    run(pipeline(`sed 's/_1//g' ./output/hybridLambda/gt-$N-$i._coal_unit`, stdout="./output/hybridLambda/out$N-$i.txt"));
    treelist = readMultiTopology("./output/hybridLambda/out$N-$i.txt");
    listCF = readTrees2CF(treelist,CFfile="./output/obsCF/tableCF$N-$i.txt");
    net1_1 = net1_2 = net1_3 = net1_4 = deepcopy(net1)    
    result1 = ticr!(net1_1,listCF,false,model="dirichlet",quartetstat="maxCF")
    p1 = result1[1]
    pAry_1[:,i] = result1[5]
    result2 = ticr!(net1_2,listCF,false,model="dirichlet",quartetstat="minpval")
    p2 = result2[1]
    pAry_2[:,i] = result2[5]
    result3 = ticr!(net1_3,listCF,false,quartetstat="pearson")
    p3 = result3[1]
    pAry_3[:,i] = result3[5]
    result4 = ticr!(net1_4,listCF,false)
    p4 = result4[1]
    pAry_4[:,i] = result4[5]
    pMat[i,:] = [p1,p2,p3,p4]
end

writedlm("./output/pMat/pMat-$N.txt", pMat)
writedlm("./output/pvalArray/pval-DmaxCF-$N.txt", pAry_1)
writedlm("./output/pvalArray/pval-DminP-$N.txt", pAry_2)
writedlm("./output/pvalArray/pval-BX-$N.txt", pAry_3)
writedlm("./output/pvalArray/pval-BQ-$N.txt", pAry_4)
