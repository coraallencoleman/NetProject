#run locally with nohup with sudo nohup julia scripts/snaq_local.jl


Pkg.add("PhyloNetworks")
Pkg.add("PhyloPlots")
using PhyloNetworks
using PhyloPlots

besttrees = readMultiTopology("/Users/cora/git_repos/NetProject/data/Cui_etal/raxml_1183genes/besttrees.tre");
#plot(besttrees[1], :RCall, useEdgeLength=true, showEdgeNumber=true)
# raxmlCF = readTrees2CF(besttrees, CFfile="tableCF.txt"); #creates table, slow
# readTrees2CF(besttrees, whichQ="rand", numQ=10, CFfile="tableCF10.txt") #takes random sample of 10 to speed things up

raxmlCF =  readTableCF("tableCF.txt") #reads table
#run SNaQ
net0 = snaq!(besttrees[1],raxmlCF, hmax=0, filename="net0", seed=1234)

writeTopology(net0, "results/bestnets_RAxMLSNaQ.tre")


cd("/Users/cora/git_repos/NetProject/data/Cui_etal/snaq")
goldNet = readMultiTopology("bestnets_calibrated_cleanNames.tre")
net0read = readTopology("results/bestnets_RAxMLSNaQ.tre")
dist = hardwiredClusterDistance(goldNet[1], net0read, false)
print("hardwired cluster distance distance: ")
println(dist)


