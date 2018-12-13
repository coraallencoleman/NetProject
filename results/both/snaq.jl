
#run on biostat server
#scp /Users/cora/git_repos/NetProject/scripts/snaq.jl allencoleman@adhara.biostat.wisc.edu:/ua/allencoleman/Phylo/scripts
#move data to server
#scp -r /Users/cora/git_repos/NetProject/data/ allencoleman@adhara.biostat.wisc.edu:/ua/allencoleman/Phylo/data
#open julia with julia0.6
#run with nohup julia0.6 snaq.jl > snaq.out 2> snaq.err in scripts/

Pkg.add("PhyloNetworks")
#kg.add("PhyloPlots") #doesnt work on server bc RCall
using PhyloNetworks
#using PhyloPlots
cd("/ua/allencoleman/Phylo/results/")

besttrees = readMultiTopology("/ua/allencoleman/Phylo/data/data/Cui_etal/raxml_1183genes/besttrees.tre");
# plot(besttrees[1], :RCall, useEdgeLength=true, showEdgeNumber=true)
raxmlCF = readTrees2CF(besttrees, CFfile="tableCF.txt"); #creates table, slow
# readTrees2CF(besttrees, whichQ="rand", numQ=10, CFfile="tableCF10.txt") #takes random sample of 10 to speed things up

# #run SNaQ

@time net0 = snaq!(besttrees[1],raxmlCF, hmax=0, filename="net0", seed=1234)

writeTopology(net0, "bestnet_RAxMLSNaQ.tre")

goldNet = readMultiTopology("/ua/allencoleman/Phylo/data/data/Cui_etal/snaq/bestnets_calibrated_cleanNames.tre")
#net0read = readTopology("bestnet_RAxMLSNaQ.tre")
dist = hardwiredClusterDistance(goldNet[1], net0, false)
print("hardwired cluster distance distance: ")
println(dist)

