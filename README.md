# Can we Improve the Speed of Likelihood-Based Phylogenetic Network Methods without Sacrificing Accuracy?
Project for Computational Networks Course (University of Wisconsin-Madison)
Julia

## SNaQ
see tutorial: http://crsl4.github.io/PhyloNetworks.jl/latest/man/inputdata/

read in MrBayes (according to Claudia's code snaq/calibration.jmd)
```julia
using PhyloNetworks
using PhyloPlots
net0 = readTopology("/Users/cora/git_repos/NetProject/data/Cui_etal/snaq/best0.tre");
net1 = readTopology("/Users/cora/git_repos/NetProject/data/Cui_etal/snaq/best1.tre");
net3 = readTopology("/Users/cora/git_repos/NetProject/data/Cui_etal/snaq/best3.tre");
net4 = readTopology("/Users/cora/git_repos/NetProject/data/Cui_etal/snaq/best4.tre");
net5 = readTopology("/Users/cora/git_repos/NetProject/data/Cui_etal/snaq/best5.tre");
plot(net0, :RCall, useEdgeLength=true, showEdgeNumber=true)
rootonedge!(net0, 15); # moves root edge
rootonedge!(net1, 43);
rootonedge!(net3, 42);
plot(net0, :RCall, useEdgeLength=true)
tree = readTopology("/Users/cora/git_repos/NetProject/data/Cui_etal/raxml_1183genes/besttrees.tre")
```

read in RAxML files
```julia
using PhyloNetworks
using PhyloPlots
besttrees = readMultiTopology("/Users/cora/git_repos/NetProject/data/Cui_etal/raxml_1183genes/besttrees.tre");
plot(besttrees[1], :RCall, useEdgeLength=true, showEdgeNumber=true)
raxmlCF = readTrees2CF(raxmltrees, CFfile="tableCF.txt");

#SNaQ
net0 = snaq!(besttrees[1],raxmlCF, hmax=0, filename="net0", seed=1234)
```


### Figures
```julia
using PhyloPlots
```

## Questions for 11 27 18 Meeting with Cécile
Should I use MrBayes or RAxML? 


