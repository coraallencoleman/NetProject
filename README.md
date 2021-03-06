# Can We Improve the Speed of Phylogenetic Network Methods without Sacrificing Accuracy?
Computational Networks Course (BMI 826) Semester Project (University of Wisconsin-Madison)
Specifically, can we remove the complex and time-consuming incomplete lineage sorting piece of the model?

## Network Method incomplete Lineage Sorting with PseudoLikelihood SNaQ 
see tutorial: http://crsl4.github.io/PhyloNetworks.jl/latest/man/inputdata/

### MrBayes + Bucky + SNaQ
(according to Claudia's code snaq/calibration.jmd)
This is most time-consuming way to calculate the network, but also the most accurate. We'll consider this (data/Cui_etal/snaq/best0.tre) the **gold standard**.
output in `/snaq/` (from 2017-07 when Claudia re-ran SNaQ, but excluding "nezy".)

```julia
using PhyloPlots
cd("/Users/cora/git_repos/NetProject/data/Cui_etal/snaq")
goldNet = readTopology("best0.tre")
mkpath("../assets/figures")
R"svg(name('gold-fixed-net.svg'), width=4, height=4)"; # hide
plot(net1, :R)
R"dev.off"(); # hide
```

### RAxML + SNaQ
for comparison with soft parsimony

read in RAxML files
```julia
using PhyloNetworks
using PhyloPlots
besttrees = readMultiTopology("/Users/cora/git_repos/NetProject/data/Cui_etal/raxml_1183genes/besttrees.tre");
plot(besttrees[1], :RCall, useEdgeLength=true, showEdgeNumber=true)
raxmlCF = readTrees2CF(besttrees, CFfile="tableCF.txt"); #creates table, slow
readTrees2CF(besttrees, whichQ="rand", numQ=10, CFfile="tableCF10.txt") #takes random sample of 10 to speed things up

#run SNaQ
net0 = snaq!(besttrees[1],raxmlCF, hmax=0, filename="net0", seed=1234)
```

## Softwired Parsimony
see https://github.com/crsl4/PhyloNetworks.jl/blob/master/docs/src/man/parsimony.md

see cat_phy.jl and parsimony.jl script
   
# Notes from 11 28 18 Meeting with Cécile

## SNaQ
Should I use MrBayes (+ BuckyCF, considers gene tree error weights those tree with a lower weight (informally) post prob for each gene) or RAxML (which doesnt consider gene tree uncertainty)? Both? Or just keep it consistent?

future idea: try iqtree instead of RAxML (iqtree is new, more accurate and faster)

## Parsimony
softwired from sequences

might be slow because it doesnt summarize patterns (yet)
I started can take subset of genes (random 100)


## Metrics
I compare these SNaQ and softwired parsimony on two metrics:  
time https://github.com/schmrlng/CPUTime.jl

```julia
Pkg.add("CPUTime")
using CPUTime
@time @CPUtime net0 = snaq!(besttrees[1],raxmlCF, hmax=0, filename="net0", seed=1234)
```

Robinson-Foulds distance
http://crsl4.github.io/PhyloNetworks.jl/latest/man/dist_reroot/

### Figures
```julia
using PhyloPlots
```