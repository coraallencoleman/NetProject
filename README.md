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
raxmlCF = readTrees2CF(besttrees, CFfile="tableCF.txt");
#readTrees2CF(besttrees, whichQ="rand", numQ=10, CFfile="tableCF10.txt") #takes random sample of 10 to speed things up

#SNaQ
net0 = snaq!(besttrees[1],raxmlCF, hmax=0, filename="net0", seed=1234)
```


### Figures
```julia
using PhyloPlots
```

### Metrics
```julia
Pkg.add("CPUTime")
using CPUTime
@time @CPUtime net0 = snaq!(besttrees[1],raxmlCF, hmax=0, filename="net0", seed=1234)
```

# Hard-wired Parsimony
see https://github.com/crsl4/PhyloNetworks.jl/blob/master/docs/src/man/parsimony.md

which data source works for parsimony?

adapted from manual
```juila
using PhyloNetworks
mkpath("../assets/figures")
using RCall
R"name <- function(x) file.path('..', 'assets', 'figures', x)"
# net1 = readTopology(joinpath(Pkg.dir("PhyloNetworks"),"examples","swadesh.out"))
# we would get net1 from analyzing the complete data, but not available with the package #? What does it mean?
using CSV, DataFrames
csvfile = joinpath(Pkg.dir("PhyloNetworks"),"examples","Swadesh.csv");
dat = CSV.read(csvfile);
head(dat) # head shows the first 6 rows only
species, traits = PhyloNetworks.readCSVtoArray(dat);
species
traits
net = readTopology("(Spanish,((English)#H1,(Norwegian,(German,#H1))));");
using PhyloPlots, RCall
R"svg(name('parsimony-fixed-net.svg'), width=4, height=4)"; # hide
R"par"(mar = [0,0,0,0]);
plot(net, :R, xlim=[0.8,7.5]);
R"dev.off"(); # hide
```


```julia
#score = parsimonySoftwired(net, species, traits)
score = parsimonyGF(net,species,traits,:softwired)
hardscore = parsimonyGF(net,species,traits,:hardwired) # error
starttree = readTopology("(((English,German),Norwegian),Spanish);");
net1 = maxParsimonyNet(starttree, dat, hmax=1, outgroup="Spanish", rootname="swadesh")
```

# Questions for 11 28 18 Meeting with Cécile

## SNaQ
Should I use MrBayes (+ BuckyCF, considers gene tree error weights those tree with a lower weight (informally) post prob for each gene) or RAxML (which doesnt consider gene tree uncertainty)? Both? Or just keep it consistent?

iqtree instead of RAxML (more accurate)
## Parsimony
softwired
use sequences for softwired

might be slow because it doesnt summarize patterns (yet) can take subset of genes (random 100)


## Metrics
time https://github.com/schmrlng/CPUTime.jl

Robinson-Foulds distance
http://crsl4.github.io/PhyloNetworks.jl/latest/man/dist_reroot/



