#run soft parsimony from sequences in alignments_1183genes

#to run on biostat server
#scp /Users/cora/git_repos/NetProject/scripts/parsimony.jl allencoleman@adhara.biostat.wisc.edu:/ua/allencoleman/Phylo/scripts
#to move data to server
#scp /Users/cora/git_repos/NetProject/data/Cui_etal/alignments_1183genes/full.phy allencoleman@adhara.biostat.wisc.edu:/ua/allencoleman/Phylo/data/data/Cui_etal/alignments_1183genes/
#scp allencoleman@adhara.biostat.wisc.edu:/ua/allencoleman/Phylo/results/* /Users/cora/git_repos/NetProject/results/parsimony/
#run with nohup julia0.6 parsimony.jl > parsimony.out 2> parsimony.err

cd("/ua/allencoleman/Phylo/data/data/Cui_etal/alignments_1183genes/")

using PhyloNetworks, CSV, DataFrames #,RCall, PhyloPlots
#R"name <- function(x) file.path('..', 'assets', 'figures', x)" #? need this?

# Read in Sequence Data for maxParsimonyNet(T::HybridNetwork, df::DataFrame)
#df = CSV.read("full.phy", delim = r"\s+")
df = readdlm("full.phy") 
df = convert(DataFrame, df)

#read in RAxML starting tree
besttrees = readMultiTopology("/ua/allencoleman/Phylo/data/data/Cui_etal/raxml_1183genes/besttrees.tre");
starttree = besttrees[1]; #starting tree

#Run Parsimony (outgroup from Claudia's paper)
@time  net1 = maxParsimonyNet(starttree, df, outgroup="Xmonticolus") #Priapella?
cd("/ua/allencoleman/Phylo/results/")
writeTopology(net1, "bestnets_Parsimony.tre")

# Calculate parsimony score
score = parsimonyGF(net1,species,traits,:softwired)
println(score)

#calculate distance
#How far are they from each other?
goldNet = readMultiTopology("/ua/allencoleman/Phylo/data/data/Cui_etal/snaq/bestnets_calibrated_cleanNames.tre")

# parsimonyNet = readTopology("results/bestnets_Parsimony.tre")
dist = hardwiredClusterDistance(goldNet[1], net1, false)
println(dist)

# # FIGURES
# cd("/ua/allencoleman/Phylo/results/")
# R"svg(name('parsimony-fixed-net.svg'), width=4, height=4)"; # hide
# R"par"(mar = [0,0,0,0]);
# plot(net1, :R)
# R"dev.off"(); # hide