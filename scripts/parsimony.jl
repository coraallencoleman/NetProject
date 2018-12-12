#run soft parsimony from sequences in alignments_1183genes

#to run on biostat server
#scp /Users/cora/git_repos/NetProject/scripts/parsimony.jl allencoleman@adhara.biostat.wisc.edu:/ua/allencoleman/Phylo/scripts
#to move data to server
#scp /Users/cora/git_repos/NetProject/data/Cui_etal/alignments_1183genes/full.phy allencoleman@adhara.biostat.wisc.edu:/ua/allencoleman/Phylo/data/data/Cui_etal/alignments_1183genes/
#run with nohup julia0.6 parsimony.jl > parsimony.out 2> parsimony.err

#TO MOVE BACK
#scp allencoleman@adhara.biostat.wisc.edu:/ua/allencoleman/Phylo/scripts/parsimony* /Users/cora/git_repos/NetProject/results/parsimony/error/
#scp allencoleman@adhara.biostat.wisc.edu:/ua/allencoleman/Phylo/results/* /Users/cora/git_repos/NetProject/results/

cd("/ua/allencoleman/Phylo/")

using PhyloNetworks, CSV, DataFrames #,RCall, PhyloPlots
#R"name <- function(x) file.path('..', 'assets', 'figures', x)" #? need this?

# Read in Sequence Data for maxParsimonyNet(T::HybridNetwork, df::DataFrame)
#df = CSV.read("full.phy", delim = r"\s+")
df = readdlm("/ua/allencoleman/Phylo/data/data/Cui_etal/alignments_1183genes/full.phy", skipstart=1); #skips first line of file
df = convert(DataFrame, df);

#remove 5 species not present in gold standard network
#goldSpecies = readdlm("/ua/allencoleman/Phylo/data/data/Cui_etal/goldSpecies.txt")
#goldNets = readMultiTopology("/ua/allencoleman/Phylo/data/data/Cui_etal/snaq/bestnets_calibrated_cleanNames.tre")
#goldNet = goldNets[2]

#setdiff(df[:,1], goldSpecies)
df[5,1] = "Xbirchmanni"
df[6,1] = "Xclemenciae"
df[13,1] =  "Xmalinche"
pruned_df = df[2:25,:];
#remove 20
#setdiff(pruned_df[:,1], goldSpecies)
pruned_df = deleterows!(pruned_df, 19)
#setdiff(pruned_df[:,1], goldSpecies)

#read in RAxML starting tree
besttrees = readMultiTopology("/ua/allencoleman/Phylo/data/data/Cui_etal/snaq/bestnets_calibrated_cleanNames.tre");
starttree = besttrees[1]; #starting tree

#Run Parsimony (outgroup from Claudia's paper)
#rooted with the southern swordtails outgroup clade (SS).
@time  net1 = maxParsimonyNet(starttree, hmax = 1, pruned_df, outgroup="Xhellerii") #southern swordfishes
cd("/ua/allencoleman/Phylo/results/")
writeTopology(net1, "bestnets_Parsimony1.tre")

@time  net3 = maxParsimonyNet(starttree, hmax = 3, pruned_df, outgroup="Xhellerii") #southern swordfishes
cd("/ua/allencoleman/Phylo/results/")
writeTopology(net1, "bestnets_Parsimony3.tre")

@time  net4 = maxParsimonyNet(starttree, hmax = 4, pruned_df, outgroup="Xhellerii") #southern swordfishes
cd("/ua/allencoleman/Phylo/results/")
writeTopology(net1, "bestnets_Parsimony4.tre")

@time  net5 = maxParsimonyNet(starttree, hmax = 5, pruned_df, outgroup="Xhellerii") #southern swordfishes
cd("/ua/allencoleman/Phylo/results/")
writeTopology(net1, "bestnets_Parsimony5.tre")

@time  net6 = maxParsimonyNet(starttree, hmax = 6, pruned_df, outgroup="Xhellerii") #southern swordfishes
cd("/ua/allencoleman/Phylo/results/")
writeTopology(net1, "bestnets_Parsimony6.tre")

@time  net7 = maxParsimonyNet(starttree, hmax = 7, pruned_df, outgroup="Xhellerii") #southern swordfishes
cd("/ua/allencoleman/Phylo/results/")
writeTopology(net1, "bestnets_Parsimony7.tre")

@time  net8 = maxParsimonyNet(starttree, hmax = 8, pruned_df, outgroup="Xhellerii") #southern swordfishes
cd("/ua/allencoleman/Phylo/results/")
writeTopology(net1, "bestnets_Parsimony8.tre")

@time  net9 = maxParsimonyNet(starttree, hmax = 9, pruned_df, outgroup="Xhellerii") #southern swordfishes
cd("/ua/allencoleman/Phylo/results/")
writeTopology(net1, "bestnets_Parsimony9.tre")
