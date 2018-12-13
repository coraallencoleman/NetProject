#run soft parsimony from sequences in alignments_1183genes

#to run on biostat server
#scp /Users/cora/git_repos/NetProject/scripts/parsimony.jl allencoleman@adhara.biostat.wisc.edu:/ua/allencoleman/Phylo/scripts
#to move data to server
#scp /Users/cora/git_repos/NetProject/data/Cui_etal/alignments_1183genes/full.phy allencoleman@adhara.biostat.wisc.edu:/ua/allencoleman/Phylo/data/data/Cui_etal/alignments_1183genes/
#run with nohup julia0.6 parsimony.jl > parsimony.out 2> parsimony.err

#TO MOVE BACK
#scp allencoleman@adhara.biostat.wisc.edu:/ua/allencoleman/Phylo/scripts/* /Users/cora/git_repos/NetProject/results/both/
#scp allencoleman@adhara.biostat.wisc.edu:/ua/allencoleman/Phylo/results/* /Users/cora/git_repos/NetProject/results/both/

cd("/ua/allencoleman/Phylo/")

using PhyloNetworks, CSV, DataFrames #,RCall, PhyloPlots

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
# besttrees = readMultiTopology("/ua/allencoleman/Phylo/data/data/Cui_etal/snaq/bestnets_calibrated_cleanNames.tre");
# starttree = besttrees[2]; #starting tree has 1 reticulation

#scp /Users/cora/git_repos/NetProject/randomNetwork.tre allencoleman@adhara.biostat.wisc.edu:/ua/allencoleman/Phylo/data/
starttree = readTopology("/ua/allencoleman/Phylo/data/randomNetwork.tre")
cd("/ua/allencoleman/Phylo/results/")
#Run Parsimony (outgroup from Claudia's paper)
#rooted with the southern swordtails outgroup clade (SS).

for i in 1:10
    @time  net1 = maxParsimonyNet(starttree, hmax = i, Nfail = 1000, pruned_df, outgroup="Xhellerii") #southern swordfishes
    writeTopology(net1, @sprintf("bestnets_Parsimony%02.d.tre", i))
end
