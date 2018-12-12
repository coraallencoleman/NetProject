#reads parsimony results from parsimony.jl. 
# Calculates distance
# Create figures

# Read in packages, network
using PhyloNetworks, CSV, DataFrames, RCall, PhyloPlots

cd("/Users/cora/git_repos/NetProject/")
parsimonyNet = readTopology("results/bestnets_Parsimony2.tre")

# Calculate Distance

#How far is this network from the gold standard network?
goldNets = readMultiTopology("/Users/cora/git_repos/NetProject/data/Cui_etal/snaq/bestnets_calibrated_cleanNames.tre")
goldNet = goldNets[3] #with hmax = 2
goldSpecies = readdlm("/Users/cora/git_repos/NetProject/data/Cui_etal/goldSpecies.txt")

# # FIGURES
cd("/Users/cora/git_repos/NetProject/figures")
R"postscript('parsimony2.eps',height=5,width=25, po=8)";
plot(parsimonyNet, :R)
R"dev.off"(); # hide

# Gold Standard figure
cd("/Users/cora/git_repos/NetProject/figures")
R"postscript('gold_standard1.eps',height=5,width=25, po=8)";
plot(goldNets[1], :R)
R"dev.off"(); # hide

cd("/Users/cora/git_repos/NetProject/figures")
R"postscript('gold_standard2.eps',height=5,width=25, po=8)";
plot(goldNets[2], :R)
R"dev.off"(); # hide

cd("/Users/cora/git_repos/NetProject/figures")
R"postscript('gold_standard3.eps',height=5,width=25, po=8)";
plot(goldNets[3], :R)
R"dev.off"(); # hide

cd("/Users/cora/git_repos/NetProject/figures")
R"postscript('gold_standard4.eps',height=5,width=25, po=8)";
plot(goldNets[4], :R)
R"dev.off"(); # hide

#TODO
# Calculate parsimony score
# score = parsimonyGF(net1,species,traits,:softwired) #TODO fix
# println(score)
