#reads parsimony results from parsimony.jl. 
# Calculates distance
# Create figures

# Read in packages, network
using PhyloNetworks, CSV, DataFrames, RCall, PhyloPlots

cd("/Users/cora/git_repos/NetProject/")
#parsimonyNet = readTopology("results/bestnets_Parsimony2.tre")

# Create DF
# results = DataFrame()
# method, hmax, time, distance
# Calculate Distance
#How far is this network from the gold standard network?
goldNets = readMultiTopology("/Users/cora/git_repos/NetProject/data/Cui_etal/snaq/bestnets_calibrated_cleanNames.tre")
goldNet = goldNets[3] #with hmax = 3

#Read in Parsimony Results, 
parsimonyNet01 = readTopology("results/both/bestnets_Parsimony01.tre")
parsimonyNet02 = readTopology("results/both/bestnets_Parsimony02.tre")
parsimonyNet03 = readTopology("results/both/bestnets_Parsimony03.tre")
parsimonyNet04 = readTopology("results/both/bestnets_Parsimony04.tre")
parsimonyNet05 = readTopology("results/both/bestnets_Parsimony05.tre")
parsimonyNet06 = readTopology("results/both/bestnets_Parsimony06.tre")
parsimonyNet07 = readTopology("results/both/bestnets_Parsimony07.tre")
parsimonyNet08 = readTopology("results/both/bestnets_Parsimony08.tre")
parsimonyNet09 = readTopology("results/both/bestnets_Parsimony09.tre")
parsimonyNet10 = readTopology("results/both/bestnets_Parsimony10.tre")


#reroot networks on correct edge
#plot(parsimonyNet01, :RCall, showEdgeNumber=true)
rootonedge!(parsimonyNet01, 44);
rootonedge!(parsimonyNet02, 44);
rootonedge!(parsimonyNet03, 44);
rootonedge!(parsimonyNet04, 44);
rootonedge!(parsimonyNet05, 44);
rootonedge!(parsimonyNet06, 44);
rootonedge!(parsimonyNet07, 44);
rootonedge!(parsimonyNet08, 44);
rootonedge!(parsimonyNet09, 44);
rootonedge!(parsimonyNet10, 44);

#Compute Distances
distances = Array{Int64,1}(zeros(10));
distances[1] = hardwiredClusterDistance(goldNet, parsimonyNet01, true)
distances[2] = hardwiredClusterDistance(goldNet, parsimonyNet02, true)
distances[3] = hardwiredClusterDistance(goldNet, parsimonyNet03, true)
distances[4] = hardwiredClusterDistance(goldNet, parsimonyNet04, true)
distances[5] = hardwiredClusterDistance(goldNet, parsimonyNet05, true)
distances[6] = hardwiredClusterDistance(goldNet, parsimonyNet06, true)
distances[7] = hardwiredClusterDistance(goldNet, parsimonyNet07, true)
distances[8] = hardwiredClusterDistance(goldNet, parsimonyNet08, true)
distances[9] = hardwiredClusterDistance(goldNet, parsimonyNet09, true)
distances[10] = hardwiredClusterDistance(goldNet, parsimonyNet10, true)

#Network Graphs
R"postscript('figures/parsimony01.eps',height=5,width=25, po=8)";
plot(parsimonyNet01, :R)
R"dev.off"(); # hide

R"postscript('figures/parsimony02.eps',height=5,width=25, po=8)";
plot(parsimonyNet02, :R)
R"dev.off"(); # hide

R"postscript('figures/parsimony03.eps',height=5,width=25, po=8)";
plot(parsimonyNet03, :R)
R"dev.off"(); # hide

R"postscript('figures/parsimony04.eps',height=5,width=25, po=8)";
plot(parsimonyNet04, :R)
R"dev.off"(); # hide

R"postscript('figures/parsimony05.eps',height=5,width=25, po=8)";
plot(parsimonyNet05, :R)
R"dev.off"(); # hide

R"postscript('figures/parsimony06.eps',height=5,width=25, po=8)";
plot(parsimonyNet06, :R)
R"dev.off"(); # hide

R"postscript('figures/parsimony07.eps',height=5,width=25, po=8)";
plot(parsimonyNet07, :R)
R"dev.off"(); # hide

R"postscript('figures/parsimony08.eps',height=5,width=25, po=8)";
plot(parsimonyNet08, :R)
R"dev.off"(); # hide

R"postscript('figures/parsimony09.eps',height=5,width=25, po=8)";
plot(parsimonyNet09, :R)
R"dev.off"(); # hide

R"postscript('figures/parsimony10.eps',height=5,width=25, po=8)";
plot(parsimonyNet10, :R)
R"dev.off"(); # hide




# # FIGURES
cd("/Users/cora/git_repos/NetProject/figures")
R'postscript('parsimony01.eps',height=5,width=25, po=8)';
plot(parsimonyNet, :R)
R"dev.off"(); # hide

cd("/Users/cora/git_repos/NetProject/figures")
R"postscript('parsimony02.eps',height=5,width=25, po=8)";
plot(parsimonyNet02, :R)
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

#TODO create random network
#switch tips to create random network
writeTopology(goldNet, "randomNetwork.tre")
#create random network