#figures for snaq

using PhyloNetworks, CSV, DataFrames, RCall, PhyloPlots

cd("/Users/cora/git_repos/NetProject/results")
net1 = readTopology("snaq_0_1.tre")
net2 = readTopology("snaq_0_2.tre")
net3 = readTopology("snaq_0_3.tre")

# Remove Extra Nodes, Reroot
PhyloNetworks.deleteLeaf!(net1, "Priapella")
PhyloNetworks.deleteLeaf!(net1, "Psjonesii")
PhyloNetworks.deleteLeaf!(net1, "Xbirref")
PhyloNetworks.deleteLeaf!(net1, "Xnezahuacoyotl")
rootonedge!(net1, 42);

PhyloNetworks.deleteLeaf!(net2, "Priapella")
PhyloNetworks.deleteLeaf!(net2, "Psjonesii")
PhyloNetworks.deleteLeaf!(net2, "Xbirref")
PhyloNetworks.deleteLeaf!(net2, "Xnezahuacoyotl")
rootonedge!(net2, 42);

PhyloNetworks.deleteLeaf!(net3, "Priapella")
PhyloNetworks.deleteLeaf!(net3, "Psjonesii")
PhyloNetworks.deleteLeaf!(net3, "Xbirref")
PhyloNetworks.deleteLeaf!(net3, "Xnezahuacoyotl")
rootonedge!(net3, 42);

## Figures
cd("/Users/cora/git_repos/NetProject")
R"postscript('figures/snaq01.eps',height=5,width=25, po=8)";
plot(net1, :R)
R"dev.off"(); # hide

cd("/Users/cora/git_repos/NetProject")
R"postscript('figures/snaq01.eps',height=5,width=25, po=8)";
plot(net2, :R)
R"dev.off"(); # hide

cd("/Users/cora/git_repos/NetProject")
R"postscript('figures/snaq01.eps',height=5,width=25, po=8)";
plot(net3, :R)
R"dev.off"(); # hide

#Distances
goldNets = readMultiTopology("/Users/cora/git_repos/NetProject/data/Cui_etal/snaq/bestnets_calibrated_cleanNames.tre")
goldNet = goldNets[3] #with hmax = 2
#Compute Distances
distances = Array{Int64,1}(zeros(10));
distances[1] = hardwiredClusterDistance(goldNet, net1, true)
distances[2] = hardwiredClusterDistance(goldNet, net2, true)
distances[3] = hardwiredClusterDistance(goldNet, net3, true)
