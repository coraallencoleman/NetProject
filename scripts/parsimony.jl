#run soft parsimony from sequences on alignments_1183genes

#"/Users/cora/git_repos/NetProject/data/Cui_etal/alignments_1183genes"

using PhyloNetworks, RCall, PhyloPlots, CSV, DataFrames
#R"name <- function(x) file.path('..', 'assets', 'figures', x)" #? need this?

# Read in Sequence Data for maxParsimonyNet(T::HybridNetwork, df::DataFrame)
#df = data frame containing the species names in column 1, or in a column named species or taxon
dat = CSV.read("/Users/cora/git_repos/NetProject/data/Cui_etal/alignments_1183genes/Loc-9984-Tr-1-1-Conf-1.000-Len-1772_64_1772.phy"; header=false, delim=' ', datarow=2)
dat = dat[:,[:Column1, :Column12]]

#read in all files in sequence directory
cd("$(homedir())/git_repos/NetProject/data/Cui_etal/alignments_1183genes/")
files = filter(r".phy$", readdir("/Users/cora/git_repos/NetProject/data/Cui_etal/alignments_1183genes/"))
df = vcat([CSV.read(f; header=false, delim=' ', datarow=2) for f in files])
#trim bad space-only columns
df = df[:,[:Column1, :Column12]]

# Subset DF for Testing
datsubset = df[1:100]

#divide up species and traits #? in Col2, make each char its own column?

#Run
#score = parsimonySoftwired(net, species, traits)
score = parsimonyGF(net,species,traits,:softwired)

#read in RAxML starting tree
besttrees = readMultiTopology("/Users/cora/git_repos/NetProject/data/Cui_etal/raxml_1183genes/besttrees.tre");
starttree = besttrees[1]; #starting tree

net1 = maxParsimonyNet(starttree, datsubset) #TODO subset just for testing

# FIGURES
mkpath("../assets/figures")
R"svg(name('parsimony-fixed-net.svg'), width=4, height=4)"; # hide
R"par"(mar = [0,0,0,0]);
plot(net1, :R, xlim=[0.8,7.5]);
R"dev.off"(); # hide