#run soft parsimony from sequences on alignments_1183genes

#"/Users/cora/git_repos/NetProject/data/Cui_etal/alignments_1183genes"

using PhyloNetworks, RCall, PhyloPlots, CSV, DataFrames
#R"name <- function(x) file.path('..', 'assets', 'figures', x)" #? need this?

# Read in Sequence Data for maxParsimonyNet(T::HybridNetwork, df::DataFrame)
#df = data frame containing the species names in column 1, or in a column named species or taxon
#read in all files in sequence directory
cd("$(homedir())/git_repos/NetProject/data/Cui_etal/alignments_1183genes/")
files = filter(r".phy$", readdir("/Users/cora/git_repos/NetProject/data/Cui_etal/alignments_1183genes/"))
df = vcat([CSV.read(f; header=false, delim=' ', datarow=2) for f in files])
#trim bad space-only columns
df = df[:,[:Column1, :Column12]]

# Subset DF for testing
datsubset = df[1:100]

#read in RAxML starting tree
besttrees = readMultiTopology("/Users/cora/git_repos/NetProject/data/Cui_etal/raxml_1183genes/besttrees.tre");
starttree = besttrees[1]; #starting tree

#Run Parsimony (outgroup from Claudia's paper)
net1 = maxParsimonyNet(starttree, datsubset, outgroup="Xmonticolus") #TODO subset just for testing

# Calculate parsimony score
score = parsimonyGF(net,species,traits,:softwired)

# FIGURES
mkpath("../assets/figures")
R"svg(name('parsimony-fixed-net.svg'), width=4, height=4)"; # hide
R"par"(mar = [0,0,0,0]);
plot(net1, :R)
R"dev.off"(); # hide