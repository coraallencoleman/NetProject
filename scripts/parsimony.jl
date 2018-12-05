#run soft parsimony from sequences in alignments_1183genes

#to run on biostat server
#scp /Users/cora/git_repos/NetProject/scripts/parsimony.jl allencoleman@adhara.biostat.wisc.edu:/ua/allencoleman/Phylo/scripts
#to move data to server
#scp -r /Users/cora/git_repos/NetProject/data/ allencoleman@adhara.biostat.wisc.edu:/ua/allencoleman/Phylo/data

#run with nohup julia0.7 parsimony.jl > parsimony.out 2> parsimony.err

cd("/ua/allencoleman/Phylo/results/")

using PhyloNetworks, RCall, PhyloPlots, CSV, DataFrames
#R"name <- function(x) file.path('..', 'assets', 'figures', x)" #? need this?

# Read in Sequence Data for maxParsimonyNet(T::HybridNetwork, df::DataFrame)
#df = data frame containing the species names in column 1, or in a column named species or taxon
#read in all files in sequence directory
files = filter(r".phy$", readdir("/ua/allencoleman/Phylo/data/data/Cui_etal/alignments_1183genes/"))
df = vcat([CSV.read(f; header=false, delim=' ', datarow=2) for f in files]) #TODO look for an external process_exited 

#check order of taxa in file. 
#trim empty space-only columns
df = df[:,[:Column1, :Column12]]

# Subset DF for testing
#datsubset = df[1:10]

#read in RAxML starting tree
besttrees = readMultiTopology("/ua/allencoleman/Phylo/data/data/Cui_etal/raxml_1183genes/besttrees.tre");
starttree = besttrees[1]; #starting tree

#Run Parsimony (outgroup from Claudia's paper)
@time  net1 = maxParsimonyNet(starttree, df, outgroup="Xmonticolus") #Priapella?
writeTopology(net1, "bestnets_Parsimony.tre")

# Calculate parsimony score
score = parsimonyGF(net,species,traits,:softwired)
println(score)

#calculate distance
#How far are they from each other?
goldNet = readMultiTopology("/ua/allencoleman/Phylo/data/data/Cui_etal/snaq/bestnets_calibrated_cleanNames.tre")

# parsimonyNet = readTopology("results/bestnets_Parsimony.tre")
dist = hardwiredClusterDistance(goldNet[1], net1, false)
println(dist)

# FIGURES
R"svg(name('parsimony-fixed-net.svg'), width=4, height=4)"; # hide
R"par"(mar = [0,0,0,0]);
plot(net1, :R)
R"dev.off"(); # hide