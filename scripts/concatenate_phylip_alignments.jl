# concatenate several phylip files into one file full.phy (Julia 1.0)
# by Cora Allen-Coleman, 2018-12-06

# assumes first line = header (no blank line before)
# assumes that each taxon are sampled in all alignments:
# taxa in the first alignment but missing from another alignment
# will trigger an error at the very end

cd("/Users/cora/git_repos/NetProject/data/Cui_etal/alignments_1183genes")
files = filter!(x -> occursin(r"\.phy$",x), readdir());

phyDict = Dict{String,String}()
taxa1 = String[] # taxa in the order in which they come in first file
filecounter = 0
for file in files
    # global taxa1 # for testing, when run in the REPL
    global filecounter += 1
    open(file) do f # open for reading by default
        readline(f) # read & discard header: "numTaxa numSites"
        counter=0
        for line in eachline(f)
            counter += 1
            m = match(r"\s*(\w+)\s+(\S+)\s*$", line) # strips away leading & trailing spaces
            m != nothing || continue # to next line. accommodates blank or weird lines
            taxon = m.captures[1]
            dna = m.captures[2]
            if haskey(phyDict, taxon)
                phyDict[taxon] *= dna
            elseif filecounter==1
                phyDict[taxon]  = dna
                push!(taxa1, taxon)
            else
                @warn "taxon $taxon appears in some alignment $filecounter: will be ignored"
            end
            # not handled here: taxon in file 1 but not in file i>1
            # counter<2 || break # for testing
        end
    end
end
totallength = unique(length(v) for v in values(phyDict));
length(totallength) == 1 ||
  error("concatenated sequences have variable lengths: $totallength")
totallength = totallength[1]

outfile = "full.phy"
open(outfile, "w") do f
    println(f, length(phyDict), " ", totallength)
    for tax in taxa1 # same order as in file 1
      write(f, tax, " ", phyDict[tax], "\n")
    end
end
