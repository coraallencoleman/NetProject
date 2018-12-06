files = filter(r"L.*\.phy$", readdir("/Users/cora/git_repos/NetProject/data/Cui_etal/alignments_1183genes/"));

#cd("/Users/cora/git_repos/NetProject/data/Cui_etal/alignments_1183genes")
#check order of taxa in first file. 

phyDict = Dict{String,String}()
for f in files
    println(f)
    filecounter = 0
    open(f) do f # open for reading by default
        linecounter = 1
        if filecounter == 0
            for line in eachline(f)
                if linecounter > 1 #skips header
                    line = strip(line)
                    dat = split(line," ")
                    taxon = dat[1]
                    phyDict[string(taxon)] = dat[12] #add sequence to current Dict
                end
                linecounter +=1
            end
        else
            for line in eachline(f)
                if linecounter > 1 #skips header
                    line = strip(line)
                    dat = split(line," ")
                    taxon = dat[1]
                    phyDict[string(taxon)] = string(phyDict[taxon], dat[12]) #add sequence to current Dict
                end
                linecounter +=1
            end
        end
    end
    filecounter += 1
end

#write to phy file
cd("/Users/cora/git_repos/NetProject/data/Cui_etal/alignments_1183genes/")
outfile = "full.phy"
open(outfile, "w") do f
    for key in keys(phyDict)
        #print without quotes
        #print(phyDict[key])
        write(outfile, key, " ", phyDict[key], "\n")
        #writedlm(outfile, phyDict[key])
    end
end

outfile = "full.phy"
f = open(outfile, "w")
for key in keys(phyDict)
    println(f, key, " ", phyDict[key], "\n")
end
close(f)