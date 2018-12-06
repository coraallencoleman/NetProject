#script to concatenate several phy files into one file full.phy

cd("/Users/cora/git_repos/NetProject/data/Cui_etal/alignments_1183genes")
files = filter(r"L.*\.phy$", readdir("/Users/cora/git_repos/NetProject/data/Cui_etal/alignments_1183genes/"));

phyDict = Dict{String,String}()
let filecounter = 1
for file in files
    open(file) do file # open for reading by default
        linecounter = 1
        if filecounter == 1 #for first file, creates a Dict entry for each taxon
            for line in eachline(file)
                if linecounter > 1 #skips header line
                    line = strip(line)
                    #split using the regex r”\s+” to capture all spaces so sequence = dat[2]
                    dat = split(line, r"\s+") 
                    taxon = dat[1]
                    phyDict[taxon] = dat[2] #add sequence
                end
                linecounter +=1
            end
        else #for later files, adds to current Dict entry
            for line in eachline(file)
                if linecounter > 1 #skips header
                    line = strip(line)
                    dat = split(line,r"\s+")
                    taxon = dat[1]
                    phyDict[taxon] *= dat[2]
                end
                linecounter += 1
            end
        end
    end
    filecounter += 1
end
end

outfile = "full.phy"
f = open(outfile, "w")
for key in keys(phyDict)
    println(f, key, " ", phyDict[key], "\n")
end
close(f)

#TODO check order of taxa in first file, write in that order


