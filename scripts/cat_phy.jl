files = filter(r".phy$", readdir("/Users/cora/git_repos/NetProject/data/Cui_etal/alignments_1183genes/"))
# df = vcat([CSV.read(f; header=false, delim=' ', datarow=2) for f in files]) #TODO look for an external process_exited 


#cd("/Users/cora/git_repos/NetProject/data/Cui_etal/alignments_1183genes")
#check order of taxa in first file. 

phyDict = Dict{String,String}()
for f in files
    open(f) do f # open for reading by default
        counter = 0
        for line in eachline(f)
            if counter >= 1
                line = strip(line)
                dat = split(line," ")
                # println(dat[1])
                # println(dat[12])
                taxon = dat[1]
                seq = string(dat[12], phyDict[taxon]) #add sequence to current Dict
                phyDict[taxon] = seq
            end
            counter += 1
        end
    end
end


