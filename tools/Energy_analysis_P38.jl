mutable struct FRAGS
      idx::Int
     name::String
   energy::Float64
  residue::String
end

# e.g.
#       julia Energy_analysis_P38.jl  ../P38_MFCC_folders  

target=ARGS[1]

if !isdir(target)
    println("The target suit $(target) is not exist")
    exit("Stopped. Reason: $(target) is not exist.")
else
    println("Monitor folder : ", target)
end

filelist = readdir(target)
println("")
println(filelist)
println("")

global fraglist=[]
#icount = 0
for ifile in filelist
    ifolder  = string(target,"/",ifile)
    #println("ifolder : ",ifolder)
    if isdir(ifolder)
        dirlist = readdir(ifolder)
        #println("dirlist : ",dirlist)
        for id in dirlist
            sline=split(id,".")
            #println("id(00) : ", id)
            if uppercase(sline[2])=="LOG"
                #println(" id(11) :  ",id)
                sline=split(split(id,".")[1],"-")
                #global icount = icount + 1
                # inum   = parse(Int32,sline[2])
                inum   = 0
                energy = 0
                outidx = string(ifolder,"/",id)

                open(outidx,"r") do readlog
                    while !eof(readlog)
                        line=readline(readlog)
                        sline=split(line)
                        ntmp=length(sline)
                        if ntmp > 4
                            if sline[1] == "SCF" && sline[2] == "Done:" 
                                energy = parse(Float64,sline[5])
                            end
                        end
                    end
                end

                #println(" frag(",inum,") : ",energy)

                #println(wlog,i," ",nlines)
                #println("FRAGS ",FRAGS(inum,id,energy))
                push!(fraglist,FRAGS(inum,id,energy,""))

            end
        end
    end
end

sort!(fraglist, by = x -> x.name)
#println("111",fraglist)
frag1=[]
frag2=[]
frag3=[]
frag4=[]
for i in 1:length(fraglist)
    println(fraglist[i])        
    if occursin("-ligand",fraglist[i].name)
        if occursin("CAP",fraglist[i].name)
            # println("inum-cap0", fraglist[i].name[4:7])
            inum = parse(Int32,fraglist[i].name[4:7])
            push!(frag3,FRAGS(inum,fraglist[i].name,fraglist[i].energy,""))
        end 
        if occursin("Chain",fraglist[i].name)
            # println("inum-cha0", fraglist[i].name[6:9])
            inum = parse(Int32,fraglist[i].name[6:9])
            push!(frag4,FRAGS(inum,fraglist[i].name,fraglist[i].energy,""))
        end 
    else
        if occursin("CAP",fraglist[i].name)
            # println("inum-cap1", fraglist[i].name[4:7])
            inum = parse(Int32,fraglist[i].name[4:7])
            push!(frag1,FRAGS(inum,fraglist[i].name,fraglist[i].energy,""))
        elseif occursin("Chain",fraglist[i].name)
            # println("inum-cha1", fraglist[i].name[6:9])
            inum = parse(Int32,fraglist[i].name[6:9])
            push!(frag2,FRAGS(inum,fraglist[i].name,fraglist[i].energy,""))
        elseif occursin("ligand0", fraglist[i].name)
            println(" ligand ", fraglist[i].name)
            global Eligand = fraglist[i].energy 
        end 
    end  
end

println(" The Eligand : ", Eligand)

for i in 1:length(frag1)
    println(i, "  ", frag1[i])
end
for i in 1:length(frag2)
    println(i, "  ", frag2[i])
end
for i in 1:length(frag3)
    println(i, "  ", frag3[i])
end
for i in 1:length(frag4)
    println(i, "  ", frag4[i])
end

deltaEChain=[]
for i in 1:length(frag4)
    if frag4[i].energy != 0.0 && frag2[i].energy != 0.0 
        dv = frag4[i].energy - frag2[i].energy - Eligand
        push!(deltaEChain,dv)
        println(" DeltaE(Chain) [",i,"] : ",dv)
    else
        println(i, " ===  skipped  === ")
    end 
end

deltaECAP=[]
for i in 1:length(frag3)
    if frag3[i].energy != 0.0 && frag1[i].energy != 0.0 
        dv = frag3[i].energy - frag1[i].energy - Eligand
        push!(deltaECAP,dv)
        println(" DeltaE(CAP) [",i,"] : ",dv)
    else
        println(i, " ===  skipped  === ")
    end 
end

println(" Sum :", sum(deltaEChain), sum(deltaECAP) )
dv0 = sum(deltaEChain) - sum(deltaECAP)
println(" Eint : ", dv0, " as the ", 627.5094*dv0," Kcal/mol") 


