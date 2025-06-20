mutable struct FRAGS
      idx::Int
     name::String
   energy::Float64
  residue::String
end

# e.g.
#       julia Energy_analysis_Gaussian.jl ../ACE2_Ab-Omicron_Q493R_0_0-8.0A-MFCC.info[INFO] ../COVID-10nodes-TEST-DFT-TRY4[FOLDER]

info=ARGS[1]
target=ARGS[2]

if !isfile(info)
    println("The information file $(info) is not exist")
    exit("Stopped. Reason: $(info) is not exist.")
else
    println("Information file : ", target )
end

if !isdir(target)
    println("The target suit $(target) is not exist")
    exit("Stopped. Reason: $(target) is not exist.")
else
    println("Monitor folder : ", target )
end

filelist = readdir(target)
println("")
println(filelist)
println("")

global fraglist=[]
#icount = 0
for ifile in filelist
    ifolder  = string(target,"/",ifile)
    println("ifolder : ",ifolder)
    if isdir(ifolder)
        dirlist = readdir(ifolder)                     
        #println("dirlist : ",dirlist)
        for id in dirlist
            sline=split(id,".")
            #println("id(00) : ", id) 
            if length(sline)>1
            if uppercase(sline[2])=="LOG" 
                #println(" id(11) :  ",id)
                sline=split(split(id,".")[1],"-")
                #global icount = icount + 1 
                inum   = parse(Int32,sline[2])
                energy = 0 
                outidx = string(ifolder,"/",id)                

                open(outidx,"r") do readlog
                    while !eof(readlog)
                        line=readline(readlog)
                        sline=split(line)
                        ntmp=length(sline)
                        if ntmp > 4
                            if sline[1] == "SCF" && sline[2] == "Done:" && sline[4] == "="
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
end  

#println("000",fraglist)
sort!(fraglist, by = x -> x.idx)
#println("111",fraglist)

open(info,"r") do readinfo
    line=readline(readinfo)
    idx = 0
    while !eof(readinfo)
        idx = idx + 1 
        line=readline(readinfo)
        sline=split(line)
        ssline=replace(sline[10],"\""=>"",","=>"")
        fraglist[idx].residue=ssline
    end  
end

# For matching
IDfrags = Array{String}(undef, length(fraglist))
for i in 1:length(fraglist)
    println("fraglist[",i,"] : ",fraglist[i].energy, " ",fraglist[i].name," ",fraglist[i].residue)
    IDfrags[i] = fraglist[i].residue
end 

#println(IDfrags)
global res=[]
for i in 1:length(fraglist)
    sline=split(fraglist[i].residue,"_")
    ntmp = length(sline)
    if ntmp == 1 
        ssline=split(sline[1],"-") 
        ntmp = length(ssline)
        if ntmp > 1            
            vec = []
            println("Interaction between spike's residue ",sline, " and hACE")
            for j in 1:length(fraglist)
                if occursin(sline[1],IDfrags[j])
                    push!(vec,j)
                end 
            end 
            push!(res,vec)
        end 
    end 
end

println(res)

Eints=[]
for i in 1:length(res)
    if length(res) == 1
        println("Interaction between spike's residue (",fraglist[res[i][1]].residue, ") and hACE is ", 0.0)
        push!(Eints,0.0)  
    else 
        Ex      = fraglist[res[i][1]].energy
        Ecap    = 0.0
        Echain  = 0.0
        Ecapx   = 0.0
        Echainx = 0.0
        Eint    = 0.0
        Ncap    = 0     # count
        Nchain  = 0
        Ncapx   = 0
        Nchainx = 0
        for j in 2:length(res[i])
            sline = split(IDfrags[res[i][j]],"_")
            ntmp  = length(sline) 
            if ntmp == 2 # dimer case
                ssline = sline[2]
                if uppercase(ssline[1:4])=="CAPA"
                    Ecapx   = Ecapx   + fraglist[res[i][j]].energy
                    Ncapx   = Ncapx   + 1
                    idxpos  = findfirst(==(ssline),IDfrags)
                    if idxpos !=  0
                        Ecap   = Ecap   + fraglist[idxpos].energy
                        Ncap   = Ncap   + 1
                    else
                        println("No $(ssline) CAP ??? ") 
                    end 
                elseif uppercase(ssline[1:4])=="CHAI"
                    Echainx = Echainx + fraglist[res[i][j]].energy 
                    Nchainx = Nchainx + 1 
                    idxpos  = findfirst(==(ssline),IDfrags)
                    if idxpos !=  0
                        Echain   = Echain   + fraglist[idxpos].energy
                        Nchain   = Nchain   + 1
                    else
                        println("No $(ssline) Chain ??? ") 
                    end 
                end 
            end
        end

        if Nchainx != Ncapx
            println("The edge case, need add code to handel this !!! ")
        else
            E1   =  Echainx - Echain - Nchain * Ex
            E2   =  Ecapx   - Ecap   - Ncap   * Ex
            Eint =  E1 - E2
        end

        push!(Eints,Eint)  
        println("Interaction between spike's residue (",fraglist[res[i][1]].residue, ") and hACE is ", Eint)

    end 
end

sline=split(info,"/")
ntmp=length(sline)
ssline=split(sline[ntmp],".")

outlog=string(ssline[1],"_result.log")
open(outlog,"w") do wlog
    for i in 1:length(res)
        println(wlog, i," ",fraglist[res[i][1]].residue, " ",Eints[i]," ",627.5094*Eints[i])
    end
    for i in 1:length(res)
        println(wlog,"Interaction between spike's residue (",fraglist[res[i][1]].residue, ") and hACE is ", Eints[i],"A.U. equally as ",627.5094*Eints[i]," Kcal/mol")
    end
end 




