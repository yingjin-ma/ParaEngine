mutable struct FRAGS
      idx::Int
     name::String
   energy::Float64
  residue::String
end

mutable struct DRUGS
     name::String
 elements::String
end

# e.g.
#       julia Assemble2casinfo.jl  ../original_document/allsdfs_withH.sdf  ../Test-test100-opt   drugfile  istart  ifinish

info=ARGS[1]
target=ARGS[2]
drugfile=ARGS[3]
istart=parse(Int32,ARGS[4])
ifinish=parse(Int32,ARGS[5])

if !isfile(info)
    println("The information SDF file $(info) is not exist")
    exit("Stopped. Reason: $(info) is not exist.")
else
    println("Information SDF file : ", target )
end

if !isdir(target)
    println("The target folder $(target) is not exist")
    exit("Stopped. Reason: $(target) folder is not exist.")
else
    println("Calculated results folder : ", target )
end

# ============================================================
#         Obtain the allsdfs_withH, reduce to string
# ============================================================

IDlist=[]
open(info,"r") do readlog
    while !eof(readlog)
        line=readline(readlog)
        sline=split(line)
        ntmp=length(sline)
        if ntmp > 0
            if sline[1] == "\$\$\$\$"
                line=readline(readlog)
                sline=split(line,"/")
                ntmp2=length(sline) 
                IDdrug=split(sline[ntmp2],".")[1]
                line=readline(readlog)
                line=readline(readlog)
                line=readline(readlog)
                StringElems=""
                while !eof(readlog) 
                    line2=readline(readlog)
                    sline2=split(line2)
                    ntmp2=length(sline2)             
                    #println("ntmp2 ", ntmp2)
                    if ntmp2 > 15   
                        #println(line2)
                        StringElems=string(StringElems,sline2[4])
                    else
                        break
                    end 
                end
#               println("StringElems : ", StringElems)
                push!(IDlist,DRUGS(IDdrug,StringElems))
            end
        end 
    end
end

# For matching
IDdrugs = Array{String}(undef, length(IDlist))
for i in 1:length(IDlist)
    IDdrugs[i] = IDlist[i].elements
end

#println("IDdrugs", IDdrugs)

# ==========================================================================
#     Obtain the calculated results, extract the elem-string string, too
# ==========================================================================

filelist = readdir(target)
println("")
println(filelist)
println("")

global fraglist=[]
global icount = 0
open(drugfile,"w") do wlog
for ifile in filelist
    ifolder  = string(target,"/",ifile)
    println("ifolder : ",ifolder)
    if isdir(ifolder)
        dirlist = readdir(ifolder)                     
#        #println("dirlist : ",dirlist)
        for id in dirlist
            sline=split(id,".")
            println("id(00) : ", id) 
            if uppercase(sline[2])=="LOG"
                println(" id(11) :  ",id)
                sline=split(split(id,".")[1],"-")
#                #global icount = icount + 1 
#                #inum   = parse(Int32,sline[2])
                energy = 0 
                outidx = string(ifolder,"/",id)               
 
                println(" outidx :  ",outidx)

##############################################################################
#         Localte the index first 
##############################################################################
                idxdrug=""
                open(outidx,"r") do readlog
                    while !eof(readlog)
                        line=readline(readlog)
                        sline=split(line)
                        ntmp=length(sline)
                        if ntmp > 0
                            if occursin("l101.exe",line)
                                line2=readline(readlog)
                                line2=readline(readlog)
                                sline2=split(line2,".")
                                ntmp2=length(sline2)
                                if uppercase(sline2[ntmp2]) == "PDB"
                                    sline3=split(sline2[ntmp2-1],"/")
                                    ntmp3=length(sline3)
                                    idxdrug=sline3[ntmp3]
                                else
                                    sline2=split(line2,"test") 
                                    # println("sline2 : ", sline2)
                                    itest=parse(Int32,sline2[2])  
                                    # println("itest : ", itest)

                                    for i in 1:3
                                        line2=readline(readlog)
                                    end 

                                    StringElems=""
                                    while !eof(readlog)  
                                        line2=readline(readlog)
                                        sline2=split(line2)
                                        ntmp2=length(sline2)

                                        # println("ntmp2",ntmp2)
                                        if ntmp2 == 4
                                            StringElems=string(StringElems,sline2[1])
                                        else
                                            break
                                        end 
                                    end
                                    # println("StringElems", StringElems)

                                    # for i in 1:length(IDdrugs)
                                    for i in istart:ifinish
                                        if StringElems == IDdrugs[i]
                                            println(StringElems, " - ", IDdrugs[i])
                                            if itest == i
                                                idxdrug = IDlist[i]
                                                break
                                            end                                             
                                            idxdrug = IDlist[i]
                                        end 
                                    end 
                                end
                            end

##############################################################################
#         Obtain the properities from the output
##############################################################################                           

                            if occursin("Full point group",line)
                                print(wlog,idxdrug.name,"  ", idxdrug.elements, "  ")
                                print(wlog," [point group : ", sline[4], " ] ")
                            end 

                            if occursin("SCF Done:",line)
                                print(wlog," [Energy(B3LYP/6-31G*) : ", sline[5], " ] ")
                            end

                            ncore=0  
                            if occursin("NAE=",line)
                                print(wlog," [N alpha electrons     : ", sline[4], " ] ")
                                print(wlog," [N  beta electrons     : ", sline[6], " ] ")
                                print(wlog," [Num. of core orbitals : ", sline[8], " ] ")
                                ncore=parse(Int32,sline[8])
                            end

                            if occursin("NOA=",line)
                                nhomo=parse(Int32,sline[4])
                                nhomo=nhomo+ncore  
                                print(wlog," [HOMO orbital (alpha)  : ", nhomo, " ] ")
                                print(wlog," [LUMO orbital (alpha)  : ", nhomo+1, " ] ")
                            end

                            if occursin("Excited State   1",line)
                                print(wlog," [Excited State   1 : ", sline[5], " eV, ", sline[7]," nm, ",  sline[9], ", ", sline[10], " ] ")
                                print(wlog," [ Excitation-manner 1 : ")
                                while !eof(readlog)
                                    line2=readline(readlog)
                                    sline2=split(line2)
                                    ntmp2=length(sline2)
                                    if occursin("->",line2)
                                        print(wlog," ( ")
                                        for i in 1:ntmp2
                                            print(wlog," ",sline2[i]," ")
                                        end
                                        print(wlog," ) ")
                                    else
                                        print(wlog," ] ")
                                        break
                                    end
                                end 
                            end

                            if occursin("Excited State   2",line)
                                print(wlog," [Excited State   2 : ", sline[5], " eV, ", sline[7]," nm, ",  sline[9], ", ", sline[10], " ] ")
                                print(wlog," [ Excitation-manner 2 : ")
                                while !eof(readlog)
                                    line2=readline(readlog)
                                    sline2=split(line2)
                                    ntmp2=length(sline2)
                                    if occursin("->",line2)
                                        print(wlog," ( ")
                                        for i in 1:ntmp2
                                            print(wlog," ",sline2[i]," ")
                                        end
                                        print(wlog," ) ")
                                    else
                                        print(wlog," ] ")
                                        break
                                    end
                                end
                            end

                            if occursin("Excited State   3",line)
                                print(wlog," [Excited State   3 : ", sline[5], " eV, ", sline[7]," nm, ",  sline[9], ", ", sline[10], " ] ")
                                print(wlog," [ Excitation-manner 3 : ")
                                while !eof(readlog)
                                    line2=readline(readlog)
                                    sline2=split(line2)
                                    ntmp2=length(sline2)
                                    if occursin("->",line2)
                                        print(wlog," ( ")
                                        for i in 1:ntmp2
                                            print(wlog," ",sline2[i]," ")
                                        end
                                        print(wlog," ) ")
                                    else
                                        print(wlog," ] ")
                                        break
                                    end
                                end
                            end

                            if occursin("Alpha  occ. eigenvalues",line)
                                dtmp1=0.0
                                while !eof(readlog)
                                    line2=readline(readlog)
                                    if occursin("Alpha virt. eigenvalues",line2)
                                        break
                                    end 
                                    sline2=split(line2)
                                    ntmp2=length(sline2)
                                    dtmp1=parse(Float64,sline2[ntmp2])
                                end
                                print(wlog," [ HOMO energy : ",dtmp1," hartree ] ")
#                               line2=readline(readlog)
                                sline2=split(line2)
                                ntmp2=length(sline2)
                                print(wlog," [ LUMO energy : ",sline2[5]," hartree ] ")
                                dtmp2=parse(Float64,sline2[5])
                                print(wlog," [ HOMO - LUMO gap : ",dtmp2-dtmp1," hartree ] ")
                            end

                            if occursin("Dipole moment",line) 
                                line2=readline(readlog) 
                                sline2=split(line2)
                                ntmp2=length(sline2)
                                print(wlog," [ Dipole moment : ",sline2[ntmp2]," Debye ] ")
                                println(wlog," ")
                            end

                        end 
                    end
                end
            end
        end  
    end
end 
end 

exit(0)

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




