mutable struct FRAGS
      idx::Int
    iatom::Int
   natoms::Int
  icharge::Int
 multiple::Int
     name::String
   energy::Float64
end

mutable struct EXCITON
    exmo1::Int
    exmo2::Int
    exci ::Float64
end

#mutable struct MOAO{T<:EXCITON}
mutable struct MOAO
    idx::Int
   nele::Int
   MOen::Any
      C::Any
      S::Any
end

#mutable struct EXSTATE{T<:EXCITON}
mutable struct EXSTATE
       idx::Int
    ex_ene::Float64
    wavenm::Float64
    oscstr::Float64
    S2    ::Float64
    exst  ::Any
end

mutable struct DATAEX
    info :: Any
      Ex :: Any
      MO :: Any
end


function readINFO(info)

    global fraglist=[]
    open(info,"r") do rinfo 
        for line in eachline(rinfo)
            sline=split(line)
            ntmp=length(sline)
            if ntmp > 6
                # idx,istart,natoms_frg,icharge,1,replace(modelpdb[imodel][i-1].residue," "=>"-"),0.0
                idx        =  parse(Int32,split(split(sline[4],("("))[2],(","))[1])
                iatom      =  parse(Int32,split(sline[5],",")[1])
                natoms     =  parse(Int32,split(sline[6],",")[1])
                icharge    =  parse(Int32,split(sline[7],",")[1])
                multiple   =  parse(Int32,split(sline[8],",")[1])
                name       =  split(sline[9],"\"")[2]
                energy     =  parse(Float64,split(sline[10],")")[1])               
                push!(fraglist,FRAGS(idx,iatom,natoms,icharge,multiple,name,energy)) 
            end 
        end 
    end 
end 

function readGAU_base(ilog,ifrag,idx)

    open(ilog,"r") do readlog
        while !eof(readlog)
            line=readline(readlog)
            sline=split(line)
            ntmp=length(sline)
            if ntmp > 0
                #println("line", line)
                if occursin("l101.exe",line)
                    line2=readline(readlog)
                    line2=readline(readlog)
                    line2=readline(readlog)
                    line2=readline(readlog)
                    line2=readline(readlog)
                    sline2=split(line2)
                    ntmp2=length(sline2)
                    # println(" line2 : ", line2 , " ntmp2 ", ntmp2, " : ",  sline2)
                    if ntmp2 == 6
                        ifrag.icharge  = parse(Int32,sline2[3])
                        ifrag.multiple = parse(Int32,sline2[6])
                    else
                        println(" Error in reading charge and multi ")
                        exit(1)
                    end    
                    # ===== Obtain the ground state energy =====                
                    while !eof(readlog)
                        line=readline(readlog)
                        sline=split(line)
                        ntmp=length(sline)
                        if ntmp > 4
                            if sline[1] == "Total" && sline[2] == "DFT" && sline[3] == "energy" 
                                ifrag.energy = parse(Float64,sline[5])
                                break
                            end 
                            if sline[1] == "SCF" && sline[2] == "Done:"
                                ifrag.energy = parse(Float64,sline[5])
                                break
                            end 
                        end
                    end

                end
            end 
        end 
    end
end 

function readGAU_exciton(ilog,idx,if_rem,if_dimer)

    istate = 0
    stateEX = []
    open(ilog,"r") do readlog
        while !eof(readlog)
            line=readline(readlog)
            if occursin("Excited State   ",line)
    
                sline   = split(line)
                ex_ene  = parse(Float64,sline[5])
                wavenm  = parse(Float64,sline[7])
                oscstr  = parse(Float64,split(sline[9],"=")[2])
                S2      = parse(Float64,split(sline[10],"=")[2])

                exst=[] 
                while !eof(readlog)
                    line2=readline(readlog)
                    sline2=split(line2)
                    ntmp2=length(sline2)
                    if occursin("->",line2)
                        exmo1 = parse(Int32,sline2[1])
                        exmo2 = parse(Int32,sline2[3])
                        exci  = parse(Float64,sline2[4])
                    else
                        break  
                    end
                    push!(exst,EXCITON(exmo1,exmo2,exci))
                end

                push!(stateEX,EXSTATE(idx,ex_ene,wavenm,oscstr,S2,exst))
                istate = istate + 1

                if if_dimer && if_rem
                    if istate == 2
                        break
                    end
                else 
                    if if_rem
                        break
                    end 
                end 

            end
        end                 
    end 
  
    exst=[]
    push!(exst,EXCITON(0,0,1.0))
    pushfirst!(stateEX,EXSTATE(idx,0.0,0.0,0.0,0.0,exst))

    return stateEX 

end 

function readGAU_MOsOVLP(ilog,idx)

    ifAoS = false
    ifMO  = false
    stateMOAO=[]
    open(ilog,"r") do readlog
        while !eof(readlog)
            line=readline(readlog)
            if occursin("Full point group",line)
                sline=split(line)
#                if uppercase(sline[4]) != "C1"
#                    println(" [point group : ", sline[4], " ] ",", but only C1 symmerty is supportted ")
#                    println(" Ignore this warning if nosym keyword is used! ")
#                    exit(1)
#                end 
            end
            if occursin("primitive gaussians",line)
                sline=split(line)
                global Norb = parse(Int32,sline[1])
                global MO  = zeros(Float64, Norb, Norb)
                global AoS = zeros(Float64, Norb, Norb)
                #AoS = Array{Float64,2}(undef, Norb, Norb)
            end 

            nocc   = 0
            MOen   = []            
            if occursin("*** Overlap ***",line)
                ifAoS = true 
                #println("overlap matrix matched - dim ", Norb)
                #println("AoS",AoS)
                icount=0
                while !eof(readlog)
                    line  = readline(readlog)
                    sline = split(line)
                    nx    = length(sline)
                    ix    = []  
                    if nx < 6 && line[1:7] == "       "
                        for i in 1:nx 
                            push!(ix, parse(Int32,sline[i]))
                        end
                        #println("ix : ", ix)
                        for i in 1:(Norb-icount*5)
                            line  = readline(readlog)
                            sline = split(line)
                            iy    = parse(Int32,sline[1])
                            # println("iy : ",iy)
                            if nx > i
                                jmax = i
                            else
                                jmax = nx 
                            end
                            for j in 1:jmax
                                #println(replace((sline[j+1]),"D"=>"E"))
                                dv = parse(Float64,replace(sline[j+1],"D"=>"E"))
                                AoS[ix[j],iy] = dv
                            end
                        end 
                        icount = icount + 1
                        if icount*5 >= Norb
                            break
                        end 
                    end
                end
                for i in 1:Norb-1
                    for j in i+1:Norb
                        AoS[j,i] = AoS[i,j]
                    end
                end
                # println("AoS",AoS)
                # into structure 
                # exit(1)
            end 

            if occursin("Molecular Orbital Coefficients",line)
                ifMO  = true
                # println("MO coeff. matched - dim ", Norb)
                icount = 0
                while !eof(readlog)
                    line  = readline(readlog)
                    sline = split(line)
                    nx    = length(sline)
                    imo    = []
                    for i in 1:nx 
                        push!(imo, parse(Int32,sline[i]))
                    end

                    line  = readline(readlog)
                    sline = split(line)
                    for i in 1:length(sline)
                        if sline[i]=="O"
                            nocc = nocc + 1
                        end
                    end 
                    line  = readline(readlog)
                    ioffset = 0
                    for i in 1:length(sline)
                        dv = parse(Float64,line[22+ioffset:31+ioffset])
                        push!(MOen,dv)
                        ioffset = ioffset + 10
                    end

                    for i in 1:Norb
                        line  = readline(readlog)
                        #println("line : ",line)
                        sline = split(line[22:22+nx*10-1])
                        for j in 1:nx
                            dv = parse(Float64,(sline[j]))
                            MO[i,imo[j]] = dv 
                        end 
                    end 
                    icount = icount + 1
                    if icount*5 >= Norb
                        break
                    end 
                end
                #println("MO",MO)
                #exit(1)
            end 
 
            if ifMO && ifAoS
                push!(stateMOAO,MOAO(idx,nocc*2,MOen,MO,AoS))
                break 
            end 
 
        end 
    end    

    return stateMOAO  

end 

function readELST(rdir,str_pop="",str_aos="",if_rem=false)

    idebug = -1000

    println("results : ",   rdir)
    println("str_pop : ",str_pop)
    println("str_aos : ",str_aos)

    filelist = readdir(rdir)
    println(filelist)     

    # idx = 0
    EXdata = []
    for ifile in filelist
        println(ifile)
        ifolder = string(rdir,"/",ifile) 
        if isdir(ifolder)
            dirlist = readdir(ifolder)
            for id in dirlist
                sline=split(id,".")
                #println("sline : ", sline)
                if uppercase(sline[2])=="GJF"
                    idx = parse(Int32,split(sline[1],"-")[2])
                    #println(sline, "  ",  idx)
                    outidx = string(ifolder,"/",sline[1],".log")                   
                    if idebug > 0
                        println(" outidx :  ",outidx)              
                    end 
 
                    readGAU_base(outidx,fraglist[idx],idx)    
                    if idebug > 0
                        println("fraglist[ ", idx, " ] : " , fraglist[idx])
                    end 

                    if_dimer = false
                    if occursin("_",fraglist[idx].name)
                        if_dimer = true
                    end
                    stateEX = readGAU_exciton(outidx,idx,if_rem,if_dimer)
                    if idebug > 0
                        for i in 1:length(stateEX)
                            println("fraglist-EX[ ", idx, " ] : " , stateEX)
                        end
                    end

                    stateMOAO = readGAU_MOsOVLP(outidx,idx)
                    if idebug > 0
                        println("fraglist-MOAO[ ", idx, " ] : " , stateMOAO[1].idx)
                    end 

                    push!(EXdata,DATAEX(fraglist[idx],stateEX,stateMOAO))

                end
            end 
        end 
    end  

    return EXdata

end



