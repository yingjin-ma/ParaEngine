using Distributed
using BenchmarkTools

@everywhere function read_last(file)
    r=""
    open(file) do io
        seekend(io)
        seek(io, position(io) - 2)
        while Char(peek(io)) != '\n'
            seek(io, position(io) - 1)
        end
        Base.read(io, Char)
        r=Base.read(io, String)
        return r
    end
end

@everywhere function read_last_abnormal(file)
    r=""
    open(file) do io
        seekend(io)       
        skip(io, -1024)  
        r=read(io,String)
        return r
    end
end

# For further details see manual section:

#=
@everywhere function check_last_NWChem_OLD(out,IOrecord,IFDONE)

    iflag = -1
    # println("out : ",out,"  IOrecord : ",IOrecord)
    sline1=rsplit(out,"/"; limit=2)
    # println("sline1 : ", sline1)
    record = open(IOrecord,"a+")
        lines=readlines(record)
        # println("lines : ",lines)
        if length(lines) > 0
            for line in lines
                sline2=rsplit(line,"/"; limit=2)
                # println("line : ", line," sline2: ",sline2," sline1: ",sline1)
                if occursin(sline1[2],sline2[2])
                    iflag = 1
                    break
                end
            end
            if iflag == 1
              println(record,out)
              # println("out : ",out)
            end
        else
            println(record,out)
            # println("out : ",out)
            iflag = 1
        end
    close(record)

    if IFDONE
        println("IFDONE is TRUE")
    else
        if iflag == 1
            while true
                sleep(5)
                tailout=read_last(out)
                println(tailout)
                if occursin("Total times  cpu:",tailout)
                    break
                end
            end
        end
    end

end
=#

@everywhere function check_last_NWChem(out,IOrecord,IFDONE)

    if IFDONE
        println("  ")
        println("IFDONE is TRUE; Skip this task ... ")
        println("  ")
    else    
        ii=0   
        while true
            ii = ii + 1
            sleep(5)
            tailout  = read_last(out)
            println(tailout)

            if occursin("Total times  cpu:",tailout)
                break
            end

#            if mod(ii,5) == 0
            println("")
            abnormal = read_last_abnormal(out)
            println(" abnormal check :", abnormal)
            println("")
            if occursin("For further details see manual",abnormal)
                break
            end 
            if occursin("Total times  cpu:",abnormal)
                break
            end 
#            end
        end
    end
    flush(stdout)
end


function nwchemtask(task,frag,atoms,par)

    # Generating NWChem headers
   
    fragidx=lpad(frag.idx,8,"0")
       name=frag.name
    #run(`cp Template.nwchem $(fragidx).nw`) 
 
    if frag.name!="unnamed"
        infile=string(workdir,"/$(name).nw")
    else 
        infile=string(workdir,"/Frag-$(fragidx).nw")
    end 

    open(infile,"w") do nwinp
        txt = read("Template.nwchem",String)
    #    # txt = replace(txt,"FRAG-i" => "Frag-$(fragidx)","X-library-basis" => "* library 6-31g","dftxc" => "xc m06-2x", "CHARGE" => "charge $(frag.icharge)" )
        txt = replace(txt,"SubTask-i" => "SubTask-$(fragidx)", "CHARGE" => "charge $(frag.icharge)" )
        print(nwinp,txt)
    end  

    inXYZ=string(workdir,"/XYZ.tmp")    
    open(inXYZ,"w") do xyztmp
        for iatom in atoms
            println(xyztmp,iatom.elem," ",iatom.coord[1]," ",iatom.coord[2]," ",iatom.coord[3])
        end
    end

    run(`sed -i "7r $inXYZ" $infile`)

    open(infile,"a") do nwinp
        for i in 1:length(QCpara)
            println(nwinp,QCpara[i])
        end 
    end
   
    #println(QCpara)
    #exit()

    if frag.name!="unnamed"
        task.infile  = "$(name).nw"
        task.outfile = "$(name).out" 
    else
        task.infile  = "Frag-$(fragidx).nw"
        task.outfile = "Frag-$(fragidx).out"
    end 

    flush(stdout)

end 

function gentask()

    @everywhere global tasklist=[]

    println("")
    if isdir(workdir)
        println("WorkDir : ",workdir)
    else 
        mkdir(workdir)
        println("WorkDir : ",workdir)
    end

#    @everywhere global IOrecord=string(workdir,"/IOrecord")
#    println("IOrecord : ",IOrecord)
#    io=open(IOrecord,"w") 
#    close(io)

    if qcdriver == "NWCHEM"
        println("Generating the NWChem tasks ... ")
        println("")        
        for i in 1:total_frags
             istart = fraglist[i].iatom
            ifinish = fraglist[i].iatom+fraglist[i].natoms-1
            print(" frags-",i," : atomlist ",istart,"-",ifinish)           
            push!(tasklist,TASKS(i,fraglist[i].idx,workdir,"","","NWCHEM",1)) 
            nwchemtask(tasklist[i],fraglist[i],atomlist[istart:ifinish],"DEFAULT")
            println(" ... done ")
        end  
        println("")
    end
 
    # println(tasklist)
    flush(stdout)

end 

function distributingtask(NSLURM)

    nnodes = -1
    iflag  = -1
    icount =  0 
    global pos_node=[]
    println("")
    println("Distributing the QC tasks ... ")
    if LBfile != ""        
        print("Assigning tasks with --> ")
        print(LBfile)
        println(" <-- ")
        # locate the overload tasks  1/2
        OLvec=[]
        # Read in the loadbalance file 
        open(LBfile,"r") do LBread
            while !eof(LBread)
                line=readline(LBread)
                sline=split(line)
                ntmp=length(sline)
                if ntmp > 0                   
                    if uppercase(sline[1]) == "NODES"
                        nnodes = parse(Int32,sline[2])                        
                        push!(pos_node,nnodes)
                        LBmat=[] 
                        for i in 1:nnodes
                            line=readline(LBread)
                            sline=split(line)
                            ntmp=length(sline)
                            LBvec=[]
                            for i in 1:ntmp
                                push!(LBvec,parse(Int32,sline[i]))
                                # If overload
                                if i>(1.0-overload)*ntmp && nnodes==NSLURM
                                    #println(" i , ntmp", i, ntmp)
                                    push!(OLvec,parse(Int32,sline[i]))
                                end   
                            end 
                            push!(LBmat,LBvec)
                        end  
                        push!(LBcube,LBmat)
                    end 
                    if uppercase(sline[1]) == "MULTINODES"
                        mnodes = parse(Int32,sline[2])
                        line=readline(LBread)
                        sline=split(line)
                        ntmp=length(sline) 
                        MNvec=[]
                        MNmat=[]
                        for i in 1:ntmp
                            push!(MNvec,parse(Int32,sline[i]))
                        end
                        push!(MNmat,mnodes)
                        push!(MNmat,MNvec)
                        push!(MNcube,MNmat)
                        #println("MNcube : ",MNcube)
                    end 
                end
            end
            #println("OLvec : ",OLvec)
            if nrepeat > 0
                OLvec2=[]
                for i in 1:length(OLvec)
                    if OLvec[i]>0
                         push!(OLvec2,OLvec[i])
                    end
                end
                nOLvec2=length(OLvec2)

                nredo = NSLURM/nrepeat
                neffe = floor(Int, nOLvec2/nredo + 1)
                println("nOLvec2 : ", nOLvec2 ,"OLvec2: ",OLvec2)

                div = nOLvec2/NSLURM
                nmod = mod(nOLvec2,NSLURM)

                for i in 1:length(pos_node)
                    if pos_node[i] == NSLURM
                        global ipos = i
                    end
                end

                iitv=0
                for i in 1:NSLURM
                    i1= iitv + 1
                    i2= iitv + neffe  
                    if length(OLvec2)-i2 <= 0
                        ix = i2 - length(OLvec2)
                        #println("OLvec2[i1:i2-ix] ",OLvec2[i1:i2-ix])
                        LBcube[ipos][i]=vcat(LBcube[ipos][i],OLvec2[i1:i2-ix])
                        if ix != 0 
                            #println("OLvec2[1:ix] ",OLvec2[1:ix])
                            LBcube[ipos][i]=vcat(LBcube[ipos][i],OLvec2[1:ix])
                        end 
                        iitv = ix
                    else
                        #println("OLvec2[i1:i2] ",OLvec2[i1:i2])
                        LBcube[ipos][i]=vcat(LBcube[ipos][i],OLvec2[i1:i2])
                    end
                    #println(i," -OLvec2 ", i1, " ", i2, " ")
                    iitv = iitv + neffe 
                    println("LBcube[",ipos,"][",i,"] : ",LBcube[ipos][i])
                end 

            end
            #exit()
        end
    else
        println("Assigning tasks with --> Default <-- ")       

        ndiv = total_frags/NSLURM
        nmod = mod(total_frags,NSLURM)

        # locate the overload tasks  1/2
        OLvec=[]
        for i in 1:NSLURM
            for j in 0:ndiv
                if j*NSLURM+i > total_frags*(1.0-overload)
                    if j*NSLURM+i <= total_frags
                        push!(OLvec,Int(j*NSLURM+i))
                    end  
                end 
            end                 
        end

        if nrepeat > 0  
           nredo = NSLURM/nrepeat 
           neffe = floor(Int, total_frags*overload/nredo + 1)
           #nintv = floor(Int, length(OLvec)/neffe)
           #println("nredo : ", nredo, " neffe : ",neffe," nintv ",nintv) 
           #println("OLvec: ",OLvec)
        end 

        # assign all the tasks      2/2
        iitv=0
        LBmat=[]
        for i in 1:NSLURM
            LBvec=[]
            for j in 0:ndiv
                if j*NSLURM+i <= total_frags
                    push!(LBvec,Int(j*NSLURM+i))
                end 
            end                
            if nrepeat > 0
                i1= iitv + 1
                i2= iitv + neffe
                #println("OLvec ", i1, " ", i2, " ")                 
                if length(OLvec)-i2 <= 0
                    ix = i2 - length(OLvec)
                    LBvec=vcat(LBvec,OLvec[i1:i2-ix])
                    if ix != 0
                        LBvec=vcat(LBvec,OLvec[1:ix])
                    end 
                    iitv = ix
                else
                    LBvec=vcat(LBvec,OLvec[i1:i2])
                end   
                #println(LBvec)
                push!(LBmat,LBvec)
                iitv = iitv + neffe
            else 
                LBvec=vcat(LBvec,OLvec)
                push!(LBmat,LBvec)
            end 
        end 
        push!(LBcube,LBmat)
 
    end 
    println() 
    println("")
    flush(stdout)

    # Multi-nodes case
    if length(MNcube)>0
        MNelems = Array{Int}(undef, total_frags)
        for i in 1:total_frags
            MNelems[i] = tasklist[i].fragid
        end 
        for i in 1:length(MNcube)
            for j in 1:length(MNcube[i][2])
                MNpos=findfirst(==(MNcube[i][2][j]),MNelems)
                #println("MNpos : ",MNpos)
                if sizeof(MNpos) > 0
                    #println("MNpos : ",MNpos, " i : " , i, " j : ", j)
                    #println("MNcube[i] : ", MNcube[i])
                    #println("MNcube[i][1] : ", MNcube[i][1])
                    ipos=MNpos[1]
                    tasklist[ipos].nnodes=MNcube[i][1][1]
                end 
            end 
        end
    end 

    println(LBcube)
    #exit()

end


@everywhere function NWChemRUN_SPAWNAT(tvec,tlist,id,snodes)
 
    println("taskvec  : ",tvec) 
    println("tasklist : ",tlist) 
    for itask in tvec
        if itask != 0
            hname = gethostname() 
            print("id",id,gethostname()," itask ",itask," ",(tlist[itask].infile)," ",(tlist[itask].outfile))

            hostlock0=tlist[itask].folder
            hostlock = string(hostlock0,"/","hostlock")
            tlist[itask].folder=string(tlist[itask].folder,"/",gethostname())
            println(" ",tlist[itask].folder)
            if !isdir(tlist[itask].folder)
                mkpath(tlist[itask].folder)
            end

            if isfile(hostlock)
                ilock = -1
                while ilock != 1
                    #lockcheck=open(hostlock,"r")
                    lines=readlines(hostlock)
                    if length(lines) > 0
                        for line in lines
                            if occursin(hname,line)
                                sleep(5)
                                println(" .. waiting for the jobs in this host (",hname,") ")
                            elseif length(snodes)-length(lines) < tlist[itask].nnodes
                                sleep(5)
                                println(" .. waiting for the jobs out of host (",hname,") ")
                            else
                                ilock = 1
                            end
                        end
                    #    close(lockcheck)
                    else
                        ilock = 1
                        println(" .. no hostlock, goon running ")
                    end
                end
            else
                println(" .. no hostlock, goon running ")
            end

            nmpi=5*tlist[itask].nnodes

            if isfile(hostlock)
                iolock=open(hostlock,"a+")
            else
                iolock=open(hostlock,"w")
            end

            lockvec  = []
            if tlist[itask].nnodes > 1
                hostflag = "-hostfile"
                hostfile = "nodes$(id)"

                ftmp = string(tlist[itask].folder,"/",hostfile)
                io=open(ftmp,"w")
                println(io, hname," slots=5")
                for inn in 1:tlist[itask].nnodes
                    if snodes[inn] != hname
                        println("inn ", inn, " snodes[inn] : ", snodes[inn])
                        println(io, snodes[inn]," slots=5")
                        println(iolock, snodes[inn])
                        push!(lockvec,snodes[inn])
                    end
                end
                close(io)
            else
                hostflag = " "
                hostfile = " "

                println(iolock, hname)
                push!(lockvec,hname)

            end

            close(iolock)

            println("hostflag : ", hostflag, " hostfile : ", hostfile)
            flush(stdout)

            RUNXX=`time ../../NWChemRUN $(nmpi) $hostflag $hostfile  $(tlist[itask].infile) $(tlist[itask].outfile) `

            try  
                run(`mv $(hostlock0)"/"$(tlist[itask].infile)  $(tlist[itask].folder)`)
                global IFDONE = false                
            catch err
                global IFDONE = true 
            end

            if IFDONE
                println("JOB $(tlist[itask].infile) already done in another Worker")                
            else 
                run(Cmd(RUNXX,dir=tlist[itask].folder,detach=true,ignorestatus=true))
            end 
            flush(stdout)
            flush(stderr)

            ccheck=string(tlist[itask].folder,"/",tlist[itask].outfile)
            rdlast=@spawn check_last_NWChem(ccheck,IOrecord,IFDONE)
            fetch(rdlast)

            for i in length(lockvec)
                println("lockvec[",i,"] : ", lockvec[i])
                run(`sed -i "/$(lockvec[i])/d" $(hostlock)`)
            end

        end 
    end 
    println("Done in NWChemRUN_SPAWNAT with inode-",id)
    flush(stdout)
end

@everywhere function NWChemRUN_SPAWN(tvec,tlist,id,snodes,iffifo)
 
    #println("taskvec  : ",tvec) 
    #println("tasklist : ",tlist) 
 
    tbreak=id
    while tbreak > 1
       tbreak = tbreak/10
    end
    sleep(tbreak)

    println(" iffifo : ", iffifo)
    flush(stdout)
    #exit(0)
  

    for itask in tvec
        if itask != 0
            hname = gethostname() 
            print("id ",id," ",gethostname()," itask ",itask," ",(tlist[itask].infile)," ",(tlist[itask].outfile))

            hostlock0=tlist[itask].folder
            hostlock = string(hostlock0,"/","hostlock")
            #tlist[itask].folder=string(tlist[itask].folder,"/",gethostname())
            folder = string(tlist[itask].folder,"/",gethostname())
            println(" ",folder) 
            if !isdir(folder)
                mkpath(folder)
            end 

            if isfile(hostlock)
                ilock = -1 
                while ilock != 1 
                    #lockcheck=open(hostlock,"r")                

                    global istamp = -1 
                    while istamp != 1
                        try
                            global lineslock=readlines(hostlock)
                            istamp = 1
                        catch err
                            sleep(0.1)
                            println("==> Re-try the hostlock read")
                        end
                    end

                    #lines=readlines(hostlock)
                    nlines=length(lineslock)

                    if length(lineslock) > 0
                        EXnodes = Array{String}(undef, nlines)
                        i=0
                        for line in lineslock
                            i=i+1
                            EXnodes[i] = line
                        end
                        ipos=findfirst(==(hname),EXnodes)
                        if sizeof(ipos)>0 
                            sleep(5)
                            println(" .. waiting for the jobs in this host (",hname,") ")               
                        #elseif length(snodes)-length(lineslock) < tlist[itask].nnodes 
                        #    sleep(5)
                        #    println(" .. waiting for the jobs out of host (",hname,") ")               
                        else
                            ilock = 1
                        end 

                    else
                        ilock = 1
                        println(" .. no host is locked, goon running ")
                    end  
                    #close(lockcheck) 
                end 
            else
                println(" .. no hostlock, goon running ")
            end 

            nmpi=5*tlist[itask].nnodes

            if isfile(hostlock)
                iolock=open(hostlock,"a+",lock = true)
                println("iolock was opened with a+")
            else
                iolock=open(hostlock,"w",lock = true)
                println("iolock was opened with w")
            end 
            flush(stdout)

            global lockvec  = [] 
            if tlist[itask].nnodes > 1
                hostflag = "-hostfile"
                hostfile = "nodes$(id)"

                ftmp = string(folder,"/",hostfile)
                io=open(ftmp,"w")
                println(io, hname," slots=5")

                ihfile = -1
                while ihfile != 1
                    hlockid=string(hostlock,"$(id)")

                    global idlock = -1
                    while idlock != 1                                         
                       try 
                          run(`cp $(hostlock) $(hlockid)`)                   
                          sleep(0.01) 
                          idlock = 1
                       catch err
                          sleep(0.1)
                          println("==> Re-try the hostlock read and cp as lockid")
                       end
                    end
                   
                    #lockcheck=open(hlockid,"r")
                    lines=readlines(hlockid)
                    nlines=length(lines) 
                    if nlines > 0
                        if length(snodes)-length(lines) < tlist[itask].nnodes
                            sleep(5)
                            println(" .. waiting for the jobs out of host (",hname,") for multi-nodes")
                        else                        
                            EXnodes = Array{String}(undef, nlines)
                            i=0
                            for line in lines
                                i=i+1  
                                EXnodes[i] = line
                            end

                            inodes=1
                            for inn in 1:length(snodes)
                                if snodes[inn] != hname
                                    ipos=findfirst(==(snodes[inn]),EXnodes)
                                    if sizeof(ipos) > 0
                                        println(">> snodes[",inn,"] : ", snodes[inn]," is not available in multi-nodes case ")                                                                            
                                    else 
                                        inodes=inodes+1
                                        println(">> snodes[",inn,"] : ", snodes[inn]," is available in the multi-nodes case ") 
                                        println(io, snodes[inn]," slots=5")
                                        println(iolock, snodes[inn])
                                        push!(lockvec,snodes[inn])
                                        if inodes >= tlist[itask].nnodes
                                            break
                                        end   
                                    end 
                                end 
                            end  
  
                            #println("Out of the for loop for multi-nodes, nline > 0 case")
                            println(iolock, hname)
                            push!(lockvec,hname)
                            ihfile = 1
                            #println("iolock was about to closed in multi-nodes, nline > 0 case")
                            flush(stdout)
                            try 
                                close(iolock)
                                println("iolock was closed in multi-nodes, nline > 0 case")
                            catch err
                                println("iolock was closed in multi-nodes (previously by break?)")             
                            end 
                            flush(stdout)

                        end                        
                    else
                        
                        inodes=1
                        for inn in 1:length(snodes)
                            if snodes[inn] != hname
                                inodes=inodes+1
                                println(">> snodes[",inn,"] : ", snodes[inn]," is available in the multi-nodes case ")
                                println(io, snodes[inn]," slots=5")
                                println(iolock, snodes[inn])
                                push!(lockvec,snodes[inn])
                                if inodes >= tlist[itask].nnodes
                                    break
                                end      
                            end
                        end

                        println(iolock, hname)
                        push!(lockvec,hname)
                        ihfile = 1 
                        println(" .. hosts are availavle, goon running ")
                        close(iolock)
                        println("iolock was closed in multi-nodes, nline = 0 case")
                        flush(stdout)

                    end
                    #close(lockcheck) 
                end
                close(io)
            else        
                hostflag = " "
                hostfile = " "
                println(iolock, hname)
                push!(lockvec,hname)
                try
                    close(iolock)
                    println("iolock was closed in single node case")
                catch err
                    println("iolock was closed in single node (previously by break?)")
                end
                sleep(0.1) 
                flush(stdout)
            end

            #supplement = filter!(ss->occursin(r"\.nw", ss),readdir(hostlock0))
            #println("000 hostlock0 : ", hostlock0, " 000 supplement : ",supplement)

            println("hostflag : ", hostflag, " hostfile : ", hostfile)
            flush(stdout)              

            try
                run(`mv $(hostlock0)"/"$(tlist[itask].infile)  $(folder)`)
                global IFDONE = false
                global infile = tlist[itask].infile
                global outfile = tlist[itask].outfile
            catch err
                global IFDONE = true
                global infile = tlist[itask].infile
                global outfile = tlist[itask].outfile
                if iffifo 
                   #println("when the fifo is activated, NWChem (.nw) for testing ")
                   supplement = filter!(ss->occursin(r"\.nw", ss),readdir(hostlock0)) 
                   println("hostlock0 : ", hostlock0, "  supplement : ",supplement) 
                   if length(supplement) > 0
                      try
                         run(`mv $(hostlock0)"/"$(supplement[1])  $(folder)`)
                         global IFDONE = false
                         global infile  = supplement[1]
                         global outfile = string(splitext(infile)[1],".out")
                         println("In FIFO, infile : ", infile, " outfile : ", outfile)
                      catch err
                         println("Conflicted by other worker")                         
                      end
                   end  
                end  
            end

            RUNXX=`time ../../NWChemRUN $(nmpi) $hostflag $hostfile  $infile $outfile `

            if IFDONE
                println("JOB $(infile) already done in another Worker")                
            else
                run(Cmd(RUNXX,dir=folder,detach=true,ignorestatus=true))
            end
            flush(stdout)
            flush(stderr)

            #ccheck=string(folder,"/",tlist[itask].outfile)
            ccheck=string(folder,"/",outfile)
            rdlast=@spawnat id check_last_NWChem(ccheck,IOrecord,IFDONE)
            fetch(rdlast)

            for i in 1:length(lockvec)
                println("lockvec[",i,"] : ", lockvec[i])
                global ised = -1
                while ised != 1
                    try 
                        run(`sed -i "/$(lockvec[i])/d" $(hostlock)`)
                        println("==> Shell's sed done")                
                        ised = 1
                    catch err
                        sleep(0.1) 
                        println("==> Re-try the shell's sed operations")                
                    end
                end   
            end 

        end 
    end 
    println("Done in NWChemRUN_SPAWN with inode-",id)
    flush(stdout)

end


function runtask(NNSLURM, IFSLURM, ISPAWN)

    println("")
    for i in 1:length(LBcube) 
        if i == NNSLURM || LBfile == ""
            println("Matched the number of NODES between SLURM and Task-Pre-Assignment ") 
            println("    NSLURM_NODES : ", NNSLURM, "        NPRE_NODES : ", i)  
            println("  ------------                ------------           ")  

            for j in 1:length(LBcube[i])
                println("  Tasks on the node-",j," :  ", LBcube[i][j])
            end  

            if length(MNcube) > 0
                println("  ------------                ------------           ")  
                for j in 1:length(MNcube)
                    println("  Multi-nodes task-ID :  ", MNcube[j][2], " with ",MNcube[j][1]," nodes ")
                end 
            end 

            println("  ------------                ------------           ")  
        end  
    end
    println("Running the QC tasks ... ")
    println("    Number of   cores : ", nprocs())
    println("    Number of workers : ", nworkers())
    if IFSLURM 
        SNODES = Array{String}(undef, NNSLURM)
        i=0 
        open("nodelist","r") do stream
            for line in eachline(stream)
                println(line)
                i = i + 1
                SNODES[i] = line
            end
        end 
    end 
    flush(stdout)

    for i in workers()
        id, pid, host = fetch(@spawnat i (myid(), getpid(), gethostname()))
        println( i, " ", id, " ", pid, " ", host, " done " )
    end
 
    @everywhere println("    process: $(myid()) on host $(gethostname())")
    println("")
    println("")
    flush(stdout)

    #println("LBcube : ",LBcube)
    #println("  ")
    println("pos_node : ", pos_node)
    println("  ")

    for i in 1:length(pos_node)
        if pos_node[i] == NNSLURM
            global ipos = i
        end 
    end  

    run=[]
    for i in 1:NNSLURM
        sleep(0.1)
        if IFSLURM
          inode = i+1
        else
          inode = i
        end
        
        LBtmp = [] 

        if LBfile != ""
            println(" LBcube[",ipos,"][",i,"] : ",LBcube[ipos][i])
            LBtmp=LBcube[ipos][i]
        else
            LBtmp=LBcube[1][i]
        end

        if ISPAWN == 1 
            push!(run,@spawn NWChemRUN_SPAWN(LBtmp,tasklist,inode,SNODES,iffifo))
        elseif ISPAWN == 2
            push!(run,@spawnat inode NWChemRUN_SPAWNAT(LBtmp,tasklist,inode,SNODES))
        end
 
    end 

    for i in 1:NNSLURM
        fetch(run[i])
    end  

end

 
