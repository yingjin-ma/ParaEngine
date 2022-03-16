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
        while true
            sleep(5)
            tailout=read_last(out)
            println(tailout)
            if occursin("Total times  cpu:",tailout)
                break
            end
        end
    end
    flush(stdout) 

end

function nwchemtask(task,frag,atoms,par)

    # Generating NWChem headers
   
    fragidx=lpad(frag.idx,8,"0")
    #run(`cp Template.nwchem $(fragidx).nw`) 
 
    infile=string(workdir,"/Frag-$(fragidx).nw")
    open(infile,"w") do nwinp
        txt = read("Template.nwchem",String)
        txt = replace(txt,"FRAG-i" => "Frag-$(fragidx)","X-library-basis" => "* library 6-31g","dftxc" => "xc m06-2x")
        print(nwinp,txt)
    end  

    inXYZ=string(workdir,"/XYZ.tmp")    
    open(inXYZ,"w") do xyztmp
        for iatom in atoms
            println(xyztmp,iatom.elem," ",iatom.coord[1]," ",iatom.coord[2]," ",iatom.coord[3])
        end
    end

    run(`sed -i "4r $inXYZ" $infile`)

    task.infile  = "Frag-$(fragidx).nw"
    task.outfile = "Frag-$(fragidx).out"

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

function distributingtask()

    nnodes = -1
    iflag  = -1
    icount =  0 

    println("")
    println("Distributing the QC tasks ... ")
    if LBfile != ""        
        print("Assigning tasks with --> ")
        print(LBfile)
        println(" <-- ")
        # Read in the loadbalance file 
        open(LBfile,"r") do LBread
            while !eof(LBread)
                line=readline(LBread)
                sline=split(line)
                ntmp=length(sline)
                if ntmp > 0                   
                    if uppercase(sline[1]) == "NODES"
                        nnodes = parse(Int32,sline[2])
                        LBmat=[] 
                        for i in 1:nnodes
                            line=readline(LBread)
                            sline=split(line)
                            ntmp=length(sline)
                            LBvec=[]
                            for i in 1:ntmp
                                push!(LBvec,parse(Int32,sline[i]))
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
                    end 
                end
            end
        end
    else
        println("Assigning tasks with --> Default <-- ")       
        LBmat=[]
        LBvec=[] 
        for i in 1:total_frags
            push!(LBvec,i)
        end 
        push!(LBmat,LBvec)
        push!(LBcube,LBmat) 
    end  
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
                println("MNpos : ",MNpos)
                if sizeof(MNpos) > 0
                    ipos=MNpos[1]
                    tasklist[ipos].nnodes=MNcube[i][1][j]
                end 
            end 
        end
    end 

    #println(tasklist)
    #stop()

end


@everywhere function NWChemRUN_SPAWNAT(tvec,tlist,id,snodes)
 
    println("taskvec  : ",tvec) 
    println("tasklist : ",tlist) 
    for itask in tvec
        if itask != 0
            print("id",id,gethostname()," itask ",itask," ",(tlist[itask].infile)," ",(tlist[itask].outfile))
            RUNXX=`time ../../NWChemRUN 5 $(tlist[itask].infile) $(tlist[itask].outfile) `
            tlist[itask].folder=string(tlist[itask].folder,"/",gethostname())
            println(" ",tlist[itask].folder) 
            if !isdir(tlist[itask].folder)
                mkpath(tlist[itask].folder)
            end

            try  
                run(`mv $(workdir)"/"$(tlist[itask].infile)  $(tlist[itask].folder)`)
                global IFDONE = false                
            catch err
                global IFDONE = true 
            end

            if IFDONE
                println("JOB $(tlist[itask].infile) already done in another Worker")                
            else 
                run(Cmd(RUNXX,dir=tlist[itask].folder,detach=true))
            end 
            flush(stdout)
            flush(stderr)

            ccheck=string(tlist[itask].folder,"/",tlist[itask].outfile)
            rdlast=@spawnat id check_last_NWChem(ccheck,IOrecord,IFDONE)
            fetch(rdlast)

        end 
    end 
    println("Done in NWChemRUN_SPAWNAT with inode-",id)
    flush(stdout)
end

@everywhere function NWChemRUN_SPAWN(tvec,tlist,id,snodes)
 
    println("taskvec  : ",tvec) 
    println("tasklist : ",tlist) 
    for itask in tvec
        if itask != 0
            hname = gethostname() 
            print("id",id,gethostname()," itask ",itask," ",(tlist[itask].infile)," ",(tlist[itask].outfile))

            tlist[itask].folder=string(tlist[itask].folder,"/",gethostname())
            println(" ",tlist[itask].folder) 
            if !isdir(tlist[itask].folder)
                mkpath(tlist[itask].folder)
            end 

            hostflag = " "
            hostfile = " "
            nmpi=5*tlist[itask].nnodes
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
                        end 
                    end 
                close(io)
            end
            println("hostflag : ", hostflag, " hostfile : ", hostfile)
            flush(stdout)

            RUNXX=`time ../../NWChemRUN $(nmpi) $hostflag $hostfile  $(tlist[itask].infile) $(tlist[itask].outfile) `

            try
                run(`mv $(workdir)"/"$(tlist[itask].infile)  $(tlist[itask].folder)`)
                global IFDONE = false
            catch err
                global IFDONE = true
            end

            if IFDONE
                println("JOB $(tlist[itask].infile) already done in another Worker")                
            else
                run(Cmd(RUNXX,dir=tlist[itask].folder,detach=true))
            end
            flush(stdout)
            flush(stderr)

            ccheck=string(tlist[itask].folder,"/",tlist[itask].outfile)
            rdlast=@spawn check_last_NWChem(ccheck,IOrecord,IFDONE)
            fetch(rdlast)

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

            if LBfile != ""
                for j in 1:length(LBcube[i])
                    println("  Tasks on the node-",j," :  ", LBcube[i][j])
                end  
            else
                println("  Tasks on ALL-nodes :  ", LBcube[1][1])
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

    run=[]
    for i in 1:NNSLURM
        if IFSLURM
          inode = i+1
        else
          inode = i
        end
        
        if LBfile != ""
            LBtmp=LBcube[NNSLURM][i]
        else
            LBtmp=LBcube[1][1]
        end   
 
        if ISPAWN == 1 
            push!(run,@spawn NWChemRUN_SPAWN(LBtmp,tasklist,inode,SNODES))
        elseif ISPAWN == 2
            push!(run,@spawnat inode NWChemRUN_SPAWNAT(LBtmp,tasklist,inode,SNODES))
        end 
    end 

    for i in 1:NNSLURM
        fetch(run[i])
    end  

end

 
